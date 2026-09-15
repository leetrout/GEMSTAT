#include <gsl/gsl_math.h>

#include "ExprPredictor.h"
#include "ExprPar.h"

#include <map>
#include <sstream>
#include <stdexcept>

//#define DEBUG

/*****************************************************
 * ParamSlots
 ******************************************************/

void ParamSlots::flatten_paths( const gsparams::DictList& d, const std::string& prefix, vector< std::string >& out )
{
    if ( d.my_type == gsparams::undecided ) return;
    if ( d.my_type == gsparams::primitive ) { out.push_back( prefix ); return; }
    int n = d.list_storage.size();
    for ( int i = 0; i < n; i++ )
    {
        std::string key;
        if ( d.my_type == gsparams::dict ) key = d.map_key_storage.at(i);
        else { std::ostringstream ss; ss << i; key = ss.str(); }
        flatten_paths( d.list_storage.at(i), prefix.empty() ? key : prefix + "/" + key, out );
    }
}

ParamSlots ParamSlots::build( const ExprPar& par, const vector< string >& motifNames )
{
    vector< std::string > paths;
    flatten_paths( par.my_pars, "", paths );
    std::map< std::string, int > index;
    for ( size_t i = 0; i < paths.size(); i++ ) index[ paths[i] ] = (int)i;
    struct Find {
        const std::map< std::string, int >& index;
        int operator()( const std::string& p ) const { std::map< std::string, int >::const_iterator it = index.find( p ); return it == index.end() ? -1 : it->second; }
    } find = { index };

    ParamSlots s;
    s.n_pars = paths.size();
    int nF = motifNames.size();
    s.maxbind.resize( nF ); s.alpha_a.resize( nF ); s.alpha_r.resize( nF );
    for ( int f = 0; f < nF; f++ )
    {
        s.maxbind[f] = find( "tfs/" + motifNames[f] + "/maxbind" );
        s.alpha_a[f] = find( "tfs/" + motifNames[f] + "/alpha_a" );
        s.alpha_r[f] = find( "tfs/" + motifNames[f] + "/alpha_r" );
    }
    s.inter.assign( nF, vector< int >( nF, -1 ) );
    for ( int a = 0; a < nF; a++ )
    {
        for ( int b = a; b < nF; b++ )
        {
            int slot = find( "inter/" + motifNames[a] + ":" + motifNames[b] );
            if ( slot < 0 ) slot = find( "inter/" + motifNames[b] + ":" + motifNames[a] );
            s.inter[a][b] = s.inter[b][a] = slot;
        }
    }
    for ( int k = 0; ; k++ )
    {
        std::ostringstream ss; ss << "qbtm/" << k;
        int slot = find( ss.str() );
        if ( slot < 0 ) break;
        s.qbtm.push_back( slot );
    }
    for ( int k = 0; ; k++ )
    {
        std::ostringstream ss; ss << "enh/" << k;
        int beta = find( ss.str() + "/beta" );
        if ( beta < 0 ) break;
        s.beta.push_back( beta );
        s.pi.push_back( find( ss.str() + "/pi" ) );
    }
    return s;
}

/*****************************************************
 * ExprFunc
 ******************************************************/

ExprFunc::ExprFunc( const ExprModel* _model, const ExprPar& _par , const SiteVec& sites_, const int seq_len, const int seq_num): expr_model(_model), motifs( _model->motifs ), actIndicators( _model->actIndicators ), maxContact( _model->maxContact ), repIndicators( _model->repIndicators ), repressionMat( _model->repressionMat ), repressionDistThr( _model->repressionDistThr ), par(_par), factorIntMat(_model->motifs.size(),_model->motifs.size(),1.0)
{
    //par = _par;//NOTE: made this const, and that solved a memory leak.

    int nFactors = par.nFactors();
    if(motifs.size() != nFactors){
        //std::cerr << " MOTFIS " << motifs.size() << " : factors par " << nFactors << std::endl;
        throw std::runtime_error(" Motifs and factors did not match");
    }
    assert( motifs.size() == nFactors );
    assert( actIndicators.size() == nFactors );
    assert( repIndicators.size() == nFactors );
    assert( repressionMat.isSquare() && repressionMat.nRows() == nFactors );
    assert( maxContact >= 0 );

    n_sites = sites_.size();//number of non-psudo sites.
    seq_number = seq_num;
    seq_length = seq_len;

    //setup legacy parameters
    maxBindingWts.clear();
    maxBindingWts.assign(nFactors,1.0);
    //maxBindingWts = vector < GEMSTAT_PAR_FLOAT_T >(nFactors);          // binding weight of the strongest site for each TF: K(S_max) [TF_max]
    txpEffects.clear();
    txpEffects.resize(nFactors,1.0);
    repEffects.clear();
    repEffects.resize(nFactors,1.0);

    std::map<std::string, int> tf_names_to_ids;
    for(int i = 0;i<expr_model->motifnames.size();i++){
        tf_names_to_ids[expr_model->motifnames.at(i)] = i;
    }

    //Do not assume that the tfs dictionary is in our internal order.
    for(int i = 0;i<nFactors;i++){
        const std::string &which_tf = expr_model->motifnames.at(i);
        maxBindingWts[i] = ((gsparams::DictList&)par.my_pars)["tfs"][which_tf]["maxbind"] ;
        txpEffects[i] = ((gsparams::DictList&)par.my_pars)["tfs"][which_tf]["alpha_a"] ;
        repEffects[i] = ((gsparams::DictList&)par.my_pars)["tfs"][which_tf]["alpha_r"] ;
    }

    //factorIntMat = Matrix(nFactors, nFactors);                      // (maximum) interactions between pairs of factors: omega(f,f')
    //Order doesn't matter for factor interactions, positions will be looked up.
    for(int k = 0;k<((gsparams::DictList&)par.my_pars)["inter"].size();k++){
        //need to split the name.
        std::string key = ((gsparams::DictList&)par.my_pars)["inter"].map_key_storage.at(k);
        double value = ((gsparams::DictList&)par.my_pars)["inter"].list_storage.at(k);

        int split_pos = key.find(":");
        int i = tf_names_to_ids[key.substr(0,split_pos)];//TODO: more defensive here. Make sure the value actually existed.
        int j = tf_names_to_ids[key.substr(split_pos+1,key.size())];
        factorIntMat.setElement(i,j,value);
        factorIntMat.setElement(j,i,value);
    }


    this->setupSitesAndBoundaries(sites_,seq_length, seq_num);
    this->buildLeftNeighbours( true, left_nbrs );

    // the values the recurrences work on for an ordinary prediction
    plain.maxBindingWts.assign( maxBindingWts.begin(), maxBindingWts.end() );
    plain.txpEffects.assign( txpEffects.begin(), txpEffects.end() );
    plain.repEffects.assign( repEffects.begin(), repEffects.end() );
    plain.basal = par.getPromoterData( seq_number ).basal_trans;
    this->fillPlainWeights();
}

void ExprFunc::buildLeftNeighbours( bool within_boundaries, vector< vector< SiteInteraction > >& out ) const
{
    int n = n_sites;
    out.assign( n + 1, vector< SiteInteraction >() );
    for ( int i = 1; i <= n; i++ )
    {
        int j_begin = within_boundaries ? boundaries[i] + 1 : 1;
        for ( int j = j_begin; j < i; j++ )
        {
            if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;
            SiteInteraction si;
            si.j = j;
            si.dist = sites[i].start - sites[j].start;
            si.w_ij = compFactorInt( sites[ i ], sites[ j ] );
            si.w_ji = compFactorInt( sites[ j ], sites[ i ] );
            si.rep_ij = testRepression( sites[ i ], sites[ j ] );
            si.rep_ji = testRepression( sites[ j ], sites[ i ] );
            out[i].push_back( si );
        }
    }
}

void ExprFunc::fillPlainWeights()
{
    plain.w_ij.assign( left_nbrs.size(), vector< gemstat_dp_t >() );
    plain.w_ji.assign( left_nbrs.size(), vector< gemstat_dp_t >() );
    for ( size_t i = 0; i < left_nbrs.size(); i++ )
    {
        plain.w_ij[i].resize( left_nbrs[i].size() );
        plain.w_ji[i].resize( left_nbrs[i].size() );
        for ( size_t k = 0; k < left_nbrs[i].size(); k++ )
        {
            plain.w_ij[i][k] = left_nbrs[i][k].w_ij;
            plain.w_ji[i][k] = left_nbrs[i][k].w_ji;
        }
    }
    plain.all_w_ji.assign( all_left_nbrs.size(), vector< gemstat_dp_t >() );
    for ( size_t i = 0; i < all_left_nbrs.size(); i++ )
    {
        plain.all_w_ji[i].resize( all_left_nbrs[i].size() );
        for ( size_t k = 0; k < all_left_nbrs[i].size(); k++ ) plain.all_w_ji[i][k] = all_left_nbrs[i][k].w_ji;
    }
    plain.right_w_ij.assign( right_nbrs.size(), vector< gemstat_dp_t >() );
    for ( size_t i = 0; i < right_nbrs.size(); i++ )
    {
        plain.right_w_ij[i].resize( right_nbrs[i].size() );
        for ( size_t k = 0; k < right_nbrs[i].size(); k++ ) plain.right_w_ij[i][k] = right_nbrs[i][k].w_ij;
    }
}

ThermoVals< gemstat_ad_t > ExprFunc::makeVals( const vector< gemstat_ad_t >& flat, const ParamSlots& slots ) const
{
    typedef gemstat_ad_t V;
    // The Logistic model works on ENERGY_SPACE parameters (see
    // ExprModel::createNewExprFunc); the flat vector is PROB_SPACE.
    const bool energy = ( expr_model->modelOption == LOGISTIC );
    struct Get {
        const vector< V >& flat; bool energy;
        V operator()( int slot, double dflt ) const
        {
            if ( slot < 0 ) return V( (gemstat_dp_t)dflt );
            return energy ? log( flat[slot] ) : flat[slot];
        }
    } get = { flat, energy };

    ThermoVals< V > v;
    int nF = motifs.size();
    v.maxBindingWts.resize( nF ); v.txpEffects.resize( nF ); v.repEffects.resize( nF );
    for ( int f = 0; f < nF; f++ )
    {
        v.maxBindingWts[f] = get( slots.maxbind[f], 1.0 );
        v.txpEffects[f] = get( slots.alpha_a[f], 1.0 );
        v.repEffects[f] = get( slots.alpha_r[f], 1.0 );
    }

    v.w_ij.assign( left_nbrs.size(), vector< V >() );
    v.w_ji.assign( left_nbrs.size(), vector< V >() );
    for ( size_t i = 0; i < left_nbrs.size(); i++ )
    {
        v.w_ij[i].resize( left_nbrs[i].size() );
        v.w_ji[i].resize( left_nbrs[i].size() );
        for ( size_t k = 0; k < left_nbrs[i].size(); k++ )
        {
            int j = left_nbrs[i][k].j;
            V normalInt = get( slots.inter[ sites[i].factorIdx ][ sites[j].factorIdx ], 1.0 );
            v.w_ij[i][k] = factorIntAffine( sites[i], sites[j], normalInt );
            v.w_ji[i][k] = factorIntAffine( sites[j], sites[i], normalInt );
        }
    }
    v.all_w_ji.assign( all_left_nbrs.size(), vector< V >() );
    for ( size_t i = 0; i < all_left_nbrs.size(); i++ )
    {
        v.all_w_ji[i].resize( all_left_nbrs[i].size() );
        for ( size_t k = 0; k < all_left_nbrs[i].size(); k++ )
        {
            int j = all_left_nbrs[i][k].j;
            V normalInt = get( slots.inter[ sites[i].factorIdx ][ sites[j].factorIdx ], 1.0 );
            v.all_w_ji[i][k] = factorIntAffine( sites[j], sites[i], normalInt );
        }
    }
    v.right_w_ij.assign( right_nbrs.size(), vector< V >() );
    for ( size_t i = 0; i < right_nbrs.size(); i++ )
    {
        v.right_w_ij[i].resize( right_nbrs[i].size() );
        for ( size_t k = 0; k < right_nbrs[i].size(); k++ )
        {
            int j = right_nbrs[i][k].j;
            V normalInt = get( slots.inter[ sites[i].factorIdx ][ sites[j].factorIdx ], 1.0 );
            v.right_w_ij[i][k] = factorIntAffine( sites[i], sites[j], normalInt );
        }
    }

    // basal transcription slot, as ExprPar::getPromoterData() picks it
    int use_enhancerID = expr_model->shared_scaling ? 0 : seq_number;
    int use_basal = expr_model->one_qbtm_per_crm ? use_enhancerID : 0;
    if ( use_basal >= (int)slots.qbtm.size() ) throw std::logic_error( "ExprFunc::makeVals: no qbtm slot for this sequence" );
    v.basal = get( slots.qbtm[ use_basal ], 1.0 );
    return v;
}

void ExprFunc::setupSitesAndBoundaries(const SiteVec& _sites, int length, int seq_num){
  #ifdef DEBUG
  cerr << "running ExprFunc::setupSitesAndBoundaries(...)" << endl;
  #endif

	
  int n = _sites.size();
  sites = SiteVec(_sites);

	//Prevents crashes if there are no annotated sites in the sequence.
	//This can happen, and it should not cause an error.
	int pseudo_start = -1000;
	int pseudo_end = length+1000;
	if(sites.size() > 0){
		pseudo_start = _sites[0].start-1000;
		pseudo_end = _sites[_sites.size()-1].start+1000;
	}

	sites.insert( sites.begin(), Site(pseudo_start,pseudo_start,true,-1,0.0,1.0) );        // start with a pseudo-site at position 0
	sites.push_back( Site(pseudo_end,pseudo_end,true,-1,0.0,1.0) );       //and another pseudo-site at the end

  // store the sequence

  boundaries.resize(n_sites+2);
  boundaries[0] = 0;//value for starting pseudosite
  double range = max( (double)expr_model->get_longest_coop_thr(), (double)repressionDistThr );
  for ( int i = 1; i <= n; i++ )
  {
      int j;
      for ( j = i - 1; j >= 1; j-- )
      {
          if ( ( sites[i].start - sites[j].start ) > range ) break;
      }//If the loop never broke, j will leak with value 0.
      int boundary = j;
      boundaries[i] = ( boundary );
  }
  boundaries[n_sites+1] = n_sites;//last-non-psudoe-site position boundary for ending pseudosite?

}

void ExprFunc::setupBindingWeights(const vector< double >& factorConcs){
  vector< gemstat_dp_t >& bindingWts = plain.bindingWts;
  bindingWts.resize(sites.size());
  bindingWts[0] = 1.0;  //for first pseudosite
  bindingWts[bindingWts.size()-1] = 1.0; //might or might not be a pseudosite, gets overwritten if not
  int n = n_sites;
  for ( int i = 1; i <= n; i++ )
  {
      bindingWts[i] = maxBindingWts[ sites[i].factorIdx ] * factorConcs[sites[i].factorIdx] * sites[i].prior_probability * sites[i].wtRatio ;
  }
}

double ExprFunc::predictExpr(const Condition& in_condition){
  double to_return = this->predictExpr( in_condition.concs );
  if(to_return < 0.0 || to_return != to_return){
	cerr << "returning a nonsense expression prediction " << to_return << endl;
	exit(1);
  }
  assert(to_return >= 0.0); //Positive expression
  assert(!(to_return != to_return)); //not NaN
  return to_return;
}

double ExprFunc::predictExpr( const vector< double >& factorConcs )
{
    // compute the Boltzman weights of binding for all sites
    setupBindingWeights(factorConcs);

    // Thermodynamic models: Direct, Quenching, ChrMod_Unlimited and ChrMod_Limited
    // compute the partition functions
    gemstat_dp_t Z_off = compPartFuncOff();
    gemstat_dp_t Z_on = compPartFuncOn();

    // compute the expression (promoter occupancy)
    gemstat_dp_t promoterOcc = thermoOccupancy( Z_off, Z_on, plain.basal );
    #ifdef DEBUG
    if(promoterOcc < 0.0 || promoterOcc != promoterOcc){
	cerr << "Ridiculous in Direct!" << endl;
	cerr << "Z_on" << Z_on << endl;
	cerr << "Z_off" << Z_off << endl;
	cerr << "basal " << plain.basal << endl;
	cerr << "=====" << endl;
	}
    #else
    assert(promoterOcc >= 0.0 && promoterOcc == promoterOcc);
    #endif
    return promoterOcc;
}

gemstat_ad_t ExprFunc::predictExprAD( ThermoVals< gemstat_ad_t >& vals, const vector< double >& factorConcs ) const
{
    setupBindingWeightsT( vals, factorConcs );
    gemstat_ad_t Z_off = compPartFuncOffAD( vals );
    gemstat_ad_t Z_on = compPartFuncOnAD( vals );
    return thermoOccupancy( Z_off, Z_on, vals.basal );
}

/*****************************************************
 * Logistic
 ******************************************************/

double Logistic_ExprFunc::predictExpr( const vector< double >& factorConcs ){
  setupBindingWeights(factorConcs);
  return kernelLogistic( plain );
}

gemstat_ad_t Logistic_ExprFunc::predictExprAD( ThermoVals< gemstat_ad_t >& vals, const vector< double >& factorConcs ) const
{
    setupBindingWeightsT( vals, factorConcs );
    return kernelLogistic( vals );
}

/*****************************************************
 * Markov
 ******************************************************/

Markov_ExprFunc::Markov_ExprFunc( const ExprModel* _model, const ExprPar& _par , const SiteVec& sites_, const int seq_len, const int seq_num) : ExprFunc( _model, _par , sites_, seq_len, seq_num){

    //Additional setup for a Markov_ExprFunc
    #ifdef DEBUG
    cerr << "running Markov_ExprFunc::setupSitesAndBoundaries(...)" << endl;
    #endif

    double range = max( (double)expr_model->get_longest_coop_thr(), (double)repressionDistThr );
    rev_bounds.resize(n_sites+2);
    rev_bounds[n_sites] = n_sites+1; //last true site points to pseudosite
    rev_bounds[0] = 1; // first pseudosite has reverse boundary pointing to first true site
    for ( int i = n_sites; i > 0; i-- )
    {
        int j;
        for ( j = i+1; j <= n_sites; j++ )
        {
            if ( ( sites[j].start - sites[i].start) > range ) break;
        }//If the loop never broke, j will leak with value n+1
        int boundary = j;
        rev_bounds[i] = ( boundary );
    }
    assert(rev_bounds[n_sites] == n_sites+1);

    right_nbrs.assign( n_sites + 1, vector< SiteInteraction >() );
    for ( int i = 1; i <= n_sites; i++ )
    {
        for ( int j = rev_bounds[i] - 1; j > i; j-- )
        {
            if ( siteOverlap( sites[ i ], sites[ j ], motifs ) ) continue;
            SiteInteraction si;
            si.j = j;
            si.dist = sites[j].start - sites[i].start;
            si.w_ij = compFactorInt( sites[ i ], sites[ j ] );
            si.w_ji = compFactorInt( sites[ j ], sites[ i ] );
            si.rep_ij = testRepression( sites[ i ], sites[ j ] );
            si.rep_ji = testRepression( sites[ j ], sites[ i ] );
            right_nbrs[i].push_back( si );
        }
    }
    this->fillPlainWeights();
}

double Markov_ExprFunc::predictExpr( const vector< double >& factorConcs )
{
    setupBindingWeights(factorConcs);
    return kernelMarkov( plain );
}

gemstat_ad_t Markov_ExprFunc::predictExprAD( ThermoVals< gemstat_ad_t >& vals, const vector< double >& factorConcs ) const
{
    setupBindingWeightsT( vals, factorConcs );
    return kernelMarkov( vals );
}

/*****************************************************
 * Partition functions
 ******************************************************/

gemstat_dp_t ExprFunc::compPartFuncOff() const { return kernelOffBasic( plain ); }
gemstat_ad_t ExprFunc::compPartFuncOffAD( const ThermoVals< gemstat_ad_t >& v ) const { return kernelOffBasic( v ); }

gemstat_dp_t ExprFunc::compPartFuncOn() const
{
    throw std::logic_error( "ExprFunc::compPartFuncOn: no on-state partition function for this model" );
}
gemstat_ad_t ExprFunc::compPartFuncOnAD( const ThermoVals< gemstat_ad_t >& v ) const
{
    throw std::logic_error( "ExprFunc::compPartFuncOnAD: no on-state partition function for this model" );
}

gemstat_dp_t Direct_ExprFunc::compPartFuncOn() const { return kernelOnDirect( plain ); }
gemstat_ad_t Direct_ExprFunc::compPartFuncOnAD( const ThermoVals< gemstat_ad_t >& v ) const { return kernelOnDirect( v ); }

Quenching_ExprFunc::Quenching_ExprFunc( const ExprModel* _model, const ExprPar& _par , const SiteVec& sites_, const int seq_len, const int seq_num) : ExprFunc( _model, _par , sites_, seq_len, seq_num)
{
    buildLeftNeighbours( false, all_left_nbrs );   // the on-state recurrence is not window-limited
    this->fillPlainWeights();
}
gemstat_dp_t Quenching_ExprFunc::compPartFuncOn() const { return kernelOnQuenching( plain ); }
gemstat_ad_t Quenching_ExprFunc::compPartFuncOnAD( const ThermoVals< gemstat_ad_t >& v ) const { return kernelOnQuenching( v ); }

gemstat_dp_t ChrMod_ExprFunc::compPartFuncOff() const { return kernelOffChrMod( plain ); }
gemstat_ad_t ChrMod_ExprFunc::compPartFuncOffAD( const ThermoVals< gemstat_ad_t >& v ) const { return kernelOffChrMod( v ); }

gemstat_dp_t ChrModUnlimited_ExprFunc::compPartFuncOn() const { return kernelOnChrModUnlimited( plain ); }
gemstat_ad_t ChrModUnlimited_ExprFunc::compPartFuncOnAD( const ThermoVals< gemstat_ad_t >& v ) const { return kernelOnChrModUnlimited( v ); }

gemstat_dp_t ChrModLimited_ExprFunc::compPartFuncOn() const { return kernelOnChrModLimited( plain ); }
gemstat_ad_t ChrModLimited_ExprFunc::compPartFuncOnAD( const ThermoVals< gemstat_ad_t >& v ) const { return kernelOnChrModLimited( v ); }

/*****************************************************
 * Pairwise terms
 ******************************************************/

double ExprFunc::compFactorInt( const Site& a, const Site& b ) const
{
    // 	assert( !siteOverlap( a, b, motifs ) );
    double maxInt = factorIntMat( a.factorIdx, b.factorIdx );
    double dist = abs( b.start - a.start );
    //bool orientation = ( a.strand == b.strand );

    FactorIntFunc* an_int_func = expr_model->coop_setup->coop_func_for(a.factorIdx, b.factorIdx);
    return an_int_func->compFactorInt( maxInt, dist, a.strand, b.strand );
}


bool ExprFunc::testRepression( const Site& a, const Site& b ) const
{
    // 	assert( !siteOverlap( a, b, motifs ) );

    double dist = abs( b.start - a.start  );
    return repressionMat( a.factorIdx, b.factorIdx ) && ( dist <= repressionDistThr );
}
