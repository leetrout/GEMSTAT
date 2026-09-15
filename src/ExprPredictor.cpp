#include <typeinfo>


#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>
#include <stdexcept>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <nlopt.hpp>

#include "ExprPredictor.h"
#include "ExprPar.h"
#include "ExprFunc.h"

double nlopt_obj_func( const vector<double> &x, vector<double> &grad, void* f_data);

/**
 * Run an NLopt optimizer and make sure the caller gets a consistent result.
 *
 * NLopt reports some stopping conditions as C++ exceptions.  roundoff_limited
 * means it could not make further progress because of floating point
 * precision; the point reached is still useful, so we keep it and only note
 * the condition.  Any other exception is printed with its reason instead of
 * being discarded.  In every case free_pars is whatever NLopt left in it and
 * obj_result is the objective evaluated at that point, so a stale value from a
 * previous stage is never returned.
 */
static void run_nlopt( nlopt::opt& optimizer, vector<double>& free_pars, double& obj_result, const char* stage, ExprPredictor* predictor )
{
    try{
        optimizer.optimize( free_pars, obj_result );
        obj_result = optimizer.last_optimum_value();
        return;
    }catch( const nlopt::roundoff_limited& ){
        cerr << stage << ": NLopt stopped because roundoff errors limited further progress; keeping the point reached." << endl;
    }catch( const std::exception& e ){
        cerr << stage << ": NLopt stopped abnormally: " << e.what() << endl;
    }
    vector<double> no_grad;
    obj_result = nlopt_obj_func( free_pars, no_grad, predictor );
}

ExprPredictor::ExprPredictor( const vector <Sequence>& _seqs, const vector< SiteVec >& _seqSites, const vector< int >& _seqLengths, TrainingDataset* _training_data, const vector< Motif >& _motifs, const ExprModel& _expr_model,
		const vector < bool >& _indicator_bool, const vector <string>& _motifNames) : TrainingAware(), seqs(_seqs), seqSites( _seqSites ), seqLengths( _seqLengths ), training_data( _training_data ),
	expr_model( _expr_model),
	indicator_bool ( _indicator_bool ), motifNames ( _motifNames ),
	search_option(UNCONSTRAINED)
{
    //TODO: Move appropriate lines from this block to the ExprModel class.
	cerr << "exprData size: " << training_data->n_rows_output() << "  " << nSeqs() << endl;
    assert( training_data->n_rows_output() == nSeqs() );
    //assert( training_data->factorExprData.nRows() == nFactors() && training_data->factorExprData.nCols() == nConds() );
    //assert( expr_model.coopMat.isSquare() && expr_model.coopMat.isSymmetric() && expr_model.coopMat.nRows() == nFactors() );
    assert( expr_model.actIndicators.size() == nFactors() );
    assert( expr_model.maxContact > 0 );
    assert( expr_model.repIndicators.size() == nFactors() );
    assert( expr_model.repressionMat.isSquare() && expr_model.repressionMat.nRows() == nFactors() );
    assert( expr_model.repressionDistThr >= 0 );

	//****** DEFAULT VALUES *********
	objOption = SSE;

	n_alternations = 4;
	n_random_starts = 5;


	max_simplex_iterations = 200;
	max_gradient_iterations = 50;


    //gene_crm_fout.open( "gene_crm_fout.txt" );

    // set the model option for ExprPar and ExprFunc
    //ExprPar::modelOption = expr_model.modelOption;//TODO: Remove both of these.
    //ExprFunc::modelOption = expr_model.modelOption;

    // set the values of the parameter range according to the model option
    if ( expr_model.modelOption != LOGISTIC && expr_model.modelOption != DIRECT )
    {
        //ExprPar::min_effect_Thermo = 0.99;
        //ExprPar::min_interaction = 0.99;
    }

    //expr_model was already initialized. Setup the parameter factory.
    param_factory = new ParFactory(expr_model, nSeqs());

	trainingObjective = NULL;

	gradient_method = GRADIENT_AD;
	par_index = param_factory->create_index_par();

	maxShift = 5;
	shiftPenalty = 0.8;
	min_delta_f_SSE = 1.0E-8;
	min_delta_f_Corr = 1.0E-8;
	min_delta_f_CrossCorr = 1.0E-8;
	min_delta_f_PGP = 1.0E-8;

	set_objective_option(objOption);

    /* DEBUG
    cout << setprecision(10);
    ExprPar foo = param_factory->createDefaultMinMax(true);
    cout << " MAXIMUMS " << endl;
    printPar(foo);
    foo = param_factory->createDefaultMinMax(false);
    cout << " MINIMUMS " << endl;
    printPar(foo);
    */
}

ExprPredictor::~ExprPredictor()
{
  delete param_factory;
  delete trainingObjective;
  delete training_data;
}

void ExprPredictor::set_objective_option( ObjType in_obj_option ){
	//TODO: Move this to the front-end or something?
    //Maybe make it have a default SSE score objective, but anything else gets specified in the front-end.
	objOption = in_obj_option;

	if(NULL != trainingObjective){
		delete trainingObjective;
		trainingObjective = NULL;
	}

    switch(in_obj_option){
      case CORR:
        trainingObjective = new AvgCorrObjFunc();
        break;
      case PGP:
        trainingObjective = new PGPObjFunc();
        break;
      case CROSS_CORR:
        trainingObjective = new AvgCrossCorrObjFunc(maxShift, shiftPenalty);
        break;
      case LOGISTIC_REGRESSION:
        trainingObjective = new LogisticRegressionObjFunc();
        break;
	case PEAK_WEIGHTED:
		trainingObjective = new PeakWeightedObjFunc();
		break;
      case SSE:
      default:
        trainingObjective = new RMSEObjFunc();
        break;
    }
}

double ExprPredictor::objFunc( const ExprPar& par )
{
    double objective_value = evalObjective( par );

    return objective_value;
}


int ExprPredictor::train( const ExprPar& par_init )
{
    par_model = par_init;

    cout << "*** Diagnostic printing BEFORE adjust() ***" << endl;
    cout << "Parameters: " << endl;
    printPar( par_model );
    cout << endl;
    cout << "Objective function value: " << objFunc( par_model ) << endl;
    cout << "*******************************************" << endl << endl;

    if ( n_alternations > 0 && this->search_option == CONSTRAINED ){
      par_model = param_factory->truncateToBounds(par_model, indicator_bool);

    }
    obj_model = objFunc( par_model );

    cout << "*** Diagnostic printing AFTER adjust() ***" << endl;
    cout << "Parameters: " << endl;
    printPar( par_model );
    cout << endl;
    cout << "Objective function value: " << objFunc( par_model ) << endl;
    cout << "*******************************************" << endl << endl;

    if ( n_alternations == 0 ) return 0;

    // alternate between two different methods
    ExprPar par_result = param_factory->create_expr_par();
    double obj_result;
	this->start_training();

    for ( int i = 1; i <= n_alternations; i++ )
    {
		this->begin_epoch(i);
        simplex_minimize( par_result, obj_result );
        par_model = par_result;

        gradient_minimize( par_result, obj_result );
        par_model = par_result;
    }

    #ifdef BETAOPTBROKEN
    optimize_beta( par_model, obj_result );
    #endif

	this->end_training();

    // commit the parameters and the value of the objective function
    //par_model = par_result;
    obj_model = obj_result;

    return 0;
}


int ExprPredictor::train( const ExprPar& par_init, const gsl_rng* rng )
{
    /*
        //for random starts:
        ExprPar par_rand_start = par_init;
        par_rand_start = param_factor->randSamplePar( rng );
        train( par_rand_start );*/
    // training using the initial values
    train( par_init );

    cout << "Initial training:\tParameters = "; printPar( par_model );
    cout << "\tObjective = " << setprecision( 5 ) << obj_model << endl;

    // training with random starts
    ExprPar par_best = par_model;
    double obj_best = obj_model;
    for ( int i = 0; i < n_random_starts; i++ )
    {
        ExprPar par_curr = par_init;
        par_curr = param_factory->randSamplePar( rng );
        train( par_curr );
        cout << "Random start " << i + 1 << ":\tParameters = "; printPar( par_model );
        cout << "\tObjective = " << setprecision( 5 ) << obj_model << endl;
        if ( obj_model < obj_best )
        {
            par_best = par_model;
            obj_best = obj_model;
        }
    }

    // training using the best parameters so far
    if ( n_random_starts ) train( par_best );
    cout << "Final training:\tParameters = "; printPar( par_model );
    cout << "\tObjective = " << setprecision( 5 ) << obj_model << endl;

    //gene_crm_fout.close();

    return 0;
}


int ExprPredictor::train()
{
    // random number generator
    gsl_rng* rng;
    gsl_rng_env_setup();
    const gsl_rng_type * T = gsl_rng_default;     // create rng type
    rng = gsl_rng_alloc( T );
    gsl_rng_set( rng, time( 0 ) );                // set the seed equal to simulTime(0)

    // training using the default initial values with random starts
    ExprPar par_default( nFactors(), nSeqs() );
    train( par_default, rng );

    return 0;
}

//Called fairly rarely, don't worry about optimality.
int ExprPredictor::predict( const SiteVec& targetSites_, int targetSeqLength, vector< double >& targetExprs, int seq_num) const
{
	return this->predict(par_model,targetSites_,targetSeqLength,targetExprs,seq_num);
}

/**
In training mode, will skip zero-weighted bins.
*/
int ExprPredictor::predict( const ExprPar& par, const SiteVec& targetSites_, int targetSeqLength, vector< double >& targetExprs, int seq_num ) const
{
	// predict the expression


		//Code for skipping during training BEGIN_SKIPPING
		/*TODO: dynamic_cast is slow, maybe it would be better to move this code that decides
		which bins to predict out to some pre-epoch place so it only gets called once.

		TODO: probably the best thing would be to make the trainingObjective request which sites it wants predicted.
		The reason I don't do that at the moment is that its the setup of the ExprFunc that takes so much time, and we currently cache and reuse that.

		For now, we value correctness above efficiency.
		*/
		Matrix *weights = NULL;
		Weighted_ObjFunc_Mixin* tmp_weighted = NULL;
		tmp_weighted = dynamic_cast<Weighted_ObjFunc_Mixin*>(this->trainingObjective);
		if( NULL != tmp_weighted ) {
			weights = tmp_weighted->get_weights();
		}

		//End of skipping code.	END_SKIPPING


    ExprFunc* func = createExprFunc( par , targetSites_, targetSeqLength, seq_num);
		targetExprs.resize(nConds());
    for ( int j = 0; j < nConds(); j++ )
    {
				//Code for skipping during training BEGIN_SKIPPING
				if( this->is_training() && weights != NULL && weights->getElement(seq_num,j) <= 0.0){
					//cerr << "TEMPORARY DEBUG CODE, skipping unweighted bin (" << seq_num << "," << j << ")." << endl;
					targetExprs[j] = 0.0;
					continue;
				}
				//cerr << "Unskipped bin, making prediction." << endl;
				//End of skipping code. END_SKIPPING


				Condition concs = training_data->getCondition( j , par );
        double predicted = func->predictExpr( concs );
        targetExprs[j] = ( predicted );
    }

    delete func;
    return 0;
}

/**
While in training mode (private in_training variable == true), this method will skip the prediction of zero-weighted positions.
*/
int ExprPredictor::predict_all( const ExprPar& par , vector< vector< double > > &targetExprs ) const
{
	vector< int > seqLengths( seqs.size() );

    for( int i = 0; i < seqs.size(); i++ ){
      seqLengths[i] = seqs[i].size();
    }


    #ifdef REANNOTATE_EACH_PREDICTION
    vector< SiteVec > seqSites( seqs.size() );
    SeqAnnotator ann( expr_model.motifs, par.energyThrFactors );
    for ( int i = 0; i < seqs.size(); i++ ) {
       	ann.annot( seqs[ i ], seqSites[ i ] );
    }
    #else
    const vector< SiteVec >& seqSites = this->seqSites;
    #endif

    //Create predictions for every sequence and condition.
    //Sequences are independent: each iteration only reads shared state (the
    //parameters, the dataset, the model) and writes its own slot of targetExprs,
    //so the result does not depend on the thread schedule.
    const int n = nSeqs();
    targetExprs.assign( n, vector< double >() );
    std::string first_error;
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for ( int i = 0; i < n; i++ ) {
        try {
            vector<double> one_seq_predictions(nConds());
            this->predict(par, seqSites[i], seqLengths[i], one_seq_predictions, i );
            targetExprs[i].swap( one_seq_predictions );
        } catch ( const std::exception& e ) {
            #ifdef _OPENMP
            #pragma omp critical(gemstat_predict_all_error)
            #endif
            if ( first_error.empty() ) first_error = e.what();
        }
    }
    if ( !first_error.empty() ) throw std::runtime_error( first_error );

	return 0;
}




void ExprPredictor::printPar( const ExprPar& par ) const
{
    cout.setf( ios::fixed );
    cout.precision( 8 );
    //     cout.width( 8 );
	cout << par.my_pars;
    cout << flush;
}


ExprFunc* ExprPredictor::createExprFunc( const ExprPar& par, const SiteVec& sites_, const int seq_length, const int seq_num ) const
{

    return expr_model.createNewExprFunc( par, sites_, seq_length, seq_num );
}


int indices_of_crm_in_gene[] =
{
};

double ExprPredictor::evalObjective( const ExprPar& par )
{
	vector<vector<double> > ground_truths;
	vector<vector<double> > predictions;

	for(int i = 0;i< nSeqs();i++){//Populate ground truths
		ground_truths.push_back(training_data->get_output_row(i));
	}

	this->predict_all(par, predictions);

    //Evaluate the objective function on that.
    double ret_val = trainingObjective->eval(ground_truths, predictions, &par);

    return ret_val;

}

int ExprPredictor::simplex_minimize( ExprPar& par_result, double& obj_result )
{
    // 	cout << "Start minimization" << endl;
    // extract initial parameters
    vector < double > pars;

    //ExprPar tmp_par_model = param_factory->changeSpace(par_model, ExprPar::searchOption == CONSTRAINED ? CONSTRAINED_SPACE : ENERGY_SPACE);
    ExprPar tmp_par_model = param_factory->changeSpace(par_model, ENERGY_SPACE);
    param_factory->separateParams(tmp_par_model, free_pars, fix_pars, indicator_bool );

    pars.clear();
    pars = free_pars;

    //SIMPLEX MINIMIZATION with NLOPT
    nlopt::opt optimizer(nlopt::LN_NELDERMEAD, pars.size());
    optimizer.set_min_objective(nlopt_obj_func, this);
    optimizer.set_initial_step(1.0);//TODO: enforce simplex starting size.
	if(max_simplex_iterations > -1){ optimizer.set_maxeval(max_simplex_iterations); }

    if(this->search_option == CONSTRAINED){
      vector<double> free_mins;
      vector<double> fix_mins;

      param_factory->separateParams(param_factory->getMinimums(), free_mins, fix_mins, indicator_bool);
      optimizer.set_lower_bounds(free_mins);

      param_factory->separateParams(param_factory->getMaximums(), free_mins, fix_mins, indicator_bool);
      optimizer.set_upper_bounds(free_mins);
    }

    run_nlopt( optimizer, free_pars, obj_result, "simplex", this );
    //Done Minimizing

    param_factory->joinParams(free_pars, fix_pars, pars, indicator_bool);
    //tmp_par_model = param_factory->create_expr_par(pars, ExprPar::searchOption == CONSTRAINED ? CONSTRAINED_SPACE : ENERGY_SPACE);
    tmp_par_model = param_factory->create_expr_par(pars, ENERGY_SPACE);
    par_result = param_factory->changeSpace(tmp_par_model, PROB_SPACE);

    printPar( par_result );
    cout << endl;

    return 0;
}


int ExprPredictor::gradient_minimize( ExprPar& par_result, double& obj_result )
{
    // 	cout << "Start minimization" << endl;
    // extract initial parameters
    vector< double > pars;
    //cout << "DEBUG: in getFreePars()" << endl;
    //par_model.getFreePars( pars, expr_model.coopMat, expr_model.actIndicators, expr_model.repIndicators );
    //cout << "DEBUG: out getFreePars()" << endl;
    //ExprPar tmp_par_model = param_factory->changeSpace(par_model, ExprPar::searchOption == CONSTRAINED ? CONSTRAINED_SPACE : ENERGY_SPACE);
    ExprPar tmp_par_model = param_factory->changeSpace(par_model, ENERGY_SPACE);

    param_factory->separateParams(tmp_par_model, free_pars, fix_pars, indicator_bool );

    pars.clear();
    pars = free_pars;

    //GRADIENT MINIMIZATION with NLOPT
    nlopt::opt optimizer(nlopt::LD_LBFGS, pars.size());
    optimizer.set_min_objective(nlopt_obj_func, this);

    //TODO: Move this to a nice lookup table or something.
    //Set the stopping criterion
    double ftol;
    switch(objOption){
      case SSE:
        ftol = min_delta_f_SSE;
        break;
      case CORR:
        ftol = min_delta_f_Corr;
        break;
      case CROSS_CORR:
        ftol = min_delta_f_CrossCorr;
        break;
      case PGP:
        ftol = min_delta_f_PGP;
        break;
      default:
        ftol = 1e-5;
        break;
    }
    optimizer.set_ftol_abs(ftol);

    if(this->search_option == CONSTRAINED){
      vector<double> free_mins;
      vector<double> fix_mins;

      param_factory->separateParams(param_factory->getMinimums(), free_mins, fix_mins, indicator_bool);
      optimizer.set_lower_bounds(free_mins);

      param_factory->separateParams(param_factory->getMaximums(), free_mins, fix_mins, indicator_bool);
      optimizer.set_upper_bounds(free_mins);
    }


    //TODO: enforce nGradientIters
	if(max_gradient_iterations > -1){ optimizer.set_maxeval(max_gradient_iterations); }

    run_nlopt( optimizer, free_pars, obj_result, "gradient descent", this );

    //Done Minimizing
    //pars now contains the optimal parameters

    param_factory->joinParams(free_pars, fix_pars, pars, indicator_bool);
    //tmp_par_model = param_factory->create_expr_par(pars, ExprPar::searchOption == CONSTRAINED ? CONSTRAINED_SPACE : ENERGY_SPACE);
    tmp_par_model = param_factory->create_expr_par(pars, ENERGY_SPACE);

    par_result = param_factory->changeSpace(tmp_par_model, PROB_SPACE);
    return 0;
}


double nlopt_obj_func( const vector<double> &x, vector<double> &grad, void* f_data){
		ExprPredictor* predictor = (ExprPredictor*)f_data;
		predictor->begin_batch();

        gsl_vector *xv = vector2gsl(x); //TODO: Ugly, remove (Make all objective functions use native STL vectors)
        double objective = gsl_obj_f(xv, f_data);

        if(!grad.empty()){
                gsl_vector *dxv = vector2gsl(grad);

                gsl_obj_df(xv,f_data,dxv,objective); // objective is f(x): don't evaluate it a second time

                for(int i = 0;i< grad.size();i++){
                        grad[i] = dxv->data[i];
                }
                cerr << " obj called, derivative : " << grad << endl;
                gsl_vector_free(dxv);
        }

        cerr << " obj called " << objective << endl;

        gsl_vector_free(xv);
        return objective;
}

double gsl_obj_f( const gsl_vector* v, void* params )
{
    // the ExprPredictor object
    ExprPredictor* predictor = (ExprPredictor*)params;

    // parse the variables (parameters to be optimized)
    //     vector< double > expv;
    //     for ( int i = 0; i < v->size; i++ ) expv.push_back( exp( gsl_vector_get( v, i ) ) );
    vector <double> temp_free_pars = gsl2vector(v);
    vector < double > all_pars;

    predictor->param_factory->joinParams(temp_free_pars, predictor->fix_pars, all_pars, predictor->indicator_bool);
    //ExprPar par = predictor->param_factory->create_expr_par(all_pars, ExprPar::searchOption == CONSTRAINED ? CONSTRAINED_SPACE : ENERGY_SPACE);
    ExprPar par = predictor->param_factory->create_expr_par(all_pars, ENERGY_SPACE);
    par = predictor->param_factory->changeSpace(par, PROB_SPACE); //TODO: WTF? This shouldn't be required because it's done in the createExprFunc method. Stack corruption or something?


    // call the ExprPredictor object to evaluate the objective function
    double obj = predictor->objFunc( par );
    return obj;
}

void gsl_obj_df( const gsl_vector* v, void* params, gsl_vector* grad )
{
    gsl_obj_df( v, params, grad, gsl_obj_f( v, params ) );
}

void gsl_obj_df( const gsl_vector* v, void* params, gsl_vector* grad, double f_val )
{
    ExprPredictor* predictor = (ExprPredictor*)params;
    if ( predictor->gradient_method == ExprPredictor::GRADIENT_AD ) gsl_obj_df_ad( v, params, grad );
    else gsl_obj_df_fd( v, params, grad, f_val );
}

void gsl_obj_df_fd( const gsl_vector* v, void* params, gsl_vector* grad, double f_val )
{
    double step = 1.0E-6;
    numeric_deriv( grad, gsl_obj_f, v, params, step, f_val );
}

void gsl_obj_df_ad( const gsl_vector* v, void* params, gsl_vector* grad )
{
    ExprPredictor* predictor = (ExprPredictor*)params;

    // the same parameter object gsl_obj_f evaluates
    vector< double > temp_free_pars = gsl2vector( v );
    vector< double > all_pars;
    predictor->param_factory->joinParams( temp_free_pars, predictor->fix_pars, all_pars, predictor->indicator_bool );
    ExprPar par = predictor->param_factory->create_expr_par( all_pars, ENERGY_SPACE );
    par = predictor->param_factory->changeSpace( par, PROB_SPACE );

    vector< double > g_prob;
    predictor->gradient_prob( par, g_prob );

    // the optimizer works in ENERGY_SPACE: p_k = exp( e_k ), so d/de_k = p_k d/dp_k;
    // and only on the free parameters
    vector< double > flat_prob;
    par.getRawPars( flat_prob );
    int c = 0;
    for ( size_t k = 0; k < flat_prob.size(); k++ )
    {
        if ( predictor->indicator_bool[k] ) gsl_vector_set( grad, c++, g_prob[k] * flat_prob[k] );
    }
}

void ExprPredictor::gradient_prob( const ExprPar& par, vector< double >& grad ) const
{
    if ( par.my_space != PROB_SPACE ) throw std::invalid_argument( "ExprPredictor::gradient_prob: par must be in PROB_SPACE" );
    typedef gemstat_ad_t V;
    typedef gemstat_ad::Tape< gemstat_dp_t > TapeT;

    vector< double > flat;
    par.getRawPars( flat );
    const int n_pars = flat.size();
    const int n = nSeqs();
    const int nc = nConds();

    // the objective's own derivatives, at the plain predictions
    vector< vector< double > > ground_truths( n ), predictions;
    for ( int i = 0; i < n; i++ ) ground_truths[i] = training_data->get_output_row( i );
    this->predict_all( par, predictions );
    vector< vector< double > > d_pred;
    vector< double > d_pars( n_pars, 0.0 );
    trainingObjective->gradient( ground_truths, predictions, &par, &par_index, d_pred, d_pars );

    // zero-weighted bins are not predicted during training (see predict())
    Matrix* weights = NULL;
    Weighted_ObjFunc_Mixin* tmp_weighted = dynamic_cast< Weighted_ObjFunc_Mixin* >( this->trainingObjective );
    if ( NULL != tmp_weighted ) weights = tmp_weighted->get_weights();

    vector< int > seqLengths( n );
    for ( int i = 0; i < n; i++ ) seqLengths[i] = seqs[i].size();

    // one reverse pass per sequence, seeded with d objective / d prediction
    vector< vector< double > > seq_grad( n, vector< double >( n_pars, 0.0 ) );
    std::string first_error;
    #ifdef _OPENMP
    #pragma omp parallel
    #endif
    {
        TapeT tape;
        tape.activate();
        #ifdef _OPENMP
        #pragma omp for schedule(dynamic)
        #endif
        for ( int i = 0; i < n; i++ )
        {
            try {
                tape.clear();
                vector< V > vars( n_pars );
                for ( int k = 0; k < n_pars; k++ ) vars[k] = V::input( flat[k] );

                ExprFunc* func = createExprFunc( par, seqSites[i], seqLengths[i], i );
                ThermoVals< V > vals = func->makeVals( vars, par_index );
                vector< int > out_index( nc, -1 );
                for ( int j = 0; j < nc; j++ )
                {
                    if ( this->is_training() && weights != NULL && weights->getElement( i, j ) <= 0.0 ) continue;
                    Condition concs = training_data->getCondition( j, par );
                    V p = func->predictExprAD( vals, concs.concs );
                    out_index[j] = p.index();
                }
                delete func;

                vector< gemstat_dp_t > adjoint( tape.size(), (gemstat_dp_t)0 );
                for ( int j = 0; j < nc; j++ )
                {
                    if ( out_index[j] >= 0 ) adjoint[ out_index[j] ] += d_pred[i][j];
                }
                tape.backward( adjoint );
                for ( int k = 0; k < n_pars; k++ ) seq_grad[i][k] = (double)adjoint[ vars[k].index() ];
            } catch ( const std::exception& e ) {
                #ifdef _OPENMP
                #pragma omp critical(gemstat_gradient_error)
                #endif
                if ( first_error.empty() ) first_error = e.what();
            }
        }
    }
    if ( !first_error.empty() ) throw std::runtime_error( first_error );

    // sum in sequence order so the result does not depend on the thread count
    grad = d_pars;
    for ( int i = 0; i < n; i++ )
        for ( int k = 0; k < n_pars; k++ ) grad[k] += seq_grad[i][k];
}

bool ExprPredictor::checkGradient( const ExprPar& par_init, ostream& os, double tol )
{
    par_model = par_init;
    ExprPar tmp_par_model = param_factory->changeSpace( par_model, ENERGY_SPACE );
    param_factory->separateParams( tmp_par_model, free_pars, fix_pars, indicator_bool );
    int n = free_pars.size();
    gsl_vector* v = vector2gsl( free_pars );
    gsl_vector* g_ad = gsl_vector_alloc( n );
    gsl_vector* g_fd = gsl_vector_alloc( n );

    gsl_obj_df_ad( v, this, g_ad );

    // central differences in ENERGY_SPACE
    gsl_vector* dv = gsl_vector_alloc( n );
    for ( int k = 0; k < n; k++ )
    {
        double x = gsl_vector_get( v, k );
        double h = 1.0e-5 * ( fabs( x ) > 1.0 ? fabs( x ) : 1.0 );
        gsl_vector_memcpy( dv, v );
        gsl_vector_set( dv, k, x + h ); double fp = gsl_obj_f( dv, this );
        gsl_vector_set( dv, k, x - h ); double fm = gsl_obj_f( dv, this );
        gsl_vector_set( g_fd, k, ( fp - fm ) / ( 2.0 * h ) );
    }

    // names of the free parameters, for the report
    vector< std::string > free_names;
    {
        size_t k = 0;
        for ( gsparams::DictList::iterator itr = par_model.my_pars.begin(); itr != par_model.my_pars.end(); ++itr, ++k )
            if ( indicator_bool[k] ) free_names.push_back( itr.get_path() );
    }

    bool ok = true;
    double worst = 0.0;
    os << "GRADIENT CHECK: " << n << " free parameters, objective " << getObjOptionStr( objOption ) << endl;
    os << "parameter	autodiff	central_difference	relative_difference" << endl;
    for ( int k = 0; k < n; k++ )
    {
        double a = gsl_vector_get( g_ad, k ), f = gsl_vector_get( g_fd, k );
        double scale = fabs( a ) > fabs( f ) ? fabs( a ) : fabs( f );
        double rel = scale > tol ? fabs( a - f ) / scale : 0.0;   // both tiny: agree
        if ( rel > worst ) worst = rel;
        bool this_ok = rel <= tol;
        if ( !this_ok ) ok = false;
        os << free_names[k] << "	" << setprecision(10) << a << "	" << f << "	" << setprecision(3) << rel << ( this_ok ? "" : "	MISMATCH" ) << endl;
    }
    os << "GRADIENT CHECK " << ( ok ? "PASSED" : "FAILED" ) << endl;
    os << "GRADIENT CHECK worst relative difference " << setprecision(3) << worst << " (tolerance " << tol << ")" << endl;

    gsl_vector_free( v ); gsl_vector_free( dv ); gsl_vector_free( g_ad ); gsl_vector_free( g_fd );
    return ok;
}


void gsl_obj_fdf( const gsl_vector* v, void* params, double* result, gsl_vector* grad )
{
    *result = gsl_obj_f( v, params );
    gsl_obj_df( v, params, grad );
}
