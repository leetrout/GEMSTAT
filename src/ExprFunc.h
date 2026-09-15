#ifndef EXPR_FUNC_H
#define EXPR_FUNC_H

#include "ExprModel.h"
#include "FactorIntFunc.h"
#include "ExprPar.h"
#include "DataSet.h"
#include "tools/ReverseAD.h"
#include "ParamSlots.h"

/*****************************************************
 * Expression Model and Parameters
 ******************************************************/

typedef long double gemstat_dp_t;
typedef gemstat_ad::Var< gemstat_dp_t > gemstat_ad_t;   // the same computations, recorded for reverse-mode differentiation


/*
 * The values the recurrences work on, for one sequence, in scalar type T.
 * T is gemstat_dp_t for an ordinary prediction and gemstat_ad_t to record
 * the computation for differentiation.
 */
template< class T >
struct ThermoVals
{
    vector< T > maxBindingWts;          // per factor: K(S_max)[TF_max]
    vector< T > txpEffects;             // per factor: alpha (activation)
    vector< T > repEffects;             // per factor: beta (repression) under ChrMod, quenching efficiency etc.
    vector< T > bindingWts;             // per site incl. pseudo-sites; depends on the condition
    vector< vector< T > > w_ij;         // per ExprFunc::left_nbrs[i][k]: interaction( sites[i], sites[j] )
    vector< vector< T > > w_ji;         // per ExprFunc::left_nbrs[i][k]: interaction( sites[j], sites[i] )
    vector< vector< T > > all_w_ji;     // per ExprFunc::all_left_nbrs[i][k] (Quenching on-state)
    vector< vector< T > > right_w_ij;   // per ExprFunc::right_nbrs[i][k] (Markov backward pass)
    T basal;                            // q_btm for this sequence
};

/* ExprFunc class: predict the expression (promoter occupancy) of an enhancer sequence */
class ExprFunc
{
    public:
        // constructors
        ExprFunc( const ExprModel* _model, const ExprPar& _par , const SiteVec& sites_, const int seq_len, const int seq_num);
        virtual ~ExprFunc() {}  // ExprFuncs are created by ExprModel and deleted through this base pointer

        // access methods
        const vector< Motif >& getMotifs() const
        {
            return motifs;
        }

        // predict the expression value of a given sequence (its site representation, sorted by the start positions) under given TF concentrations
        virtual double predictExpr( const vector< double >& factorConcs );
        virtual double predictExpr( const Condition& in_condition );
        const ExprPar& getPar() const { return par; }

        /*
         * Differentiable prediction.  flat_pars are the PROB_SPACE parameters of
         * this sequence's ExprPar as gemstat_ad_t inputs of the active tape (see
         * ParamSlots); the returned value is recorded on that tape.  makeVals()
         * derives the per-sequence quantities once; predictExprAD() is then
         * called per condition.
         */
        ThermoVals< gemstat_ad_t > makeVals( const vector< gemstat_ad_t >& flat_pars, const ParamSlots& slots ) const;
        virtual gemstat_ad_t predictExprAD( ThermoVals< gemstat_ad_t >& vals, const vector< double >& factorConcs ) const;

        //static ModelType modelOption;             // model option
    protected:
        //setup functions that may be useful to subclasses
        virtual void setupSitesAndBoundaries(const SiteVec& _sites, int length, int seq_num);
        void setupBindingWeights(const vector< double >& factorConcs);
        // TF binding motifs

        const ExprModel* expr_model;
        const vector< Motif >& motifs;

        // control parameters
//        const FactorIntFunc* intFunc;             // function to compute distance-dependent TF-TF interactions
        const vector< bool >& actIndicators;      // 1 if the TF is in the activator set
        int maxContact;                           // the maximum contact
        const vector< bool >& repIndicators;      // 1 if the TF is in the repressor set
        const IntMatrix& repressionMat;           // repression matrix: R(f,f') = 1 if f can repress f'
        double repressionDistThr;                 // distance threshold for repression: d_R
        int seq_number;
        int seq_length;

        // model parameters
        ExprPar par;//NOTE: Removing "const" here caused the copy constructor to be called. Thus the ExprFunc gets its own copy that will not have problems when the original par is changed.
                          //NOTE: (Additional) put const back, copying is slow, and par should be constant for an ExprFunc, this also fixed a memory leak.

        // the sequence whose expression is to be predicted
        SiteVec sites;
        int n_sites;  //Useful because the sitevec with pseudosites etc. might change, this should be the number of true sites.
        vector< int > boundaries;                 // left boundary of each site beyond which there is no interaction

        /*
         * Per-sequence cache of everything in the DP inner loops that does not
         * depend on the condition (the TF concentrations): which earlier sites a
         * site can interact with (no overlap, inside the interaction window), the
         * pairwise interaction weight in both argument orders, the distance and
         * both repression flags.  The DP recurrences below are evaluated once per
         * condition, so without this cache all of that was recomputed for every
         * condition even though it only changes with the parameters.  Lists keep
         * the iteration order of the original loops so results are bit-identical.
         */
        struct SiteInteraction {
            int j;          // index of the other site
            int dist;       // |sites[i].start - sites[j].start|
            double w_ij;    // compFactorInt( sites[i], sites[j] )
            double w_ji;    // compFactorInt( sites[j], sites[i] )
            bool rep_ij;    // testRepression( sites[i], sites[j] )
            bool rep_ji;    // testRepression( sites[j], sites[i] )
        };
        vector< vector< SiteInteraction > > left_nbrs;      // left_nbrs[i]: sites j in (boundaries[i], i), increasing j
        vector< vector< SiteInteraction > > all_left_nbrs;  // every non-overlapping site j < i (Quenching on-state; empty otherwise)
        vector< vector< SiteInteraction > > right_nbrs;     // sites j > i in decreasing j (Markov backward pass; empty otherwise)
        vector< int > rev_bounds;                           // right boundary of each site (Markov backward pass; empty otherwise)
        // Fill out[i] with the non-overlapping sites j < i, in increasing j.
        // within_boundaries limits j to (boundaries[i], i) as most recurrences do;
        // the Quenching on-state recurrence considers every earlier site.
        void buildLeftNeighbours( bool within_boundaries, vector< vector< SiteInteraction > >& out ) const;

        // the values the recurrences use for an ordinary (non-differentiated) prediction
        ThermoVals< gemstat_dp_t > plain;
        void fillPlainWeights();   // copy the cached interaction weights of every neighbour list into plain

        // compute the TF-TF interaction between two occupied sites
        double compFactorInt( const Site& a, const Site& b ) const;
        // the same interaction as an affine function of the pair's interaction parameter, in scalar type T
        template< class T > T factorIntAffine( const Site& a, const Site& b, const T& normalInt ) const;

        double compDen() const;
        double compNum( int TFid ) const;

        // test if one site represses another site
        bool testRepression( const Site& a, const Site& b ) const;

        // compute the partition function when the BTM is bound
        virtual gemstat_dp_t compPartFuncOn() const;
        virtual gemstat_ad_t compPartFuncOnAD( const ThermoVals< gemstat_ad_t >& v ) const;
        // compute the partition function when the basal transcriptional machinery (BTM) is not bound
        virtual gemstat_dp_t compPartFuncOff() const;
        virtual gemstat_ad_t compPartFuncOffAD( const ThermoVals< gemstat_ad_t >& v ) const;

        /*
         * The recurrences themselves, written once for any scalar type.  Each
         * model's compPartFunc*() / compPartFunc*AD() pair just calls its kernel.
         */
        template< class T > T kernelOffBasic( const ThermoVals< T >& v ) const;
        template< class T > T kernelOffChrMod( const ThermoVals< T >& v ) const;
        template< class T > T kernelOnDirect( const ThermoVals< T >& v ) const;
        template< class T > T kernelOnQuenching( const ThermoVals< T >& v ) const;
        template< class T > T kernelOnChrModUnlimited( const ThermoVals< T >& v ) const;
        template< class T > T kernelOnChrModLimited( const ThermoVals< T >& v ) const;
        template< class T > T kernelLogistic( const ThermoVals< T >& v ) const;
        template< class T > T kernelMarkov( const ThermoVals< T >& v ) const;
        template< class T > static T thermoOccupancy( const T& Z_off, const T& Z_on, const T& basal );
        template< class T > void setupBindingWeightsT( ThermoVals< T >& v, const vector< double >& factorConcs ) const;

        /*
        At setup time, this will get populated from the SNOT object.
        */
        // parameters
        vector < GEMSTAT_PAR_FLOAT_T > maxBindingWts;          // binding weight of the strongest site for each TF: K(S_max) [TF_max]
        Matrix factorIntMat;                      // (maximum) interactions between pairs of factors: omega(f,f')
        vector < GEMSTAT_PAR_FLOAT_T > txpEffects;             // transcriptional effects: alpha for Direct and Quenching model, exp(alpha) for Logistic model (so that the same default values can be used). Equal to 1 if a TF is not an activator under the Quenching model
        vector < GEMSTAT_PAR_FLOAT_T > repEffects;             // repression effects: beta under ChrMod models (the equlibrium constant of nucleosome association with chromatin). Equal to 0 if a TF is not a repressor.
        //vector < GEMSTAT_PAR_FLOAT_T > pis;
        //     double expRatio; 		// constant factor of measurement to prediction

        vector < GEMSTAT_PAR_FLOAT_T > betas;
        //vector < GEMSTAT_PAR_FLOAT_T > energyThrFactors;

    private:
};

class Logistic_ExprFunc : public ExprFunc {
  public:
      // constructors
      Logistic_ExprFunc( const ExprModel* _model, const ExprPar& _par , const SiteVec& sites_, const int seq_len, const int seq_num) : ExprFunc( _model, _par , sites_, seq_len, seq_num){} ;

      double predictExpr( const vector< double >& factorConcs );
      gemstat_ad_t predictExprAD( ThermoVals< gemstat_ad_t >& vals, const vector< double >& factorConcs ) const;
};

class Markov_ExprFunc : public ExprFunc {
  public:
      Markov_ExprFunc( const ExprModel* _model, const ExprPar& _par , const SiteVec& sites_, const int seq_len, const int seq_num);// : ExprFunc( _model, _par , sites_, seq_len, seq_num);
      double predictExpr( const vector< double >& factorConcs );
      gemstat_ad_t predictExprAD( ThermoVals< gemstat_ad_t >& vals, const vector< double >& factorConcs ) const;
};

class Direct_ExprFunc : public ExprFunc {
  public:
      // constructors
      Direct_ExprFunc( const ExprModel* _model, const ExprPar& _par , const SiteVec& sites_, const int seq_len, const int seq_num) : ExprFunc( _model, _par , sites_, seq_len, seq_num){} ;
  protected:
    // compute the partition function when the BTM is bound
    gemstat_dp_t compPartFuncOn() const;
    gemstat_ad_t compPartFuncOnAD( const ThermoVals< gemstat_ad_t >& v ) const;

};

class Quenching_ExprFunc : public ExprFunc {
  public:
      // constructors
      Quenching_ExprFunc( const ExprModel* _model, const ExprPar& _par , const SiteVec& sites_, const int seq_len, const int seq_num);
  protected:
    // compute the partition function when the BTM is bound
    gemstat_dp_t compPartFuncOn() const;
    gemstat_ad_t compPartFuncOnAD( const ThermoVals< gemstat_ad_t >& v ) const;
};

class ChrMod_ExprFunc : public ExprFunc {
  public:
      // constructors
      ChrMod_ExprFunc( const ExprModel* _model, const ExprPar& _par , const SiteVec& sites_, const int seq_len, const int seq_num) : ExprFunc( _model, _par , sites_, seq_len, seq_num){} ;
  protected:
    // compute the partition function when the BTM is bound
    virtual gemstat_dp_t compPartFuncOn() const = 0;
    virtual gemstat_ad_t compPartFuncOnAD( const ThermoVals< gemstat_ad_t >& v ) const = 0;
    // compute the partition function when the basal transcriptional machinery (BTM) is not bound
    gemstat_dp_t compPartFuncOff() const;
    gemstat_ad_t compPartFuncOffAD( const ThermoVals< gemstat_ad_t >& v ) const;

};

class ChrModUnlimited_ExprFunc : public ChrMod_ExprFunc {
  public:
      // constructors
      ChrModUnlimited_ExprFunc( const ExprModel* _model, const ExprPar& _par , const SiteVec& sites_, const int seq_len, const int seq_num) : ChrMod_ExprFunc( _model, _par , sites_, seq_len, seq_num){} ;
  protected:
    // compute the partition function when the BTM is bound
    gemstat_dp_t compPartFuncOn() const;
    gemstat_ad_t compPartFuncOnAD( const ThermoVals< gemstat_ad_t >& v ) const;
};

class ChrModLimited_ExprFunc : public ChrMod_ExprFunc {
  public:
      // constructors
      ChrModLimited_ExprFunc( const ExprModel* _model, const ExprPar& _par , const SiteVec& sites_, const int seq_len, const int seq_num) : ChrMod_ExprFunc( _model, _par , sites_, seq_len, seq_num){} ;
  protected:
    // compute the partition function when the BTM is bound
    gemstat_dp_t compPartFuncOn() const;
    gemstat_ad_t compPartFuncOnAD( const ThermoVals< gemstat_ad_t >& v ) const;
};

#include "ExprFuncKernels.hpp"

#endif
