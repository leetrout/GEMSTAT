#ifndef GEMSTAT_PARAM_SLOTS_H
#define GEMSTAT_PARAM_SLOTS_H

#include "ExprPar.h"

/*
 * Where each ExprFunc-level parameter lives in the flat (traversal order)
 * PROB_SPACE parameter vector of an ExprPar.  Built once by
 * ParamSlots::build() from the parameter dictionary; -1 means "not present"
 * (e.g. no interaction parameter for a pair of factors: the interaction is
 * then the constant 1).
 */
struct ParamSlots
{
    vector< int > maxbind;              // per factor: tfs[f].maxbind
    vector< int > alpha_a;              // per factor: tfs[f].alpha_a
    vector< int > alpha_r;              // per factor: tfs[f].alpha_r
    vector< vector< int > > inter;      // per factor pair (symmetric): inter[a:b]
    vector< int > qbtm;                 // per basal-transcription slot: qbtm[i]
    vector< int > beta;                 // per enhancer slot: enh[i].beta
    vector< int > pi;                   // per enhancer slot: enh[i].pi
    int n_pars;                         // total length of the flat vector

    static ParamSlots build( const ExprPar& par, const vector< string >& motifNames );
    static void flatten_paths( const gsparams::DictList& d, const std::string& prefix, vector< std::string >& out );
};

#endif
