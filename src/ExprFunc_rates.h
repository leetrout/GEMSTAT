#ifndef EXPR_FUNC_RATES_H
#define EXPR_FUNC_RATES_H

#include "ExprFunc.h"

/*
 * Rates model.  Four partition functions are computed at once, sharing the
 * loop over the sites (about a third faster than four separate recurrences):
 *   O  : sites bound, no arc effect
 *   A  : every bound site carries its alpha_a
 *   B  : every bound site carries its alpha_r
 *   AB : every bound site carries both
 * q_btm scales arc A, the "pi" promoter parameter scales arc B, and the
 * expression is the harmonic combination 2 pA pB / (pA + pB) of the two arc
 * probabilities.
 */
class Rates_ExprFunc : public ExprFunc {
  public:
      Rates_ExprFunc( const ExprModel* _model, const ExprPar& _par , const SiteVec& sites_, const int seq_len, const int seq_num) : ExprFunc( _model, _par , sites_, seq_len, seq_num){} ;
      double predictExpr( const vector< double >& factorConcs );
      gemstat_ad_t predictExprAD( ThermoVals< gemstat_ad_t >& vals, const vector< double >& factorConcs ) const;
  protected:
      template< class T > T kernelRates( const ThermoVals< T >& v ) const;
};

template< class T >
T Rates_ExprFunc::kernelRates( const ThermoVals< T >& v ) const
{
    int n = n_sites;
    struct Parts { T O, A, B, AB; };
    vector< Parts > Z( n + 1 ), Zt( n + 1 );
    Z[0].O = 1.0; Z[0].A = 1.0; Z[0].B = 1.0; Z[0].AB = 1.0;
    Zt[0] = Z[0];

    for ( int i = 1; i <= n; i++ )
    {
        const T& alpha_a = v.txpEffects[ sites[ i ].factorIdx ];
        const T& alpha_r = v.repEffects[ sites[ i ].factorIdx ];
        const vector< SiteInteraction >& nbrs = left_nbrs[i];
        const vector< T >& w = v.w[i];

        const Parts& start = Zt[ boundaries[i] ];
        gemstat_ad::Accumulator< T > accO( start.O ), accA( start.A ), accB( start.B ), accAB( start.AB );
        for ( size_t k = 0; k < nbrs.size(); k++ )
        {
            const Parts& zj = Z[ nbrs[k].j ];
            accO.add( w[k], zj.O );
            accA.add( w[k], zj.A );
            accB.add( w[k], zj.B );
            accAB.add( w[k], zj.AB );
        }
        Z[i].O  = v.bindingWts[i] * accO.result();
        Z[i].A  = v.bindingWts[i] * accA.result()  * alpha_a;
        Z[i].B  = v.bindingWts[i] * accB.result()  * alpha_r;
        Z[i].AB = v.bindingWts[i] * accAB.result() * alpha_a * alpha_r;

        Zt[i].O  = Z[i].O  + Zt[i - 1].O;
        Zt[i].A  = Z[i].A  + Zt[i - 1].A;
        Zt[i].B  = Z[i].B  + Zt[i - 1].B;
        Zt[i].AB = Z[i].AB + Zt[i - 1].AB;
    }

    T Z_O  = Zt[n].O;
    T Z_A  = Zt[n].A  * v.basal;
    T Z_B  = Zt[n].B  * v.pi;
    T Z_AB = Zt[n].AB * ( v.basal * v.pi );
    T Z_total = Z_O + Z_A + Z_B + Z_AB;
    T prob_A = ( Z_A + Z_AB ) / Z_total;
    T prob_B = ( Z_B + Z_AB ) / Z_total;
    return 2.0 * ( prob_A * prob_B ) / ( prob_A + prob_B );
}

#endif //EXPR_FUNC_RATES_H
