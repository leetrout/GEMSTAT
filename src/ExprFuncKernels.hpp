/*
 * The partition-function recurrences of every GEMSTAT model, written once for
 * any scalar type T.  Included at the end of ExprFunc.h.
 *
 * T is gemstat_dp_t (long double) for an ordinary prediction and gemstat_ad_t
 * (a recorded long double, see tools/ReverseAD.h) when the gradient is wanted.
 * The loops, the order of the floating point operations and the NaN checks
 * are those of the original per-model functions.
 */
#ifndef EXPR_FUNC_KERNELS_HPP
#define EXPR_FUNC_KERNELS_HPP

#include <iostream>
#include <cmath>
#include <cstdlib>

template< class T >
inline T gemstat_logistic( const T& x )
{
    return 1.0 / ( 1.0 + exp( -x ) );
}

template< class T >
void ExprFunc::setupBindingWeightsT( ThermoVals< T >& v, const vector< double >& factorConcs ) const
{
    v.bindingWts.resize( sites.size() );
    v.bindingWts[0] = 1.0;                       // first pseudo-site
    v.bindingWts[ v.bindingWts.size() - 1 ] = 1.0;
    for ( int i = 1; i <= n_sites; i++ )
    {
        const Site& s = sites[i];
        v.bindingWts[i] = v.maxBindingWts[ s.factorIdx ] * factorConcs[ s.factorIdx ] * s.prior_probability * s.wtRatio;
    }
}

template< class T >
T ExprFunc::thermoOccupancy( const T& Z_off, const T& Z_on, const T& basal )
{
    T efficiency = Z_on / Z_off;
    return efficiency * basal / ( 1.0 + efficiency * basal );
}

template< class T >
T ExprFunc::kernelOffBasic( const ThermoVals< T >& v ) const
{
    int n = n_sites;
    vector< T > Z( n + 1 );
    Z[0] = 1.0;
    vector< T > Zt( n + 1 );
    Zt[0] = 1.0;

    for ( int i = 1; i <= n; i++ )
    {
        T sum = Zt[ boundaries[i] ];
        if ( sum != sum )
        {
            cout << "DEBUG: sum nan" << "\t" << gemstat_ad::value_of( Zt[ boundaries[i] ] ) << endl;
            exit(1);
        }
        const vector< SiteInteraction >& nbrs = left_nbrs[i];
        const vector< T >& w = v.w_ij[i];
        gemstat_ad::Accumulator< T > acc( sum );
        for ( size_t k = 0; k < nbrs.size(); k++ )
        {
            acc.add( w[k], Z[ nbrs[k].j ] );
        }
        sum = acc.result();
        if ( sum != sum || isinf( sum ) )
        {
            cout << "DEBUG: sum nan/inf at site " << i << " (factor " << sites[ i ].factorIdx << ")\t" << gemstat_ad::value_of( sum ) << endl;
            for ( size_t k = 0; k < nbrs.size(); k++ )
                cout << "  j=" << nbrs[k].j << "\tfactor " << sites[ nbrs[k].j ].factorIdx << "\tw=" << gemstat_ad::value_of( w[k] ) << "\tZ[j]=" << gemstat_ad::value_of( Z[ nbrs[k].j ] ) << endl;
            exit(1);
        }

        Z[i] = v.bindingWts[ i ] * sum;
        if ( Z[i] != Z[i] )
        {
            cout << "DEBUG: Z bindingWts[i]: " << sites[i].factorIdx << "\t" << gemstat_ad::value_of( v.bindingWts[ i ] ) << "\t" << gemstat_ad::value_of( sum ) << endl;
            exit(1);
        }
        Zt[i] = Z[i] + Zt[i - 1];
    }
    return Zt[n];
}

template< class T >
T ExprFunc::kernelOffChrMod( const ThermoVals< T >& v ) const
{
    int n = n_sites;
    vector< T > Z0( n + 1 );
    Z0[0] = 1.0;
    vector< T > Z1( n + 1 );
    Z1[0] = 1.0;
    vector< T > Zt( n + 1 );
    Zt[0] = 1.0;

    for ( int i = 1; i <= n; i++ )
    {
        T sum = Zt[ boundaries[i] ];
        const vector< SiteInteraction >& nbrs = left_nbrs[i];
        const vector< T >& w = v.w_ji[i];
        const bool rep_i = repIndicators[ sites[i].factorIdx ];
        gemstat_ad::Accumulator< T > acc0( sum );
        for ( size_t k = 0; k < nbrs.size(); k++ )
        {
            int j = nbrs[k].j;
            acc0.add( w[k], Z0[j] );
            if ( nbrs[k].dist > repressionDistThr ) acc0.add( Z1[j] );
        }
        T sum0 = acc0.result();
        T sum1 = sum;
        if ( rep_i )
        {
            gemstat_ad::Accumulator< T > acc1( sum );
            for ( size_t k = 0; k < nbrs.size(); k++ )
            {
                int j = nbrs[k].j;
                acc1.add( w[k], Z1[j] );
                if ( nbrs[k].dist > repressionDistThr ) acc1.add( Z0[j] );
            }
            sum1 = acc1.result();
        }
        Z0[i] = v.bindingWts[i] * sum0;
        if ( rep_i ) Z1[i] = v.bindingWts[i] * v.repEffects[ sites[i].factorIdx ] * sum1;
        else Z1[i] = 0.0;
        Zt[i] = Z0[i] + Z1[i] + Zt[i - 1];
    }
    return Zt[n];
}

template< class T >
T ExprFunc::kernelOnDirect( const ThermoVals< T >& v ) const
{
    int n = n_sites;
    vector< T > Z( n + 1 );
    Z[0] = 1.0;
    vector< T > Zt( n + 1 );
    Zt[0] = 1.0;

    for ( int i = 1; i <= n; i++ )
    {
        const vector< SiteInteraction >& nbrs = left_nbrs[i];
        const vector< T >& w = v.w_ji[i];
        gemstat_ad::Accumulator< T > acc( Zt[ boundaries[i] ] );
        for ( size_t k = 0; k < nbrs.size(); k++ )
        {
            acc.add( w[k], Z[ nbrs[k].j ] );
        }
        T sum = acc.result();
        if ( actIndicators[ sites[ i ].factorIdx ] )
        {
            Z[ i ] = v.bindingWts[ i ] * v.txpEffects[ sites[ i ].factorIdx ] * sum;
        }
        if ( repIndicators[ sites[ i ].factorIdx ] )
        {
            Z[ i ] = v.bindingWts[ i ] * v.repEffects[ sites[ i ].factorIdx ] * sum;
        }
        Zt[i] = Z[i] + Zt[i - 1];
    }
    return Zt[n];
}

template< class T >
T ExprFunc::kernelOnQuenching( const ThermoVals< T >& v ) const
{
    int n = n_sites;
    int N0 = maxContact;
    const int W = N0 + 1;
    vector< T > Z1( ( n + 1 ) * W );
    vector< T > Z0( ( n + 1 ) * W );
    #define GS_Q(M, i, k) M[ (i) * W + (k) ]

    // k = 0
    for ( int i = 0; i <= n; i++ )
    {
        const vector< SiteInteraction >& nbrs = all_left_nbrs[i];
        const vector< T >& w = v.all_w_ji[i];
        gemstat_ad::Accumulator< T > acc1( T( 1.0 ) ), acc0( T( 0.0 ) );
        for ( size_t kk = 0; kk < nbrs.size(); kk++ )
        {
            int j = nbrs[kk].j;
            T z = GS_Q(Z1, j, 0) + GS_Q(Z0, j, 0);
            if ( nbrs[kk].rep_ji ) acc0.add( w[kk], z );
            else acc1.add( w[kk], z );
        }
        GS_Q(Z1, i, 0) = v.bindingWts[i] * acc1.result();
        GS_Q(Z0, i, 0) = v.bindingWts[i] * acc0.result();
    }

    // k >= 1
    for ( int k = 1; k <= N0; k++ )
    {
        for ( int i = 0; i <= n; i++ )
        {
            if ( i < k )
            {
                GS_Q(Z1, i, k) = 0.0;
                GS_Q(Z0, i, k) = 0.0;
                continue;
            }
            const vector< SiteInteraction >& nbrs = all_left_nbrs[i];
            const vector< T >& w = v.all_w_ji[i];
            gemstat_ad::Accumulator< T > acc1( T( 0.0 ) ), acc0( T( 0.0 ) );
            for ( size_t kk = 0; kk < nbrs.size(); kk++ )
            {
                int j = nbrs[kk].j;
                int act_not_rep = ( actIndicators[ sites[j].factorIdx ] ? 1 : 0 ) * ( 1 - ( nbrs[kk].rep_ij ? 1 : 0 ) );
                T effect = act_not_rep * GS_Q(Z1, j, k - 1) * v.txpEffects[ sites[j].factorIdx ];
                T z = GS_Q(Z1, j, k) + GS_Q(Z0, j, k) + effect;
                if ( nbrs[kk].rep_ji ) acc0.add( w[kk], z );
                else acc1.add( w[kk], z );
            }
            GS_Q(Z1, i, k) = v.bindingWts[i] * acc1.result();
            GS_Q(Z0, i, k) = v.bindingWts[i] * acc0.result();
        }
    }

    T Z_on( 1.0 );
    for ( int i = 1; i <= n; i++ )
    {
        for ( int k = 0; k <= N0; k++ )
        {
            T term = GS_Q(Z1, i, k) + GS_Q(Z0, i, k);
            Z_on += term;
        }
        int act = actIndicators[ sites[i].factorIdx ] ? 1 : 0;
        for ( int k = 0; k <= N0 - 1; k++ )
        {
            Z_on += act * GS_Q(Z1, i, k) * v.txpEffects[ sites[i].factorIdx ];
        }
    }
    #undef GS_Q
    return Z_on;
}

template< class T >
T ExprFunc::kernelOnChrModUnlimited( const ThermoVals< T >& v ) const
{
    int n = n_sites;
    vector< T > Z0( n + 1 );
    Z0[0] = 1.0;
    vector< T > Z1( n + 1 );
    Z1[0] = 1.0;
    vector< T > Zt( n + 1 );
    Zt[0] = 1.0;

    for ( int i = 1; i <= n; i++ )
    {
        T sum = Zt[ boundaries[i] ];
        const vector< SiteInteraction >& nbrs = left_nbrs[i];
        const vector< T >& w = v.w_ji[i];
        const bool rep_i = repIndicators[ sites[i].factorIdx ];
        gemstat_ad::Accumulator< T > acc0( sum );
        for ( size_t k = 0; k < nbrs.size(); k++ )
        {
            int j = nbrs[k].j;
            acc0.add( w[k], Z0[j] );
            if ( nbrs[k].dist > repressionDistThr ) acc0.add( Z1[j] );
        }
        T sum0 = acc0.result();
        T sum1 = sum;
        if ( rep_i )
        {
            gemstat_ad::Accumulator< T > acc1( sum );
            for ( size_t k = 0; k < nbrs.size(); k++ )
            {
                int j = nbrs[k].j;
                acc1.add( w[k], Z1[j] );
                if ( nbrs[k].dist > repressionDistThr ) acc1.add( Z0[j] );
            }
            sum1 = acc1.result();
        }
        Z0[i] = v.bindingWts[i] * v.txpEffects[ sites[i].factorIdx ] * sum0;
        if ( repIndicators[ sites[i].factorIdx ] ) Z1[i] = v.bindingWts[i] * v.repEffects[ sites[i].factorIdx ] * sum1;
        else Z1[i] = 0.0;
        Zt[i] = Z0[i] + Z1[i] + Zt[i - 1];
    }
    return Zt[n];
}

template< class T >
T ExprFunc::kernelOnChrModLimited( const ThermoVals< T >& v ) const
{
    int n = n_sites;
    int N0 = maxContact;
    const int W = N0 + 1;
    vector< T > Z0( ( n + 1 ) * W );
    vector< T > Z1( ( n + 1 ) * W );
    vector< T > Zt( ( n + 1 ) * W );
    #define GS_Q(M, i, k) M[ (i) * W + (k) ]
    GS_Q(Z0, 0, 0) = 0.0;
    GS_Q(Z1, 0, 0) = 0.0;
    GS_Q(Zt, 0, 0) = 1.0;
    for ( int k = 1; k <= N0; k++ )
    {
        GS_Q(Z0, 0, k) = 0.0;
        GS_Q(Z1, 0, k) = 0.0;
        GS_Q(Zt, 0, k) = 0.0;
    }

    for ( int k = 0; k <= N0; k++ )
    {
        for ( int i = 1; i <= n; i++ )
        {
            const vector< SiteInteraction >& nbrs = left_nbrs[i];
            const vector< T >& w = v.w_ji[i];
            const bool rep_i = repIndicators[ sites[i].factorIdx ];

            T sum0, sum0A( 0.0 ), sum1;
            {
                gemstat_ad::Accumulator< T > acc( GS_Q(Zt, boundaries[i], k) );
                for ( size_t kk = 0; kk < nbrs.size(); kk++ )
                {
                    int j = nbrs[kk].j;
                    acc.add( w[kk], GS_Q(Z0, j, k) );
                    if ( nbrs[kk].dist > repressionDistThr ) acc.add( GS_Q(Z1, j, k) );
                }
                sum0 = acc.result();
            }
            if ( k > 0 )
            {
                gemstat_ad::Accumulator< T > acc( GS_Q(Zt, boundaries[i], k - 1) );
                for ( size_t kk = 0; kk < nbrs.size(); kk++ )
                {
                    int j = nbrs[kk].j;
                    acc.add( w[kk], GS_Q(Z0, j, k - 1) );
                    if ( nbrs[kk].dist > repressionDistThr ) acc.add( GS_Q(Z1, j, k - 1) );
                }
                sum0A = acc.result();
            }
            sum1 = GS_Q(Zt, boundaries[i], k);
            if ( rep_i )
            {
                gemstat_ad::Accumulator< T > acc( GS_Q(Zt, boundaries[i], k) );
                for ( size_t kk = 0; kk < nbrs.size(); kk++ )
                {
                    int j = nbrs[kk].j;
                    acc.add( w[kk], GS_Q(Z1, j, k) );
                    if ( nbrs[kk].dist > repressionDistThr ) acc.add( GS_Q(Z0, j, k) );
                }
                sum1 = acc.result();
            }
            GS_Q(Z0, i, k) = v.bindingWts[i] * sum0;
            if ( actIndicators[ sites[i].factorIdx ] && k > 0 ) GS_Q(Z0, i, k) += v.bindingWts[i] * v.txpEffects[ sites[i].factorIdx ] * sum0A;
            if ( rep_i ) GS_Q(Z1, i, k) = v.bindingWts[i] * v.repEffects[ sites[i].factorIdx ] * sum1;
            else GS_Q(Z1, i, k) = 0.0;
            GS_Q(Zt, i, k) = GS_Q(Z0, i, k) + GS_Q(Z1, i, k) + GS_Q(Zt, i - 1, k);
        }
    }

    T total( 0.0 );
    for ( int k = 0; k <= N0; k++ ) total += GS_Q(Zt, n, k);
    #undef GS_Q
    return total;
}

template< class T >
T ExprFunc::kernelLogistic( const ThermoVals< T >& v ) const
{
    // total occupancy of each factor
    vector< T > factorOcc( motifs.size() );
    for ( int i = 1; i <= n_sites; i++ )
    {
        factorOcc[ sites[i].factorIdx ] += v.bindingWts[i] / ( 1.0 + v.bindingWts[i] );
    }
    T totalEffect( 0.0 );
    for ( size_t i = 0; i < motifs.size(); i++ )
    {
        T effect = v.txpEffects[i] * factorOcc[i];
        totalEffect += effect;
    }
    return gemstat_logistic( v.basal + totalEffect );
}

template< class T >
T ExprFunc::kernelMarkov( const ThermoVals< T >& v ) const
{
    int n = n_sites;

    vector< T > Z( n + 2 );
    Z[0] = 1.0;
    vector< T > Zt( n + 2 );
    Zt[0] = 1.0;

    vector< T > backward_Z( n + 2 );
    backward_Z[ backward_Z.size() - 1 ] = 1.0;
    vector< T > backward_Z_sum( n + 1 );
    vector< T > backward_Zt( n + 2 );
    backward_Zt[ backward_Zt.size() - 1 ] = 1.0;

    // recurrence forward
    for ( int i = 1; i <= n; i++ )
    {
        const vector< SiteInteraction >& nbrs = left_nbrs[i];
        const vector< T >& w = v.w_ij[i];
        gemstat_ad::Accumulator< T > acc( Zt[ boundaries[i] ] );
        for ( size_t k = 0; k < nbrs.size(); k++ )
        {
            acc.add( w[k], Z[ nbrs[k].j ] );
        }
        T sum = acc.result();
        Z[ i ] = v.bindingWts[ i ] * sum;
        Zt[i] = Z[i] + Zt[i - 1];
    }

    // recurrence backward (rev_bounds and right_nbrs are built by the Markov_ExprFunc constructor)
    for ( int i = n; i >= 1; i-- )
    {
        const vector< SiteInteraction >& nbrs = right_nbrs[i];
        const vector< T >& w = v.right_w_ij[i];
        gemstat_ad::Accumulator< T > acc( backward_Zt[ rev_bounds[i] ] );
        for ( size_t k = 0; k < nbrs.size(); k++ )
        {
            acc.add( w[k], backward_Z[ nbrs[k].j ] );
        }
        T sum = acc.result();
        backward_Z_sum[i] = sum;
        backward_Z[ i ] = sum * v.bindingWts[i];
        backward_Zt[i] = backward_Z[i] + backward_Zt[i + 1];
    }

    // marginal binding probability of every site, then the expression it implies
    T sum_total( 0.0 );
    for ( int i = 1; i <= n; i++ )
    {
        T one_final_Z = Z[i] * backward_Z_sum[i];
        T bindprob = one_final_Z / backward_Zt[1];
        assert( bindprob >= 0.0 );
        assert( bindprob <= 1.0 );

        T log_effect( 0.0 );
        if ( actIndicators[ sites[ i ].factorIdx ] )
        {
            log_effect = log( v.txpEffects[ sites[ i ].factorIdx ] );
        }
        if ( repIndicators[ sites[ i ].factorIdx ] )
        {
            log_effect = log( v.repEffects[ sites[ i ].factorIdx ] );
        }
        sum_total += log_effect * bindprob;
    }
    T Z_on = exp( sum_total ) * v.basal;
    return Z_on / ( 1.0 + Z_on );
}

template< class T >
T ExprFunc::factorIntAffine( const Site& a, const Site& b, const T& normalInt ) const
{
    double dist = abs( b.start - a.start );
    FactorIntFunc* an_int_func = expr_model->coop_setup->coop_func_for( a.factorIdx, b.factorIdx );
    double scale, offset, post;
    bool clamp_to_one;
    an_int_func->affineForm( dist, a.strand, b.strand, scale, offset, clamp_to_one, post );
    T w = scale * normalInt + offset;
    if ( clamp_to_one ) w = max( 1.0, w );
    return post * w;
}

#endif
