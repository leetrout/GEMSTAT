#include "ExprFunc_rates.h"

double Rates_ExprFunc::predictExpr( const vector< double >& factorConcs )
{
    setupBindingWeights( factorConcs );
    return kernelRates( plain );
}

gemstat_ad_t Rates_ExprFunc::predictExprAD( ThermoVals< gemstat_ad_t >& vals, const vector< double >& factorConcs ) const
{
    setupBindingWeightsT( vals, factorConcs );
    return kernelRates( vals );
}
