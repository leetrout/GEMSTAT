#include "ObjFunc.h"

double RMSEObjFunc::eval(const vector<vector<double> >& ground_truth, const vector<vector<double> >& prediction,
  const ExprPar* par){

    assert(ground_truth.size() == prediction.size());
    int nSeqs = ground_truth.size();
    int nConds = ground_truth[0].size();
    double squaredErr = 0.0;

    for(int i = 0;i<ground_truth.size();i++){
      double beta = 1.0;
      #ifdef BETAOPTTOGETHER
        if(NULL != par)
          beta = par->getBetaForSeq(i);
        squaredErr += least_square( prediction[i], ground_truth[i], beta, true );
      #else
        squaredErr += least_square( prediction[i], ground_truth[i], beta );
      #endif
  }

    double rmse = sqrt( squaredErr / ( nSeqs * nConds ) );
    return rmse;
}

double AvgCorrObjFunc::eval(const vector<vector<double> >& ground_truth, const vector<vector<double> >& prediction,
  const ExprPar* par){

    assert(ground_truth.size() == prediction.size());
    int nSeqs = ground_truth.size();
    int nConds = ground_truth[0].size();
    double totalSim = 0.0;

    for(int i = 0;i<ground_truth.size();i++){

      totalSim += corr(  prediction[i], ground_truth[i] );
  }

    return -totalSim/nSeqs;
  }

double PGPObjFunc::eval(const vector<vector<double> >& ground_truth, const vector<vector<double> >& prediction,
  const ExprPar* par){

        assert(ground_truth.size() == prediction.size());
        int nSeqs = ground_truth.size();
        int nConds = ground_truth[0].size();
        double totalPGP = 0.0;

        for(int i = 0;i<ground_truth.size();i++){
          double beta = 1.0;
          #ifdef BETAOPTTOGETHER
        	beta = par->getBetaForSeq(i);
                totalPGP += pgp(  prediction[i], ground_truth[i], beta, true);
        	#else
        	totalPGP += pgp(  prediction[i], ground_truth[i], beta );
        	#endif
      }

      return totalPGP / nSeqs;
  }

double AvgCrossCorrObjFunc::exprSimCrossCorr( const vector< double >& x, const vector< double >& y )
  {
      vector< int > shifts;
      for ( int s = -maxShift; s <= maxShift; s++ )
      {
          shifts.push_back( s );
      }

      vector< double > cov;
      vector< double > corr;
      cross_corr( x, y, shifts, cov, corr );
      double result = 0, weightSum = 0;
      //     result = corr[maxShift];
      result = *max_element( corr.begin(), corr.end() );
      //     for ( int i = 0; i < shifts.size(); i++ ) {
      //         double weight = pow( shiftPenalty, abs( shifts[i] ) );
      //         weightSum += weight;
      //         result += weight * corr[i];
      //     }
      //     result /= weightSum;

      return result;
  }

double AvgCrossCorrObjFunc::eval(const vector<vector<double> >& ground_truth, const vector<vector<double> >& prediction,
  const ExprPar* par){

    assert(ground_truth.size() == prediction.size());
    int nSeqs = ground_truth.size();
    int nConds = ground_truth[0].size();
    double totalSim = 0.0;

    for(int i = 0;i<ground_truth.size();i++){
        totalSim += exprSimCrossCorr( prediction[i], ground_truth[i] );
    }

  return -totalSim / nSeqs;
  }


double LogisticRegressionObjFunc::eval(const vector<vector<double> >& ground_truth, const vector<vector<double> >& prediction,
  const ExprPar* par){
    assert(ground_truth.size() == prediction.size());
    int nSeqs = ground_truth.size();
    int nConds = ground_truth[0].size();
    double totalLL = 0.0;

    for(int i = 0;i<ground_truth.size();i++){
      vector<double> Y = ground_truth[i];
      vector<double> Ypred = prediction[i];

      double one_sequence_LL = 0.0;

      for(int j = 0;j<Y.size();j++){
        double one_gt = Y[i];
        double pred_prob = logistic(w*(Ypred[j] - bias));
        double singleLL = one_gt*log(pred_prob) + (1.0 - one_gt)*log(1.0 - pred_prob);
        one_sequence_LL += singleLL;
      }

      totalLL += one_sequence_LL;
    }
    return -totalLL;
}

RegularizedObjFunc::RegularizedObjFunc(ObjFunc* wrapped_obj_func, const ExprPar& centers, const ExprPar& l1, const ExprPar& l2)
{
  my_wrapped_obj_func = wrapped_obj_func;
  ExprPar tmp_energy_space = centers.my_factory->changeSpace(centers, ENERGY_SPACE);
  tmp_energy_space.getRawPars(my_centers );

  //It doesn't matter what space these are in, they are just storage for values.
  l1.getRawPars(lambda1 );
  l2.getRawPars(lambda2 );
  cache_pars = vector<double>(my_centers.size(),0.0);
  //cache_diffs(my_centers.size(),0.0);
  //cache_sq_diffs(my_centers.size(),0.0);
}

double RegularizedObjFunc::eval(const vector<vector<double> >& ground_truth, const vector<vector<double> >& prediction, const ExprPar* par){

  double objective_value = my_wrapped_obj_func->eval( ground_truth, prediction, par );



  double l1_running_total = 0.0;
  double l2_running_total = 0.0;

  ExprPar tmp_energy_space = par->my_factory->changeSpace(*par, ENERGY_SPACE);
  tmp_energy_space.getRawPars(cache_pars );

  for(int i = 0;i<cache_pars.size();i++){
    double the_diff = abs(cache_pars[i] - my_centers[i]);
    l1_running_total += lambda1[i]*the_diff;
    l2_running_total += lambda2[i]*pow(the_diff,2.0);
  }

  objective_value += l1_running_total + l2_running_total;

  return objective_value;
}

double PeakWeightedObjFunc::eval(const vector<vector<double> >& ground_truth, const vector<vector<double> >& prediction,
  const ExprPar* par){
    double on_threshold = 0.5;
    assert(ground_truth.size() == prediction.size());
    int nSeqs = ground_truth.size();
    int nConds = ground_truth[0].size();
    double squaredErr = 0.0;

    for(int i = 0;i<ground_truth.size();i++){
      double beta = 1.0;
      #ifdef BETAOPTTOGETHER
        if(NULL != par)
          beta = par->getBetaForSeq(i);
        squaredErr += wted_least_square( prediction[i], ground_truth[i], beta, on_threshold, true );
      #else
        squaredErr += wted_least_square( prediction[i], ground_truth[i], beta, on_threshold );
      #endif
  }

    double rmse = sqrt( squaredErr / ( nSeqs * nConds ) );
    return rmse;
}

double Weighted_RMSEObjFunc::eval(const vector<vector<double> >& ground_truth, const vector<vector<double> >& prediction,
  const ExprPar* par){
    #ifndef BETAOPTTOGETHER
        assert(false);
    #endif

    assert(ground_truth.size() == prediction.size());
    int nSeqs = ground_truth.size();
    int nConds = ground_truth[0].size();
    double squaredErr = 0.0;

    for(int i = 0;i<ground_truth.size();i++){
      double beta = 1.0;
      if(NULL != par){ beta = par->getBetaForSeq(i); }

      for(int j = 0;j<nConds;j++){
          double single_sqr_error = (beta*prediction[i][j] - ground_truth[i][j]);
          single_sqr_error = weights->getElement(i,j)*pow(single_sqr_error,2);
          squaredErr += single_sqr_error;
      }
    }

    double rmse = sqrt( squaredErr / total_weight );
    return rmse;
}

void Weighted_ObjFunc_Mixin::set_weights(Matrix *in_weights){
    if(NULL != weights){delete weights;}
    weights = in_weights;

    //Caculate the total weight.
    total_weight = 0.0;
    for(int i = 0;i<weights->nRows();i++){
        for(int j = 0;j<weights->nCols();j++){
            total_weight+=weights->getElement(i,j);
        }
    }
}


/*****************************************************
 * Gradients (see ObjFunc::gradient in ObjFunc.h)
 ******************************************************/

// the position, in the flat parameter vector, of the beta that scales sequence i
static int beta_slot_for( const ExprPar* par_index, int i )
{
    if ( NULL == par_index ) return -1;
    return (int)par_index->getBetaForSeq( i );
}

void ObjFunc::gradient(const vector<vector<double> >& ground_truth, const vector<vector<double> >& prediction, const ExprPar* par, const ExprPar* par_index,
                       vector<vector<double> >& d_prediction, vector<double>& d_pars)
{
    // central differences: eval() is cheap compared with a prediction pass
    vector<vector<double> > pred = prediction;
    d_prediction.assign( pred.size(), vector<double>() );
    for ( size_t i = 0; i < pred.size(); i++ )
    {
        d_prediction[i].assign( pred[i].size(), 0.0 );
        for ( size_t j = 0; j < pred[i].size(); j++ )
        {
            double x = pred[i][j];
            double h = 1.0e-6 * ( fabs( x ) > 1.0 ? fabs( x ) : 1.0 );
            pred[i][j] = x + h; double fp = eval( ground_truth, pred, par );
            pred[i][j] = x - h; double fm = eval( ground_truth, pred, par );
            pred[i][j] = x;
            d_prediction[i][j] = ( fp - fm ) / ( 2.0 * h );
        }
    }

    if ( NULL == par ) return;
    vector<double> flat;
    par->getRawPars( flat );
    for ( size_t k = 0; k < flat.size(); k++ )
    {
        double x = flat[k];
        double h = 1.0e-6 * ( fabs( x ) > 1.0 ? fabs( x ) : 1.0 );
        flat[k] = x + h; ExprPar pp = par->my_factory->create_expr_par( flat, PROB_SPACE ); double fp = eval( ground_truth, prediction, &pp );
        flat[k] = x - h; ExprPar pm = par->my_factory->create_expr_par( flat, PROB_SPACE ); double fm = eval( ground_truth, prediction, &pm );
        flat[k] = x;
        d_pars[k] += ( fp - fm ) / ( 2.0 * h );
    }
}

void RMSEObjFunc::gradient(const vector<vector<double> >& ground_truth, const vector<vector<double> >& prediction, const ExprPar* par, const ExprPar* par_index,
                           vector<vector<double> >& d_prediction, vector<double>& d_pars)
{
    // rmse = sqrt( sum_ij ( beta_i p_ij - g_ij )^2 / N )
    int nSeqs = ground_truth.size();
    int nConds = ground_truth[0].size();
    double N = (double)nSeqs * nConds;
    double rmse = eval( ground_truth, prediction, par );
    d_prediction.assign( nSeqs, vector<double>( nConds, 0.0 ) );
    if ( rmse <= 0.0 ) return;   // exact fit: the derivative is not defined, use 0

    for ( int i = 0; i < nSeqs; i++ )
    {
        double beta = 1.0;
        int bslot = -1;
        #ifdef BETAOPTTOGETHER
        if ( NULL != par ) { beta = par->getBetaForSeq(i); bslot = beta_slot_for( par_index, i ); }
        #else
        // the objective solves for the best beta itself; by the envelope theorem
        // the derivative with respect to the predictions is the partial at that beta
        { double num = 0, den = 0;
          for ( int j = 0; j < nConds; j++ ) { num += prediction[i][j] * ground_truth[i][j]; den += prediction[i][j] * prediction[i][j]; }
          beta = num / den; }
        #endif
        double dbeta = 0.0;
        for ( int j = 0; j < nConds; j++ )
        {
            double r = beta * prediction[i][j] - ground_truth[i][j];
            d_prediction[i][j] = beta * r / ( N * rmse );
            dbeta += prediction[i][j] * r / ( N * rmse );
        }
        if ( bslot >= 0 ) d_pars[bslot] += dbeta;
    }
}

void Weighted_RMSEObjFunc::gradient(const vector<vector<double> >& ground_truth, const vector<vector<double> >& prediction, const ExprPar* par, const ExprPar* par_index,
                                    vector<vector<double> >& d_prediction, vector<double>& d_pars)
{
    // rmse = sqrt( sum_ij w_ij ( beta_i p_ij - g_ij )^2 / total_weight )
    int nSeqs = ground_truth.size();
    int nConds = ground_truth[0].size();
    double rmse = eval( ground_truth, prediction, par );
    d_prediction.assign( nSeqs, vector<double>( nConds, 0.0 ) );
    if ( rmse <= 0.0 ) return;

    for ( int i = 0; i < nSeqs; i++ )
    {
        double beta = 1.0;
        int bslot = -1;
        if ( NULL != par ) { beta = par->getBetaForSeq(i); bslot = beta_slot_for( par_index, i ); }
        double dbeta = 0.0;
        for ( int j = 0; j < nConds; j++ )
        {
            double w = weights->getElement(i,j);
            double r = beta * prediction[i][j] - ground_truth[i][j];
            d_prediction[i][j] = w * beta * r / ( total_weight * rmse );
            dbeta += w * prediction[i][j] * r / ( total_weight * rmse );
        }
        if ( bslot >= 0 ) d_pars[bslot] += dbeta;
    }
}

void RegularizedObjFunc::gradient(const vector<vector<double> >& ground_truth, const vector<vector<double> >& prediction, const ExprPar* par, const ExprPar* par_index,
                                  vector<vector<double> >& d_prediction, vector<double>& d_pars)
{
    my_wrapped_obj_func->gradient( ground_truth, prediction, par, par_index, d_prediction, d_pars );

    // penalties are on the ENERGY_SPACE values e_k = log p_k
    vector<double> flat;
    par->getRawPars( flat );
    for ( size_t k = 0; k < flat.size(); k++ )
    {
        double diff = log( flat[k] ) - my_centers[k];
        double sign = diff > 0 ? 1.0 : ( diff < 0 ? -1.0 : 0.0 );
        d_pars[k] += ( lambda1[k] * sign + 2.0 * lambda2[k] * diff ) / flat[k];
    }
}

void AvgCorrObjFunc::gradient(const vector<vector<double> >& ground_truth, const vector<vector<double> >& prediction, const ExprPar* par, const ExprPar* par_index,
                              vector<vector<double> >& d_prediction, vector<double>& d_pars)
{
    // objective = - mean_i corr( prediction_i, ground_truth_i ); no direct parameter dependence
    int nSeqs = ground_truth.size();
    d_prediction.assign( nSeqs, vector<double>() );
    for ( int i = 0; i < nSeqs; i++ )
    {
        const vector<double>& x = prediction[i];
        const vector<double>& y = ground_truth[i];
        int n = x.size();
        d_prediction[i].assign( n, 0.0 );
        double x_bar = mean( x ), y_bar = mean( y );
        double sxx = 0, syy = 0, sxy = 0;
        for ( int j = 0; j < n; j++ )
        {
            sxx += ( x[j] - x_bar ) * ( x[j] - x_bar );
            syy += ( y[j] - y_bar ) * ( y[j] - y_bar );
            sxy += ( x[j] - x_bar ) * ( y[j] - y_bar );
        }
        if ( sxx <= 0.0 || syy <= 0.0 ) continue;   // corr undefined: eval() returns NaN there anyway
        double denom = sqrt( sxx * syy );
        double corr_xy = sxy / denom;
        for ( int j = 0; j < n; j++ )
        {
            double dcorr = ( y[j] - y_bar ) / denom - corr_xy * ( x[j] - x_bar ) / sxx;
            d_prediction[i][j] = -dcorr / nSeqs;
        }
    }
}
