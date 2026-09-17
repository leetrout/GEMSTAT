#ifndef OBJFUNC_H
#define OBJFUNC_H

#include "ExprPar.h"

class Weighted_ObjFunc_Mixin {
  public:
    Weighted_ObjFunc_Mixin() : weights(NULL) {}
    ~Weighted_ObjFunc_Mixin(){if(NULL != weights){delete weights;}}

    virtual void set_weights(Matrix *in_weights);//Some child classes might autocreate weights
    virtual Matrix* get_weights(){ return weights;}
    virtual double get_total_weight(){return total_weight;}
  protected:
      Matrix *weights;
      double total_weight;
};

class ObjFunc {
public:
  ObjFunc(){}
  virtual ~ObjFunc(){}

public:
  //pure virtual function
  virtual double eval(const vector<vector<double> >& ground_truth, const vector<vector<double> >& prediction, const ExprPar* par){return -1;}

  /*
   * Derivatives of eval() with respect to its inputs, for the reverse-mode
   * gradient of the whole objective (ExprPredictor::gradient_prob):
   *   d_prediction[i][j] = d eval / d prediction[i][j]
   *   d_pars[k]         += d eval / d p_k   for the *direct* dependence on the
   *                        flat PROB_SPACE parameter vector p of par (the
   *                        per-sequence scale beta, regularisation terms);
   *                        the dependence through the predictions is not
   *                        included.  d_pars must already have par's length.
   *   par_index is the index parameter set (ParFactory::create_index_par())
   *                        used to find a parameter's position in d_pars.
   * The default is a central-difference approximation; objectives with a
   * simple closed form override it.
   */
  virtual void gradient(const vector<vector<double> >& ground_truth, const vector<vector<double> >& prediction, const ExprPar* par, const ExprPar* par_index,
                        vector<vector<double> >& d_prediction, vector<double>& d_pars);
};


class RMSEObjFunc: public ObjFunc {
public:
  ~RMSEObjFunc(){}
  double eval(const vector<vector<double> >& ground_truth, const vector<vector<double> >& prediction, const ExprPar* par);
  void gradient(const vector<vector<double> >& ground_truth, const vector<vector<double> >& prediction, const ExprPar* par, const ExprPar* par_index,
                vector<vector<double> >& d_prediction, vector<double>& d_pars);
};

class AvgCorrObjFunc: public ObjFunc {
public:
  ~AvgCorrObjFunc(){}
  double eval(const vector<vector<double> >& ground_truth, const vector<vector<double> >& prediction, const ExprPar* par);
  void gradient(const vector<vector<double> >& ground_truth, const vector<vector<double> >& prediction, const ExprPar* par, const ExprPar* par_index,
                vector<vector<double> >& d_prediction, vector<double>& d_pars);
};

class PGPObjFunc: public ObjFunc {
public:
  ~PGPObjFunc(){}
  double eval(const vector<vector<double> >& ground_truth, const vector<vector<double> >& prediction, const ExprPar* par);
};

class AvgCrossCorrObjFunc: public ObjFunc {
public:
  AvgCrossCorrObjFunc(int _maxShift, double _shiftPenalty){  maxShift = _maxShift; shiftPenalty = _shiftPenalty; }
  ~AvgCrossCorrObjFunc(){}
  double exprSimCrossCorr( const vector< double >& x, const vector< double >& y );
  double eval(const vector<vector<double> >& ground_truth, const vector<vector<double> >& prediction, const ExprPar* par);
private:
  int maxShift;
  double shiftPenalty;
};

class LogisticRegressionObjFunc: public ObjFunc {
public:
  LogisticRegressionObjFunc(){ w = 1.0; bias = 0.5;}
  LogisticRegressionObjFunc(double _w, double _bias){w = _w; bias = _bias;}
  ~LogisticRegressionObjFunc(){}
  double eval(const vector<vector<double> >& ground_truth, const vector<vector<double> >& prediction, const ExprPar* par);
private:
  double w;
  double bias;
};

//TODO: Make this class also inherit from Weighted_ObjFunc_Mixin and explicitly calculate a weight matrix.
class PeakWeightedObjFunc: public ObjFunc {
public:
  ~PeakWeightedObjFunc(){}
  double eval(const vector<vector<double> >& ground_truth, const vector<vector<double> >& prediction, const ExprPar* par);
};

class Weighted_RMSEObjFunc: public RMSEObjFunc, public Weighted_ObjFunc_Mixin {
public:
    Weighted_RMSEObjFunc() : RMSEObjFunc(), Weighted_ObjFunc_Mixin() {}
    ~Weighted_RMSEObjFunc(){};
  double eval(const vector<vector<double> >& ground_truth, const vector<vector<double> >& prediction, const ExprPar* par);
  void gradient(const vector<vector<double> >& ground_truth, const vector<vector<double> >& prediction, const ExprPar* par, const ExprPar* par_index,
                vector<vector<double> >& d_prediction, vector<double>& d_pars);

};

#endif
