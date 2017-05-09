#ifndef ___DATAIP___
#define ___DATAIP___

#include "TROOT.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraph.h"

using namespace std;

class DataInterpolator {
 public:
  DataInterpolator (int, double*, double*);    //construtor

  double* CubicSplineCurvatures();
  TF1* CubicSpline(double*);                   //retorna TF1 com a funçao interpoladora por cubic spline
  TF1* CubicSplineSegment(double*, double);    //retorna o segmento que contem um dado x
  double CSEvaluate(double*,double);           //retorna o valor de f(x) atraves do cubic spline
  double CubicSplineDeriv(double*, double);    //retorna a derivada em x do cubic spline
  double CubicSplineDerivN(double*, double);   //retorna a derivada numerica em x do cubic spline
  TF1* Polynomial();                           //retorna TF1 com o polinomio interpolador de lagrange
  double PEvaluate(double);                    //retona o valor de p(x) atraves do polinomio de lagrange
  double PolynomialDerivN(double);             //retorna a derivada numerica em x do polinomio de lagrange
  TGraph* Draw();                              //grafico com os pontos
  TGraph* CubicDerivDif(double*);              //grafico com os pontos da diferença entre as derivadas do cubic spline
  TGraph* DerivNDif(double*);                  //grafico com os pontos da diferença entre a derivada numerica do cubic spline e do pol. de lagrange
  TGraph* InterpolDif(double*);                //grafico com os pontos da diferença entre os dois interpoladores

  ~DataInterpolator() {;}                      //destrutor

 private:
  int N;           //numero de pontos
  double *x, *y;   //array com pontos x e outro com pontos y 

};

#endif
