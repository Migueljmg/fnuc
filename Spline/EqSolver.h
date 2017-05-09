#ifndef __EEQSOLVER__
#define __EEQSOLVER__

#include <iostream>
#include <vector>
#include "Vec.h"
using namespace std;

class EqSolver {
 public:
  EqSolver();

  //eliminação de Gauss: matriz M e vector de constantes B
  //double* GaussElimination(double**, double*, int);
Vec GaussElimination(double**, double*, int);
  //decomposição LU com |L|=1
  void LUdecomposition(double**, int , double**, double**);
  Vec LUdecompositionSolver(double** , double *, int ,double **, double**);
  double** LUdecomposition3(double*, double*, double*, int);
  void LUsolve3(double*, double*, double*, double*, int, double*);
  ~EqSolver();

 private:
  double** m; //matriz (ou outra forma)
  double* b; //vector de constantes (ou outra forma)
  int nn;

};
#endif
