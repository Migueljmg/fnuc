#include "EqSolver.h"
#include <iostream>
#include <vector>
#include <cmath>
#include "Vec.h"

using namespace std;
EqSolver::EqSolver() {
  nn=0;
}

//eliminação de Gauss: matriz M e vector de constantes B e resolução do sistema conjunta
//escolheu-se alterar o retorno da função de modo a ser possivel resolver o sistema aquando da eliminção de gauss

//double* EqSolver::GaussElimination(double** M, double* B, int n) {
Vec EqSolver::GaussElimination(double** M, double* B, int n) {
 //Determinant Expansion by Minors
  double **A =new double*[n];
  double det =0;
  Vec s(n);

  for(int i=0; i<n-1; i++) {
    A[i]=new double [n];
    for(int j=0; j<n; j++) 
      A[i][j]=M[i+1][j];
    det=det + pow(-1,i)*M[i][0]*A[i][0];
  }
  
  if(det == 0)
    cout << "Não tem solução única o sistema considerado" << endl;
  else {
 //Ver se pivot e nulo ou seja mudança de linha necessária

    /*      double p = M[i][i];
      int indice= i;
	
      for(int l=i+1; l<n; l++) {
	if(p > abs(M[l][i])) {
	  indice = l;
	  p = M[l][i];
	}
      }
      //Row exchange
      for(int l = 0; l<n; l++) {
	double t = M[i][l];
	M[i][l] = M[indice][l];
	M[indice][l] = t;
      }
    */

	//Elimination
    for(int i=0; i<n-1; i++) {
      for(int j=i+1; j<n; j++) {
	double c = -(M[j][i]/M[i][i]);
	for(int l = i; l<n; l++) {
	  M[j][l]=M[j][l]+c*M[i][l];
	}
	B[j] = B[j]+ c*B[i];
      }
    }
    
    //Guardar na classe, a matriz de gauss e o vetor b apos eliminaçao de gauss
    nn=n;
    b= new double[n];
    m=new double*[n];
    for(int i=0; i<n; i++)  {
      b[i]=B[i];
      m[i]=new double[n];
      for(int j=0; j<n; j++) {
	m[i][j]=M[i][j];
      }
    }
    
    double* x=new double[n];

    x[n-1]=B[n-1]/M[n-1][n-1];
    double soma = 0;
    for(int i = n-2; i>=0; i--) {
      soma = 0;
      for(int j = i; j<n; j++) 
	soma = soma+ x[j]*M[i][j];
      x[i]=(B[i]-soma)/M[i][i];
    } 
    
    for(int i=0; i<n; i++) {
      if(i==0)
	cout << "X[" << n << "] = \t|" << x[i] << "|\t" << "B[" << n << "] = \t|" << B[i] <<  "|\t ";
      else
	cout << "\t|" << x[i] << "|\t\t|" << B[i] << "|\t ";

      for(int j=0; j<n;j++) {
	if(i==0 && j==0) 
	  cout << " M[" << n << "][" << n << "] = \t|";
	else if (j==0)
	  cout << "\t\t|";
	if((j+1)%4)
	  cout << M[i][j] << "   ";
	else 
	  cout << M[i][j] << "|\n";
      }
    }
    cout << "\n" ;
    s.SEntries(n, x);
    }

  return s;
}
 
//decomposição LU com |L|=1
void EqSolver::LUdecomposition(double** M, int n, double** U, double** L){

  for(int i=0; i<n; i++) {
    L[i][i]=1;
    U[0][i]=M[0][i];
  }

    for(int i=0; i<n-1; i++) {
      for(int j=i+1; j<n; j++) {
	L[j][i]=(M[j][i]/M[i][i]);
	for(int l = i; l<n; l++) {
	  M[j][l]=M[j][l]-L[j][i]*M[i][l];
	  U[j][l]=M[j][l];
	}
      }
    }
    
    for(int i=0; i<n; i++) {
      for(int j=0; j<n;j++) {
	if(i==0 && j==0) 
	  cout << "U[" << n << "][" << n << "] = \t|";
	else if (j==0)
	  cout << "\t\t|";
	if((j+1)%4)
	  cout << U[i][j] << "   ";
	else 
	  cout << U[i][j] << "|\n";
      }
    }

    for(int i=0; i<n; i++) {
      for(int j=0; j<n;j++) {
	if(i==0 && j==0) 
	  cout << " \nL[" << n << "][" << n << "] = \t|";
	else if (j==0)
	  cout << "\t\t|";
	if((j+1)%4)
	  cout << L[i][j] << "   ";
	else 
	  cout << L[i][j] << "|\n";
      }
    }
    
}

Vec EqSolver::LUdecompositionSolver(double** M, double *B, int n, double** U, double** L){
  
  LUdecomposition(M, n, U,L);
  double* y0= new double[n];
  double soma=0;

  y0[0]=B[0];
  for(int i=1;i<n; i++) {
    y0[i]=B[i];
    for(int j=0; j<i; j++) 
      y0[i]-=L[i][j]*y0[j];
    }

  double* x= new double[n];

    x[n-1]=y0[n-1]/U[n-1][n-1];
    soma = 0;
    for(int i = n-2; i>=0; i--) {
      soma = 0;
      for(int j = i; j<n; j++) 
	soma = soma+ x[j]*U[i][j];
      x[i]=(y0[i]-soma)/U[i][i];
    } 

    Vec s(n,x);
    return s;
}

double** EqSolver::LUdecomposition3(double* dp, double* du, double* dl, int n){
  double* v = new double[n];
  double* l = new double[n];
  double** ar = new double*[3];

  v[0]=dp[0];
  
  for(int i=1; i<n ; i++) {
    l[i]=dl[i]/v[i-1];
    v[i]=dp[i]-l[i]*du[i];
  }
  ar[0]=l;
  ar[1]=v;
  ar[2]=du;

  return ar;
}  

void EqSolver::LUsolve3(double* dp, double* du, double* dl, double* f, int n, double* x){
  //Au=f    A=LU
  //Solve Ly=f
  double** ar = new double*[3];

  for(int i=0; i<3; i++)
    ar[i]=new double[n];

  double* y = new double[n];

  ar=LUdecomposition3(dp, du, dl, n);

  y[0]=f[0];
  for(int i=1; i<n; i++) {
    y[i]=f[i]-ar[0][i]*y[i-1]; //i=1...n
  }

  //Solve Uu=y

  x[n-1]=y[n-1]/ar[1][n-1]; 

  /*
    for(int i=(n-2); i>=0; i--) 
    x[i]=(y[i]-ar[2][i]*x[i+1])/ar[1][i]; //i=0...n
  */

  for(int i=(n-2); i>=0; i--) 
    x[i]=(y[i]-ar[2][i+1]*x[i+1])/ar[1][i]; //i=0...n-1
}

EqSolver::~EqSolver(){;}
  
