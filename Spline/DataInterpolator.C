#include "DataInterpolator.h"
#include <iostream>
#include "TROOT.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraph.h"
#include <sstream>
#include "Vec.h"
#include "EqSolver.h"
#include "TMath.h"

using namespace std;

DataInterpolator::DataInterpolator( int np, double *xp, double *yp)
{//associa os dados introduzidos aos data members da classe
  N=np;                          
  double* xaux=new double[N];
  double* yaux=new double[N];
  for (int i=0;i<N;++i)
      xaux[i]=xp[i];
  for (int j=0;j<N;++j)
      yaux[j]=yp[j];

  x=xaux;
  y=yaux;

}

double* DataInterpolator::CubicSplineCurvatures()
{//calcula os k's
  double* a=new double[N-2];    //diagonal inferior da matriz dos k's
  double* b=new double[N-2];    //diagonal principal da matriz dos k's
  double* c=new double[N-2];    //diagonal superior da matriz dos k's
  double* g=new double[N-2];    //vector coluna com as constantes
  double* k1=new double[N-2];   //array com os k's excluindo o k[0] e k[N-1], que sao forçados a 0

  a[0]=0;
  c[0]=0;
  //obtençao dos coeficientes da matriz
  for (int i=0;i<N-2;++i)
    {
      b[i]=2*(x[i]-x[i+2]);
      g[i]=6*((y[i]-y[i+1])/(x[i]-x[i+1]) - (y[i+1]-y[i+2])/(x[i+1]-x[i+2]));
    }

  for (int i=1;i<N-2;++i)
    {
      a[i]=x[i]-x[i+1];
      c[i]=a[i];
    }

  EqSolver Sol;
  //resolve a matriz e insere os resultados no array k1
  Sol.LUsolve3(b,c,a,g,N-2,k1);
  
  double* k=new double[N];

  k[0]=0;
  k[N-1]=0;
  for (int i=1;i<N-1;++i)
    k[i]=k1[i-1];                //gera array k com todos os valores

  return k;
}

double CubicAux(double* xp,double* par)
{
  Float_t xx = xp[0];
  double f;
  int N=par[0];
  for(int i=1;i<N;++i)
    {
      if(xx<=par[i+1] && xx>=par[i]) 
	{   
	  f=(par[2*N+i]/6.)*(((xx-par[i+1])*(xx-par[i+1])*(xx-par[i+1])/(par[i]-par[i+1]))-(xx-par[i+1])*(par[i]-par[i+1]))-(par[2*N+i+1]/6.)*(((xx-par[i])*(xx-par[i])*(xx-par[i])/(par[i]-par[i+1]))-(xx-par[i])*(par[i]-par[i+1]))+ ((par[N+i]*(xx-par[i+1])-par[N+i+1]*(xx-par[i]))/(par[i]-par[i+1]));
	}
    }
  return f;
}

TF1* DataInterpolator::CubicSpline(double* k)
{
  double param[3*N+1];
  param[0]=N;                 //iguala par[0] ao numero de pontos

  for(int i=0;i<N;++i)
    param[i+1]=x[i];          //parametros com valores iguais aos de x

  for(int j=0;j<N;++j)
    param[N+j+1]=y[j];       //parametros com valores iguais aos de y

  for(int j=0;j<N;++j)
    param[2*N+j+1]=k[j];     //parametros com valores iguais aos de k

  TF1* func = new TF1("CubicSpline",CubicAux,x[0],x[N-1],3*N+1);
  func->SetParameters(param);
  func->SetTitle("Cubic Spline");

  return func;
}

TF1* DataInterpolator::CubicSplineSegment(double* k, double xp)
{
  double* kp=new double[2];
  int i=0;

  if(xp>x[N-1]||xp<x[0])
    {
      cout << "\nERRO! O ponto introduzido nao pertence ao dominio da funçao\n" << endl;
      cout << "Sera criada a funçao completa!\n" << endl;
      TF1* normal = DataInterpolator::CubicSpline(k);
      return normal;
    }
  else
    {//vai-se procurar o troço a que pertence xp
      for (int j=0;j<=N-1;++j)
	{
	  if(xp==x[j] && xp!=x[N-1])
	    {
	      i=j;      
	      break;
	    }
	  else if(xp==x[N-1])
	    {
	      i=N-2;
	      break;
	    }
	  else if(xp>x[j] && xp<x[j+1])
	    {
	      i=j;
	      break;
	    }
	  else continue;
	}
    
      kp[0]=k[i];
      kp[1]=k[i+1];

      //passagem de doubles para strings
      kp[0]=kp[0]/6;
      std::ostringstream kis;
      kis << kp[0];
      std::string ki= kis.str();

      kp[1]=kp[1]/6;
      std::ostringstream kiups;
      kiups << kp[1];
      std::string kiup = kiups.str();

      std::ostringstream xis;
      xis << x[i];
      std::string xi= xis.str();

      std::ostringstream xiups;
      xiups << x[i+1];
      std::string xiup= xiups.str();

      std::ostringstream yis;
      yis << y[i];
      std::string yi= yis.str();

      std::ostringstream yiups;
      yiups << y[i+1];
      std::string yiup= yiups.str();

      double ddif;
      ddif=x[i]-x[i+1];
      std::ostringstream difs;
      difs << ddif;
      std::string dif= difs.str();

      //string com o segmento da funçao 
      string function = ki + "*(" + "(((x-" + xiup + ")^3)/(" + dif + "))" + "-(x-" + xiup + ")*(" + dif + ")) - " +  kiup + "*(" + "(((x-" + xi + ")^3)/(" + dif + "))" + "-(x-" + xi + ")*(" + dif + ")) + ((" + yi + "*(x-" + xiup + ")-" + yiup + "*(x-" + xi + "))/(" + dif + "))";
      
      //conversao para char*
      char* cfunction = new char[function.size() + 1];
      std::copy(function.begin(), function.end(), cfunction);
      cfunction[function.size()] = '\0';

      TF1* f = new TF1("SFunction", cfunction, x[i],x[i+1]);
 
      f->SetTitle("Spline Segment");

      return f;
    }
}

double DataInterpolator::CSEvaluate(double*k, double z)
{
  int i=0;
  double cse;

  //procura-se o troço a que pertence z
  for (int j=0;j<N-1;++j)
    {
      if(z>x[j] && z<x[j+1] || z==x[j])
	{
	  i=j;
	  break;
	}
      else if(z==x[j+1])
	{
	  i=j+1;
	  break;
	}
      else continue;
    }

  //calculo do valor no ponto
  cse=(k[i]/6)*((z-x[i+1])*(z-x[i+1])*(z-x[i+1])/(x[i]-x[i+1]) - (x[i]-x[i+1])*(z-x[i+1])) - (k[i+1]/6)*((z-x[i])*(z-x[i])*(z-x[i])/(x[i]-x[i+1]) - (z-x[i])*(x[i]-x[i+1])) + (y[i]*(z-x[i+1])-y[i+1]*(z-x[i]))/(x[i]-x[i+1]);

  return cse;
}

double DataInterpolator::CubicSplineDeriv(double* k, double a)
{
  int i=0;
  double df;
  
  if(a>x[N-1] || a<x[0])
    {
      df=0;
      cout << "\nO ponto introduzido nao pertence ao dominio da funçao representada" << endl;
      return df;
    }
  else
    {
      //procura-se o segmento
      for (int j=0;j<=N-1;++j)
	{
	  if(a==x[j] && a!=x[N-1])
	    {
	      i=j;
	      break;
	    }
	  else if(a==x[N-1])
	    {
	      i=N-2;
	      break;
	    }
	  else if(a>x[j] && a<x[j+1])
	    {
	      i=j;
	      break;
	    }
	  else continue;
	}
      //calculo da derivada no ponto
      df=(k[i]/6)*(3*(a-x[i+1])*(a-x[i+1])/(x[i]-x[i+1]) - (x[i]-x[i+1])) - (k[i+1]/6)*(3*(a-x[i])*(a-x[i])/(x[i]-x[i+1]) - (x[i]-x[i+1])) + (y[i]-y[i+1])/(x[i]-x[i+1]);

      return df;
    }
}
double DataInterpolator::CubicSplineDerivN(double* k, double a)
{
  double dfN=0;
  double h=0.00001;
  double fa=DataInterpolator::CSEvaluate(k,a);
  double fuph=DataInterpolator::CSEvaluate(k,a+h);
  double fup2h=DataInterpolator::CSEvaluate(k,a+2*h);
  double fdwh=DataInterpolator::CSEvaluate(k,a-h);
  double fdw2h=DataInterpolator::CSEvaluate(k,a-2*h);

  if(a>x[N-1] || a<x[0])
    {
      dfN=0;
      cout << "\nO ponto introduzido nao pertence ao dominio da funçao representada" << endl;
    }
  else
    {//forward derivative
      if(a==x[0]||(a>x[0] && a<(x[0]+2*h)))
	{
	  dfN=(1/(2*h))*(-1*fup2h + 4*fuph -3*fa);
	}
      else if(a==x[N-1]||(a<x[N-1]&&a>(x[N-1]-2*h)))
	{//backward derivative
	  dfN=(1/(2*h))*(fdw2h -4*fdwh + 3*fa);
	}
      else 
	{//central derivative
	  //dfN=(1/(12*h))*((fdw2h + 8*fuph) - (8*fdwh + fup2h));
	  dfN=(1/(2*h))*(fuph-fdwh);
	}
    }

 return dfN;
}

TF1* DataInterpolator::Polynomial()
{
  string li[N], lagrange;

  for(int i=0;i<N;++i)
    {
      //conversao de doubles em strings
      std::ostringstream xis;
      xis << x[i];
      std::string xi= xis.str();

      std::ostringstream yis;
      yis << y[i];
      std::string yi= yis.str();

      for(int j=0;j<N;++j)
	{
	  std::ostringstream xjs;
	  xjs << x[j];
	  std::string xj= xjs.str();

	  double dqc=x[i]-x[j];
	  std::ostringstream qcs;
	  qcs << dqc;
	  std::string qc= qcs.str();

	 
	  if(j==i) continue;
	  
	  li[i] = "*((x-" + xj + ")/(" + qc + "))"+ li[i]; 
	}

      if(i==0) 
	lagrange = yi + li[i];
      else  
	lagrange = yi + li[i] + '+' + lagrange;
    }
  //conversao para char*
  char* lfunction = new char[lagrange.size() + 1];
  std::copy(lagrange.begin(), lagrange.end(), lfunction);
  lfunction[lagrange.size()] = '\0';

  TF1* fl= new TF1("Function", lfunction, x[0],x[N-1]);
  fl->SetTitle("Polynomial");
  return fl;
}

double DataInterpolator::PEvaluate(double z)
{
  double L=0;
  double *l=new double[N];

  for(int i=0;i<N;++i)
    {
      l[i]=1;
    }
  
  for(int i=0;i<N;++i)
    {
      for(int j=0;j<N;++j)
	{
	  if(j==i) continue;
	    
	  l[i]=((z-x[j])/(x[i]-x[j]))*l[i];
	}
      L=y[i]*l[i]+L;
    }

  return L;
}

double DataInterpolator::PolynomialDerivN(double a)
{
  double PdfN=0;
  double h=0.00001;
  double fa=DataInterpolator::PEvaluate(a);
  double fuph=DataInterpolator::PEvaluate(a+h);
  double fup2h=DataInterpolator::PEvaluate(a+2*h);
  double fdwh=DataInterpolator::PEvaluate(a-h);
  double fdw2h=DataInterpolator::PEvaluate(a-2*h);

  if(a>x[N-1] || a<x[0])
    {
      PdfN=0;
      cout << "\nO ponto introduzido nao pertence ao dominio da funçao representada" << endl;
    }
  else
    {
      if(a==x[0]||(a>x[0] && a<(x[0]+2*h)))
	{//forward derivative
	  PdfN=(1/(2*h))*(-1*fup2h + 4*fuph -3*fa);
	}
      else if(a==x[N-1]||(a<x[N-1] && a<(x[N-1]-2*h)))
	{//backward derivative
	  PdfN=(1/(2*h))*(fdw2h -4*fdwh + 3*fa);
	}
      else 
	{//central diference
	  //PdfN=(1/(12*h))*((fdw2h + 8*fuph) - (8*fdwh + fup2h));
	  PdfN=(1/(2*h))*(fuph-fdwh);
	}
    }

 return PdfN;
}

TGraph* DataInterpolator::Draw()
{
  TGraph* g =new TGraph(N,x,y);
  g->SetTitle("Dots");

  return g;
}

TGraph* DataInterpolator::CubicDerivDif(double* k)
{
  double a=x[0];
  double* ydif=new double[201];
  double* xdif=new double[201];
  //calcula a diferença em cada ponto
  for (int i=0;i<201;++i)
    {
      xdif[i]=a;
      ydif[i]=DataInterpolator::CubicSplineDerivN(k,a)-DataInterpolator::CubicSplineDeriv(k,a);
      a=a+1;
    }

  TGraph* gdif =new TGraph(201,xdif,ydif);
  gdif->SetTitle("Numerical-Analytical CS Derivative Difference");

  return gdif;
}


TGraph* DataInterpolator::DerivNDif(double* k)
{
  double a=x[0];
  double* ydif=new double[201];
  double* xdif=new double[201];
  //calcula a diferença em cada ponto
  for (int i=0;i<201;++i)
    {
      xdif[i]=a;
      ydif[i]=DataInterpolator::CubicSplineDerivN(k,a)-DataInterpolator::PolynomialDerivN(a);
      a=a+1;
    }

  TGraph* gNdif =new TGraph(201,xdif,ydif);
  gNdif->SetTitle("Numerical Derivatives Difference");

  return gNdif;
}

TGraph* DataInterpolator::InterpolDif(double* k)
{
  double a=x[0];
  double* ydif=new double[201];
  double* xdif=new double[201];
  //calcula a difenreça em cada ponto
  for (int i=0;i<201;++i)
    {
      xdif[i]=a;
      ydif[i]=DataInterpolator::CSEvaluate(k,a)-DataInterpolator::PEvaluate(a);
      a=a+1;
    }

  TGraph* gdif =new TGraph(201,xdif,ydif);
  gdif->SetTitle("Interpolation Difference");

  return gdif;
}

