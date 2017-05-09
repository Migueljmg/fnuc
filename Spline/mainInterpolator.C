#include <iostream>
#include <string>
#include "TROOT.h"
#include "TH1F.h"
#include "TF1.h"
#include "DataInterpolator.h"
#include "Vec.h"
#include "EqSolver.h"
#include "TGraph.h"
#include "TLatex.h"

#include "TSpline.h"
#include "TVirtualPad.h"
#include "TSystem.h"
#include "Riostream.h"
#include "TClass.h"
#include "TMath.h"

#include "TFrame.h"  //gClient
#include "TCanvas.h"  //Ecra principal dos graficos
#include "TPad.h"  //"Subcanvas"
#include "TFrame.h"  //gClient
#include "TCanvas.h"  //Ecra principal dos graficos
#include "TApplication.h"  //Janela
#include "TAxis.h"  //Eixos dos graficos
#include "TMultiGraph.h"  //Varios graficos sobrepostos
#include "TList.h"

using namespace std;

int main(int argc, char **argv)
{
  int N=13;
  double x[13];
  double y[13];

  //Janela
  TApplication theApp("App", &argc, argv);
  theApp.InitializeGraphics();
  
  TCanvas *c1 = new TCanvas("c1","Spline",200,10,1280,720);
  c1->SetFillColor(0);
  //c1->SetGrid();
  c1->GetFrame()->SetFillColor(21);
  c1->GetFrame()->SetBorderSize(12);

  //Pontos
  y[0]=1700;
  y[1]=1750;
  y[2]=1800;
  y[3]=1850;
  y[4]=1900;
  y[5]=1950;
  y[6]=2000;
  y[7]=2050;
  y[8]=2100;
  y[9]=2150;
  y[10]=2200;
  y[11]=2250;
  y[12]=2300;

  x[0]=4379.87;
  x[1]=4442.02;
  x[2]=4503.29;
  x[3]=4563.72;
  x[4]=4623.33;
  x[5]=4682.17;
  x[6]=4740.25;
  x[7]=4797.62;
  x[8]=4854.28;
  x[9]=4910.28;
  x[10]=4965.63;
  x[11]=5020.36;
  x[12]=5074.48;

  DataInterpolator A((int)13,x,y);

  double* k=new double[N];
  k=A.CubicSplineCurvatures();  
  /* 
  //alinea a)
  double ce;
  ce=A.CSEvaluate(k,4.5);
 
  cout << "\n Valor da interpolaçao por Cubic Spline em E= MeV: " << ce << " (mbarn)"<< endl;
  */

  //alinea b)
  TGraph* g =A.Draw();
  g->SetMarkerStyle(7);
  g->SetTitle("Energia em funcao de B");
  g->GetYaxis()->SetTitle("E (keV)");
  g->GetXaxis()->SetTitle("B (T)");
  g->Draw("AP");
  c1->Update();
  
  TF1* function = A.CubicSpline(k);
  function->SetLineColor(kCyan+3);
  function->SetLineWidth(2);
  function->Draw("SAME");
  c1->Update();
  
  cout << "\n Valor da funçao eval: " << function->Eval(4742.8) << endl;

  c1->Update();
  c1->Modified();
  c1->Draw();
  c1->Print("EvsB.pdf");
  getchar();

  theApp.Terminate();

  delete [] k;

  return 0;
}
