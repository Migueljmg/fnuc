#include "Vec.h"
#include <iostream>
using namespace std;

Vec::Vec() {        //default constructor
 N=0; 
 entries= new double[1];
 entries[0]=0;
}

Vec::Vec(int p) {     //construtor de p elementos todos a 0
  N=p;
  entries=new double[p];
  for(int j=0; j<p; j++)
    entries[j]=0; 
}

Vec::Vec(int p, double a) {    //construtor de p elementos todos a a
  N=p;
  entries=new double[p];
  for(int j=0; j<p; j++)
    entries[j]=a;
}

Vec::Vec(int p, double* a) {   //construtor de p elementos a partir de array a
  N=p;
  entries=new double[p];
  for(int i=0; i<p; i++) {
  entries[i]=a[i];
  }
}

//copy constructors

Vec::Vec(const Vec& a) {
  entries=new double[a.N];
  N=a.N;
  for(int j=0; j<N; j++) {
  entries[j]=a.entries[j];
  }
}

Vec::Vec(const Vec* a) {
  entries=new double[a->N];
  N=a->N;
  for(int j=0; j<N; j++) {
  entries[j]=a->entries[j];
  }
}

void Vec::SetSize(int siz) {
   N=siz;
}

//mudar entries a partir de um array de doubles
void Vec::SEntries(int siz, double* vc) {
  N=siz;
  entries=new double[N];
  for(int j=0; j<N; j++) 
    entries[j]=vc[j];
}

 //mudar entries a partir de um certo valor
void Vec::SEntries(int siz, double vc) { 
    N=siz;
    entries=new double[N];
    for(int j=0; j<siz; j++) 
      entries[j]=vc;
}

 //mudar entries[nm] para um certo valor
void Vec::SEntrie(int nm, double vc) {
  entries[nm]=vc;
}

 //mudar entries a partir de um vector
void Vec::SEntries(int siz,const Vec& vc) {
    N=vc.N;
    entries=new double[N];
    for(int j=0; j<N; j++) 
      entries[j]=vc.entries[j];
}

double* Vec::GetEntries() {          //retorna os valores do vector
  static double* gp;
  gp=new double[N];
  for(int j=0; j<N ; j++) 
    gp[j]=entries[j];

  return gp;
}

double Vec::GetEntrie(int n) {        //retorna a entrada pretendida
  if(n<N)
    return entries[n];
  else {
    cout << "Entrada invalida" << endl;
    return 0;
  }
}

void Vec::Print() {          //Print do vector
  cout <<"<";
  for(int j=0; j<N; j++){
    if((j+1)==((j+1)/2))
      cout << "\n" << endl;
    if(j!=(N-1))
    cout << " " << entries[j]<< " , ";
    else
      cout << "  " << entries[j] << " ";
    }
  cout <<">\n " <<endl ;
}

const Vec& Vec::operator=(const Vec& cv) { //copy assigment
  if (this != &cv) {
    N=cv.N;
    entries=new double[N+1];
    for(int i=0; i<=N ; i++)
      entries[i]=cv.entries[i];
  }
  return *this;
}

const Vec& Vec::operator-(const Vec& cv) {  //Subtraçao de vetores
static Vec cp;
  cp.N=N;
  cp = *this;
  cp -= cv;
  return cp;
}

const Vec& Vec::operator-() { //Negativo de um vector

 for(int i=0;i<N;i++)
   entries[i]=-entries[i];

  return *this;
}

const Vec& Vec::operator-=(const Vec& cv) {
 
if(N<cv.N) {
   int* pc =new int[cv.N+1];
   N=cv.N;
   //copia do vetor extendido
   for(int i=0; i<=cv.N ; i++) {
     if(i<N)
       pc[i]=entries[i];
     else
       pc[i]=0;
   }
   //alocaçao de memoria
   entries =new double[cv.N+1];
   for(int i=0; i<=cv.N ; i++) {
     entries[i]=pc[i];        //transferencia do valor correto para entries
     if(i>N)
	entries[i]=0;
      entries[i] = entries[i]-cv.entries[i];
    }
    return *this;
 }
 else {
   for(int i=0; i<=cv.N ; i++) {
       entries[i]-=cv.entries[i];
   } 
    return *this;
  }
}

const Vec& Vec::operator+=(const Vec& cv) {
  if(N<cv.N) {
    int* pc =new int[cv.N];
    N=cv.N;
    //copia do vetor estendido
   for(int i=0; i<=cv.N ; i++) {
     if(i<=N)
       pc[i]=entries[i];
     else
       pc[i]=0;
   }
   //alocaçao de memoria
   entries =new double[cv.N+1];
   for(int i=0; i<=cv.N ; i++) {
     entries[i]=pc[i];                //transferencia valor certo para entries
      if(i>N)
	entries[i]=0;
      entries[i] = entries[i]+cv.entries[i];
    }
   delete [] pc;
    return *this;
 }
 else {
   for(int i=0; i<=cv.N ; i++) {
       entries[i]+=cv.entries[i];
   } 
    return *this;
  }
}

const Vec& Vec::operator+(const Vec& cv) {  //operaçao soma de vetores
  static Vec cp;
  cp.N=this->N;
  cp.entries=new double[(this->N)+1];
  cp = *this;
  cp += cv;
  return cp;
}

const Vec& Vec::operator+(){               //operador positivo
  static Vec cp;
  cp.N=N;
  cp = *this;
  return cp;
}

const Vec& Vec::operator*(const Vec& cv) { //multiplicaçao por outro vetor
  static Vec cp;
  cp = *this;
  cp *= cv;                        
  return cp;
}


const Vec& Vec::operator*=(const Vec& cv) {   //Definiçao do operador *=
  //caso o vetor seja menor do que o vetor pelo qual multiplica
  if(N<cv.N) {
    int* pc =new int[cv.N];
    N=cv.N;
    //Guardar o valor de entries
    for(int i=0; i<=cv.N ; i++) {
      if(i<=N)
	pc[i]=entries[i];
      else
	pc[i]=0;
    }
    //alocaçao de memoria
    entries =new double[cv.N+1];
    for(int i=0; i<=cv.N ; i++) {
      entries[i]=pc[i];       //transferencia de valor para entries
      if(i>N)
	entries[i]=0;
      entries[i] = entries[i]*cv.entries[i];
    }
 delete [] pc;
    return *this;
 }
 else {                                    
   for(int i=0; i<=cv.N ; i++) {
       entries[i]*=cv.entries[i];
   } 
    return *this;
  }
}

const double& Vec::operator[](const int& s) {
  return entries[s];
}

int Vec::Size() {                        //Retorna o tamanho do vector
  return N;
}

void Vec::Scale(int a) {                 //Multiplicaçao por um escalar
  for (int i=0; i<N ; i++) {
    entries[i]=a*entries[i];
  }
}

double Vec::Dot(Vec p) {                 //Produto interno com um vector p
  double sum;

//Se o numero de elementos dos vectores or diferente nao e possivel realizar produto interno
  if(N != p.N) {                           
   cout << "Operaçao Invalida" << endl;
   return 0;
  }
  else {
    for (int i=0; i<N ; i++) {
      sum=p.entries[i]*entries[i]+sum;
    }
  }
  return sum;
}

Vec:: ~Vec() {
  //delete [] entries;
}
