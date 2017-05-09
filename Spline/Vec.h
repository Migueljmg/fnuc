#ifndef __VEC__
#define __VEC__

#include <iostream>

using namespace std;

class Vec {
 public:
  Vec();
  Vec(int);
  Vec(int, double);
  Vec(int, double*);
  Vec(const Vec& );
  Vec(const Vec* );

  void SetSize(int);
  void SEntries(int, double*);
  void SEntries(int, double);
  void SEntrie(int, double);
  void SEntries(int, const Vec&);
  double* GetEntries();
  double GetEntrie(int);
  void Print(); 

  const Vec& operator+(const Vec&); 
  const Vec& operator-(const Vec&);
  const Vec& operator-();
  const Vec& operator+();
  const Vec& operator*=(const Vec&);
  const Vec& operator=(const Vec&);
  const Vec& operator-=(const Vec&);
  const Vec& operator+=(const Vec&);
  const double& operator[](const int&);
  const Vec& operator*(const Vec&);

  int Size();
  void Scale(int);
  double Dot(Vec);
  ~Vec();

 private:
  int N; //number of elements
  double *entries; // pointer to array of doubles
};
#endif
