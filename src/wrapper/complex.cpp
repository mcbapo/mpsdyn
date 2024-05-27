#include "complex.h"
#include <cmath>

using namespace std;


complex_t operator+(complex_t a,complex_t b){
  complex_t c;
  c.re=a.re+b.re;c.im=a.im+b.im;
  return c;
}

complex_t& operator+=(complex_t& a,complex_t b){
  a.re+=b.re;
  a.im+=b.im;
  return a;
}

complex_t operator-(complex_t a){
  complex_t c={-a.re,-a.im};
  return c;
}
complex_t operator-(complex_t a,complex_t b){
  complex_t c={a.re-b.re,a.im-b.im};
  return c;
}
complex_t operator*(complex_t a,complex_t b){
  complex_t c={a.re*b.re-a.im*b.im,a.re*b.im+a.im*b.re};
  return c;
}

complex_t& operator*=(complex_t& a,complex_t b){
  double rA=a.re;
  a.re=a.re*b.re-a.im*b.im;
  a.im=rA*b.im+a.im*b.re;
  return a;
}


complex_t operator*(double x,complex_t a){
  complex_t c={x*a.re,x*a.im};
  return c;
}

complex_t operator*(complex_t a,double x){
  return x*a;
}

bool operator==(complex_t a,complex_t b){
  return (a.re==b.re)&&(a.im==b.im);
}

bool operator!=(complex_t a,complex_t b){
  return (a.re!=b.re)||(a.im!=b.im);
}

complex_t conjugate(complex_t a){
  complex_t c={a.re,-a.im};
  return c;
}

complex_t operator/(complex_t a,complex_t b){
  return (1./(b.re*b.re+b.im*b.im))*(a*conjugate(b) );
}

complex_t operator/(complex_t a,double x){
  return (1./x)*a;
}

#include <iostream>
ostream& operator<<(ostream& os,complex_t c){
  os<<c.re;
  if(c.im>=0) os<<"+"<<c.im<<"i";
  else os<<"-"<<-c.im<<"i";
  return os;
}

#include <cstdlib>

istream& operator>>(istream& is,complex_t& c){
  string aux;
  is>>aux;
  // Now I have to parse it!
  int epos=aux.find_first_of("e");
  int ipos=aux.find_first_of("i");
  int signpos=aux.substr(1).find_first_of("+-") +1;
  if (signpos == epos+1){
    signpos=aux.substr(epos+2).find_first_of("+-") + epos+2;
  }
  // The real part goes from 0 to signpos-1, and the imag, from
  // signpos until ipos-1 (including the sign symbol)
  string real=aux.substr(0,signpos-1+1);
  string imag=aux.substr(signpos,ipos-1-signpos+1);
  //cout<<"Read from "<<aux<<", real part "<<real<<" and imag "<<imag<<endl;
  c.re=atof(real.c_str());
  c.im=atof(imag.c_str());
  //cout<<"Now in numbers, c is "<<c<<endl;
  // cout << "String=" << aux << " epos=" << epos << " signpos=" << signpos << " ipos=" << ipos << " real=" << c.re << " imag=" << c.im << endl;
  return is;
}

double abs(complex_t c){
  return sqrt(c.re*c.re+c.im*c.im);
}

double phase(complex_t c){
  //cout<<"Phase of c "<<c<<", i.e re="<<c.re<<", im="<<c.im<<endl;
  if(c.im==0)
    if(c.re>=0) return 0;
    else return M_PIl;
  if(c.re==0)
    if(c.im>0) return .5*M_PIl;
    else return -.5*M_PIl;
  // none is zero
  double ph=atan(c.im/c.re);
  if(c.re<0)
    if(c.im>0) return ph+M_PIl;
    else return ph-M_PIl;
  return ph;
  // WRONG!!
  //  return c.im==0?0:(c.re==0?(c.im>0?.5*M_PIl:(c.im<0?-.5*M_PIl:0)):atan(c.im/c.re));
}

double real(complex_t c){
  return c.re;
}
double imag(complex_t c){
  return c.im;
}

complex_t sqrt(complex_t c){
  double ph=phase(c);
  complex_t d={cos(ph/2),sin(ph/2)};
  return sqrt(abs(c))*d;
}

complex_t cos(complex_t c){
  complex_t d={cos(c.re)*cosh(c.im),-sin(c.re)*sinh(c.im)};
  return d;
}

complex_t sin(complex_t c){
  complex_t d={sin(c.re)*cosh(c.im),cos(c.re)*sinh(c.im)};
  return d;
}

complex_t cosh(complex_t c){
  complex_t d={cos(c.im)*cosh(c.re),+sin(c.im)*sinh(c.re)};
  return d;
}

complex_t sinh(complex_t c){
  complex_t d={cos(c.im)*sinh(c.re),sin(c.im)*cosh(c.re)};
  return d;
}

complex_t exp(complex_t c){
  double modul=exp(c.re);
  complex_t d={modul*cos(c.im),modul*sin(c.im)};
  return d;
}

complex_t power(complex_t c,double x){
  double z=abs(c);
  double phi=phase(c);
  z=pow(z,x);
  phi=phi*x;
  complex_t res;
  res.re=z*cos(phi);res.im=z*sin(phi);
  return res;
}
