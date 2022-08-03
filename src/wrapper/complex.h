#ifndef COMPLEX_H
#define COMPLEX_H

#include <fstream>

// Imaginary unit
#ifndef STRUCTARRAYPROB
#define I_c (complex_t){0.,1.}
#define ONE_c (complex_t){1.,0.}
#define ZERO_c (complex_t){0.,0.}
#endif

// From standard glibc headers
#ifndef M_PIl
# define M_PIl		3.1415926535897932384626433832795029L  /* pi */
#endif


/** 
Encapsulates real an imaginary part as a simple double array.
*/
typedef struct complex_t{
  double re;
  double im;
} complex_t;

complex_t operator+(complex_t a,complex_t b);
complex_t& operator+=(complex_t& a,complex_t b);
complex_t operator-(complex_t a);
complex_t operator-(complex_t a,complex_t b);
complex_t operator*(complex_t a,complex_t b);
complex_t& operator*=(complex_t& a,complex_t b);
complex_t operator*(double x,complex_t a);
complex_t operator*(complex_t a,double x);
complex_t operator/(complex_t a,complex_t b);
complex_t operator/(complex_t a,double x);
complex_t conjugate(complex_t a);
bool operator==(complex_t a,complex_t b);
bool operator!=(complex_t a,complex_t b);
std::ostream& operator<<(std::ostream& os,complex_t c);
std::istream& operator>>(std::istream& is,complex_t& c);
double abs(complex_t c);
double phase(complex_t c);
double real(complex_t c);
double imag(complex_t c);
complex_t sqrt(complex_t c);
complex_t cos(complex_t c);
complex_t sin(complex_t c);
complex_t cosh(complex_t c);
complex_t sinh(complex_t c);
complex_t exp(complex_t c);
complex_t power(complex_t c,double k);

#ifdef STRUCTARRAYPROB
static complex_t I_c={0.,1.};
static complex_t ONE_c={1.,0.};
static complex_t ZERO_c={0.,0.};
#endif

#endif // COMPLEX_H
