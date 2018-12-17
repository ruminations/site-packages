#cython: nonecheck=False, boundscheck=False
"""
Package: algebra
Module : mathx.pyx
Website: https://github.com/ruminations/site-packages
License: https://github.com/ruminations/Licenses#design-license
Initial Copyright 2018

A math utility in the spirit of the builtin 'math' module extending
available functionality in various ways.

Fast implementation of some useful auxilliary math functions:
   ferf(x)    : a fast approximation of the error function
   fcdf(x)    : a fast approximation of the cumulative distribution function
   cbrt(x)    : cube root function
   qntrt(x)   : fifth root function

Several simple convenience functions:
   even(n)    : return 'True' if integer 'n' is even
   odd(n)     : return 'True' if integer 'n' is odd
   sign(x,n)  : return 'copysign(1.0,x)', forced to 0.0 for -1e-n < x < 1e-n
   zero(x,n)  : return 'x', forced to 0.0 for -1e-n < x < 1e-n
   clean(x,n) : return 'x' rounded to 'n' significant decimal digits
   len(obj)   : return 'obj.__len__()', permitting lengths unsanitized by python

Some combinatorial functions :
   factorial(n)     : n!, limited to base architecture integers
   permutation(n,k) : n!/k!, limited, but efficiently computed
   combination(n,k) : n!/(k!(n-k)!),  limited, but efficiently computed

Two booleans applicable to any objects implementing arithmetic:
   delta_is_zero(a,b,n)       : return 'True' if
                                  0 <= abs(a-b) <= 1e-n
   delta_ratio_is_zero(a,b,n) : return 'True' if
                                  0 <= abs(a-b)/max(abs(a),abs(b)) <= 1e-n

The error function algorithm is documented at:
  https://github.com/ruminations/Essays/Approximating_erf(x).pdf
The algorithm implemented as a C extension module is theoretically about
20x as fast as 'math.erf'; for comparison, the Bürmann approximation is
theoretically about 10x as fast.  Implemented as here in cython, in the
typical use context involving python invocation overhead, the algorithm
is about 3x as fast as 'math.erf'.  Accuracy is about 10x that of the
Bürmann approximation, but still many orders of magnitude less than
'math.erf', which has double precision accuracy.

The Kullback-Leibler divergence of the base distribution of 'ferf' from
the base distribution of 'math.erf' calculated over the domain:
    [-6.0000 .. 6.0001 (dx == .0001)]
is about:
    6.413910846374e-05
The divergence of the base distribution of 'math.erf' from the base
distribution of 'ferf' over the same domain is:
    5.694665020225e-05
The approximation is therefore quite accurate.

The n-th root algorithm is:
  1) guess root == 1.0
  2) improve guess by replacing with mean of n values (for nth root)
     with the loop invariant: <product of those n values> == x
This converges slightly faster than Newton's method.
"""
__version__=20180422
__version_history__=(20180325,20180404)
__all__=['ferf','fcdf','cbrt','qntrt', \
         'even','odd','len','sign','clean','zero', \
         'factorial','permutation','combination', \
         'delta_is_zero','delta_ratio_is_zero']
cdef extern from "math.h":
  ##cdef double erf( double ) ## excised after testing
  cdef double fabs( double )
  cdef double ceil( double )
  cdef double sqrt( double )
  cdef double tanh( double )
  cdef double log10( double )
  cdef double scalbn( double, int )        #   scalbn(x,n) == x*2.0^n
  cdef double copysign( double, double )   # copysign(x,y) == x with sign(y)
  double M_PI       # pi           == 3.14159265358979323846
  double M_2_SQRTPI # 2.0/sqrt(pi) == 1.12837916709551257390
  double M_SQRT1_2  # 1.0/sqrt(2)  == 0.70710678118654752440
cdef extern from "float.h":
  double DBL_MIN # double precision minimum value
  double DBL_EPSILON # double precision minimum value
# constants

cdef double NM_EPSILON  = 5.0*DBL_EPSILON         # Newton's method termination
cdef double A           = 0.1009107485756095      # for ferf, empirically determined
cdef double LOG_ARG_MIN = 2.0*DBL_MIN*DBL_EPSILON # avoid infinite logarithm

# C functions

cdef double c_cbrt(double x):
  cdef double a
  cdef double b
  cdef double c
  if x==0.0: return 0.0
  a=1.0; b=c=fabs(x) # compute positive, return copysign
  while fabs(a-b)>NM_EPSILON:
    b=(a+scalbn(sqrt(c/a),1))/3.0 # (a+2.0*sqrt(abs(x)/a))/3.0
    a,b=b,a
  return copysign(a,x)

cdef double c_qntrt(double x):
  cdef double a
  cdef double b
  cdef double c
  if x==0.0: return 0.0
  a=1.0; b=c=fabs(x) # compute positive, return copysign
  while fabs(a-b)>NM_EPSILON:
    b=(a+scalbn(sqrt(sqrt(c/a)),2))/5.0 # (a+4.0*sqrt(sqrt(abs(x)/a)))/5.0
    a,b=b,a
  return copysign(a,x)

cdef double c_ferf(double x):
  return tanh(x*(A*x*x+M_2_SQRTPI)) # tanh( A*x^3 + x*2/sqrt(pi) )

cdef double c_fcdf(double x):
  return scalbn(1.0+c_ferf(M_SQRT1_2*x),-1) # .5*(1.0+c_ferf(sqrt(2)*x))

cdef int c_factorial(int n):
  cdef int f
  cdef int i
  f=1
  if n<0: return -1 # bad arguments
  for i in range(n,0,-1): f*=i
  return f

cdef int c_permutation(int n, int k):
  cdef int d
  cdef int p
  cdef int i
  d,p=n-k,1
  if k<0 or n<0 or d<0: return -1 # bad arguments
  for i in range(n,d,-1): p*=i
  return p

cdef int c_combination(int n, int k):
  cdef int p
  cdef int i
  p=c_permutation(n,k)
  if p<0: return -1 # bad arguments
  for i in range(k,0,-1): p/=i
  return p

cdef double c_zero(double x, int p):
  cdef double y
  y=fabs(x)
  if y <= LOG_ARG_MIN: return 0.0
  if p+int(ceil(log10(y))) <= 0: return 0.0
  return x

## excised after testing:
"""
cdef c_execute_erf(double x,int iterations):
  cdef double r
  for i in xrange(iterations): r=erf(x)
cdef c_execute_ferf(double x,int iterations):
  cdef double r
  for i in xrange(iterations): r=c_ferf(x)
#"""

# python wrappers

def cbrt(double x):
  """Return the cube root of 'x'."""
  return c_cbrt(float(x))

def qntrt(double x):
  """Return the fifth root of 'x'."""
  return c_qntrt(float(x))

def ferf(double x):
  """Return the error function of 'x'.  The algorithm is about twice as
     fast as the Bürmann approximation and is about an order of magnitude
     more accurate:
       Absolute error:
         max(abs(algebra.ferf-math.erf)) == 0.000357842621864 at
         abs(x) == 0.8783  and   abs(x) == 1.83295
       Relative error:
         max(abs( (algebra.ferf-math.erf)/math.erf )) == .000466 at
         abs(x) == 0.8"""
  return c_ferf(float(x))

def fcdf(double x):
  """Return the cummulative distribution function of 'x', utilizing
     the C implementation 'c_ferf' of 'ferf'."""
  return c_fcdf(float(x))

def factorial(int n):
  """Compute integer 'n' factorial.  Return 'None' for negative argument."""
  f=c_factorial(int(n))
  return f if f>0 else None

def permutation(int n,int k):
  """Choose 'k' elements from a set of 'n', n>=k, order matters.
     Return integer or 'None' for negative arguments or k>n."""
  p=c_permutation(int(n),int(k))
  return p if p>0 else None

def combination(int n,int k):
  """Choose 'k' elements from a set of 'n', n>=k, order is unimportant.
     Return integer or 'None' for negative arguments or k>n."""
  c=c_combination(int(n),int(k))
  return c if c>0 else None

def zero(double x,int precision=15):
  """Return 0.0 when 0.0 <= abs(x) <= 10**(-precision) else return 'x'."""
  return c_zero(x,precision)

def even(int i):
  """True when integer 'i' is even."""
  return i&1==0

def odd(int i):
  """True when integer 'i' is odd."""
  return i&1==1

def sign(double x,int precision=15):
  """Return: -1.0 => x < 0; 0.0 => x==0; 1.0 => x > 0.  'sign' will
     return 0.0 fox 'x' within 'precision' decimal places of zero."""
  if c_zero(x,precision)==0.0: return 0.0
  return copysign(1.0,x)

def clean(double x, int n=15):
  """Round 'x' to 'n' significant decimal digits, e.g.
       clean( 0.012345, 2) ==       0.012
       clean( 0.012345, 4) ==       0.01235
       clean(-0.012345, 4) ==      -0.01235
       clean(1234500.0, 4) == 1235000.0"""
  cdef double y
  y=fabs(x)
  if y < LOG_ARG_MIN: return 0.0
  return round( x, int( n-ceil(log10(y)) ) )

def len(object obj not None):
  """Returns 'obj.__len__()', omitting sanitizing checks imposed by
     '__builtins__.len(obj)'.  Unlike 'len', return values may be negative,
     floating point, or any arbitrary object in cases where such values
     define a meaningful 'length', accomodating specifically geometric cases."""
  return obj.__len__()

def delta_is_zero(object a not None, object b not None, int precision=15):
  """Return 'round(abs(a-b),precision)==0.0'.  Use when values are <= unity
     i.e. floating point noise is already beyond precision decimal places.
     Works with any type defining 'abs' that returns a 'float'."""
  return c_zero(abs(a-b),precision)==0.0

def delta_ratio_is_zero(object a not None, object b not None, int precision=15):
  """Return 'round(abs(a-b)/max(abs(a),abs(b)),precision)==0.0'.  Use when
     values may be very large: the division pushes floating point noise
     beyond precision decimal places.  Works with any type defining 'abs'
     that returns a 'float'."""
  return c_zero(abs(a-b)/max(abs(a),abs(b)),precision)==0.0

## excised after testing:
"""
from time import time
def erf_exec_time_compare(x=1.0,iterations=1000000):
  t=time(); c_execute_erf(x,iterations); t=time()-t
  s=time(); c_execute_ferf(x,iterations); s=time()-s
  return t/s
# This code shows the cython ferf can be as much as 3x as fast as
# the libm implementation of erf, but is typically about 2.2x.  The libm
# implementation, unlike python's math.erf uses piecewise approximation,
# typically calculating two eight degree polynomials and an exponential,
# plus a few floating point additions and multiplications.  It is about
# 3x as fast as math.erf; thus cython still introduces considerable
# overhead that a pure C extension module would eliminate.  Optimally,
# ferf should be about 20x math.erf and 6x the libm erf.
#
# In typical use, however, called via the python wrappers, this code is
# 3x faster than the math.erf implementation.  Thus a pure C extension
# module will be necessary.
#"""
