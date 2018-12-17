#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Package: algebra
Module : calculus.py
Website: https://github.com/ruminations/site-packages
License: https://github.com/ruminations/Licenses#design-license
Initial Copyright 2018

Defines a set of classes whose '__call__' and '__getitem__' interface
return scalar functions of derivatives, compatible with the interface
expected by 'geometry.parametric.Path'.

The '__call__' and '__getitem__' methods return anonymous functions of a
scalar variable.  A 'geometry.parametric.Path' may then be constructed,
representing a 'Function' by a trivial parameterization as a path in the
plane for calculating arc length parameterized differential quantities.

Tools are provided to compose these classes, and when such a composition
is passed to 'geometry.parametric.Path', the result is a powerful general
tool for calculating differential geometric quantities.

The '__repr__' and '__str__' methods follow the convention that 'repr(f)'
for 'f' a 'Function' object instance returns a string that is valid Python
source text, equivalent to the defining source.  Calling 'str(f)', in
contrast yields an ad-hoc string representation, more or less consistent
with typical ascii text expressions of mathematics.  If this text were
parsed and used to reproduce 'f', the result would not be exact, for
floating point quantities are rounded to four decimal places.  Moreover,
the syntax is not valid python text.
"""
import math # standard 'math' functions are prefixed for clarity: e.g. math.sqrt
from math import copysign
from operator import __mul__ as mul
from algebra.mathx import factorial,permutation,combination,clean
from geometry.gauss import C,Angle
from data.catalog import Physics as phy
__all__=['ConvergenceError','Domain','Function','Monomial','Polynomial', \
         'Sum','Partition','Exp','Log','Trig','HypTrig','Integral','Bruno', \
         'kld','Bose','Maxwell','Normal','LogNormal','BlackBody','BBpower']
__version__=201801215 # clean up doc strings
__version_history__=(20180326,20180403,20180424,20180430,20180504)
## TODO:
# 0) as always, writing a doctest suite.
#---------------+---------------+---------------+---------------#
#---------------+--- Implementation Utilities --+---------------#
#---------------+---------------+---------------+---------------#

class ConvergenceError(Exception):
  """Use to record data representing the final state of a numerical
     calculation that fails to converge.  This can be used for
     debugging or programmatically restarting a calculation with
     better initial values based on the failed values."""
  def __init__(self,value,intermediates,iteration):
    """Raise the error with this data - e.g for Newton's method:
         value         : x - the estimate of the root in the current cycle
         intermediates : (f(x),f_prime(x)) - the function & derivative
                         at x whose roots are sought
         iteration     : the cycle that failed
       Other calculations may provide data in any suitable way, including
       setting the values to arbitrary objects."""
    self._s = value,intermediates,iteration
    super(ConvergenceError,self).__init__(str(self))

  def __repr__(self): return "ConvergenceError({!r},{!r},{!r})".format(*self._s)

  def __str__(self):
    return "value={}, intermediates={},  iteration={}".format(*self._s)

  @property
  def intermediates(self):
    """Intermediate values of the numerical calculation when aborted."""
    return self._s[1]

  @property
  def iteration(self):
    """The numerical calculation iteration when aborted."""
    return self._s[2]

  @property
  def value(self):
    """The numerical calculation approximation when aborted."""
    return self._s[0]

class Bruno(object):
  """Compute the n-th derivative of the composition of two functions
     'f' and 'g', using Faà di Bruno's algorithm via the Bell polynomial
     of two vectors 'u' and 'v'.

     Calling an instance 'b=Bruno(f,g,n)' on a value 'x', i.e. 'b(x)',
     yields the value of:  d^n(f(g(x))) / dx^n

     The algorithm is based on the mathematics described in:
        http://dx.doi.org/10.1080/23311835.2016.1220670

     If 'f' and 'g' are instances of a subclass of 'Function' defined in
     this module, then 'f[i]' and 'g[i]' yield functions calculating the
     i-th respective derivative.  Invoking a 'Bruno' instance evaluates
     the Bell polynomial for the two vectors 'u=[h(x) for h in f[1:n+1]]'
     and 'v=[h(x) for h in g[1:n+1]]', returning a floating point value.

     The Bell polynomial is thus a combinatorial inner product, useful for
     computing the n-th order derivative chain rule, and other quantities.

     This implementation is purely recursive and does not compute partial
     Bell polynomials.  Translating indices in the reference to Python
     was tricky enough for this simpler algorithm.  Code readers beware,
     this code can look wrong, but is actually correct.  If you are getting
     incorrect results, there is probably an error in your derivatives."""
  __slots__=['_f','_g','_n']
  def __init__(self,f,g,n):
    """'f' and 'g' are iterables over functions and their derivatives; 'n'
       is the order of derivative of f(g(x)) evaluated by 'Bruno(f,g,n)(x)'."""
    self._f,self._g,self._n = f,g,n
  def __call__(self,x):
    """Evaluate the n-th derivative of f(g(x)) at 'x'."""
    gx,l = self._g[0](x),self._n+1
    return self._y([h(gx) for h in self._f[1:l]],[h(x) for h in self._g[1:l]])
  def _y(self,u,v):
    """The recursive Bell algorithm for Yn, where 'n==len(u)==len(v)>=1'"""
    # u[0] defines default Y0 at this recursion level
    # separating out of sum avoids null lists in '_y' recursions
    # 'sum( ... range(0)) == 0', so recursion terminates when 'n == 0'.
    n=len(u)-1; return u[0]*v[n]+ \
      sum(combination(n,k)*self._y(u[1:n-k+1],v[0:n-k])*v[k] for k in range(n))

def kld(p,q,d):
  """Kullback-Leibler divergence of distribution 'q' from 'p' \
     using values in Domain object 'd' to approximate the integral.

     Loosely, this is a measure of 'surprise' when using 'q' in a context
     where 'p' is expected.  When 'q' is 'p', the divergence is zero.
     Larger divergences imply larger surprises.  If 'p' represents
     empirical data, and 'q' is an analytic function, a smaller divergence
     implies 'q' is a better approximation of 'p' than some other function
     that produces a larger divergence.

     More generally, 'kld' serves as a measure of the global fidelity of
     a function 'q' used to approximate a function 'p'.  In those cases,
     'p' and 'q' are not normalized to unity, and 'kld' may be negative.

     That negative value implies 'q', on average, over estimates 'p' where
     a positive value implies an average under estimate.  Larger absolute
     magnitude implies a worse approximation, but in this generalized
     application, small magnitude does not necessarily imply a good
     approximation, it just implies a balance of errors."""
  s,i,it,lim = 0.0,0,iter(d),d.count
  v0=next(it); p0,q0 = p(v0),q(v0); f0=p0*math.log(p0/q0)
  while True: # trapezoidal integration of p*log(p/q) over domain
    if i>=lim: break
    v1=next(it); p1,q1 = p(v1),q(v1); f1=p1*math.log(p1/q1)
    s += (v1-v0)*(f0+f1)/2. # accomodate geometric intervals
    v0,f0 = v1,f1; i+=1
  return s

#---------------+---------------+---------------+---------------#
#---------------+--------- Base Tools ----------+---------------#
#---------------+---------------+---------------+---------------#

class Domain(object):
  """Defines a scalar domain as an iterable over which a function or a
     parameterized path may be evaluated."""
  __slots__=['_s','_e','_c','_g','_l','_h',
            '_i','_v','_d','_si','_ei'] # iterator vars this line

  def __init__(self, start, end, count=32,geometric=False):
    """'start' is a scalar at the beginning of the interval defined by
       this 'Domain'; 'end' is a scalar at the end of the interval;
       'count' specifies how many sub-intervals the interval is to be
       divided into.

       If start>end, the iterator is implicitly reversed, i.e. return
       values proceed from 'end' to 'start', i.e. return values ALWAYS
       proceed from more negative to more positive.

       Sub-intervals are arithmetic by default, with size calculated as
       dx=(end-start)/float(count).  If 'geometric==True', sub-intervals
       are multiplicative with ratio dx=abs(end/start)**(1/float(count)).
       A 'ValueError' is raised for geometric intervals that include or
       span the origin.

       Normal geometric intervals increase as the iterator proceeds from
       start to end; intervals decrease when abs(start)>abs(end), but
       produce an identical set of values, either due to reversal or to
       negative bounds.  Arithmetic return sets are also identical, but
       differ from geometric values by definition."""
    self._s,self._e = float(start),float(end)
    self._c=round(abs(count)); self._g=bool(geometric)
    self._i=self._d=self._v=0 # iterator variables
    self._rectify_bounds() # sets self._l and self._h

  def __contains__(self,v):
    """Return 'True' if 'v' is in the closed interval defining this
       'Domain'.  'v' need not be a return value of the iterable."""
    return v>=self._l and v<=self._h

  def __iter__(self):
    self._i=-1; l,c = abs(self.length),float(self.count)
    s,e = self._s,self._e
    if self._g: self._v,self._d = 1.0,(e/s)**(1.0/c) # v is cumulative
    else: self._v,self._d = 0.0,l/c # d is difference or ratio of v(i),v(i+1)
    if s>e: self._si,self._ei = self._e,self._s # reverse return
    else:   self._si,self._ei = self._s,self._e # sequence
    return self

  def __len__(self):
    """Return number of intervals in this domain."""
    return self.count

  def next(self):
    self._i+=1
    if self._i>self._c:  raise StopIteration
    if self._i==0:       return self._si # guard against
    if self._i==self._c: return self._ei # cumulative error
    if self._g: self._v*=self._d; return self._si*self._v
    else:       self._v+=self._d; return self._si+self._v

  def __repr__(self): return "Domain({},{},count={},geometric={})".format( \
                                self._s,self._e,self._c,self._g)

  def _rectify_bounds(self):
    """Private method to check for valid geometric bounds and rectify
       data after property setting."""
    s,e = self._s,self._e
    if s>e: self._l,self._h = e,s # for containment
    else:   self._l,self._h = s,e # and bounds test
    if not self._g: return
    if e==0.0 or s==0.0 or (copysign(1.0,e)!=copysign(1.0,s)): raise ValueError, \
      "Geometric interval includes origin: start == {} and end == {}".format(s,e)

  @property
  def delta(self):
    """Get the arithmetic differential interval length or
       the geometric interval ratio."""
    if self._g: return (self._e/self._s)**(1.0/self._c)
    return abs(self.length)/float(self._c)

  @property
  def count(self):
    """Get or set the number of intervals in the 'Domain'."""
    return self._c
  @count.setter
  def count(self,integer): self._c=int(abs(integer))

  @property
  def start(self):
    """Get or set the start point of the 'Domain'."""
    return self._s
  @start.setter
  def start(self,scalar): self._s=float(scalar); self._rectify_bounds()

  @property
  def end(self):
    """Get or set the end point of the 'Domain'."""
    return self._e
  @end.setter
  def end(self,scalar): self._e=float(scalar); self._rectify_bounds()

  @property
  def length(self):
    """Get the difference 'end - start' of the 'Domain'.  This may be
       negative if 'start' is greater than 'end'."""
    return self._e-self._s

  @property
  def span(self):
    """Get a 'C(l,h)' where l < h are the terminals of the 'Domain'."""
    return C(self._l,self._h)

  @property
  def geometric(self):
    """Get or set a boolean enabling geometric progression across the
       'Domain'.  The default is 'False'.  When 'True', the intervals
       intervals are more closely spaced nearer the origin, increasing
       geometrically to the final return value."""
    return self._g
  @geometric.setter
  def geometric(self,boolean): self._g=bool(boolean); self._rectify_bounds()


class Function(object):
  """An infinite iterator class for objects that define and provide
     access to a sequence of anonymous functions, starting with a base
     function and subsequently yielding higher order derivatives.

     'Function' object instances may be composed in an intuitive way via
     the arithmetic operator overrides defined by this class.  For example,
         'normal=Exp(1.0/sqrt(pi),-1.0)( Monomial(2.0) )'
     Defines the standard normal distribution with mu=0.0 sigma=1.0 and
         'maxwell=Monomial(2.0)*normal'
     defines an unnormalized related Maxwellian distribution."""
  __slots__=['_i','_j','_k']

  def __init__(self):
    self._i,self._j,self._k = 0,None,1 # iterator slice indices

  def __call__(self,x,derivative=0):
    """Evaluate the 'derivative' at 'x', where the zero derivative is the
       defining function.  If 'isinstance(x,Function)==True' return a
       'Function' instance representing the composition self(x(<parameter>))."""
    if isinstance(x,Function): return _Chain(self,x)
    return self[derivative](x)

  def __add__(self,g): return Sum(self,g)

  def __radd__(self,g): return Sum(g,self)

  def __div__(self,g): return _Product(self,~g)

  def __rdiv__(self,g): return _Product(g,~self)

  def __getitem__(self,n):
    """Override this method to return the n-th derivative of a user
       'Function' definition sub-class.  The default defined here is
       'None' which raises a glaring error if you forget to override."""
    return None

  def __getslice__(self,*s):
    """Get a finite slice of the iterator as a tuple, or get an infinite
       iterator:
         function[:3]    yields the function and the first two derivatives.
         function[2:4]   yields the second and third derivatives.
         function[2:5:2] yields the second and fourth derivatives.
         function[1::2]  yields an iterator over all odd derivatives.
         function[1:]    yields an iterator over all derivatives.
         function[:]     yields an iterator over all derivatives, starting
                         with the function or zero order derivative."""
    s=slice(*s)
    self._i = 0 if s.start is None else s.start
    self._j = s.stop # 'None' if omitted => infinite
    self._k = 1 if s.step  is None else s.step
    if self._j is None: return iter(self)
    return tuple(self[i] for i in range(self._i,self._j,self._k))

  def __invert__(self): return _Chain(Monomial(-1.0),self)

  def __iter__(self):
    self._i,self._j,self._k = 0,None,1
    return self

  def __mul__(self,g):
    """Return a 'Function' instance representing self(x)*g(x)."""
    return _Product(self,g)

  def __rmul__(self,g): return _Product(g,self)

  def __neg__(self): return _Negate(self)

  def __or__(self,g):
    """Return an 'Integral' instance representing the numerical inner
       product integral of self(x)*g(x)."""
    return Integral(self,g)

  def __sub__(self,g): return Sum(self,-g)

  def __rsub__(self,g): return Sum(g,-self)

  def next(self):
    """Get the next derivative in this iterator's sequence."""
    if self._i==self._j: raise StopIteration
    f=self[self._i]; self._i+=self._k; return f

  def _nm(self,f,fp,g=1.0):
    """Analytic Newton's method."""
    delta,epsilon = 1.0,.0000000000000003 # 15.7 decimal places
    i,limit = 0,20
    while abs(delta)>epsilon:
      try:delta=f(g)/fp(g)
      except ZeroDivisionError: raise ConvergenceError(g,(f(g),fp(g),delta),i)
      g-=delta; i+=1
      if i==limit: raise ConvergenceError(g,(f(g),fp(g),delta),i)
    return g

  def extremum(self,guess):
    """Return a value near 'guess' where the first derivative of the
       function evaluates to zero."""
    fp,f2p = self[1:3]; return self._nm(fp,f2p,guess)

  def inflection(self,guess):
    """Return a value near 'guess' where the second derivative of the
       function evaluates to zero."""
    f2p,f3p = self[2:4]; return self._nm(f2p,f3p,guess)

  def root(self,guess):
    """Return a value near 'guess' where the function evaluates to zero."""
    f,fp = self[0:2]; return self._nm(f,fp,guess)

class _Chain(Function):
  """Private class for implementing 'Function.__call__'."""
  __slots__=['_f','_g']
  def __init__(self,f,g):
    """Create a function implementing the composition of 'Function' objects
       'f' and 'g': p(x)=f(g(x)) and all of its derivatives."""
    if not all(isinstance(h,Function) for h in (f,g)): raise ValueError, \
      "{!r} and {!r} must be 'Function' object instances".format(f,g)
    super(_Chain,self).__init__(); self._f,self._g = f,g
  def __getitem__(self,n):
    """Return a static function that evaluates the n-th derivative of
       this 'Function' at a parameter 'x' when called on a value 'x'.

       This is implemented through the recursive function 'Bruno(f,g,n)'."""
    if n==0: return lambda x: self._f[0](self._g[0](x))
    return Bruno(self._f,self._g,n)
  def __repr__(self): return "{!r}({!r})".format(self._f,self._g)
  def __str__(self):  return str(self._f).replace('x',"({!s})".format(self._g))

class _Product(Function):
  """Private class for implementing 'Function.__mul__'."""
  __slots__=['_f','_g']
  def __init__(self,f,g):
    """Create a function implementing the product of 'Function' objects
       'f' and 'g': p(x)=f(x)*g(x) and all of its derivatives."""
    if not all(isinstance(h,Function) for h in (f,g)): raise ValueError, \
      "{!r} and {!g} must be 'Function' object instances".format(f,g)
    super(_Product,self).__init__(); self._f,self._g = f,g
  def __getitem__(self,n):
    """Return a static function that evaluates the n-th derivative of
       this 'Function' at a parameter 'x' when called on a value 'x'."""
    return lambda x: sum(combination(n,k)*self._f[n-k](x)*self._g[k](x) \
                          for k in range(0,n+1))
  def __repr__(self): return "{!r}*{!r}".format(self._f,self._g)
  def __str__(self):  return "({!s})*({!s})".format(self._f,self._g)

class _Negate(Function):
  """Private class for implementing 'Function.__neg__'."""
  __slots__=['_f']
  def __init__(self,f):
    """Create a function implementing the negation of the 'Function'
       object 'f': r(x)=-f(x) and all of its derivatives."""
    if not isinstance(f,Function): raise ValueError, \
      "{!r} must be a 'Function' object instance".format(f)
    super(_Negate,self).__init__(); self._f=f
  def __getitem__(self,n):
    """Return a static function that evaluates the n-th derivative of
       this 'Function' at a parameter 'x' when called on a value 'x'."""
    return lambda x: -self._f[n](x)
  def __repr__(self): return "-{!r}".format(self._f)
  def __str__(self): return "-({!s})".format(self._f)

class Partition(Function):
  """Class for implementing 'Function.__getslice__[:]'.
     May be used independently."""
  __slots__=['_f','_n']
  def __init__(self,function,n=1):
    """Create a function implementing a partition of 'Function' object
       'function' beginning with its 'n'-th derivative.  The default
       base function is the first derivative of 'function'."""
    if not isinstance(f,Function): raise ValueError, \
      "{!r} must be a 'Function' object instance".format(f)
    super(Partition,self).__init__(); self._f,self._n = function,n
  def __getitem__(self,n):
    """Return a static function that evaluates the n-th derivative of
       this 'Partition' at a parameter 'x' when called on a value 'x'."""
    return self._f[n+self._n]
  def __repr__(self): return "Partition({!r},{!r})".format(self._f,self._n)
  def __str__(self):  return "D{}({!s})".format(self._n,self._f)

class Sum(Function):
  """Class for summing an iterable of 'Function' object instances."""
  __slots__=['_fs']
  def __init__(self,*functions):
    """Create a function implementing the sum of 'Function' objects
       contained in 'functions': s(x)=f0(x)+f1(x)+...
       and all of its derivatives of 's'."""
    super(Sum,self).__init__(); self._fs=tuple(functions)
    if not all(isinstance(f,Function) for f in self._fs): raise ValueError, \
      "All arguments must be 'Function' object instances"
  def __getitem__(self,n):
    """Return a static function that evaluates the n-th derivative of
       this 'Sum' at a parameter 'x' when called on a value 'x'."""
    return lambda x: sum(f[n](x) for f in self._fs)
  def __repr__(self): return '+'.join(  "{!r}".format(f) for f in self._fs)
  def __str__(self):  return '+'.join("({!s})".format(f) for f in self._fs)

#---------------+---------------+---------------+---------------#
#---------------+---- Function sub-classes -----+---------------#
#---------------+---------------+---------------+---------------#

class Monomial(Function):
  """The most basic scalar function."""
  __slots__=['_c','_e']

  class _f(object):
    """Anonymous monomial."""
    __slots__=['_c','_e','_g']
    def __init__(self,c,e):
      self._c,self._e = c,e
      if round(e)==e: # improve accuracy for integral exponents
        e=self._e=int(e)
        if 0> e: self._g=lambda x: self._c/reduce(mul,(x,)*-self._e)
        else   : self._g=lambda x: self._c*reduce(mul,(x,)* self._e,1.0)
      else     : self._g=lambda x: self._c*(x**self._e)
    def __call__(self,x): return self._g(x)

  def __init__(self,exponent=1.0,coefficient=1.0):
    """'exponent' and 'coefficient' define the 'Monomial':
         f(x) = coefficient*x^(exponent) and all of its derivatives.
       The exponent need neither be positive nor integral - fractional
       exponents are interpreted correctly facilitating root functions.
       The default is the identity monomial f(x)=x"""
    super(Monomial,self).__init__(); self._c,self._e = coefficient,exponent
  def __getitem__(self,n):
    """Return a static function that evaluates the n-th derivative of
       this 'Monomial' at a parameter 'x' when called on a value 'x'."""
    if n==0: return lambda x: self._c*(x**self._e)
    e=self._e-n; c=self._c*reduce(mul,[self._e-i for i in range(n)])
    if c==0.0: return lambda x: 0.0 # integral exponent, n>self._e
    return self._f(c,e)
  def __repr__(self): return "Monomial({!r},{!r})".format(self._e,self._c)
  def __str__(self):
    s=chr(ord('+')+(int(self._c<0.0)<<1))
    return "{}{:.4f}*x^{:.4f}".format(s,abs(self._c),self._e)

class Polynomial(Function):
  """The second most basic scalar function."""
  __slots__=['_c']

  class _f(object):
    """Anonymous polynomial."""
    __slots__=['_c']
    def __init__(self,c): self._c=tuple(c)
    def __call__(self,x):
      r=0.0
      for c in self._c: r=x*r+c
      return r

  def __init__(self,coefficients=(1.0,1.0)):
    """'coefficients' is an iterable yielding the constant coefficient
       for each monomial term that defines the 'Polynomial'.  For
       example, Polynomial([3.0,2.0,1.0]) defines the quadratic
       f(x) = 3.0*x^2 + 2.0*x + 1.0  and all of its derivatives.
       The default is the linear binomial f(x)=x+1.0"""
    super(Polynomial,self).__init__(); self._c=tuple(coefficients)
  def __getitem__(self,n):
    """Return a static function that evaluates the n-th derivative of
       this 'Polynomial' at a parameter 'x' when called on a value 'x'."""
    l=len(self._c)
    if n>=l: return lambda x: 0.0
    perm=(permutation((l-1)-i,n) for i in range(l-n))
    return self._f(c*p for c,p in zip(self._c[:l-n],perm) )
  def __len__(self):
    """Given 'p' a 'Polynomial' its degree is 'len(p) == len(coefficients)-1'."""
    return len(self._c)-1
  def __repr__(self): return "Polynomial({!r})".format(self._c)
  def __str__(self):
    l=[(chr(ord('+')+(int(c<0.0)<<1)),abs(c),len(self._c)-1-i) \
        for i,c in enumerate(self._c)]
    return ''.join(["{}{:.4f}*x^{}".format(*t) for t in l[:-1]]) +\
                    "{}{:.4f}".format(*l[-1][:-1]) # constant term

class Exp(Function):
  """The third most basic scalar function - the exponential to the base 'e'."""
  __slots__=['_a','_b','_c']
  def __init__(self,a=1.0,b=1.0,c=0.0):
    """'a' and 'b' are convenience parameters so that common exponential
       forms  do not require composition to define.  This specifies
       the function: f(x) = b*exp(a*x)+c and all of its derivatives.
       The default is the simple exponential f(x)=e^x"""
    super(Exp,self).__init__(); self._a,self._b,self._c = map(float,(a,b,c))
  def __getitem__(self,n):
    """Return a static function that evaluates the n-th derivative of
       this 'Exp' at a parameter 'x' when called on a value 'x'."""
    if n==0: return lambda x:self._b*math.exp(self._a*x)+self._c
    return lambda x: self._b*(self._a**n)*math.exp(self._a*x)
  def __repr__(self): return "Exp({!r},{!r},{!r})".format(self._a,self._b,self._c)
  def __str__(self):
    s='' if (self._c==0.0) else "{}{:.4f}".format( \
       chr(ord('+')+(int(self._c.rad<0.0)<<1)),abs(self._c) )
    return "{:.4f}*e^({:.4f}*x){}".format(self._b,self._a,s)

class Log(Function):
  """The fourth most basic scalar function - the logarithm to the base 'e'."""
  __slots__=['_a','_b','_c','_m']
  def __init__(self,a=1.0,b=1.0,c=0.0):
    """'a', 'b', and 'c' are convenience parameters so that common logarithmic
       forms  do not require composition to define.  This specifies the function:
       f(x) = a*log(b*x+c) and all of its derivatives.
       The default is the simple logarithm f(x)=log(x)"""
    super(Log,self).__init__(); self._a,self._b,self._c = map(float,(a,b,c))
    self._m=lambda x: self._b*x+self._c
  def __getitem__(self,n):
    """Return a static function that evaluates the n-th derivative of
       this 'Log' at a parameter 'x' when called on a value 'x'."""
    if n==0: return lambda x:self._a*math.log(self._m(x))
    return lambda x: -self._a*factorial(n-1)*((-self._b/self._m(x))**n)
  def __repr__(self): return "log({!r},{!r},{!r})".format(self._a,self._b,self._c)
  def __str__(self):
    s='' if (self._c==0.0) else "{}{:.4f}".format( \
       chr(ord('+')+(int(self._c.rad<0.0)<<1)),abs(self._c) )
    return "{:.4f}*log({:.4f}*x{})".format(self._a,self._b,s)

class Trig(Function):
  """The elementary trigonometic functions - cosine and sine."""
  __slots__=['_a','_b','_c','_f','_fp','_s']
  def __init__(self,a=1.0,b=1.0,c=Angle(0.0),sine=False):
    """'a', 'b', and 'c' are convenience parameters so that linear forms
       do not require composition to define.  This specifies the function:
       f(x) = a*cos(bx+c) and all of its derivatives.

       The default is the simple cosine f(x)=cos(x).  Using 'sine==True'
       creates a 'Trig' instance using 'sin' instead of 'cos' using the
       same constants.  The same effect can be produced with the default
       calculation function 'cos' by adding Angle(pi/2.0,'rad') to the
       value of 'c' that is passed.  The value 'c' may be a floating point
       number.  In that case, it is interpreted in units of radians."""
    super(Trig,self).__init__()
    self._a,self._b = map(float,(a,b)); self._c,self._s=Angle(c,'rad'),bool(sine)
    if sine: self._f,self._fp = math.sin,math.cos
    else:    self._f,self._fp = math.cos,math.sin
  def __getitem__(self,n):
    """Return a static function that evaluates the n-th derivative of
       this 'Trig' at a parameter 'x' when called on a value 'x'."""
    if self._s:
      neg,f = (1.0,1.0,-1.0,-1.0)[n%4]*self._a,(self._f,self._fp)[n%2]
    else:
      neg,f = (1.0,-1.0,-1.0,1.0)[n%4]*self._a,(self._f,self._fp)[n%2]
    return lambda x: neg*(self._b**n)*f(self._b*x+self._c)
  def __repr__(self): return "Trig({!r},{!r},{!r})".format(self._a,self._b,self._c)
  def __str__(self):
    s='' if (self._c.hemi==0.0) else "{}{:.4f}".format( \
       chr(ord('+')+(int(self._c.rad<0.0)<<1)),abs(self._c.rad) )
    if self._s: return "{:.4f}*sin({:.4f}*x{})".format(self._a,self._b,s)
    return "{:.4f}*cos({:.4f}*x{})".format(self._a,self._b,s)

class HypTrig(Function):
  """The elementary hyperbolic trigonometic functions - cosh and sinh."""
  __slots__=['_a','_b','_c','_f','_fp','_s']
  def __init__(self,a=1.0,b=1.0,c=0.0,sinh=False):
    """'a', 'b', and 'c' are convenience parameters so that linear forms
       do not require composition to define.  This specifies the function:
       f(x) = a*cosh(bx+c) and all of its derivatives.

       The default is the simple hyperbolic cosine f(x)=cos(x).  Using
       'sinh==True' creates a 'HypTrig' instance using 'sinh' instead of
       'cosh' using the same constants."""
    super(HypTrig,self).__init__()
    self._a,self._b,self._c = map(float,(a,b,c)); self._s=bool(sinh)
    if sinh: self._f,self._fp = math.sinh,math.cosh
    else:    self._f,self._fp = math.cosh,math.sinh
  def __getitem__(self,n):
    """Return a static function that evaluates the n-th derivative of
       this 'HypTrig' at a parameter 'x' when called on a value 'x'."""
    f=(self._f,self._fp)[n%2]
    return lambda x: (self._b**n)*f(self._b*x+self._c)
  def __repr__(self): return "Trig({!r},{!r},{!r})".format(self._a,self._b,self._c)
  def __str__(self):
    s='' if (self._c==0.0) else "{}{:.4f}".format( \
       chr(ord('+')+(int(self._c<0.0)<<1)),abs(self._c) )
    if self._s: return "{:.4f}*sinh({:.4f}*x{})".format(self._a,self._b,s)
    return "{:.4f}*cosh({:.4f}*x{})".format(self._a,self._b,s)

class Integral(Function):
  """A numerical integral of a 'Function' object instance over a 'Domain'
     instance.  Analytical derivatives are represented by the integrand."""
  __slots__=['_f','_l','_u','_dx','_dy','_lim','_g']
  def __init__(self,function):
    """'function' is the integrand.  By default, the function independent
       variable is the integral's upper limit and the lower limit is zero.
       These and other behaviors may be changed via properties."""
    super(Integral,self).__init__(); self._f=function
    self._l,self._u,self._dx,     self._dy,self._lim,self._g = \
        0.0,   None,      .1,.000000000001,     1000,  False
  def __call__(self,x,derivative=0):
    """Evaluate the 'derivative' at 'x', where the zero derivative is
       defined by the numerical integral of the defining integrand.  If
       'isinstance(x,Function)==True' return a 'Function' instance
       representing the composition self(x(<parameter>)).

       'x==None' indicates an infinite integral.  This calculation will
       terminate for integrands known to quickly asymptotically reach
       zero: the integration stops when the summand becomes less than
       the tolerance set via properties."""
    if isinstance(x,Function): return _Chain(self,x)
    if derivative>0: return self._f[derivative-1](x)
    return self._trapezoid(x)
  def __getitem__(self,n):
    """Return a static function that evaluates the n-th derivative of
       this 'Integral' at a parameter 'x' when called on a value 'x'."""
    if n==0: return lambda x: self._trapezoid
    return self._f[n-1]
  def __repr__(self): return "Integral({!r})".format(self._f)
  def _trapezoid(self,x):
    """Evaluate the integral subject to the current parameters using a
       numerical trapezoidal approximation."""
    f,l,u,dx,g = self._f,self._l,self._u,self._dx,self._g
    if l is None and u is None: raise ValueError, \
      "Integral.lower and Integral.upper can not both be 'None'"
    c = (lambda a,b: round(abs(math.log(b/a)/math.log(dx+1.0)))) if g else \
        (lambda a,b: round(abs((b-a)/dx)))
    if x is not None: # finite integral
      if l is None: d,s = Domain(x,u,c(x,u),g),copysign(1.0,u-x)
      if u is None: d,s = Domain(l,x,c(l,x),g),copysign(1.0,x-l)
      else:         d,s = Domain(l,u,c(l,u),g),copysign(1.0,u-l)
      xs=list(d); dxs=[b-a for a,b in zip(xs[:-1],xs[1:]) ]
      ys=tuple(f(v) for v in xs)
      means=tuple((yl+yu)/2.0 for yl,yu in zip(ys[:-1],ys[1:]))
      return s*sum(m*dx for m,dx in zip(means,dxs))
    dy,lim = self._dy,self._lim # infinite integral
    if l is None:
      if u>=0.0: raise ValueError, \
        "Upper bound == {} on integral to -infinity must be negative.".format(u)
      if g: v,dx = u,(1.0+dx)
      else: v,dx = u,  -dx
    else:
      if l<0.0: raise ValueError, \
        "Lower bound == {} on integral to +infinity must be positive.".format(l)
      if g: v,dx = l,(1.0+dx)
      else: v,dx = l,   dx
    if g: yl,yh = f(v),f(v*dx)
    else: yl,yh = f(v),f(v+dx)
    sm=0.0; i=0
    while i<lim:
      du=dx*(yl+yh)/2.0; sm+=du
      if du<dy: return copysign(1.0,v)*sm
      v = v*dx if g else v+dx
      yl,yh = yh,f(v)
    raise ConvergenceError(sm,(v,du,yl,yh),i)

  @property
  def geometric(self):
    """Get or set a boolean enabling geometric progression across the
       integral's 'Domain'.  The default is 'False'.  When 'True', the
       intervals are more closely spaced nearer the origin, increasing
       geometrically to 'abs(upper)'."""
    return self._g
  @geometric.setter
  def geometric(self,boolean): self._g=bool(boolean)

  @property
  def lower(self):
    """Get or set the lower limit of the integral.  Setting to 'None'
       implies the lower limit is the integral's independent variable."""
    return self._l
  @lower.setter
  def lower(self,limit): self._l = limit if limit is None else float(limit)

  @property
  def upper(self):
    """Get or set the upper limit of the integral.  Setting to 'None'
       implies the upper limit is the integral's independent variable."""
    return self._u
  @upper.setter
  def upper(self,limit): self._u = limit if limit is None else float(limit)

  @property
  def delta(self):
    """Get or set the discrete differential of the integral.  The default
       is .1, but appropriate sizes are highly application dependent.  A
       nominal value of about 1.0/100.0 times the largest interval to be
       integrated may work in many cases.

       For geometric domains, 'delta' is the differential increase of the
       unit interval in the geometric progression, that is:
          delta = ratio - 1.0   or   ratio = delta + 1.0
       so the ratio is assumed to be greater than unity, leading to
       increasing intervals in progression away from the origin."""
    return self._dx
  @delta.setter
  def delta(self,dx): self._dx=abs(float(dx))

  @property
  def limit(self):
    """Get or set the iteration termination limit for infinite integrals.
       The default is 1000.  Larger values may be required for slowly
       converging integrals.  A 'ConvergenceError' is raised if the limit
       is exceeded when evaluating infinite integrals."""
    return self._lim
  @limit.setter
  def limit(self,iterations): self._lim=abs(int(iterations))

  @property
  def tolerance(self):
    """Get or set the summand termination tolerance for infinite integrals.
       The default is .000000000001 (12 decimal places).  Larger values
       may be required for slowly converging integrals.  Values smaller
       than 14 decimal places or larger than .1 raise a 'ValueError'."""
    return self._dy
  @tolerance.setter
  def tolerance(self,dy):
    self._dy=abs(float(dy))
    if self._dy<.00000000000001 or self._dy>.1: raise ValueError, \
      "Integration tolerance {} must be in range [1e-14 : .1]".format(dy)

#---------------+---------------+---------------+---------------#
#---------------+--- Distribution sub-classes --+---------------#
#---------------+---------------+---------------+---------------#

class Bose(Function):
  """An abstraction of the Bose-Einstein distribution that may be used
     for, among other puposes, expressing the black body radiation
     distribution by composing this distribution with a particular
     parameterization and degeneracy function."""
  __slots__=['_a','_b','_c','_f']
  def __init__(self,a=1.0,b=1.0,c=-1.0):
    """'a' and 'b' are convenience parameters so that common
       forms  do not require composition to define.  This specifies
       the function: f(x) = 1.0/(b*exp(ax)+c) and all of its derivatives.
       The default is the simple  f(x)=1.0/(e^x-1.0); if c is positive,
       the function is a Fermi-Dirac distribution."""
    super(Bose,self).__init__(); self._a,self._b,self._c = map(float,(a,b,c))
    self._f=~Exp(a,b,c)
  def __getitem__(self,n): return self._f[n]
  def __repr__(self): return "Bose({!r},{!r},{!r})".format(self._a,self._b,self._c)
  def __str__(self):  return "1.0/({:.4f}e^({:.4f}*x)-{:.4f})".format(self._b,self._a,self._c)

class Normal(Function):
  """The normal or Gaussian distribution function."""
  __slots__=['_s','_m','_f']
  def __init__(self,sigma=1.0,mu=0.0):
    """'sigma' is the standard deviation and 'mu' is the mean specifying:
          f(x) = exp(-((x-mu)^2)/(2*sigma^2)) / sqrt(2*pi*sigma^2)
       and all of its derivatives. The default is the standard distribution:
          f(x) = e^(-(x^2)/2)/sqrt(2*pi)"""
    super(Normal,self).__init__(); self._s,self._m = map(float,(sigma,mu))
    t=1.0/abs(self._s) # Stigler's formulation
    self._f=Exp(1.0,t/math.sqrt(2.0*math.pi))
    self._f=self._f(Monomial(2.0,-0.5)(Polynomial((t,-self._m*t))))
  def __getitem__(self,n): return self._f[n]
  def __repr__(self): return "Normal({!r},{!r})".format(self._s,self._m)
  def __str__(self):  return "{:.4f}e^(-((x-{:.4f})^2)/{:.4f})".format( \
    1.0/math.sqrt(2.0*self._v*math.pi),self._m,2.0*self._v)

class LogNormal(Function):
  """The log-normal function."""
  __slots__=['_s','_m','_v','_f']
  def __init__(self,sigma=1.0,mu=0.0):
    """'sigma' is the standard deviation and 'mu' is the mean specifying:
          f(x) = exp(-((log(x)-mu)^2)/(2*sigma^2)) / (x*sqrt(2*pi*sigma^2))
       and all of its derivatives. The default is the standard distribution:
          f(x) = e^(-.5*(log(x))^2) / (x*sqrt(2*pi))"""
    super(LogNormal,self).__init__(); self._s,self._m = map(float,(sigma,mu))
    self._f=Monomial(-1.0)*( Normal(self._s,self._m)(Log()) )
  def __getitem__(self,n): return self._f[n]
  def __repr__(self): return "LogNormal({!r},{!r})".format(self._s,self._m)
  def __str__(self):  return str(self._f)

class Maxwell(Function):
  """The Maxwellian distribution function."""
  __slots__=['_a','_f']
  def __init__(self,a=1.0):
    """'a' is typically sqrt(m/(kT)).  This specifies the function:
           f(x) = (a^3)*sqrt(2.0/pi)*(x^2)*exp(-.5*(ax)^2)
       and all of its derivatives. The default is the distribution:
           f(x) = sqrt(2/pi)(e^(-.5*x^2))."""
    super(Maxwell,self).__init__(); self._a=float(a); a2=self._a*self._a
    self._f=Monomial(2.0,self._a*a2*math.sqrt(2.0/math.pi))* \
            Exp()(Monomial(2.0,-.5*a2))
  def __getitem__(self,n): return self._f[n]
  def __repr__(self): return "Maxwell({!r})".format(self._a)
  def __str__(self):  return "{:.4f}*(x^2)*e^(-({:.4f}*x)^2)".format( \
       self._a*self._a*self._a*math.sqrt(2.0/math.pi),sqrt(.5)*self._a)

class BlackBody(Function):
  """Computes the black body radiation power distribution function
     at a fixed temperature according to various parameterizations:
         wavelength, frequency, energy

     The return value is in SI units of
       Joules/square-meter/steradian/<differential unit>"""
  __slots__=['_u','_t','_l','_f']
  def __init__(self,temperature=6504.0,unit='eV',log=False):
    """The  default 'temperature' corresponds to a CIE-D65 illuminant.
       The differential unit parameterization is given by 'unit' as follows:
          'unit'   meaning        parameterization
         ------------------------------------------
          'nm'    nano-meters       wavelength
          'Thz'   Tera-Hertz        frequency
          'eV'    electron-Volts    energy
        Setting 'log=True' composes the logarithm of the specified unit."""
    c,h,k,q = map(phy.value,('c','h','k','q'))
    self._u,self._t,self._l = str(unit),float(temperature),bool(log)
    kt=k*temperature; hc=h*c; giga,tera = 1.0e9,1.0e12; u=self._u.lower()
    if   u=='ev' :
         t       = q/hc
         a,b,e,i = 2.0*c*q*t*t*t,q/kt, 3.0, 1.0
    elif u=='nm' :
         s,t     = hc*giga,giga
         a,b,e,i = 2.0*c*s*t*t*t,s/kt,-5.0,-1.0
    elif u=='thz':
         s,t     = h*tera,tera/c
         a,b,e,i = 2.0*c*s*t*t*t,s/kt, 3.0, 1.0
    else: raise ValueError, "Unrecognized units: {!r}".format(unit)
    super(BlackBody,self).__init__()
    self._f=(Bose()(Monomial(i,b)))*Monomial(e,a)
    if self._l: self._f=self._f(Log())
  def __getitem__(self,n): return self._f[n]
  @property
  def unit(self):
    """Get string representing the parameter units."""
    return self._u
  @property
  def temperature(self):
    """Get numerical value of the black body temperature."""
    return self._t

class BBpower(object):
  """Implements a function via the '__call__' method that yields the
     total radiant power in kW/sqcm of a black body as a function of
     temperature over a >>fixed<< and >>finite<< range of integration."""
  __slots__=['_l','_u','_un','_bb','_d','_cnv']
  _conv_ev={'nm':1e9*phy.value('h')*phy.value('c')/phy.value('q'), \
           'thz':1e12*phy.value('h')/phy.value('q'), \
            'ev':1.0}
  def __init__(self,lower=890.0,upper=390.0,unit='nm',delta='eV'):
    """The default parameters are convenient for working with CIE-LMS data.
       'lower' and 'upper' are the limits of integration.  'unit' indicates
       the units in which these limits are expressed.  'delta' indicates
       the differential parameter to use for the power distribution function
       in the black body integrand.  Values are converted automatically to
       be compatible with the chosen integrand.  Unit values may be:
          'unit'   meaning        parameterization
         ------------------------------------------
          'nm'    nano-meters       wavelength
          'Thz'   Tera-Hertz        frequency
          'eV'    electron-Volts    energy"""
    if upper is None: # Assume eV units
      assert unit==delta=='eV', "Infinite integral must be expressed in 'eV'."
      self._l,self._u,self._un,self._bb,self._d = float(lower),upper,'eV','eV',1.0/200.0
      if self._l==0.0: self._l=self._d
    else:
      self._l,self._u,self._un,self._bb = float(lower),float(upper),str(unit),str(delta)
      # convert limit units to eV and then to distribution units
      self._cnv = self._conv_ev[self._un.lower()]/self._conv_ev[self._bb.lower()]
      if self._un.lower()=='nm': self._l,self._u = 1.0/self._l,1.0/self._u
      self._l,self._u = self._l*self._cnv,self._u*self._cnv
      if self._bb.lower()=='nm': self._u,self._l = clean(1.0/self._l),clean(1.0/self._u)
      self._d=(self._u-self._l)/1000.0

  def __call__(self,t):
    """'t' is the temperature in degrees Kelvin.  Return the total radiant
       power in kW/sqcm for this black body at this temperature over the
       specified spectral limits."""
    i=Integral(BlackBody(t,self._bb)); i.lower=self._l; i.delta=self._d
    return i(self._u)*math.pi*1.0e-7 # 1e-3 => kW; 1e-4 => sqcm from W/sqm

if __name__=="__main__":
  print "===Beginning module calculus.py test suite==="
  print "\nTest 'Polynomial' against explicit function definition."
  def g(x): return ((1.0*x+3.0)*x-2.0)*x-1.0
  def gp(x): return (1.0*3.0*x+3.0*2.0)*x-2.0
  f=Polynomial([1.0,3.0,-2.0,-1.0])
  print "\n Cubic and critical points:"
  print "    str: {!s}\n   repr: {!r}".format(f,f)
  for i in range(-35,15):
    x=i/10.
    print "    {:4.1f}  {:6.3f}  {:6.3f}  {:6.3f}  {:6.3f}  {:6.3f}  {:6.3f}".format( \
               x,g(x),f(x),gp(x),f(x,1),f(x,2),f(x,3))
  print "\n  Roots, extrema, and inflection:"
  print "    ",f.root(-3),f.root(-.5),f.root(1)
  print "    ",f.extremum(-2),f.extremum(.5)
  print "    ",f.inflection(0)
  print "\n\nExponential of logarithm (net y=x):"
  g=lambda x,d: round(Exp()(Log())(x,d),15)
  for i in range(1,20):
    x=i/10.
    print "    {:4.1f}  {:7.4f}  {:7.4f}  {:7.4f}  {:7.4f}".format( \
               x,g(x,0),g(x,1),g(x,2),g(x,3))
  print "\n\nLogarithm of exponential (net y=x):"
  g=lambda x,d: round(Log()(Exp())(x,d),15)
  for i in range(1,20):
    x=i/10.
    print "    {:4.1f}  {:7.4f}  {:7.4f}  {:7.4f}  {:7.4f}".format( \
               x,g(x,0),g(x,1),g(x,2),g(x,3))
  print "\n\nComposition of Monomial sqrt and Monomial square (net y=x):"
  g=Monomial(1.0/2.0)(Monomial(2.0))
  for i in range(1,15):
    x=i/10.
    print "    {:4.1f}  {:7.4f}  {:7.4f}  {:7.4f}  {:7.4f}".format( \
               x,g(x,0),g(x,1),g(x,2),g(x,3))
  print "\n\nComposition of Monomial cube and Monomial cbrt (net y=x):"
  g=Monomial(3.0)(Monomial(1.0/3.0))
  for i in range(1,15):
    x=i/10.
    print "    {:4.1f}  {:7.4f}  {:7.4f}  {:7.4f}  {:7.4f}".format( \
               x,g(x,0),g(x,1),g(x,2),g(x,3))
  print "\n\nCube root:"
  g=Monomial(1.0/3.0) # cube root
  for i in range(1,20):
    x=i/2.
    print "    {:4.1f}  {:7.4f}  {:7.4f}  {:7.4f}  {:7.4f}".format( \
               x,g(x,0),g(x,1),g(x,2),g(x,3))
  print "\n\n1st Partition of cube root:"
  g=Partition(g,1) # 1.0/(3.0*cbrt^2)
  for i in range(1,20):
    x=i/2.
    print "    {:4.1f}  {:7.4f}  {:7.4f}  {:7.4f}  {:7.4f}".format( \
               x,g(x,0),g(x,1),g(x,2),g(x,3))
  print "\n\nSquare root minus cube root:"
  g=Monomial(1.0/2.0)-Monomial(1.0/3.0)
  for i in range(1,10):
    x=i/1.
    print "    {:4.1f}  {:7.4f}  {:7.4f}  {:7.4f}  {:7.4f}".format( \
               x,g(x,0),g(x,1),g(x,2),g(x,3))
  print "\n\nSimple reciprocal:"
  g=Monomial(-1.0)
  for i in range(1,20):
    x=i/2.
    print "    {:4.1f}  {:7.4f}  {:7.4f}  {:7.4f}  {:7.4f}".format( \
               x,g(x,0),g(x,1),g(x,2),g(x,3))
  print "\n\nInverse square:"
  g=Monomial(-2.0)
  for i in range(1,20):
    x=i/2.
    print "    {:4.1f}  {:7.4f}  {:7.4f}  {:7.4f}  {:7.4f}".format( \
               x,g(x,0),g(x,1),g(x,2),g(x,3))
  print "\n\nIter test:"
  print "\n  Inverse square derivatives accessed by iterator:"
  g=Monomial(-2.0)
  for i in range(1,20):
    x=i/2.; v=[x]; it=iter(g)
    while len(v)<5: v+=[next(it)(x)]
    print "    {:4.1f}  {:7.4f}  {:7.4f}  {:7.4f}  {:7.4f}".format(*v)
  print "\n\nQuadratic as a product of binomials:"
  g=Polynomial((1.0,1.0))*Polynomial((1.0,1.0))
  for i in range(-15,15):
    x=i/10.
    print "    {:4.1f}  {:7.4f}  {:7.4f}  {:7.4f}  {:7.4f}".format( \
               x,g(x,0),g(x,1),g(x,2),g(x,3))
  print "\n\nQuadratic as a composition of monomial and binomial:"
  g=(Monomial(2.0))(Polynomial((1.0,1.0)))
  for i in range(-15,15):
    x=i/10.
    print "    {:4.1f}  {:7.4f}  {:7.4f}  {:7.4f}  {:7.4f}".format( \
               x,g(x,0),g(x,1),g(x,2),g(x,3))
  print "\n\nTrig cosine ; sine:"
  g,f = Trig(),Trig(sine=True)
  for i in range(-16,16):
    x=i/10.
    print "    {:4.1f}  {:6.3f}  {:6.3f}  {:6.3f}  {:6.3f} ; {:6.3f}  {:6.3f}  {:6.3f}  {:6.3f}".format( \
               x,g(x,0),g(x,1),g(x,2),g(x,3),f(x,0),f(x,1),f(x,2),f(x,3))
  print "\n\nNormal distribution sigma=2, mean=10:"
  g=Normal(2.0,10.0)
  for i in range(1,40):
    x=i/2.
    print "    {:4.1f}  {:7.4f}  {:7.4f}  {:7.4f}  {:7.4f}".format( \
               x,g(x,0),g(x,1),g(x,2),g(x,3))
  print "\n\nLogNormal distribution mode == 1.0:"
  g=LogNormal(1.0/3.0,1.0/9.0)
  for i in range(1,40):
    x=i/10.
    print "    {:4.1f}  {:7.4f}  {:7.4f}  {:7.4f}  {:7.4f}".format( \
               x,g(x,0),g(x,1),g(x,2),g(x,3))
  print "\n\nMaxwell distribution:"
  g=Maxwell(1.0/3.0)
  for i in range(1,40):
    x=i/2.
    print "    {:4.1f}  {:7.4f}  {:7.4f}  {:7.4f}  {:7.4f}".format( \
               x,g(x,0),g(x,1),g(x,2),g(x,3))
  print "\n\nBose distribution times cubic monomial:"
  g=Bose()*Monomial(3.0) # canonical black body
  for i in range(1,40):
    x=i/10.
    print "    {:4.1f}  {:7.4f}  {:7.4f}  {:7.4f}  {:7.4f}".format( \
               x,g(x,0),g(x,1),g(x,2),g(x,3))
  print "\n\nCIE-D65 BlackBody distribution MW/sr/sqm parameterized by eV:"
  g=BlackBody()
  for i in range(1,40):
    x=i/10.; sc=1e-6
    print "    {:4.3f}  {:7.4f}  {:7.4f}  {:7.4f}  {:7.4f}".format( \
               x,g(x,0)*sc,g(x,1)*sc,g(x,2)*sc,g(x,3)*sc)
  print "\n\nVisible black body power by temperature:"
  print "  kW/sqcm parameterized by eV:"
  g=BBpower()
  for i in range(500,10001,500):
    x=float(i)
    print "    {:> 8.1f}  {:9.4f}".format(x,g(x))
  print "  kW/sqcm parameterized by Thz:"
  g=BBpower(delta='Thz')
  for i in range(500,10001,500):
    x=float(i)
    print "    {:> 8.1f}  {:9.4f}".format(x,g(x))
  print "  kW/sqcm parameterized by nm:"
  g=BBpower(delta='nm')
  for i in range(500,10001,500):
    x=float(i)
    print "    {:> 8.1f}  {:9.4f}".format(x,g(x))
  print "  kW/sqcm parameterized by eV approximating 0 .. infinity:"
  g=BBpower(0,None,'eV'); s=phy.value('σ')
  def sb(T): return s*T*T*T*T*1e-7 # Stefan-Boltzmann, 0 to infinity
  for i in range(500,10001,500):
    x=float(i)
    print "    {:> 8.1f}  {:9.4f}  {:9.4f}".format(x,g(x),sb(x))
  print "\n\nIntegral of y=x:"
  g=Integral(Monomial(1.0))
  for i in range(1,10):
    x=i/1.
    print "    {:4.1f}  {:7.4f}  {:7.4f}  {:7.4f}  {:7.4f}".format( \
               x,g(x,0),g(x,1),g(x,2),g(x,3))
  print "\n\nIntegral of normal distribution:"
  g=Integral(Exp()(Monomial(2.0,-1.0)))
  for i in range(1,10):
    x=i/5.
    print "    {:4.1f}  {:7.4f}  {:7.4f}  {:7.4f}  {:7.4f}".format( \
               x,g(x,0),g(x,1),g(x,2),g(x,3))
  g.tolerance=.00000000000001; g.lower=g.delta=.001
  print "\n\nvarious 'infinite' Integral~s of exp(-x^2)"
  print "\n   lower == delta == .001      linear:",g(None)
  #g.lower=g.delta=.0001
  #print "\n   lower == delta == .0001     linear:",g(None)
  # g.delta=.0001  yields 0.886226925439 in about 4 seconds
  #g.lower=g.delta=.00001
  #print "\n   lower == delta == .00001    linear:",g(None)
  # g.delta=.00001  yields 0.886236925345 in about 8 seconds
  # g.delta=.000001 yields 0.886227924325 in about 2 minutes
  # exact: sqrt(pi)/2.0 == 0.8862269254527579
  g.geometric=True; g.lower=g.delta=.001
  print "\n   lower == delta == .001   geometric:",g(10)
  #g.geometric=True; g.lower=g.delta=.0001
  #print "\n   lower == delta == .0001  geometric:",g(10)
  # g.delta=.0001  yields 0.88612692693 in about 8 seconds
  #g.geometric=True; g.lower=g.delta=.00001
  #print "\n   lower == delta == .00001 geometric:",g(10)
  # g.delta=.00001  yields 0.886216925467 in about 20 seconds
  print "\n\nKullback-Leibler divergence of distributions underlying\n  'mathx.ferf' and 'math.erf':"
  m,a = 1.12837916709551257390,0.1009107485756095
  # simple direct implementation
  def q(x): csh=math.cosh(x*(a*x*x+m)); return (3.0*a*x*x+m)/csh/csh
  # constructive implementation: works
  #dhtn=Monomial(-2.0,4.0)(Exp(1.0)+Exp(-1.0))
  #q=Polynomial((3.0*a,0.0,m))*dhtn( Polynomial((a,0.0,m,0.0)) )
  #
  # direct implementation
  #h=Polynomial((a,0.0,m,0.0))
  #q=Polynomial((3.0*a,0.0,m))*Monomial(-2.0,1.0)(HypTrig()(h))
  def p(x): return 2.0*math.exp(-x*x)/math.sqrt(math.pi)# direct implementation is faster
  #p=lambda x: 2.0*(Normal(math.sqrt(.5))(x))           #  and more precise in last 4 digits
  d=Domain(-6.0,6.1,120)                                # kld integral muddies last 3 digits
  print "    In domain [-6.0, 6.1    by .1   ]: {!r}".format(kld(p,q,d))
  d=Domain(-6.0,6.01,1200)
  print "    In domain [-6.0, 6.01   by .01  ]: {!r}".format(kld(p,q,d))
  d=Domain(-6.0,6.001,12000)
  print "    In domain [-6.0, 6.001  by .001 ]: {!r}".format(kld(p,q,d))
  d=Domain(-6.0,6.0005,24000)
  print "    In domain [-6.0, 6.0005 by .0005]: {!r}".format(kld(p,q,d))
  #d=Domain(-6.0,6.0001,120000)
  #print "    In domain [-6.0, 6.0001 by .0001]: {!r}".format(kld(p,q,d))
  #for i in range(-60,60): print i/10.,p(i/10.),q(i/10.)
  print "\n\n===All tests complete==="
  # Conclusion here: limited spectrum has a complex, non-simple polynomial form
  # unlike the Stefan-Boltzmann law.  Stefan-Boltzmann corroborated within accuracy
  # of the discrete infinite integral approximation.
  '''
  from geometry.space import V,SqM
  g=BBpower()
  dom=[500.,5000.]#[500.,1000.,2000.,4000.,8000.,16000.,32000.]
  bbp=[(float(i),g(float(i))) for i in dom]
  #bbp=[(float(i),sb(float(i))) for i in dom] # sb == Stefan-Boltzmann
  #bbp=[(0.,0.)]+bbp; dom=[0.]+dom
  print repr(dom)
  print repr(bbp)
  v=V([math.log(t[1]) for t in bbp])
  print repr(v)
  #m=SqM([V([e**i for i in range(len(dom))]) for e in dom])
  m=SqM([V([1.0,math.log(e)]) for e in dom])
  print repr(m)
  u=m.inverse|v
  print repr(u)
  print map(repr,[math.exp(u[0]),u[1]])
  '''
