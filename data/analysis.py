#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Package: data
Module : analysis.py
Website: https://github.com/ruminations/site-packages
License: https://github.com/ruminations/Licenses#design-license
Initial Copyright 2018

This module is a repository for discrete data smoothing operations providing
two callable classes:

Function     : calculate critical points using discrete numerical
               derivatives derived from a continuous function.

Interpolant  : linearly interpolate a function defined by a data set.
"""
from geometry.gauss import C,Angle
from algebra.calculus import ConvergenceError
__all__=('Function','Interpolant')
__version__=20181216
__version_history__=(20180506,)
###TODO:
# 0) as always, writing a doctest suite.

class Function(object):
  """Instances of this class evaluate roots, extrema, and inflectia
     near a guess using Newton's method and numerically evaluated
     derivatives.  Roots have about 8 digits of accuracy and are
     rounded to 8 decimal places.  Optionally, analytic derivatives
     may be supplied at instantiation, yielding double precision results.

     Bad guesses that infinite loop raise 'ConvergenceError' after the
     property 'limit' iterations which defaults to 40.

     Infinite loops result when a guess is at an extremum of the primary
     function corresponding to the type of critical point that is sought.
       e.g. an extremum of the base function when a root is sought, or
            an extremum of the first derivative when a base function
               extremum is sought, or
            an extremum of the second derivative when a base function
               inflection is sought.

     Infinite looping can occur for some guesses for functions that are
     symmetric.  For periodic functions, guesses near an extremum will
     converge after some 10's of iterations to some angular value well
     outside the periodic interval.  It is more accurate to calculate
     the value from some guess that converges inside the interval, and
     that can be done by using the first result, modulo the interval, as
     a better guess.  This must be implemented by the user, since neither
     periodicity nor its interval can be infered.

     Infinite looping can also occur for functions that are asymptotic
     to zero at either infinity, and may raise a 'ConvergenceError' due
     to a 'ZeroDivisionError', depending on how rapidly the computation
     diverges.

     A 'ConvergenceError' caused by a 'ZeroDivisionError' can also mean
     the termination tolerance is set too small.  Nota Bene: incorrect
     derivatives will certainly raise errors.

     The class is convenient when the functional form is complicated
     and the derivatives, in consequence, are extremely complicated
     and therefore inconvenient to formulate.  Curvature and torsion
     calculations typically qualify.  It is also useful for continuous
     function definitions for which no closed form derivatives exist.

     When computing the extrema and inflectia of the scalar curvature of
     a 'parametric.Path', Newton's method can converge much more slowly
     than the usual doubling of accuracy with each iteration.  This is
     partly due to the gradual variations of curvature, especially for
     closed paths.  Another reason is that positions underlying the scalar
     curvature are in a space of two or more dimensions, and Newton's
     method is, in general, fractal in such spaces.  Reducing the spatial
     vector to a scalar thus introduces some floating point rounding noise.

     It is possible to get double precision accuracy with these algorithms
     using the module 'decimal' and 60 decimal places of precision, but
     this generally requires the passed in function to be evaluated in
     conformance to the expectations of the 'decimal' module.

     This more limited implementation gives results sufficiently
     accurate to determine reasonable spliced arc approximations to a
     curve for the purposes of display."""

  def __init__(self,f,fp=None,f2p=None,f3p=None,f4p=None,f5p=None):
    """'f' is a scalar valued function of a scalar variable for which
       it is desired to numerically calculate critical points.

       Optional parameters are first, second, third, fourth and fifth
       derivative functions, which may be supplied if they are known.
       If 'fp' is supplied,
         'Function.root' has double precision accuracy.
       If 'fp' and 'f2p' are supplied,
         'Function.extremum' has double precision accuracy.
       If 'f2p' and 'f3p' are supplied,
         'Function.inflection' has double precision accuracy.
       Derivatives beyond those are necessary to exploit the full benefits
       of a 'parametric.Path' object formed using a 'Function' instance."""
    self._dfs=(f,fp,f2p,f3p,f4p,f5p); self._max_loop=40; self._p=8
    self._t=self._dx=.5**26; self._d_5=self._dx*.5; self._d1_5=self._dx+self._d_5
    self._at=.0000000000000003 # 15.7 decimal places

  def __call__(self,x,derivative=0):
    """Return the complex position of the trivial parameterization of the
       indicated derivative of the function represented by this instance
       at the parameter 'x', e.g.:
         derivative    return value
         ----------    ------------
              0         C(x  ,  f(x))
              1         C(1.0, fp(x))
              >1        C(0.0,f<d>p(x)), where <d>==derivative"""
    if derivative==0: return C(x,self._dfs[0](x))
    if derivative==1: return C(1.0,self._dfs[1](x))
    return C(0.0,self._dfs[derivative](x))

  def __getitem__(self,i):
    """Return a static function that evaluates the i-th trivially
       parameterized complex derivative value of the path as a 'C' object
       when passed a abcissa value 'x'. 'i in (0,1,2,3,4,5)'."""
    if i==0: return lambda x: C(x,self._dfs[0](x))
    if i==1: return lambda x: C(1.0,self._dfs[1](x))
    return lambda x: C(0.0,self._dfs[i](x))

  def __getslice__(self,*s):
    s=slice(*s)
    i = 0 if s.start is None else s.start
    j = 5 if s.stop  is None else s.stop
    k = 1 if s.step  is None else s.step
    return (self[l] for l in range(i,j,k))

  def _df(self,t):
    """Discrete differential of function 'f' at 't'."""
    f,d_5 = self._dfs[0],self._d_5
    return f(t+d_5)-f(t-d_5)

  def _d2f(self,t):
    """Second order discrete differential of function 'f' at 't'."""
    f,dx = self._dfs[0],self._dx
    return f(t+dx)+f(t-dx)-2.0*f(t) # (f(t+dx)-f(t))-(f(t)-f(t-dx))

  def _d3f(self,t):
    """Third order discrete differential of function 'f' at 't'."""
    f,d_5,d1_5 = self._dfs[0],self._d_5,self._d1_5
    return (f(t+d1_5)-f(t-d1_5))+3.0*(f(t-d_5)-f(t+d_5))

  def _nm(self,f,fp,g=1.0):
    """Analytic Newton's method."""
    delta,epsilon = 1.0,self._at
    i,limit = 0,self._max_loop
    while abs(delta)>epsilon:
      try: delta=f(g)/fp(g)
      except ZeroDivisionError:
        try: ex=ConvergenceError(g,(f(g),fp(g)),i)
        except ZeroDivisionError:
          print "g,f(g),i; fp(g)->ZeroDivision:",g,f(g),i
          raise
        print ex.args,repr(ex.message),"done printing"
        raise ex
      g-=delta; i+=1
      if i==limit: raise ConvergenceError(g,(f(g),fp(g),delta),i)
    return g

  def _dnm(self,df,ddf,g=1.0):
    """Numerical Newton's method."""
    delta,i,limit,dt = 1.0,0,self._max_loop,self._t
    while abs(delta)>dx:
      try: delta=( dx/ddf(g) )*df(g) # order avoids some under/overflow
      except ZeroDivisionError: raise ConvergenceError(g,((df(g),dx),(ddf(g),dx)),i)
      g-=delta; i+=1
      if i==limit: raise ConvergenceError(g,((df(g),dx),(ddf(g),dx)),i)
    return round(g,self._p)

  def extremum(self,guess):
    """Return a value near 'guess' where the first derivative of the
       function evaluates to zero."""
    fp,f2p = self._dfs[1:3]
    if (fp is not None) and (f2p is not None):
      return self._nm(fp,f2p,guess)
    return self._dnm(self._df,self._d2f,guess)

  def inflection(self,guess):
    """Return a value near 'guess' where the second derivative of the
       function evaluates to zero."""
    f2p,f3p = self._dfs[2:4]
    if (f2p is not None) and (f3p is not None):
      return self._nm(f2p,f3p,guess)
    return self._dnm(self._d2f,self._d3f,guess)

  def root(self,guess):
    """Return a value near 'guess' where the function evaluates to zero."""
    f,fp = self._dfs[:2]
    if fp is not None: return self._nm(f,fp,guess)
    return self._dnm(f,self._df,guess)

  @property
  def differential(self):
    """Get or set the discrete first derivative differential.  This
       automatically adjusts the differentials used for higher order
       derivatives.  Defaults to .5**26."""
    return self._dx
  @differential.setter
  def differential(self,dx):
    self._dx=float(dx); self._d_5=self._dx*.5; self._d1_5=self._dx+self._d_5

  @property
  def limit(self):
    """Get or set the maximum number of iterations allowed before
       raising a 'ConvergenceError'.  Defaults to 40."""
    return self._max_loop
  @limit.setter
  def limit(self,iterations): self._max_loop=int(iterations)

  @property
  def rounding(self):
    """Get or set the number of decimal places return values from
       discrete calculations are rounded to.  Defaults to 8."""
    return self._p
  @rounding.setter
  def rounding(self,places): self._p=int(places)

  @property
  def discrete_tolerance(self):
    """Get or set the tolerance used to terminate iteration using
       discrete derivatives.  Defaults to dx."""
    return self._t
  @discrete_tolerance.setter
  def discrete_tolerance(self,dt): self._t=float(dt)

  @property
  def tolerance(self):
    """Get or set the tolerance used to terminate iteration using
       analytic derivatives.  Defaults to .0000000000000003 (15.7
       decimal places)."""
    return self._at
  @tolerance.setter
  def tolerance(self,at): self._at=float(at)


class Interpolant(object):
  """Represent a discrete data set as a C_0 function using linear interpolation.
     The dataset abcissae need not be uniformly spaced."""
  def __init__(self,dataset,default=None):
    """'dataset' is an iterable yielding pairs
          (<scalar_abcissa>,<function_value>)
       where <function_value> need not be scalar, but must implement
       arithmetic operations on the class instances.  The pairs are
       sorted on <scalar_abscissa> prior to use.  'default' is the value
       returned for abscissae outside the range of abcissae in 'dataset'."""
    self._ds=list(dataset); self._ds.sort(key=lambda v:v[0]); self._d=default
    if len(self._ds)<2: raise ValueError( \
      "Dataset length {} must be at least 2.".format(len(self._ds)))
    if any(len(e)!=2 for e in self._ds): raise ValueError( \
      "Initialization dataset must be an iterable yielding pairs.")
    self._ds=tuple((float(a),o) for a,o in self._ds) # sanitize
    self._minabcissa,self._maxabcissa = self._ds[0][0],self._ds[-1][0]
  def __call__(self,t):
    """Evaluate the 'Interpolant' at abcissa parameter 't'.  Values of 't'
       outside the domain of the dataset abcissa return initialization 'default'."""
    if not t in self: return self._d
    a,b = self._pair(t); xa,ya = a; xb,yb = b
    return ya+(yb-ya)*((t-xa)/(xb-xa)) # ya + dy*(dt/dx)
  def __contains__(self,t):
    """Return 'True' if 't' is >= the minimum and <= the maximum of the abscissae."""
    return t>=self._minabcissa and t<=self._maxabcissa
  def __len__(self):
    """Length of the initial dataset."""
    return len(self._ds)
  def _pair(self,t):
    """Get the two sequential dataset pairs defining an abcissa interval
       containing 't' using binary search."""
    i,j = 0,len(self._ds)-1
    while i+1!=j:
      k=i+j>>1
      if   self._ds[k][0]>t: j=k
      elif self._ds[k][0]<t: i=k
      else: # equal
        if   self._ds[j][0]>t: return self._ds[k],self._ds[k+1]
        elif self._ds[i][0]<t: return self._ds[k-1],self._ds[k]
    return self._ds[i],self._ds[j]
  @property
  def domain(self):
    """Get a tuple '(<minimum_abscissa>,<maximum_abcissa>)'"""
    return self._minabcissa,self._maxabcissa


if __name__=="__main__":
  from data.cie import lms_10deg_rectified as lms
  i=Interpolant(lms)
  l,m = map(int,i.domain)
  print all(abs(a-b)==0.0 for a,b in zip([y for x,y in lms],[i(x) for x in range(l,m+1)]))
