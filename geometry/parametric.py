#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Package: geometry
Module : parametric.py
Website: https://github.com/ruminations/site-packages
License: https://github.com/ruminations/Licenses#design-license
Initial Copyright 2017

Provides a class 'Path' to encapsulate differential geometry calculations
for arbitrarily parameterized paths.  A useful mechanism to provide the
necessary derivatives is the class 'algebra.calculus.Function', as well
as other classes in the 'algebra.calculus' module.  Discrete data may be
employed via 'data.analysis.Function' and 'data.analysis.Interpolant'.

Defines class 'Ovoid' to implement some of the calculations necessary to
render a general ovoid as an SVG file.  This serves to provide a useful
non-trivial test framework for 'Path' and 'data.analysis.Function'.
"""
###TODO:
# 1) add 4th & 5th scalar curvature derivatives to 'Path'
# 2) finish lift & tilt properties of 'Path'
# 3) adapt to use 'Domain' class and automate extrema search in 'Domain'
# 4) add 'Ovoid' property that returns arc approximation
# 5) write some serious unit tests
from math import pi,sin,cos,sqrt
from geometry.gauss import C,Angle
from algebra.calculus import ConvergenceError
from data.analysis import Function # Ovoid extremum & tests use this line 339 e.g.

__all__=('Path','Ovoid')
__version__=20181209 # Finally moved discrete 'Function' code out
__version_history__=(20170208,20171224,20180309)

class Path(object):
  """Container for the differential curvature sequence of a path.
     Arbitrarily parameterized initialization function classes produce
     arc length parameterized attribute values.

     The user function parameter is truly arbitrary and is simply passed
     through to the user function definitions.  For example, if the user
     function expects an 'Angle' object, it is passed through so the
     user function may calculate its return value by accessing suitable
     attributes.

     In a more complex example, the arbitrary parameter may be a 'Q'
     (quaternion) object tracing a trajectory on S3, the 4-D sphere.
     The function may then calculate a differential 'D' on S2, the 3-D
     sphere that is the base space of a Hopf fibration.  The arc length
     trajectory and its differential properties thus describe a path on
     S2 embedded in R3, representative of the path traced on S3 by the
     actual parameter."""

  def __init__(self,p):
    """'p' is a class whose call method returns a position on a path by
       any convenient parameterization and whose index method returns
       a conveniently parameterized function yielding a derivative of
       the path, or a slice of such functions.  Indexing should retrieve
       [p,p',p'',p''',p'''',p'''''] in order to compute all attributes.
       Fewer derivatives imply fewer computable 'Path' attributes.

       NOTA BENE: derivatives MUST be formulated correctly.  Failing
       convergence of 'Function' object's Newton's method with accurate
       guesses is a clear indication of a faulty derivative formulation."""
    self._p=p

  def __call__(self,t,curvature=0):
    """Return the specified arc length parameterized curvature
       differential for the path corresponding to non-arc-length
       parameter 't'.  The curvatures are:
         'curvature'     differential
         -----------     ------------
           0              position
           1              tangent vector
           2              curvature vector
           3              torsion vector (unimplemented)
           4              tilt vector    (unimplemented)
        In two dimensions torsion and tilt are zero.  In three dimensions
        tilt is zero."""
    return self[curvature](t)

  def __getitem__(self,i):
    """Return a static function that evaluates the i-th arc length
       parameterized vector derivative value of the path when passed a
       user defined arbitrary parameterization value. 'i in (0,1,2,3,4,5)'."""
    return lambda t: self._curvatures[i](self,t)

  def __getslice__(self,*s):
    return (lambda t:f(self,t) for f in self._curvatures[slice(*s)])

  def _a0(self,t,order=0):
    """Return the weighted acceleration vector corresponding to
       parameter 't'.  This is A described in appendix II of the
       LyX document."""
    p1,p2 = (d(t) for d in self._p[1:3])
    return p2/(p1|p1)

  def _a1(self,t):
    """Return the first derivative of the weighted acceleration
       vector corresponding to parameter 't'.  This is A' described
       in appendix II of the LyX document."""
    p1,p2,p3 = (d(t) for d in self._p[1:4])
    n1,p21 = (p1|p1),(p2|p1)
    return p3/n1-((2.0/n1/n1)*p21)*p2

  def _a2(self,t):
    """Return the second derivative of the weighted acceleration
       vector corresponding to parameter 't'.  This is A'' described
       in appendix II of the LyX document."""
    p1,p2,p3,p4 = (d(t) for d in self._p[1:5])
    n1,n2,p21,p31 = (p1|p1),(p2|p2),(p2|p1),(p3|p1)
    return p4/n1-(2.0/n1/n1)*(2.0*p21*p3+(p31-4.0*p21*p21/n1+n2)*p2)

  def _a3(self,t):
    """Return the third derivative of the weighted acceleration
       vector corresponding to parameter 't'.  This is A''' described
       in appendix II of the LyX document."""
    p1,p2,p3,p4,p5 = (d(t) for d in self._p[1:6])
    n1,n2,p21,p31,p41,p32 = (p1|p1),(p2|p2),(p2|p1),(p3|p1),(p4|p1),(p3|p2)
    r=p21/n1; rp=p21*r
    return p5/n1-(2.0/n1/n1)*( ((3.0*p21)*p4+(3.0*(p31-4.0*rp+n2))*p3) \
                              +(p41+3.0*p32-12.0*r*(p31-2.0*rp+n2))*p2 )

  #def _u0(self,t): return self._k1(t)

  def _u1(self,t):
    """Return the first derivative of the unit tangent vector,
       corresponding to parameter 't'.  This is u0' described in
       appendix II of the LyX document."""
    return abs(self._p[1](t))*self._k2(t)

  def _u2(self,t):
    """Return the second derivative of the unit tangent vector,
       corresponding to parameter 't'.  This is u0'' described in
       appendix II of the LyX document."""
    u0,u1 = self._k1(t),self._u1(t)
    a0,a1 = self._a0(t),self._a1(t)
    return abs(self._p[1](t))*(a1-((a1|u0)+(a0|u1))*u0)

  def _u3(self,t):
    """Return the third derivative of the unit tangent vector,
       corresponding to parameter 't'.  This is u0''' described in
       appendix II of the LyX document."""
    u0,u1,u2 = self._k1(t),self._u1(t),self._u2(t)
    a0,a1,a2 = self._a0(t),self._a1(t),self._a2(t)
    return abs(self._p[1](t))*( a2 - ((a2|u0)+2.0*(a1|u1)+(a0|u2))*u0 \
                                  - ((a1|u0)+(a0|u1))*u1 + (a0|u0)*u2 )

  def _k0(self,t):
    """Return the position on the path corresponding to parameter 't'."""
    return self._p(t)

  def _k1(self,t):
    """Return the unit tangent to the path corresponding to parameter 't'."""
    return self._p[1](t).unit

  def _k2(self,t):
    """Return the curvature vector of the path corresponding to parameter 't'."""
    a,u0 = self._a0(t),self._k1(t)
    return a-u0*(a|u0)

  def _k3(self,t):
    """Return the torsion vector of the path corresponding to parameter 't'."""
    raise NotImplementedError

  def _k4(self,t):
    """Return the tilt vector of the path corresponding to parameter 't'."""
    raise NotImplementedError

  _curvatures=(_k0,_k1,_k2,_k3,_k4)

  @property
  def kappa(self):
    """Get a static function to evaluate the scalar curvature.
       Roots => geometric inflection points on the path."""
    return lambda t: abs(self._k2(t))

  def _kp1(self,t):
    u0,u1 = self._k1(t),self._u1(t)
    a0,a1,ks  = self._a0(t),self._a1(t),self.kappa(t)
    unscaled=((a1|a0)-(u0|a0)*((u1|a0)+(u0|a1))) # check 0.0/0.0 case
    return unscaled/ks if unscaled else 0.0      # return 0.0 for 0.0/0.0
    # version restricted to two dimensions:
    #xp,yp = (self._p[1](t)).xy; xpp,ypp = (self._p[2](t)).xy
    #x3p,y3p = (self._p[3](t)).xy
    #n=xp*xp+yp*yp
    #return (xp*y3p-yp*x3p)/n/sqrt(n)-3.0*self.kappa(t)*(xp*xpp+yp*ypp)/n

  @property
  def kappa_p1(self):
    """Get a static function to evaluate the first derivative of
       the  scalar curvature.  Roots => curvature extrema."""
    return lambda t:self._kp1(t)

  def _kp2(self,t):
    u0,u1,u2 = self._k1(t),self._u1(t),self._u2(t)
    a0,a1,a2 = self._a0(t),self._a1(t),self._a2(t)
    ks,kps,s = self.kappa(t),self._kp1(t),(u1|a0)+(u0|a1)
    return ((a2|a0)+(a1|a1)-s*s-(u0|a0)*((u2|a0)+2.0*(u1|a1)+(u0|a2))-kps*kps)/ks
    # version restricted to two dimensions:
    #xp,yp = (self._p[1](t)).xy; xpp,ypp = (self._p[2](t)).xy
    #x3p,y3p = (self._p[3](t)).xy; x4p,y4p = (self._p[4](t)).xy
    #n=xp*xp+yp*yp; n2=sqrt(n); m=(xp*xpp+yp*ypp)/n
    #return (xpp*y3p+xp*y4p-ypp*x3p-yp*x4p)/n/n2 - \
    #       3.0*( 2.0*self.kappa_p1(t)*m + self.kappa(t)*( \
    #             (xpp*xpp+ypp*ypp+xp*x3p+yp*y3p)/n+m*m  ) )

  @property
  def kappa_p2(self):
    """Get a static function to evaluate the second derivative of
       the  scalar curvature.  Roots => curvature inflectia."""
    return lambda t:self._kp2(t)

  def _kp3(self,t):
    u0,u1,u2,u3 = self._k1(t),self._u1(t),self._u2(t),self._u3(t)
    a0,a1,a2,a3 = self._a0(t),self._a1(t),self._a2(t),self._a3(t)
    ks,kps,kpps = self.kappa(t),self._kp1(t),self._kp2(t)
    return ( (a3|a0) + 3.0*( (a2|a1) - kpps*kps \
                - ((u1|a0)+(u0|a1))*((u2|a0)+2.0*(u1|a1)+(u0|a2)) ) \
                - (u0|a0)*((u3|a0)+3.0*((u2|a1)+(u1|a2))+(u0|a3)) )/ks
    # version restricted to two dimensions never implemented

  @property
  def kappa_p3(self):
    """Get a static function to evaluate the third derivative of
       the  scalar curvature."""
    return lambda t:self._kp3(t)

class Ovoid(object):
  """A generalized ellipse:
     1) generated from two circles of radius 'a' and 'b' with centers
        separated by 'c' such that 0 <= c <= a and 0 < b <= a.
     2) The x-coordinate is determined by the intersection of a pencil
        of lines through the center of the circle with radius 'b' with
        the circle of radius 'a'.
     3) The corresponding y-coordinate is determined by the intersection
        of the same pencil with the circle of radius 'b'.
     4) If, in addition, a = b + c, then the two circles are tangent
        at one point and ovoid is called an ooid.
    The resulting forms include ordinary circles and ellipses as well
    as a diversity of egg shapes that are continuously differentiable.
    Parameterization is via 'Angle' objects, but functions will
    correctly evaluate for radian valued scalars."""
  __slots__=['_a','_b','_c','_path']; _p=14

  def __init__(self,a,b,c):
    """'a' is the radius of the large circle; 'b' is the radius of the
        small circle; 'c' is the distance between their centers."""
    if not(0.0<=c and c<=a and 0.0<b and b<=a): raise ValueError, \
      "Must have 0 <= (c=={}) <= (a=={}) and 0 < (b=={}) <= (a=={}).".format(c,a,b,a)
    self._a,self._b,self._c = a,b,round(c,self._p) # clarify divide by zero case
    self._path=Path(self)

  def __call__(self,theta,derivative=0):
    """Return the complex position of the ovoid (default) for angular
       parameter 'theta' as a 'C' object, or one of its complex
       derivatives (1,2,3,4).  The Complex values are interpreted as
       2-dimensional vectors."""
    return self[derivative](theta)

  def __getitem__(self,i):
    """Return a static function that evaluates the i-th 'Angle' object
       parameterized complex derivative value of the 'Ovoid' path in the
       complex plane.  'i in (0,1,2,3,4,5)' is 'True'."""
    return lambda t: self._derivatives[i](self,t)

  def __getslice__(self,*s):
    return (lambda t:f(self,t) for f in self._derivatives[slice(*s)])

  @staticmethod
  def _rad(t):
    """Extract floating point radian value from an 'Angle' object."""
    try: return t.rad
    except AttributeError: return float(t)

  def _aux(self,theta,values=1):
    """Return the tuple of sin and the auxilliary function
         s,f = sin(theta),sqrt(a*a-c*c*sin(theta)*sin(theta))
       and sequential derivatives if values > 1 : fp, fpp, f3p, and f4p.
       Appends alternately sin(2*theta) or cos(2*theta) to the tuple."""
    s=sin(theta); a,c = self._a,self._c
    f=sqrt(a*a-c*c*s*s); values=int(values)-1
    if not values: return s,f
    norm=not(round(f,15)==0.0)
    theta=2.0*theta; s2=sin(theta)
    fp = -c*c*s2/f/2.0 if norm else 0.0; values-=1
    if not values: return s,f,fp,s2
    c2=cos(theta)
    fpp = -(c*c*c2+fp*fp)/f if norm else 0.0; values-=1
    if not values: return s,f,fp,fpp,c2
    fpp_f = fpp/f if norm else 0.0
    sub_ex = -3.0*fpp_f-4.0 if norm else 0.0
    f3p=fp*sub_ex; values-=1
    if not values: return s,f,fp,fpp,f3p,s2
    fp_f = fp/f if norm else 0.0;
    f4p=fpp*sub_ex + 12.0*fp_f*fp*(fpp_f+1.0); values-=1
    if not values: return s,f,fp,fpp,f3p,f4p,c2
    f5p=16.0*fp+10.0*fpp_f*f3p - 5.0*fp_f*f4p; values-=1
    if not values: return s,f,fp,fpp,f3p,f4p,f5p,s2
    raise ValueError,"{} values requested; only 6 defined.".format(values+4)

  def _d0(self,theta):
    """Return the complex plane position for angular parameter 'theta'."""
    theta=self._rad(theta); c=cos(theta); s,f = self._aux(theta)
    return C(c*(self._c*c+f),self._b*s)

  def _d1(self,theta):
    """Return the complex first derivative for angular parameter 'theta'."""
    theta=self._rad(theta); c=cos(theta); s,f,fp,s2 = self._aux(theta,2)
    return C(fp*c-f*s-s2*self._c,self._b*c)

  def _d2(self,theta):
    """Return the complex second derivative for angular parameter 'theta'."""
    theta=self._rad(theta); c=cos(theta); s,f,fp,fpp,c2 = self._aux(theta,3)
    return C((fpp-f)*c-2.0*fp*s-2.0*c2*self._c,-self._b*s)

  def _d3(self,theta):
    """Return the complex third derivative for angular parameter 'theta'."""
    theta=self._rad(theta); c=cos(theta); s,f,fp,fpp,f3p,s2 = self._aux(theta,4)
    return C(4.0*s2*self._c+(f3p-3.0*fp)*c-(3.0*fpp-f)*s,-self._b*c)

  def _d4(self,theta):
    """Return the complex fourth derivative for angular parameter 'theta'."""
    theta=self._rad(theta); c=cos(theta); s,f,fp,fpp,f3p,f4p,c2 = self._aux(theta,5)
    return C(8.0*c2*self._c+(f4p-6.0*fpp+f)*c-4.0*(f3p-fp)*s,self._b*s)

  def _d5(self,theta):
    """Return the complex fifth derivative for angular parameter 'theta'."""
    theta=self._rad(theta); c=cos(theta); s,f,fp,fpp,f3p,f4p,f5p,s2 = self._aux(theta,6)
    return C(-16.0*s2*self._c+(f5p-10.0*f3p+5.0*fp)*c-(5.0*f4p-10.0*fpp+f)*s,self._b*c)

  _derivatives=(_d0,_d1,_d2,_d3,_d4,_d5)

  @property
  def curvature_minima(self):
    """Return a sorted tuple of 'Angle' objects representing the
       parameters where the 'Ovoid' curvature has a local minimum."""
    angles=[Angle(i,'deg') for i in xrange(-170,185,10)]
    k=self._path.kappa; f=Function(self._path.kappa_p1,self._path.kappa_p2)
    positions=[round(k(a),15) for a in angles]
    minima=[angles[i] for i,p in enumerate(positions) \
              if p<positions[i-1] and p<positions[(i+1)%36]]
    return sorted(Angle(f.root(g.rad),'rad') for g in minima)

  @property
  def curvature_maxima(self):
    """Return a sorted tuple of 'Angle' objects representing the
       parameters where the 'Ovoid' curvature has a local maximum."""
    angles=[Angle(i,'deg') for i in xrange(-170,185,10)]
    k=self._path.kappa; f=Function(self._path.kappa_p1,self._path.kappa_p2)
    positions=[round(k(a),15) for a in angles]
    maxima=[angles[i] for i,p in enumerate(positions) \
              if p>positions[i-1] and p>positions[(i+1)%36]]
    return sorted(Angle(f.root(g.rad),'rad') for g in maxima)

  @property
  def curvature_extrema(self):
    """Return a sorted tuple of 'Angle' objects representing the
       parameters where the 'Ovoid' curvature has a local extremum."""
    return sorted(self.curvature_minima+self.curvature_maxima)

  @property
  def curvature_inflectia(self):
    """Return a sorted tuple of 'Angle' objects representing the
       parameters where the 'Ovoid' curvature has an inflection point."""
    angles=[Angle(i,'deg') for i in xrange(-170,185,10)]
    k=self._path.kappa_p1; f=Function(self._path.kappa_p2,self._path.kappa_p3)
    positions=[round(k(a),15) for a in angles]
    minima=[angles[i] for i,p in enumerate(positions) \
              if p<positions[i-1] and p<positions[(i+1)%36]]
    maxima=[angles[i] for i,p in enumerate(positions) \
              if p>positions[i-1] and p>positions[(i+1)%36]]
    return sorted(Angle(f.root(g.rad),'rad') for g in minima+maxima)

  @property
  def curvature_critica(self):
    """Return a sorted tuple of tuples '(Angle,type)' objects representing
       the parameters where the 'Ovoid' curvature has a critical point.
         type value    critical point type
         ----------    -------------------
             -1             minimum
              0          inflection point
              1             maximum"""
    return sorted([(c,-1) for c in self.curvature_minima]+ \
                  [(c, 1) for c in self.curvature_maxima]+ \
                  [(c, 0) for c in self.curvature_inflectia])

if __name__=="__main__":
  o=Ovoid(2.0,1.0,1.0); p=Path(o) # ooid
  #o=Ovoid(2.0,1.0,0.0); p=Path(o) # ellipse
  '''
  for i in range(-10,11):
    a=Angle(float(i)/10.0)
    print a.deg
    print "    p",o(a),p(a)
    print "   dp",o[1](a),p[1](a),abs(p[1](a))
    print "  dpp",o[2](a),p[2](a),1/abs(p[2](a))
    print "  d3p",o[3](a)
    print "  d4p",o[4](a)
  scalar_k=Function(p.kappa)
  scalar_kp=Function(p.kappa_p1,p.kappa_p2)
  scalar_kpp=Function(p.kappa_p2)
  #for d in range(-180,185,5): print d,p.kappa_p1(Angle(d,'deg'))
  #for d in range(-17,17): print d,p.kappa_p1(Angle(d,'deg'))
  print
  #print Angle(Function(p.kappa,p.kappa_p1,p.kappa_p2).extremum(2.3),'rad').deg
  #extrema=[scalar_kp.root(Angle(d*pi/180).rad) for d in range(-180,185,5)]
  #extrema=sorted([deg(v) for v in extrema if v is not None])
  #inflectia=[scalar_kpp.root(Angle(d*pi/180).rad) for d in range(-180,185,5)]
  #inflectia=sorted([deg(v) for v in extrema if v is not None])
  ''''''#0.1846793
  print "extrema"
  for i in range(len(extrema)):
    print "   ",extrema[i]
  print "inflectia"
  for i in range(len(inflectia)):
    print "   ",inflectia[i]
  #'''
  f=Function(lambda t: p._u2(t).x)#p._kp2)#
  print map(lambda a: a.deg,o.curvature_extrema)
  print map(lambda a: a.deg,o.curvature_inflectia)
  for t in o.curvature_critica: print "type: {1:3d}  location: {0!r}".format(*t)
  '''
  for d in range(-175,185,5):
     #print d,p.kappa(Angle(d,'deg')),p._kp1(Angle(d,'deg')), \
     #        p._kp2(Angle(d,'deg')),p._kp3(Angle(d,'deg')),f._df((Angle(d,'deg')))/f._dx
    print d,p._u3(Angle(d,'deg')).x,f._df(Angle(d,'deg'))/f._dx,p._u2(Angle(d,'deg')).x
  '''
  from math import exp
  p=Path(Function(exp,exp,exp,exp,exp,exp))
  f=Function(p.kappa,p.kappa_p1,p.kappa_p2,p.kappa_p3)
  print "exponential curvature extremum:   {!r})".format(f(f.extremum(0.0)))
  print "exponential y-value at {!r}: {!r}".format(f.extremum(0.0),exp(f.extremum(0.0)))
  print
  print "exponential curvature inflection: {!r})".format(f(f.inflection(-1.0)))
  print "exponential y-value at {!r}: {!r}".format(f.inflection(-1.0),exp(f.inflection(-1.0)))
  print
  print "exponential curvature inflection: {!r})".format(f(f.inflection(.7)))
  print "exponential y-value at {!r}: {!r}".format(f.inflection(.7),exp(f.inflection(.7)))
  print
  from math import degrees as deg
  a,w = 1.0,1.0; w2=w*w; w3=w2*w; w4=w2*w2; w5=w3*w2; g=.5/a
  p=Path(Function(lambda t:  a*   cos(w*t), lambda t: -a*w* sin(w*t), \
                  lambda t: -a*w2*cos(w*t), lambda t:  a*w3*sin(w*t), \
                  lambda t:  a*w4*cos(w*t), lambda t: -a*w5*sin(w*t)))
  f=Function(p.kappa,p.kappa_p1,p.kappa_p2,p.kappa_p3)
  print "cosine curvature extremum:   {!r})".format(f(f.extremum(0.0)))
  print "cosine y-value at {!r} (degrees): {!r}".format(deg(f.extremum(0.0)),cos(f.extremum(0.0)))
  print
  print "cosine curvature inflection: {!r})".format(f(f.inflection(-g)))
  print "cosine y-value at {!r} (degrees): {!r}".format( \
           deg(f.inflection(-g)),a*cos(w*f.inflection(-g)))
  print
  print "cosine curvature inflection: {!r})".format(f(f.inflection(g)))
  print "cosine y-value at {!r} (degrees): {!r}".format( \
           deg(f.inflection(g)),a*cos(w*f.inflection(g)))
  print
  o=Ovoid(1.0,1.0,1.0) # pathological edge case
  print o.curvature_critica
  print
  from math import cosh,sinh
  p=Path(Function(sinh,cosh,sinh,cosh,sinh,cosh))
  f=Function(p.kappa,p.kappa_p1,p.kappa_p2,p.kappa_p3)
  print "sinh curvature extremum:   {!r})".format(f(f.extremum(1.0)))
  print "sinh y-value at {!r}: {!r}".format(f.extremum(1.0),sinh(f.extremum(1.0)))
  print
  print "sinh curvature inflection: {!r})".format(f(f.inflection(1.0)))
  print "sinh y-value at {!r}: {!r}".format(f.inflection(1.0),sinh(f.inflection(1.0)))
  print
  print "sinh curvature extremum: {!r})".format(f(f.extremum(-1.0)))
  print "sinh y-value at {!r}: {!r}".format(f.extremum(-1.0),sinh(f.extremum(-1.0)))
  print
  print "sinh curvature inflection: {!r})".format(f(f.inflection(-1.0)))
  print "sinh y-value at {!r}: {!r}".format(f.inflection(-1.0),sinh(f.inflection(-1.0)))
  print
  #print "sinh curvature extremum: {!r})".format(f(f.extremum(0.0)))
  #print "sinh y-value at {!r}: {!r}".format(f.extremum(0.0),sinh(f.extremum(0.0)))
  # zero division - issue is sharp corner
  print
  #for i in range(-30,30): print i/100.,p.kappa_p1(i/100.)
