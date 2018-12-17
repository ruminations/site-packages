#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Package: geometry
Module : point.py
Website: https://github.com/ruminations/site-packages
License: https://github.com/ruminations/Licenses#design-license
Initial Copyright 2013

Module point defines three sub-types of 'space.V' that are four
dimensional. The components may be thought of as labelled [ w x y z ].
"""
__version__=20181216 # change imports to package imports
__version_history__=(20131105,20131222,20140703,20171211,20180128)
__all__=['D','P','Q','null','origin']
from math import sqrt
from geometry.gauss import C
try:
  #raise ImportError # Uncomment to use pure python
  from geometry.space_4D import V,SqM
  print "Using cython compiled space_4D.so in module point"
except ImportError:
  from geometry.space import V,SqM
  print "Using pure python space.py in module point"

class D(V):
  """Defines an affine space: the set of differentials on a space of
     absolute positions.  The positions implied are points of either
     type 'P' or 'Q' defined in this module.

     The object is represented as a quadruple with 0.0 for the first
     component.  Thus a 'D' instance represents a relative point in
     3-dimensional space, corresponding to a point on the plane at
     infinity for 'P' objects, and a pure quaternion for 'Q' objects.

     The first 0.0 component is largely invisible in the application
     interface: addition and subtraction of 'D' objects work as expected
     and iteration over a 'D' object yields 3 elements.

     The model for infinite elements returned by operations on the set
     of differentials is Gaussian.  That is, the infinite element is a
     singleton represented by 'None'.  The corresponding inversion is
     the zero differential predefined in this module as 'null'.

     Differentials 'D' may be thought of as relative positions used to
     define the bulk of a particular geometric configuration refered
     to salient positions of type 'P' which serve to position and orient
     the geometry in space.

     Implementation Notes:
     --------------------
     V.__init__ is over-ridden to enforce proper formation.  V.__add__
     and V.__sub__ are over-ridden to implement P.__radd__ and
     P.__rsub__ in a more secure way.  V.__iter__ is over-ridden to
     omit the first 0.0 element.

     Unary minus returns the negative differential as expected.
     Unary positive overides V.__pos__ and returns the implicit
     sum 'origin+self', which is a position of type 'P'.  Thus the
     idiom +d is useful for converting differentials to positions.

     Finally, V.__xor__ is over-ridden to implement the cross product
     of two differentials, which is well-defined only in this three
     dimensional affine space.
  """
  __slots__=[]

  def __init__(self,*u):
    """Actual parameters may be:
        0)  3 numeric values
        1)  a 3 element tuple or list or vector 'V' or a sub-type
       The value '0.0' is pre-pended to form a quadruple.
    """
    if not u: # class method '_' is caller: return empty instance
      super(D,self).__init__(); return # Realize superclass __slots__
    if len(u)==1: u=tuple(*u) # u[0] is a tuple, list, or V
    if len(u)==3: super(D,self).__init__(float(e) for e in ((0.0,)+tuple(u)))
    else: raise ValueError,"Differentials must have 3 spatial components."

  def __add__(self,u):
    if isinstance(u,D): return super(D,self).__add__(u)
    return u+self # __radd__

  def __iter__(self): return iter(self.p)

  def __lshift__(self,i):
    """Left shifts only the three dimensional relative position 'D.p' by 'i'."""
    v=self.p; l=len(v) # speedup: local lookup
    return D(v[i:l]+v[0:i])

  def __pos__(self):
    """Notational idiom defined for convenience: '(+d)' (unary positive)
       Return the position of type 'P' == 'origin+d', where 'd' is 'self'."""
    return P(self)

  def __rshift__(self,i):
    """Right shifts only the three dimensional relative position 'D.p' by 'i'."""
    v=self.p; l=len(v) # speedup: local lookup
    return D(v[l-i:l]+v[0:l-i])

  def __sub__(self,u):
    if isinstance(u,D): return super(D,self).__sub__(u)
    return -u+self # __rsub__

  def __xor__(self,d):
    if not isinstance(d,D): raise TypeError,\
      "The cross product is not defined for type: %s" % (type(d),)
    u0,u1,u2=self.p; v0,v1,v2=d.p
    return D(u1*v2-u2*v1,u2*v0-u0*v2,u0*v1-u1*v0)

  @property
  def i(self):
    """The x coordinate of the differential."""
    return self[1]

  @property
  def j(self):
    """The y coordinate of the differential."""
    return self[2]

  @property
  def k(self):
    """The z coordinate of the differential."""
    return self[3]

  @property
  def p(self):
    """The tuple representation of the three dimensional relative
       position of this differential 'D'."""
    return self.v[1:]

  @property
  def unit(self):
    """This differential normalized to unit magnitude.
       Returns 'None' for the zero differential."""
    try: return self/abs(self)
    except ZeroDivisionError: return None

null=D(0.0,0.0,0.0)

class P(V):
  """Defines a homogeneous representation [ w x y z ] of the projective
     three dimensional coordinate [ x/w y/w z/w ].  Points on the plane
     at infinity are represented by w==0.

     Positions 'P' may be thought of as absolute positions in world
     space used to define a reference frame for placing instances of
     objects specifying geometric configurations defined using
     differentials 'D'.  Positions 'P' thus define salient or special
     positions in world space.  Among those, the choice of origin is
     fundamental, and so that position is defined in this module and
     made available as 'origin'.  The implementation's underlying
     orthonormal Cartesian units implicitly define orientation.

     Implementation Notes:
     --------------------
     Operator over-rides implement the algebraic variation introduced
     by the homogeneous redundancy.  Most often, 'w==1' is 'True', but
     geometric interpretation of an operation can force renormalization.
     The '__invert__' method corresponding to '~' is defined to return
     a normalized representation as a 3-dimensional 'V' object where
     w==1 is implicit.

     Unary idioms are used to convert positions of type 'P' to
     differentials of type 'D':
       1) d=-p # set d == origin-p
       2) d=+p # set d == p-origin
     V.__neg__ and V.__pos__ are over-ridden to implement.
  """
  __slots__=['__p']

  def __init__(self,*u):
    """Actual parameters may be:
        0)  3 numeric values
        1)  a 3 element tuple or list or vector 'V' or sub-type
        2)  4 numeric values
        3)  a 4 element tuple or list or vector 'V' or sub-type
        4)  a scalar and a 3 element vector V or sub-type
        5)  a scalar and a 3 element tuple or list
       When a total of 3 elements are presented, w==1.0 is assumed.
    """
    self.__p=None # Realize __slots__
    if not u: # class method '_' is caller: return empty instance
      super(P,self).__init__(); return # Realize superclass __slots__
    if len(u)==1: u=tuple(*u) # u is a tuple, list, or V
    if len(u)==3: u=(1.0,)+tuple(u)
    elif len(u)==2:
      s,v=u # s is scalar; v is a tuple, list, or V
      u=(s,)+tuple(v)
    if len(u)==4:
      if round(u[0],15)==0.0: u=(0.0,)+u[1:] # unambiguous infinity
      super(P,self).__init__(float(e) for e in u)
    else: raise ValueError,"Homogeneous points must be 4-dimensional."

  def __abs__(self):
    """Return 'None' for a point at infinity."""
    if self.infinite: return None
    return abs(~self)

  def __add__(self,u):
    """If 'self' is infinite, self==return_value."""
    if isinstance(u,D): return super(P,self).__add__(self[0]*u)
    raise TypeError,"Point addition with %s is undefined." % (type(u),)

  def __eq__(self,u): return ~self==~u

  def __invert__(self):
    """Fold the homogeneous redundancy and return the resulting 3-vector.
       If the representation is of a point at infinity, return the
       3-vector representing the direction.
    """
    if self.__p is None:
      if self.infinite: self.__p=V(self.v[1:])
      else:
        self.__p=V((self/self[0]).v[1:])
    return self.__p

  def __ne__(self,u): return ~self!=~u

  def __neg__(self):
    """Notational idiom defined for convenience: '(-p)' (unary negative)
       Return the differential implied by the subtraction 'origin-self'.
       Infinite points return 'None': the differential representation
       of infinity."""
    if self.infinite: return None
    return D(-(~self))

  def __pos__(self):
    """Notational idiom defined for convenience: '(+p)' (unary positive)
       Return the differential implied by the subtraction 'self-origin'.
       Infinite points return 'None': the differential representation
       of infinity."""
    if self.infinite: return None
    return D(~self)

  def __sub__(self,u):
    """1) Subtracting a differential returns a point 'P'.
          If 'self' is infinite, self==return_value.
       2) Point subtraction returns a differential 'D'.
          One infinite point causes 'None' to be returned, i.e.
          the differential representation of infinity.
          Both infinite returns the differential, i.e. the differential
          representing the distance between the two points in the plane
          at infinity."""
    if isinstance(u,D): return super(P,self).__sub__(self[0]*u)
    if isinstance(u,P):
      if self.infinite ^ u.infinite: return None
      return D((~self)-(~u))
    raise TypeError,"Point subtraction with %s is undefined." % (type(u),)

  def __xor__(self,u):
    """Angles are well defined for all points, finite and infinite."""
    return (~self)^(~u)

  @property
  def p(self):
    """The tuple representation of the three dimensional position of
       this homogeneous point 'P'.  Returns a 4-tuple (0,x,y,z) for
       points at infinity."""
    if self.infinite: return self.v
    return (~self).v

  @property
  def infinite(self):
    """Returns 'True' for points on the plane at infinity."""
    return self[0]==0.0

origin=P(1.0,0.0,0.0,0.0)

class Q(V):
  """Defines a quaternion [ w x y z] where w,x,y,z are the coefficients
     of the corresponding h (scalar), i, j, k components.

     'Q' represents a four dimensional member of the set of quaternions.
     A unit magnitude quaternion, [ cos(t) a*sin(t) b*sin(t) c*sin(t) ]
     may be thought of as representing a rotation operator on a three
     dimensional differential 'D' about the unit magnitude differential
     [ a b c ] by an angle 2*t.

     Implementation Notes:
     --------------------
     Operator over-rides implement the algebraic variation introduced
     by the geometric interpretation.  A differential 'D' may be thought
     of as a pure quaternion where 'w==0' is 'True'.

     'V.__mul__' and 'V.__div__' are over-ridden since quaternion products
     and quotients are well defined.  The '__invert__' method implements
     the quaternion reciprocal or inverse.  Finally, the '__call__'
     method implements conjugation by a  quaternion: 'q(d)' applies 'q'
     as a rotation operator to the differential 'd' <==> 'q*d*~q'.
  """
  __slots__=['__d','__inv']

  def __init__(self,*u):
    """Actual parameters may be:
        0)  3 numeric values
        1)  a 3 element tuple or list or vector 'V' or sub-type
        2)  4 numeric values
        3)  a 4 element tuple or list or vector 'V' or sub-type
        4)  a scalar and a 3 element vector V or sub-type
        5)  a scalar and a 3 element tuple or list
       When a total of 3 elements are presented, w==0.0 is assumed.
    """
    self.__d=self.__inv=None # Realize __slots__
    if not u: # class method '_' is caller: return empty instance
      super(Q,self).__init__(); return # Realize superclass __slots__
    if len(u)==1: u=tuple(*u) # u is a tuple, list, or V
    if len(u)==3: u=(0.0,)+tuple(u)
    elif len(u)==2:
      s,v=u # s is scalar; v is a tuple, list, or V
      u=(s,)+tuple(v)
    if len(u)==4: super(Q,self).__init__(float(e) for e in u)
    else: raise ValueError,"Quaternions must be 4-dimensional."

  def __call__(self,d):
    """Apply 'self' as a rotation operator to the differential 'd'.
       That is, compute and return 'self*d*(~self)' as a differential."""
    if isinstance(d,D): return (self*d*(~self)).d
    raise TypeError, \
      "Quaternion rotation operand must be a differential type 'D'"

  def __div__(self,q):
    if isinstance(q,Q): return self*~q
    return super(Q,self).__div__(q)

  def __invert__(self):
    if self.__inv is None:
      self.__inv=super(Q,Q(self.h,-self.d)).__div__(self|self)
      self.__inv.__inv=self
    return self.__inv

  def __mul__(self,q):
    if isinstance(q,Q): return \
      Q(self.h*q.h-(self.d|q.d),q.d*self.h+self.d*q.h+(self.d^q.d))
    if isinstance(q,D): return \
      Q(-(self.d|q),q*self.h+(self.d^q))
    return super(Q,self).__mul__(q)

  __rmul__=__mul__

  @property
  def h(self):
    """The scalar component of the quaternion."""
    return self[0]

  @property
  def i(self):
    """The x coordinate of the vector component of the quaternion."""
    return self[1]

  @property
  def j(self):
    """The y coordinate of the vector component of the quaternion."""
    return self[2]

  @property
  def k(self):
    """The z coordinate of the vector component of the quaternion."""
    return self[3]

  @property
  def d(self):
    """The vector differential component of the quaternion."""
    if self.__d is None: self.__d=D(self.v[1:])
    return self.__d

  @property
  def conj(self):
    """The quaternion conjugate of this quaternion."""
    return Q(self.h,-self.d)

  @property
  def unit(self):
    """This quaternion normalized to unit magnitude.
       Returns 'None' for the zero quaternion."""
    try: return self/abs(self)
    except ZeroDivisionError: return None

  @property
  def sqrt(self):
    """Return the quaternion operator that, when applied twice, is
       equivalent to this operator applied once.  That is, for a unit,
       return a rotation operator about the same axis by half the angle.

       Returns 'None' for the zero quaternion or the gimbal lock
       quaternion: [-1 0 0 0]."""
    q=self.unit
    if q is None: return None # zero quaternion
    d=q.d; z=C(q[0]+1.,abs(d)).unit; d=d.unit # C(cos,sin)+C(1) unitized
    if z is None: return None # gimbal lock quaternion: [-1 0 0 0]
    if d is None: return Q(sqrt(abs(self)),q.d) # scalar sqrt
    return Q(z.x,d*z.y)*sqrt(abs(self))

  @property
  def SqM(self):
    """This quaternion converted to a homogeneous 4x4 matrix."""
    h,i,j,k=self.v; h2,i2,j2,k2=(e*e for e in self.v)
    ij,jk,ik,hi,hj,hk=i*j,j*k,i*k,h*i,h*j,h*k
    return SqM( V(1.,          0.,          0.,          0.), \
                V(0., h2+i2-j2-k2,  2.*(ij+hk),  2.*(ik-hj)), \
                V(0.,  2.*(ij-hk), h2-i2+j2-k2,  2.*(jk+hi)), \
                V(0.,  2.*(ik+hj),  2.*(jk-hi), h2-i2-j2+k2))

if __name__=='__main__':
  print "===Beginning module point.py test suite==="
  print "origin=",origin
  a=P(1,2,3); b=P((2,3,4)); c=P([3,4,5]); d=P(V(4,5,6))
  print tuple(repr(e) for e in (a,b,c,d))
  a=P(1,1,2,3); b=P((1,2,3,4)); c=P([1,3,4,5]); d=P(V(1,4,5,6))
  print tuple(repr(e) for e in (a,b,c,d))
  b=P(1,(2,3,4)); c=P(2.0,[3,4,5]); d=P(0.0,V(4,5,6))
  print tuple(repr(e) for e in (b,c,d))
  from sys import exc_info as ex
  try: a=P(1,2,3,4,5)
  except: print 'Caught:',ex()[1:2]
  try: a=P(1)
  except: print 'Caught:',ex()[1:2]
  print "===__init__ tests complete==="
  print repr(a)
  print repr(-a),repr(a+(a-b)),abs(a)
  print repr(a-b),repr(a-(a-b))
  print repr((a-b)-a),repr((a-b)+a)
  print a==a, a==b, a!=b, a!=a
  print a.p,a.infinite,a^b
  print repr(~c),repr(~d),c.p,d.p
  print c.infinite,d.infinite,c^d
  print "===class 'P' tests complete==="
  q=Q(1,2,3); i=Q(1,0,0)
  print q,i,'\n',q*i,i*q
  print i(q.d),q*q,q.d|q.d,'\n',q(i.d)
  print ~q,'\n',q*~q,~q*q
  print q.sqrt, q.sqrt*q.sqrt
  print i.sqrt, i.sqrt*i.sqrt
  print i.sqrt(D(0,1,0)), i(D(0,1,0))
  print q.SqM, q.unit.SqM|i.d,q.unit(i.d)
  try: a=Q(1,2,3,4,5)
  except: print 'Caught:',ex()[1:2]
  try: a=q(a)
  except: print 'Caught:',ex()[1:2]
  print "===class 'Q' tests complete==="
  diff=a-b
  print type(diff),repr(diff), abs(diff)
  print "===class 'D' tests complete==="
