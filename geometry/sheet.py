#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Package: geometry
Module : sheet.py
Website: https://github.com/ruminations/site-packages
License: https://github.com/ruminations/Licenses#design-license
Initial Copyright 2013

Module 'sheet' defines three dimensional surface objects in homogeneous
projective space together with containment operators.

In projective geometry, there is the plane at infinity, which every line
outside that plane intersects in a single point.  The 'two' infinite
points on either end of a line are identified with the intersection.
This preserves projective duality: every point is a pencil of lines, and
every line is a pencil of points, even in the plane at infinity.

In classic projective geometry, a plane can not divide space, for it is
always possible to go around through infinity.  All space is contained
in the plane at infinity, so all planes are equivalent: they all contain
all of space.

Space can, however, be divided by two planes.  These two planes always
intersect in a line (possibly on the plane at infinity).  For example,
two parallel planes divide space into the region between them, and the
region outside them (which is connected through the plane at infinity.

Here, we adopt the convention that the plane at infinity is considered
to be implicitly part of every surface.  In consequence, a plane divides
space into two parts, consistent with Euclidean intuition.  The same
convention is applied to any surface, finite or infinite.  One side of
an object of infinite extent is arbitrarily chosen as an 'interior' and
the other side is 'exterior'.  Of course, closed surfaces always divide
space into an interior and an exterior, and which is deemed interior is
equally arbitrary.  Usually the finite volume is called 'interior'.

Using these ideas, every sheet defines the method '__contains__', so for
a position p and a sheet s 'p in s' is a boolean that returns 'True' for
homogeneous coordinates that are part of the strictly interior volume.

Every sheet also defines the '__call__' method so 's(p)' returns zero
when 'p' is on 's', '-1' when 'p' is an interior point, and '1' when 'p'
is an exterior point.

Finally, every sheet defines the '__pos__' method (unary positive), so
'+s' causes the sheet to invert its current interior and exterior
conventions and ++s leaves them unchanged.
"""
__version__=20181216 # change imports to package imports
__version_history__=(20131106,20131204,20140613,20171205,20180410)
__all__=['Plane','Quadric']
try:
  #raise ImportError # Uncomment to use pure python
  from geometry.space_4D import V,SqM
  print "Using cython compiled space_4D.so in module sheet"
except ImportError:
  from geometry.space import V,SqM
  print "Using pure python space.py in module sheet"
from algebra.mathx import sign
from geometry.point import P,D,origin
## following flags error when importing geometry.figure
import geometry.figure as figure# avert circular import: used in
#try: from geometry.figure import Line,Circle
#except ImportError: pass # skip circular import second pass

def flip_sign(x):
  """Return: 1 => x < 0; 0 => x==0; -1 => x > 0"""
  return -sign(x)

class Plane(V):
  """Defines a homogeneous representation [ d a b c ] of the algebraic
     definition of a plane: ax+by+cz+d=0.  Alternatively, if p is a
     homogeneous position 'P' and s is a 'Plane' then 'p|s==0' means
     p is a point on the plane s.
  """
  # local x,y,z Cartesian coordinate frame axis units; containment sense
  __slots__=['__u','__v','__w','__s']

  def __init__(self,*p):
    """Actual parameters may be:
        1)  3 homogeneous coordinate positions 'P' on the 'Plane'
        2)  a single homogeneous position 'P': the 'Plane' comprises all
            points orthogonal to 'P' translated by 'P'
        3)  a magnitude and differential D specifying normal direction.
            This is useful for specifying planes through the origin,
            where method 2) would yield no direction orthogonal.

       Internally, a plane is stored as a canonical 4-vector:
         [ -m a b c]
       where m is the orthogonal distance to the plane from the origin
       and [a b c] is the unit differential vector orthogonal to the
       plane that points toward the plane.

       The origin is arranged to be interior by default.  The unary
       positive operator '+' must be used to reorient a sheet with
       exterior origin, while '++' leaves orientation unchanged.

       For planes through the origin, [a b c] is one of the two unit
       orthogonals, the specific differential depending on the actual
       parameters that were passed.  Naturally, the origin is on the plane.
    """
    self.__u=self.__v=self.__w=None; self.__s=sign # Realize __slots__
    if not p: # class method '_' is caller: return empty instance
      super(Plane,self).__init__(); return # Realize superclass __slots__
    if len(p)==1: # given the orthogonal position
      p=p[0]; m=abs(p); w=(p-origin).unit
      x,y,z=w.p; u=(D(z,-x,y)^w).unit # guarantee non-zero normal
    elif len(p)==2: # given magnitude and direction
      m,w=p; w=w.unit
      x,y,z=w.p; u=(D(z,-x,y)^w).unit # guarantee non-zero normal
    elif len(p)==3: # calculate the orthogonal position
      if not(all(isinstance(e,P) for e in p)): raise TypeError,\
        "'Plane' instances are defined by type 'P' points."
      a,b,c=p; u,v=(a-b).unit,(b-c).unit; w=(u^v).unit; m=w|a
      if m<0.0: m,w = -m,-w
      if round(abs(w),15)==0.0: raise ValueError,\
        "Points in 'Plane' definition are collinear."
    else: raise ValueError,\
      "1 or 3 points define 'Plane' instances: %s received" % (len(p),)
    if w is None: raise ValueError, \
      "Zero length directional orthogonal for plane: %r" % (p,)
    v=(u^w).unit; self.__u,self.__v,self.__w=u,v,w
    super(Plane,self).__init__((-m,)+w.p)

  def __call__(self,p): return self.__s(round(p|self,15))

  def __contains__(self,p): return self(p) == -1

  def __pos__(self):
    """Flip the interpretation of interior and exterior points for this
       sheet.  Returns a reference to this sheet."""
    if self.__s==sign: self.__s=flip_sign
    else: self.__s=sign
    return self

  def __xor__(self,*s):
    """Return the intersection of surfaces '*s' and 'self':
         1) single plane       : return a 'Line'
         2) tuple of two planes: return a position 'P'
         3) sphere             : return s^self
       Returns 'self' for coincident planes.
       Returns a 'Line' for two unique planes and one redundant plane."""
    if type(s[0])==tuple: p0,p1=s[0]
    else:p0,p1 = s[0],self
    if type(p0)==Sphere: return p0^self
    if self.same(p0): p0=None
    if self.same(p1): p1=None
    if p0 is None: p0,p1=p1,p0
    if p0 is None: return self # coincident planes
    if p1 is None: # two unique planes => return a 'Line'
      d=(self.W^p0.W).unit; ref=Plane(0.0,d)
      v=P((SqM(origin,ref,self,p0).inverse|origin)[1:])
      return figure.Line(v,v+d)
    else: # three unique planes => return a position 'P'
      return P(round(e,15) for e in \
                (SqM(origin,p0,self,p1).inverse|origin)[1:])

  def same(self,p):
    """Return 'True' when plane 'p' is equivalent to self."""
    test=round(abs(self-p),15)
    return test==0.0 or test==2.0 and \
            (round(self[0],15)==0.0 and round(p[0],15)==0.0)

  @property
  def distance(self):
    """Return the orthogonal distance from the origin to the plane."""
    return -self[0]

  @property
  def W(self):
    """Return the unit orthogonal differential of the plane."""
    return self.__w

  @property
  def U(self):
    """Return the real axis unit differential in the plane."""
    return self.__u

  @property
  def V(self):
    """Return the imaginary axis unit differential in the plane."""
    return self.__v

  @property
  def origin(self):
    """Return the nearest position on the plane to the origin.
       This is the virtual origin used for planar figures."""
    return P(-self[0]*self.__w)

class Quadric(SqM):
  """Defines a homogeneous 4x4 matrix representation of the algebraic
     definition of a quadric surface: if p is a homogeneous position 'P'
     and s is a 'Quadric' then 'p|s|p==0' means p is a point on the
     quadric surface 's'.  Spheres, ellipsoids, paraboloids, and
     hyperboloids of one and two sheets are all quadric surfaces.

     The criterion for a point p being on s is equivalent to the
     algebraic quadratic form that defines the surface.  The matrix
     diagonal holds the coefficients for the constant, x^2, y^2, and z^2
     terms respectively.  The matrix is usually symmetric with off diagonal
     components representing half of a linear or cross term coefficient.

     For example s[3,4]+s[4,3] is the y*z coefficient.  Since linear and
     cross terms represent rotations and translations of a canonical
     quadric surface symmetrically positioned at the origin and aligned
     with the coordinate axes, we can simplify the definition by only
     providing a 4-tuple representing the diagonal.

     A more complex quadric with cross terms can then be created by
     multiplying the quadric by a transformation matrix representing
     a sequence of rotations, scalings, and translations.

     This process excludes forms that require an independent linear
     term: paraboloids and linear hyperboloids.  For those objects, an
     entire matrix must be provided.

     It is nevertheless possible to provide geometrically meaningless
     parameters, for example (1,1,1,1) which represents a sphere with
     imaginary radius.  This is preliminary code: not all pathological
     parameters are screened out.
  """
  __slots__=['__s']

  def __init__(self,*u):
    """Actual parameters may be:
        1) 4 numeric values
        2) A tuple or list of values
        3) A SqM
       Tuple values are [d a b c] where d is the constant term, and a,b,c
       are the coefficients of the x,y,z square terms respectively.

       For example:
        1) 'Quadric(-r*r,1,1,1)' creates a sphere with radius 'r'.
        2) 'Quadric(-1,1/(a*a),1/(b*b),1/(c*c))' creates an ellipsoid
           with x,y,z semi-axes of a,b,c respectively.
        3) 'Quadric(-1,1/(r*r),1/(r*r),-1)' creates a hyperboloid of
           one sheet with a radius 'r' waist in the x-y plane.
    """
    self.__s=sign # Realize __slots__
    if not u: # class method '_' is caller: return empty instance
      super(Quadric,self).__init__(); return # Realize superclass __slots__
    if len(u)==1: u=tuple(*u) # u is a tuple, list, or V or SqM
    if len(u)==4:
      if not isinstance(u[0],V): # SqM was NOT passed in
        u=(V(float(e),0.0,0.0,0.0)>>i for i,e in enumerate(u))
      super(Quadric,self).__init__(u)
    else: raise ValueError,\
      "A'Quadric' surface definition must have four components."

  def __call__(self,p): return self.__s(round(p|self|p,15))

  def __contains__(self,p): return self(p) == -1

  def __pos__(self):
    """Flip the interpretation of interior and exterior points for this
       sheet.  Returns a reference to this sheet."""
    if self.__s==sign: self.__s=flip_sign
    else: self.__s=sign
    return self

class Sphere(Quadric):
  """A quadric surface restricted to the special case of a sphere."""
  __slots__=['__c','__r']
  def __init__(self,c=origin,r=1.0):
    """Create a 'Sphere' of radius 'r' centered at 'c', a position 'P'."""
    self.__c,self.__r = P(c),float(r); d=-self.__c; x,y,z = d
    super(Sphere,self).__init__(
      map(V,[((d|d)-r*r,x,y,z),(x,1.,0.,0.),(y,0.,1.,0.),(z,0.,0.,1.)]))

  def __xor__(self,s):
    """Return the intersection of plane or sphere 's' and 'self':
         1) disjoint     : return 'None'
         2) Coincident   : return 'self'
         3) tangent      : return a 'P'
         3) intersecting : return a 'Circle'"""
    if type(s)==Sphere:
      p0,p1 = self.center,s.center; r0,r1 = self.radius,s.radius
      u=(p1-p0).unit; d,dmax = round(abs(u),15),r0+r1
      if d>dmax: return None # disjoint
      if d==0.0 and r0==r1: return self # coincident
      if d==dmax: return p0+r0*u # tangent
      l=figure.Line(p0,p1) # line through the two sphere's centers
      d0,d1 = r0*l.U,r1*l.U
      c0=figure.Circle(p0+d0,p0+r0*l.V,p0-d0) # coplanar great circles
      c1=figure.Circle(p1+d1,p1+r1*l.V,p1-d1)
      i0,i1 = c0^c1; v=i1-i0 # two points on the sphere intersection circle
      r=abs(v)/2.0; p3=i0+v/2.0 # radius,center of intersection
      return figure.Circle(i0,p3+r*l.W,i1)
    if type(s)==Plane:
      p0,w = self.center,s.W; r0,d = self.radius,s.distance
      d0=d-(+p0|w); r1=round(abs(d0),15)
      if r1>r0: return None # disjoint
      if r1==r0: return p0+r0*w # tangent
      c,r = p0+d0*w,sqrt(r0*r0-d0*d0); u=r*s.U
      return figure.Circle(c+u,c+r*s.V,c-u)
    raise ValueError, \
      "Intersection of 'Sphere' and {} is undefined.".format(repr(s))

  @property
  def radius(self):
    """Return the radius of this 'Sphere'."""
    return self.__r

  @property
  def center(self):
    """Return center of this 'Sphere' as a 'P'."""
    return self.__c

if __name__=='__main__':
  print "===Beginning module sheet.py test suite==="
  from math import sqrt
  r=sqrt(2); o=origin
  a=P(1,r/2,1,0); b=P(1,r/2,-1,0)
  c=P(1,-r/2,0,1); d=P(1,-r/2,0,-1)
  print "tetrahedron vertices"
  for p in (a,b,c,d): print repr(p)
  faces=[Plane(*f) for f in ((a,b,c),(a,b,d),(c,a,d),(b,c,d))]
  print "tetrahedron face plane representations"
  for f in faces: print f
  print "a tetrahedron edge as the intersection of two faces"
  print faces[0]^faces[1]
  print "a tetrahedron vertex as the intersection of three faces"
  print faces[0]^(faces[1],faces[2])
  print "tetrahedron face unit normals"
  for f in faces: print repr(f.W), abs(f.W)
  print "tetrahedron face origin containment"
  print tuple(o in f for f in faces)
  exteriors=[+f for f in faces]; s=faces[0]
  print "tetrahedron origin inverted containment"
  print tuple(o in f for f in exteriors)
  print "tetrahedron vertice etc, containment"
  for p in (a,b,c,d,o,a+(a-d)): print s(p),
  print
  print "tetrahedron face representation"
  for f in faces: print repr(f)
  sphere=Quadric(-4,1,1,1); hyperboloid=Quadric(-1,1,1,-1)
  ellipsoid=Quadric(-1,1,.25,.1)
  paraboloid=Quadric(SqM(V(-1,0,0,-1),V(0,1,0,0),V(0,0,1,0),V(-1,0,0,0)))
  print "quadric surface representations"
  print repr(sphere),repr(hyperboloid),repr(ellipsoid),repr(paraboloid)
  print "Origin containment:",\
    o in sphere, o in hyperboloid, o in ellipsoid, o in paraboloid
  x=P(3,0,0); y=P(0,3,0); z=P(0,0,3); w=P(0,0,10)
  print "x outside containment:",\
    x in sphere, x in hyperboloid, x in ellipsoid, x in paraboloid
  print "y outside containment:",\
    y in sphere, y in hyperboloid, y in ellipsoid, y in paraboloid
  print "z outside containment:",\
    z in sphere, z in hyperboloid, z in ellipsoid, z in paraboloid
  print "w outside containment:",\
    w in sphere, w in hyperboloid, w in ellipsoid, w in paraboloid
  s,p = Sphere(),Plane(.5,D(1.0,0.0,0.0))
  print "intersection of unit sphere and yz plane at x==.5"
  print p^s
