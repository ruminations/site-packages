#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Package: geometry
Module : figure.py
Website: https://github.com/ruminations/site-packages
License: https://github.com/ruminations/Licenses#design-license
Initial Copyright 2013

Module 'figure' defines two dimensional planar figures embedable in a
homogeneous projective space of three dimensions.  The module is tested
only for three dimensional homogeneous positions and differentials,
but the ideas apply for suitably defined higher dimensional spaces.

The key for higher dimensions is a meaningful cross-product, implemented
via a general wedge product.  This entails formalizing differential
geometry's tangent, normal and binormal of space curves for the degenerate
case of a line, using the generalized Frenet equations.  That is beyond
the scope of the effort here.

There is still an issue with circular imports.
If you import geometry.sheet first, then importing geometry.figure works.
"""
__all__=['Line','Edge','Triangle','Circle','Arc','Frame']
__version__=20181216 # change imports to package imports
__version_history__=(20131113,20140602,20150929,20171213,20180411)
### TODO:
# 0) Edge/Line containment testing
# 1) Fix Triangle.barycenter to work with neg coordinates; trilinear coords
# 2) Continue refactor to use complex numbers for planar calculations - accuracy
# 3) check and fix various ## marks
from math import sqrt,acos
from algebra.mathx import len,odd,sign,delta_ratio_is_zero
from numbers import Real
from geometry.gauss import C,Angle,Cline
from geometry.point import P,D,Q,origin
try: from geometry.sheet import Plane
except ImportError: pass # skip circular import second pass

class Frame(object):
  """A reference frame translator: reference a local complex (C objects)
     plane to three given spatial positions.  Facilitate translating
     arbitrary in plane spatial 'P' objects to equivalent C objects and
     translating C objects to equivalent 'P' objects."""
  __slots__=['_u','_v','_defs','_p0','_dp1','_dp2']
  def __init__(self,p0,p1,p2):
    """'pi' are non-collinear positions 'P'.  'p0' corresponds to the
       complex plane origin.  'p1' is in the direction of the real axis.
       'p2' is used to find the orthogonal imaginary axis."""
    self._p0,self._dp1,self._dp2 = p0,p1-p0,p2-p0
    self._u=self._dp1.unit; self._v=self._dp2.unit # v as a temporary
    c=(self._u|self._v); self._v=(self._v-c*self._u).unit # u,v now orthogonal, in plane
    self._defs=tuple( (p0, p1, p2, C(0.0), C(abs(self._dp1)), \
                     C(abs(self._dp2)*c,(self._dp2|self._v))) )

  def __call__(self,zp):
    """If 'zp' is a 'C' object, return the equivalent 'P' object.
       If 'zp' is a 'P' object, return the equivalent 'C' object."""
    if isinstance(zp,C): return self._p0 + self._u*zp.x + self._v*zp.y
    if isinstance(zp,P):
      delta=zp-self._p0; return C(delta|self._u,delta|self._v)
    raise ValueError,"{!r} must be a 'C' or 'P' object.".format(zp)

  def __repr__(self): return "_Frame({!r},{!r},{!r})".format( \
                     self._p0,self._p0+self._dp1,self._p0+self._dp2)

  def __getitem__(self,i):
    """Get the 'i'-th defining point:
         0 <= i <= 2   as a 'P' value
         3 <= i <= 5   as a 'C' value"""
    return self._defs[i]

  def __getslice__(self,*args): return self._defs[slice(*args)]

  @property
  def U(self):
    """Get the unit 'D' defining the real axis."""
    return self._u

  @property
  def V(self):
    """Get the unit 'D' defining the imaginary axis."""
    return self._v


class Line(object):
  """Represents an infinite line in three dimensions."""
  __slots__=['_p','_u','_n','_args']
  # plane containing unit segment, unit normal
  # binormal = (u^n).unit points from origin to plane.origin
  # _u,_n,binormal form an orthonormal Cartesian basis
  def __init__(self,p,dp):
    """'p' and 'dp' are distinct positions 'P' on the line, or
       'p' is a position 'P' and dp is a differential 'D' from
       'p' to another point on the line."""
    self._args=(p,dp)
    if isinstance(dp,D): u,p1,p0 = dp.unit,p+dp,p
    else: u,p1,p0 = (dp-p).unit,dp,p
    if u is None: raise ValueError, \
      "Indistinct points defining line: %s,%s" % (p0,p1)
    d=p0-origin; du=d.unit
    if du is None or round(abs(du)-abs(du|u),15)==0.0: # 0.0 <= inner product <= 1.0
      # 'origin' on the line: the binormal is arbitrary: guarantee normal component
      if any(e==0.0 for e in u.p): # planar: just rotate
        i=0;t=u
        while t[-1]!=0.0: i,t = i+1,t<<1
        b=D(-t[2],t[1],0.0)>>i
      else:
        b=D(e if odd(i) else -e for i,e in enumerate(u<<1)).unit # b|u != 0.0
        b=(b-(b|u)*u).unit # orthogonalize
    else:
      b=(d-(d|u)*u).unit # binormal to the line at point nearest to origin
    self._p=Plane((p0|b),b)
    self._u=u; self._n=(b^u).unit # normal to the line in plane

  def __contains__(self,p):
    """Return 'True' for a point on the 'Line'."""
    print "\n in contains",self._p(p),(p|self._p)
    if self._p(p) != 0: return False # not in containing plane
    print "in plane"
    return round((p-self._p.origin|self._n),15)==0.0

  def __repr__(self):
    """Return a string showing internal representation."""
    return self.__class__.__name__+"({!r},{!r})".format( \
                        self.origin.clean, (self.origin+self._u).clean )

  def __str__(self):
    """Return a string showing intial arguments."""
    return self.__class__.__name__+"({!r},{!r})".format(*self._args)

  def __xor__(self,l):
    """Return the minimal line segment separating 'self' and line 'l'.

       For coplanar lines, return intersection  position 'P'(possibly in
       the plane at infinity).  Raises a 'ValueError' for coincident lines.

       If 'isinstance(l,Circle)', return '(l|self)' - the intersection(s)
       of this line with that 'Circle'."""
    if isinstance(l,Circle): return l^self
    p0=P(self._p.origin); u0=self._u; p1=P(l._p.origin); u1=l._u
    plane_n=u1^u0; s=abs(plane_n) # mutually normal plane, sin(angle)
    n=plane_n.unit
    if n is None: # parallel => sin == 0.0
      if round(abs(p0-p1),15)==0.0: # coincident
        raise ValueError, \
              "Attempt to intersect coincident lines: %r, %r" % (self,l)
      return P(0.0,u0.p) # point at infinity on the lines
    r=((p1-p0)^n)/s # directed law of sines: length/sin (complement projection)
    b=p0+u0*(r|u1); e=p1+u1*(r|u0) # law of sines (complement projection)
    if delta_ratio_is_zero(e,b,14): return b # intersecting
    return Edge(b,e)

  @property
  def distance(self):
    """Return the line's distance from the origin."""
    return round(self._p.distance,15)

  @property
  def plane(self):
    """Return the line's plane."""
    return self._p

  @property
  def origin(self):
    """Return the nearest position on the line to the coordinate 'origin'.
       This is the virtual origin used for linear figures in the plane
       orthogonal to 'self.W'."""
    return self._p.origin

  @property
  def U(self):
    """Return the line's unit differential."""
    return self._u

  @property
  def V(self):
    """Return the line's in plane unit normal differential."""
    return self._n

  @property
  def W(self):
    """Return the unit binormal differential orthogonal to the defining
       plane containing the line.  The binormal points from the coordinate
       'origin' to 'self.origin'"""
    return self._p.W

class Edge(Line):
  """Represents a finite segment on a line in three dimensions."""
  __slots__=['_b','_e','_l']
  def __init__(self,p0,p1):
    """'p0' and 'p1' are line segment endpoint spatial positions."""
    self._b,self._e = p0,p1; self._l=abs(p1-p0)
    super(Edge,self).__init__(p0,p1)

  def __contains__(self,p):
    """Return 'True' for positions on the line and
       strictly interior to the segment."""
    if super(Edge,self).__contains__(p):
      return 0 < abs(p-self._b) < len(self)
    return False

  def __len__(self):
    """Return the length of this line segment."""
    return self._l

  def __repr__(self):
    return self.__class__.__name__+"({!r},{!r})".format(*self.vertices)

  @property
  def vertices(self):
    """Return segment's end points as a tuple: '(p0,p1)'."""
    return (self._b,self._e)

class Circle(object):
  """Represents a planar circle in three dimensions."""
  __slots__=['_c','_r','_p','_u','_v','_f','_ref','_args']
  #        center,radius,plane,uv Ds, distance, init args
  def __init__(self,*p):
    """Initialization arguments may be in any order as:
         r,c,p     scalar radius 'r', a 'P' center 'c', and a 'Plane' 'p'
         p0,p1,p2  three non-collinear positions 'P' on the circle
         d0,p0,c
         p0,d0,p1  correllated tangent 'D' and positions 'P' for two
         p0,p1,d1  positions  on the circle or a position and a center.
                   The specification is for an arc from 'p0' to 'p1', with
                   the tangent differential to its position."""
    a,b,c = self._args = p
    if any(isinstance(v,Plane) for v in p): # radius, center, plane
      if isinstance(c,Plane): a,c = c,a
      if isinstance(b,Plane): a,b = b,a
      if isinstance(b,P): b,c = c,b # a plane, b radius, c center
      assert isinstance(c,P) and isinstance(b,Real), \
             "Invalid circle parameters {!r}, {!r}, {!r}.".format(a,b,c)
      self._p,self._r,self._c = a,float(b),c
      self._u,self._v   = self._p.U,self._p.V
      self._f,self._ref = Frame(self._c,self._u,self._v),self._c+self._u*self._r
      return
    d=d0=d1=None
    if any(isinstance(v,D) for v in (a,b,c)): # differential and positions
      if isinstance(a,D): d,p0,p1 = a,b,c
      if isinstance(b,D): p0,d0,p1 = a,b,c
      else:               p0,p1,d1 = a,b,c
      assert isinstance(p0,P) and isinstance(p1,P), \
             "Invalid circle parameters {!r}, {!r}, {!r}.".format(a,b,c)
      if d is not None:
        self._c,u = p1,p0-p1; self._r=abs(u)
        d=(d-(d|u)*u).unit; p1=self._c+d*self._r
      else:
        if d0 is None: # tangent applies to p1
          d=d1.unit; self._f=Frame(p0,p1,p1+d) # p1+d NOT on circle
          b1,b2 = self._f[-2]/2.0,p1 ##; u1,n1 = b1.unit,b1.ortho
          n1,n2 = b1.ortho,self._f(p0+d)*C(1j) # tangent as a 'C': p0 is origin
          ##n2=((d|u1)*n1-(d|n1)*u1).unit # rotate d 90 deg
        else:          # tangent applies to p0
          d=d0.unit; self._f=Frame(p0,p0+d,p1) # p0+d NOT on circle
          b2,b1 = self._f[-1]/2.0,p0 ##; u2,n2 = b2.unit,b2.ortho
          n2,n1 = b2.ortho,self._f(p0+d)*C(1j) # tangent as a 'C': p0 is origin
          ##n1=((d|u2)*n2-(d|n2)*u2).unit # rotate d 90 deg
        self._c=self._f(Cline(b1,b1+n1)^Cline(b2,b2+n2)); self._r=abs(p0-self._c)
      ##b=(p1-p0)/2.0; u=b.unit; cs=(d|u) ### reference until tested
      ##n1=(d-ip*u).unit; n2=(ip*n1-(d|n1)*u).unit
      ##if round(ip,14)==0.0: self._c=p0+b # differential orthogonal to p1-p0
      ##else: self.__c=Line(p0+b,n1)^Line(p1,n2)
      self._p,self._f = Plane(self._c,p0,p1),Frame(self._c,p0,p1) # reframe to center
      self._u,self._v = self._f.U,self._f.V; self._ref=p0
      return
    assert all((isinstance(a,P),isinstance(b,P),isinstance(c,P))), \
           "Invalid circle parameters {}, {}, {}.".format(a,b,c)
    self._f=Frame(a,b,c);
    try: self._p=Plane(a,b,c)
    except TypeError: raise ValueError, \
      "Collinear points defining circle: {}, {}, {}".format(a,b,c)
    c1,c2=self._f[-2:]; b1,b2 = c1/2.0,c2/2.0
    self._c=self._f(Cline(b1,b1+b1.ortho)^Cline(b2,b2+b2.ortho))
    self._r=abs(a-self._c); self._f=Frame(self._c,a,c) # reframe to center
    self._u,self._v = self._f.U,self._f.V; self._ref=a

  def __call__(self,p):
    """Return 'None' for non-coplanar 'p'.
       Return:
        -1 : 'p' inside the circle
         0 : 'p' on the circle
         1 : 'p' outside the circle"""
    if round(self._p(p),15)!=0.0: return None
    return sign(round(abs(p-self._c)-self._r,15))

  def __contains__(self,p):
    """Return 'True' for a point on the circle."""
    return self(p)==0

  def __repr__(self):
    """Return a string showing internal representation."""
    d,p,c = self.tangent(Angle(0.0)),self._ref,self._c
    return self.__class__.__name__+"({!r},{!r},{!r})".format( \
                              d.clean,self._ref.clean,c.clean )

  def __str__(self):
    """Return a string showing initial arguments."""
    return self.__class__.__name__+"({!r},{!r},{!r})".format(*self._args)

  def __xor__(self,f): ## refer to LyX doc 'Circle_Intersection.lyx'
    """Compute the intersection of this circle with figure 'f'.
       Supported figures: 'Circle' and 'Line'"""
    sc,sr,sp = self.center,self.radius,self.plane
    if isinstance(f,Line):
      fu,fv,fo = f.U,f.W,f.origin # line unit, orthogonal, reference point
      if sp(fo) == 0 and sp(fo+fu) == 0: # coplanar: fo==sp.origin==self.origin
        do=(sc-fo); xv=(do|fu)*fu; yv=do-xv; d=round(abs(yv),15)
        if d > sr: return None # disjoint
        if d == sr: return ((sc-yv).clean,) # tangent
        dxv=sqrt(sr*sr-d*d)*fu
        return (fo+xv+dxv,fo+xv-dxv) # intersecting
      else: # intersects plane in a point
        c=(sp.W|fu) # cosine
        if round(c,15) == 0.0: return None # point at infinity
        d=sp.distance-(sp.W|(fo-origin)); pint=fo+fu*d/c
        if round(abs(sc-pint)-sr,15)==0.0: return (pint,)
        return None
    if isinstance(f,Circle):
      fc,fr,fp = f.center,f.radius,f.plane
      rs=sr+fr; rd=fr-sr; u=fc-sc; d,u=abs(u),u.unit
      if sp.same(fp): # coplanar
        if d > rs: return None # separated
        if d < abs(rd): return None # containment
        if rd == 0.0 and u is None: return None # coincident
        if d == rs: return (sc+u*sr,) # separate tangent
        if d == abs(rd): return (sc+sr*(-u if fr>sr else u),) # containment tangent
        v=u^sp.W; g2=rs*rd; m=(g2/d-d)/2.0 # g2 is signed due to rd !
        x,y = (d-g2/d)/2.0,sqrt(sr*sr-m*m); xv,yv=sc+u*x,v*y
        return (xv+yv,xv-yv) # intersection pair
      else:
        l=(sp^fp); si,fi = self^l,f^l
        if si is None or fi is None: return None # disjoint
        pts=tuple(a for a in si for b in fi if round(abs(a-b),12) == 0.0)
        if pts == (): return None # disjoint or linked
        return pts # length 1 =>pivoted length 2 => hinged
    raise ValueError,"Type({!r}) == {} is not a 'Line' or a 'Circle'".format(f,type(f))

  def tangent(self,angle=Angle(0.0)):
    """Return the unit tangent at 'Angle' object 'angle' from the reference
       position self.vertex in the direction of 'p1' used to create the
       circle."""
    return self._f(C(angle+Angle(.5)))-self.center

  @property
  def center(self):
    """Return the circle's center position."""
    return self._c

  @property
  def origin(self):
    """Return the nearest position in 'self.plane', orthogonal to 'self.W',
       to the coordinate 'origin'.  This is the virtual origin used for
       operations on circles."""
    return self._p.origin

  @property
  def radius(self):
    """Return the circle's radius."""
    return self._r

  @property
  def plane(self):
    """Return the circle's plane."""
    return self._p

  @property
  def vertex(self):
    """Get the reference vertex on the circle.  This is 'p0' used to
       create the circle, or P(self._center+self.radius*self.U,0.0,0.0)
       if created with a plane and a center."""
    return self._ref

  @property
  def U(self):
    """Get the circle's in plane unit reference differential.
       This differential points from 'self.center' to 'self.vertex'."""
    return self._u

  @property
  def V(self):
    """Get the circle's in plane unit normal differential."""
    return self._v

  @property
  def W(self):
    """Get the unit binormal differential orthogonal to the defining
       plane containing the circle."""
    return self._p.W

class Arc(Circle):
  """Represents a planar Arc three dimensions."""
  __slots__=['_a']

  def __init__(self,*p):
    """'p0', 'p1', and 'p2' are non-collinear spatial positions.  The
       arc extends from 'p0' through 'p1' to 'p2'."""
    if not all(isinstance(e,P) for e in p) or len(p)!=3: raise TypeError, \
      "'Arc' is defined by three 'geometry.point.P' positions."
    super(Arc,self).__init__(*p); self._a=None

  def __call__(self,alpha):
    """Interpolate a position on the arc:
          'alpha'     return
        --------------------
           -1.0       'p2'
       -1.0> & >0.0   position on circular arc excluding 'p1',
                      linear in anglular deviation.
            0.0       'p0'
        0.0< & <1.0   position on circular arc containing 'p1',
                      linear in angular deviation.
            1.0       'p2'
       Values outside the range are converted to the corresponding
       cyclic modulus inside the range."""
    a=Angle(float(alpha)) # sanitize
    c,da = self.center,t*self._a/2.0; q=Q(cos(da.rad),sin(da.rad)*self.W).unit
    return c+q(self._args[0]-c)

  def __contains__(self,p):
    """Return 'True' for a point on the arc."""
    if super(Arc,self).__call__(p) != 0.0: return None
    return ((p-self.center)|self.V) >= self.chord

  @property
  def angle(self):
    """Get an 'Angle' object representing the angle subtended by the 'Arc'."""
    if self._a is not None: return self._a
    p0,p1,p2 = self._args; c=self.center
    self._a=Angle(acos(p0-c|p2-c),'rad')


class Triangle(Frame):
  """Represents a planar triangle in three dimensions."""
  __slots__=['_v','_e','_centroid','_area','_circle','_edges']
  def __init__(self,p0,p1,p2):
    """'p0', 'p1', and 'p2' are non-collinear spatial positions 'P' for
       the triangle's vertices."""
    self._v=(p0,p1,p2)
    if not all(isinstance(p,P) for p in self._v): raise TypeError, \
      "Triangle is defined by three positions of type 'P'."
    super(Triangle,self).__init__(self,p0,p1,p2)
    self._edges=tuple(map(Edge,zip(v,v<<1)))
    self._e=(p1-p0,p2-p1,p0-p2) # sum(self._e)==origin
    self._area=self._heron(self._e)
    self._centroid=P(sum(~p for p in self._v)/3.0) # mean of vertices
    self._circle=Circle(p0,p1,p2)

  def __call__(self,p): ## will not work until barycenter is generalized
    """Return 'None' for non-coplanar 'p'.
       Return:
        -1 : 'p' inside the triangle
         0 : 'p' on the triangle
         1 : 'p' outside the triangle"""
    if not p in self.__circle.plane: return None
    b=~self.barycentrer(p) # normalize by folding the homogeneous coordinate
    if all((e>=0.0 and e<=1.0) for e in b):
      if any((e==0.0) for e in b): return 0
      else: return -1
    return 1

  def __contains__(self,p): return self(p)==-1

  def __repr__(self): return "Triangle({!r},{!r},{!r})".format(*self._v)

  def __xor__(self,l):
    """Return the intersection position(s) of 'self' with Line 'l'.  This
       can be any position in the triangle's plane for 'l' not in the
       triangle's plane.  For 'l' parallel, this will be a point at infinity.

       For 'l' in the triangle's plane, return a tuple of the intersections
       with the edges.  The tuple can have length zero, one or two.  If the
       line coincides with an edge, return 'None'."""
    if not self._circle.plane.same(l.plane): return (self._circle.plane^l.plane)^l
    try: ints=tuple(l^e for e in self._edges)
    except ValueError: return None # conincident with an edge
    return tuple(i for i in ints if any(i in e for e in self._edges))

  def _heron(self,*d): ## need to generalize to include neg areas
    """Compute the area of a triangle via Heron's formula for a triple
       of 'D' differentials representing the sides of the triangle."""
    sq=D(e|e for e in d); sq_s=sum(sq); q_s=sum(sq|sq)
    return sqrt(sq_s*sq_s-2*q_s)/4.0

  def barycenter(self,p): ### need to generalize to neg coordinates
    """Return a 'P' object whose components are the homogeneous planar
       barycentric coordinates of position 'p' in this triangle.

       The redundant homogeneous coordinate contains the area of this
       triangle, so the normalized barycenter"""
    if not isinstance(p,P): raise TypeError,"Expected a 'P' object."
    if not p in self._circle.plane: raise ValueError, \
      "Point {!r} is not in the same plane as this triangle".format(p)
    trgls=tuple(tuple(b-a,p-b,a-p) for a,b in zip(self._v,self.v<<1))
    return P(self._area,tuple(self._heron(t) for t in trgls))

  @property
  def area(self):
    """Return the triangle's area via using Heron's formula"""
    return self._area

  @property
  def centroid(self):
    """Return the triangle's centroid position 'G'."""
    return self._centroid

  @property
  def circumcenter(self):
    """Return the triangle's circumcenter position 'C'."""
    return self._circle.center

if __name__=='__main__':
  z0=D(0,2,0)
  z1=D(2,0,0)
  z2=D(2,2,0)
  p=[P(z) for z in [z0,z1,z2]]
  edge=Edge(p[0],p[1])
  l0=Line(p[2],z2)   # y=x
  l1=Line(p[1],p[0]) # y=-x+1
  c0=Circle(*p)
  p=[P(-z+z2*(1-1/sqrt(2))) for z in [z0,z1,z2]]
  c1=Circle(*p)
  print "\n\n=== Line tests: ===\n"
  print l0,'\n',l1
  print repr(l0),'\n',repr(l1)
  try:
    l=Line(+z2,+z2); print "Line(p0,p0) error raising failed"
  except ValueError as e: print "successfully raised error:\n  ",e
  try: print l0^l0, "Coincident line intersection error raising failed"
  except ValueError as e: print "successfully raised error:\n  ",e
  print "\nl0,l1 intersection == P(1.0,(1.0, 1.0, 0.0)):",repr((l0^l1).clean)
  print "\n\n=== Edge tests: ===\n"
  print +z2/2.0,repr(+z2/2.0)
  print edge,len(edge),repr(edge._p),+z2/2.0 in edge
  print "\n\n=== Circle tests: ===\n"
  print c0,'\n',c1
  import inspect
  d=dict((k,v) for k,v in inspect.getmembers(\
          Circle, predicate=inspect.isdatadescriptor)\
          if not k.startswith('_'))
  for k,v in sorted(d.iteritems()):
    print "  %-20r : %r" % (k,v.__get__(c0))
    print "  %-20r : %r" % (k,v.__get__(c1))
  print c0.vertex in c0, c0(c0.vertex)
  print c0^c1
  print c1^c0
  print
  l=Line(origin,P(0.0,0.0,1.0))
  print l,c0^l # single point, orthogonal intersection
  l=Line(P(0.0,2.0,0.0),D(1.0,1.0,0.0))
  print l,c0^l # single point, coplanar tangent intersection
  print "\n\n=== Practical tests: ===\n"
  l=Line(P(1.0, 389.0, 0.0, 0.0),P(1.0, 389.0, 1.0, 0.0))
  c=Circle( P(1.0, 390.0, 0.00295241992118,  0.0), \
          P(1.0, 391.0, 0.00357727459669, 0.0), \
          P(1.0, 392.0, 0.00433214623759, 0.0)  )
  #l=Line(P(1.0, 831.0, 0.0, 0.0),P(1.0, 831.0, 1.0, 0.0))
  #c=Circle( P(1.0, 828.0, 1.77308261662e-06,  0.0), \
  #        P(1.0, 829.0, 1.67308571253e-06, 0.0), \
  #        P(1.0, 830.0, 1.57919944105e-06, 0.0)  )
  print l,'\n',repr(l)
  print c,'\n',repr(c)
  print l^c
  print
