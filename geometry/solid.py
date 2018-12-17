#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Package: geometry
Module : solid.py
Website: https://github.com/ruminations/site-packages
License: https://github.com/ruminations/Licenses#design-license
Initial Copyright 2013

This module defines objects useful for dividing space into distinct
volumes, and determining containment and intersection with other objects.
"""
__all__=['Box','Closure']
__version__=20181216 # change imports to package imports
__version_history__=(20131107,)

from geometry.point import P,D
from geometry.sheet import Plane,sign,flip_sign

class Closure(object):
  """A 'Closure' is a container class of 2-dimensional surface classes
     in Euclidean 3-space.  It divides space into interior and exterior
     partitions.  A 'Closure' thus defines a geometrically more complex
     partition than its component surfaces, and may be a component of
     another 'Closure' or may contain another 'Closure' as a component.

     One or both partitions defined by a 'Closure' may be infinite in
     extent.  A plane is perhaps the simplist component, which by itself
     defines two infinite partitions.

     A 'Closure' uses its surfaces to create a well defined separation
     of an interior volume and exterior volume via boolean composition
     of the component containments.

     The property 'union' specifies the way the partitions will be
     combined to form the closure.  The default is intersection.
     Property 'union' values are:
       'False' :  Set intersection - a position is in the closure only
                  if it is inside all component partitions.
       'True'  :  Set union - a position is inside the closure if it is
                  in any of the component partitions.

     So, for example, if you wanted to place a round chamfer on the
     edge of a box (form a set difference), you could do this:
       c=Closure(+cyl,p1,p2,+p3); b=Closure(box,+c)
     where 'box' is an instance of 'Box', 'cyl' is a sheet defining a
     cylinder, 'p1' and 'p2' are planes from the box that are tangent
     to the cylinder, and 'p3' joins the lines of tangency.

     The first closure defines an interior point as one outside of the
     cylinder and inside the prism defined by the three planes.  The
     second closure makes sure a point is both inside the box and
     outside of the first closure.
  """
  __slots__=['__sheets','__union','__s']

  def __init__(self,*s):
    """'*s' specifies class instances representing a partition of space.
       The instance objects may be sheets or other closures.
    """
    self.__sheets=tuple(*s); self.__union=False; self.__s=sign

  def __call__(self,p):
    """Return, with respect to orientation convention in force:
         -1 : for an interior point
          0 : for a surface point
          1 : for an exterior point"""
    tests=[s(p) for s in self.__sheets]
    if self.union:
      if any(r==-1 for r in tests): return self.__s(-1)
      if any(r==0 for r in tests): return 0
    else: # intersection
      if all(r==-1 for r in tests): return self.__s(-1)
      tests=[r for r in tests if r!=0]
      if all(r==-1 for r in tests): return 0
    return self.__s(1)

  def __contains__(self,p):
    """Returns, with respect to orientation convention in force:
         'True'  : 'p' is a position 'P' strictly inside the closure.
         'False' : 'p' is a position 'P' on or outside the closure.
    """
    return self(p)==-1

  def __pos__(self):
    """Flip the interpretation of interior and exterior points for this
       'Closure'.  Returns a reference to this 'Closure'."""
    if self.__s==sign: self.__s=flip_sign
    else: self.__s=sign
    return self

  @property
  def union(self):
    """Return 'True' if partitions combine via set union,
       'False' if they combine via set intersection."""
    return self.__union

  @union.setter
  def union(self,value):
    if isinstance(value,bool): self.__union=value
    else: raise TypeError, \
      "Property 'union' must be boolean: %s received." % type(value)
  '''
  @property
  def extent(self):
    """Return a reference to the closure's bounding box."""
    return self.__surf.extent

  @property
  def position(self):
    """Return the vector offset of the closure's location."""
    return self.__surf.position
  '''
class Box(Closure):
  """'Box' is defined by the corner points of a rectangular prism in
     Euclidean 3-space aligned with the coordinate axes."""
  __slots__=['__llf','__urr','__closure']
  def __init__(self,lower_left_front,upper_right_rear):
    """Specify the box by specifying positions of type 'P' defining
       opposite points on a diagonal.  Specific orientation is moot."""
    self.__llf=lower_left_front; self.__urr=upper_right_rear; planes=[]
    for i,t in enumerate((self.x_span,self.y_span,self.z_span)):
      vmin,vmax=t
      planes.append(Plane(abs(vmin),D(vmin-vmax,0,0)>>i))
      planes.append(Plane(abs(vmax),D(vmax-vmin,0,0)>>i))
    super(Box,self).__init__(planes)

  @property
  def x_span(self):
    """Return the tuple (x_min,x_max)."""
    return self.x_min,self.x_max

  @property
  def x_min(self):
    """Return the least value of the x coordinate."""
    return min(self.__urr[1],self.__llf[1])

  @property
  def x_max(self):
    """Return the greatest value of the x coordinate."""
    return max(self.__urr[1],self.__llf[1])

  @property
  def y_span(self):
    """Return the tuple (y_min,y_max)."""
    return self.y_min,self.y_max

  @property
  def y_min(self):
    """Return the least value of the y coordinate."""
    return min(self.__urr[2],self.__llf[2])

  @property
  def y_max(self):
    """Return the greatest value of the y coordinate."""
    return max(self.__urr[2],self.__llf[2])

  @property
  def z_span(self):
    """Return the tuple (z_min,z_max)."""
    return self.z_min,self.z_max

  @property
  def z_min(self):
    """Return the least value of the z coordinate."""
    return min(self.__urr[3],self.__llf[3])

  @property
  def z_max(self):
    """Return the greatest value of the z coordinate."""
    return max(self.__urr[3],self.__llf[3])

  @property
  def extent(self):
    """Return a reference to the closure's bounding box."""
    return self

class Tetrahedron(Closure):
  """'Tetrahedron' is defined by the four vertices in Euclidean 3-space."""
  __slots__=['__vertices','__closure','__box','__volume']
  def __init__(self,p0,p1,p2,p3):
    """Specify the tetrahedron by specifying positions of type 'P'
       defining the vertices."""
    self.__vertices=(p0,p1,p2,p3); planes=[]
    planes.append(Plane(p0,p1,p2))
    planes.append(Plane(p0,p1,p3))
    planes.append(Plane(p0,p2,p3))
    planes.append(Plane(p1,p2,p3))
    super(Tetrahedron,self).__init__(planes)
    self.__box=Box()
    self.__volume=None

  @property
  def centroid(self):
    """Return a reference to the closure's bounding box."""
    return self.__centroid

  @property
  def volume(self):
    """Return a reference to the closure's bounding box."""
    return self.__volume

  @property
  def extent(self):
    """Return a reference to the closure's bounding box."""
    return self.__box

if __name__=='__main__':
  print "===Beginning module solid.py test suite==="
  from point import origin
  p,q=P(1,1,1),P(-1,-1,-1)
  b=Box(p,q)
  print "origin in box", origin in b, b(origin)
  print "corner in box", p in b, b(p)
  print "corner in box", q in b, b(q)
  r=P(2,2,2)
  print "Outlier", r in b, b(r)
  r=P(1,1,2)
  print "Outlier in edge", r in b, b(r)
  r=P(1,0,2)
  print "Outlier in surface", r in b, b(r)
