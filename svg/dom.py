#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Package: svg
Module : dom.py
Website: https://github.com/ruminations/site-packages
License: https://github.com/ruminations/Licenses#design-license
Initial Copyright 2015

Provides a direct low level interface to generate SVG text files
conforming to the SVG specification.
"""
###TODO:
# 0) as always, writing a doctest suite.
# 1) look into using a meta-class for key customization
from __future__ import unicode_literals
from math import pi,sqrt,cos,sin,tan,atan2
from geometry.gauss import C,Angle,O
__all__=['Key','SVG', 'Group', 'Transform', 'Defs', 'Use', 'ClipPath',\
         'Circle', 'Ellipse', 'Rectangle', 'Path', 'Text']
__version__=20180609
# added 'xform' property,tuple slots,'None' nop for __add__; span property
__version_history__=(20150430,20150617,20151027,20160105,20170315,20171226)

#---------------+---------------+---------------+---------------#
#---------------+---------- Utilities ----------+---------------#
#---------------+---------------+---------------+---------------#

class _(float):
  """Minimally print precision controlled floats.
      _(1.00000)  => 1.0
      _(1.000012) => 1.00001
      _(1.000015) => 1.00002
      _(1.1)      => 1.10
  """
  precision=5
  def __str__(self):
    s="{{:.{0}f}}".format(self.precision)
    s=s.format(self)
    if s.endswith('0'): s=s.rstrip('0')+'0'
    return s

class Matrix(object):
  """Rudimentary 2D homogeneous matrix functionality."""
  from operator import __mul__ as _s_mul
  __slots__=['_m','_inv','_inv_ok']

  def __init__(self,row0=(1.0,0.0,0.0),row1=(0.0,1.0,0.0),row2=(0.0,0.0,1.0)):
    """A 2D homogeneous matrix is 3x3 - 3 rows of 3, identity by default. """
    self._m=(tuple(row0),tuple(row1),tuple(row2)); self._inv=None; self._inv_ok=False

  def __call__(self,z):
    """Apply self to the complex number 'z', folding in the homogeneous scaling."""
    x,y = C(z); v_mtrc=self._v_metric # speedup: local lookup
    t=tuple(v_mtrc(r,(x,y,1.0)) for r in self._m)
    return C(t[0]/t[2],t[1]/t[2])

  def __invert__(self):
    """Return the transpose of 'self'."""
    return Matrix(*zip(*self._m))

  def __getitem__(self,i):
    """Return the ith row as a tuple."""
    return self._m[i]

  def __mul__(self,n):
    """Return self*n for compatible matrices 'm' and 'n'."""
    n=~n; v_mtrc=self._v_metric # speedup: local lookup
    return Matrix(*[[v_mtrc(r,c) for c in n._m] for r in self._m])

  def __rmul__(self,n): return n*self

  def __repr__(self): return "Matrix({0},{1},{2})".format(*self._m)

  @staticmethod   # private vector function
  def _v_metric(u,v):
    """Return the inner product of commensurate vectors 'u' and 'v'.
      '_v_metric' is the Euclidean metric tensor, hence the name."""
    return sum(map(Matrix._s_mul,u,v))

  @property
  def inverse(self): # algorithm documented in 'matrix_inverse.py'
    """Return the inverse of this 'Matrix'; 'None' if singular."""
    if self._inv_ok: return self._inv
    m=[e for t in self._m  for e in t]
    ct02 = m[4]*m[8]-m[5]*m[7], m[2]*m[7]-m[1]*m[8], m[1]*m[5]-m[2]*m[4]
    ct35 = m[5]*m[6]-m[3]*m[8], m[0]*m[8]-m[2]*m[6], m[2]*m[3]-m[0]*m[5]
    ct68 = m[3]*m[7]-m[4]*m[6], m[1]*m[6]-m[0]*m[7], m[0]*m[4]-m[1]*m[3]
    try:
      d = 1.0/(m[0]*ct02[0]+m[1]*ct35[0]+m[2]*ct68[0])
      self._inv=Matrix(
        *( map(Matrix._s_mul,(d,d,d),v) for v in (ct02,ct35,ct68) )   )
    except ZeroDivisionError: self._inv=None
    finally: self._inv_ok=True
    return self._inv


class Key(object):
  """Infinite iterable to generate sequential SVG id strings of the form
       'id="<name>_<#>"' e.g. 'id="group_0"', 'id="group_1"' etcetera.
     DOM definition instances obtain an identity by calling 'next(key)'.

     The default id string may be overridden on an element by supplying
     an instance of this class for the 'key' keyword argument.  This
     is useful for labeling custom sub-classes of standard elements."""
  __slots__=['_nm','_id']
  def __init__(self,name): self._nm,self._id = name,-1
  def __iter__(self): return self
  def next(self): self._id+=1; return 'id="{}_{}"'.format(self._nm,self._id)

#---------------+---------------+---------------+---------------#
#---------------+------- SVG DOM classes -------+---------------#
#---------------+---------------+---------------+---------------#


### NOTE TO SELF: I'd like to add a facility for xml comment and the
# <metadata> tag.  Both these are not 'Element' objects as defined here:
# they have no transformation.
#
# A comment could appear anywhere, but there is really only one that is
# important for an automatically generated file: a header comment
# that follows the entity declarations.  This should name the software
# and version used to generate the file, and might also reference
# a license.
#
# A <metadata> tag really only belongs in one place - after the svg
# preamble but before the actual geometry definitions or descriptions.
# The simplest method may be making them properties of the 'SVG' object
# appearing in these two fixed places in the output.
class SVG(object):
  """'SVG' encapsulates the rudimentary standalone DOM model for
     .svg files to facilitate generating those files programmatically."""

  _XML_header= """<?xml version="1.0" encoding="UTF-8"?>"""
  _entity_preamble= \
"""
<!DOCTYPE svg [
"""
# Simplified DOCTYPE declaration Reference:
#   https://jwatt.org/svg/authoring/  (a Mozilla employee)
#   https://www.w3.org/TR/SVGTiny12/intro.html#defining
  _entity_template= """<!ENTITY {0} "{1}">\n"""
  _entity_close= "]>"
  _SVG_preamble= \
"""
<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:xlink="http://www.w3.org/1999/xlink"
     enable-background="new 0 0 &x_view; &y_view;"
     width="&x_extent;cm" height="&y_extent;cm"
     viewBox="0 0 &x_view; &y_view;" x="0px"  y="0px">
"""
  _close_SVG="</svg>\n"

  def __init__(self,extent=(10,10),view=(512,512),origin=(0,0)):
    """Create a new .svg file proxy for a display image contained in
         'extent[0]' cm of x-coordinate and
         'extent[1]' cm of y-coordinate where
       'Graphic' objects are defined using the traditional coordinates
       of mathematics consistent with
         the positive x-axis extends to the right
         the positive y-axis extends upward
       The output file generates data for each 'Graphic' with resolution
         'view[0]' x pixels wide by
         'view[1]' y pixels tall in
       accordance with the SVG specification with the upper left
       corner as origin and data consistent with
         the positive x-axis extending to the right
         the positive y-axis extending downward
       while maintaining the visual orientation of each 'Graphic' in the
       view consistent with traditional mathematical convention.
       The output view data is automatically mapped and interpolated to
       the display extent at display resolution by the user agent."""
    self._indent=0; self._tab="  " # two spaces
    self._entities={} # entity name => value
    self._group=g=Group()
    g.xform=Transform( \
              translate=(view[0]/2.0-origin[0],view[1]/2.0+origin[1]), \
                  scale=(1.0,-1.0))
    self._entities.update(zip(("x_view","y_view"),map(str,view)))
    self._entities.update(zip(("x_extent","y_extent"),map(str,extent)))

  def __add__(self,e):
    """Add another SVG element 'e' to this image."""
    if e is None: return self
    if not issubclass(type(e),Element): raise ValueError, \
      "Object e==%r is not a sub-class of 'Element'." % (e,)
    self._group+e
    return self

  def _entity_update(self):
    for e in self._group:
      if issubclass(type(e),Graphic): self._entities.update(iter(e))

  def indent(self,n=1): self._indent+=int(n)

  def dedent(self,n=1): self._indent-=int(n)

  @property
  def group(self):
    """Return a reference to this svg's outer group."""
    return self._group

  @property
  def tab(self):
    """Return a string of spaces == "  "*<# tabs>."""
    return self._tab*self._indent

  @tab.setter
  def tab(self,n):
    "'n' is the number of tabs (default == two spaces)."
    self._indent=int(n)

  @property
  def xml(self):
    """Yield an html/xml/svg conformant text string representing the graphic."""
    # header and preamble
    self.tab=0
    s=self.tab+self._XML_header+self._entity_preamble
    # entity declaration
    self.indent(); self._entity_update()
    for t in sorted(self._entities.items()):
      s+=self.tab+self._entity_template.format(*t)
    s+=self.tab+self._entity_close
    self.dedent()
    # svg body text
    s+=self.tab+self._SVG_preamble; self.indent()
    s+=self._group.xml(self)
    self.dedent()
    s+=self.tab+self._close_SVG
    return s

  '''
  def Gradient(self,*args): return self._Gradient(self,*args)
  class _Gradient(_Element):
    """ """
    def __init__(self):
      """"""
      pass

    @property
    def entities(self): return []

    @property
    def xml(self):
      pass
  '''

class Element(object):
  """Base class for all SVG element types."""
  __slots__=['_clp','_transform','_merge','_id','_ll','_ur','_parent']
  def __init__(self):
    self._clp=None; self._transform=None; self._merge=False
    self._ll=self._ur=O; self._id=''; self._parent=None

  def __str__(self):
    return ' '.join([s for s in [self.clip, self.transform] if s])

  @property
  def span(self):
    """Return a pair of 'C' values representing the lower left corner and
       upper right corner of user unit space spanned by all components in
       this 'Element' in its transformed display space."""
    m = self.Matrix() if self._transform is None else self._transform.merge
    xs,ys = zip(self._ll,self._ur)  # bounding box corner coordinates
    xsys = zip(*map(m,[C(u,v) for u in xs for v in ys])) # transform corners
    return C(map(min,xsys)),C(map(max,xsys)) # transformed bounding box

  @property
  def clip(self):
    """Return xml attribute string invoking any defined clipping path."""
    t='clip-path="url(#{})"'
    return t.format(self._clp.key.split('"')[-2]) if self._clp else ''

  @clip.setter
  def clip(self,clip_path): self._clp=clip_path

  @property
  def transform(self):
    """Return xml attribute string invoking any defined transformations."""
    if not self._transform: return ''
    t = self._transform if not self._merge else self._transform.merge
    return 'transform="{}"'.format(t)

  @property
  def xform(self):
    """Get or set 'Transform' object instance on this element."""
    return self._transform
  @xform.setter
  def xform(self,transform): self._transform=transform

  @property
  def merge(self):
    """True when this instance collapses the transformation stack to
       a single matrix for xml output."""
    return self._merge

  @merge.setter
  def merge(self,boolean): self._merge=bool(boolean)

  @property
  def key(self): return self._id


class Graphic(Element):
  """Base class for SVG image graphic element classes.
     Class properties are used to represent element attributes."""
  __slots__=['_entities','_f','_s','_sw','_slj']
  def __init__(self):
    """Create a new renderable element."""
    super(Graphic,self).__init__(); self._entities={}
    self._f="White"; self._s="Black"; self._sw=0.0
    self._slj="miter"

  def __add__(self,name,value):
    """Add an SVG entity definition tuple (name,value) to this element."""
    self._entities[name]=value
    return self

  def __iter__(self): return iter(self._entities.items())

  def __str__(self): # '^' indicates insert newline and tab
    s0=super(Graphic,self).__str__()
    l=[self.fill, self.stroke, self.stroke_width, self.stroke_linejoin]
    s1=' '.join([a for a in l if a])
    return '^'.join([s for s in [s0,s1] if s])

  @property
  def fill(self):
    return 'fill="{}"'.format(self._f) if self._f else ''

  @fill.setter
  def fill(self,c): self._f=str(c)

  @property
  def stroke(self): # spec default is "none": also null when width==0
    include=(self._s and (self._s!="none")) # >lower case< for svg
    return 'stroke="{}"'.format(self._s) if (include and self._sw) else ''

  @stroke.setter
  def stroke(self,c): self._s=str(c)

  @property
  def stroke_width(self): # spec default==1 null when stroke is null:
    include=(self._s and (self._s!="none")) # >lower case< for svg
    return 'stroke-width="{}"'.format(_(self._sw)) if include else ''

  @stroke_width.setter
  def stroke_width(self,v): self._sw=float(v)

  @property
  def stroke_linejoin(self):
    return 'stroke-linejoin="{}"'.format(self._slj) if self._slj else ''

  @stroke_linejoin.setter
  def stroke_linejoin(self,join="round"): self._slj=join


class Container(Element):
  """Base class for SVG image container element classes."""
  __slots__=['_elements']
  def __init__(self):
    """Setup list of contained elements."""
    super(Container,self).__init__(); self._elements=[]

  def __add__(self,e):
    """Add another SVG sub-element 'e' to this element."""
    if e is None: return self
    if not issubclass(type(e),Element): raise ValueError, \
      "Object e=={!r} is not a sub-class of 'Element'.".format(e)
    if e._parent is not None: raise ReferenceError, \
      "Element {!r} is already contained by element {!r}".format(e,e._parent)
    e._parent=self; self._elements.append(e)
    return self

  def __iter__(self): return iter(self._elements)

  def __str__(self):
    s=super(Container,self).__str__()
    return ' '+s if s else ''

  def empty(self):
    """Remove all elements in this container."""
    self._elements=[]

  def remove(self,old):
    """Remove element copies of 'old' from this container. NOP if not found."""
    while True:
      try: self._elements.remove(old)
      except ValueError: return

  def replace(self,old,new):
    """Find the element 'old' in this container and replace it with 'new'.
       If 'old' is not found, append 'new'."""
    try:
      i=self._elements.index(old); self._elements[i]=new
    except ValueError: self._elements.append(new)
    return new

  @property
  def span(self):
    """Return a pair of 'C' values representing the lower left corner and
       upper right corner of user unit space spanned by all components in
       this 'Container' in their transformed display space."""
    ll,ur = zip(*(e.span for e in self))
    xmin,ymin = map(min,zip(*ll)); xmax,ymax = map(max,zip(*ur))
    self._ll,self._ur = C(xmin,ymin),C(xmax,ymax)
    return Element.span # computes transformed bounding box


class Group(Container):
  """A class to encapsulate an svg group container."""
  __slots__=[]
  _key=Key("group") # class attribute to sequence instance ids

  def __init__(self,key=None):
    """Create a new svg group container."""
    super(Group,self).__init__()
    if key is None: self._id=next(Group._key)
    else: self._id=next(key)

  def __str__(self): return self.key+super(Group,self).__str__()

  def xml(self,svg):
    attributes=str(self).replace('^','\n    '+svg.tab)
    s=svg.tab+'<g {}>\n'.format(attributes)
    svg.indent()
    for e in self: s+=e.xml(svg)
    svg.dedent()
    s+=svg.tab+"</g>\n"
    return s


class Defs(Container):
  """A class to encapsulate an svg definition container."""
  __slots__=[]
  _key=Key("defs") # class attribute to sequence instance ids

  def __init__(self,key=None):
    """Create a new svg definition container."""
    super(Defs,self).__init__()
    if key is None: self._id=next(Defs._key)
    else: self._id=next(key)

  def __str__(self): return self.key+super(Defs,self).__str__()

  def xml(self,svg):
    attributes=str(self).replace('^','\n    '+svg.tab)
    s=svg.tab+'<defs {}>\n'.format(attributes)
    svg.indent()
    for e in self: s+=e.xml(svg)
    svg.dedent()
    s+=svg.tab+"</defs>\n"
    return s


class ClipPath(Container):
  """A class to encapsulate an svg clipping path container."""
  __slots__=[]
  _key=Key("clip_path") # class attribute to sequence instance ids

  def __init__(self,key=None):
    """Create a new svg clipping path container."""
    super(ClipPath,self).__init__()
    if key is None: self._id=next(ClipPath._key)
    else: self._id=next(key)

  def __str__(self): return self.key+super(ClipPath,self).__str__()

  def xml(self,svg):
    attributes=str(self).replace('^','\n    '+svg.tab)
    s=svg.tab+'<clipPath {}>\n'.format(attributes)
    svg.indent()
    for e in self: s+=e.xml(svg)
    svg.dedent()
    s+=svg.tab+"</clipPath>\n"
    return s


class Use(Graphic):
  """A class to encapsulate an svg use def reference."""
  __slots__=['_def']
  _key=Key("use") # class attribute to sequence instance ids

  def __init__(self,definition,key=None):
    """'definition' is an svg defs element that this use reference will render."""
    super(Use,self).__init__()
    if not isinstance(definition._parent,Defs): raise ReferenceError, \
      "Definition 'Element' {!r} is not contained by {!r}, not a 'Defs' 'Element'". \
         format(definition,definition._parent)
    self._def=definition; self.fill=''; self.stroke=''; self.stroke_linejoin=''
    if key is None: self._id=next(Use._key)
    else: self._id=next(key)

  def __str__(self):
    return ' '.join([self.key, self.link, super(Use,self).__str__()])

  @property
  def link(self): return 'xlink:href="#{}"'.format(self._def.key.split('"')[1])

  @property
  def span(self):
    """Return a pair of 'C' values representing the lower left corner and
       upper right corner of user unit space spanned by all components in
       this 'Container' in their transformed display space."""
    ll,ur = zip(*(e.span for e in self._def))
    xmin,ymin = map(min,zip(*ll)); xmax,ymax = map(max,zip(*ur))
    self._ll,self._ur = C(xmin,ymin),C(xmax,ymax)
    return Element.span # computes transformed bounding box

  def xml(self,svg):
    attributes=str(self).replace('^','\n    '+svg.tab)
    return svg.tab+'<use {0}/>\n'.format(attributes)


class Transform(object):
  """Container for an SVG element transformation attribute string.
     A transform is a list of operations applied to a graphic element
     from last to first.  The transform object is managed as a stack
     where the last item is always the first to be applied."""
  __slots__=['_stack','_id']
  _key=Key("xform") # class attribute to sequence instance ids

  _templates={
    't': "translate({0},{1})",
    's': "scale({0},{1})",
    'r': "rotate({0})",
    'x': "skewX({0})",
    'y': "skewY({0})",
    'm': "matrix({0},{1},{2},{3},{4},{5})"
    }

  def __init__(self, translate=None, rotate=None, scale=None, key=None):
    """Parameters are:.
         'translate' : tuple of x,y coordinates
         'rotate'    : angle of rotation in degrees
         'scale'     : tuple of x,y multiplicative scaling
       The scale, if present, is applied first; the rotation, if
       present, next; and finally the translation, if present.
       More complex transformations may be composed with the methods.
    """
    self._stack=[]
    if key is None: self._id=next(Transform._key)
    else: self._id=next(key)
    if translate: self.translate(*translate)
    if rotate: self.rotate(rotate)
    if scale: self.scale(*scale)

  def __str__(self): # '^' indicates insert newline and tab
    l=[self._templates[t[0]].format(*map(_,t[1:-1])) for t in self._stack]
    l=[l[max(0,s):s+3] for s in range(-1,len(l),3)] # 3 per line after first 2
    return '^  '.join([' '.join([s for s in sl if s]) for sl in l])

  def __len__(self):
    """Return the number of stacked transformations."""
    return len(self._stack)

  def clear(self):
    """Empty the transformation stack."""
    self._stack=[]

  def pop(self):
    """Remove the last transformation from the transform list."""
    self._stack.pop()

  # stored tuple format: (template_key,..svg_args..,matrix_rep)
  def translate(self,tx=0.0,ty=0.0):
    """Push a translation of x by tx and y by ty on the transform list."""
    self._stack.append( ('t',tx,ty,Matrix((1.0,0.0,tx),(0.0,1.0,ty)) ) )

  def scale(self,sx,sy):
    """Push a scale in x by sx and in y by sy on the transform list."""
    self._stack.append( ('s',sx,sy,Matrix((sx,0.0,0.0),(0.0,sy,0.0)) ) )

  def rotate(self,a):
    """Push a rotation by angle 'a' in degrees on the transform list."""
    a=Angle(a,'deg'); z=a.unit; c,s = z
    self._stack.append( \
        ('r',a.deg,Matrix((c,-s,0.0),(s,c,0.0)) ) )

  def skewX(self,a):
    """Push a skew of the x coordinate by angle 'a' in degrees
       on the transform list."""
    a=Angle(a,'deg'); t=tan(a.rad)
    self._stack.append( \
        ('x',a.deg,Matrix(row0=(1.0,t,0.0)) ) )

  def skewY(self,a):
    """Push a skew of the x coordinate by angle 'a' in degrees
       on the transform list."""
    a=Angle(a,'deg'); t=tan(a.rad)
    self._stack.append( \
        ('y',a.deg,Matrix(row1=(t,1.0,0.0)) ) )

  def matrix(self,sx,kx,tx,ky,sy,ty):
    """Push a general matrix transformation, where for each 'x' or 'y'
       component the prefix 't' is the translation, 's' is the scale,
       and 'k' is the skew component.""" # ordering is correct if strange
    self._stack.append( \
        ('m',sx,ky,kx,sy,tx,ty,Matrix((sx,kx,tx),(ky,sy,ty)) ) )

  @property
  def merge(self):
    """Calculate the single matrix equivalent of the current stack.
       Each access returns a new 'Transform' object and leaves the
       stack of this object unchanged."""
    m=reduce(lambda m,n: m*n, [t[-1] for t in self._stack], Matrix())
    t=Transform(); t.matrix(*(e for v in m._m[0:2] for e in v))
    return t

  @property
  def key(self): return self._id


class Circle(Graphic):
  """A class to encapsulate an svg circle element definition."""
  __slots__=['_r','_c']
  _key=Key("disk") # class attribute to sequence instance ids

  def __init__(self,radius=1.0,center=C(0.0,0.0),key=None):
    """Create a new circle at 'center' with size 'radius'."""
    super(Circle,self).__init__()
    self._r=radius; self._c=C(center); self.stroke_linejoin=''
    dvabs=C(abs(radius),abs(radius))
    self._ll,self._ur = self._c-dvabs,self._c+dvabs
    if key is None: self._id=next(Circle._key)
    else: self._id=next(key)

  def __str__(self): # '^' indicates insert newline and tab
    l=[self.key, self._cxcy, self._radius]
    s0,s1 = ' '.join([a for a in l]),super(Circle,self).__str__()
    return s0+'^'+s1 if s1 else s0

  @property
  def _cxcy(self): return 'cx="{0}" cy="{1}"'.format(*map(_,self._c))

  @property
  def center(self): return self._c

  @center.setter
  def center(self,c=C(0.0,0.0)): self._c=C(c)

  @property
  def _radius(self): return 'r="{}"'.format(_(self._r))

  @property
  def radius(self): return self._r

  @radius.setter
  def radius(self,r): self._r=float(r)

  def xml(self,svg):
    attributes=str(self).replace('^','\n  '+svg.tab)
    return svg.tab+'<circle {}/>\n'.format(attributes)


class Ellipse(Graphic):
  """A class to encapsulate an svg ellipse element definition."""
  __slots__=['_v','_c']
  _key=Key("oval") # class attribute to sequence instance ids

  def __init__(self,vertex=C(1.0,.5),center=C(0.0,0.0),key=None):
    """Create a new ellipse tangent to the edges of a rectangle at
       'center' with corner 'vertex'."""
    super(Ellipse,self).__init__()
    self._v=C(vertex); self._c=C(center); self.stroke_linejoin=''
    dvabs=C(map(abs,self._v))-self._c
    self._ll,self._ur = self._c-dvabs,self._c+dvabs
    if key is None: self._id=next(Ellipse._key)
    else: self._id=next(key)

  def __str__(self): # '^' indicates insert newline and tab
    l=[self.key, self._cxcy, self._rxry]
    s0,s1 = ' '.join([a for a in l]),super(Ellipse,self).__str__()
    return s0+'^'+s1 if s1 else s0

  @property
  def _cxcy(self): return 'cx="{0}" cy="{1}"'.format(*map(_,self._c))

  @property
  def center(self): return self._c

  @center.setter
  def center(self,c=C(0.0,0.0)): self._c=C(c)

  @property
  def _rxry(self): return 'rx="{0}" ry="{1}"'.format(*map(_,self._v))

  @property
  def vertex(self): return self._v

  @vertex.setter
  def vertex(self,v): self._v=C(v)

  def xml(self,svg):
    attributes=str(self).replace('^','\n  '+svg.tab)
    return svg.tab+'<ellipse {}/>\n'.format(attributes)


class Rectangle(Graphic):
  """A class to encapsulate an svg rectangle element definition."""
  __slots__=['_v','_c']
  _key=Key("rect") # class attribute to sequence instance ids

  def __init__(self,vertex=C(.5,.5),center=C(0.0,0.0),key=None):
    """Create a new rectangle at 'center' defined by one 'vertex'."""
    super(Rectangle,self).__init__()
    self._v=C(vertex); self._c=C(center); dvabs=C(map(abs,self._v))-self._c
    self._ll,self._ur = self._c-dvabs,self._c+dvabs
    if key is None: self._id=next(Rectangle._key)
    else: self._id=next(key)

  def __str__(self): # '^' indicates insert newline and tab
    l=[self.key, self._xy, self._dxdy]
    s0,s1 = ' '.join([a for a in l]),super(Rectangle,self).__str__()
    return s0+'^'+s1 if s1 else s0

  @property
  def _xy(self):
    halfdiag=C(map(abs,self._c-self._v))
    return 'x="{0}" y="{1}"'.format(*map(_,self._c-halfdiag))

  @property
  def center(self): return self._c

  @center.setter
  def center(self,c=C(0.0,0.0)): self._c=C(c)

  @property
  def _dxdy(self):
    diag=2.0*C(map(abs,self._c-self._v))
    return 'width="{0}" height="{1}"'.format(*map(_,diag))

  @property
  def vertex(self): return self._v

  @vertex.setter
  def vertex(self,v=C(.5,.5)): self._v=C(v)

  def xml(self,svg):
    attributes=str(self).replace('^','\n  '+svg.tab)
    return svg.tab+'<rect {}/>\n'.format(attributes)


class Text(Graphic):
  """A class to encapsulate an svg text element definition."""
  __slots__=['_anch','_s','_v','_t','_fnt','_domb']
  _key=Key("text") # class attribute to sequence instance ids

  def __init__(self,size=5, vertex=C(-2.0,-2.0), \
                    text="none", font="Arial", key=None):
    """Create a new text element with lower left corner 'vertex', 'size'
       pixels tall.  Use 'font' when rendering the string passed in 'text'.
       The default text color is '#404040'."""
    super(Text,self).__init__()
    self._anch="start"; self._s=size; self._v=C(vertex); self._t=text
    self._fnt=font; self.fill="#404040"; self.stroke_linejoin=''
    self._domb="auto"; self.xform=Transform() # empty transformation
    if key is None: self._id=next(Text._key)
    else: self._id=next(key)

  def __str__(self): # '^' indicates insert newline and tab
    s0=' '.join([a for a in [self.key, self.anchor]])
    if self._domb!="auto": s0=s0+'^'+self.dominant_baseline
    s1=' '.join([a for a in [self.font, self.size]])
    s0=s0+'^'+s1 if s1 else s0
    # correct for Cartesian coordinate mirror and initially position text
    self._transform.translate(*self._v); self._transform.scale(1.0,-1.0)
    s1=super(Text,self).__str__(); self._transform.pop(); self._transform.pop()
    return s0+'^'+s1 if s1 else s0

  @property
  def anchor(self): return 'text-anchor="{}"'.format(self._anch)

  @anchor.setter
  def anchor(self,position="start"): self._anch=position

  @property
  def dominant_baseline(self): return 'dominant-baseline="{}"'.format(self._domb)

  @dominant_baseline.setter
  def dominant_baseline(self,position="auto"): self._domb=position

  @property
  def font(self): return 'font-family="{}"'.format(self._fnt)

  @font.setter
  def font(self,font="Arial"): self._fnt=font

  @property
  def size(self): return 'font-size="{}"'.format(_(self._s))

  @size.setter
  def size(self,size=50): self._s=size

  @property
  def text(self): return self._t

  @text.setter
  def text(self,text="none"): self._t=text

  @property
  def vertex(self): return self._v

  @vertex.setter
  def vertex(self,v=C(-2.0,-2.0)): self._v=C(v)

  @property
  def span(self):
    """Return a pair of 'C' values representing the lower left corner and
       upper right corner of user unit space spanned by all characters in
       this 'Text' in their transformed display space.  This is an estimate
       assuming an 'em' box size x 2*size (max-width x (ascent+descent))
       for all characters."""
    ll,ur = self._a-C(0.0,self._s),self._a+C(len(self._t*self._s,self._s))
    # correct for Cartesian coordinate mirror and initially position text
    self._transform.translate(*self._v); self._transform.scale(1.0,-1.0)
    sp=Element.span; self._transform.pop(); self._transform.pop()
    return sp

  def xml(self,svg):
    self.merge = len(self._transform)>0 # more than orientation and position
    attributes=str(self).replace('^','\n  '+svg.tab)
    return svg.tab+'<text {0}>{1}</text>\n'.format(attributes,self.text)


class Path(Graphic):
  """A class to encapsulate an svg path element definition.
     Currently absolute and relative move, line and arc are supported.
     Horizontal and vertical lines are specified as ordinary lines.
     Due to automated precision control, the special commands would
     only save 4 characters, and would harm readability.
  """
  __slots__=['_p','_start','_pos','_pts']
  _key=Key("path") # class attribute to sequence instance ids

  _templates={
    'M': "^  M {0},{1}",
    'L': "^  L {0},{1}",
    'A': "^  A {0},{1} {2} {3:.0f} {4:.0f} {5},{6}",
    'm': "^  m {0},{1}",
    'l': "^  l {0},{1}",
    'a': "^  a {0},{1} {2} {3:.0f} {4:.0f} {5},{6}",
    'z': "^  z"
    }
  def __init__(self,key=None):
    """Create a new empty path element."""
    super(Path,self).__init__()
    self._p=[]; self._start=self._pos=C(0.0); self._pts=[]
    if key is None: self._id=next(Path._key)
    else: self._id=next(key)

  def __str__(self): # '^' indicates insert newline and tab
    s0=self.key+' d="'
    for c in self._p:
      s0+=self._templates[c[0]].format(*map(_,c[1]))
    s0+='^"'; s1=super(Path,self).__str__()
    return s0+'^'+s1 if s1 else s0

  def move(self,z=C(0.0)):
    z=self._pos=self._start=C(z)
    self._p.append(('M',z)); self._pts.append(z)

  def line(self,z=C(1.0)):
    z=self._pos=C(z); self._p.append(('L',z)); self._pts.append(z)

  def arc(self,radius=1.0, end=C(1.0), \
            clockwise=True, bigarc=False, orient=0):
    """An arc is defined from the current position to the end position
       along an elliptical path.  There are four such paths.  'big' and
       'clockwise' select which is rendered.  When 'radius' is a scalar,
       a circle is indicated; when complex, an ellipse.  'orient' is the
       relative angle in degrees of the axes of the path defining ellipse."""
    end=self._pos=C(end); rx,ry = C(radius)
    if ry==0.0: ry=rx
    cw = 1 if clockwise else 0; b = 1 if bigarc else 0
    self._p.append(('A',(rx,ry,orient,b,cw,end.x,end.y))); self._pts.append(end)

  def rmove(self,z=C(0.0)):
    z=C(z); self._pos+=z; self._start+=z
    self._p.append(('m',z)); self._pts.append(self._pos)

  def rline(self,z=C(1.0)):
    z=C(z); self._pos+=z
    self._p.append(('l',z)); self._pts.append(self._pos)

  def rarc(self,radius=1.0, end=C(1.0), \
            clockwise=True, bigarc=False, orient=0):
    """The endpoint 'end' is relative; otherwise the same as 'arc'."""
    end=C(end); self._pos+=end; rx,ry = C(radius)
    if ry==0.0: ry=rx
    cw = 1 if clockwise else 0; b = 1 if bigarc else 0
    self._p.append(('a',(rx,ry,orient,b,cw,end.x,end.y))); self._pts.append(self._pos)

  def close(self):
    self._pos=self._start; self._p.append(('z',())); self._pts.append(self._pos)

  @property
  def cursor(self):
    """The current absolute position of the path end."""
    return self._pos

  @property
  def span(self):
    """Return a pair of 'C' values representing the lower left corner and
       upper right corner of user unit space spanned by all vertices in
       this 'Path' in their transformed display space."""
    xs,ys = zip(*(e for e in self._pts))
    xmin,ymin = map(min,zip(*ll)); xmax,ymax = map(max,zip(*ur))
    self._ll,self._ur = C(xmin,ymin),C(xmax,ymax)
    return Element.span # computes transformed bounding box

  def xml(self,svg):
    if not self._p: return ""
    attributes=str(self).replace('^','\n  '+svg.tab)
    return svg.tab+'<path {0}/>\n'.format(attributes)


if __name__=="__main__":
  svg=SVG()
  clp=ClipPath()
  clp+=Circle(70,complex(75,0)); svg+=clp
  g=Group(); g.clip=clp; svg+=g
  g+=Circle(300)
  e=Ellipse(complex(100,50))
  e.fill="Yellow"; g+=e
  r=Rectangle(complex(50,50))
  r.fill="Blue"; g+=r
  p=Path()
  v=sqrt(150)
  z0,z1,z2 = complex(50),complex(-25,v/2.0),complex(-25,-v/2.0)
  p.move(z0);p.arc(50,z1);p.line(z2);p.line(z0);p.close(); g+=p
  t=Text(text="Greetings"); t.fill="Orange"; g+=t
  t.xform=Transform(rotate=35.26,translate=(40,0),scale=(5,5))
  d=Defs();
  r=Rectangle(complex(50,50)); r.fill="Green";
  t=Transform(rotate=10); t.skewX(12); t.skewY(24)
  t.rotate(45); t.translate(100,0)
  r.xform=t.merge
  d+=r; u=Use(r); u.clip=clp; g+=d; g+=u
  print svg.xml
  m=Matrix((2.,1.,0.),(-1.,2.,0.),(0.,0.,2.))
  print "Matrix inversion test:\n", m,' x \n',m.inverse,' = \n',m*m.inverse
