#cython: nonecheck=False, boundscheck=False
"""
Package: geometry
Module : gauss.pyx
Website: https://github.com/ruminations/site-packages
License: https://github.com/ruminations/Licenses#design-license
Initial Copyright 2017

An implementation of the Gaussian geometry of the complex plane.  This
differs from the builtin type 'complex' by defining several properties
and operators that are geometrically useful.  It is NOT a sub-type of
the builtin, but will interoperate with the builtin, coercing the builtin
to the types defined here.

The class 'C' implements complex number arithmetic, operations, and
properties.  The values 'U', 'V', and 'O' defined here are the basis
vectors C(1.0,0.0), C(0.0,1.0), and origin C(0.0,0.0) respectively.

The class 'Angle' is a convenience class for building rotation operators
with 'C' that encapsulates unit conversions and supplement calculation.

The class 'Moebius' is an operator constructed from and applicable to
'C' objects.  It facilitates a unified representation of general
transformations of the Gaussian plane. *** Only implemented in gauss.py ***

Implementation notes:

This module is a 2-dimensional cython optimization of the space_4D.py
class 'V' renamed here to 'C' that implements complex number arithmetic.
The value of this over the builtin type 'complex' is that the geometric
and vector properties of complex numbers may be manipulated as easily as
the algebraic field properties.  The module preserves the application
programming interface of space.py to the extent possible, and provides
all of the functionality of the builtin complex number interface in a
geometrically and algbraically more intuitive style.

In particular, the package maintains support of arithmetic notation for
operations on 2-vectors (complex numbers).  The notation maintains
augmentation by the Dirac notation for vector inner products, e.g.
(u|v) computes the inner product of vectors 'u' and 'v'.  The inversion
operator '~' is used for complex conjugation: (~C(a,b)) == C(a,-b).

In addition the operators '&' and '//' implement component-wise
multiplication and division of 'C' values. Thus, '(u|v) == sum(u&v)' and
'(v//v) == C(1.0,1.0) == U+V', whereas 'v/v == C(1.0) == U'.

On the other hand, slicing a two element object makes little sense.  In
addition, the altered semantics of python 3 introduce quite a bit of
overhead when __getitem__ supports slices.  Therefore slicing is omitted.
"""
#### TODO:
# 1) complete coverage in test suite
# 2) incorporate 'Moebius' class when it is complete
# 3) ultimately this module needs to be made into a bonafide c-extension
# 4) there are some accuracy failures without multiply accumulate -
#       the fma function does not work well for multiple accumulations
#       and C's rounding of intermediaries introduces cumulative error
#       Possible solution: change rounding mode to toward zero
__version__=20180413
__version_history__=(20170207,20171225,20180325)
__all__=['Angle','C','U','V','O','Cline','Ccircle']
from algebra.mathx import sign,zero,clean,delta_ratio_is_zero

cdef extern from "math.h":
  cdef double fabs( double )
  cdef double sqrt( double )
  cdef double atan2( double, double )
  cdef double cos( double )
  cdef double sin( double )
  double M_PI

# comparison operator dispatch table
cdef short _lt_(double a, double b): return a< b
cdef short _le_(double a, double b): return a<=b
cdef short _eq_(double a, double b): return a==b
cdef short _ne_(double a, double b): return a!=b
cdef short _gt_(double a, double b): return a> b
cdef short _ge_(double a, double b): return a>=b
ctypedef short (*cmp_op)(double,double)
cdef cmp_op _op_[6]
_op_[:] = [&_lt_,&_le_,&_eq_,&_ne_,&_gt_,&_ge_]
# angular unit indices
cdef enum: hemi=0
cdef enum: turn=1
cdef enum: deg =2
cdef enum: rad =3
# angular conversion constants indexed by indices
cdef double _conv[4]
_conv[:] = [1.0,0.5,180.0,M_PI]

cdef class Angle(object):
  """
  The 'Angle' object has four attributes:
    'hemi':0  The representation of the angle in hemicycles
    'turn':1  The representation of the angle in turns
    'deg' :2  The representation of the angle in degrees
    'rad' :3  The representation of the angle in radians

  A hemicycle is the angle expressed as a fraction of a semicircle.
  For example, the hemicycles -1, -.5, -.25, 0, .25, .5, 1 correspond,
  respectively, to the degree measures -180, -90, -45, 0, 45, 90, 180.
  This is useful in conjunction with frequency calculations.

  'deg' and 'rad' properties are normalized to +-180 and +-pi; 'hemi' and
  'turn' properties are unnormalized.  The unnormalized value is used for
  arithmetic operations so that 'turn' units accumulate.  Equivalent
  angles that differ in turn measure will test not equal, and comparisons
  are also in an absolute sense.

  'Angle' behaves essentially like a float with properties that provide
  convenient access to angular unit conversions.  It is more minimal -
  only basic arithmetic is supported.  In addition to the properties,
  the invert operator returns the normalized supplementary 'Angle':
    Angle(2/3) == ~Angle(1/3); Angle(-2/3) == ~Angle(-1/3)
  """
  cdef double __rad
  cdef double __hemi

  def __cinit__(self, object angle not None, unit=0):
    cdef double a
    if isinstance(angle,Angle): unit=0
    elif isinstance(unit,str):
      if unit=="hemi": unit=hemi
      if unit=="turn": unit=turn
      if unit=="deg": unit=deg
      if unit=="rad": unit=rad
    elif unit<0: unit=hemi
    elif unit>3: unit=rad
    a=angle.__float__()/_conv[unit]; self.__hemi=a
    while a>1.0: a-=2.0
    while a<=-1.0: a+=2.0
    self.__rad=a*M_PI

  def __init__(self,*arg):
    """Cython signature: (object angle not None, unit=0)
       Create an angle 'angle' using units from enumeration 'unit' in
         ('hemi', 'turn', 'deg', 'rad') == (0,1,2,3).
       'hemi' is the default and is used internally. 'angle' may also be
       of type 'Angle' in order to make a copy.  'unit' may be a string."""
    pass # used to display __cinit__ docstring under help(...)

  def __abs__(self): return Angle(fabs(self.__hemi))

  def __add__(object u not None, object v not None): # cython has no reversed ops
    cdef double a, b
    a=u.__float__(); b=v.__float__()
    return Angle(a+b)

  def __div__(object u not None, object v not None): # cython has no reversed ops
    cdef double a, b
    a=u.__float__(); b=v.__float__()
    return Angle(a/b)

  def __float__(self): return self.__hemi

  def __invert__(self):
    """Return the supplement of this 'Angle'."""
    if self.__rad<0.0: return Angle(-M_PI-self.__rad,3)
    return Angle(M_PI-self.__rad,3)

  def __mul__(object u not None, object v not None): # cython has no reversed ops
    cdef double a, b
    a=u.__float__(); b=v.__float__()
    return Angle(a*b)

  def __neg__(self): return Angle(-self.__hemi)

  def __pos__(self): return self

  def __sub__(object u not None, object v not None): # cython has no reversed ops
    cdef double a, b
    a=u.__float__(); b=v.__float__()
    return Angle(a-b)

  def __richcmp__(object u not None, object v not None, short op):
    """cython op|cmp: 0|<, 1|<=, 2|==, 3|!=, 4|>, 5>= """
    cdef double a, b
    a=u.__float__(); b=v.__float__()
    return bool(_op_[op](a,b))

  def __hash__(self): return hash(self.__hemi)

  def __repr__(self): return self.__class__.__name__+"({!r})".format(self.__hemi)

  def __str__(self): return str(self.deg)

  @property
  def deg(self):
    """Return the normalized angle in degree measure."""
    return self.__rad*180.0/M_PI

  @property
  def hemi(self):
    """Return the unnormalized angle in hemi-turn measure."""
    return self.__hemi

  @property
  def rad(self):
    """Return the normalized angle in radian measure."""
    return self.__rad

  @property
  def turn(self):
    """Return the unnormalized angle in turn measure."""
    return self.__hemi/2.0

  @property
  def unit(self):
    """Return a unit 'C' operator that rotates by this 'Angle' object."""
    return C(self)

cdef double _temp[3] # acts as class attribute for 'C'
cdef double* _tmp_cplx= &_temp[0] # pointer to temporary builtin 'complex' data
cdef double* _tmp_sclr= &_temp[1] # pointer to temporary builtin scalar data
_temp[2] = 0.0 # scalar complex component is always zero

cdef C _coerce_C(e):
  """Convert a scalar or 'complex' value to a 'C'."""
  if isinstance(e,C): return e
  return C(e)

cdef double* _coerce_data(e):
  """Convert a scalar or 'complex' or 'C' object to its data."""
  if isinstance(e,C): return C.__get_v(e)
  if isinstance(e,complex):
    _temp[0]=e.real; _temp[1]=e.imag
    return _tmp_cplx
  _temp[1]=e.__float__()
  return _tmp_sclr

cdef C _new(double* data): # @classmethod replacement for 'C'
  """Internal init call speedup avoids integrity checks.
     'data' is a double array pointer."""
  cdef double* that
  inst=C(); that=C.__get_v(inst)
  that[0],that[1] = data[0],data[1]
  return inst

cdef class C(object):
  """Class 'C' represents a complex number as a Euclidean 2-vector.

     A vector 'C' may be constructed from one or two numbers, a complex
     number, a list, or a tuple:
       a=C(1,2); b=C(complex(-2,-1)); c=C([-2.0,-1.0]); d=C((1,0)); e=C(3)
     A unit rotation operator is returned for an 'Angle':
       r=C(Angle(pi/3.0))
     An argument of type 'C' is also accepted and makes a copy.
     'C' is an immutable type, and so may be used as a dictionary key.

     Operator overloading:
       '+'   : Vector or complex addition    <=> a+b == C(-1.0,1.0)
       '-'   : Vector or complex subtraction <=> a-b == C(3.0,3.0)
       '*'   : Multiplication <=>  pi*c , c*pi, a*b
       '/'   : Division <=>  c/pi, pi/c, a/b
       '-'   : Reflection through the origin <=>  -a == C(-1.0,-2.0)
       '~'   : Complex conjugation <=> ~a == C(1.0,-2.0)
       '|'   : Inner product <=>  a|b == -4.0
       '&'   : Component-wise product <=>  a&b  == C(a.x*b.x,a.y*b.y)
       '//'  : Component-wise quotient <=> a//b == C(a.x/b.x,a.y/b.y)
       'abs' : Euclidean magnitude = sqrt(a|a)
       '=='  : Component-wise equality   <=> (a==a,a==b) == True,False
       '!='  : Component-wise inequality <=> (a!=a,a!=b) == False,True
       '<<'  : Rotate components left  <=>  a<<3 == C(2,1)
       '>>'  : Rotate components right <=>  a>>3 == C(2,1)
       'in'  : Containment or iteration <=>  2 in a == True; for e in a
       '[i]' : Extract a component <=>  a[1]==2
       '^'   : a^b is the unit complex rotation operator that aligns a
               with b, i.e.
                 b ==abs(b)/abs(a)*(a^b)*a, where a^b == b.unit/a.unit
               Notice abs(a) == ans((a^b)*a).
               Returns 'None' if abs(a) or abs(b) is zero.
  """
  cdef double __v[2] # ctype data values for this vector
  cdef double __abs # the Euclidean magnitude of this vector
  cdef short __absflag # flag magnitude as valid

  cdef double* __get_v(self): return &self.__v[0]

  def __cinit__(self,*u): # *arg,**kwarg are implicit or can be explicit
    cdef double* this
    self.__abs=0.0;  self.__absflag=0
    if not u: return # class method '_new' is caller: return empty instance
    if len(u)==1: # single argument: scalar, tuple, list, or error
      u=u[0] # floats,ints,longs and bools have .real and .imag properties defined
      if isinstance(u,Angle): u=(cos(u.rad),sin(u.rad)) # rotation operator
      elif isinstance(u,C): u=u.xy # copy a 'C' object
      else:
        try: u=(u.real,u.imag)
        except AttributeError: pass # tuple, list, or error
    try: u=tuple(map(float,u)) # multi-argument or tuple from previous, or error
    except TypeError: raise TypeError,"'C' elements must be numeric scalars."
    if len(u)!=2:
      raise ValueError,"Only two scalar elements form a 'C'; {} given.".format(len(u))
    this=C.__get_v(self)
    this[0],this[1] = u[0],u[1]

  def __abs__(self):
    if self.__absflag: return self.__abs
    self.__abs=sqrt(self|self); self.__absflag=1; return self.__abs

  def __add__(object u not None, object v not None): # cython has no reversed ops
    cdef double rval[2]
    cdef double* this
    cdef double* that
    this,that = _coerce_data(u),_coerce_data(v)
    rval[0],rval[1] = this[0]+that[0],this[1]+that[1]
    return _new(rval)

  def __and__(object u not None, object v not None): # cython has no reversed ops
    """Return component wise product of complex values treated as a
       vector.  Scalars are treated as 'C' objects with zero in the
       imaginary component.  This product is equivalent to treating
       each value as a diagonal matrix, and is useful in some
       2-dimensional graphic applications."""
    cdef double rval[2]
    cdef double* this
    cdef double* that
    this,that = _coerce_data(u),_coerce_data(v)
    rval[0],rval[1] = this[0]*that[0],this[1]*that[1]
    return _new(rval)

  def __call__(self, m): ##### Moebius m not None): # currently unimplemented
    """Return the application of Moebius transform 'm' to 'self'.  If 'c'
       is a 'C' instance and 'm0', 'm1', 'm2', and 'm3' are Moebius instances,
       c(m0)(m1)(m2)(m3) == m3(m2(m1(m0(c))))"""
    raise NotImplementedError, \
      "'gauss.pyx' does not currently support 'Moebius' objects."
    #return m(self)

  def __contains__(self,double e):
    cdef double* this
    this=C.__get_v(self)
    if this[0]==e: return True
    if this[1]==e: return True
    return False

  def __richcmp__(object u not None, object v not None, short op):
    cdef object t,f
    cdef double* this
    cdef double* that
    this,that = _coerce_data(u),_coerce_data(v)
    if op == 2: t,f=True,False # __eq__
    elif op == 3: t,f=False,True # __ne__
    else: raise NotImplementedError,"4-D rich comparison opcode: %s" % (op,)
    if this[0]!=that[0]: return f
    if this[1]!=that[1]: return f
    return t

  def __div__(object u not None, object v not None): # cython has no reversed ops
    cdef double rval[2]
    cdef double n
    cdef double* this
    cdef double* that
    this,that = _coerce_data(u),_coerce_data(v)
    n=that[0]*that[0]+that[1]*that[1] # (v|v)
    rval[0],rval[1] = \
      (this[0]*that[0]+this[1]*that[1])/n,(this[1]*that[0]-this[0]*that[1])/n
    return _new(rval)

  def __floordiv__(object u not None, object v not None): # cython has no reversed ops
    """Return component wise quotient of complex values treated as a
       vector.  Scalars are treated as 'C' objects with zero in the
       imaginary component and therefore raise a 'ZeroDivisionError'.
       This quotient is equivalent to treating each value as a diagonal
       matrix, and is useful in some 2-dimensional graphic applications."""
    cdef double rval[2]
    cdef double* this
    cdef double* that
    this,that = _coerce_data(u),_coerce_data(v)
    rval[0],rval[1] = this[0]/that[0],this[1]/that[1]
    return _new(rval)

  def __getitem__(self, short i):
    # python 3 slice code: need to modify parameter type so slicing is not supported:
    # if isinstance(i,slice): return tuple(self[j] for j in range(*i.indices(2)))
    if i<0: i+=2
    if i>=2 or i<0: raise IndexError,"Index == %r out of range; 'C' length == 2" % (i,)
    return C.__get_v(self)[i]

  def __hash__(self): return hash(self.xy)

  def __invert__(self):
    """Return the complex conjugate of this 'C'."""
    cdef double rval[2]
    cdef double* this
    this=C.__get_v(self)
    rval[0],rval[1] = this[0],-this[1]
    return _new(rval)

  def __iter__(self): return iter(self.xy)

  def __len__(self): return 2

  def __lshift__(self,short index):
    cdef double rval[2]
    cdef double* this
    cdef short t
    this=C.__get_v(self); t=index%2
    if t==1: rval[0],rval[1] = this[1],this[0]
    else: rval[0],rval[1] = this[0],this[1]
    return _new(rval)

  def __mul__(object u not None, object v not None): # cython has no reversed ops
    cdef double rval[2]
    cdef double* this
    cdef double* that
    this,that = _coerce_data(u),_coerce_data(v)
    rval[0],rval[1] = this[0]*that[0]-this[1]*that[1],this[0]*that[1]+this[1]*that[0]
    return _new(rval)

  def __neg__(self):
    cdef double rval[2]
    cdef double* this
    this=C.__get_v(self)
    rval[0],rval[1] = -this[0],-this[1]
    return _new(rval)

  def __or__(object u not None, object v not None): # cython has no reversed ops
    """Return inner product of complex values treated as a vector.  Scalars
       are treated as 'C' objects with zero in the imaginary component."""
    cdef double* this
    cdef double* that
    this,that = _coerce_data(u),_coerce_data(v)
    return this[0]*that[0]+this[1]*that[1]

  def __repr__(self):
    return self.__class__.__name__+"({!r}, {!r})".format(*self)

  def __rshift__(self, short index):
    cdef double rval[2]
    cdef double* this
    cdef short t
    this=C.__get_v(self); t=index%2
    if t==1: rval[0],rval[1] = this[1],this[0]
    else: rval[0],rval[1] = this[0],this[1]
    return _new(rval)

  def __str__(self):
    """'C' components with self.tidy applied."""
    return "[ {!r}  {!r}j ]".format(*self.tidy)

  def __sub__(object u not None, object v not None): # cython has no reversed ops
    cdef double rval[2]
    cdef double* this
    cdef double* that
    this,that = _coerce_data(u),_coerce_data(v)
    rval[0],rval[1] = this[0]-that[0],this[1]-that[1]
    return _new(rval)

  def __xor__(object u not None, object v not None): # cython has no reversed ops
    if isinstance(u,C): v=_coerce_C(v) # called with arguments
    else: u=_coerce_C(u)            # in original order on right
    try: return v.unit/u.unit
    except TypeError: return None # u or v is zero => unit == None

  @property
  def arg(self):
    """Return the polar coordinate angle of this 'C' as an 'Angle'."""
    cdef double* this
    this=C.__get_v(self)
    return Angle(atan2(this[1],this[0]),3)

  @property
  def mag(self):
    """Return the polar coordinate radius == abs(self) of this 'C'."""
    return abs(self)

  @property
  def polar(self):
    """Return the polar coordinate radius,theta tuple (self.mag,self.arg)."""
    return self.mag,self.arg

  @property
  def sqrt(self):
    """Return the complex operator that, when applied twice, is
       equivalent to this operator applied once.  That is, for a unit,
       return a rotation operator about the origin by half the angle.
       Returns 'None' for the origin."""
    cdef double rval[2]
    cdef double m
    cdef double* this
    cdef C u
    u=self.unit
    if u is None: return None
    this=C.__get_v(u); m=sqrt(abs(self)/2.0)
    rval[0]=m*sqrt(this[0]+1.0)
    if this[1]<0.0: rval[1]=-m*sqrt(1.0-this[0])
    else: rval[1]=m*sqrt(1.0-this[0])
    return _new(rval)

  @property
  def ortho(self):
    """Return a unit magnitude orthogonal to this 'C'.
       Return 'None' when this 'C' is the origin.
       Equivalent to a positive 90 degree rotation of the unit."""
    cdef double rval[2]
    cdef double m
    cdef double* this
    this=C.__get_v(self); m=abs(self)
    try:
      rval[0],rval[1] = -this[1]/m,this[0]/m
      return _new(rval)
    except ZeroDivisionError: return None

  @property
  def unit(self):
    """Return this 'C' normalized to unit magnitude.
       Return 'None' when this 'C' is the origin."""
    cdef double rval[2]
    cdef double m
    cdef double* this
    this=C.__get_v(self); m=abs(self)
    try:
      rval[0],rval[1] = this[0]/m,this[1]/m
      return _new(rval)
    except ZeroDivisionError: return None

  @property
  def xy(self):
    """Return the Cartesian x,y tuple representation of this 'C'."""
    cdef double* this
    this=C.__get_v(self)
    return this[0],this[1]

  @property
  def x(self):
    """Return the Cartesian x coordinate or real part of this 'C'."""
    cdef double* this
    this=C.__get_v(self)
    return this[0]

  @property
  def y(self):
    """Return the Cartesian y coordinate or imaginary part of this 'C'."""
    cdef double* this
    this=C.__get_v(self)
    return this[1]

  @property
  def clean(self):
    """This 'C' with components rounded to 15 significant digits.
       Converts, for example:
         2.220446000000013e-16 -> 2.22044600000001e-16
         2.220446000000003e-16 -> 2.220446e-16
         0.5000000000000001    -> 0.5
         0.7499999999999999    -> 0.75
         5000000000000001.      -> 5000000000000000."""
    cdef double rval[2]
    cdef double* this
    this=C.__get_v(self)
    rval[0],rval[1] = clean(this[0],15),clean(this[1],15)
    return _new(rval)

  @property
  def tidy(self):
    """This 'C' with C.clean first applied to components.  Then component
       elements 'e' such that 0.0 <= abs(e) <= 10**(-15) are replaced
       with 0.0 exactly.

       This is useful for printing values clearly and eliminating differences
       that are actually null.  CAUTION: when legitimate data is very small,
       it WILL be corrupted.  In such cases, use C.clean directly."""
    cdef double rval[2]
    cdef double* this
    this=C.__get_v(self)
    rval[0],rval[1] = zero(clean(this[0],15),15),zero(clean(this[1],15),15)
    return _new(rval)

U,V,O = C(1),C(1j),C(0)

cdef class Cline(object):
  """Representation of a line in the complex plane."""
  cdef C _u
  cdef C _v
  cdef C _get_u(self): return self._u
  cdef C _get_v(self): return self._v

  def __cinit__(self, C w, C z):
    self._u=(z-w).unit
    if self._u is None: raise ValueError, \
      "Coincident points defining line:\n    {!r}, {!r}".format(w,z)
    self._v=~self._u<<1 #.ortho
    self._v=self._v*(self._v|w)

  def __init__(self,C w,C z):
    """Cython signature: (C w, C z)
       Given 'C' values 'w' and 'z', create 'C' values u,s such that a
       point 'x' on the line is given by x=u*t+s, where u is a unit, t is
       a scalar, abs(s) is minimal, and both 'w' and 'z' are on the line.

       This is a Plücker coordinate representation where w => s, the start,
       and z-w => u, specifies the direction."""
    pass # used to display __cinit__ docstring under help(...)

  def __call__(self,double t):
    """Calculate and return the complex value on the line corresponding
       to the scalar 't'."""
    return Cline._get_u(self)*t+Cline._get_v(self)

  def __xor__(object k,object l): # cython has no reversed ops
    """Return the unique intersection of self and line 'l'.  If the two
       lines are parallel, return 'None'."""
    cdef C uk
    cdef C ul
    cdef C vk
    cdef C vl
    cdef C c
    cdef C ret
    cdef ulo
    cdef double det
    if isinstance(k,Ccircle): return k^l # no reversed ops rectify
    if isinstance(l,Ccircle): return l^k
    uk,ul = Cline._get_u(k),Cline._get_u(l); ulo=~ul<<1; det=(uk|ulo)
    if zero(det,15)==0.0: return None # parallel
    vk,vl = Cline._get_v(k),Cline._get_v(l); c=C(vk|~uk<<1,vl|ulo)
    ret= C( (C(-(ul.x),uk.x)|c)/det, (C(-(ul.y),uk.y)|c)/det )
    return ret

  def __repr__(self):
    cdef C v
    v = Cline._get_v(self)
    return "Cline({!r},{!r})".format(v,v+Cline._get_u(self))

  def __str__(self):
    cdef C v
    v = (Cline._get_v(self)).tidy
    return "Cline({!r},{!r})".format(v,(v+Cline._get_u(self)).tidy)

  @property
  def U(self):
    """Get the unit defining the direction of the line'."""
    return Cline._get_u(self)

  @property
  def origin(self):
    """Get the starting point (t==0) of the line'."""
    return Cline._get_v(self)

class Ccircle(object):
  """Represents a circle in the complex plane."""
  __slots__=['_c','_r','_u','_v','_p0','_p1','_p2']
  #        center,radius,  uv   , p0   ,init args
  def __init__(self,p0,p1,p2,t=None):
    """Initialization arguments are 'C' objects:
        d       semantics
       ----     -----------
       None      p0,p1,p2 are distinct points on the circle
        0        p0,p1 points on the circle, p2 tangent at p0
        1        p0,p1 points on the circle, p2 tangent at p1
        2        p0 a point on circle, p1 the center, p2 ignored"""
    assert isinstance(p0,C) and isinstance(p1,C), \
      "Invalid circle parameters p0 or p1:\n  {!r},\n  {!r}.".format(p0,p1)
    self._p0,self._p1,self._p2 = p0,p1,p2
    if t==2: # position,center
      self._u,self._c = p0-p1,p1
      self._u,self._r = self._u.unit,abs(self._u)
      try: self._v=self._u.ortho
      except AttributeError: raise ValueError, \
        "Indistinct positions defining circle center and point:\n  {!r}\n == {!r}". \
          format(p1,p0)
      self._p1,self._p2 = self._c+self._r*self._v,self._c-self._r*self._u
      return
    assert isinstance(p2,C), \
      "Invalid circle parameter p2:\n  {!r}.".format(p2)
    dp,br = p1-p0,(p1+p0)/2.0; nr=dp.ortho
    assert nr is not None, \
      "Indistinct positions p0 and p1:\n  {!r}\n  {!r}".format(p0,p1)
    if t is not None: # differential and positions
      nt=p2.ortho
      assert nt is not None, "Tangent p2 can not be zero:\n  {!r}".format(p2)
      if   t==0: bt=p0 # tangent applies to p0
      elif t==1: bt=p1 # tangent applies to p1
      else: raise ValueError,"Circle tangent specifier out of range: {!r}".format(t)
    else: # three positions on the circle
      bt,nt = (p2+p1)/2.0,(p2-p1).ortho
    self._c=Cline(br,br+nr)^Cline(bt,bt+nt)
    if self._c is None: raise ValueError, \
    "Collinear points or tangent defining circle:\n  {!r},\n  {!r},\n  {!r}". \
                                                     format(p0,p1,p2)
    self._u=p0-self._c; self._u,self._r = self._u.unit,abs(self._u)
    assert self._u is not None, \
      "Bug: implies p0==calculated center - can not happen here."
    self._v=self._u.ortho
    if t is None: # p2 is a position
      if (self._v|p1-self._c)<0.0: self._v=-self._v # orient correctly
    else:         # p2 is a differential
      if (self._v|p2)<0.0:         self._v=-self._v # orient correctly
      dp=dp.unit; assert dp is not None,"Bug: implies p0==p1 - can not happen here."
      if (self._v|dp)<0.0: dp= dp.ortho   # dp bisects p1-p0,
      else:                dp=-dp.ortho   #    oriented correctly
      dp=C((self._u|dp),(self._v|dp))     # dp rotates u to bisector
      dp=dp*dp; dp=(dp*dp).unit           # dp rotates u by 2*<angle p1-p0>
      self._p2=self._c+self._u*dp*self._r # p2 is p0 mirrored across p1

  def __call__(self,z):
    """Return:
        -1 : 'z' inside the circle
         0 : 'z' on the circle
         1 : 'z' outside the circle"""
    t,r = abs(z-self._c),self._r
    if delta_ratio_is_zero(t,r,15): return 0.0
    return sign(t-r)

  def __contains__(self,z):
    """Return 'True' for a point on the circle."""
    return self(z)==0.0

  def __repr__(self):
    """Return a string showing internal representation."""
    return self.__class__.__name__+"({!r},{!r},{!r},t=0)".format( \
                                 self._p0,self._p1,self._p2 )

  def __xor__(self,f): ## refer to LyX doc 'Circle_Intersection.lyx'
    """Compute the intersection of this circle with figure 'f'.
       Supported figures: 'Circle' and 'Line'"""
    if isinstance(self,Cline): self,f = f,self # no reversed ops rectify
    sc,sr = self.center,self.radius
    if isinstance(f,Cline):
      fu,fo = f.U,f.origin # line unit, reference point
      do=(sc-fo); xv=(do|fu)*fu; yv=do-xv; d=round(abs(yv),15)
      if d > sr: return None # disjoint
      if d == sr: return (sc-yv,) # tangent
      dxv=sqrt(sr*sr-d*d)*fu
      return (fo+xv+dxv,fo+xv-dxv) # intersecting
    if isinstance(f,Ccircle):
      fc,fr = f.center,f.radius
      rs=sr+fr; rd=fr-sr; u=fc-sc; d,u=abs(u),u.unit
      if d > rs: return None # separated
      if d < abs(rd): return None # containment
      if rd == 0.0 and u is None: return None # coincident
      if d == rs: return (sc+u*sr,) # separate tangent
      if d == abs(rd): return (sc+sr*(-u if fr>sr else u),) # containment tangent
      v=u.ortho; g2=rs*rd; m=(g2/d-d)/2.0 # g2 is signed due to rd !
      x,y = (d-g2/d)/2.0,sqrt(sr*sr-m*m); xv,yv=sc+u*x,v*y
      return (xv+yv,xv-yv) # intersection pair
    raise ValueError,"Type({!r}) == {} is not a 'Line' or a 'Circle'".format(f,type(f))

  def tangent(self,angle=Angle(0.0)):
    """Return the unit tangent at 'Angle' object 'angle' from the reference
       position 'self.vertex' in the direction of 'self.V'."""
    return C(angle)*self._v

  @property
  def center(self):
    """Return the circle's center position."""
    return self._c

  @property
  def radius(self):
    """Return the circle's radius."""
    return self._r

  @property
  def vertex(self):
    """Get the reference vertex on the circle.  This is initialization 'p0'."""
    return self._p0

  @property
  def U(self):
    """Get the circle's reference unit.  This differential points
       from 'self.center' to 'self.vertex'."""
    return self._u

  @property
  def V(self):
    """Get the circle's unit normal differential.  This differential gives
       the direction positive angles move from 'self.vertex'."""
    return self._v
