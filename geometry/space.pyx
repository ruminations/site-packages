#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Package: geometry
Module : space.pyx  *** see conversion notes at module end ***
Website: https://github.com/ruminations/site-packages
License: https://github.com/ruminations/Licenses#design-license
Initial Copyright 2013

Provides a linear algebra package supporting arithmetic notation for
operations on arbitrary vectors and matrices.  The notation is augmented
by the Dirac notation for matrix multiplication, e.g. (u|M|v) computes
the conjugation of matrix 'M' by vectors 'u' and 'v'.

Implementation notes:
  0) Numerical variations have been benchmarked using Python 2.5.4 and
     this version is believed to optimize speed.

  1) The array module is not used because matrices can not be represented
     using them without a <<lot>> of index calculations, and vectors only
     benefit for slicing operations.
"""
__version__=20180413
__version_history__=(20130209,20131106)
__all__=['V','M','SqM']
cdef extern from "math.h":
    cdef double sqrt( double )
from algebra.mathx import zero,clean
from itertools import imap,izip # yields C code speedup w/ operators
from operator import mul as _m, eq as _eq, ne as _ne
from operator import neg as _n, add as _a, sub as _s

class V(object):
  """Class 'V' represents a Euclidean vector.

     A vector 'V' may be constructed from numbers, a list, or a tuple:
       a=V(1,2,3,4); b=V([-4.0,-3.0,-2.0,-1.0]); c=V((-1,0,1))
     'V' is an immutable type, and so may be used as a dictionary key.

     Operator overloading:
       '+'   : Vector addition    <=>  a+b
       '-'   : Vector subtraction <=>  a-b
       '*'   : Scalar multiplication <=>  pi*c , c*pi
       '/'   : Scalar division <=>  c/pi
       '-'   : Reflection through the origin <=>  -a
       '|'   : Inner product <=>  a|b or a|m , m|a and m is a matrix
       'abs' : Euclidean magnitude = sqrt(a|a)
       '=='  : Component-wise equality   <=>  a==a,a==b == True,False
       '!='  : Component-wise inequality <=>  a!=a,a!=b == False,True
       '<<'  : Rotate components left  <=>  a<<3 == V(4,1,2,3)
       '>>'  : Rotate components right <=>  a>>3 == V(2,3,4,1)
       'in'  : Containment or iteration <=>  3 in a == True; for e in a
       '[i]' : Extract a component <=>  a[2]==3
       '[i:j]' : Slicing <=>  a[1:3] == V(2,3); a[::2] == V(1,3)
       '^'   : a^b is the angle between a and b represented as the
               equivalent unit complex number in the plane of a and b.
                 a^b    = complex(cos(t),sin(t)) where
                 cos(t) = (a|b)/(abs(a)*abs(b))
                 sin(t) = sqrt(1-cos(t)*cos(t))
               Returns 'None' if abs(a) or abs(b) is zero.
  """
  __slots__=['__v','__abs']

  @staticmethod
  def numeric(e):
    """Return 'True' if 'e' is an 'int', 'float', or 'complex' instance.
       Override in subclasses permitting extended numeric types."""
    return isinstance(e,int) or isinstance(e,float) or \
           isinstance(e,complex)

  def __init__(self,*u):
    self.__v=self.__abs=None # Realize __slots__
    if not u: return # class method '_' is caller: return empty instance
    num=self.numeric # speedup: local lookup
    if not num(u[0]): u=tuple(*u)
    if all(num(e) for e in u): self.__v=u
    else: raise ValueError,"Vector elements must be numeric."

  @classmethod
  def _(cls,*u):
    """Internal init call speedup avoids integrity checks.
      'u' is a tuple, list, or generator of numeric values."""
    inst=cls(); inst.__v=tuple(*u)
    if len(inst.__v) > 0: return inst
    raise ValueError, "Report bug: null vector definition encountered."

  def __abs__(self):
    if self.__abs is None: self.__abs=sqrt(self|self)
    return self.__abs

  def __add__(self,u): return self._(imap(_a,self.__v,u.v))

  def __contains__(self,e): return e in self.__v

  def __richcmp__(self,u,op): # cython adaptation replacing python comparison ops
    if op == 2: return all(imap(_eq,self.__v,u.v)) # __eq__
    return any(imap(_ne,self.__v,u.v)) # __ne__; op == 3

  def __div__(self,s): return self._(a/s for a in self.__v)

  def __getitem__(self,i): return self.__v[i]

  def __getslice__(self,i,j): return self._(self.__v[i:j])

  def __hash__(self): return hash(self.__v)

  def __iter__(self): return iter(self.__v)

  def __len__(self): return len(self.__v)

  def __lshift__(self,i):
    v=self.__v; l=len(v) # speedup: local lookup
    return self._(v[i:l]+v[0:i])

  def __mul__(self,s): return self._(a*s for a in self.__v)

  __rmul__=__mul__

  def __neg__(self): return self._(imap(_n,self.__v))

  def __or__(self,u):
    if isinstance(u,V): return sum(imap(_m,self.__v,u.v))
    if not isinstance(u,M): return NotImplemented
    if len(self)!=u.rank[0]: raise ValueError, \
      "Multiplicatively incommensurate vector and matrix: %s with %s."\
      % (len(self),u.rank)
    m=(~u).m; v_mtrc=V.__or__; v_=V._ # speedup: local lookup
    return v_(v_mtrc(self,r) for r in m)

  def __ror__(self,u):
    if not isinstance(u,M): return NotImplemented
    if len(self)!=u.rank[1]: raise ValueError, \
      "Multiplicatively incommensurate vector and matrix: %s with %s."\
      % (len(self),u.rank)
    m=u.m; v_mtrc=V.__or__; v_=V._ # speedup: local lookup
    return v_(v_mtrc(self,r) for r in m)

  def __repr__(self):
    return self.__class__.__name__+"(%s)" % \
           (', '.join(repr(e) for e in self.__v),)

  def __rshift__(self,i):
    v=self.__v; l=len(v) # speedup: local lookup
    return self._(v[l-i:l]+v[0:l-i])

  def __str__(self):
    return "[ %s ]" % ('  '.join(str(e) for e in self.__v),)

  def __sub__(self,u): return self._(imap(_s,self.__v,u.v))

  def __xor__(self,u):
    try:
      c=(self|u)/(abs(self)*abs(u))
      return complex(c,sqrt(1-c*c))
    except ZeroDivisionError: return None

  @property
  def clean(self):
    """This vector 'V' with components rounded to 15 significant digits.
       Converts, for example:
         2.220446000000013e-16 -> 2.22044600000001e-16
         2.220446000000003e-16 -> 2.220446e-16
         0.5000000000000001    -> 0.5
         0.7499999999999999    -> 0.75
         5000000000000001.      -> 5000000000000000."""
    return self._(clean(e,15) for e in self.__v)

  @property
  def tidy(self):
    """This vector 'V' with V.clean first applied to components.  Then
       component elements 'e' such that 0.0 <= abs(e) <= 10**(-15) are
       replaced with 0.0 exactly.

       This is useful for printing values clearly and eliminating differences
       that are actually null.  CAUTION: when legitimate data is very small,
       it WILL be corrupted.  In such cases, use V.clean directly."""
    return self._(zero(clean(e,15),15) for e in self.__v)

  @property
  def v(self):
    """The tuple representation of this vector 'V'."""
    return self.__v

class M(object):
  """Class 'M' represents a Euclidean matrix.

     An 'M' is constructed from a sequence, list, or tuple of vectors 'V':
       a=M(V(1,2,3,4), V(4,3,2,1)); b=M([V(-1,-2,-3,-4), V(-4,-3,-2,-1)])
     'M' is an immutable type, and so may be used as a dictionary key.

     Operator overloading:
       '+'   : Matrix addition    <=>  a+b
       '-'   : Matrix subtraction <=>  a-b
       '*'   : Scalar multiplication <=>  pi*a , a*pi
       '/'   : Scalar division <=>  a/pi
       '-'   : Matrix negation <=>  -a
       '~'   : Matrix transposition <=>  ~a
       '|'   : Matrix product <=>  a|b
       'abs' : Volume of prism defined by rows = sqrt(det(a|~a))
       '=='  : Component-wise equality   <=>  a==a,a==b = True,False
       '!='  : Component-wise inequality <=>  a!=a,a!=b = False,True
       '<<'  : Rotate rows up   <=>  a<<1 == M(V(4,3,2,1),V(1,2,3,4))
       '>>'  : Rotate rows down <=>  a>>1 == M(V(4,3,2,1),V(1,2,3,4))
       'in'  : Containment or iteration <=>  v in a ; for v in a
       '[i]' : Extract a vector component <=>  a[1] == V(4,3,2,1)
       '[i,j]' : Extract a component <=>  a[1,2]==2 ; a[1]==V(4,3,2,1)
       '[(i,j):(k,l):(m,n)]' : Slicing returns the sub-matrix between
               (i,j) and (k,l) striking the (m,n) row and column <=>
                 a[::(1,2)] == M(V(1,2,4)); a[(1,2)::] == M(V(2,1))
                 a[:(1,2):] == a[:(1,2)] == M(V(1,2))
  """
  __slots__=['__m','__rank','__vol']

  def __init__(self,*n):
    self.__m=self.__rank=self.__vol=None # Realize __slots__
    if not n: return # class method '_' is caller: return empty instance
    if not isinstance(n[0],V): n=tuple(*n)
    if not all(isinstance(v,V) for v in n):
      raise ValueError,"Matrix components must be vectors of type V."
    l=len(n[0]); self.__rank=(len(n),l)
    if all(len(v)==l for v in n): self.__m=n
    else: raise ValueError,"Incommensurate vectors in matrix definition."

  @classmethod
  def _(cls,*n):
    """Internal init call speedup avoids integrity checks.
      'n' is a tuple, list, or generator of vector values."""
    inst=cls(); inst.__m=n=tuple(*n); inst.__rank=(len(n),len(n[0]))
    if len(inst.__m) > 0: return inst
    raise ValueError, "Report bug: null matrix definition encountered."

  def __abs__(self):
    if self.__vol is None: self.__vol=sqrt(abs(self.__mp(self)))
    return self.__vol

  def __add__(self,n):
    v_add=V.__add__ # speedup: local lookup
    return self._(v_add(u,v) for u,v in izip(self.__m,n.m))

  def __contains__(self,e): return e in self.__m

  def __div__(self,s):
    v_sc=V.__div__ # speedup: local lookup
    return self._(v_sc(v,s) for v in self.__m)

  def __richcmp__(self,n,op): # cython adaptation replacing python comparison ops
    v_rc=V.__richcmp__ # speedup: local lookup
    if op == 2: return all(v_rc(u,v,2) for u,v in izip(self.__m,n.m)) # __eq__
    return any(v_rc(u,v,3) for u,v in izip(self.__m,n.m)) # __ne__; op == 3

  def __getitem__(self,e):
    if isinstance(e,slice): return self.__getslice__(e)
    if not isinstance(e,tuple): return self.__m[e]
    try:
      i,j=e; return self.__m[i][j]
    except: raise ValueError, "At most two indices in a subscript."

  def __getslice__(self,s):
    ul,lr,mnr=s.start,s.stop,s.step
    m_=self._; v_=V._; m=self.__m # speedup: local lookup
    if ul is None: bi,bj=0,0
    else: bi,bj=ul
    if lr is None: ei,ej=self.rank # property lookup for subclasses
    else: ei,ej=lr
    if mnr is None: si,sj=ei,ej
    else: si,sj=mnr
    return m_(v_(u.v[bj:sj]+u.v[sj+1:ej]) for u in m[bi:si]+m[si+1:ei])

  def __hash__(self): return hash(self.__m)

  def __invert__(self):
    v_=V._ # speedup: local lookup
    return self._(v_(v) for v in izip(*self.__m))

  def __iter__(self): return iter(self.__m)

  def __len__(self): return self.__rank[0]

  def __lshift__(self,i):
    m=self.__m; l=len(m) # speedup: local lookup
    return self._(m[i:l]+m[0:i])

  def __mul__(self,s):
    v_sc=V.__mul__ # speedup: local lookup
    return self._(v_sc(v,s) for v in self.__m)

  __rmul__=__mul__

  def __neg__(self): return self._(-v for v in self.__m)

  def __mp(self,n): # internal speedup helper
    if self.rank[1]!=n.rank[1]: raise ValueError, \
      "Multiplicatively incommensurate matrices: %s with %s."\
      % (self.rank,(n.rank[1],n.rank[0]))
    if self.rank[0]==n.rank[0]: m_=SqM._
    else: m_=M._
    n=n.m; v_mtrc=V.__or__; v_=V._ # speedup: local lookup
    return m_(v_(v_mtrc(r,c) for c in n) for r in self.__m)

  def __or__(self,n):
    if not isinstance(n,M): return NotImplemented
    n=~n; return self.__mp(n)

  def __repr__(self):
    name=self.__class__.__name__; l=len(name)+1
    return '\n'+name+'(%s)\n' % \
           ((',\n'+l*' ').join(repr(v) for v in self.m),)

  def __rshift__(self,i):
    m=self.__m; l=len(m) # speedup: local lookup
    return self._(m[l-i:l]+m[0:l-i])

  def __str__(self): return '\n[ %s ]\n' % \
                      ('\n  '.join(str(v) for v in self.__m),)

  def __sub__(self,n):
    v_sub=V.__sub__ # speedup: local lookup
    return self._(v_sub(u,v) for u,v in izip(self.__m,n.m))

  @property
  def m(self):
    """The tuple value of the matrix represented by this instance."""
    return self.__m

  @property
  def rank(self):
    """The tuple (rows,columns) for this matrix instance."""
    return self.__rank

  @property
  def volume(self):
    """The volume of the prism defined by the rows of this matrix
       instance. If the matrix is square, this is a signed volume
       with volume=determinant(m)."""
    return abs(self)

  @property
  def clean(self):
    """This matrix 'M' with components rounded to 15 significant digits.
       Converts, for example:
         2.220446000000013e-16 -> 2.22044600000001e-16
         2.220446000000003e-16 -> 2.220446e-16
         0.5000000000000001    -> 0.5
         0.7499999999999999    -> 0.75
         5000000000000001.      -> 5000000000000000."""
    return self._(v.clean for v in self.__m)

  @property
  def tidy(self):
    """This matrix 'M' with M.clean first applied to components.  Then
       component elements 'e' such that 0.0 <= abs(e) <= 10**(-15) are
       replaced with 0.0 exactly.

       This is useful for printing values clearly and eliminating differences
       that are actually null.  CAUTION: when legitimate data is very small,
       it WILL be corrupted.  In such cases, use M.clean directly."""
    return self._(v.tidy for v in self.__m)

class SqM(M):
  """A 'SqM' is a square matrix comprising n row vectors of length n.

     Linearly independent vectors form a basis spanning a Euclidean
     n-space.  Dependent vectors represent a projection of the n-space
     into an m-dimensional sub-space.

     Operator overriding of super class M:
       'abs' : Signed volume of prism defined by rows = det(inst.m)
  """
  __slots__=['__inv']

  def __init__(self,*n):
    super(SqM,self).__init__(*n); self.__inv=None # Realize __slots__
    if not n: return # class method '_' is caller: return empty instance
    if self.rank[0]!=self.rank[1]: raise ValueError, \
      "A SqM must be a square matrix, self.rank is: %s" % (self.rank,)

  @classmethod
  def _(cls,*n):
    """Internal init call speedup avoids integrity checks.
      'n' is a tuple, list, or generator of vector values."""
    inst=cls(); inst._M__m=n=tuple(*n)
    inst._M__rank=(len(n),len(n[0]))
    if len(inst._M__m) > 0: return inst
    raise ValueError, \
      "Report bug: null square matrix definition encountered."

  def __alt(self,s,i): # alternating sum helper
    """Return -s if 'i' is odd or s if 'i' is even."""
    if i & 1: return -s # alt(s,i) call 2x faster than inline
    else: return s      # (cmp(1,i&1)-cmp(i&1,0))*s

  def __abs__(self):
    """Return the determinant of square matrix 'm'."""
    if self._M__vol is None:
      row=self.m[0]
      if self.rank==(1,1): self._M__vol=row[0]
      else:
        a,det=self.__alt,abs # speedup: local lookup
        self._M__vol=sum(a(e*det(self[::(0,c)]),c) \
                         for c,e in enumerate(row) if e!=0)
    return self._M__vol

  def __coftr(self): # inverse helper
    """Return the cofactor matrix transpose of square matrix 'm'."""
    a,det,v_ = self.__alt,abs,V._
    xr=xrange(self.rank[0]) # speedup: local lookup
    return self._(v_(a(det(self[::(i,j)]),i-j) for i in xr) for j in xr)

  @property
  def inverse(self):
    """The inverse of 'self.m' or 'None' if abs(m)==0."""
    vol=self._M__vol; inv=self.__inv # speedup: local lookup
    if vol is not None:
      if vol==0: return inv
      if inv is not None: return inv
    c=self.__coftr() # speedup: local lookup
    self._M__vol=vol=self.m[0]|(V._(v[0] for v in c))
    if vol != 0:
      vol=1./vol; self.__inv=inv=c*vol
      inv._SqM__inv=self; inv._M__vol=vol
    return inv

# Conversion notes: 20180406
#
# 'geometry.space.pyx' (this code) is identical to 'geometry.space.py'
# with the following exceptions:
#
# 1) 'sqrt' is imported from the 'math.h' C library
# 2) Comparison operator overloads conform to cython:
#       __eq__ and __ne__ are implemented by __richcmp__
# 3) Module tests have been moved to 'geometry.setup_space.py' since
#       this module must be compiled first.  The tests are identical
#       except they are prefaced with an import statements from the
#       newly created space.so file.
#
# In summary, the only optimization is 'sqrt' and that which may be
# derived by simplistic cythoning of a python module.
#
# It is probably worthwhile making a more refined translation ... someday.
