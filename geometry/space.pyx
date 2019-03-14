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
__version__=20190309 # add mutable types & maintenance
__version_history__=(20130209,20131106,20180413)
__all__=['V','U','M','N','SqM','SqN']
cdef extern from "math.h":
    cdef double sqrt( double )
from algebra.mathx import zero,clean
from itertools import imap,izip # yields C code speedup w/ operators
from operator import mul as _m, eq as _eq, ne as _ne
from operator import neg as _n, add as _a, sub as _s
## TODO:
# 0) write meta-class that instantiates a square 'M' or 'N' as 'SqM' or 'SqN'

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
  __slots__=['_v','_abs']

  @staticmethod
  def numeric(e):
    """Return 'True' if 'e' is an 'int', 'float', or 'complex' instance.
       Override in subclasses permitting extended numeric types."""
    return isinstance(e,int) or isinstance(e,float) or \
           isinstance(e,complex)

  def __init__(self,*u):
    self._v=self._abs=None # Realize __slots__
    if not u: return # class method '_' is caller: return empty instance
    num=self.numeric # speedup: local lookup
    if not num(u[0]): u=tuple(*u)
    if all(num(e) for e in u): self._v=tuple(u)
    else: raise ValueError,"Vector elements must be numeric."

  @classmethod
  def _(cls,*u):
    """Internal init call speedup avoids integrity checks.
      'u' is a tuple, list, or generator of numeric values."""
    t= list if issubclass(cls,N) else tuple
    inst=cls(); inst._v=t(*u)
    if len(inst._v) > 0: return inst
    raise ValueError, "Report bug: null vector definition encountered."

  @classmethod
  def zero(cls,n):
    """Return a new 'n' dimensional zero vector.  Call on the class."""
    return cls([0.0]*n)

  @classmethod
  def basis(cls,i,n):
    """Return a new 'n' dimensional orthonormal basis vector with the
       'i'-th element set to unity.  Call on the class."""
    return cls(1.0 if i==j else 0.0 for j in range(n))

  def __abs__(self):
    if self._abs is None: self._abs=sqrt(self|self)
    return self._abs

  def __add__(self,u): return self._(imap(_a,self._v,u.v))

  def __contains__(self,e): return e in self._v

  def __richcmp__(self,u,op): # cython adaptation replacing python comparison ops
    if op == 2: return all(imap(_eq,self._v,u.v)) # __eq__
    return any(imap(_ne,self._v,u.v)) # __ne__; op == 3

  def __div__(self,s): return self._(a/s for a in self._v)

  def __getitem__(self,e): # all py3 & newer py2
    """Get a slice as a 'V' or a component numeric scalar."""
    if isinstance(e,slice): return self._(self._v[e]) # py3 or py2 slice [b:e:i]
    return self._v[e] # otherwise return a scalar component

  def __hash__(self): return hash(self._v)

  def __iter__(self): return iter(self._v)

  def __len__(self): return len(self._v)

  def __lshift__(self,i):
    v=self._v; l=len(v) # speedup: local lookup
    return self._(v[i:l]+v[0:i])

  def __mul__(self,s): return self._(a*s for a in self._v)

  __rmul__=__mul__

  def __neg__(self): return self._(imap(_n,self._v))

  def __or__(self,u):
    if isinstance(u,V): return sum(imap(_m,self._v,u.v))
    if not isinstance(u,M): return NotImplemented
    if len(self)!=u.rank[0]: raise ValueError, \
      "Multiplicatively incommensurate vector and matrix: {!r} with {!r}.". \
                                                  format(len(self),u.rank)
    m=(~u).m; v_mtrc=V.__or__; v_=self._ # speedup: local lookup
    return v_(v_mtrc(self,r) for r in m)

  def __ror__(self,u):
    if not isinstance(u,M): return NotImplemented
    if len(self)!=u.rank[1]: raise ValueError, \
      "Multiplicatively incommensurate vector and matrix: {!r} with {!r}.". \
                                                  format(len(self),u.rank)
    m=u.m; v_mtrc=V.__or__; v_=self._ # speedup: local lookup
    return v_(v_mtrc(self,r) for r in m)

  def __repr__(self):
    return self.__class__.__name__+"(%s)" % \
           (', '.join(repr(e) for e in self._v),)

  def __rshift__(self,i):
    v=self._v; l=len(v) # speedup: local lookup
    return self._(v[l-i:l]+v[0:l-i])

  def __str__(self):
    return "[ %s ]" % ('  '.join(str(e) for e in self._v),)

  def __sub__(self,u): return self._(imap(_s,self._v,u.v))

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
    return self._(clean(e,15) for e in self._v)

  @property
  def tidy(self):
    """This vector 'V' with V.clean first applied to components.  Then
       component elements 'e' such that 0.0 <= abs(e) <= 10**(-15) are
       replaced with 0.0 exactly.

       This is useful for printing values clearly and eliminating differences
       that are actually null.  CAUTION: when legitimate data is very small,
       it WILL be corrupted.  In such cases, use V.clean directly."""
    return self._(zero(clean(e,15),15) for e in self._v)

  @property
  def v(self):
    """The tuple representation of this vector 'V'."""
    return self._v


class U(V):
  """Updateable sub-class of vector 'V'.  'U' is a mutable version of 'V'
     much as a 'list' is a mutable version of 'tuple'.  The main consequence
     is that U.__abs__ computes the magnitude on every invocation."""

  def __init__(self,*u):
    super(U,self).__init__(*u)
    if not u: return # class method '_' is caller: return empty instance
    self._v=list(self._v)

  def __abs__(self):
    self._abs=sqrt(self|self)
    return self._abs

  def __setitem__(self,e,v): # all py3 & newer py2
    """Set a slice to an iterable 'v' or a component to a numeric 'v'."""
    num=self.numeric # speedup: local lookup
    if isinstance(e,slice): b=all(num(e) for e in v)
    else: b=num(v)
    if not b: raise ValueError,"Vector elements must be numeric."
    self._v[e]=v # py3 or single element or py2 slice [b:e:i]


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
     '[i:j]' : Extract a slice of rows <=> a[0:1] == M[V(1,2,3,4)]
     '[i,j]' : Extract a component <=>  a[1,2]==2
     '[(i,j):(k,l):(m,n)]' : Slicing returns the sub-matrix between
               (i,j) and (k,l) striking the (m,n) row and column <=>
                 a[::(1,2)] == M(V(1,2,4)); a[(1,2)::] == M(V(2,1))
                 a[:(1,2):] == a[:(1,2)] == M(V(1,2))"""
  __slots__=['_m','_rank','_vol']

  def __init__(self,*n):
    self._m=self._rank=self._vol=None # Realize __slots__
    if not n: return # class method '_' is caller: return empty instance
    if not isinstance(n[0],V): n=tuple(*n)
    if not all(isinstance(v,V) for v in n):
      raise ValueError,"Matrix components must be vectors of type V."
    l=len(n[0]); self._rank=(len(n),l)
    if all(len(v)==l for v in n): self._m=n
    else: raise ValueError,"Incommensurate vectors in matrix definition."

  @classmethod
  def _(cls,*n):
    """Internal init call speedup avoids integrity checks.
      'n' is a tuple, list, or generator of vector values."""
    inst=cls(); inst._m=n=tuple(*n); inst._rank=(len(n),len(n[0]))
    if len(inst._m) > 0: return inst
    raise ValueError, "Report bug: null matrix definition encountered."

  @classmethod
  def zero(cls,n,m):
    """Return a new 'n' row by 'm' column zero matrix.  Call on the class."""
    R=U if issubclass(cls,N) else V
    return cls(R.zero(m) for i in range(n))

  def __abs__(self):
    if self._vol is None: self._vol=sqrt(abs(self._mp(self))) # sqrt(abs(SqM))
    return self._vol

  def __add__(self,n):
    v_add=V.__add__ # speedup: local lookup
    return self._(v_add(u,v) for u,v in izip(self._m,n.m))

  def __contains__(self,e): return e in self._m

  def __div__(self,s):
    v_sc=V.__div__ # speedup: local lookup
    return self._(v_sc(v,s) for v in self._m)

  def __richcmp__(self,n,op): # cython adaptation replacing python comparison ops
    v_rc=V.__richcmp__ # speedup: local lookup
    if op == 2: return all(v_rc(u,v,2) for u,v in izip(self._m,n.m)) # __eq__
    return any(v_rc(u,v,3) for u,v in izip(self._m,n.m)) # __ne__; op == 3

  def __getitem__(self,e): # all py3 & newer py2
    if isinstance(e,slice): return self._getslice(e) # get a sub-array
    if isinstance(e,int): return self._m[e] # get a row vector
    try: # tuple => get a scalar component
      i,j=e; return self._m[i][j]
    except ValueError: raise ValueError, "Scalar subscripts require 2 integer indices."
    except TypeError: raise ValueError, "A scalar subscript is a 2 integer iterable."
    except IndexError: raise IndexError, \
      "Subscript {!r} exceeds mutable matrix rank: {!r}.".format(e,self.rank)

  def _getslice(self,s): # custom slicing helper
    t = s.start,s.stop,s.step # int => row slice; tuple => sub-matrix
    if any(isinstance(i,int) for i in t): # get slice of rows
      i,j,k = t; r,c = self.rank
      i = 0 if (i is None) else i
      j = r if (j is None) else j
      k = 1 if (k is None) else k
      q=j-i
      if q<=0:   p = 0
      elif k>=q: p = 1
      else:      p = q-(q//k)*(k-1) # len(range(i,j,k))
      P,Q = (SqM,SqN) if (p==c) else (M,N) # return proper type: sq/rect
      R = Q if isinstance(self,N) else P   # return proper type: mutable/immutable
      return R._(self._m[s])
    v_=U._ if isinstance(self,N) else V._
    ul,lr,mnr = t; m=self._m # speedup: local lookup
    bi,bj = (0,0)     if (ul is None) else ul
    ei,ej = self.rank if (lr is None) else lr # property lookup for subclasses
    r,c = ei-bi,ej-bj
    P,Q = (SqM,SqN) if (r==c) else (M,N) # return proper type: sq/rect
    R = Q if isinstance(self,N) else P   # return proper type: mutable/immutable
    if mnr is None: return R(v_(u.v[bj:ej]) for u in m[bi:ei]) # get sub-matrix
    si,sj = mnr # get sub-matrix strike row & column
    return R(v_(u.v[bj:sj]+u.v[sj+1:ej]) for u in m[bi:si]+m[si+1:ei])

  def __hash__(self): return hash(self._m)

  def __invert__(self):
    v_=V._ # speedup: local lookup
    return self._(v_(v) for v in izip(*self._m))

  def __iter__(self): return iter(self._m)

  def __len__(self): return self._rank[0]

  def __lshift__(self,i):
    m=self._m; l=len(m) # speedup: local lookup
    return self._(m[i:l]+m[0:i])

  def __mul__(self,s):
    v_sc=V.__mul__ # speedup: local lookup
    return self._(v_sc(v,s) for v in self._m)

  __rmul__=__mul__

  def __neg__(self): return self._(-v for v in self._m)

  def _mp(self,n): # internal speedup helper
    if self.rank[1]!=n.rank[1]: raise ValueError, \
      "Multiplicatively incommensurate matrices: {!r} with {!r}.".format( \
                                        self.rank,(n.rank[1],n.rank[0]) )
    mut = isinstance(self,N) or isinstance(n,N)
    if self.rank[0]==n.rank[0]: m_ = SqN._ if mut else SqM._
    else: m_ = N._ if mut else M._
    v_ = U._ if mut else V._ ; n=n.m; v_mtrc=V.__or__ # speedup: local lookup
    return m_(v_(v_mtrc(r,c) for c in n) for r in self._m)

  def __or__(self,n): # matrix|vector handled by V.__ror__
    if not isinstance(n,M): return NotImplemented
    n=~n; return self._mp(n)

  def __repr__(self):
    name=self.__class__.__name__; l=len(name)+1
    return '\n'+name+'(%s)\n' % \
           ((',\n'+l*' ').join(repr(v) for v in self.m),)

  def __rshift__(self,i):
    m=self._m; l=len(m) # speedup: local lookup
    return self._(m[l-i:l]+m[0:l-i])

  def __str__(self): return '\n[ %s ]\n' % \
                      ('\n  '.join(str(v) for v in self._m),)

  def __sub__(self,n):
    v_sub=V.__sub__ # speedup: local lookup
    return self._(v_sub(u,v) for u,v in izip(self._m,n.m))

  @property
  def m(self):
    """The tuple value of the matrix represented by this instance."""
    return self._m

  @property
  def rank(self):
    """The tuple (rows,columns) for this matrix instance."""
    return self._rank

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
    return self._(v.clean for v in self._m)

  @property
  def tidy(self):
    """This matrix 'M' with M.clean first applied to components.  Then
       component elements 'e' such that 0.0 <= abs(e) <= 10**(-15) are
       replaced with 0.0 exactly.

       This is useful for printing values clearly and eliminating differences
       that are actually null.  CAUTION: when legitimate data is very small,
       it WILL be corrupted.  In such cases, use M.clean directly."""
    return self._(v.tidy for v in self._m)


class N(M):
  """Updateable sub-class of matrix 'M'.  'N' is a mutable version of 'M'
     much as a 'list' is a mutable version of 'tuple'.  The main consequence
     is that N.__abs__ computes the volume on every invocation and 'N' is
     constructed from 'U' objects.."""

  def __init__(self,*n):
    super(N,self).__init__() # Realize __slots__
    if not n: return # class method '_' is caller: return empty instance
    if not isinstance(n[0],U): n=list(*n)
    if not all(isinstance(v,U) for v in n):
      raise ValueError,"Mutable matrix components must be vectors of type U."
    l=len(n[0]); self._rank=(len(n),l)
    if all(len(v)==l for v in n): self._m=n
    else: raise ValueError,"Incommensurate vectors in matrix definition."

  def __abs__(self):
    self._vol=sqrt(abs(self._mp(self))) # sqrt(abs(SqM))
    return self._vol

  def __setitem__(self,e,v):
    if isinstance(e,slice): self._setslice(e,v); return # set a sub-array
    if isinstance(e,int): # set a row vector
      if not isinstance(v,U): raise ValueError, \
        "Mutable matrix rows must be vectors of type 'U'."
      if len(v)!=self.rank[1]: raise ValueError, \
        "Set vector length: {} incommensurate with matrix row length: {}.". \
          format(len(v),self.rank[1])
      self._m[e]=v; return
    if not V.numeric(v): raise ValueError,"Mutable matrix entries must be numeric."
    try: # tuple => set a scalar component
      i,j=e; self._m[i][j]=v
    except ValueError: raise ValueError, "Scalar subscripts require 2 integer indices."
    except TypeError: raise ValueError, "A scalar subscript is a 2 integer iterable."
    except IndexError: raise IndexError, \
      "Subscript {!r} exceeds mutable matrix rank: {!r}.".format(e,self.rank)

  def _setslice(self,s,v): # custom slicing helper
    if not isinstance(v,M): raise ValueError,"Setting matrix must be sub-type of 'M'."
    ul,lr,mnr=s.start,s.stop,s.step # int => row slice; tuple => sub-matrix
    if isinstance(ul,int): # set slice of rows
      if not v.rank[1]==self.rank[1]: raise ValueError, \
        "Set matrix row length: {} incommensurate with matrix row length: {}.". \
          format(v.rank[1],self.rank[1])
      self._m[s]=v; return
    v_=U._; m_=self._; m=self._m # speedup: local lookup
    bi,bj = (0,0)     if (ul is None) else ul
    ei,ej = self.rank if (lr is None) else lr # property lookup for subclasses
    if mnr is None: # set sub-matrix
      for i in range(bi,ei): m[i][bj:ej]=v[i][bj:ej]; return # sub-matrix
    si,sj = mnr #set sub-matrix omit row & column
    for i in range(bi,ei): m[i][bj:sj],m[i][sj+1:ej]=v[i][bj:sj],v[i][sj+1:ej]; return


class SqM(M):
  """A 'SqM' is a square matrix comprising n row vectors of length n.

     Linearly independent vectors form a basis spanning a Euclidean
     n-space.  Dependent vectors represent a projection of the n-space
     into an m-dimensional sub-space.

     Operator overriding of super class M:
       'abs' : Signed volume of prism defined by rows = det(inst.m)"""
  __slots__=['_inv']

  def __init__(self,*n):
    super(SqM,self).__init__(*n); self._inv=None # Realize __slots__
    if not n: return # class method '_' is caller: return empty instance
    if self.rank[0]!=self.rank[1]: raise ValueError, \
      "A SqM must be a square matrix, self.rank is: %s" % (self.rank,)

  @classmethod
  def _(cls,*n):
    """Internal init call speedup avoids integrity checks.
      'n' is a tuple, list, or generator of vector values."""
    inst=cls(); inst._m=n=tuple(*n)
    inst._rank=(len(n),len(n[0]))
    if len(inst._m) > 0: return inst
    raise ValueError, \
      "Report bug: null square matrix definition encountered."

  @classmethod
  def zero(cls,n):
    """Return a new 'n' by 'n' zero matrix.  Call on the class."""
    R=U if issubclass(cls,SqN) else V
    return cls(R.zero(n) for i in range(n))

  @classmethod
  def identity(cls,n):
    """Return a new 'n' dimensional identity matrix.  Call on the class."""
    R=U if issubclass(cls,SqN) else V
    return cls(R.basis(i,n) for i in range(n))

  def _alt(self,s,i): # alternating sum helper
    """Return -s if 'i' is odd or s if 'i' is even."""
    if i & 1: return -s # alt(s,i) call 2x faster than inline
    else: return s      # (cmp(1,i&1)-cmp(i&1,0))*s

  def __abs__(self):
    """Return the determinant of square matrix 'm'."""
    if self._vol is not None: return self._vol
    row=self._m[0]
    if self.rank==(1,1): self._vol=row[0]; return self._vol
    a,det=self._alt,abs # speedup: local lookup
    self._vol=sum(a(e*det(self[::(0,c)]),c) for c,e in enumerate(row) if e!=0)
    return self._vol

  def _coftr(self): # inverse helper
    """Return the cofactor matrix transpose of square matrix 'm'."""
    a,det,v_ = self._alt,abs,V._
    xr=xrange(self.rank[0]) # speedup: local lookup
    return self._(v_(a(det(self[::(i,j)]),i-j) for i in xr) for j in xr)

  @property
  def inverse(self):
    """The inverse of 'self.m' or 'None' if abs(m)==0."""
    vol=self._vol; inv=self._inv # speedup: local lookup
    if vol is not None:
      if vol==0: return inv
      if inv is not None: return inv
    c=self._coftr() # speedup: local lookup
    self._vol=vol=self.m[0]|(V._(v[0] for v in c))
    if vol != 0:
      vol=1./vol; self._inv=inv=c*vol
      inv._inv=self; inv._vol=vol
    return inv


class SqN(SqM,N):
  """Updateable sub-class of 'SqM'.  'SqN' is a mutable version of 'SqM'
     much as a 'list' is a mutable version of 'tuple'.  The main consequence
     is that N.__abs__ computes the determinant on every invocation and 'SqN'
     is constructed from 'U' objects."""

  def __init__(self,*n):
    super(SqN,self).__init__(*n); self.__inv=None # Realize __slots__
    if not n: return # class method '_' is caller: return empty instance
    if self.rank[0]!=self.rank[1]: raise ValueError, \
      "A SqN must be a square matrix, self.rank is: %s" % (self.rank,)

  def __abs__(self):
    """Return the determinant of square matrix 'm'."""
    row=self._m[0]
    if self.rank==(1,1): self._vol=row[0]; return self._vol
    a,det=self._alt,abs # speedup: local lookup
    self._vol=sum(a(e*det(self[::(0,c)]),c) for c,e in enumerate(row) if e!=0)
    return self._vol

# Implementation notes: 20190312
#
# There are now significant differences between this .pyx module and the
# corresponding .py module: primarily in the form of new mutable content.
#
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
