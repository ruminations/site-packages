#cython: nonecheck=False, boundscheck=False
"""
Package: geometry
Module : space_4D.pyx
Website: https://github.com/ruminations/site-packages
License: https://github.com/ruminations/Licenses#design-license
Initial Copyright 2013

This module is a cython optimization of space.py for the 4-dimensional
space that is so essential to graphics applications.  It preserves the
application programming interface of space.py to the extent possible.

 1) In preparation for conversion of the geometric subclass interfaces
    to a unified Grassman-Clifford algebra representation, vectors V,
    the basic mathematical tool, are limited to six components or
    dimensions, which is the number of components in a 4-D bivector.

 2) Consequently, only the nxn 'SqM' sub-type of 'M' is supported with n <= 4.

The package maintains support of arithmetic notation for operations on
4-vectors and 4x4 matrices.  The notation maintains augmentation by the
Dirac notation for matrix multiplication, e.g. (u|M|v) computes the
conjugation of matrix 'M' by vectors 'u' and 'v'.

Implementation notes:
  0) This is a complete reimplementation of the space.py interface as
     a cython generated C-language extension module.

  1) Data is internally represented by statically allocated arrays of
     the C floating point double type.
"""
__version__=20181211 # sanitize V initialization input
__version_history__=(20131106,20140526,20140530,20171207,20180413)
__all__=['V','SqM']
from algebra.mathx import zero,clean

cdef extern from "math.h":
  cdef double sqrt( double )

# Internal use utility functions

cdef short numeric(e): # internal numeric type check
  """Return 'True' if 'e' is an 'int', 'float', or 'complex' instance.
     Override in subclasses permitting extended numeric types."""
  if isinstance(e,float): return 1
  if isinstance(e,int): return 1
  if isinstance(e,complex): return 1
  return 0

cdef short commensurate(a,b):
  if len(a)==len(b): return 1
  return 0

cdef class V(object):
  """Class 'V' represents a Euclidean 4-vector.

     A vector 'V' may be constructed from four numbers, a list, or a tuple:
       a=V(1,2,3,4); b=V([-4.0,-3.0,-2.0,-1.0]); c=V((-1,0,1,0))
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
       '[i:j]' : Slicing <=>  a[1:3] == (2,3); a[::2] == (1,3)
       '^'   : a^b is the angle between a and b represented as the
               equivalent unit complex number in the plane of a and b.
                 a^b    = complex(cos(t),sin(t)) where
                 cos(t) = (a|b)/(abs(a)*abs(b))
                 sin(t) = sqrt(1-cos(t)*cos(t))
               Returns 'None' if abs(a) or abs(b) is zero.
  """
  cdef double __v[6] # ctype data values for this vector
  cdef double __abs # the Euclidean magnitude of this vector
  cdef short __l # number of valid components (dimension)
  cdef short __absflag # flag magnitude as valid
  cdef short __tflag # flag tuple as valid
  cdef tuple __t # python tuple type of the same values

  cdef double* __get_v(self): return &self.__v[0]
  cdef void*   __set_l(self,i): self.__l=i

  def __cinit__(self): #,*arg,**kwarg):# implicit
    cdef double* this
    cdef short i
    self.__abs=0.0;self.__l=0;  self.__absflag=0; self.__tflag=0
    this=V.__get_v(self)
    for i in range(6): this[i]=0.0

  def __init__(self,*u):
    cdef double* this
    cdef short i,t
    if not u: return # class method '_' is caller: return empty instance
    if not numeric(u[0]): u=tuple(*u)
    if len(u)>6: raise ValueError, \
      "A 4-D V is limited to 6 (bivector) components, len(self) == %s" % (len(u),)
    self.__t=u; self.__tflag=1; t=1; self.__l=len(u)
    for i in range(self.__l): t &= numeric(u[i])
    if t: # all numeric
      this=V.__get_v(self)
      for i in range(self.__l): this[i]=u[i]
    else: raise ValueError,"Vector elements must be numeric."

  cdef V _(self, double* u): #@classmethod replacement
    """Internal init call speedup avoids integrity checks.
       'u' is a double array pointer."""
    cdef double* that
    cdef short i
    inst=self.__class__(); that=V.__get_v(inst); V.__set_l(inst,self.__l)
    for i in range(self.__l): that[i]=u[i]
    return inst

  def __abs__(self):
    if self.__absflag: return self.__abs
    self.__abs=sqrt(self|self); self.__absflag=1; return self.__abs

  def __add__(self, V u):
    cdef double rval[6]
    cdef double* this
    cdef double* that
    cdef short i
    this=V.__get_v(self); that=V.__get_v(u)
    for i in range(len(self)): rval[i]=this[i]+that[i]
    if commensurate(self,u): return V._(self,rval)
    raise ArithmeticError, \
          "Addition of incommensurate vectors: %r,%r" % (self,u)

  def __contains__(self,double e):
    cdef double* this
    cdef short i
    this=V.__get_v(self)
    for i in range(len(self)):
      if this[i]==e: return True
    return False

  def __richcmp__(self, V u, short op):
    cdef object t,f
    cdef double* this
    cdef double* that
    cdef short i
    this=V.__get_v(self); that=V.__get_v(u)
    if op == 2: t,f=True,False # __eq__
    elif op == 3: t,f=False,True # __ne__
    else: raise NotImplementedError,"4-D rich comparison opcode: %s" % (op,)
    if not commensurate(self,u): return f
    for i in range(len(self)):
      if this[i] != that[i]: return f
    return t

  def __div__(self, double s):
    cdef double rval[6]
    cdef double* this
    cdef short i
    this=V.__get_v(self)
    for i in range(len(self)): rval[i]=this[i]/s
    return V._(self,rval)

  def __getitem__(self, short i):
    if i<0: i=i+len(self) # valid indices are -len <= i < len
    if i>=0 and i<len(self): return V.__get_v(self)[i]
    raise IndexError, "Index out of range; vector length == %r".format(len(self))

  def __getslice__(self, short i, int j):
    # j must be int to handle self[n:] => i=n, j=2147483647=0x7FFFFFFF
    return self.v[i:j]

  def __hash__(self): return hash(self.v)

  def __iter__(self): return iter(self.v)

  def __len__(self): return self.__l

  def __lshift__(self,short index):
    cdef double rval[6]
    cdef double* this
    cdef short i,t
    this=V.__get_v(self); t=len(self)-index
    for i in range(t): rval[i]=this[index+i]
    for i in range(index): rval[t+i]=this[i]
    return V._(self,rval)

  def __mul__(this, that):
    if isinstance(this,V): return V.__mul(this,that)
    if isinstance(that,V): return V.__mul(that,this) # cython defect
    raise TypeError,"Invalid scalar multiplication operand types: %r,%r" %\
                      (type(this),type(that))

  cdef V __mul(V self, double s): # __mul__ helper
    cdef double rval[6]
    cdef double* this
    cdef short i
    this=V.__get_v(self)
    for i in range(len(self)): rval[i]=this[i]*s
    return V._(self,rval)

  def __neg__(self):
    cdef double rval[6]
    cdef double* this
    cdef short i
    this=V.__get_v(self)
    for i in range(len(self)): rval[i] = -this[i]
    return V._(self,rval)

  def __or__(self,u):
    cdef double rval[6]
    cdef double t
    cdef double* this
    cdef double* that
    cdef short i,j
    this=V.__get_v(self)
    if isinstance(u,V):
      t=0.0; that=V.__get_v(u)
      for i in range(len(self)): t += this[i]*that[i]
      return t
    if not isinstance(u,SqM): raise TypeError, \
      "Invalid 4-D inner product operand type: %s" % (type(u),)
    that=SqM.__get_m(u)
    for j in range(len(self)):
      t=0.0;
      for i in range(u.rank[0]): t += this[i]*that[(i<<2)+j]
      rval[j]=t
    return V._(self,rval)

  cdef V __ror(self, SqM n): # __or__ helper
    cdef double rval[6]
    cdef double t
    cdef double* this
    cdef double* that
    cdef short i,j,row
    this=V.__get_v(self); that=SqM.__get_m(n)
    for i in range(n.rank[1]):
      t=0.0; row=i<<2
      for j in range(len(self)): t += that[row+j]*this[j]
      rval[i]=t
    return V._(self,rval)

  def __repr__(self):
    return self.__class__.__name__+"(%s)" % \
           (', '.join(repr(e) for e in self.v),)

  def __rshift__(self, short index):
    cdef double rval[6]
    cdef double* this
    cdef short i,t
    this=V.__get_v(self); t=len(self)-index
    for i in range(index): rval[i]=this[t+i]
    for i in range(t): rval[index+i]=this[i]
    return V._(self,rval)

  def __str__(self):
    return "[ %s ]" % ('  '.join(str(e) for e in self.v),)

  def __sub__(self, V u):
    cdef double rval[6]
    cdef double* this
    cdef double* that
    cdef short i
    this=V.__get_v(self); that=V.__get_v(u)
    for i in range(len(self)): rval[i]=this[i]-that[i]
    return V._(self,rval)

  def __xor__(self, V u):
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
    cdef double rval[6]
    cdef double* this
    cdef short i
    this=V.__get_v(self)
    for i in range(6): rval[i]=clean(this[i],15)
    return V._(self,rval)

  @property
  def tidy(self):
    """This vector 'V' with V.clean first applied to components.  Then
       component elements 'e' such that 0.0 <= abs(e) <= 10**(-15) are
       replaced with 0.0 exactly.

       This is useful for printing values clearly and eliminating differences
       that are actually null.  CAUTION: when legitimate data is very small,
       it WILL be corrupted.  In such cases, use V.clean directly."""
    cdef double rval[6]
    cdef double* this
    cdef short i
    this=V.__get_v(self)
    for i in range(6): rval[i]=zero(clean(this[i],15),15)
    return V._(self,rval)

  property v:
    def __get__(self):
      """The tuple representation of this vector 'V'."""
      cdef double* this
      cdef short i
      if self.__tflag: return self.__t
      l=[]; this=V.__get_v(self)
      for i in range(len(self)): l.append(this[i])
      self.__t=tuple(l); self.__tflag=1; return self.__t

cdef class SqM(object):
  """Class 'SqM' represents a 4x4 Euclidean matrix.

     Linearly independent vectors form a basis spanning a Euclidean
     4-space.  Dependent vectors represent a projection of the 4-space
     into an m-dimensional sub-space.

     An 'M' is constructed from a sequence, list, or tuple of vectors 'V':
       a=M(V(1,2,3,4), V(4,3,2,1),V(2,1,4,3), V(3,4,1,2))
       b=M([V(-1,-2,-3,-4), V(-4,-3,-2,-1),V(-2,-1,-4,-3), V(-3,-4,-1,-2)])
     'M' is an immutable type, and so may be used as a dictionary key.

     Operator overloading:
       '+'   : Matrix addition    <=>  a+b
       '-'   : Matrix subtraction <=>  a-b
       '*'   : Scalar multiplication <=>  pi*a , a*pi
       '/'   : Scalar division <=>  a/pi
       '-'   : Matrix negation <=>  -a
       '~'   : Matrix transposition <=>  ~a
       '|'   : Matrix product <=>  a|b
       'abs' : Signed 4-volume of prism defined by rows <=> det(a)
       '=='  : Component-wise equality   <=>  a==a,a==b = True,False
       '!='  : Component-wise inequality <=>  a!=a,a!=b = False,True
       '<<'  : Rotate rows up   <=>  a<<1 == M(V(4,3,2,1), ... ,V(1,2,3,4))
       '>>'  : Rotate rows down <=>  a>>1 == M(V(3,4,1,2), ... ,V(2,1,4,3))
       'in'  : Containment or iteration <=>  v in a ; for v in a
       '[i]' : Extract a vector component <=>  a[1] == V(4,3,2,1)
       '[i,j]' : Extract a component <=>  a[1,2]==2 ; a[1]==V(4,3,2,1)
       '[(i,j):(k,l):(m,n)]' : Slicing returns the sub-matrix between
               (i,j) and (k,l) striking the (m,n) row and column <=>
                 a[::(1,2)] == M(V(1,2,4)); a[(1,2)::] == M(V(2,1))
                 a[:(1,2):] == a[:(1,2)] == M(V(1,2))
  """

  cdef double __m[16] # ctype data values for this matrix
  cdef double __det # signed volume of the parallelepiped
  cdef SqM __inv # the inverse of this matrix
  cdef short __detflag # ==0 => invalid ==1 => valid
  cdef short __invflag # ==0 => invalid ==1 => valid ==2 => singular
  cdef short __tflag # flag tuple as valid
  cdef tuple __t # python tuple type of the same values

  cdef double* __get_m(self): return &self.__m[0]

  def __cinit__(self): #,*arg,**kwarg):# implicit
    self.__det=0.0; self.__detflag=0; self.__invflag=0; self.__tflag=0

  def __init__(self,*n):
    cdef short i,j,row
    if not n: return # class method '_' is caller: return empty instance
    if not isinstance(n[0],V): n=tuple(*n) # list or tuple argument
    if not all(isinstance(v,V) for v in n):
      raise ValueError,"SqM components must be 4-vectors of type V."
    if len(n)!=4: raise ValueError, \
      "A 4-D SqM must be 4x4, self.rank is: (%s,%s)" % (4,len(n))
    self.__t=n; self.__tflag=1
    for i in range(4):
      row=i<<2
      for j in range(4): self.__m[row+j]=n[i][j]

  cdef SqM _(self, double* n): #@classmethod replacement
    """Internal init call speedup avoids integrity checks.
       'n' is a double array pointer."""
    cdef double* that
    cdef short i
    inst=self.__class__(); that=SqM.__get_m(inst)
    for i in range(16): that[i]=n[i]
    return inst

  def __abs__(self):
    if self.__detflag: return self.__det
    self.inverse; return self.__det

  def __add__(self, SqM n):
    cdef double rval[16]
    cdef double* this
    cdef double* that
    cdef short i
    this=SqM.__get_m(self); that=SqM.__get_m(n)
    for i in range(16): rval[i]=this[i]+that[i]
    return SqM._(self,rval)

  def __contains__(self, V v):
    cdef short i
    for i in range(4):
      if self.m[i]==v: return True
    return False

  def __div__(self, double s):
    cdef double rval[16]
    cdef double* this
    cdef short i
    this=SqM.__get_m(self)
    for i in range(16): rval[i]=this[i]/s
    return SqM._(self,rval)

  def __richcmp__(self, SqM n, int op):
    cdef object t,f
    cdef double* this
    cdef double* that
    cdef short i
    this=SqM.__get_m(self); that=SqM.__get_m(n)
    if op == 2: t,f=True,False # __eq__
    elif op == 3: t,f=False,True # __ne__
    else: raise NotImplementedError,"4-D rich comparison opcode: %s" % (op,)
    for i in range(16):
      if this[i] != that[i]: return f
    return t

  def __getitem__(self,e):###
    if isinstance(e,slice): raise NotImplementedError #return self.__getslice__(e)
    if not isinstance(e,tuple): return self.m[e]
    try:
      i,j=e; return self.m[i][j]
    except: raise ValueError, "At most two indices in a subscript."
  """
  def __getslice__(self,s):###
    ul,lr,mnr=s.start,s.stop,s.step
    m_=self._; v_=V._; m=self.__m # speedup: local lookup
    if ul is None: bi,bj=0,0
    else: bi,bj=ul
    if lr is None: ei,ej=self.rank # property lookup for subclasses
    else: ei,ej=lr
    if mnr is None: si,sj=ei,ej
    else: si,sj=mnr
    return m_(v_(u.v[bj:sj]+u.v[sj+1:ej]) for u in m[bi:si]+m[si+1:ei])
  """
  def __hash__(self): return hash(self.m)

  def __invert__(self): # matrix transpose
    cdef double n[16]
    cdef double* this
    cdef short i,j,row
    this=SqM.__get_m(self)
    for j in range(4):
      row=j<<2
      for i in range(4): n[row+i]=this[(i<<2)+j]
    return SqM._(self,n)

  def __iter__(self): return iter(self.m)

  def __len__(self): return 4

  def __lshift__(self, short index):
    cdef double rval[16]
    cdef double* this
    cdef short i,t
    this=SqM.__get_m(self); index=index<<2; t=16-index
    for i in range(t): rval[i]=this[index+i]
    for i in range(index): rval[t+i]=this[i]
    return SqM._(self,rval)

  def __mul__(this, that):
    if type(this)==SqM: return SqM.__mul(this,that)
    if type(that)==SqM: return SqM.__mul(that,this) # cython defect
    raise TypeError,"Invalid scalar multiplication operand types: %r,%r" %\
                      (type(this),type(that))

  cdef SqM __mul(SqM self, double s): # __mul__ helper
    cdef double rval[16]
    cdef double* this
    cdef short i
    this=SqM.__get_m(self)
    for i in range(16): rval[i]=this[i]*s
    return SqM._(self,rval)

  def __neg__(self):
    cdef double rval[16]
    cdef double* this
    cdef short i
    this=SqM.__get_m(self)
    for i in range(16): rval[i] = -this[i]
    return SqM._(self,rval)

  def __or__(self,n):
    cdef double rval[16]
    cdef double t
    cdef double* this
    cdef double* that
    cdef short i,j,k,row
    if isinstance(n,V): return V.__ror(n,self)
    if not isinstance(n,SqM): raise TypeError, \
      "Invalid 4-D matrix product operand type: %s" % (type(n),)
    this=SqM.__get_m(self); that=SqM.__get_m(n)
    for i in range(4):
      row=i<<2
      for j in range(4):
        t=0.0
        for k in range(4): t+=this[row+k]*that[(k<<2)+j]
        rval[row+j]=t
    return SqM._(self,rval)

  def __repr__(self):
    name=self.__class__.__name__; l=len(name)+1
    return '\n'+name+'(%s)\n' % \
           ((',\n'+l*' ').join(repr(v) for v in self.m),)

  def __rshift__(self, short index):
    cdef double rval[16]
    cdef double* this
    cdef short i,t
    this=SqM.__get_m(self); index=index<<2; t=16-index
    for i in range(index): rval[i]=this[t+i]
    for i in range(t): rval[index+i]=this[i]
    return SqM._(self,rval)

  def __str__(self): return '\n[ %s ]\n' % \
                      ('\n  '.join(str(v) for v in self.m),)

  def __sub__(self, SqM n):
    cdef double rval[16]
    cdef double* this
    cdef double* that
    cdef short i
    this=SqM.__get_m(self); that=SqM.__get_m(n)
    for i in range(16): rval[i]=this[i]-that[i]
    return SqM._(self,rval)

  @property
  def clean(self):
    """This matrix 'SqM' with components rounded to 15 significant digits.
       Converts, for example:
         2.220446000000013e-16 -> 2.22044600000001e-16
         2.220446000000003e-16 -> 2.220446e-16
         0.5000000000000001    -> 0.5
         0.7499999999999999    -> 0.75
         5000000000000001.      -> 5000000000000000."""
    cdef double rval[16]
    cdef double* this
    cdef short i
    this=SqM.__get_m(self)
    for i in range(16): rval[i] = clean(this[i],15)
    return SqM._(self,rval)

  @property
  def tidy(self):
    """This matrix 'SqM' with SqM.clean first applied to components.  Then
       component elements 'e' such that 0.0 <= abs(e) <= 10**(-15) are
       replaced with 0.0 exactly.

       This is useful for printing values clearly and eliminating differences
       that are actually null.  CAUTION: when legitimate data is very small,
       it WILL be corrupted.  In such cases, use SqM.clean directly."""
    cdef double rval[16]
    cdef double* this
    cdef short i
    this=SqM.__get_m(self)
    for i in range(16): rval[i]=zero(clean(this[i],15),15)
    return SqM._(self,rval)

  property m:
    def __get__(self):
      """The tuple value of the matrix represented by this instance."""
      if self.__tflag: return self.__t
      inst=V(0.,0.,0.,0.)
      self.__t=(V._(inst,&self.__m[0]),V._(inst,&self.__m[4]), \
                V._(inst,&self.__m[8]),V._(inst,&self.__m[12]))
      self.__tflag=1; return self.__t

  property rank:
    def __get__(self):
      """The tuple (rows,columns) for this matrix instance."""
      return (4,4)

  property volume:
    def __get__(self):
      """The volume of the prism defined by the rows of this matrix
         instance. If the matrix is square, this is a signed volume
         with volume=determinant(m)."""
      return abs(self)

  property inverse:
    def __get__(self):
      if self.__invflag==1: return self.__inv
      if self.__invflag==2: return None # singular
      cdef double coftr[16]
      cdef double t0,t1,t2,t3,t4
      cdef double* m
      cdef short i
      m=SqM.__get_m(self)
      # compute 2x2 determinants in cofactor matrix
      # (36 products and 18 differences) interleaved with
      # assembly of the cofactor matrix transpose from 3x3 determinants
      # (48 products and 32 sum/differences)
      t0 = m[10]*m[15] - m[11]*m[14] # d10151114
      t1 = m[ 9]*m[15] - m[11]*m[13] # d_9151113
      t2 = m[ 9]*m[14] - m[10]*m[13] # d_9141013
      coftr[ 0]= m[5]*t0-m[6]*t1+m[7]*t2 # 4,4,4 : redundant usage counts
      coftr[ 1]= m[2]*t1-m[3]*t2-m[1]*t0 # 4,4,4
      t3 = m[ 8]*m[15] - m[11]*m[12] # d_8151112
      t4 = m[ 8]*m[14] - m[10]*m[12] # d_8141012
      coftr[ 4]= m[6]*t3-m[7]*t4-m[4]*t0 # 4,4,4
      coftr[ 5]= m[0]*t0-m[2]*t3+m[3]*t4 # 4,4,4
      # 5 temporaries required to here, 4 in use
      t0 = m[ 8]*m[13] - m[ 9]*m[12] # d_813_912
      coftr[ 8]= m[4]*t1-m[5]*t3+m[7]*t0 # 4,4,4
      coftr[ 9]= m[1]*t3-m[3]*t0-m[0]*t1 # 4,4,4
      # 5 temporaries required to here, 3 in use
      coftr[12]= m[5]*t4-m[6]*t0-m[4]*t2 # 4,4,4
      coftr[13]= m[0]*t2-m[1]*t4+m[2]*t0 # 4,4,4
      # 3 temporaries required to here, 0 in use
      t0 = m[ 6]*m[15] - m[ 7]*m[14] # d_615_714
      t1 = m[ 5]*m[15] - m[ 7]*m[13] # d_515_713
      t2 = m[ 5]*m[14] - m[ 6]*m[13] # d_514_613
      coftr[ 2]= m[1]*t0-m[2]*t1+m[3]*t2 # 2,2,2
      t3 = m[ 4]*m[15] - m[ 7]*m[12] # d_415_712
      t4 = m[ 4]*m[14] - m[ 6]*m[12] # d_414_612
      coftr[ 6]= m[2]*t3-m[3]*t4-m[0]*t0 # 2,2,2
      # 5 temporaries required to here, 4 in use
      t0 = m[ 4]*m[13] - m[ 5]*m[12] # d_413_512
      coftr[10]= m[0]*t1-m[1]*t3+m[3]*t0 # 2,2,2
      # 5 temporaries required to here, 3 in use
      coftr[14]= m[1]*t4-m[2]*t0-m[0]*t2 # 2,2,2
      # 3 temporaries required to here, 0 in use
      t0 = m[ 6]*m[11] - m[ 7]*m[10] # d_611_710
      t1 = m[ 5]*m[11] - m[ 7]*m[ 9] # d_511_7_9
      t2 = m[ 5]*m[10] - m[ 6]*m[ 9] # d_510_6_9
      coftr[ 3]= m[2]*t1-m[3]*t2-m[1]*t0 # 2,2,2
      t3 = m[ 4]*m[11] - m[ 7]*m[ 8] # d_411_7_8
      t4 = m[ 4]*m[10] - m[ 6]*m[ 8] # d_410_6_8
      coftr[ 7]= m[0]*t0-m[2]*t3+m[3]*t4 # 2,2,2
      # 5 temporaries required to here, 4 in use
      t0 = m[ 4]*m[ 9] - m[ 5]*m[ 8] # d_4_9_5_8
      coftr[11]= m[1]*t3-m[3]*t0-m[0]*t1 # 2,2,2
      # 5 temporaries required to here, 3 in use
      coftr[15]= m[0]*t2-m[1]*t4+m[2]*t0 # 2,2,2
      # 3 temporaries required to here, 0 in use
      # compute determinant of m
      t0=0.0
      for i in range(4): t0 += m[i]*coftr[i<<2]
      self.__det=t0; self.__detflag=1
      # check for a singular matrix
      if round(t0,15)==0.0:
        self.__invflag=2; return None
      # return the inverse
      t1=1.0/t0
      for i in range(16): coftr[i] *= t1
      self.__inv=SqM._(self,coftr)
      self.__inv.__inv=self; self.__inv.__det=t1
      self.__inv.__detflag=1; self.__inv.__invflag=1; self.__invflag=1
      return self.__inv
