"""
Package: data
Module : sort.py
Website: https://github.com/ruminations/site-packages
License: https://github.com/ruminations/Licenses#design-license
Initial Copyright 2013

Class implementation of a list stored binary tree priority queue data
structure 'PQ', sometimes called a 'heap'.  The term 'heap', however, can
be confused with a dynamically allocated memory pool.  A 'PQ' algorithm is
equivalent to an n*log(n) sorting algorithm computationally, often called
'heapsort'.  The terminology here clarifies the use cases.

Elements are accessed via the root node.  The root node always holds an
extreme of all items currently in the 'PQ' under its comparison.  The
default comparison is '<=', and in that case, the extreme is the <least>
item.  Other 'PQ' nodes are inaccessible to the user.

Insertion and removal are 'stable', that is, equal priority items are
removed in the order they are inserted: first in, first out.

Usage:

q=PQ()        : Create an empty 'PQ' with comparison '<='.
q.set(iter)   : Clear the 'PQ' and O(n) set from iterable.
len(q)        : Gives the number of data items in the 'PQ'.
str(q)        : Returns the string representation of the 'PQ' array.
repr(q)       : Returns the repr() representation of the 'PQ'.
q.root=item   : Push an item onto the 'PQ' and shift it into place.
item=q.root   : Pull the root and shift a replacement into the root.
                Iteration is supported via this property and __iter__.
item=q.head   : Examine the root item (<least>) without altering the 'PQ'.
item=q.tail   : Examine the <greatest> item without altering the 'PQ'.
item=q(item)  : Size invariant 'PQ' push followed by pull.  NOP when
                the pushed item is <least>.
item=q[item]  : Size invariant 'PQ' pull followed by push.  Guarantees
                returning new data to process.
q=PQ(relop)   : Create a 'PQ' using 'relop' as the comparison.  'relop'
                takes two 'PQ' items as operands and returns a boolean:
                  relop(op1,op2)==True  => op1 moves closer to the root
                  relop(op1,op2)==False => op2 moves closer to the root
                The operator module is a useful source for comparison
                functions, or a custom function can be coded.
bool(q)       : 'False' when empty, 'True' otherwise.  Application of
                'bool' is automatic in, for example, 'while q: q.root',
                which empties 'q' without an 'IndexError'.

Acknowledgement:
  inspired by Kevin O'Connor's python standard library 'heapq.py' module.

Implementation Notes:

  0) Comments assume '<=' is the relop.  The term <least> is used where,
     with a different ordering, it would be replaced by another term,
     e.g. <greatest> when '>=' is the relop.

  1) A 'PQ' is an array representation of a binary tree.  The two nodes
     stemming from the index 'i' node are stored at indices 2*i+1 and
     2*i+2.  As is usual in python, indices range from 0 to n-1, where
     n is the number of nodes in the tree (array).  It can be seen there
     will be no conflicting storage indices.

  2) By introducing a comparison relation, for example '<=' and ensuring
     a heap invariant, e.g. h[i] <= h[2*i+1] and h[i] <= h[2*i+2], is
     True for all 'i', an ordered data structure is created.  O(log2(n))
     operations maintain the invariant for each insertion or removal of
     data.  Index references beyond the end of the array are treated as
     an extreme of the comparison, e.g. in the default case as infinite.

  3) User defined data may be complex - e.g. a class with properties - it
     is simply a python reference, which in the simplest case refers to a
     python numeric value.  The user may therefore easily define a 'relop'
     and easily prioritize (sort) arbitrarily complex data.

  4) The shifting code is very brief.  Therefore it is replicated
     directly in the public methods.  This avoids function call overhead
     and redundant fetching of data the method already has available.

  5) Shifting after a pull shifts the entire chain of <least> values
     forward, then looks for a place to put the former list end.  This
     will usually be near the end of the <least> chain.  Additionally,
     the operation moves <least> values up, optimizing future pulls.

  6) *2 and /2 operations are implemented with bit shifts <<1 and >>1.

  7) 'PQ.set(iterable)' requires O(n) time - faster than the O(n*log2(n))
     that successive root accesses requires.
"""
__all__=('PQ','RelOp')
__version__=20181125 # prepare for publication
__version_history__=(20130212,20130611)
## STATUS
# 0) implemented: PQ
## TODO:
# 0) as always, writing a doctest suite.
# 1) look into adding inckey and deckey: followed by shiftup/shiftdown respectively
# 2) look into adding remove: deckey(-inf) followed by pull root
# 3) look into heap on top PQ

class RelOp(object):
  """Container class for priority calculations for a 'PQ'."""
  __slots__=('_k','_r') # key function, relation operator
  from operator import le
  def __init__(self,key=(lambda d:d), relop=le):
    """'key(data)' returns a value that relop compares to determine
       the relative priority of two data items:
         relop( key(op1),key(op2) ) == True  => op1 moves closer to root
         relop( key(op1),key(op2) ) == False => op2 moves closer to root
       By default, the data is compared directly using '<='."""
    self._k,self._r = key,relop
  def __call__(self,op1,op2):
    """Return 'relop(key(op1),key(op2))'."""
    return self._r( self._k(op1),self._k(op2) )
  def __repr__(self): return "RelOp({!r},{!r})".format(self._k,self._r)
  @property
  def key(self):
    """Get or set the priority of this 'PQentry'."""
    return self._k(self._d)
  @key.setter
  def key(self,priority): return self._m


class PQ(object):
  """Wraps a binary tree representation stored in a list ordered by a
     comparison operation.  The default comparison is '<='."""
  __slots__=('_a','_relop')

  def __init__(self,relop=RelOp()):
    """Create an empty 'PQ' and assign 'relop' as for priority management."""
    self._relop=relop; self._a=[]

  def __call__(self,item):
    """Efficient push followed by a pull.  'PQ' size remains invariant.
       NOP when the pushed 'item' is <least>."""
    a=self._a; op=self._relop; n=len(a) # speedup: local lookup
    if n==0: return item
    if op(item,a[0]): return item
    ret=a[0]; stem,i=1,0
    while stem < n: # shift chain of <least> values
      twin=stem+1
      if twin<n and op(a[twin],a[stem]): stem=twin
      a[i],i,stem = a[stem],stem,(stem<<1)+1
    while i > 0: # shift back along chain of root nodes
      node=(i-1)>>1; value=a[node]
      if op(value,item): break
      a[i],i = value,node # item is <least>
    a[i]=item # drop target in its place
    return ret

  def __contains__(self,item): raise NotImplementedError

  def __getitem__(self,item):
    """Efficient pull followed by a push.  'PQ' size remains invariant.
       Guarantees getting new data.  Useful for restoring a processed
       'item' to the 'PQ' and returning the next 'item' to process."""
    a=self._a; op=self._relop; n=len(a) # speedup: local lookup
    if n==0: raise IndexError, "Pull from empty 'PQ': {!r}.".format(self)
    ret=a[0]; stem,i=1,0
    while stem < n: # shift chain of <least> values
      twin=stem+1
      if twin<n and op(a[twin],a[stem]): stem=twin
      a[i],i,stem = a[stem],stem,(stem<<1)+1
    while i > 0: # shift back along chain of root nodes
      node=(i-1)>>1; value=a[node]
      if op(value,item): break
      a[i],i = value,node # item is <least>
    a[i]=item # drop target in its place
    return ret

  def __delitem__(self,item): #### needs work
    """Asynchronously remove the first instance of 'item' and re-shuffle
       the 'PQ'.  If 'item' is not in the 'PQ', raise 'ValueError'."""
    # conceptually, set priority -infinity, shuffle up, self.root
    # no need to actually set priority: simple omit comparisons in shuffle
    try: i=item.index(); stem=(i<<1)+1
    except ValueError: raise PriorityError, \
      "Can not delete: {!r} is not in 'PQ' {!r}.".format(item,self)
    a=self._a; op=self._relop; n=len(a) # speedup: local lookup
    while i > 0: # shift back along chain of root nodes
      node=(i-1)>>1; value=a[node]
      if op(value,item): break
      a[i],i = value,node # item is <least>
    del(self.root)

  def __iter__(self): return self

  def __len__(self): return len(self._a)

  def __repr__(self): return "PQ({!r})".format(self._relop)

  def __str__(self): return str(self._a)

  def next(self):
    """Pull the next <least> item."""
    if bool(self): return self.root
    raise StopIteration

  def _pull(self):
    """Pull root of the 'PQ', maintaining the 'PQ' invariant."""
    a=self._a; op=self._relop; n=len(a) # speedup: local lookup
    if n==0: raise IndexError, "Pull from empty 'PQ': {!r}.".format(self)
    if n==1: return a.pop()
    ret,item = a[0],a.pop()
    n-=1; stem,i=1,0
    while stem < n: # shift chain of <least> values
      twin=stem+1
      if twin<n and op(a[twin],a[stem]): stem=twin
      a[i],i,stem = a[stem],stem,(stem<<1)+1
    while i > 0: # shift back along chain of root nodes
      node=(i-1)>>1; value=a[node]
      if op(value,item): break
      a[i],i = value,node # item is <least>
    a[i]=item # drop target in its place
    return ret

  def _push(self,item=None):
    """Push 'item' onto the 'PQ', maintaining the 'PQ' invariant."""
    a=self._a; op=self._relop; i=len(a) # speedup: local lookup
    a.append(None) # make space in the list
    while i > 0: # shift back along chain of root nodes
      node=(i-1)>>1; value=a[node]
      if op(value,item): break
      a[i],i = value,node # item is <least>
    a[i]=item # drop item in its place

  root=property(_pull,_push,None,"Root node access.")

  def set(self,it):
    """Clear the 'PQ' and set with iterable 'it' values in O(n) time."""
    a=self._a=list(it); op=self._relop; n=len(a)
    for base in reversed(xrange(n>>1)): # impose heap invariant
      i=base; item=a[base]; stem=(base<<1)+1
      while stem < n: # shift chain of <least> values
        twin=stem+1
        if twin<n and op(a[twin],a[stem]): stem=twin
        a[i],i,stem = a[stem],stem,(stem<<1)+1
      while i > base: # shift back along chain of root nodes
        node=(i-1)>>1; value=a[node]
        if op(value,item): break
        a[i],i = value,node # item is <least>
      a[i]=item # drop target in its place

  @property
  def head(self):
    """Examine the root item (<least>) without altering the 'PQ'.
       Returns 'None' for an empty 'PQ'.  O(1) time."""
    try: return self._a[0]
    except IndexError: return None

  @property
  def tail(self):
    """Examine the <greatest> item without altering the 'PQ'.
       Returns 'None' for an empty 'PQ'.  O(log2(n)) time."""
    if not self: return None


  @property
  def relop(self):
    """The priority management class used to order this 'PQ'."""
    return self._relop

if __name__ == "__main__":
  print "No syntax errors."
  print "Module tests:"
  q=PQ()
  print "'PQ' relop,root:",q.relop,',',q.head
  data=[1, 3, 5, 7, 9, 2, 4, 6, 8, 0]
  for item in data:
    q.root=item
    print "'PQ':",q
  print "'PQ' relop,root:",q.relop,',',q.head
  sort=[]
  while len(q):
    sort.append(q.root)
    print "sort:",sort,"  'PQ':  ",q
  data=[9, 8, 7, 7, 9, 2, 4, 6, 8, 0]
  for item in data:
    q.root=item
    print "'PQ':",q
  for i in xrange(10,20): print q(i),
  print "'PQ':",q
  print "'PQ' relop,root:",q.relop,',',q.head
  from operator import ge
  p=PQ(ge)
  while q.head:
    p.root=q.root
    print "p:",p,"  q: ",q
  print "'PQ' relop,root:",p.relop,',',q.head
  print bool(p),bool(q) # calls len()
  q.set(data)
  print "q after setting:",q
  print bool(p),bool(q) # calls len()
  for i in xrange(5): print q[7],
  print "q:",q
  for e in q: print e,
  print "\np:",repr(p)
  for e in p:
    sort.append(e)
    print "sort:",sort
    print "p:",p
  try:
    print q.root # raise IndexError
    print "Failed to raise IndexError"
  except IndexError: print "Successfully raised IndexError."
  try:
    print p, 5 in p, p # raise NotImplementedError
    print "Failed to raise NotImplementedError"
  except NotImplementedError:
    print "Successfully raised NotImplementedError"
  print "Tests completed successfully."
