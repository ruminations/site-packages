# Elaborate cython setup version that works with libraries
# Run from the command line as:
#     python setup_mathx.py build_ext --inplace
# builds an extension module .so file, overwriting any existing file.
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules=[Extension("mathx", ["mathx.pyx"], libraries =["m"])]
setup(name        = "mathx_cython",\
      cmdclass    = {"build_ext" : build_ext },\
      ext_modules = ext_modules)
if __name__ == '__main__':
  try: import mathx
  except ImportError:
    raise ImportError,"Extension failed to compile - omitting unit tests"
  print "version:",mathx.__version__
  from mathx import cbrt,qntrt,ferf,fcdf ##,erf_exec_time_compare # excised after testing
  print "===Beginning module mathx.pyx test suite==="
  from time import time
  print "===Beginning numeric performance tests==="
  t=time();
  print "\nCube root tests:"
  for i in range(-9,10):
    r=cbrt(i)
    print "  cbrt({!r}) == {!r}; cubed == {!r}".format(i,r,r*r*r)
  r=cbrt(100)
  print "  cbrt({!r}) == {!r}; cubed == {!r}".format(100,r,r*r*r)
  print "elapsed time:",time()-t
  t=time()
  print "\nFifth root tests:"
  for i in range(-9,10):
    r=qntrt(i)
    print "  qntrt({!r}) == {!r}; ^5 == {!r}".format(i,r,r*r*r*r*r)
  r=qntrt(100)
  print "  qntrt({!r}) == {!r}; ^5 == {!r}".format(100,r,r*r*r*r*r)
  print "elapsed time:",time()-t
  from math import erf
  t=time()
  print "\nError function tests:"
  for i in range(-9,10): print "  ",i,ferf(i), ferf(i)-erf(i)
  print "elapsed time:",time()-t
  t=time()
  print "\nCumulative distribution function tests:"
  for i in range(-9,10): print "  ",i,fcdf(i)
  print "elapsed time:",time()-t
  ## excised after testing
  ##print "\n1,000,000 iterations time(c_erf)/time(ferf): {}".format(erf_exec_time_compare(2.0))
  print "\nfactorial, permutation, and combination tests:"
  from mathx import factorial,permutation,combination
  print map(repr,[factorial(0),factorial(1),factorial(3),factorial(-1)])
  print map(repr,[permutation(0,0),permutation(3,0),permutation(0,3)])
  print map(repr,[permutation(-1,0),permutation(3,-1),permutation(-1,-3)])
  print map(repr,[permutation(3,1),permutation(3,2),permutation(3,3)])
  print map(repr,[combination(0,0),combination(3,0),combination(0,3)])
  print map(repr,[combination(-1,0),combination(3,-1),combination(-1,-3)])
  print map(repr,[combination(3,1),combination(3,2),combination(3,3)])
  from mathx import len,even,odd
  class len_test(object):
    """test arbitrary object return for mathx.len"""
    def __len__(self): return self
  print "mathx.len test:", len(len_test())
  print "even odd test on (3,4)",map(even,(3,4)),map(odd,(3,4))
  print "\n===All tests complete==="
