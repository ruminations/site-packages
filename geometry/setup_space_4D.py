# Elaborate cython setup version that works with libraries
# Run from the command line as:
#     python setup_space_4D.py build_ext --inplace
# builds an extension module .so file, overwriting any existing file.
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules=[Extension("space_4D", ["space_4D.pyx"], libraries =["m"])]
setup(name        = "space_4D_cython",\
      cmdclass    = {"build_ext" : build_ext },\
      ext_modules = ext_modules)
if __name__ == '__main__':
  try: import space_4D
  except ImportError:
    raise ImportError,"Extension failed to compile - omitting unit tests"
  print
  print "===Beginning module space_4D.pyx test suite==="
  from time import time
  from math import pi
  from space_4D import V,SqM
  t=time()
  a=V(1,2,3,4); b=V([4,3,2,1]); c=V((-1,0,1,0))
  print a.v,type(a.v)
  print repr(a),b,c,-c, len(c)
  print type(a),type(b),type(c),type(-c)
  print repr(a+b), repr(a-b)
  print type(a+b), type(a-b)
  print a|b, abs(c)
  print type(a|b), type(abs(c))
  print pi*c, c*pi
  print type(pi*c), type(c*pi)
  #print repr(c[1]), repr(c[:-1]), a[::2], repr(c[1:2])
  #print type(c[1]), type(c[:-1]), type(a[::2]), type(c[1:2]), type(a<<1)
  print a==b,(a<<1)==(a>>3),a<<1, a!=b,a!=a
  print a^b, b^a
  print "Iterator test:", tuple('<'+str(e)+'>' for e in a)
  from sys import exc_info as ex
  try: a+1
  except: print 'Caught:',ex()[1:2]
  try: a-1
  except: print 'Caught:',ex()[1:2]
  try: a+c
  except: print 'Caught:',ex()[1:2]
  try: a-c
  except: print 'Caught:',ex()[1:2]
  try: a*b
  except: print 'Caught:',ex()[1:2]
  try: a/b
  except: print 'Caught:',ex()[1:2]
  try: a|c
  except: print 'Caught:',ex()[1:2]
  print "===Vector tests complete==="
  g=SqM(V(1, 2, 3, 4), V([4, 3, 2, 1]),-a,-b)
  h=SqM([V(-1, -2, -3, -4), V([-4, -3, -2, -1]),a,b])
  d=SqM(a,b,-a,-b)
  print d.m, g.rank, h.rank
  print g,h,-h,~g, repr(g+h), repr(g-h)
  print type(g),type(h),type(-h),type(~g), type(g+h), type(g-h)
  print g<<1,g>>1
  print g==-h,g==h,g!=-h,g!=h,(g<<1)==-(g>>1)
  print g|~h, abs(g|~h)
  print type(g|~h), type(abs(g|~h))
  print pi*g,h*pi,len(g),len(~h)
  print type(pi*g),type(h*pi)
  #print repr(g[(1,2)]), h[1,2],type(g[(1,2)])
  #print g[:(1,2):],g[:(1,2)],g[(1,2):],g[(1,2)::],g[::(1,2)]
  #print map(type,(g[:(1,2):],g[:(1,2)],g[(1,2):],g[(1,2)::],g[::(1,2)]))
  print g|~g, abs(g)
  print type(g|~g), type(abs(g))
  print repr(g[1]),repr(b),type(g[1])
  print (g|d),(d|~g)
  print type(g|d),type(d|~g)
  print "===Singular Square Matrix tests complete==="
  s=SqM(V(3,-2,-1,4),V(2,1,2,-1),V(1,2,-3,-4),V(-2,1,-2,-1))
  print s, abs(s), s.inverse
  print type(s), type(abs(s)), type(s.inverse)
  print s|s.inverse, s.inverse.inverse, abs(s.inverse)
  print type(s|s.inverse), type(s.inverse.inverse), type(abs(s.inverse))
  print s|s
  print c|s, s|c, (c|(s|s)|c)
  print type(c|s), type(s|c), type((c|(s|s)|c))
  print c/.1, type(c/.1)
  for i in g: print type(i),
  print
  for i in s: print type(i),
  print
  print repr(s),s.m
  print "===Square Matrix tests complete==="
  d={a:'a',b:'b',c:'c'}
  for k,v in d.items(): print k,':',v
  print "===Dictionary hash tests complete==="
  print "elapsed time:",time()-t
  print "===Beginning numeric performance test==="
  t=time(); z=SqM(V(1.0,0.0,0.0,0.0),V(0.0,1.0,0.0,0.0), \
                  V(0.0,0.0,1.0,0.0),V(0.0,0.0,0.0,1.0))
  print z, type(z)
  for i in xrange(50): z=(z|s)
  print z,type(z)
  print "{!r}.clean ->".format(V(2.220446049250313e-16,0.5000000000000001,0.7499999999999999))
  print "    {!r}".format(V(2.220446049250313e-16,0.5000000000000001,0.7499999999999999).clean)
  print "{!r}.tidy  ->".format(V(2.220446049250313e-16,0.5000000000000001,0.7499999999999999))
  print "    {!r}".format(V(2.220446049250313e-16,0.5000000000000001,0.7499999999999999).tidy)
  print "elapsed time:",time()-t
  print "===All tests complete==="
