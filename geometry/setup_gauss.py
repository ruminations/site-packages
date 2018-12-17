# Elaborate cython setup version that works with libraries
# Run from the command line as:
#     python setup_gauss.py build_ext --inplace
# builds an extension module .so file, overwriting any existing file.
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules=[Extension("gauss", ["gauss.pyx"], libraries =["m"])]
setup(name        = "gauss_cython",\
      cmdclass    = {"build_ext" : build_ext },\
      ext_modules = ext_modules)
if __name__ == '__main__':
  try: import gauss
  except ImportError:
    raise ImportError,"Extension failed to compile - omitting unit tests"
  print
  print "===Beginning module gauss.pyx test suite==="
  from time import time
  from math import pi,sqrt
  from gauss import Angle,C,U,V,O,Cline,Ccircle ###,Moebius # currently unimplemented
  t=time()
  print "init tests:\n"
  print "  predefined U,V,O:",map(repr,(U,V,O))
  a=C(1,2); b=C([2,1]); c=C((-1,1)); d=C(3); e=C(4.0); f=C(complex(4,5)); g=C(Angle(pi/3.0))
  print "  multi,list,tuple,int,float,complex repr:\n    ",map(repr,(a,b,c,d,e,f,g))
  print "  multi,list,tuple,int,float,complex type:\n    ",map(type,(a,b,c,d,e,f,g))
  print "\ninit error tests:\n"
  from sys import exc_info as ex
  print "  text argument:"
  try: C("text")
  except: print '    Caught:',ex()[1:2]
  print "  two char arguments:"
  try: C("tx")
  except: print '    Caught:',ex()[1:2]
  print "  boolean arguments:"
  try: C(True,False)
  except: print '    Caught:',ex()[1:2]
  print "  three arguments:"
  try: C(1,2,3)
  except: print '    Caught:',ex()[1:2]
  print "  single complex arguments:"
  try: C(1,complex(3,4))
  except: print '    Caught:',ex()[1:2]
  print "  multi complex arguments:"
  try: C(complex(1,2),complex(3,4))
  except: print '    Caught:',ex()[1:2]
  print "\ninner product tests:\n"
  print "  (a|b),(b|c):    ", map(repr,((a|b),(b|c)))
  z=complex(-1,1)
  print "  (z|b),(b|z) z==complex(-1,1):    ", map(repr,((z|b),(b|z)))
  print "  multi,list,tuple,int,float,complex abs:\n    ",map(abs,(a,b,c,d,e,f))
  print "\narithmetic tests:\n"
  s=['a+b','z+a','a+z','a-b','z-a','a-z','a*b','z*a','a*z','a/b','z/a','a/z']
  print "  a+b,z+a,a+z,a-b,z-a,a-z,a*b,z*a,a*z,a/b,z/a,a/z:\n    "
  for l,v in zip(s,map(repr,(a+b,z+a,a+z,a-b,z-a,a-z,a*b,z*a,a*z,a/b,z/a,a/z))):
    print "    ",l,v
  print "\nunary tests:\n"
  print "  (-a),(~a):    ", map(repr,(-a,~a))
  print "\nproperty tests:\n"
  print "  a.polar,a.xy,a.x,a.y:\n    ", map(repr,(a.polar,a.xy,a.x,a.y))
  print "  a.sqrt,a.unit,a.ortho,a.y:\n    ", map(repr,(a.sqrt,a.unit,a.ortho))
  print "\nAngle tests:"
  a=Angle(pi/3.0,"rad"); b=Angle(pi/6.0,"rad"); c=1.0/6.0; d=1.0/3.0
  print "  rad :  ",a.rad,b.rad
  print "  deg :  ",a.deg,b.deg
  print "  turn:  ",a.turn,b.turn
  print "  hemi:  ",a.hemi,b.hemi
  print "  unit:  ",a.unit
  print "  sum tests  :  ",(a+b).deg,(c+a).deg,(a+c).deg
  print "  diff tests :  ",(a-b).deg,(d-b).deg,(a-c).deg
  print "  prod tests :  ",(a*2).deg,(2*a).deg,(a*2.0).deg,(2.0*a).deg
  print "  div tests  :  ",(b/c).deg,(c/b).deg,(b/b).deg
  print "  unary tests:  ",(-a).deg,(~a).deg,(~(-a)).deg
  print "  a< b:  ",a< b
  print "  a<=c:  ",a<=c
  print "  b==c:  ",b==c
  print "  d!=a:  ",d!=a
  print "  a> b:  ",a> b
  print "  a>=c:  ",a>=c
  print "  copy a:",Angle(a).deg
  print "Cline tests:"
  k,l = Cline(C(1.0,1.0),C(2.0,2.0)),Cline(C(2.0,0.0),C(0.0,2.0))
  print "Clines k and l (repr):\n  {!r}\n  {!r}".format(k,l)
  print "Clines k and l (str) :\n  {!s}\n  {!s}".format(k,l)
  print "Their intersection (repr):\n  {!r}".format(k^l)
  print "Cline l's origin, U, and l(1.0):\n  {!r}\n  {!r}\n  {!r}". \
           format(l.origin,l.U,l(1.0))
  print "\n\n=== regression tests from 'geometry.figure.Circle': ==="
  # setup
  z0=C(0,2)
  z1=C(2,0)
  z2=C(2,2)
  p=[z0,z1,z2]
  l0=Cline(p[2],p[2]+z2)   # y=x
  l1=Cline(p[1],p[0]) # y=-x+1
  c0=Ccircle(*p)
  p=[C(-z+z2*(1-1/sqrt(2))) for z in [z0,z1,z2]]
  c1=Ccircle(*p)
  print "\n=== Circle tests: ===\n"
  print c0,'\n',c1
  def print_properties(cls,obj):
    import inspect
    d=dict( (k,v) for k,v in inspect.getmembers( \
        cls, predicate=inspect.isdatadescriptor) if not k.startswith('_') )
    for k,v in sorted(d.iteritems()):
      print "  {!r:<20} : {!r}".format(k,v.__get__(obj))
  print "c0 properties:"
  print_properties(Ccircle,c0)
  print "c1 properties:"
  print_properties(Ccircle,c1)
  print "test Ccircle.__contains__ and __call__:",c0.vertex in c0, c0(c0.vertex)
  print c0^c1
  print c1^c0
  print
  l=Cline(C(0.0,2.0),C(1.0,1.0))
  print l,c0^l # single point, coplanar tangent intersection
  print "\n\n=== Practical tests: ===\n"
  l=Cline(C(389.0, 0.0),C(389.0, 1.0))
  c=Ccircle( C(390.0, 0.00295241992118), \
            C(391.0, 0.00357727459669), \
            C(392.0, 0.00433214623759)  )
  #l=Cline(C(831.0, 0.0),C(831.0, 1.0,))
  #c=Ccircle( P(1.0, 828.0, 1.77308261662e-06), \
  #          C(829.0, 1.67308571253e-06), \
  #          C(830.0, 1.57919944105e-06)  )
  print l,'\n',repr(l)
  print c,'\n',repr(c)
  print "c properties:"
  print_properties(Ccircle,c)
  print
  print
  print l^c
  print
  print "magnitudes (intersections-center), radius:\n  ", \
         [abs(c.center-i) for i in l^c],repr(c.radius)
  print
  print
  print "repr C(1.0,2.220446049250313e-16), clean, tidy:\n  ",
  print repr(C(1.0,2.220446049250313e-16)),"\n  ",
  print repr(C(1.0,2.220446049250313e-16).clean),"\n  ",
  print repr(C(1.0,2.220446049250313e-16).tidy)
  print "str C(1.0,2.220446049250313e-16), clean, tidy:\n  ",
  print str(C(1.0,2.220446049250313e-16)),"\n  ",
  print str(C(1.0,2.220446049250313e-16).clean),"\n  ",
  print str(C(1.0,2.220446049250313e-16).tidy)
  print
  '''# currently unimplemented
  print "Moebius tests:\n"
  z,o,c = C(0,1),C(1),C(1,1)
  m=Moebius(z,o,c)
  print "Mapping of defining points zero, one, center:\n  {!r}, {!r}, {!r}". \
                                    format( repr(z(m)),repr(o(m)),repr(c(m)) )
  print "Mapping of points C(1,.5), C(.5,1), C(1.5,1):\n  {!r}, {!r}, {!r}". \
                 format( repr(m(C(1,.5))),repr(m(C(.5,1))),repr(m(C(1.5,1))) )
  '''
  print "===All tests complete==="
