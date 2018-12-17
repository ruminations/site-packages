#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Package: data
Module : catalog.py
Website: https://github.com/ruminations/site-packages
License: https://github.com/ruminations/Licenses#design-license
Initial Copyright 2017

Defines:
1) 'DataTable'
  A metaclass enforcing an immutable singleton simplified dictionary.
  This provides a degree of type safety against accidental alteration of
  'constant' data tables by client modules.  It also prevents redundant
  memory allocation when accessing invariant data.
2) 'Catalog'
  A base class implementing an empty 'DataTable' object type.
3) 'AWG'
  A readonly 'Catalog' of electrical wire gauge sizes and properties.
4) 'Paper'
  A readonly 'Catalog' of standard paper and display sizes.
5) 'SFM'
  A readonly 'Catalog' of typical material surface feet per minute.
6) 'Drill'
  A readonly 'Catalog' of standard drill bit sizes.
7) 'Greek'
  A readonly 'Catalog' of Greek letter names.
8) 'Physics'
  A readonly 'Catalog' of SI unit physical constants.
9) 'WhitePoint'
  A readonly 'Catalog' of standard illumination CIE-XYZ constants.
10)'WebNameToColor'
  A readonly 'Catalog' of standard web color names to SRGB values.
11)'WebColorToName'
  A readonly 'Catalog' of standard web SRGB color values to names.
12)'WebNameToCamel'
  A readonly 'Catalog' of standard web color names to camel case.

Usage:

0)  Primary data values may be accessed calling the meta-method 'value'
    on the class:
      AWG.value('12')    yields  floating point 0.0808 == 12 gauge wire size
      Greek.value('α')   yields  'alpha': the name of the unicode character
      Physics.value('c') yields  299792458.0: the speed of light in SI units

1)  An entire record may be accessed by the meta-method '__getitem__'
    on the class:
      AWG['12']    yields  0.0808
      Greek['α']   yields  'alpha'
      Physics['c'] yields  (299792458.,'v',"speed of light")

2)  A reference to the singleton is returned by the meta-property 'table'
    called on the class, and it may be used like an ordinary dictionary
    instance, except any attempt to alter data will raise an error:
      AWG.table['12']    yields  0.0808
      Greek.table['α']   yields  'alpha'
      Physics.table['c'] yields  (299792458.,'v',"speed of light")

3)  The meta-method '__call__' also returns the singleton reference:
      AWG()['12']    yields  0.0808
      Greek()['α']   yields  'alpha'
      Physics()['c'] yields  (299792458.,'v',"speed of light")

4)  The meta property 'keys' returns a tuple of valid key values.

5)  The meta-property 'choices' returns a sorted list of tuples
      (<descriptive string>,<key>)
    when the data record is a tuple; otherwise returns a sorted key tuple.

6)  The above access conventions also work on the instance, e.g.:
      AWG.table.value('12')    yields  0.0808
      Greek.table.value('α')   yields  'alpha'
      Physics.table.value('c') yields  299792458.0

7)  User defined 'Catalog' sub-class methods and properties are accessible
    only on the instance e.g. AWG.table.ampacity(key) or AWG().ampacity(key)

  When constructing a new 'Catalog' whose items are tuples, the index 0
  of each item should be the primary data value, and the index -1 of each
  item should be a string describing the table entry in order to be
  consistent with the above descibed usage conventions.

  For a simple 'Catalog', each item may be a scalar, where a single string
  value is considered a scalar, omitting any descriptive string.  The
  'Catalog' doc-string gives instructions for 'Catalog' construction.
  Also see the 'DataTable' doc-string.
"""
#from __future__ import unicode_literals # metaclass attribute name must be str not unicode
from math import pi
from geometry.space import V
__all__=['DataTable','Catalog','AWG','Paper','SFM','Drill','Greek', \
         'Physics','WhitePoint','WebNameToColor','WebColorToName', \
         'WebNameToCamel']
__version__=20181216 # correct doc-string errata
__version_history__=(20170331,20171224,20180223,20180329,20181127)
# super calls avoid attribute lookup enforcement conflicts
class DataTable(type):
  """Metaclass for an immutable singleton simplified data dictionary.
     'DataTable' type classes provide only basic operations:
       0) catalog retrieval:     'instance=UserCatalog()'
       1) key access:            'instance[key]' or 'UserCatalog[key]'
       2) iteration over keys:   'for key in instance: ...' or 'in UserCatalog'
       3) key containment:       'if key in instance: ...' or 'in UserCatalog'
       4) meta-property 'table': 'UserCatalog.table is instance' (True)
       5) optional UserCatalog.method and UserCatalog.property definitions
          that use only these limited operations and set no attributes.
       6) meta-method 'value':    'UserCatalog.value(key) or instance.value(key)
       7) meta-property 'keys':   'UserCatalog.keys' or 'instance.keys'
       8) meta-property 'choices':'UserCatalog'.choices' or instance.choices'"""
  def __new__(cls,name,bases,dictionary):
    """Setup access controls on singleton."""
    dictionary["_inst"]=None; dictionary["_dict"]={}
    dictionary["__slots__"]=[] # forbids altering instance attributes
    dictionary["__getattribute__"]=DataTable._getattribute_proxy
    dictionary["__setitem__"]=classmethod(DataTable._setitem_proxy)
    dictionary["__getitem__"]=classmethod(lambda cls,key: \
      super(DataTable,cls).__getattribute__("_dict")[key] )
    dictionary["__iter__"]=classmethod(lambda cls: \
      iter(super(DataTable,cls).__getattribute__("_dict")) )
    dictionary["__len__"]=classmethod(lambda cls: \
      len(super(DataTable,cls).__getattribute__("_dict")) )
    return super(DataTable,cls).__new__(cls,name,bases,dictionary)
  def __init__(cls,name,bases,dictionary):
    """Set the singleton instance then delete 'cls.__init__'."""
    super(DataTable,cls).__setattr__("_inst",super(DataTable, cls).__call__())
    super(DataTable,cls).__delattr__("__init__")
  def __call__(cls):
    """Return the singleton of this 'DataTable' class."""
    return cls.table
  def __delattr__(cls,name):
    """Forbid deletion of class attributes."""
    raise AttributeError,"'{}' forbids deleting class attributes.".format(cls.__name__)
  def __getattribute__(cls,name):
    """Prevent indirect corruption of private data via methods on the
       data by denying access to the private data."""
    if name in ("__init__","_dict","_inst"): return None # this is a compromise:
    # raise AttributeError, ... causes help(cls) to print repr(cls)
    return super(DataTable,cls).__getattribute__(name)
  def __setattr__(cls,name,val):
    """Forbid class attribute alteration."""
    raise AttributeError,"'{}' forbids altering class attributes.".format(cls.__name__)
  # instance __setattr__ and __delattr__ are forbidden by empty __slots__
  @staticmethod
  def _getattribute_proxy(self,name):
    """Method to intercept instance.__getattribute__(...), forbidding
       private data corruption but permitting method and property access."""
    if name in ("__init__","_dict","_inst"): raise AttributeError, \
      "'{}' instance has no readable attribute '{}'.".format(self.__class__.__name__,name)
    if name in ("value","choices","keys","table","__len__"):
      return self.__metaclass__.__getattribute__(self.__class__,name)
    return super(Catalog,self).__getattribute__(name)
  @staticmethod # converted to classmethod when 'dictionary[__setitem__]' is set
  def _setitem_proxy(cls,key,val): # using @classmethod directly is incorrect
    """Method to intercept instance.__setitem__(...), permitting initialization
       but forbidding subsequent alteration."""
    if cls.table is None: super(DataTable,cls).__getattribute__("_dict")[key]=val
    else: raise KeyError, \
      "'{}' instance is a readonly 'DataTable' type.".format(cls.__name__)
  def __getitem__(cls,key):
    """Meta-method returning the item indexed by 'key' in the singleton
       of this 'DataTable' class object."""
    return cls.table[key]
  def __iter__(cls):
    """Meta-method returning iterator over the singleton of this
       'DataTable' class object."""
    return iter(cls.table)
  def __len__(cls):
    """Meta-method returning the length of the singleton of this
       'DataTable' class object."""
    return len(cls.table)
  def value(cls,key,default=None):
    """Meta-method returning the index zero value of the item indexed by 'key'
       in the singleton of this 'DataTable' class or the item if scalar.  A
       single string is considered a 'scalar', i.e. the whole string is returned.
       Returns 'default' on 'KeyError' - does not raise."""
    try:
      t=cls.table[key]
      if isinstance(t,str): return t
      try: return t[0]
      except TypeError: return t
    except KeyError: return default
  @property
  def choices(cls):
    """Meta-property returning a sorted tuple of '(description,key)' tuples.
       If table items are scalars, returns a sorted tuple of keys."""
    t,keys = cls.table,cls.keys
    if isinstance(t[keys[0]],str): return tuple(sorted(keys))
    try: return tuple(sorted((t[k][-1],k) for k in keys))
    except TypeError: return tuple(sorted(keys))
  @property
  def keys(cls):
    """Meta-property returning sorted keys, using 'repr(key)' comparison
       to avoid unicode/str issues."""
    # Note: Direct sorting raises 'UnicodeDecodeError' because of differing types.
    return tuple(sorted(cls.table,key=repr))
  @property
  def table(cls):
    """Meta-property returning the singleton of this 'DataTable' class."""
    return super(DataTable,cls).__getattribute__("_inst")


class Catalog(object):
  """Base class for a readonly simplified singleton dictionary.
     Usage:
       1) create a sub-class 'class UserCatalog(Catalog):'.
       2) write a parameterless __init__ method setting 'self[key]=value'
          for each key,value pair in the immutable.
       3) optionally write properties and methods that return computed
          attributes without setting any new instance or class attributes."""
  __metaclass__=DataTable
  def __init__(self): pass


class AWG(Catalog):
  """'Catalog' of American wire gauge sizes.  Keys are 'int(<gauge_number>)'
     where '-1','-2','-3' are used instead of '00','000','0000'.  Value
     units are in inches."""
  def __init__(self):
    for g in range(-3,41): self[g]=float("{:5.4f}".format( .005*(92**((36-g)/39.0)) ))
  def ampacity(self,key):
    """Calculate the current carrying capacity of wire gauge 'key' in amps."""
    return self.cross_section(key)/(700.0*pi*.0005*.0005) # 700 circular mils/amp
  def cross_section(self,key):
    """Calculate the cross sectional area of wire gauge 'key' in
       square inch units."""
    return pi*self[key]*self[key]/4.0


class Paper(Catalog):
  """'Catalog' of standard paper and display parameters. The key values
     are more or less defacto standard usage.  The value tuple is:
       (width,height), units, descriptive_string
     where 'units' is in ('mm','in','px')."""
  def __init__(self):
    self['ha'] = ( 5.5,8.5),'in',"ANSI_HA_portrait"
    self['a'] = ( 8.5,11.0),'in',"ANSI_A_portrait"
    self['b'] = (11.0,17.0),'in',"ANSI_B_portrait"
    self['c'] = (17.0,22.0),'in',"ANSI_C_portrait"
    self['d'] = (22.0,34.0),'in',"ANSI_D_portrait"
    self['e'] = (34.0,44.0),'in',"ANSI_E_portrait"
    for k in tuple(self): self[k.upper()] = \
      tuple(reversed(self[k][0])),'in',"ANSI_{}_landscape".format(k.upper())
    self['cga']    = ( 320.0, 200.0),'px',"cga_320x200"
    self['qvga']   = ( 320.0, 240.0),'px',"qvga_320x240"
    self['vga']    = ( 640.0, 480.0),'px',"vga_640x480"
    self['ntsc']   = ( 720.0, 480.0),'px',"ntsc_720x480"
    self['hd450']  = ( 854.0, 450.0),'px',"hd450_854x450"
    self['svga']   = ( 800.0, 600.0),'px',"svga_800x600"
    self['pal']    = ( 768.0, 576.0),'px',"pal_768x576"
    self['wsvga']  = (1024.0, 600.0),'px',"wsvga_1024x600"
    self['xga']    = (1024.0, 768.0),'px',"xga_1024x768"
    self['wntsc']  = (1152.0, 768.0),'px',"wntsc_1152x768"
    self['wxga']   = (1280.0, 768.0),'px',"wxga_1280x768"
    self['qcga']   = (1280.0, 800.0),'px',"qcga_1280x800"
    self['sxga']   = (1280.0,1024.0),'px',"sxga_1280x1024"
    self['hd720']  = (1280.0, 720.0),'px',"hd720_1280x720"
    self['hd768']  = (1366.0, 768.0),'px',"hd768_1366x768"
    self['sntsc']  = (1440.0, 960.0),'px',"sntsc_1440x960"
    self['sxga+']  = (1400.0,1050.0),'px',"sxga+_1400x1050"
    self['wsxga']  = (1680.0,1050.0),'px',"wsxga_1680x1050"
    self['uxga']   = (1600.0,1200.0),'px',"uxga_1600x1200"
    self['hd1080'] = (1920.0,1080.0),'px',"hd1080_1920x1080"
    self['wuxga']  = (1920.0,1200.0),'px',"wuxga_1920x1200"
    self['qxga']   = (2048.0,1536.0),'px',"qxga_2048x1536"
    self['wqxga']  = (2560.0,1600.0),'px',"wqxga_2560x1600"
    self['qsxga']  = (2560.0,2048.0),'px',"qsxga_2560x2048"
    self['a0'] = ( 841.0,1189.0),'mm',"METRIC_A0_portrait"
    self['a1'] = ( 594.0, 841.0),'mm',"METRIC_A1_portrait"
    self['a2'] = ( 420.0, 594.0),'mm',"METRIC_A2_portrait"
    self['a3'] = ( 297.0, 420.0),'mm',"METRIC_A3_portrait"
    self['a4'] = ( 210.0, 297.0),'mm',"METRIC_A4_portrait"
    self['a5'] = ( 148.0, 210.0),'mm',"METRIC_A5_portrait"
    self['a6'] = ( 105.0, 148.0),'mm',"METRIC_A6_portrait"
    self['a7'] = (  74.0, 105.0),'mm',"METRIC_A7_portrait"
    self['a8'] = (  52.0,  74.0),'mm',"METRIC_A8_portrait"
    self['a9'] = (  37.0,  52.0),'mm',"METRIC_A9_portrait"
    self['a10']= (  26.0,  37.0),'mm',"METRIC_A10_portrait"
    for k in ('a{}'.format(i) for i in range(11)): self[k.upper()] = \
      tuple(reversed(self[k][0])),'mm',"METRIC_{}_landscape".format(k.upper())
  def aspect_ratio(self,key):
    """Calculate the aspect ratio of 'Paper.table[key]'."""
    x,y = self[key][0]
    return x/y


class SFM(Catalog):
  """'Catalog' translating material names to surface feet per minute (SFM)."""
  def __init__(self):
    self['Al']       = 250.0
    self['Cu_alloy'] = 125.0
    self['Cu']       =  85.0
    self['Ti']       =  40.0
    self['SS']       =  50.0
    self['Fe_cast']  =  70.0
    self['Fe_hard']  =  35.0
    self['Fe_alloy'] =  75.0
    self['Fe_tool']  =  40.0
    self['Wood']     = 330.0
    self['Plastic']  = 450.0


class Drill(Catalog):
  """'Catalog' of standard drill sizes.  The key values are formatted as:
       letter drills:     letter.upper() e.g. 'A'
       fractional drills: 'num/denom'    e.g. '1/2' or '1/1'
       metric drills:     str(float)     e.g. '0.5' or '1.0' or '0.75'
       number drills:     str(integer)   e.g. '80'
     The value is a float(<diameter in inches>) to four places."""
  def __init__(self):
    self['A'] = 0.2340
    self['B'] = 0.2380
    self['C'] = 0.2420
    self['D'] = 0.2460
    self['E'] = 0.2500
    self['F'] = 0.2570
    self['G'] = 0.2610
    self['H'] = 0.2660
    self['I'] = 0.2720
    self['J'] = 0.2770
    self['K'] = 0.2810
    self['L'] = 0.2900
    self['M'] = 0.2950
    self['N'] = 0.3020
    self['O'] = 0.3160
    self['P'] = 0.3230
    self['Q'] = 0.3320
    self['R'] = 0.3390
    self['S'] = 0.3480
    self['T'] = 0.3580
    self['U'] = 0.3680
    self['V'] = 0.3770
    self['W'] = 0.3860
    self['X'] = 0.3970
    self['Y'] = 0.4040
    self['Z'] = 0.4130
    self['1/64']  = 0.0156
    self['1/32']  = 0.0312
    self['3/64']  = 0.0469
    self['1/16']  = 0.0625
    self['5/64']  = 0.0781
    self['3/32']  = 0.0938
    self['7/64']  = 0.1094
    self['1/8']   = 0.1250
    self['9/64']  = 0.1406
    self['5/32']  = 0.1562
    self['11/64'] = 0.1719
    self['3/16']  = 0.1875
    self['13/64'] = 0.2031
    self['7/32']  = 0.2188
    self['15/64'] = 0.2344
    self['1/4']   = 0.2500
    self['17/64'] = 0.2656
    self['9/32']  = 0.2812
    self['19/64'] = 0.2969
    self['5/16']  = 0.3125
    self['21/64'] = 0.3281
    self['11/32'] = 0.3438
    self['23/64'] = 0.3594
    self['3/8']   = 0.3750
    self['25/64'] = 0.3906
    self['13/32'] = 0.4063
    self['27/64'] = 0.4219
    self['7/16']  = 0.4375
    self['29/64'] = 0.4531
    self['15/32'] = 0.4688
    self['31/64'] = 0.4844
    self['1/2']   = 0.5000
    self['33/64'] = 0.5156
    self['17/32'] = 0.5312
    self['35/64'] = 0.5469
    self['9/16']  = 0.5625
    self['37/64'] = 0.5781
    self['19/32'] = 0.5938
    self['39/64'] = 0.6094
    self['5/8']   = 0.6250
    self['41/64'] = 0.6406
    self['21/32'] = 0.6562
    self['43/64'] = 0.6719
    self['11/16'] = 0.6875
    self['45/64'] = 0.7031
    self['23/32'] = 0.7188
    self['47/64'] = 0.7344
    self['3/4']   = 0.7500
    self['49/64'] = 0.7656
    self['25/32'] = 0.7813
    self['51/64'] = 0.7969
    self['13/16'] = 0.8125
    self['53/64'] = 0.8281
    self['27/32'] = 0.8438
    self['55/64'] = 0.8594
    self['7/8']   = 0.8750
    self['57/64'] = 0.8906
    self['29/32'] = 0.9062
    self['59/64'] = 0.9219
    self['15/16'] = 0.9375
    self['61/64'] = 0.9531
    self['31/32'] = 0.9688
    self['63/64'] = 0.9844
    self['1/1']   = 1.0000
    self['0.5']  = 0.0197
    self['0.75'] = 0.0295
    self['1.0']  = 0.0394
    self['1.25'] = 0.0492
    self['1.5']  = 0.0591
    self['1.75'] = 0.0689
    self['2.0']  = 0.0787
    self['2.25'] = 0.0886
    self['2.5']  = 0.0984
    self['2.75'] = 0.1083
    self['3.0']  = 0.1181
    self['3.25'] = 0.1280
    self['3.5']  = 0.1378
    self['3.75'] = 0.1476
    self['4.0']  = 0.1575
    self['4.25'] = 0.1673
    self['4.5']  = 0.1772
    self['4.75'] = 0.1870
    self['5.0']  = 0.1969
    self['5.25'] = 0.2067
    self['5.5']  = 0.2165
    self['5.75'] = 0.2264
    self['6.0']  = 0.2362
    self['6.25'] = 0.2461
    self['6.5']  = 0.2559
    self['6.75'] = 0.2657
    self['7.0']  = 0.2756
    self['7.25'] = 0.2854
    self['7.5']  = 0.2953
    self['7.75'] = 0.3051
    self['8.0']  = 0.3150
    self['8.25'] = 0.3248
    self['8.5']  = 0.3346
    self['8.75'] = 0.3445
    self['9.0']  = 0.3543
    self['9.25'] = 0.3642
    self['9.5']  = 0.3740
    self['9.75'] = 0.3839
    self['10.0'] = 0.3937
    self['10.5'] = 0.4134
    self['11.0'] = 0.4331
    self['11.5'] = 0.4528
    self['12.0'] = 0.4724
    self['12.5'] = 0.4921
    self['13.0'] = 0.5118
    self['13.5'] = 0.5315
    self['14.0'] = 0.5512
    self['14.5'] = 0.5709
    self['15.0'] = 0.5906
    self['15.5'] = 0.6102
    self['16.0'] = 0.6299
    self['16.5'] = 0.6496
    self['17.0'] = 0.6693
    self['17.5'] = 0.6890
    self['18.0'] = 0.7087
    self['18.5'] = 0.7283
    self['19.0'] = 0.7480
    self['19.5'] = 0.7677
    self['20.0'] = 0.7874
    self['20.5'] = 0.8071
    self['21.0'] = 0.8268
    self['21.5'] = 0.8465
    self['22.0'] = 0.8661
    self['22.5'] = 0.8858
    self['23.0'] = 0.9055
    self['23.5'] = 0.9252
    self['24.0'] = 0.9449
    self['24.5'] = 0.9646
    self['25.0'] = 0.9843
    self['80'] = 0.0135
    self['79'] = 0.0145
    self['78'] = 0.0160
    self['77'] = 0.0180
    self['76'] = 0.0200
    self['75'] = 0.0210
    self['74'] = 0.0225
    self['73'] = 0.0240
    self['72'] = 0.0250
    self['71'] = 0.0260
    self['70'] = 0.0280
    self['69'] = 0.0292
    self['68'] = 0.0310
    self['67'] = 0.0320
    self['66'] = 0.0330
    self['65'] = 0.0350
    self['64'] = 0.0360
    self['63'] = 0.0370
    self['62'] = 0.0380
    self['61'] = 0.0390
    self['60'] = 0.0400
    self['59'] = 0.0410
    self['58'] = 0.0420
    self['57'] = 0.0430
    self['56'] = 0.0465
    self['55'] = 0.0520
    self['54'] = 0.0550
    self['53'] = 0.0595
    self['52'] = 0.0635
    self['51'] = 0.0670
    self['50'] = 0.0700
    self['49'] = 0.0730
    self['48'] = 0.0760
    self['47'] = 0.0785
    self['46'] = 0.0810
    self['45'] = 0.0820
    self['44'] = 0.0860
    self['43'] = 0.0890
    self['42'] = 0.0935
    self['41'] = 0.0960
    self['40'] = 0.0980
    self['39'] = 0.0995
    self['38'] = 0.1015
    self['37'] = 0.1040
    self['36'] = 0.1065
    self['35'] = 0.1100
    self['34'] = 0.1110
    self['33'] = 0.1130
    self['32'] = 0.1160
    self['31'] = 0.1200
    self['30'] = 0.1285
    self['29'] = 0.1360
    self['28'] = 0.1405
    self['27'] = 0.1440
    self['26'] = 0.1470
    self['25'] = 0.1495
    self['24'] = 0.1520
    self['23'] = 0.1540
    self['22'] = 0.1570
    self['21'] = 0.1590
    self['20'] = 0.1610
    self['19'] = 0.1660
    self['18'] = 0.1695
    self['17'] = 0.1730
    self['16'] = 0.1770
    self['15'] = 0.1800
    self['14'] = 0.1820
    self['13'] = 0.1850
    self['12'] = 0.1890
    self['11'] = 0.1910
    self['10'] = 0.1935
    self['9']  = 0.1960
    self['8']  = 0.1990
    self['7']  = 0.2010
    self['6']  = 0.2040
    self['5']  = 0.2055
    self['4']  = 0.2090
    self['3']  = 0.2130
    self['2']  = 0.2210
    self['1']  = 0.2280
  def rpm(self,key,sfm):
    """Calculate the drilling rpm to use for 'Drill.table[key]'
       with a material requiring 'sfm' surface feet per minute."""
    return (sfm*12.0)/(pi*self[key])


# These are standard HTML/CSS color names and RGB values
# There are 139 unique names; 141 with two alias names.
# The SVG specification additionally recognizes
#   aqua           : alias for cyan
#   darkgrey       : alias for darkgray
#   darkslategrey  : alias for darkslategray
#   dimgrey        : alias for dimgray
#   fuchsia        : alias for magenta
#   grey           : alias for gray
#   lightgrey      : alias for lightgray
#   lightslategrey : alias for lightslategray
#   slategrey      : alias for slategray
# that are not currently included in 'color_value'.
# In addition, 'RebeccaPurple' is not in the SVG color list.
class WebNameToColor(Catalog):
  """'Catalog' translating camel case HTML/CSS color names to hexadecimal
     SRGB strings e.g. 'DarkGreen'.  Values are expressed as '#rrggbb'.
     'Aqua' and 'Fuchsia' are aliases for 'Cyan' and 'Magenta'"""
  def __init__(self):
    self['Black']                = '#000000'
    self['Navy']                 = '#000080'
    self['DarkBlue']             = '#00008b'
    self['MediumBlue']           = '#0000cd'
    self['Blue']                 = '#0000ff'
    self['DarkGreen']            = '#006400'
    self['Green']                = '#008000'
    self['Teal']                 = '#008080'
    self['DarkCyan']             = '#008b8b'
    self['DeepSkyBlue']          = '#00bfff'
    self['DarkTurquoise']        = '#00cde1'
    self['MediumSpringGreen']    = '#00fa9a'
    self['Lime']                 = '#00ff00'
    self['SpringGreen']          = '#00ff7f'
    self['Cyan']                 = '#00ffff'
    self['MidnightBlue']         = '#191970'
    self['DodgerBlue']           = '#1e90ff'
    self['LightSeaGreen']        = '#20b2aa'
    self['ForestGreen']          = '#228b22'
    self['SeaGreen']             = '#2e8b57'
    self['DarkSlateGray']        = '#2f4f4f'
    self['LimeGreen']            = '#32cd32'
    self['MediumSeaGreen']       = '#3cb371'
    self['Turqoise']             = '#40e0d0'
    self['RoyalBlue']            = '#4169e1'
    self['SteelBlue']            = '#4682b4'
    self['DarkSlateBlue']        = '#483d8b'
    self['MediumTurqoise']       = '#48d1cc'
    self['Indigo']               = '#4b0082'
    self['DarkOliveGreen']       = '#556b2f'
    self['CadetBlue']            = '#5f9ea0'
    self['CornflowerBlue']       = '#6495ed'
    self['RebeccaPurple']        = '#663399'
    self['MediumAquaMarine']     = '#66cdaa'
    self['DimGray']              = '#696969'
    self['SlateBlue']            = '#6a5acd'
    self['OliveDrab']            = '#6b8e23'
    self['SlateGray']            = '#708090'
    self['LightSlateGray']       = '#778899'
    self['MediumSlateBlue']      = '#7b68ee'
    self['LawnGreen']            = '#7cfc00'
    self['Chartreuse']           = '#7fff00'
    self['AquaMarine']           = '#7fffd4'
    self['Maroon']               = '#800000'
    self['Purple']               = '#800080'
    self['Olive']                = '#808000'
    self['Gray']                 = '#808080'
    self['SkyBlue']              = '#87ceeb'
    self['LightSkyBlue']         = '#87cefa'
    self['BlueViolet']           = '#8a2be2'
    self['DarkRed']              = '#8b0000'
    self['DarkMagenta']          = '#8b008b'
    self['SaddleBrown']          = '#8b4513'
    self['DarkSeaGreen']         = '#8fbc8f'
    self['LightGreen']           = '#90ee90'
    self['MediumPurple']         = '#9370db'
    self['DarkViolet']           = '#9400d3'
    self['PaleGreen']            = '#98fb98'
    self['DarkOrchid']           = '#9932cc'
    self['YellowGreen']          = '#9acd32'
    self['Sienna']               = '#a0522d'
    self['Brown']                = '#a52a2a'
    self['DarkGray']             = '#a9a9a9'
    self['LightBlue']            = '#add8e6'
    self['GreenYellow']          = '#adff2f'
    self['PaleTurquoise']        = '#afeeee'
    self['LightSteelBlue']       = '#b0c4de'
    self['PowderBlue']           = '#b0e0e6'
    self['FireBrick']            = '#b22222'
    self['DarkGoldenRod']        = '#b8860b'
    self['MediumOrchid']         = '#ba55d3'
    self['RosyBrown']            = '#bc8f8f'
    self['DarkKhaki']            = '#bdb76b'
    self['Silver']               = '#c0c0c0'
    self['MediumVioletRed']      = '#c71585'
    self['IndianRed']            = '#cd5c5c'
    self['Peru']                 = '#cd853f'
    self['Chocolate']            = '#d2691e'
    self['Tan']                  = '#d2b48c'
    self['LightGray']            = '#d3d3d3'
    self['Thistle']              = '#d8bfd8'
    self['Orchid']               = '#da70d6'
    self['GoldenRod']            = '#daa520'
    self['PaleVioletRed']        = '#db7093'
    self['Crimson']              = '#dc143c'
    self['Plum']                 = '#dda0dd'
    self['Gainsboro']            = '#dcdcdc'
    self['BurlyWood']            = '#deb887'
    self['LightCyan']            = '#e0ffff'
    self['Lavender']             = '#e6e6fa'
    self['DarkSalmon']           = '#e9967a'
    self['Violet']               = '#ee82ee'
    self['PaleGoldenRod']        = '#eee8aa'
    self['LightCoral']           = '#f08080'
    self['Khaki']                = '#f0e68c'
    self['AliceBlue']            = '#f0f8ff'
    self['HoneyDew']             = '#f0fff0'
    self['Azure']                = '#f0ffff'
    self['SandyBrown']           = '#f4a460'
    self['Wheat']                = '#f5deb3'
    self['Beige']                = '#f5f5dc'
    self['WhiteSmoke']           = '#f5f5f5'
    self['MintCream']            = '#f5fffa'
    self['GhostWhite']           = '#f8f8ff'
    self['Salmon']               = '#fa8072'
    self['AntiqueWhite']         = '#faebd7'
    self['Linen']                = '#faf0e6'
    self['LightGoldenRodYellow'] = '#fafad2'
    self['OldLace']              = '#fdf5e6'
    self['Red']                  = '#ff0000'
    self['Magenta']              = '#ff00ff'
    self['DeepPink']             = '#ff1493'
    self['OrangeRed']            = '#ff4500'
    self['Tomato']               = '#ff6347'
    self['HotPink']              = '#ff69b4'
    self['Coral']                = '#ff7f50'
    self['DarkOrange']           = '#ff8c00'
    self['LightSalmon']          = '#ffa07a'
    self['Orange']               = '#ffa500'
    self['LightPink']            = '#ffb6c1'
    self['Pink']                 = '#ffc0cb'
    self['Gold']                 = '#ffd700'
    self['PeachPuff']            = '#ffdab9'
    self['NavajoWhite']          = '#ffdead'
    self['Moccasin']             = '#ffe4b5'
    self['Bisque']               = '#ffe4c4'
    self['MistyRose']            = '#ffe4e1'
    self['BlanchedAlmond']       = '#ffebcd'
    self['PapayaWhip']           = '#ffefd5'
    self['LavenderBlush']        = '#fff0f5'
    self['SeaShell']             = '#fff5ee'
    self['CornSilk']             = '#fff8dc'
    self['LemonChiffon']         = '#fffacd'
    self['FloralWhite']          = '#fffaf0'
    self['Snow']                 = '#fffafa'
    self['Yellow']               = '#ffff00'
    self['LightYellow']          = '#ffffe0'
    self['Ivory']                = '#fffff0'
    self['White']                = '#ffffff'
    self['Aqua']                 = '#00ffff'
    self['Fuchsia']              = '#ff00ff'
class WebColorToName(Catalog):
  """'Catalog' translating web color strings in the form '#rrggbb' to
     camel case web color names."""
  def __init__(self):
    for k in WebNameToColor: self[WebNameToColor[k]]=k
class WebNameToCamel(Catalog):
  """'Catalog' translating lower case web color names to camel case."""
  def __init__(self):
    for k in WebNameToColor: self[k.lower()]=k

class Greek(Catalog):
  """'Catalog' translating Greek character names to the utf-8
     representation and translating each utf-8 code point to the
     corresponding name.  Capital versions are accessed via key
     capitalization e.g. 'Alpha' instead of 'alpha'.  Hex utf-8 byte
     string keys are also acceptable.  So, for example:

        Greek()['alpha']  => unichr(0x3b1) (python type 'unicode')

        Greek()[u'u03b1']      => 'alpha' | decoded unicode:
        Greek()[u'α']          => 'alpha' | key type 'unicode'
        Greek()[unichr(0x3b1)] => 'alpha' |

        Greek()['α']           => 'alpha' | two bytes - encoded:
        Greek()['\\xce\\xb1']    => 'alpha' | key type 'str'

     'encoded' means the byte string stored in memory - may be utf-8,
     utf-16, utf-32, or ascii - requires mapping with 'repr' or 'bytearray'
     to produce unified comparison types for sorting.
     'decoded' means converted to python's internal representation which
     can vary, but is always uniformly handled.

     Notice the variant 'ϕ' encoded as '\\xcf\\x95' is not in the utf-8
     contiguous sequence, but 'φ' encoded as '\\xcf\\x86' is."""
     # NOTA BENE: '\\xce\\xb1' in doc strings prints as '\xce\b1' which is the
     #            correct usage entry text.  Single escapes print as 'α'.
  def __init__(self):
    names=("Alpha", "Beta", "Gamma", "Delta", "Epsilon", "Zeta", "Eta", \
           "Theta", "Iota", "Kappa", "Lamda", "Mu", "Nu", "Xi", "Omicron", \
           "Pi", "Rho", "Sigma", "Tau", "Upsilon", "Phi", "Chi", "Psi", "Omega")
    code_pts=range(0x391,0x3a2)+range(0x3a3,0x3aa) # 0x3a2 and 0x3c2 not used
    for n in range(len(names)):
      uchr,nm = unichr(code_pts[n]),names[n]
      self[nm]=uchr; self[uchr]=nm; self[uchr.encode('utf-8')]=nm
      uchr,nm = unichr(code_pts[n]+0x20),names[n].lower()
      self[nm]=uchr; self[uchr]=nm; self[uchr.encode('utf-8')]=nm
  def py2(self,key):
    """Used to convert the utf-8 representation to the python 2.x string
       used when a source utf-8 symbol name indexes a dictionary.
       For example:  Physics.table.[Greek.table.py2('alpha')]
       properly accesses the tuple defined by Physics.table['α'], which
       in python 2.x has an actual key == '\\ce\\b1', i.e. the string
       representation of the bytearray encoding utf-8 in memory.
       In python 3.x this is not an issue, for all strings are unicode,
       so the actual key implied by source text 'α' is u'α'."""
    return str(bytearray(self[key],'utf8'))


class Physics(Catalog):
  """'Catalog' of physical constants in SI units.  Values are tuples of
     the form:
       (numerical_constant, base_unit_string, description_string)
     The base_unit_string is symbolic as follows:
       'l' :  length          (meters)
       'j' :  energy          (Joules)
       'q' :  electric charge (Coulombs)
       't' :  time            (seconds)
       'k' :  temperature     (Kelvins)
     In contrast to MKSA, mass (kilograms) and current (amperes) are the
     derived units 'j*t*t/(l*l)' and 'q/t' respectively.

     In electrodynamics, the units of frequency and velocity can often
     be more basic than those of length and time.  Moreover, we adopt
     the convention that angular measure has units and use the base unit
     'cycle' == 2*pi radians.  Additionally, in plasma physics, it can
     be useful to count integral objects, to which we assign the base
     unit 'count'.  Finally, for symmetry, it is useful to define the
     unit of magnetic charge.  Therefore, where possible, units are
     abbreviated with:
       'n' :  count     (integral # of items, e.g. atoms)
       'a' :  turns     ( 2*pi radians - traditionally dimensionless)
       'v' :  velocity  ('l/t' or meters/second)
       'f' :  frequency ('a/t' or cycles/second)
       'w' :  magnetic charge (Weber)
     making the kilogram 'j/(v*v)'.  Note also w*q == j*t, for conversions.
     Also notice α == 4.0*(ϕ/q)*sqrt(ε/μ)=4.0*(ϕ/q)/Z0 and units of w/q == Ω."""
  def __init__(self):
    self['q']  = (1.6021766208e-19,'q',"electron charge") # NIST standard
    self['φ']  = (2.067833831170082e-15,'w',"magnetic flux quantum") # calculated: h/(2q)
    self['qm'] = (-1.758820024e11,'q*v*v/j',"electron charge to mass ratio") # NIST standard
    self['c']  = (299792458.,'v',"speed of light") # NIST standard
    self['μ']  = (pi*4.0e-7,'w*w/j',"vacuum permeability") # defined
    self['ε']  = (8.854187817620389e-12,'q*q/j',"vacuum permittivity") # calculated: 1/(c*c*μ)
    self['Z0'] = (376.73031346177066,'w/q',"vacuum radiation impedance") # calculated: c*μ or sqrt(μ/ε)
    self['h']  = (6.626070040e-34,'j*t/a',"Planck constant") # NIST standard
    self['k']  = (1.38064852e-23,'j/(k*n)',"Boltzmann constant") # NIST standard
    self['Av'] = (6.022140857e23,'n',"Avogadro's number") # NIST standard
    self['α']  = (137.03599914178824,'n',"fine structure constant reciprocal") # calculated: 4.0*ϕ*sqrt(μ/ε)/q
    self['σ']  = (5.670366816083269e-08,'j/(t*l*l*k*k*k*k)',"Stefan-Boltzmann constant") # calculated: (2./15.)*(pi^5 * k^4)/(c^2 * h^3)
    self['b']  = (2.897772914526216e-3,'l*k',"Wien's displacement constant") # calculated: h*c/(k*x), x root of (x-5)*exp(x)+5 = 0
    self['Lsol'] =(3.828e26,'j/t',"solar luminosity") # wikipedia            # see https://en.wikipedia.org/wiki/Wien%27s_displacement_law
  def units(self,key):
    """Get the SI base unit string of the constant 'key'."""
    return self[key][1]


class WhitePoint(Catalog):
  """'Catalog' of CIE standard illuminants.  Values are tuples of the form:
       (<CIE-XYZ 'V' object>,<color temperature, Kelvin>,<description>)
     where 'V' is a three dimensional 'geometry.space.V' vector with floating
     point coordinates in the range 0.0 < value < 1.27 , the temperature is a
     float in the range 0.0 < value < 10000.0 and the description is text."""
     # source: https://www.mathworks.com/help/images/ref/whitepoint.html
     # &       http://www.brucelindbloom.com/index.html?Eqn_ChromAdapt.html
  def __init__(self):
    self['A']  =(V(1.09850,1.00000,0.35585),2856.0,"CIE illuminant A: tungsten lamp")
    self['B']  =(V(0.99072,1.00000,0.85223),5000.0,"CIE illuminant B: <no reference specified>") # color temp is guess
    self['C']  =(V(0.98074,1.00000,1.18232),6774.0,"CIE illuminant C: daylight - north sky")
    self['E']  =(V(1.00000,1.00000,1.00000),10000.0,"CIE illuminant E: equal energy radiator") # color temp is guess
    self['F2'] =(V(0.99186,1.00000,0.67393),4100.0,"CIE illuminant F2: cool white fluorescent lamp")
    self['F7'] =(V(0.95041,1.00000,1.08747),6500.0,"CIE illuminant F7: broadband fluorescent lamp")
    self['F11']=(V(1.00962,1.00000,0.64350),4000.0,"CIE illuminant F11: triband fluorescent lamp")
    self['ICC']=(V(31595.0,32768.0, 27030.0)/32768.0,5000,"ICC illuminant: approximates (0.9642,1.000,0.8249) by fixed point") # color temp is guess
    self['D50']=(V(0.96422,1.00000,0.82521),5003.0,"CIE illuminant D65: daylight - sunrise")
    self['D55']=(V(0.95682,1.00000,0.92149),5500.0,"CIE illuminant D55: daylight - morning")
    self['D65']=(V(0.95047,1.00000,1.08883),6504.0,"CIE illuminant D65: daylight - noon")
    self['D75']=(V(0.94972,1.00000,1.22638),7500.0,"CIE illuminant D75: daylight - overcast")


if __name__=="__main__":
  print "Check default instance is singleton:",Paper.table is Paper()
  paper=Paper();
  print "Check 'created' instance is singleton:",paper is Paper.table
  try:
    paper.__init__(); print "instance.__init__() error raising failed"
  except AttributeError as e: print "successfully raised error:\n  ",e
  try:
    del(Paper._inst); print "del(class._inst) error raising failed"
  except AttributeError as e: print "successfully raised error:\n  ",e
  try:
    if Paper._dict is None: print "private class attribute access yielded 'None'"
    else: print "private class attribute access denial failed"
  except AttributeError as e: print "private class attribute access raised error:\n  ",e
  try:
    Paper._inst=None; print "class attribute set error raising failed"
  except AttributeError as e: print "successfully raised error:\n  ",e
  try:
    paper["some_key"]="some_nonsense"; print "key,value setting error raising failed"
  except KeyError as e: print "successfully raised error:\n  ",e
  try:
    del(paper._dict); print "del(instance._dict) error raising failed"
  except AttributeError as e: print "successfully raised error:\n  ",e
  try:
    paper._inst=None; print "instance attribute set error raising failed"
  except AttributeError as e: print "successfully raised error:\n  ",e
  try:
    paper._inst; print "instance attribute access error raising failed"
  except AttributeError as e: print "successfully raised error:\n  ",e
  try:
    paper.some_attr="some_nonsense"; print "instance attribute creation error raising failed"
  except AttributeError as e: print "successfully raised error:\n  ",e
  print "Starting iteration over property test:"
  for t in paper.choices: print t
  print "Starting iteration over self and method test:"
  for k in paper: print "{:7}  {:8.4f}".format(k,paper.aspect_ratio(k))
  print "Default value for missing key test: default=={!r}".format(paper.value('some_key'))
  print "Containment test: 'A' in paper == {}; 'some_key' in paper == {}".format( \
         'A' in paper, 'some_key' in paper)
  print "Testing supercall access and setting of forbidden attribute:"
  super(DataTable,Paper).__getattribute__("_dict")["some_key"]="some_nonsense"
  print "  Succeeded: key = '{}' ; value = '{}'".format("some_key",paper["some_key"])
  print SFM.table.keys
  print Drill.table.rpm('1/4',SFM.table['Ti'])
  for k in AWG:
    print k,
    print AWG.table[k],
    print AWG.table.ampacity(k)
  print AWG.table.choices
  print Physics.choices
  pc=Physics
  print "Test meta-method key access:",pc['c']
  from math import sqrt
  print "speed of light {} == {} (1/sqrt(mu*q0))".format( \
            pc.value('c'),1./sqrt(pc.value('μ')*pc.value('ε')))
  print repr(1.0/pc.value('c')/pc.value('c')/pc.value('μ'))
  for k in Greek().keys: print repr(k),k,':',Greek.value(k),Greek.table[k]
  for k in ('alpha','Alpha','α',Greek.table['Alpha']): print Greek.table[k]
  print WhitePoint().choices
  print Physics.table[Greek.table.py2('sigma')]
  print len(WebColorToName), "Fuchsia" in WebColorToName
  help(DataTable)
