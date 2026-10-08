from setuptools import setup,Extension
from Cython.Build import cythonize
from pathlib import Path
import sys,ctypes.util
# Linux is the recorded platform. The conventional Homebrew locations also
# permit a source build on macOS; that platform is not claimed as tested here.
include=[];library=[]
for prefix in ['/opt/homebrew','/usr/local']:
 if (Path(prefix)/'include/gmpxx.h').exists():
  include.append(prefix+'/include');library.append(prefix+'/lib')
if sys.platform.startswith('linux'):
 soname=ctypes.util.find_library('mpfr') or 'libmpfr.so.6'
 extra=['-l:'+soname,'-lgmpxx','-lgmp']
else:extra=['-lmpfr','-lgmpxx','-lgmp']
setup(name='kg-independent-mpfr-endpoints',ext_modules=cythonize([
 Extension('flint',['flint.pyx'],language='c++',include_dirs=include,
  library_dirs=library,extra_compile_args=['-O3','-std=c++17'],extra_link_args=extra)
],compiler_directives={'language_level':3}))
