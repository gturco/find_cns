from distutils.core import setup, Extension
from Cython.Distutils import build_ext

setup(
        name="blastmisc"
        ,ext_modules = [
            Extension("blast_misc", ["blast_misc.pyx"]) # ["blast_misc.pyx"])
         ]
        ,cmdclass = {'build_ext':build_ext}
)


