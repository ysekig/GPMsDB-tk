# setup.py
import os

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext


def version():
    setup_dir = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(setup_dir, 'GPMsDB_tk', 'VERSION'), 'r') as f:
        return f.readline().strip()


ext_modules = [
    Extension('GPMsDB_tk.calc', sources=['GPMsDB_tk/calc.pyx'])
]

setup(
    cmdclass={'build_ext': build_ext},
    ext_modules=ext_modules
)

setup(
    name='GPMsDB_tk',
    python_requires='>=3.6',
    version=version(),
    author='Yuji Sekiguchi',
    author_email='y.sekiguchi@aist.go.jp',
    packages=['GPMsDB_tk'],
    scripts=['bin/GPMsDB_tk'],
    package_data={'GPMsDB_tk': ['VERSION']},
    url='',
    mdclass={'build_ext': build_ext},
    ext_modules=ext_modules,
    description='Toolkit for bacterial and archaeal identification based on MALDI-TOF MS peak lists.'
)
