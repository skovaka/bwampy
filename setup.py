from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext
from glob import glob
import os
import subprocess
import sys

try:
    from pybind11.setup_helpers import Pybind11Extension, ParallelCompile, naive_recompile
    ParallelCompile("NPY_NUM_BUILD_JOBS", needs_recompile=naive_recompile).install()
except:
    from setuptools import Extension as Pybind11Extension

about = {}
with open("_bwampy/__about__.py") as fp:
    exec(fp.read(), about)

ROOT_DIR = os.getcwd()

#SUBMOD_DIR = os.path.join(ROOT_DIR, "submods")
SUBMOD_DIR = ROOT_DIR
SUBMODS = [
    "bwa", 
]

class get_pybind_include(object):
    """Helper class to determine the pybind11 include path
    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked. """

    def __str__(self):
        import pybind11
        from pybind11.setup_helpers import Pybind11Extension
        return pybind11.get_include()

class pre_build(build_ext):
    def run(self):
        submods_loaded = True
        for submod in SUBMODS:
            if not os.path.exists(os.path.join(SUBMOD_DIR, submod)):
                submods_loaded = False
                break

        if not submods_loaded:
            sys.stderr.write("Downloading submodules\n")
            subprocess.check_call([
                "git", "submodule", "update", "--init"
            ])
        else:
            sys.stderr.write("All submodules present\n")

        if os.path.exists("./bwa/libbwa.a"):
            sys.stderr.write("Found libbwa.a\n")
        else:
            sys.stderr.write("building libbwa\n")

            subprocess.check_call([
                "make", 
                 "-C", "./bwa", 
                 "-f", "../src/Makefile_bwa"
            ])

            os.chdir(ROOT_DIR)

        build_ext.run(self)

bwampy = Pybind11Extension(
    "bwampy",

    sources = [ #glob("src/**/*.cpp", recursive=True),
       "src/bwampy.cpp",
    ],

    include_dirs = [
        "./",
        get_pybind_include()
    ],

    library_dirs = [
        "./bwa"
    ],
    
    libraries = ["bwa", "z", "dl", "m"],

    extra_compile_args = ["-std=c++11", "-O3", "-g"],

    define_macros = [("PYBIND", None)]#, ("PYDEBUG", None)]
)

requires=[
    'pybind11>=2.6.0', 
    #'read-until==3.0.0',
    #'pandas>=1.1.5',
    #'plotly>=5.0.0',
    #'dash>=2.0.0',
    #'scipy>=1.5.4',
    #'toml>=0.10.2',
    #'ont_fast5_api',
    #'pysam'
],

setup(
    name = about["__title__"],
    version = about["__version__"],
    description = about["__summary__"],
    author = about["__author__"],
    author_email = about["__email__"],
    url = about["__uri__"],

    classifiers=[
      "Programming Language :: Python :: 3"
    ],

    python_requires=">=3.8",

    setup_requires=['pybind11>=2.6.0'],
    install_requires=requires,

    packages=find_packages(),
    include_package_data=True,
    ext_modules = [bwampy],
    cmdclass={'build_ext': pre_build},
)
