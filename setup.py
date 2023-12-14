import sys
from distutils.core import setup
from setuptools.command.test import test as TestCommand

#class PyTest(TestCommand):
#    def finalize_options(self):
#        TestCommand.finalize_options(self)
#        self.test_args = []
#        self.test_suite = True
#
#    def run_tests(self):
#        import pytest
#        errno = pytest.main(self.test_args)
#        sys.exit(errno)
#
#long_description = open('README.md').read()

setup(
    name='euler-plate-solver',
    version='0.0.1',
    author='Fred Dux',
    author_email='duxfrederic@gmail.com',
    description='A plate solver',
#    long_description=long_description,
#    long_description_content_type='text/markdown',
#    packages=['starred', 'starred.deconvolution', 'starred.plots', 'starred.psf', 'starred.utils'],
#    requires=['astropy', 'dill', 'jax', 'jaxlib', 'jaxopt', 'matplotlib', 'numpy', 'scipy', 'optax', 'tqdm', 'emcee',
#              'pyregion'],
#    cmdclass={'test': PyTest}
)
