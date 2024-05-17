from setuptools import setup

setup(name='ldsc',
      version='2.0b',
      description='LD Score Regression (LDSC)',
      url='http://github.com/abrantesas/ldsc_py3',
      author='Anthony Abrantes, Brendan Bulik-Sullivan and Hilary Finucane',
      author_email='antshaabr@gmail.com',
      license='GPLv3',
      packages=['ldscore'],
      scripts=['ldsc.py', 'munge_sumstats.py'],
      install_requires = [
            'bitarray>=2.5,<2.6',
            'scipy>=1.11,<1.12',
            'numpy>=1.24,<1.25',
            'pandas>=2.0,<2.1',
            'pyranges>=0.0.129'
      ]
)
