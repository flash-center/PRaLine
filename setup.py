from setuptools import setup

setup(name='lin-prad',
      version='0.1.0',
      packages=['lin_prad'],
      entry_points={
          'console_scripts': ['lin-reconstruct = lin_prad.reconstruct:prad_wrap',
                              'lin-analyze = lin_prad.analysis:prad_wrap'],
                },
      )
