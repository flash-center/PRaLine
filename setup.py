from setuptools import setup

setup(name='praline',
      version='0.0.6',
      packages=['praline'],
      entry_points={
          'console_scripts': ['lin-reconstruct = praline.reconstruct:prad_wrap',
                              'lin-analyze = praline.analysis:prad_wrap'],
                },
      )
