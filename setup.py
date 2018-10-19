from setuptools import setup, find_packages

setup(
      name='rnaloc',
      version='0.1.0',
      description='Python code to analyze RNA localization data.',
      url='https://github.com/muellerflorian/rna-loc',
      author='Florian MUELLER',
      author_email='muellerf.research@gmail.com',
      license='MIT',
      packages=find_packages(),
      include_package_data=True,
      install_requires=[
          'numpy',
          'scikit-image',
          'nested_lookup',
          'read-roi',
          'scipy',
          'matplotlib',
          'Pillow',
      ],
      zip_safe=False)
