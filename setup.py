from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()

setup(
    name='npmine',
    version='0.0.1',
    description='Retrieve chemical information from scientific literature',
    classifiers=[
    'Development Status :: 1 - Planning',
    'Environment :: Web Environment',
    'Intended Audience :: Science/Research',
    'License :: Free for non-commercial use',
    'Natural Language :: English',
    'Programming Language :: Python :: 3.7',
    'Topic :: Software Development :: Libraries :: Python Modules'],
    license='BSD-3-Clause',
    packages=['npmine'],
    author='Ana Carolina L. Coelho',
    author_email='ana.lunardello.coelho@usp.br'
    install_requires=[
          'markdown',
          'beautifulsoup4',
          'requests',
          'json',
          'subprocess',
          'itertools',
    ], 
    dependency_links=['https://www.rdkit.org/docs/Install.html']
)
