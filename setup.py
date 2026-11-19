from setuptools import setup, find_packages

setup(
    name='Markov_DNA',
    version='2.0.0',
    packages=["Markov_DNA"],
    package_dir={'':'src'},
    install_requires=['numpy'],
    description='A Markov Model DNA sequence generator to generate pseudo-replicate sequences based on an input sequence.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='Erik Blázquez',
    author_email='erikblazfer@outlook.es',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.9',
)