from setuptools import setup, find_packages

setup(
    name='flowsig',
    version='0.1.2',
    packages=find_packages(exclude=['tests*']),
    description='Uses graphical modeling to infer communication-driven intercellular flows from single-cell RNA-sequencing and spatial transcriptomics data.',
    url='https://github.com/axelalmet/flowsig',
    author='Axel A. Almet',
    author_email='axelalmet@gmail.com',
    include_package_data=True,
    package_data={'': ['data/*.csv', 'data/*.txt']}
)
