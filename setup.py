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
    package_data={'': ['data/*.csv', 'data/*.txt']},
    python_requires='>=3.8',
    install_requires=[
        'numpy',
        'scipy',
        'anndata',
        'pandas',
        'scanpy',
        'spatial-factorization @ git+https://github.com/willtownes/spatial-factorization-py.git@c7c7fbb22a7ed9abccca94b759e4601b78d874b',
        'matplotlib',
        'seaborn',
        'networkx',
        'causaldag',
        'graphical_models',
        'joblib',
        'squidpy'
    ]
)
