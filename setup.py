from setuptools import setup, find_packages

setup(
    name="pybridizer",
    version="0.2.0",
    packages=find_packages(),
    install_requires=[
        "biopython>=1.79",
        "pandas>=1.3.0",
        "biothings_client>=0.2.6",
        "matplotlib>=3.4.0",
        "seaborn>=0.11.0",
        "numpy>=1.20.0"
    ],
    author="Chintan Trivedi",
    author_email="c.trivedi@ucl.ac.uk",
    description="A tool for designing HCR (Hybridization Chain Reaction) v3 probes",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/ctucl/pybridizer",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    python_requires=">=3.6",
)