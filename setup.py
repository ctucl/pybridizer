from setuptools import setup, find_packages

setup(
    name="pybridizer",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "biopython",
        "pandas",
        "biothings_client",
        "matplotlib",
        "seaborn",
        "numpy"
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
    ],
    python_requires=">=3.6",
)