import os, sys
from setuptools import setup, find_packages


def read_requirements():
    """Parse requirements from requirements.txt."""
    reqs_path = os.path.join('.', 'requirements.txt')
    with open(reqs_path, 'r') as f:
        requirements = [line.rstrip() for line in f]
    return requirements

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="omics-dashboard",
    version="1.0",
    author="AC",
    author_email="acristof@mit.edu",
    description="DGE, annotation, and visualization",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="http://web.mit.edu/foundry/",
    packages=find_packages(),
    license="MIT",
    install_requires=read_requirements(),
    classifiers=[
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering :: Bio-Informatics'],
    include_package_data=True
)


