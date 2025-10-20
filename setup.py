#!/usr/bin/env python
"""
Setup script for Quantum Biology Framework
"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="fHEOM",
    version="1.0.0",
    author="Roberto Ignacio Henriquez-Perozo",
    author_email="roberto@henriquezperozo.com",
    description="Low-rank factorization for efficient HEOM simulations of quantum coherence in biological systems",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/rihp/fHEOM",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Physics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.8",
    install_requires=[
        "qutip>=5.2.0",
        "numpy>=1.24.0",
        "scipy>=1.10.0",
        "matplotlib>=3.7.0",
        "seaborn>=0.12.0",
        "pandas>=2.0.0",
        "tqdm>=4.65.0",
    ],
    extras_require={
        "gpu": ["cupy-cuda12x>=13.0.0"],
        "dev": ["pytest>=7.3.0", "black>=23.0.0", "flake8>=6.0.0"],
    },
)

