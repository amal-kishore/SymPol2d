"""
Setup script for SYMPOL2D
"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="sympol2d",
    version="0.1.0",
    author="SYMPOL2D Development Team",
    description="SYMmetry-based prediction of POLarity in 2D bilayers",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/sympol2d",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
    python_requires=">=3.7",
    install_requires=[
        "numpy>=1.19.0",
        "ase>=3.20.0",
    ],
    entry_points={
        "console_scripts": [
            "sympol2d=sympol2d.cli:main",
        ],
    },
)