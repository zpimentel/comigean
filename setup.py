from setuptools import setup, find_packages
import pkg_resources


with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="comigean",
    version="0.0.38",
    author="Zachary T. Pimentel",
    author_email="zpimentel@uri.edu",
    description="Enables comparitive phylogenomic analyses of microbial genomes",
    long_description=long_description,
    long_description_content_type="text/markdown",
    #url="",
    packages=['comigean',],
    #packages = find_packages(),
    entry_points = {'console_scripts': ['comigean = comigean.cli:commands_entry']},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires = ['argh']
)
