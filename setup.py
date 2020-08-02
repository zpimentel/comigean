from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="MicroContext",
    version="0.0.38",
    author="Zachary T. Pimentel",
    author_email="zpimentel@uri.edu",
    description="Enables comparitive phylogenomic analyses of microbial genomes",
    long_description=long_description,
    long_description_content_type="text/markdown",
    #url="",
    packages = find_packages(),
    entry_points = {'console_scripts': ['comoidean = comigean.cli:main']},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires = ['argh']
)
