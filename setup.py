import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="PyRational",
    version="0.0.1",
    author="Alessio Benavoli",
    author_email="benavoli@gmail.com",
    description="A package for rational modelling and dual probabilisitc programming",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/PyRational/PyRational",
    packages=setuptools.find_packages(),
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: The 3-Clause BSD License",
        "Operating System :: OS Independent",
    ),
)
