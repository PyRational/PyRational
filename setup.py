import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="PyRational",
    version="0.0.4",
    author="Alessio Benavoli",
    author_email="benavoli@gmail.com",
    description="A package for rational modelling and dual probabilisitc programming",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/PyRational/PyRational",
    packages = ["PyRational","PyRational.models","PyRational.special"],
    package_dir={'PyRational': 'PyRational'},  
    include_package_data = True,
    install_requires = ['numpy>=1.7', 'scipy>=0.16', 'SymPy>=1.2'],
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
    ),
)
