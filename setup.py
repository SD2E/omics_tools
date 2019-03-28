import setuptools


with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="Omics-dashboard",
    version="1.0",
    author="AC",
    author_email="acristof@mit.edu",
    description="DGE and Annotation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="http://web.mit.edu/foundry/",
    packages=setuptools.find_packages(),
    license="MIT",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)


