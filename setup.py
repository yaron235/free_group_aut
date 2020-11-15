import setuptools

with open("README.md", "r") as f:
    long_description = f.read()

setuptools.setup(
    name="free_group_aut",
    version="0.0.1",
    author="Yaron Brodsky",
    author_email="yaron235@gmail.com",
    description="A package implementing algorithms on free groups",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yaron235/free_group_aut",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',
)
