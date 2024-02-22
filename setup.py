from setuptools import find_packages, setup

setup(
    name="td_tools",
    version="2.1.3",
    pythonrequires=">=3.9",
    packages=find_packages(),
    package_data={
        "td_tools": ["*.txt"],
    },
    description="Reaction Transform Descriptor Tools",
    author="DocMinus",
    author_email="alexander.minidis@gmail.com",
    url="https://github.com/DocMinus/RxnTransformDescriptors",
)
