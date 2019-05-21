# coding: utf-8
import os

from setuptools import find_packages, setup


def find_requires():
    dir_path = os.path.dirname(os.path.realpath(__file__))
    requirements = []
    with open(os.path.join(dir_path, "requirements.txt"), "r") as reqs:
        requirements = reqs.readlines()
    return requirements


setup(
    name="simple_eos",
    version="0.1.0",
    description="Simple python package to calculate EOS based on ASE snippet",
    packages=find_packages(),
    python_requires=">=3.6",
    install_requires=find_requires(),
    include_package_data=True,
    zip_safe=True,
    entry_points={
        "console_scripts": [
            "get_eos = simple_eos.cli:main",
        ]
    },
)
