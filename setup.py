#!/usr/bin/env python3

from setuptools import setup

setup(
    name="divscan",
    version="1.0",
    packages=[
        "divscan",
        "divscan.test"
    ],
    scripts=["bin/divscan.py"],
)