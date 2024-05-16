from setuptools import setup, find_packages
import os

with open('README.md') as f:
    long_description = f.read()

setup(
    name="bridgernadesigner",
    version='0.0.1',
    description='Software package to design bridge RNAs as described by Durrant & Perry et al. 2024',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/hsulab-arc/BridgeRNADesigner',
    author="Matt Durrant",
    author_email="matthew@arcinstitute.org",
    license="MIT",
    packages=find_packages(),
    install_requires=[
        'biopython',
        'click==7.0'
    ],
    zip_safe=False,
    entry_points = {
        'console_scripts': [
            'brna-design = bridgernadesigner.main:cli'
        ]
}
)