from setuptools import setup, find_packages
import os

def package_files(directory):
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            paths.append(os.path.join(path, filename))
    return paths

with open('README.md') as f:
    long_description = f.read()

setup(
    name="bridgernadesigner",
    version='0.0.1_dev',
    description='Software package to design bridge RNAs as described by Durrant & Perry et al.',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/hsulab-arc/BridgeRNADesigner',
    author="Matt Durrant",
    author_email="matthew@arcinstitute.org",
    license="MIT",
    packages=find_packages(),
    package_data={
        "bridgernadesigner": package_files("data/"), #+ package_files("snakemake/")
    },
    include_package_data=True,
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