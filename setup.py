#
# processGenbank setuptools script
#
#
from setuptools import setup, find_packages


def get_version():
    """
    Get version number from the genbank_utils module.
    The easiest way would be to just ``import genbank_utils ``
    """
    import os
    import sys

    sys.path.append(os.path.abspath('genbank_utils'))
    from version_info import VERSION as version # noqa
    sys.path.pop()

    return version


def get_readme():
    """
    Load README.md text for use as description.
    """
    with open('README.md', encoding='utf-8') as f:
        return f.read()


setup(
    # Module name (lowercase)
    name='genbank_utils',

    # Version
    version=get_version(),

    description='This module contains juptyter notebook and python code to process extract rodent virus sequences from genbank files for probe panel sequence selection',  # noqa

    long_description=get_readme(),

    license='',

    # author='',

    # author_email='',

    maintainer='',

    maintainer_email='',

    url='https://github.com/jnarag/processGenbank.git',

    # Packages to include
    packages=find_packages(include=('genbank_utils', 'genbank_utils.*')),
    include_package_data=True,

    # List of dependencies
    install_requires=[
        # Dependencies go here!
        'biopython'
    ],
    python_requires='>=3.10',
)
