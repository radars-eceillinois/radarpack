"""A setuptools based setup module.
See:
https://packaging.python.org/en/latest/distributing.html
"""

import setuptools
from glob import glob

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="radarpack",
    version="0.0.1",
    author="Pablo Reyes",
    author_email="pablo.reyes@sri.com",
    description="Radar Tools",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/illinoisrsssradar/radarpack",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    #install_requires=['numpy','scipy',
    #    'pyigrf @ git+https://github.com/radars-eceillinois/pyigrf.git'],
    # Optional
    #package_data={'pyigrf': ['pyigrf/*.txt'],
    #                         },
    #include_package_data=True,
)
