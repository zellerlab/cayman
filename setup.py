# coding: utf-8
from setuptools import setup, find_packages
from codecs import open
from os import path
import sys

from cayman import __version__
here = path.abspath(path.dirname("__file__"))

with open(path.join(here, "DESCRIPTION.md"), encoding="utf-8") as description:
	description = long_description = description.read()

	name="cayman"
	version = __version__

	if sys.version_info.major != 3:
		raise EnvironmentError(f"""{name} is a python module that requires python3, and is not compatible with python2.""")

	setup(
		name=name,
		version=version,
		description=description,
		long_description=long_description,
		url="https://github.com/zellerlab/cayman",
		author="Christian Schudoma",
		author_email="cschu1981@gmail.com",
		license="MIT",
		classifiers=[
			"Development Status :: 4 - Beta",
			"Topic :: Scientific/Engineering :: Bio-Informatics",
			"License :: OSI Approved :: MIT License",
			"Operating System :: POSIX :: Linux",
			"Programming Language :: Python :: 3.7",
			"Programming Language :: Python :: 3.8",
			"Programming Language :: Python :: 3.9",
			"Programming Language :: Python :: 3.10",
		],
		zip_safe=False,
		keywords="",
		packages=find_packages(exclude=["test"]),
		install_requires=[line.strip() for line in open("requirements.txt", "rt")],
		entry_points={
			"console_scripts": [
				"cayman=cayman.__main__:main",
			],
		},
		package_data={},
		include_package_data=True,
		data_files=[],
	)
