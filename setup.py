from setuptools import find_packages, setup

dependencies=''
with open("requirements.txt","r") as f:
	dependencies = f.read().splitlines()


setup(
	name="obscurae",
	version='0.0.1',
	description='Stellar line profiles distorted by a transiting planet or/and a spot crossing the disk',
	url='https://github.com/emilknudstrup/obscurae',
	author='Emil Knudstrup',
	author_email='emil.knudstrup@chalmers.se',
	packages=find_packages(where="src"),
	package_dir={"": "src"},
	include_package_data=True,
    classifiers = ["Programming Language :: Python :: 3"],
	install_requires = dependencies
)
