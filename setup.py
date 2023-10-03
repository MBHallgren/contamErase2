from setuptools import setup, find_packages

setup(
    name='NanoDecon',
    version='1.0.0',
    packages=find_packages(),
    data_files=[],
    include_package_data=True,
    url='https://https://github.com/MBHallgren/NanoDecon',
    license='',
    install_requires=(),
    author='Malte B. Hallgren',
    scripts=['bin/nanodecon'],
    author_email='malhal@food.dtu.dk',
    description='NanoDecon - Decontamination of Nanopore sequencing data'
)