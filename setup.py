from setuptools import setup, find_packages

setup(
    name='CNAsim',
    version='1.3.5',
    author='Samson Weiner',
    author_email='samson.weiner@uconn.edu',
    description='CNAsim is a software package for simulation of single-cell CNA data from tumors.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/samsonweiner/CNAsim',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent'
    ],
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'cnasim=CNAsim.main:main',
        ],
    },
    python_requires='>=3.7'
)