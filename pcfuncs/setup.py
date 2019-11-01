from setuptools import setup
import os.path as path

# Get the path to our current directory
path_here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
readme_path = path.join(path_here, "README.md")
with open(readme_path) as readme_file:
    readme_text = readme_file.read()

setup(
    name='pcfuncs',
    version="0.0",
    author="Dr. Samuel Payne",
    author_email="sam_payne@byu.edu",
    description="Helper functions for protein complex analysis",
    long_description=readme_text,
    long_description_content_type="text/markdown",
    url="https://github.com/PayneLab/ProteinComplexes_2.0/tree/master/pcfuncs",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: Apache Software License',
    ],
    python_requires='>=3.6',
	install_requires=[
        'cptac>=0.6.5',
		'numpy>=1.16.3',
		'pandas==0.25.*',
		'requests>=2.21.0',
		'beautifulsoup4>=4.7.1',
		'scipy>=1.2.1',
		'urllib3>=1.24.2',
    ],
	)
