import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="nPhase",
    version="1.0.9",
    author="Omar Abou Saada",
    author_email="oabousaada@unistra.fr",
    description="nPhase is a command line ploidy agnostic phasing pipeline and algorithm which phases samples of any ploidy with sequence alignment of long and short read data to a reference sequence.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://https://github.com/nPhasePipeline/nPhase",
    entry_points = {
        'console_scripts': [
            'nphase = bin.nPhasePipeline:main'
        ]
    },
    packages=setuptools.find_packages(),
    install_requires=["plotnine","sortedcontainers"],
    classifiers=[
        "Programming Language :: Python :: 3.8",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Development Status :: 5 - Production/Stable",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    python_requires='>=3.8',
)

