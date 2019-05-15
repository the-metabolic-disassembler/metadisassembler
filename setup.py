from setuptools import setup

with open("requirements.txt") as f:
    requirements = f.read().splitlines()

with open("README.rst") as f:
    long_desc = f.read()

def main():
    setup(
        name="metadisassembler",
        version="0.0.7",
        url="https://github.com/the-metabolic-disassembler/metadisassembler",
        author="Kohei Amano",
        author_email="amanok2167@gmail.com",
        maintainer="Masaaki Kotera",
        maintainer_email="maskot@chemsys.t.u-tokyo.ac.jp",
        description="The Metabolic Disassembler",
        long_description=long_desc,
        license="MIT",
        keywords=[
            "biosynthesis",
            "metabolic pathway",
            "natural product",
            "cheminformatics"
        ],
        classifiers=[
            "Development Status :: 5 - Production/Stable",
            "License :: OSI Approved :: MIT License",
            "Operating System :: Microsoft :: Windows :: Windows 10",
            "Operating System :: MacOS :: MacOS X",
            "Operating System :: POSIX :: Linux",
            "Programming Language :: Python :: 3.6",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            "Topic :: Scientific/Engineering :: Chemistry",
            "Topic :: Software Development :: Libraries :: Python Modules",
        ],
        packages=["metadisassembler"],
        entry_points={
            'console_scripts':[
                'metadisassembler = metadisassembler.__init__:main',
                ],
        },
        install_requires=requirements,
        zip_safe=False,
        include_package_data=True,
    )


if __name__ == "__main__":
    main()
