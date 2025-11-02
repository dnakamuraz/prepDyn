from setuptools import setup, find_packages

setup(
    name="prepdyn",
    version="0.1.0",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    python_requires=">=3.8",
    install_requires=[
        "biopython",
        # add other dependencies
    ],
    entry_points={
        "console_scripts": [
            "GB2MSA=src.GB2MSA:main",
            "addSeq=src.addSeq:main",
            "prepDyn=src.prepDyn:main",
        ],
    },
)

