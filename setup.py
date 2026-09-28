from pathlib import Path

from setuptools import setup


README = Path(__file__).with_name("README.md").read_text(encoding="utf-8")


setup(
    name="prepdyn",
    version="0.7.0",
    description="Preprocessing toolkit for dynamic homology.",
    long_description=README,
    long_description_content_type="text/markdown",
    py_modules=[
        "GB2MSA",
        "UP2AP",
        "addSeq",
        "prepDyn",
        "prepDyn_auxiliary",
    ],
    package_dir={"": "src"},
    python_requires=">=3.8",
    include_package_data=False,
    install_requires=[
        "biopython",
        "matplotlib",
        "numpy",
        "termcolor",
    ],
    entry_points={
        "console_scripts": [
            "GB2MSA=GB2MSA:main",
            "UP2AP=UP2AP:main",
            "addSeq=addSeq:main",
            "prepDyn=prepDyn:main",
        ],
    },
)
