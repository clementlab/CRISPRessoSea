from setuptools import setup, find_packages

setup(
    name="CRISPRessoSea",
    version="0.1.3",
    description="CRISPRessoSea: pooled CRISPR analysis pipeline",
    author="Kendell Clement",
    author_email="k.clement@utah.edu",
    url="https://github.com/clementlab/CRISPRessoSea",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "pandas",
        "matplotlib",
        "seaborn",
        "scipy",
        "statsmodels",
        "jinja2",
        "crispresso2",
    ],
    entry_points={
        "console_scripts": [
            "CRISPRessoSea=CRISPRessoSea.CRISPRessoSea:main"
        ]
    },
    include_package_data=True,
    python_requires=">=3.7",
)
