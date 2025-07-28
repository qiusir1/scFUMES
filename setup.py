from pathlib import Path
from setuptools import setup, find_packages

this_dir = Path(__file__).parent
README = (this_dir / "README.md").read_text(encoding="utf-8")

setup(
    name="scFUMES",
    version="1.0.0",
    description="Single cell FUnctional MEtabolite-Sensor (scFUMES)",
    long_description=README,
    long_description_content_type="text/markdown",
    author="Yunguang Qiu",
    license="MIT",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    include_package_data=True,
    python_requires=">=3.9",
    install_requires=[
        "scanpy>=1.10",
        "pandas>=2.0",
        "numpy>=1.26",
        "scipy>=1.12",
        "statsmodels>=0.14",
    ],
    entry_points={
        "console_scripts": [
            "scfumes=scfumes.cli:main",
        ],
    },
)
