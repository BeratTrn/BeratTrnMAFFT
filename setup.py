from setuptools import setup, find_packages

setup(
    name="BeratTrnMAFFT", 
    version="0.1.0",
    author="Berat",
    description="Biyoinformatik dersi MAFFT algoritması dönem projesi",
    packages=find_packages(),
    install_requires=[
        "numpy", 
    ],
    python_requires='>=3.6',
)