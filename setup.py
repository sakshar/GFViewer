from setuptools import setup, find_packages

setup(
    name="gfviewer",
    version="1.0.3",
    author="Sakshar Chakravarty",
    author_email="schak026@ucr.edu",
    description="A tool for visualizing the localization of multi-gene families across the genome of a given organism",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/sakshar/GFViewer",
    packages=find_packages(),
    install_requires=[
        "biopython",
        "pandas",
        "matplotlib",
        "reportlab",
        "pypdf2",
        "openpyxl"
    ],
    python_requires=">=3.8, <3.13",  # Restrict Python version
    entry_points={
        "console_scripts": [
            "gfviewer = gfviewer.main:main"
        ]
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
