import setuptools

setuptools.setup(
    name="scMDCF",
    version="0.0.1",
    author="Yue Cheng",
    author_email="chengyue22@mails.jlu.edu.cn",
    description="scMDCF Enables Cross-modality Cell Heterogeneity Elucidation and Interaction between Multiomic Profiles at Single-cell Resolution",
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url="https://github.com/DARKpmm/scMDCF",
    packages=['scMDCF'],
    classifiers=[
        "Programming Language :: Python :: 3.9",
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering :: Artificial Intelligence',
    ],
)