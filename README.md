# scMDCF

[![scMDCF badge](https://img.shields.io/badge/scMDCF-python-blue)](https://github.com/DARKpmm/scMDCF)
[![PyPI badge](https://img.shields.io/pypi/v/scMDCF.svg)](https://pypi.org/project/scMDCF/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

`scMDCF` is a python package containing tools for clustering single cell multi-omics data based on cross-modality contrastive learning to learn the common latent representation and assign clustering.

- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Usage](#usage)
- [Data Availability](#data-availability)
- [License](#license)

# Overview
Single-cell multi-omics (scMulti-omics) technologies have revolutionized our understanding of cellular functions and interactions by enabling the simultaneous measurement of diverse cellular modalities. However, the inherent complexity, high-dimensionality, and heterogeneity of these datasets pose substantial challenges for integration and analysis across different modalities. To address these challenges, we develop a single-cell multi-omics deep learning model (scMDCF) based on contrastive learning, tailored for the efficient characterization and integration of scMulti-omics data. scMDCF features a cross-modality contrastive learning module that harmonizes data representations across different omics types, ensuring consistency while accommodating conditional entropy to preserve data heterogeneity. Furthermore, a cross-modality feature fusion module is designed to extract common low-dimensional latent representations of scMulti-omics data, effectively balancing the characteristics of these diverse omics data. Extensive empirical studies demonstrate that scMDCF outperforms existing state-of-the-art scMulti-omics models across various types of scMulti-omics data. In particular, scMDCF exhibits progressive capability in extracting cell-type specific peak-gene associations and cis-regulatory elements from SNARE-seq data, as well as in elucidating immune regulation from CITE-seq data. Furthermore, we demonstrate that in the post-BNT162b2 mRNA SARS‐CoV‐2 vaccination dataset, scMDCF successfully annotates specific vaccine-induced B cell subpopulations through integrative and multimodal analysis, uncovering dynamic interactions and regulatory mechanisms within the immune system after vaccination.
![The framework plot of scMDCF](https://github.com/DARKpmm/scMDCF/raw/main/scMDCF.png)

# System Requirements
## Hardware requirements
`scMDCF` package requires only a standard computer with enough RAM to support the in-memory operations.

## Software requirements
### OS requirements
This package is supported for *Linux*. The package has been tested on the following systems:
* Linux: Ubuntu 18.04

### Python Dependencies
`scMDCF` mainly depends on the Python scientific stack.
    numpy
    pytorch
    scanpy
    pandas
    scikit-learn
For specific setting, please see <a href="https://github.com/DARKpmm/scMDCF/blob/main/requirements.txt">requirements</a>.

# Installation Guide
## Install from PyPi
    conda create -n scMDCF_env python=3.9.16
    conda activate scMDCF_env
    pip install scMDCF==1.1.2

# Usage
`scMDCF` is a deep embedding learning method for single-cell multi-omics data clustering, which can be used to:
* CITE-seq dataset clustering. The example can be seen in the <a href="https://github.com/DARKpmm/scMDCF/tree/main/tutorial/main_CITE.py">main_CITE.py</a>
* SNARE-seq (paired RNA-seq and ATAC-seq) dataset clustering. The example can be seen in the <a href="https://github.com/DARKpmm/scMDCF/tree/main/tutorial/main_SNARE.py">main_SNARE.py</a>

# Data Availability
The datasets we used can be download in <a href="https://github.com/DARKpmm/scMDCF/tree/main/dataset">dataset</a>

# License
This project is covered under the **MIT License**.
