# Comparing directed networks via denoising graphlet distributions

This repository contains code to perform network comparison based on graphlet distributions. This code implements the methods described in the arXiv preprint: https://arxiv.org/abs/2207.09827

It provides a python implementation of the Earth Mover's Distance (EMD), coded by the authors of [Identifying networks with common organizational principles](https://doi.org/10.1093/comnet/cny003). This repository further extends the original code by adapting it to support directed networks and introduces a denoising mechanism using linear projections, aimed at creating reducing noise in graphlet frequencies prior to the comparing their distributions using EMD.

---

# Requirements

Python package requirements are listed in the `requirements.txt` file, generated using anaconda. This code has been tested in Python 3.5 and 3.8, in a Unix and MacOS operating system.

This repository uses [ORCA](https://doi.org/10.1093/bioinformatics/btt717) to calculate the frequencies of undirected graphlets. The interface to call ORCA is implemented in R, tested with R versions 3.4 and 4.2.

To calculate the frequencies of directed graphlets, this repository uses [G-Tries](https://ieeexplore.ieee.org/abstract/document/7501518), with a modified version of the code available at [this webpage](https://www.dcc.fc.up.pt/~daparicio/software.html). G-Tries uses [nauty](https://doi.org/10.1016/j.jsc.2013.09.003) to compute graph isomorphism tests. Nauty and G-Tries require a C and C++ compiler to be used.

---

# Usage

Before using the methods implemented by the code in this repository for directed networks, G-Tries and nauty need to be compiled. To do this follow the steps below:

```
cd gtscanner
make
```

The basic command to use NetEmd is the following:

```python NetEmd.py [input directory] [graphlet size] [graph direction] [NetEmd version] [number of threads] [optional argument]```

- *Input directory*: the directory where the networks to be compared are stored.
- *Graphlet size*: the maximum graphlet size to use for orbit frequencies. If the size given is *k*, then orbit are calculated from graphlets of size *k*, *k-1*, ..., 2 (degree distribution). The allowed values are: **2** (directed networks only), **3**, **4** and **5**.
- *Graph direction*: whether the graph is directed or undirected. The allowed values are: **dir** or **undir**.
- *NetEmd version*: which version of NetEmd to use. It can be one of the following values:
  - **allorbs**: the original NetEmd formulation. Uses all orbits from the specified graphlet size.
  - **weighted**: when comparing two networks, uses all orbits that have occurrences in at least one of the networks. If the frequency of an orbit is 0 in both networks, that orbit is not included in the NetEmd average.
  - **pca**: uses Principal Component Analysis as the linear projection to denoise the graphlet frequencies.
  - **ica**: uses Independent Component Analysis as the linear projection to denoise the graphlet frequencies.
- *Number of threads*: how many threads to use in parallel for the computation. In undirected networks, it controls how many instances of ORCA are launched (how many networks are counted in parallel). In directed networks, it controls how many threads are used by G-Tries to calculate the orbit frequencies of a single network. In directed and undirected networks, it controls how many parallel instances are launched to calculate the EMD between distributions of each orbit (one instance per orbit), using Python's `multiprocessing` package.
- *Optional arguments*: additional arguments for the number of components used in PCA and ICA. If the NetEmd version is **pca**, then this argument should be a number greater than 0 and lesser than 1 that represents the minimum percentage of explained variance included in the components. If the NetEmd version is **ica**, then this argument is the number of components used by ICA, it should be a number greater than 0 and lesser than the total number of orbits. If the NetEmd version is **allorbs** or **weighted**, this argument is ignored.
