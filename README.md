# rust_quma

An implementation of the `quma` module from [**Pyllelic**](https://github.com/Paradoxdruid/pyllelic) in [Rust](https://www.rust-lang.org/).

## Goal

This implementation intends to explore the conversion of computationally intensive python code in pyllelic into Rust for performance optimization.

(Note: The python implementation uses [bio-python](https://biopython.org/) for pairwise alignment, which is [already implemented in C](https://github.com/biopython/biopython/blob/master/Bio/Align/_aligners.c))

Once cloned and dependencies resolved, you can build this module locally with: `maturin develop`

To compare the speed of the implementation, run `python test.py`

## About pyllelic

[**Pyllelic**](https://github.com/Paradoxdruid/pyllelic): a tool for detection of allelic-specific methylation variation in bisulfite DNA sequencing files.

Pyllelic documention is available at **<https://paradoxdruid.github.io/pyllelic/>**

## Authors

This software is developed as academic software by [Dr. Andrew J. Bonham](https://github.com/Paradoxdruid) at the [Metropolitan State University of Denver](https://www.msudenver.edu). It is licensed under the GPL v3.0.

This software incorporates implementation from [QUMA](http://quma.cdb.riken.jp), licensed under the GPL v3.0.
