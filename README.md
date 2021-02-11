# MJDSigGen.jl

[![Documentation for development version](https://img.shields.io/badge/docs-dev-blue.svg)](https://legend-exp.github.io/MJDSigGen.jl/dev)
[![License](http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat)](LICENSE.md)
[![Build Status](https://github.com/legend-exp/MJDSigGen.jl/workflows/CI/badge.svg?branch=master)](https://github.com/legend-exp/MJDSigGen.jl/actions?query=workflow%3ACI)
[![Codecov](https://codecov.io/gh/legend-exp/MJDSigGen.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/legend-exp/MJDSigGen.jl)


## Documentation

* [Documentation for development version](https://legend-exp.github.io/MJDSigGen.jl/dev)

MJDSigGen.jl provides a [Julia](http://julialang.org/) wrapper around
[David Radford's](http://radware.phy.ornl.gov/) field-calculation and
signal-generation package for non-segmented high-purity Germanium detectors,
[mjd_siggen](http://radware.phy.ornl.gov/MJ/mjd_siggen/).

If you refer to this software in a scientific publication, please primarily
cite the mjd_siggen package and it's author, David Radford.

The mjd_siggen C-code is included in this Julia package, with kind permission
of the author.
