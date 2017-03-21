# MJDSigGen.jl

[![Build Status](https://travis-ci.org/mppmu/MJDSigGen.jl.svg?branch=master)](https://travis-ci.org/mppmu/MJDSigGen.jl)
[![Coverage Status](https://coveralls.io/repos/github/mppmu/MJDSigGen.jl/badge.svg?branch=master)](https://coveralls.io/github/mppmu/MJDSigGen.jl?branch=master)
[![codecov.io](http://codecov.io/github/mppmu/MJDSigGen.jl/coverage.svg?branch=master)](http://codecov.io/github/mppmu/MJDSigGen.jl?branch=master)

MJDSigGen.jl provides a [Julia](http://julialang.org/) wrapper around
[David Radford's](http://radware.phy.ornl.gov/) field-calculation and
signal-generation package for non-segmented high-purity Germanium detectors,
[mjd_siggen](http://radware.phy.ornl.gov/MJ/mjd_siggen/).

If you refer to this software in a scientific publication, please primarily
cite the mjd_siggen package and it's author, David Radford.

The mjd_siggen C-code is included in this Julia package, with kind permission
of the author.
