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
[mjd_siggen](https://github.com/radforddc/icpc_siggen).

If you refer to this software in a scientific publication, please primarily
cite the mjd_siggen package and it's author, David Radford.

The mjd_siggen C-code is included in this Julia package, with kind permission
of the author.

## Basic Usage

To get started you'll need a config file to specify detector geometry and some simulation parameters. See the [example.config](examples/config/example.config) to get started. Next, you'll need a lookup table of drift velocities as in [drift_vel_torr.tab](examples/config/drift_vel_tcorr.tab). See that the *drift_name* field in the config file matches the relative path from the config to the velocity table file.

Now we can start using MJDSigGen. If you haven't installed it already, install it as follows:

```julia
julia> using Pkg; Pkg.add(url="https://github.com/legend-exp/MJDSigGen.jl.git")
```

We'll also need the weighting field, potential and electric field in the detector. To simulate these we use the `fieldgen` function as follows (This may take a while). The electric field data and the weighting potential data will be written to the files specified in the config under *field_name* and *wp_name* respectively.

```julia
julia> using MJDSigGen;

julia> fieldgen("examples/config/example.config");
```

We can now start simulating signals. Initialise a `SigGenSetup` object like so:

```julia
julia> setup = SigGenSetup("examples/config/example.config");
```

To simulate a signal produced by a drifting charge cluster starting at `(10, 10, 10)`, do:

```julia
julia> get_signal!(setup, (10, 10, 10))
q: -1.0 t: 122 n: 1 ((10.00 10.00 10.00)=>(12.67 12.67 -0.20))
q: 1.0 t: 375 n: 1 ((10.00 10.00 10.00)=>(0.80 0.80 -0.18))
800-element Array{Float32,1}:
 0.0014776316
 0.0047704997
 0.008109845
 0.011552348
 0.015177331
 ⋮
 1.0000001
 1.0000001
 1.0000001
 1.0000001
```

We can also look at the path, velocity or cluster size of the last simulation by examining the relevant properties of our `setup` object. For example, to see the velocity of the holes cluster, we can use:

```julia
julia> setup.instant_vel_h
8000×3 LinearAlgebra.Transpose{Float32,Base.ReshapedArray{Float32,2,Base.ReinterpretArray{Float32,1,MJDSigGen.CartPoint{Float32},Array{MJDSigGen.CartPoint{Float32},1}},Tuple{}}}:
 -0.0332866  -0.0332866  0.0399359
 -0.0333873  -0.0333873  0.0394825
 -0.0334852  -0.0334852  0.0390309
 -0.0335816  -0.0335816  0.0385836
 -0.0336817  -0.0336817  0.0381274
  ⋮                      
  0.0         0.0        0.0
  0.0         0.0        0.0
  0.0         0.0        0.0
  0.0         0.0        0.0
```

The three columns refer to x, y and z components respectively. Velocity is given in *step_time_calc* ns steps (see config), while the signal vector is given in *step_time_out* ns steps.
