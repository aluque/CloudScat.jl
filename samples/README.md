# Example files

This folder contains example files for CloudScat.  Each example is heavily documented and the combination of all examples should make for the lack of a proper documentation.

## File description

* `sample.jl`: The simplest example.  The cloud is modelled as a wide cylinder with a homogeneous composition.
* `scan_depths.jl`: A set of simulations for an impulsive and localized sources at varying cloud depth and droplet radius.
* `cloud_geometry.jl`: Shows how to define complex cloud geometries and how to produce Mathematica code to render the cloud geometry in 3D.
* `variable_droplet_radius.jl`: Here we use an inhomogeneous cloud composition where the droplet radius changes as a funcion of altitude.
* `rings.jl`: Shows how to define you own geometrical shapes and in particular how to create new shapes based on an implicit equation.  Here we model the cloud top as containing concentric ondulations.

## Use

After you have installed the CloudScat package you can find
the location of the `samples` folder by
```julia
> using CloudScat
> CloudScat.SAMPLES_PATH
```

It is recommended that you copy the example files to a folder of your choice and use them as templates for your own simulations.  You can run e.g. the `sample.jl` example with
```bash
> julia --color=yes sample.jl
```
or, if you prefer to do everything from the julia prompt,
```julia
julia> include("sample.jl"); run()
```

The files are designed to produce an output file based on the input filename. For example `sample.jl` produces a `sample.h5` that you can view with the provided tools.  In this manner, you can just copy any of the input files, apply any desired change and start a simulation.
