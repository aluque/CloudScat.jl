# CloudScat.jl

*Cloud-scattering code*

**CloudScat** is a Montecarlo code that predicts the signal observed by a satellite due to an optical source located inside or above a thundercloud.  It considers both Rayleigh (molecular) scattering and Mie scattering with cloud droplets.


## Install

The code is distributed as a [julia](https://julialang.org) package.  To install it, after you have downloaded and installed julia, start julia and go into the package manager by pressing `]`.  Then type

```julia
add http://gitlab.com/aluque/cloudscat
```
You can exit the package manager by pressing backspace at the beginning of the prompt.

## Use

### Provided examples

The `samples` folder contains example computations. It is recommended that you 
start by copying these files into another folder and use them as templates for 
your simulations.  After you have installed the CloudScat package you can find
the location of the `samples` folder by
```julia
> using CloudScat
> CloudScat.SAMPLES_PATH
```

The file `sample.jl` contains a simple example.  You can run it
with
```bash
> julia --color=yes sample.jl
```
or, if you prefer to do everything from the julia prompt,
```julia
julia> include("sample.jl"); run()
```

If the simulation finishes correctly it produces a file called `sample.h5` that contains the output data from the simulation. To view this output you can use one of the scripts contained in the `util` folder of the repo. For example to plot the simulated photometer data of observer one uses

```bash
> python plot_photo.py input_file.h5 --observer=1
```

To see the image recorded by the same observer use
```bash
> python plot_image input_file.h5 --observer=1
```

### Complex geometries

CloudScat can manage very complex cloud geometries, defined as combinations and
linear transformations of a set of elementary shapes. Full details are given
in the example file `cloud_geometry.jl`. If you want to define new geometrical
shapes beyond those provided by default, look at `src/geometry.jl`. A geometry 
type must define methods to compute (quickly if possible) intersections between 
the defined shape and a straight line and to test whether a point is inside the 
figure.  For an example on how to do this look at `src/rings.jl`, where the cloud
top has the shape of concentric rings.


### Parallel computations

If your machine contains several cores you can benefit from using several threads for your simulation. This can be accomplished by setting the environment variable
`JULIA_NUM_THREADS` to the number of threads.  For example, for 4 threads invoke
julia as e.g.
```bash
> JULIA_NUM_THREADS=4 julia --color=yes sample.jl
```
