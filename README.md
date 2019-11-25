# CloudScat.jl

*Cloud-scattering code*

**CloudScat** is a Montecarlo code that predicts the signal observed by a satellite due to an optical source located inside or above a thundercloud.  It considers both Rayleigh (molecular) scattering and Mie scattering with cloud droplets.


## Install

The code is distributed as a [julia](https://julialang.org) package.  To install it, after you have downloaded and installed julia, start julia and go into the package manager by pressing `]`.  Then type

```julia
# Use Mie scattering code by Olli Wilkman
add https://github.com/dronir/MieScatter.jl.git
add http://gitlab.com/aluque/cloudscat
```
You can exit the package manager by pressing backspace at the beginning of the prompt.

## Use

Look at the file `samples/sample.jl` for an example computation.  You can run it
with
```bash
> julia --color=yes samples/sample.jl
```
or, ig you prefer to do everything from the julia prompt,
```julia
julia> include("samples/sample.jl")
```

If the simulation finishes correctly it produces a file called `sample.h5` that contains the output data from the simulation. To view this output you can use one of the scripts contained in the `util` folder of the repo. For example to plot the simulated photometer data of observer one uses

```bash
> python plot_photo.py input_file.h5 --observer=1
```

To see the image recorded by the same observer use
```bash
> python plot_image input_file.h5 --observer=1
```


