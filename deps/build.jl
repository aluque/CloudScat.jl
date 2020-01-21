using Pkg

@info "Installing MieScatter"

Pkg.activate()
Pkg.add(PackageSpec(url="https://github.com/dronir/MieScatter.jl.git"))
