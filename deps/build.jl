using Pkg

@info "Installing MieScatter"

Pkg.activate("CloudScat")
Pkg.add(PackageSpec(url="https://github.com/dronir/MieScatter.jl.git"))


