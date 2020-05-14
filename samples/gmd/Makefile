ozone_absorption_*.h5: ozone_absorption.jl
	echo julia --startup-file=no $<

wavelength_2*.h5: wavelength_2.jl
	echo julia --startup-file=no $<

%.h5: %.jl
	echo julia --startup-file=no $<

all: asim_27386.h5 first_hitting.h5 first_hitting_lowerbnd.h5 \
	wavelength_2_337nm_10km_10um.h5 ozone_absorption_MidLat_Summer_Hale.h5

