import numpy as np
import scipy.constants as co
from scipy.interpolate import interp1d
import yaml

import PyMieScatt

def main():
    # Refraction index relative to the medium

    # Jursa for 377 nm; note that the convention here is to use positive
    # imaginary part for absorption
    # m = 1.345 + 8.45e-9j

    lmbd = 777 * co.nano
    
    wl, n, k = np.loadtxt("Hale.dat", unpack=True)
    wl *= co.micro

    m = interp1d(wl, n)(lmbd) + 1j * interp1d(wl, k)(lmbd)

    print(f"lmbd = {lmbd}")
    print(f"m = {m}")
    
    r = 20 * co.micro
    d = 2 * r
    
    Qext, Qsca, Qabs, g, Qpr, Qback, Qratio = PyMieScatt.MieQ(m, lmbd, d)

    print(f"Qext = {Qext}")
    print(f"Qsca = {Qsca}")
    print(f"Qabs = {Qabs}")
    print(f"g = {g}")

    # Single-scattering albedo:
    omega0 = Qsca / Qext

    print(f"ω₀ = {omega0}")
    print(f"1 - ω₀ = {1 - omega0}")

    print("=" * 24)
    print(yaml.dump({
        "parameters": {
            "Qext": float(Qext),
            "g": float(g),
            "ω₀": float(omega0)}},
                    indent=4, 
                    default_flow_style=False,
                    allow_unicode=True))
    
if __name__ == '__main__':
    main()
