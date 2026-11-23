"""
iVPM (Inverse Visco-Plastic Model)
Closed-form visco-plastic continuum for the co-seismic response of moving landslides.

Reference:
Sportoletti, M. (2026). A closed-form visco-plastic continuum for the co-seismic
response of moving landslides. Zenodo. https://doi.org
"""

import numpy as np

# Physical and time constants
G = 9.81          # Gravitational acceleration (m/s^2)
YR = 31536000.0   # Seconds in a standard 365-day year


def keff(n: float, Vc_mm_yr: float, S: float, alpha_deg: float, rho: float = 2000.0) -> float:
    """
    Compute the operational viscous consistency Keff (in Pa * s^n) calibrated
    from the observed pre-seismic creep rate. Valid only for n > 0.
    """
    if n <= 0:
        raise ValueError("keff is defined only for n > 0.")
    Vc = Vc_mm_yr / 1000.0 / YR
    return (
        rho * G * S * np.sin(np.radians(alpha_deg))
        * (n * S / ((n + 1.0) * Vc)) ** n
    )


def d_ivpm(n: float, Vc_mm_yr: float, alpha_deg: float, Ia: float, Sa: float, T: float) -> float:
    """
    Compute the closed-form co-seismic permanent displacement DiVPM (in m).
    Evaluates cleanly down to the rate-independent limit n = 0.
    """
    Vc = Vc_mm_yr / 1000.0 / YR
    f = (n + 1.0) / (2.0 * n + 1.0)
    return (
        T ** n * (Sa * Ia) ** 2 * Vc ** n
        / (2.0 * f * G * np.sin(np.radians(alpha_deg)))
    ) ** (1.0 / (n + 1.0))


def d_jibson_m(Ia: float, ky: float) -> float:
    """
    Empirical Jibson (2007) rigid-block regression for comparison (in m).
    """
    return 0.01 * 10 ** (2.401 * np.log10(Ia) - 3.481 * np.log10(ky) - 3.230)


def profile(n: float, z: np.ndarray | float, S: float) -> np.ndarray | float:
    """
    Normalized steady-state velocity profile v(z) / Vc across thickness S.
    """
    if n <= 0:
        raise ValueError("profile is defined only for n > 0.")
    return 1.0 - (1.0 - np.clip(z, 0.0, S) / S) ** ((n + 1.0) / n)


if __name__ == "__main__":
    # Benchmark check: La Sorbella strong-motion scenario (Table 2, Irpinia 1980)
    disp_m = d_ivpm(n=0.12, Vc_mm_yr=15.0, alpha_deg=8.0, Ia=1.4, Sa=1.5, T=30.0)
    print(f"La Sorbella (Irpinia 1980, Sa=1.5) -> DiVPM = {disp_m:.3f} m")
