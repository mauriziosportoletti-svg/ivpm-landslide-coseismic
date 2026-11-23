# iVPM: Closed-form visco-plastic continuum for co-seismic landslide displacement

Reference Python implementation of the **iVPM (Inverse Visco-Plastic Model)** described in:

> Sportoletti, M. (2026). *A closed-form visco-plastic continuum for the co-seismic response of moving landslides*. Zenodo. https://doi.org

## Overview
`ivpm.py` provides analytical functions to estimate permanent co-seismic displacements (\(D_{\mathrm{iVPM}}\)) for pre-existing moving landslides calibrated directly on pre-seismic creep rates (\(V_c\), e.g., from InSAR monitoring) and slope inclination (\(\alpha\)), without requiring bulk density (\(\rho\)) or sliding thickness (\(S\)) in the operative displacement expression.

## Requirements
- Python 3.9+
- `numpy`

## Quick Usage
```python
from ivpm import d_ivpm

# Example: La Sorbella landslide scenario (n=0.12, Vc=15 mm/yr, alpha=8 deg)
disp = d_ivpm(n=0.12, Vc_mm_yr=15.0, alpha_deg=8.0, Ia=1.4, Sa=1.5, T=30.0)
print(f"Co-seismic displacement: {disp:.3f} m")
```
