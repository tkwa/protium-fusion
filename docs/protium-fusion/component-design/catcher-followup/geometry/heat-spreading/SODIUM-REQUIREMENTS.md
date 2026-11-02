# ✨ Sodium spreader flow requirements

A short sodium phase-change circuit has enough thermodynamic capacity to be a credible heat-spreader direction. The calculation below establishes **flow, pressure and inventory scales**, while the k≈934 transient response still requires a complete evaporator/contact/return/condenser design.

[sodium_requirements.py](sodium_requirements.py) uses the actual **0.88115 MW/m² C13 material-plus-delayed peak** and 116.52 W/m² mean. At 1000 K, published sodium correlations give saturation pressure **19.22 kPa**, latent heat **4.015 MJ/kg**, vapor density **0.05824 kg/m³** and surface tension **0.13545 N/m**. [Argonne property documentation](https://sas-doc.nse.anl.gov/5.4/Part04/Ch12/12.13.html); [Fink–Leibowitz assessment](https://doi.org/10.2172/94649). The companion [primary evidence review](material-evidence/heat-pipes.md) records property quality and distinguishes measured loaded performance from optimistic effective-conductivity normalizations.

The full peak requires **0.21946 kg/(m²·s)** of phase flow. For 2 mm-radius vapor channels:

| Vapor-area fraction | Vapor speed | Reynolds number | Laminar-reference pressure drop | Saturation-temperature equivalent |
|---:|---:|---:|---:|---:|
| 2% | 188.4 m/s | 2348 | 1268 Pa | 5.42 K |
| 5% | 75.36 m/s | 939 | 524 Pa | 2.24 K |

The 2% case approaches transition, so its laminar pressure loss is a reference rather than a qualified result. The 5% case is a more comfortable preliminary flow screen. Pressure loss alone does not include entrainment, evaporation or condensation resistance.

The inventory example also assigns 0.5% liquid-return area, 0.2% wick-solid volume and 20 μm channel walls with selected solid density 8000 kg/m³. Restoring the original **350 kg/m² carbon** around these volumes increases total thickness from 175 mm to **179.93/185.77 mm** in the 2%/5% cases. Wall, wick, liquid and vapor together add **4.16/5.18 kg/m²**, before headers, evaporator and condenser faces, contacts, seals and fittings. This is a geometric inventory scenario, not a selected alloy or a complete bill of materials.

The peak liquid return speed is **0.0563 m/s**. A uniformly fine wick must meet capillary pressure and permeability simultaneously; declaring a small meniscus radius alone does not close it. The saved curve gives required Darcy permeability against that constraint. A separate example uses a 2 mm-radius liquid artery plus a 100 μm-thick fine evaporator over 50% of the face, with specified permeability 10⁻¹² m². Their calculated liquid drops are only about **4 Pa and 10.5 Pa**. Ideal complete wetting with 10 μm menisci would supply approximately **27.1 kPa** capillary pressure. These deliberately separated functions avoid assigning the whole 186 mm return length the fine evaporator's permeability. Actual wetting, interface geometry, dryout and entrainment remain qualification inputs.

The evaporator must collect heat across substantial front area. Putting it only on the 5% vapor-channel ends would impose **17.6 MW/m² locally**, compared with **0.881 MW/m²** over a full-area evaporator. A planar wick/manifold, dense distributed contacts or another explicit collection geometry is needed; vapor cross-section alone does not supply the required thermal conductance. The prompt flash first deposits energy locally and can feed the spreader later.

The circuit stays sealed from catcher gas and ash. Fluid inventory and additional solids must fit a named service allocation or enter the plant mass denominator. Blind/staggered channels or new coupled photon transport are needed if voids alter shielding; conserving average carbon column is not enough to inherit the previous capture calculation. The material packet provides useful high-temperature operating precedent, without certifying this pulsed receiver architecture. [Numerical requirements](sodium-requirements.json).
