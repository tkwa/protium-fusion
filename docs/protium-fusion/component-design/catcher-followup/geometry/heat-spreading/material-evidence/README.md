# ✨ High-temperature graphite evidence for radial heat spreading

**A bulk radial conductivity of a few hundred W/mK has direct experimental support at 1000–1200 K.** The companion heat-spreading calculation's constant effective **325.14 W/mK** requirement for a 1250 K front-temperature screen is a plausible assembly target using oriented pyrolytic graphite. Its **933.64 W/mK** requirement for a 1200 K screen exceeds the measured basal conductivity recovered here by about 2.4 times at 1200 K. The stronger target calls for enhanced transport, such as heat pipes, or a different thermal layout.

This packet supplies material evidence and assembly bounds. The companion finite-slab calculation supplies the temperature requirements, using the same 350 kg/m² receiver, 1000 K baseline, 150 µm surface laminate and actual dispatch. No probability interval is assigned to future assembly conductivity.

## Primary measured conductivity

Roth, Watson, Moss and Drotning measured pyrolytic graphite along its basal planes with laser-flash diffusivity, measured density and a common literature heat-capacity curve. Their compression-annealed specimen was treated at 3000°C for one hour and had density **2240 kg/m³**. The numerical appendix, printed p.84, supports the following interpolation; Figure12 and the appendix page were visually checked. [SAND88-2057, full primary report](https://www.osti.gov/servlets/purl/6188316).

| Temperature | Basal material k, W/mK | Ideal parallel-packing k at effective density 2000 kg/m³, W/mK |
|---|---:|---:|
| 1000 K | 423 | 378 |
| 1100 K | 400 | 357 |
| 1200 K | 386 | 345 |
| 1250 K | 368 | 328 |

The appendix points nearest this window are (750.9°C,404.206 W/mK), (901.0°C,395.811 W/mK) and (1052.8°C,339.491 W/mK). `data.json` preserves the full twelve-point specimen table. `calculate.py` uses linear interpolation without extrapolation. Calibration agreement within 5% and heat-capacity agreement within 3% are reported method checks; they do not establish a total-conductivity confidence interval. The large low-temperature conductivity declines substantially with temperature. [Roth et al., methods and appendix](https://www.osti.gov/servlets/purl/6188316).

The report's less-oriented grades give lower values. Separately, Moore, Graves and McElroy directly measured nearly isotropic AXF-5Q and H337 graphite at **74 and 85 W/mK at 1000 K**, respectively, with an author-stated conductivity uncertainty of 5%. The generic word “graphite” therefore cannot justify a 325 W/mK design value. [Full primary report, TableII](https://www.osti.gov/servlets/purl/4387457-kSLmOu/), [published paper](https://doi.org/10.13182/NT74-A16277).

## Orientation, density and interfaces

For a uniaxial material with basal conductivity k_a, c-axis conductivity k_c, c-axis unit vector c and radial unit vector n,

`k_radial = k_a + (k_c − k_a)(c·n)²`.

Orienting the bulk basal planes to contain the radial direction gives the favorable k_a. The separate surface laminate retains basal planes parallel to the wall for lateral gas routing. Its radial resistance remains in the thermal ledger. Radial lamellae or tiles below that laminate can provide the favorable bulk orientation; their junctions and structural connections must preserve heat flow.

At fixed 350 kg/m², dense 2240 kg/m³ graphite is **156.25 mm** thick. The archived 2000 kg/m³ model is **175 mm** thick. Keeping that thickness by using a solid fraction 2000/2240 gives the last table column under ideal parallel conduction and negligible void conduction. That column is an upper bound for this assumed packing. A constant assembled k=325 lies only about 1% below its value at 1250 K; realistic joints, misorientation or damage can consume this margin. Much of the bulk remains cooler than the front, so this hottest-point comparison is conservative, but only a temperature-dependent assembly calculation can exploit that fact.

At peak flux **0.88115 MW/m²**, an additional 10 K contact-drop allocation requires `R_contact ≤ 1.135e−5 m²K/W`. This is an added-interface requirement. The existing `150 µm / (4 W/mK) + 4e−6 m²K/W` surface-and-post resistance is already included in the companion sweep and must not be counted twice.

The present review does not establish the same specimen's c-axis conductivity at 1000–1200 K. Taylor's primary abstract covers 100–900 K; Null et al. cover room temperature and 1500–2700 K; Vardi and Hoch cover 1260–2200 K. Full Taylor and Null PDFs returned HTTP403. Their temperature coverage and access outcomes are retained in `data.json`; the archived surface-laminate k=4 assumption remains separately qualified. [Taylor1966](https://doi.org/10.1080/14786436608211993), [Null et al.1973](https://doi.org/10.1016/0008-6223(73)90058-4), [Vardi and Hoch1963](https://doi.org/10.1111/j.1151-2916.1963.tb11695.x).

## Irradiation and restoration

Moore et al. irradiated graphite at 823–1023 K to 2.6–5.0×10²⁵ fast neutrons/m², above 50 keV. Conductivity near 900 K fell by roughly a factor of 2.5 in those specimens. This demonstrates that operation near 1000 K does not by itself preserve pristine conductivity. The grade, neutron dose and defect population differ from the proposed basal material. [Moore et al., full report and Fig5](https://www.osti.gov/servlets/purl/4387457-kSLmOu/).

Swank, Rohrbaugh, Cottle and Windes tested NBG-18, IG-110 and PCEA after irradiation around 600±50°C, including 2–4.5 dpa specimens. Thermal-diffusivity recovery began near **1000°C (1273 K)**; only lower-dose specimens fully recovered by their highest **2380°C** anneal. A 1000 K bulk operating temperature is below that reported recovery onset. These experiments support scheduled restoration or replacement when dosimetry requires it; they do not determine its interval for this receiver. [Full primary INL/CON-18-45742 report](https://www.osti.gov/servlets/purl/1484477).

The [species-resolved stopping packet](../../../stopping/README.md) places primary implantation and its strongest displacement dose near the surface. Its thin-layer damage proxy should not be imposed uniformly on the entire bulk slab. Bulk neutron/secondary-radiation exposure, contact degradation and the conductivity of manufactured radial joints require separate lifetime qualification.

## Reproduction and scope

Run `python3 calculate.py` in this directory. `results.json` records 15 passing unit/interpolation/geometry checks and hashes its input and script. This does not rerun or modify the companion transient model. Source PDFs and rendered inspection images remain under `output/catcher-followup/heat-spreading-evidence/`; only factual extracts, URLs and SHA256 provenance are archived here. The scratch names `ornl-6188316.pdf` and `snead-annealing.pdf` are legacy filenames: their actual authors and institutions are correctly recorded in `data.json`.

The useful design conclusion is to keep oriented graphite around **300–400 W/mK** as a concrete candidate and to test its assembled, irradiated conductivity. The 1200 K screen requires additional heat-spreading capability. The separate heat-pipe evidence review assesses that alternative's demonstrated operating temperatures and device limits.
