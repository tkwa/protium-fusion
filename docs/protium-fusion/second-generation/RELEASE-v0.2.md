# v0.2 — Large-N15 thermal design ✨

This release selects `spherical_large_N15_16` for the website's current design:

| Quantity | Selected value |
| --- | --- |
| Net electrical output | 9.8727962 × 10²¹ W |
| Complete installed mass | 3.5554631 × 10²² kg |
| Specific output | 0.2776796 W/kg |
| Helium throughput / earlier selected reference | 16 |
| N15 mass per pulse | 1.1309734 × 10¹³ kg |
| Useful hot source / radiator | 1600 K / 400 K |
| Net heat-conversion package efficiency | 54.5715% |
| Electrical-to-prepared-pulse efficiency | 60% |
| Radiator exterior coverage | 83.5% |
| Serialized chamber occupancy | 28.2% |
| Pulses in the 370-day dispatch | 1879 |

The specific-output gain is 13.36 times the 20.784 mW/kg selected reference, rounded to 13.4. This result uses one reaction centre and thermal conversion. The mass account retains 1000 kg/m² of service equipment, complete peak-rated launch hardware, two-pulse electrical storage, maintenance, spares and installed reserves. Full availability is assumed; 90% availability gives about 0.250 W/kg averaged over the calendar.

The larger spherical N15 burn converts 77.45% at 128 radial zones and 79.13% at 256 zones. Those simulations begin with prescribed hot inflows. Cold source formation, three-dimensional assembly, capsule flight and neutralization, spectral receiver transients, irradiated material lifetime near 1600 K, and a completed thermal/recovery equipment bill remain engineering qualifications. Numerical closure supplies no whole-reactor feasibility probability.

The extended CNO account retains oxygen as recoverable catalyst inventory and sizes return work to actual modeled oxygen production. The hypothetical new oxygen feed equal to 1% of catalyst throughput is outside this design basis. Headline numbers use the closed primary quartet account with held oxygen and replacement carbon. A separately allocated oxygen-return calculation subtracts approximately 0.13% of net export before added apparatus; complete coupled mixed-source recovery remains development work.

## Reproduction and evidence

- [Publication inputs](publication-inputs.json) freeze the selected plant, cycle, N15 record, inherited transport responses, and original source hashes.
- [Publication events](publication-events.csv) contain the actual serialized focus times. [Buffer screen](publication-buffer.json) applies 16-fold pump and outlet capacity to that calendar, with 1000 K gas and 0.5 m³/m² volume; its resolved pressure range is about 0.339–2.069 Pa.
- [Research report](README.md), [oxygen account](oxygen-network.md), [source review](sources.md), and the accompanying scripts preserve the investigation and its numerical qualifications. Their dated statements describe the research run; this release note owns the publication selection.
- [Published plant data](../../../public/blog/ordinary-hydrogen-fusion/selected-plant.json) and [energy balance](../../../public/blog/ordinary-hydrogen-fusion/catcher-interface.json) share the publication-input hash.
- `python3 tools/verify_article_balance.py` checks the selected source identity and energy, mass and isotope closure using the standard library. The original study's `validate.py` supplies its more extensive scientific checks; some source histories retain research-machine paths.
- The website's [simulation page](https://tkwa.me/blog/ordinary-hydrogen-fusion/plasma-simulation/) supplies the four selected raw traces, updated plots and a compact reproduction ZIP. The new N15 trace and that ZIP are also attached to the GitHub release. Their input hashes bind the raw data and reducer.

The website changes preserve prose except Thomas's explicitly approved corrections and the requested version label. The risk table retains its inherited scores and states that their application to the selected operating conditions awaits reassessment. The v0.1 tag remains unchanged.
