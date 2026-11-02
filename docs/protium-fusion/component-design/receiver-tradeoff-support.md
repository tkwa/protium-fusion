# ✨ Support and bulk heat storage at the exact lower-capture comparison masses

This follow-up applies [receiver-prestress.md](receiver-prestress.md) to the parent's exact median-carbon128 rows with an additional 0.005 missed-hot-energy fraction and **all missed hot energy deposited in cold equipment**. It preserves the supplied plant masses, pulse rates and energy destinations. No 99% hot-capture requirement is imposed. The scalar cold-equipment energy charge still needs a physical shielding and heat-routing design.

`receiver_tradeoff_support.py` reads `/home/tkwa/src/protium-fusion/output/current-rate128-receiver-tradeoff.json` and the archived `mixture-design/current-rate128-reduced-inputs.json`. Output and source hashes are in `receiver-tradeoff-support.json`. It performs no hydro or repeated plant optimization.

| Receiver | Installed mass | Hot fraction | Required isotropic frame tension | Biaxial prestrain | Stress in 8 kg/m² equivalent cold load layer |
|---|---:|---:|---:|---:|---:|
| 500 kg/m² | 3.62306×10²² kg | 0.89143 | 702.0 kN/m | 0.562% | 175.5 MPa |
| 750 kg/m² | 4.10009×10²² kg | 0.92748 | 599.3 kN/m | 0.479% | 149.8 MPa |
| 1000 kg/m² | 4.57680×10²² kg | 0.95082 | 451.4 kN/m | 0.361% | 112.8 MPa |

The tension gives a_eff=4.4n², maintaining a positive margin for all three stream families. The stated prestrain uses E·h=10⁸ N/m and ν=.2. The equivalent load layer has density 2000 kg/m³ and leaves 12 kg/m² of the retained 20 kg/m² frame allocation for connections and related hardware. These are simultaneous material/joint requirements, not measurements of a completed lattice.

Stream mass uses 95.5%, 89.7% and 84.6% of the corresponding existing scalar support allowance. The inherited 0.1%-per-revolution loss assumption gives 0.00247, 0.00279 and 0.00312 W/m². No guide or frame package is silently removed. All three thickness cases have no growing eigenmode in the extended 9DOF affine model for each of the 30-, 100- and 300-day attitude-feedback choices. Local guide feedback and the full discrete network retain their separate qualification requirements.

The finite bulk heat-capacity check uses actual nominal128 smooth kinetic energy, endpoint internal energy, escaped prompt accumulator and the isotope ledger's subsequent non-neutrino decay energy. **All** of this heat is assigned to the hot absorber for the capacity upper screen, even though the supplied plant sends part into cold equipment. No apparent hot-buffer benefit is obtained by forgetting leaked energy elsewhere.

| Receiver | Graphite thickness | Bulk diffusion time at k=50 W/m/K | Sum-of-stages bulk temperature swing bound | Largest single-stage uniform rise |
|---|---:|---:|---:|---:|
| 500 kg/m² | 0.250 m | 1.25 h | 156.8 K | 68.8 K |
| 750 kg/m² | 0.375 m | 2.81 h | 104.5 K | 45.9 K |
| 1000 kg/m² | 0.500 m | 5.00 h | 78.4 K | 34.4 K |

Density is 2000 kg/m³ and heat capacity is 1800 J/kg/K. With each stage periodic and extraction held at its mean, its zero-mean stored energy ranges by no more than one complete pulse energy. Summing those ranges gives the conservative bulk bound ΣE_pulse/(Σ_receiver c_p). The bound concerns the uniform bulk mode; it neither adds independent surface maxima nor establishes the hottest surface temperature. A finite bulk operating swing has to be allowed. Holding the back face exactly at 1000 K instead would demand much higher short-term heat pickup, so that boundary condition cannot be used to obtain free buffering.

At the actual rates, full-retained mean incident heat is 111.6–118.5 W/m². A steady one-dimensional front-to-back gradient would be only 0.56–1.18 K. The one-day spacing allowance for large events occupies 19.1–20.2% of calendar time, a necessary capacity check rather than a completed scheduling proof. Material pulses and their local skin/contact resistance must still be handled separately; [surface-damage.md](surface-damage.md) gives that qualification. Its current128 annual exposure at the higher .99-capture-family rate conservatively covers the lower rates used here.

Neither support nor finite bulk capacity rules out 500 or 750 kg/m² in this screen. The thinner choices demand larger allowed temperature swings and a stronger prestressed frame than the 2000 kg/m² comparison. Maximum front temperature, fatigue/contact performance, transient photon and secondary-electron transport through the cold equipment, and realizable guide control remain the architecture discriminators. No complete plant or preferred receiver thickness is selected by this memo.
