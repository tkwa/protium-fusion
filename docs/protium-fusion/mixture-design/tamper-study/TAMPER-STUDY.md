# Second-generation pulse design: tamper study ✨

Add a tamper comparison to the hydrogen-enriched, staged CNO pulse design. The aim is to increase burning time and reduce pulse mass while retaining a positive whole-cycle energy and material balance. The tamper concepts below are proposed candidates; no tampered pulse has yet been simulated, and no material or layer thickness has been selected.

## Motivation and baseline

The existing second-generation mixture and component design uses staged, optically thick, prescribed converging pulses, with a 30 m reference radius for the major pulses. See [mixture design](../README.md), [pulse methods](../methods.md), and [component design](../../component-design/README.md). Radius in these calculations is a reference scale, not the outer radius of the evolving fuel.

A new 64-zone radius sensitivity used the earlier N14 + p stage, full reaction network, corrected compression-cap viscosity and the prescribed 30 m hot inflow. Each run scaled every initial radius while preserving the initial density, temperature, composition and velocity profile. The fuel input energy consequently scales as radius cubed. The model uses LTE and a constant gray opacity of 0.001 m²/kg. It is a comparison for motivating a tamper, not a run of the selected enriched fuel and not a result for a tampered capsule.

| Reference radius | Outer radius near peak temperature | Fuel mass | Capture product per initial N14 | Prompt capture nuclear energy / fuel input |
|---:|---:|---:|---:|---:|
| 30 m | 77.8 m | 1.131e13 kg | 0.5394 | 8.196 |
| 10 m | 20.8 m | 4.189e11 kg | 0.2174 | 3.304 |
| 3 m | 6.56 m | 1.131e10 kg | 0.005853 | 0.08895 |
| 1.5 m | 3.25 m | 1.414e9 kg | 0.002141 | 0.03253 |

The 3 m run still has a naive center-to-surface Thomson optical depth of approximately 4.5e6 near its maximum sampled temperature. The decline in capture yield therefore cannot be attributed to low scattering depth alone. A shorter hydrodynamic time and a changed temperature and density history are plausible reasons; the sweep does not isolate their individual contributions. The local free–free absorption depth of about 825 quoted in discussion was a half-cell diagnostic at the 30 m baseline's hottest point, not a requirement for the complete pulse.

The 30 m case is the previously saved corrected-viscosity regression. The 10 m, 3 m and 1.5 m cases were run for four reference sound times with relative tolerance 1e-6 and absolute tolerance 1e-11; all reported completion, with energy residuals below 1.5e-7 of initial energy. The 30 m regression used relative tolerance 1e-7. These single-resolution outputs do not establish convergence or a realizable cold driver. The saved [probe](radius_probe.py), JSON results and spatial traces in this directory preserve the new runs; the probe uses the repository’s pulse solver and input profile and requires the scientific Python environment documented by the research. Run `python docs/protium-fusion/mixture-design/tamper-study/radius_probe.py 3 --output /tmp/protium-radius-3` from the repository root. No energy from subsequent beta decay is included in the stated capture-stage ratio.

## Candidate configurations

Compare the following with an untamped control. Define initial mass, pressure, density, radial velocity and composition throughout each layer; assemble them before peak compression and evolve the entire fuel–tamper system. Keep the initial fuel profile fixed for a first mechanism comparison, then compare again at equal **total** delivered energy including the tamper.

1. **Low-Z inertial layer:** a hydrogen-rich or carbon-rich sacrificial outer layer. It provides mass and mechanical impedance with modest high-Z contamination. Hydrogen and carbon might be recoverable within the existing fuel cycle, but their isotopic, thermal and processing costs must be counted.
2. **Bismuth layer:** a converging high-Z bismuth layer surrounding the fuel. Evaluate its inertia, shock response and radiation absorption/re-emission using an ionization-dependent equation of state and opacity. Cold solid density and cold X-ray attenuation are preparation references, not hot-plasma properties. Bismuth recovery, activation and isotope changes need a separate material ledger.
3. **Low-Z buffer plus bismuth:** put a carbon-rich or hydrogen-rich buffer between the reacting fuel and bismuth. Sweep buffer and bismuth masses independently to see whether the buffer improves interface behavior and limits bismuth transport into the burn. Establish its effect rather than assuming it prevents mixing.
4. **Staged inflow / fuel tamper:** retain an outer, cooler CNO–hydrogen layer whose inward momentum confines the central burn and can itself later be processed as fuel. This tests whether useful fuel can supply confinement without a separate heavy-material inventory. Track its compression, burn and recovery separately from the core.

For the initial survey, use shell-to-fuel mass ratios 0, 0.03, 0.1, 0.3 and 1 as exploration points, plus at least two choices of shell inward speed and arrival timing. These are search points, not recommended specifications or probability intervals. Evaluate nitrogen capture first at 3 m and 10 m reference radii; compare with the 30 m untamped case. Then transfer promising ideas to the actual hydrogen-enriched second-generation mixtures and other reaction stages.

## Screening calculation

For intuition only, a thin shell initially at rest under constant pressure satisfies Δr ≈ p t²/(2Σ), with Σ the mass per unit area. At p=1e18 Pa, t=5 µs and Δr=2 m, Σ≈6.25e6 kg/m²; assigning that column at 6.56 m radius gives ≈3.4e9 kg, about 30% of the 3 m fuel mass. Pressure was selected from the untamped outermost resolved cell near peak temperature; neither that pressure nor a stationary thin shell describes the evolving interface. Cold bismuth at this column would not form a thin layer at this radius. This calculation is a parameter-scale prompt, not a physical mass prediction or bounds on required confinement.

An illustrative kinetic-energy scale for 3.4e9 kg at 1e6 m/s is 1.7e21 J, compared with 3.493e22 J of fuel input in the 3 m untamped control. Count the tamper's actual initial kinetic and thermal energy, compression work, driver inefficiency, subsequent product recovery, and material processing. Some outgoing energy may be recovered; do not assume a recovery fraction or count it twice.

## Required model and decision criteria

- Add distinct fuel, buffer and tamper materials to radiation hydrodynamics with material-specific EOS, opacity and electron/ion thermodynamics. Resolve interfaces and follow changing ionization and radiation transport; benchmark the coupled solver against shock-tube and energy-conservation controls. A gray LTE survey is a first screen, followed by frequency-dependent transport at consequential boundaries.
- Record fuel density and temperature histories, interface pressure and velocity, shock timing, burn fraction, escaped radiation and particle energy. Compare paired runs both at fixed fuel state and fixed total input. Repeat at multiple zoning and timestep settings and retain failed or stopped runs.
- Model or bound interface instability and material transport; high-Z admixture changes opacity, bremsstrahlung, dilution and heat capacity. The dependence in an opaque pulse must be evaluated with radiation transport, not an optically thin loss formula alone. Two- or three-dimensional instability and assembly work are required before claiming a feasible configuration.
- Compute prompt stage gain as nuclear energy less explicitly escaped weak energy divided by **fuel plus tamper preparation and delivery energy**. Account separately for captured outgoing energy, delayed beta decays, catalyst return, additional processing and electrical efficiencies in the whole-cycle ledger. Do not promote a capture-stage ratio to net electric gain.
- The exploratory 3 m control converts 0.5853% of N14 with a prompt capture nuclear/input ratio of 0.08895. If energy release per conversion remains approximately constant, approximately 6.6% conversion would yield a ratio of one against the *old fuel-only* input. The tamper adds costs, so this is a diagnostic threshold, not the plant break-even condition.
- Recompute pulse mass, frequency, collection load, chemical/isotope handling and mass replacement for the second-generation plant. Check whether a lower fuel pulse mass survives the cost of tamper inventory and a higher firing rate. Revisit the other reaction stages and their waiting-point inventories.

The immediate acceptance criterion for a candidate is a converged comparison showing a higher capture-stage gain **at equal total delivered energy** than the untamped control, followed by a closed whole-cycle energy and material account. No subjective 80% credible interval is assigned to tamper performance: there is not yet a calibrated physical model or experimental evidence from which to justify one.

## External references

- [Sandia, inertial confinement fusion and metal liners](https://www.sandia.gov/pulsed-power/inertial-confinement/).
- [NIST bismuth mass attenuation coefficients](https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z83.html) and [material density](https://physics.nist.gov/PhysRefData/XrayMassCoef/tab1.html); cold-material reference only.
- [Fuel–shell mix and fusion performance](https://www.osti.gov/servlets/purl/1172490).
