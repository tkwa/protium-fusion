# Small prescribed-hot stages and factory discriminator ✨

29 October 2026. This bounded follow-up tests geometric scaling with the existing spherical 27-isotope model. It supplies completed A/B1/C13 stages, a completed 3 m helium factory, and an explicitly incomplete B2 calculation. It does not model preparation, multidimensional assembly, collection, or engineering performance.

## Results

The reference radius is the radius of the original uniform compressed comparison sphere. It is not the initial inflow radius or the full stagnation radius. Spatially scaling each saved hot inflow preserves its temperature, density and velocity profiles while reducing mass as radius cubed and hydrodynamic time in proportion to radius.

| Stage | Reference radius | Mass, kg | Prompt deposited Q | Target conversion | Peak temperature | Status |
|---|---:|---:|---:|---:|---:|---|
| A, full network | 30 m | 1.131e13 | 8.181641 | 0.538400 | 1.149 GK | Completed |
| B1 | 30 m | 1.131e11 | 6.367205 | 0.600477 | 0.4575 GK | Completed |
| C13 | 30 m | 1.131e13 | 15.479467 | 0.813129 | 1.577 GK | Completed |
| B2, tighter tolerance | 30 m | 1.131e13 | 10.516806 | 0.819339 | 1.638 GK | Failed during late expansion; endpoint only |
| B2, radiation off | 30 m | 1.131e13 | 10.516697 | 0.819332 | 1.638 GK | Failed during late expansion; endpoint only |
| He factory | 3 m | 1.131e12 | 0.000952541 | 0.000349412 | 1.592 GK | Completed |

All completed runs reach four reference sound times. Their absolute energy residuals relative to input are at most 2.39e-5; baryon residuals are below 8e-14. Trace negative mass fractions are tiny numerical errors and are disclosed in every output. The new full-network A gain differs from the earlier reduced-network small A by only 1.13e-6, and conversion differs by 8.59e-8. This is a useful network comparison at this trajectory, not grid convergence. All new runs have 64 zones.

B1's original reference radius was 1,500 m, B2's 300 m, and C13's 100 m. Thus the same 30 m endpoint imposes different duration reductions. No common factor of ten was silently applied. The prescribed temperatures/densities remain those in the original stage inputs. B2 begins with pure C12 plus two protons per catalyst; no C13 admixture is included.

## B2 failure is numerically identifiable

The first default-tolerance B2 attempt raised a shell-crossing exception during a BDF Newton trial. The 1e-7 tolerance retry stops at 35.182 microseconds, after 2.2315 of the requested four sound times, with `Required step size is less than spacing between numbers`. All accepted endpoint shell velocities are outward. The core is about 0.792 GK, while the outermost cell has cooled to 0.301 MK.

The legacy solver writes `termination_reason=final_time` whenever no hot guard fires, even when `success=false`. This packet preserves that raw field and adds `effective_termination_reason=solver_failure`. Consumers must check `success` and the effective status. The saved endpoint is not a completed B2 pulse.

At that cold outer state, extrapolated library coefficients for F18(n,gamma)F19 and Ne21(n,gamma)Ne22 rise to approximately 4.4e68 and 3.1e68 cm³ mol⁻¹ s⁻¹. Neutron and Ne21 abundances contain roundoff negatives of order 1e-30; the solver evaluates its unclipped mass-action polynomial for Newton differentiability. Huge low-temperature extrapolations amplify those negligible numerical traces. The actual dominant cold-cell flow with nonnegative abundances is ordinary N13 decay. The exact coefficient audit is in `B2-low-temperature-rate-audit.json`.

A sensitivity evaluates network coefficients at `max(T,10 MK)` while retaining the physical temperature in the EOS and fluid equations. The weak rates in this network are temperature independent. This is a transparent numerical regularization, not a replacement measurement of cold neutron-capture rates. Charged CNO captures at 10 MK are already negligible on the microsecond timescale; neutron capture requires a physical low-energy treatment if real neutron inventories become relevant. The regularization passes the original stiff region and then raises the existing EOS-table guard when the trial outer temperature reaches about 80,000 K. It therefore does **not** supply a completed result, and no gain is invented for it. Its error and configuration are saved separately.

The accepted B2 tail is already close to capture freeze-out: O14 conversion rises from 0.819230 at 33.247 microseconds to 0.819339 at 35.182 microseconds, only 0.0133% relative. Its gain is already above ten and all shells are expanding. Continued cooling with most released energy retained in expanding material is physically plausible; completing cold-boundary transport and EOS treatment remains necessary to certify the final output. Nothing in this endpoint indicates that useful energy must disappear before collection. The saved tail history is in `B2-tail-conversion.json`.

Turning off radiation also fails in the cold tail, at 34.913 microseconds. Nevertheless, its endpoint gain and conversion agree with the radiating endpoint to roughly 1e-5 relative. This isolates a numerical/cold-boundary limitation from failure to release capture energy. It is insufficient to certify later freeze-out or the final isotope distribution.

The primary database's [rate-fitting documentation](https://reaclib.jinaweb.org/help.php?intCurrentNum=0&topic=rate_fitting) discusses avoiding nonphysical behavior outside fitted ranges. Its indexed text was readable; direct retrieval returned 403 and was not bypassed. No universal 10 MK validity threshold is asserted for every library rate.

## Conditional cycle accounting

`cycle-check.json` applies the existing isotope ledger to completed A/B1/C13, the B2 endpoint, and the existing 30 m factory. `cycle-He3-check.json` substitutes the completed 3 m factory. **Both ledgers assume the B2 endpoint can be frozen and processed.** They do not silently treat its failed solver as a completed stage or prove full-cycle engineering closure.

For the radiating B2 endpoint, 30% electrical-to-delivered efficiency and 40% thermal conversion give:

| Factory | Catalyst loss per stage | Available nuclear heat, MeV/He | Prepared input, MeV/He | Net electricity before other auxiliaries, MeV/He |
|---|---:|---:|---:|---:|
| 30 m | 0 | 24.45010 | 2.39254 | +2.76193 |
| 30 m | 1e-4 | 24.49367 | 2.75855 | +1.70572 |
| 3 m | 0 | 24.45005 | 3.22727 | +0.31337 |
| 3 m | 1e-5 | 24.45438 | 3.58968 | −0.74796 |
| 3 m | 1e-4 | 24.49345 | 6.85651 | −10.31506 |

The smaller factory reduces single-shot mass and energy by 1,000 relative to the 30 m factory, but reduces carbon conversion by another factor of about ten. Its preparation cost per produced carbon therefore increases about tenfold. At 50% driver and 50% conversion, the 3 m factory / 1e-4-loss sensitivity still gives approximately +1.962 MeV/He before other auxiliaries. Driver efficiency, isotope processing and material return must be selected together with the collector's pulse tolerance.

These ledgers retain the original approximate beta-neutrino allowance and count delayed heat separately. Their energy and baryon identities are checked by the existing ledger code. They establish an accounting sensitivity, not a prediction of the chosen future installation.

## Opacity check

`transport-check.json` independently extracts actual saved spatial trajectories using the earlier qualified microphysics diagnostic. In the burning regions, the 1st-percentile cell effective thermalization depths are about 17,014 for A, 230 for B1, 18,403 for the radiating B2 endpoint, 7,287 for C13, and 667,441 for He3. Reducing both scattering and true absorption by 100 leaves these values at 170, 2.30, 184, 72.9 and 6,674 respectively. This is a coefficient stress test, not a certified lower bound.

B1 has the smallest margin. With that severe common suppression, its 99th-percentile absorption/evolution-time ratio is 0.0795; the other completed capture stages remain below 0.0034. The bulk can remain opaque at these scales. Dilute interfaces, resolved-layer conduction, and the cooling B2 outer boundary still require appropriate transport/EOS treatment. The opacity diagnostic does not repair a failed hydro solution.

## Reproduction and provenance

`small_stages.py` reads original stage JSONs and their saved inflows; it writes only to the specified output directory. Run A, B1 and C13 at default radius 30 m; B2 additionally used `--rtol 1e-7 --tag=-tight`, `--no-radiation --tag=-norad`, and `--rate-floor 10000000 --tag=-ratefloor`. The helium run used `He --radius 3`. All stage source files remain unchanged. Python needs the repository's existing numpy/scipy/numba/pynucastro environment; this work used `/tmp/protium-python`.

Each run records source hashes. The packet manifest includes code, EOS table and output hashes. Full traces are retained in `/tmp/protium-small-stages`; explicitly downsampled `*-compact.npz` archives preserve 32 distributed checkpoints, the last eight, and the stored peak-temperature checkpoint. `compact-traces.json` records sampling counts and source hashes. Use full traces to reproduce the original transport quadratures. The saved source and compact output JSONs are sufficient to rerun the simulations. No guard stop, failure or sensitivity is promoted to a qualified physical deployment result.
