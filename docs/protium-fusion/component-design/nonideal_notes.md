# ✨ Matched-energy Coulomb pressure sensitivity

This test asks whether the initial local-pressure correction can overturn the preferred fuel mixtures. It adds a consistent, deliberately simple Coulomb free-energy term to the saved fully ionized EOS, then repeats the full 27-isotope,64-zone hydrodynamics at the same initial total prepared energy, density/temperature profiles and physical duration as the uncorrected controls. The eight cases are A-H1/A-H3, B1-H1/B1-H2, B2-H2/B2-H2.5 and C13-H1/C13-H4. Grey opacity remains 0.001 m²/kg. At most two integrations run concurrently.

For initial abundances Y_i, define the frozen charge-weighted moment M=ΣY_iZ_i^(5/3)/ΣY_iZ_i. The stress free energy per mass is

F_C=u_C=−0.9(e²/m_u)[4π/(3m_u)]^(1/3) M ρ^(1/3)Y_e^(4/3),

with e² understood as e²/(4πε₀) in SI. The implementation uses cgs consistently with FastEOS. It gives P_C=ρu_C/3, ∂u_C/∂ρ=u_C/(3ρ), ∂P_C/∂ρ=4P_C/(3ρ), and ∂u_C/∂Y_e=4u_C/(3Y_e). There is no temperature derivative, entropy or heat-capacity correction. The adiabatic exponent is recomputed from the full derivatives. The added term satisfies the thermodynamic identity P_C=ρ²∂u_C/∂ρ and exactly matches the source solver's composition-energy identity when M is frozen.

This is the ion-sphere strong-coupling scale discussed in [Potekhin and Chabrier's ion-liquid EOS work](https://doi.org/10.1103/PhysRevE.62.8554). It is too attractive at weak/moderate coupling and freezes the charge moment as nuclear composition changes. Reaction-rate screening remains unchanged. Consequently this is a controlled hydrodynamic stress test, rather than a complete evolving-composition Coulomb plasma calculation or a guaranteed conservative bound on yield. An independent reviewer confirmed the derivative/accounting construction before the runs.

The initial negative internal-energy correction is balanced by a small increase in convergent kinetic energy. New mixtures also retain the original comparison's internal-energy change. The resulting total initial energy equals the original uncorrected reference to floating-point precision. The corrected center sound speed changes only the internal time unit; the requested duration is converted so that the physical endpoint remains identical. Frozen solver, network and baseline files are unmodified: each process substitutes its own EOS factory and initial composition.

## Validation and evidence

Finite differences verify the added density-energy, density-pressure and composition derivatives to approximately 5×10⁻¹¹ relative; exact Maxwell/composition identities agree to a few 10⁻¹⁶. An independent SI ion-sphere normalization agrees to 4.6×10⁻¹⁰, reflecting the stated constants. Every runtime EOS evaluation monitors positive total pressure, heat capacity and adiabatic sound-speed squared. The independent saved-trace audit additionally checks positive isothermal compressibility.

The comparison reconstructs kinetic and corrected internal energies directly from saved shell states, and reads the three cumulative energy channels with their saved units. It checks initial/final energies, physical duration and all sampled conservation residuals rather than accepting only the reported final residual. A postprocessed instantaneous charge moment supplies a separate scale for the approximation being omitted; it does not silently change the simulated trajectory.

| Stage and H:catalyst | Control nuclear gain Q | Coulomb-stress Q | Relative change |
|---|---:|---:|---:|
| A H1 | 8.19634 | 8.23894 | +0.520% |
| A H3 | 10.28882 | 10.27940 | -0.092% |
| B1 H1 | 6.38197 | 6.50204 | +1.881% |
| B1 H2 | 7.54562 | 7.55905 | +0.178% |
| B2 H2 | 10.55503 | 10.50609 | -0.464% |
| B2 H2.5 | 10.79986 | 10.76176 | -0.353% |
| C13 H1 | 15.53898 | 15.54164 | +0.017% |
| C13 H4 | 23.65567 | 23.59034 | -0.276% |

All eight requested physical endpoints completed. The largest absolute change in stage nuclear gain is 1.881%; each selected mixture retains its stage-Q advantage over its reference. The maximum energy-conservation residual over all saved times is 3.83×10⁻⁶ of prepared energy. The eight independent audits include 64 endpoint/energy identities, alongside the six EOS algebra/normalization checks. Pressure, heat capacity, isothermal compressibility and adiabatic sound-speed squared remain positive at all saved states.

The frozen-moment omission is also quantified: substituting the instantaneous charge moment only in the postprocessed strong-coupling energy changes it by at most 0.885% of initial prepared energy during B2-H2.5's burn, and 0.089% at that endpoint. This provides the scale of one omitted term; it does not evaluate the resulting trajectory or screened rates. The observed yield changes are deterministic scenario differences. A credible interval for true yield would require an explicit EOS/rate uncertainty model.


## Closed-cycle energy comparison

`nonideal_cycle.py` reuses the generalized nuclear ledger on all 16 combinations of these stage choices, with perfect eventual H/CNO return, 60% driver efficiency, 60% heat conversion and 99% energy capture. All 32 control/stress balances are feasible and their entire ranking by pre-auxiliary electrical energy per helium remains unchanged. The A-H3/B1-H2/B2-H2.5/C13-H4 candidate changes from 12.567506 to 12.566413 MeV/He, a decrease of 0.00870%. This includes the altered full isotope output rather than interpreting the named single-stage target yield as the complete cycle.

This metric's best B2 choice is H2 in both cases, slightly ahead of H2.5; the parent's plant-specific-power candidate can favor H2.5 because pulse timing and apparatus mass also matter. The saved-ledger check has not rerun source, receiver, recovery or plant-mass optimization. Its result is that this particular thermodynamically consistent pressure stress does not overturn the nuclear-energy mixture ordering.

## Reproduction and downstream use

`nonideal_eos.py` owns the EOS correction and algebra checks. `nonideal_run.py` runs one case using the archived corrected solver. `nonideal_batch.py --out /home/tkwa/src/protium-fusion/output/nonideal-runs` runs the specified eight cases, two at a time. `nonideal_compare.py` performs the independent trace audit; its optional `--archive-dir` copies numerical result records while correcting the inherited solver's generic no-Coulomb limitation text. Raw simulation outputs remain preserved. The copied metadata explicitly records that distinction.

The runtime is the shared Python 3.10 interpreter with NumPy, SciPy, pynucastro and Numba under `/tmp/protium-python`. The batch script sets the library path, a writable Matplotlib directory, one BLAS thread per process and a separate Numba cache below the scratch output. Scripts and results record source hashes. The corrected result snapshots are in `nonideal_records/`; full traces remain in the scratch output and are hashed by the comparison.

The generalized nuclear ledger can use the saved isotope outputs and energy channels. A receiver or trace reducer that reconstructs internal energy must instantiate `NonidealEOS` with the saved frozen moment and strength. Reusing original FastEOS for those records would re-add the correction and fail the intended energy accounting. These numerical sensitivities add no new claim about cold-feed preparation, multidimensional assembly, photospheric transport or receiver hardware.
