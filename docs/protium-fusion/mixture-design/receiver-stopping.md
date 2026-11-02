# ✨ Selected-mixture receiver and stopping checkpoint

The generalized receiver reduction is in `receiver.py`, with one cached result per saved source case in `receiver-results/`. The selected material calculation is in `receiver_stopping.py` and `receiver-stopping.json`. These files use the new mixture traces; the earlier station and receiver packets remain unchanged.

The source rating reconstructs initial internal energy with the actual electron and ion abundances of each mixture. Its launch-time bins retain the earlier source-allocation convention, and an exact local half-cell kinetic allocation is also reported. Outgoing calculations retain the complete mass and distinguish the solver's nodal kinetic quadrature from a smooth, mass-conserving velocity distribution. The named radiation accumulator supplies prompt radiation without a fitted normalization. Prompt and slow-material temperature maxima occur separately and should not be added to create a simultaneous surface temperature. The surface-deposition prompt result is a mathematical reference; finite spectral absorption is assessed in the radiation packet.

For a cell with constant mass per velocity, the finite-radius mapping is

\[
 t=t_0-\frac{\Delta r}{\Delta v}+\frac{b}{v},\qquad
 b=R-r_i+v_i\frac{\Delta r}{\Delta v},\qquad
 \frac{dm}{dt}=\frac{m_i}{\Delta v}\frac{v^2}{b}.
\]

Thus the local energy and momentum fluxes are respectively \(v^2/2\) and \(v\) times the mass flux per receiver area. A constant-flux analytic heat-conduction check and an independent Gaussian-quadrature check of the ballistic mass, momentum and energy Jacobians pass. Each source's K, U and radiated energy is checked against its own record. The `mass999` schedule reserves a separate route for the slowest 0.1% of mass on each leg. The finite `full_resolved` schedule refers to shell representatives; the continuous central cell extends down to zero speed and has no finite literal last arrival.

The selected stopping cases are the fixed-prepared-energy H2 versions of A, B1, B2 and C13, plus C13-H4, at R=1e9 m. The actual endpoint isotope mixture is propagated through its flight using the network's one-body weak decays. Sixty-four Gaussian velocity packets per cell give the range quantiles below; doubling that order changes the reported quantiles by at most 0.21%. The dominant uncertainty is the stopping/material model rather than that quadrature.

| Case | Maximum MeV/u | Graphite CNO 99.9% path, µm | Aluminium, µm | Alumina, µm |
|---|---:|---:|---:|---:|
| A-H2 | 1.186 | 11.34 | 13.18 | 6.38 |
| B1-H2 | 0.548 | 5.74 | 7.28 | 3.69 |
| B2-H2 | 1.641 | 17.71 | 19.65 | 10.30 |
| C13-H2 | 2.396 | 30.65 | 33.28 | 17.35 |
| C13-H4 | 2.730 | 36.46 | 39.33 | 20.43 |

The calculation retains the hydro model’s Newtonian kinetic-energy convention. At the fastest selected endpoint, the relativistic energy correction is 0.44%, smaller than the stated stopping-model uncertainty. The assumed densities are 2000, 2300 and 3900 kg/m³. These are continuous-slowing-down path lengths; projected depth, straggling, roughness, damage and transient material state require separate qualification. The underlying electronic stopping combines the archived measured [IAEA data](https://www-nds.iaea.org/stopping/) with the previously documented interpolation and extrapolation, plus ZBL elastic stopping. Entry charge and equilibration in the target remain uncertainties. Multiplying electronic stopping by 0.5 or 2 is an explicit sensitivity test, not an experimentally established confidence interval.

Nominal 50 µm alumina covers the calculated range of every included charged nucleus in these selected cases. Under the half-electronic-stopping sensitivity, some fast H/He reaches the backing while the CNO remains within the nominal coating. This changes the recovery reservoir; it does not establish physical escape. A thick closed backing also catches particles which pass a thinner coating. Annual batch treatment remains an available extraction schedule, and incomplete extraction leaves retained stock for reprocessing. Its material loading, isotope dilution and chemistry are handled in the material-loop packet.

The nuclear screen integrates an assumed energy-dependent event cross section along the stopping trajectory:

\[
 P_j(E_0)=1-\exp\left[-\int_0^{E_0}\frac{n_t}{\rho}\,
 \frac{\sigma_j(E)}{S_j(E)}\,dE\right].
\]

A constant one-barn sensitivity gives the following CNO event probabilities in parts per million. It represents a response coefficient for an assumed cross section, and cannot be read as either a measured reaction probability or permanent catalyst destruction.

| Case | Graphite | Aluminium | Alumina |
|---|---:|---:|---:|
| A-H2 | 36.5 | 24.0 | 28.4 |
| B1-H2 | 27.5 | 18.0 | 22.0 |
| B2-H2 | 42.6 | 27.3 | 32.4 |
| C13-H2 | 69.3 | 42.7 | 49.1 |
| C13-H4 | 74.9 | 45.9 | 53.0 |

The recorded rows also give all-projectile events per initial catalyst, so a selected cycle can weight them using its actual stage frequency. They do not import the previous cycle's fixed return allowance. An inelastic excitation, a recoverable change between useful catalyst isotopes, and conversion into a nuclide outside the chosen catalyst set have different consequences. Recoil stopping and closed-system recapture also differ from external escape.

The illustrative carbon-on-carbon fusion extension from [Notani et al. (2012)](https://doi.org/10.1103/PhysRevC.85.014607) gives approximately 0.69, 5.34 and 7.04 fusion events per million initial catalysts for B2-H2, C13-H2 and C13-H4 respectively. This inherits the earlier model's measured-fit range and explicitly extrapolated low/high tails. It covers only C12/C13 on carbon. It is not a total CNO nuclear-loss calculation.

Measured N/O heavy-ion fusion data exist and can constrain the remaining calculation: [Kovar et al. (1979)](https://doi.org/10.1103/PhysRevC.20.1305) includes C12+N14, C12+N15 and C12+O16, while [Gomez del Campo et al. (1979)](https://doi.org/10.1103/PhysRevC.19.2170) reports N14+C12 fusion and direct reaction components. The faster new tails now overlap part of these measured energy ranges. This checkpoint has not digitized their excitation functions or resolved all product branches. Aluminium and oxide require their own reaction channels; the constant-barn rows preserve that uncertainty instead of treating every event as lost carbon.

This is an ion-stopping component calculation. Quasineutral plume electrons, secondary-electron escape and return, sheath currents and dielectric charging still belong in the collector heat/current circuit. The collector must account for charge balance and any required return path so that charging does not invalidate its assumed trajectories.

The present result supports continued development of a material collector. It supplies no new necessity for nearly complete upstream electromagnetic braking. The deciding engineering issues remain surface survival, chemical release, local gas capture and a closed material-routing system, with nuclear-channel closure still required for a final fresh-carbon ledger.
