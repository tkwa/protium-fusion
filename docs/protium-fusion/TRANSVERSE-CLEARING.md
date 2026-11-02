# Transverse lepton clearing and ion response ✨

The next proposed collector applies a rotating transverse electric field to move trace electrons and positrons toward the channel wall. Two calculations now constrain this idea: collisionless ion trajectories in the saved equilibrium, and the collective motion of two rigid ion profiles. The first supplies a bare kinetic response; the second shows why induced fields matter. Neither calculation supplies a collection efficiency or changes the plant's net output.

The reference remains the 1.6 MeV titanium source at 137.149 T effective field, with a 20.732 mm channel radius. The existing plant allocation leaves 1.11714 MeV per primary titanium event for extra electrical losses before it falls below 1 mW/kg.

## Guiding-center motion and capture

For an applied potential

`phi1 = -E0 r cos(theta - phase(t))`,

the leading guiding-center equations in a prescribed uniform axial magnetic field are

`dr/dt = -(E0/B) sin(psi)`,

`dpsi/dt = omega_E(r) - omega_d(t) - E0 cos(psi)/(B r)`,

where `psi = theta - phase(t)` and `omega_E = -Er/(Br)`. The leading electric drift is independent of charge sign. These equations omit collisions, finite gyroradius, radiation, induced fields and the collection boundary.

They give a useful speed limit, `|dr/dt| <= E0/B`. They also retain the phase equation: a rotating field does not automatically move every particle outward at that speed. An actual chirp must capture enough of the birth distribution, preserve the required phase, and remove particles after they approach the wall. Collisionless phase-space conservation limits assumptions about capture into a narrow resonance.

From the median reaction radius, the ideal radial crossing-time lower bounds are 0.228 s, 22.81 ms and 2.281 ms at 10, 100 and 1000 V/m, respectively. These are kinematic lower bounds, not residence predictions. Continuous births, collisions, drift back toward the core, intermittent operation and failed capture all affect the actual population.

## Frozen ion trajectories

`ion_dipole_spectrum.py` samples both equilibrium ion species and follows their exact nonlinear radial force in the common rotating frame. Length is measured in `sqrt(T mp)/(eB)` and time in `mp/(eB)`. With dimensionless mass `M`, charge `Z` and signed rotation `w`, the equation is

`r'' = -(Z/M) grad(psi0) + w(Z/M+w) r - (Z/M+2w) zhat cross r'`.

The conserved rotating energy is

`Hrot/T = M u²/2 + Z psi0 - M w(Z/M+w) r²/2`.

The integrator splits the nonlinear potential kicks around an exact magnetic helix. Initial radii follow the saved number-density cumulative distributions; angles are uniform and rotating velocities are Maxwellian. The production runs use two independent scrambled Sobol ensembles of 16,384 particles per species, 2 ps steps, 131.072 ns duration and 16 ps saved samples. A 4096-particle control halves the orbit step while preserving the sample cadence.

For `z = x + iy`, the two-sided Hann-windowed spectrum is

`S_z(delta) = <|integral g(t) z(t) exp(-i delta t) dt|²> / integral g(t)² dt`.

The bare circular-field response is

`P_rot/(length E0²) = N_line q² delta² S_z(delta)/(4 kT)`.

Here `delta = omega_d - w_physical` is the drive frequency seen in the ion frame. An independently forced harmonic oscillator verifies the factor 1/4 and polarization convention. The saved source spectra omit the ions' perturbed electric field. Their absorption cannot be assigned directly to an externally driven reactor without solving for the total field.

At laboratory angular frequency -1.51e9 rad/s, interpolation of the refined spectrum gives 3.3144e-4 and 5.9444e-5 MeV per primary event per `(V/m)²` for protons and titanium. A prescribed uniform total field of 100 V/m would therefore correspond to about 3.91 MeV of rotating-frame heating per event in this bare calculation. The coefficient is a diagnostic of overlap with ion motion; it is not an established RF operating cost.

The maximum saved rotating-energy error is 4.90e-4 T for protons and 8.61e-5 T for titanium. Halving the step in the matched 4096-particle control reduces the maximum error by a factor of four. Parseval and an independent velocity second-moment identity pass. In three broad frequency bands covering the clearing sweep, the two refined ensembles differ by 0.16–5.82%. The initial 4096-particle proton tail differed by 21.1%; the larger ensemble changes that tail by 16.3% and reduces the independent-seed difference to 4.54%. Narrow spectral bins and the low-density tails have no percent-level convergence claim. The frequency resolution is 4.794e7 rad/s.

## A collective check

`rigid_ion_dipole.py` allows the proton and titanium profiles to translate by different small transverse displacements. This closure retains mutual electrostatic forces and the conducting cylinder's dipole image field. It omits shape changes and kinetic damping.

For line charges `lambda_i`, line masses `M_i`, and laboratory cyclotron frequencies `Omega_i`, the complex centroid equation is

`M_i z_i'' + i M_i Omega_i z_i' = sum_j A_ij z_j + lambda_i E(t)`.

The force matrix is

`A = K [[1,-1],[-1,1]] + lambda lambda^T/(2 pi epsilon0 R²)`,

`K = q_H q_Ti integral n_H n_Ti d²r/(2 epsilon0)`.

The overlap coefficient follows by translating one circular density through the other's field and averaging its field gradient. A separate finite-displacement two-dimensional force integration agrees within 7.3e-7 relative at the first displacement and 1.9e-7 after halving it. The mutual forces cancel in the summed centroid equation. The cylindrical image term follows from the leading dipole solution of the grounded boundary.

For a mode proportional to `exp(i omega t)`, the roots satisfy

`det[diag(M_i omega(omega+Omega_i)) + A] = 0`.

All four roots in this closure are real. Their clockwise laboratory frequencies are 1986.354, 801.759, 226.440 and 23.138 MHz. This is stability of four rigid centroid coordinates, not full kinetic stability.

The lepton electric-drift frequency falls from 248.542 MHz near the axis to 22.731 MHz at the wall. Consequently an outward frequency sweep crosses two centroid-model poles:

| Clockwise mode frequency | Matching radius | Reaction fraction inside that radius |
| --- | ---: | ---: |
| 226.440 MHz | 5.85775 mm | 91.4209% |
| 23.138 MHz | 20.54849 mm | Approximately 100% |

The first mode mainly displaces the two species in opposite directions; the second is nearly common translation. A crossing calls for a coupled response calculation. It does not by itself determine damping, ion loss or a maximum usable drive.

The equal-charge-to-mass control reproduces the common centroid's analytic frequencies: internal forces cancel and a uniform field drives only common translation. This illustrates the limitation of summing independent orbits in a fixed self-field. Collective motion can change the response substantially even when those orbits are integrated accurately.

## Heat, angular momentum and the electrical ledger

For a single rotating perturbation, laboratory work and canonical torque obey

`P_lab = omega_d tau`,

`P_rot = P_lab - w_physical tau = delta tau`.

The review by [Danielson et al. (2015)](https://doi.org/10.1103/RevModPhys.87.247), especially sections II.D–E and the rotating-wall discussion, provides the equilibrium and thermodynamic context. Laboratory rotating-wall experiments also show mode coupling and heating; their measured performance is not transferred to this MeV source.

In the inner clearing band, `w_physical < omega_d < 0`, so positive rotating heating accompanies negative laboratory work: the drive extracts ordered energy while the plasma gains entropy. Restoring its original angular momentum requires work. For rotating heat `H`, ideal restoration supplies `R = -w_physical H/delta`, and the clearing drive receives `C = -omega_d H/delta`; `R-C=H`.

At -1.51e9 rad/s, these amounts are `R=2.03104 H` and `C=1.03104 H`. With assumed conversion efficiencies `eta_in` and `eta_out`, the electrical input is `R/eta_in - eta_out C`, before RF field losses, electronics, heat recovery and hardware. If this heating can replace existing ion heating, subtract the displaced heater input `H/eta_heater` once. At all three efficiencies 95%, the incremental input is 0.10582 H before waste-heat credit. At 99.9% torque-transfer efficiencies and a 95% displaced heater, the ideal incremental term is -0.04957 H. These are accounting examples; neither efficiency nor usable heating fraction is established for this field arrangement.

The main titanium stage already requires 8.669 MeV per primary event of thermal preparation in the saved source ledger. A useful heating credit must fit the actual species, distribution, replacement throughput and control requirements. Excessive heating cannot receive an unlimited credit. Converter mass and losses must follow the gross opposing transfers R and C as well as their difference.

## Research decision and reproduction

Radial clearing remains a conditional route. Its speed limit alone does not establish capture, and its frozen ion response alone does not establish net electrical loss. The two collective frequency crossings and the heat/torque identity are the useful results of this checkpoint. Full field-response and lepton-capture calculations should be driven by a concrete collector design and a loss allocation, rather than by further arbitrary orbit refinement.

An axial traveling field is another candidate because its phase velocity can be placed above most ion velocities. Its penetration through the cylindrical electrodes, plasma dielectric response, capture from the lepton birth distribution, RF losses and finite end barriers must be included before comparing it with radial clearing. No numerical advantage is assigned here.

The scripts and saved outputs accompany this note. Regenerate the reference spectra with `ion_dipole_spectrum.py --n 4096`; use `--dt 1e-12 --stride 16` for its step control, and `--seed 917` for the independent ensemble. The production spectra use `--n 16384` with seeds 42 and 917 and the corresponding explicit output names. Run `rigid_ion_dipole.py`, `validate_ion_dipole.py` and `plot_ion_dipole.py` in this directory. All 100 mathematical and broad-band comparison checks pass; the retained figure was visually inspected. The manifest records code, data and existing-input hashes. All jobs are complete.
