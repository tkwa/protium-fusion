# Ion rotation, encounter curvature and screening sensitivity ✨

The reference titanium plasma rotates at minus half the titanium cyclotron frequency. In its rotating frame, the magnetic and Coriolis terms acting on thermal titanium motion cancel exactly. Electric and centrifugal forces remain, and the electron guiding center follows a curved path through this frame. Evaluating those terms gives a more useful encounter test than the laboratory ion gyroradius alone.

Independent laboratory and rotating-frame orbit integrations verify the transformation and conserved rotating energy. At the median reaction radius, sampled titanium/electron relative paths for the warm probe deviate from straight motion by at most 0.59% of one static screening length over one characteristic encounter time. The cold probe gives 6.26%. These sampled comparisons support a local approximation over that interval while leaving slow encounters, density gradients and collective screening to be assessed.

A separate screening sensitivity test shows why that last qualification matters. Reducing the interaction length tenfold changes the warm state's parallel cooling from 29.37 to 18.24 MeV/s. Its perpendicular heating barely changes. Long-range interactions remain significant even after the numerical integration has converged.

## Thermal motion in the rotating equilibrium

The saved source follows the common-temperature, common-rotation Gibbs equilibrium described in [Dubin and O'Neil (1999)](https://doi.org/10.1103/RevModPhys.71.87). Let signed rigid rotation be `w < 0`, laboratory cyclotron frequency be `Omega_i = q_i B/m_i > 0`, and thermal velocity in the rotating frame be `u = v - w zhat cross r`.

Directly transforming the Lorentz equation gives

```
du/dt = (q_i/m_i) E + w(Omega_i+w) r
        - (Omega_i+2w) zhat cross u.
```

The radial force derives from

```
V_eff(r) = q_i phi(r) - m_i w(Omega_i+w) r²/2,
H_rot = m_i |u|²/2 + V_eff(r) = H_lab - w P_phi.
```

Here `P_phi` includes the magnetic contribution to canonical angular momentum. The density satisfies `d ln n_i/dr = F_eff/(kT_i)`. Both this identity and the energy invariant are checked against the saved nonlinear electric field.

For the reference field 137.149 T:

| Quantity | Angular frequency |
|---|---:|
| Rigid plasma rotation | -2.97454e9 rad/s |
| Titanium laboratory cyclotron frequency | +5.94909e9 rad/s |
| Titanium rotating-frame thermal frequency | 0 |
| Proton rotating-frame thermal frequency | +7.18820e9 rad/s |

The electron cyclotron correction from changing frames is 0.0247% in magnitude. Proton thermal motion retains substantial magnetic turning. Titanium's remaining effective force also varies with radius; its density rises in part of the core and then falls through the outer reaction region. The local effective-potential curvature is retained in the orbit check.

## Finite trajectories and spatial ordering

The diagnostic encounter time is

```
tau = b / sqrt(U² + vz,e² + vth,i²).
```

It uses one ion thermal speed, rather than a sampled relative-speed distribution. `b` is either 0.1 or 1 static screening length. The warm electron has the same parallel velocity as the 3467.60505 eV probe in [HELICAL-RESOLUTION.md](HELICAL-RESOLUTION.md). The cold probe has zero parallel velocity. The electron guiding center follows its leading circular E-cross-B orbit; axial motion is linear.

The ion's laboratory Lorentz orbit and its rotating-frame orbit are integrated independently through the nonlinear radial field, both forward and backward from the selected radius. Initial transverse thermal velocities comprise zero and eight directions at the ion thermal speed. All 58 coordinate, density-gradient and invariant checks pass. The largest normalized verification residual is about 1.02e-10.

At the median reaction radius, 4.09975 mm, the static screening length is 0.56547 mm. The maximum sampled departure of the ion/electron relative path from its initial straight tangent, divided by `b`, is:

| Ion and electron probe | b = 0.1 screening length | b = 1 screening length |
|---|---:|---:|
| Titanium, cold | 0.625% | 6.26% |
| Titanium, warm | 0.0589% | 0.589% |
| Proton, cold | 1.32% | 13.2% |
| Proton, warm | 0.603% | 6.02% |

This is a sample of unperturbed relative paths. It does not integrate an impact-parameter ensemble or bound arbitrarily slow relative encounters. The force-correlation correction also depends on the frequency and wavevector weighting, so a position percentage cannot be applied directly to a collision rate.

Spatial homogeneity presents a separate limitation. At the median radius, the linear variation `b |d ln n/dr|` across one screening length is 0.101 for Ti and 0.157 for H. At the radius containing 99% of reactions, the Ti value is 2.84; at 99.9%, it is 10.5. The outer source therefore requires an inhomogeneous treatment or an explicitly limited interaction range. Shrinking `b` tenfold reduces those linear variation measures tenfold. Reaction weighting alone does not determine the birth locations of contaminating leptons.

## Removing the uncertain distant interval consistently

The sensitivity screen increases the Yukawa screening wave number from `kappa` to `f kappa`, with `f = 1, 2, 10`. It simultaneously reduces the prescribed total logarithm from 15 to `15 - ln(f)`. This keeps the short-range contribution nearly fixed while removing distant interactions. Holding the total logarithm at 15 would instead transfer the removed interval to the unmagnetized contribution.

The exact filtered logarithm has a useful independent form. For `a = (kappa/kcut)²`,

```
L_filtered = [(1+a) exp(a) E1(a) - 1]/2.
```

It follows by substituting `y = k²/kappa²` in the radial spectral integral. Its small-a leading term changes by `-ln(f)`. At `f=10`, finite-filter corrections change the complementary logarithm by 1.16e-5, or 1.8 parts per million. The validator checks that correction against the analytic expression rather than assuming exact equality.

The matching filter is fixed at `kcut = 2 Omega_e/vth,Ti`. For the warm probe:

| Screening multiplier f | Perpendicular rate, MeV/s | Parallel rate, MeV/s | Total rate, MeV/s |
|---|---:|---:|---:|
| 1 | +3.59139 | -29.36825 | -25.77686 |
| 2 | +3.59144 | -26.01728 | -22.42584 |
| 10 | +3.59323 | -18.24476 | -14.65153 |

For an initially stationary electron, the corresponding total heating rates are 2309.11, 2244.13 and 2093.22 MeV/s. These remain initial-state rates; heating rapidly changes the electron velocity. The removed tensor contribution agrees with the independent distant guiding-center operator to within the recorded small finite-gyroradius corrections. Radial refinement from 512 to 768 nodes at `f=10` agrees to numerical precision. Covariance remains positive. All 48 screening-sensitivity checks pass.

This screen measures sensitivity to suppressing a distant logarithmic interval. A dynamical dielectric response can reshape that interval through collective modes and damping, rather than simply remove it. The table therefore supplies an uncertainty test, not an upper or lower physical bound.

## Consequence for the investigation

The laboratory ion gyroradius was too crude a discriminator for titanium in this rotating source. The transformed dynamics support the short-encounter approximation more strongly in the central warm-electron region. Density gradients, slow encounters, screening, total-logarithm choices and finite electron energy still constrain its accuracy. A source transport calculation should retain those sensitivities and the correlated tensor.

The next practical discriminator is lepton supply and collection. Nuclear radiation, residual-gas ionization, atomic charge changes, injection and wall emission determine the actual electron and positron population. Their complete removal cost, including voltage circulation and hardware, must fit within the conditional plant's remaining 1.117 MeV per primary titanium event. Positrons need a mechanism compatible with the positive end barriers. No new plant load, stationary distribution, residence time or future breakthrough date is selected by this checkpoint.

Reproduce with `rotating_collision_scales.py`, `validate_rotating_collisions.py`, `screening_sensitivity.py` and `validate_screening_sensitivity.py`. Each accepts `--root` for the saved source; output and data flags are described by `--help`. The manifest records inputs, results, hashes and terminal timings.
