# ✨ Positive-prestress support for thinner receivers

A 1000 kg/m² receiver is a live support candidate. The former negative weight in the third stream family was a limit of the **zero-prestress** construction. A cold tension frame can add the inward force needed to make all three stream families positive, within the existing 20 kg/m² frame allocation under explicit strength, stiffness and joint requirements.

`receiver_prestress.py` holds the archived median-carbon128 quartet, pulse rate and source/service budgets fixed, then asks the plant ledger for comparison masses at 500–2000 kg/m². Its 99% capture input provides comparison masses only; lower capture and radiation heating of cold services belong to the parent's separate thermal/electrical scan. No photon-capture performance is inferred here.

Let a=GM/(2R³), n²=GM_sun/(5 AU)³, Σ=M/(4πR²), and T be isotropic frame tension in N/m. The added inward acceleration is `2T/(ΣR)`, so define a_eff=a+2T/(ΣR²). At common stream angular speed ω the three surface densities are

\[
\Sigma_i={\Sigma\over2\omega^2}(a_{eff}+4n^2,\ a_{eff}-2n^2,\ a_{eff}-4n^2).
\]

Choosing a_eff≥4.4n² supplies an explicit positive margin above 4n². This is a design choice, not a minimum proven robustness margin. The script also retains at least 10,000 N/m ordinary pretension for thicker candidates.

| Receiver allocation | Comparison installed mass | Frame tension | Biaxial prestrain at E·h=10⁸ N/m, ν=.2 | Cold stress if 8 kg/m² is load-bearing material |
|---|---:|---:|---:|---:|
| 500 kg/m² | 3.627×10²² kg | 701,440 N/m | 0.561% | 175 MPa |
| 1000 kg/m² | 4.578×10²² kg | 450,765 N/m | 0.361% | 113 MPa |
| 1500 kg/m² | 5.530×10²² kg | 19,280 N/m | 0.0154% | 4.82 MPa |
| 2000 kg/m² | 6.482×10²² kg | 10,000 N/m | 0.008% | 2.50 MPa |

The last column is an equivalent distributed layer at density 2000 kg/m³, thickness 4 mm. Its effective in-plane modulus would be 25 GPa for the chosen E·h. It consumes 8 of the existing 20 kg/m² frame allowance, leaving 12 for joints, connections, load spreading and redundancy. A particular lattice must achieve the effective modulus, isotropic tension and allowable stress simultaneously; the equivalent layer does not certify a fabricated truss. For comparison, retaining only the old ideal 1 mm, 100 GPa, 2 kg/m² load path would require 451 MPa at the 1000 kg/m² point. Reassigning the frame allowance is preferable to treating that ideal minimum as a fixed mass limit.

At 1000 kg/m², added inward pressure is 0.000820 Pa. Stream family densities are 0.04853, 0.01387 and 0.002311 kg/m²; their combined mass is 9.839×10¹⁷ kg. This is 84.6% of the plant's conservative scalar support-stream allowance. At the inherited 0.1% kinetic-energy loss per revolution, the explicit streams dissipate 0.00312 W/m². No guide-package, source or processing mass is reduced merely because the streams are lighter. The frame is already inside services and is not charged twice.

The static construction alone is insufficient. The new global affine calculation includes the second variation of the prestress energy T A(F), **and** the existing frame elastic energy. For F=I+E the area expansion's quadratic coefficient is

\[
{2\operatorname{tr}(E)^2+2\|E\|^2-3\operatorname{tr}(E^2)\over15}.
\]

The normalized prestress energy is 3T/(ΣR²) times this coefficient. The actual self-gravity coefficient remains a; only the stream support uses a_eff. Independent quadrature of the nonlinear surface area, homoeoid gravity and guided-ellipse perimeter reproduces the second variation within 5×10⁻⁸ in units of a. The first variation also closes, checking static equilibrium. A separate agent independently derived the area coefficient and Laplace force.

With the same E·h=10⁸ N/m elastic frame and 30-, 100- or 300-day attitude-feedback choices, every tested 9DOF case has no growing linear eigenmode. The heaviest local stream/body ratio is 0.616 of the earlier local-control reference at 500–1500 kg/m², within its tested 0.5–2 range. The light third family needs its own gain setting; a result for the heaviest rail is not a proof for every carrier and network mode. Sensor noise, real-space control, local trim dissipation, joints, construction prestressing and fault redistribution retain their earlier qualification requirements.

These results remove a spurious hard mass floor at the zero-prestress positive-weight boundary. They do not select a photon capture fraction or prove a complete thinner plant. The existing finite-buffer and local plate checks, renewed skin/contact assessment, and the separately charged destination of leaked radiation must be applied to each thermal candidate. The script is cheap to rerun at a different actual installed mass; a changed capture/heat-routing ledger must not silently inherit the numbers in this comparison table.
