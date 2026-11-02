# B2 cold nuclear-rate audit

This independent audit diagnoses the old failed B2 tail. It changes no historical results or production runner. The first completion attempt should correct the upstream artificial-viscosity defect identified by the integration agent and retain the original network. A cold-rate replacement is warranted only if the corrected trajectory still needs one.

## Upstream defect and rate symptoms

The original scalar artificial pressure is activated by a negative radial velocity difference. In spherical geometry that can occur while the cell volume is increasing. An isotropic added pressure dissipates through `−q div(v)/rho`; applying it at positive volume divergence therefore produces cooling. Requiring negative discrete volume divergence for this scalar viscosity is a model correction, independent of any nuclear-rate changes. A tensor viscosity could include shear, but would require consistent tensor work rather than the existing scalar-pressure term.

The old failed endpoint has an outer temperature of 301,365 K. Its fitted F18(n,gamma)F19 and Ne21(n,gamma)Ne22 coefficients are respectively 4.385e68 and 3.125e68 cm³ mol⁻¹ s⁻¹. A conservative single-s-wave unitarity diagnostic at that temperature is about 5.2e11 in the same units. This comparison is not a complete evaluated cross-section bound; it makes the enormous fit divergence unambiguously unphysical.

The scan finds an additional charged-particle defect if a trajectory reaches 100,000 K: O16+C12→He4+Mg24 (network index135) overflows. Its coefficient at10MK is only2.422e-180. Modifying only neutron rates would not prevent that separate overflow.

The initial prescribed inflow already includes cells below10MK (minimum2.20MK), before they heat. It is therefore incorrect to assert that a successful unmodified run necessarily remains above10MK. The relevant checks are the magnitude and physical validity of the actual reaction flows, together with the absence of a numerically destabilizing extrapolation.

## Physical continuation if one is needed

T. Rauscher's primary [rate-fitting guidance](https://nucastro.org/forum/viewtopic.php?id=22) identifies0.01–10GK as the usual fit interval, requires well-behaved extrapolation, and describes constant neutron rates for s-wave capture or rates proportional to temperature for p-waves. This is general guidance, not a verified validity interval for every coefficient in this particular library. The documentation also requires consistency between forward and reverse reactions.

A targeted contingency is a constant `N_A<sigma v>` continuation below10MK for binary, exothermic neutron-induced channels, anchored to their value at10MK. This corresponds to the familiar s-wave1/v cross-section asymptote. The present network has39 such channels. Charged and weak rates retain their actual-temperature evaluation. For the separately diagnosed C12+O16 rate, use a matched Gamow-suppressed continuation if necessary; its tiny10MK anchor also permits an explicit negligible-flow truncation bound. Reverse partners must be adjusted consistently if their low-temperature flows matter. None of these continuations is an evaluated low-energy measurement.

The uncertainty test should span zero cold neutron-induced rates and a generous bounded upper-rate alternative, leaving all hot coefficients unchanged. A completed trajectory showing unchanged dominant isotope yields and energy across that test would support insensitivity to the missing cold microscopic detail. A global temperature floor applied to every rate is unnecessary and obscures which physics was changed.

## Inventory and numerical precision

In the old saved160-checkpoint trace, four cells ever fall below10MK. Their largest positive neutron molar abundance is6.73e-64, whereas the largest positive neutron abundance anywhere in the full trajectory is1.78e-11. The first number is far below the solver's absolute abundance tolerance of1e-11 and cannot be interpreted as physical precision.

Assigning one full absolute tolerance of neutron abundance to every cell that becomes cold gives a mass-weighted neutron budget1.7603e-12. Capturing that diagnostic inventory at10MeV per neutron would change energy by4.71e-10 of prepared input. This is a useful tolerance-scale sensitivity, not a rigorous global ODE error theorem. The cold neutron-producing flows evaluated at10MK are also negligible in the saved states. A corrected completed trace must be audited again rather than inheriting those old values.

The largest mass-weighted negative baryon inventory anywhere in the old trace is1.7735e-12. Its deliberately conservative rest-energy scale is4.42e-8 of prepared input. Finite abundance errors still need explicit reporting; clamping them silently cannot establish better accuracy.

## Positivity and energy handling

The current network uses the same reaction flows for `dY=S f`, `q=Q f`, and neutrino losses, with `A S=0` and `Q=−M S`. Those algebraic identities hold under any common rate continuation. Modifying a reaction coefficient consistently in its species and energy terms preserves them.

Keeping the differentiable mass-action polynomial during Newton iterations is acceptable only while accepted negative inventories remain tiny and bounded. Physically bounded rate coefficients remove the mechanism that turns1e-30-level traces into enormous spurious flows. A positivity projection or clipping/renormalization of accepted states must separately record its isotope correction and change thermal energy by the corresponding rest-mass difference; otherwise it silently changes the conserved energy. No such projection is needed merely to hide negligible negative printed abundances.

`b2_nuclear_rate_audit.py` reproduces the coefficient scan, trace-inventory diagnostics, hashes and four algebraic assertions in `b2-nuclear-rate-audit.json`. The full old trace is `/tmp/protium-small-stages/scaled-hot-B2-64-R30-tight.npz`; the historical endpoint is preserved in the reassessment packet. This audit does not claim a completed hydro calculation or a new gain.
