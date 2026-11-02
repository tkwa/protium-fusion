# Residual gas and incoming packets ✨

The saved graphite thermal-evaporation source passes a useful neutral-gas screening test at the existing **1.1×10⁹ m receiver radius**. This supports continuing a design with local contamination control. The calculation does not establish the total sputter, dust, desorption, or charged-plasma load, and does not by itself justify a global membrane.

The source is **0.005722405 kg/s**, or 180,585 kg/year, from [results.json](results.json). The earlier 0.005665 kg/s summary omitted the conservative late-window tail contribution. A 400-fold enhancement gives **2.289 kg/s**, compared with the existing 3,548.56 kg/s structural-carbon allocation. That comparison concerns material replacement; a contamination limit can be substantially tighter.

## Model evidence and source interpretation

[plume.py](plume.py) uses the [NIST-JANAF C](https://janaf.nist.gov/tables/C-003.html), [C₂](https://janaf.nist.gov/tables/C-113.html), and [C₃](https://janaf.nist.gov/tables/C-138.html) thermochemical tables, a unit evaporation coefficient, and the saved receiver temperature histories. The largest prompt temperature is **1,713 K**; the largest slow debris/decay temperature is **1,387 K**. Adding separate temperature maxima into a simultaneous 1,965 K surface would overestimate the evaporation of these separated pulses. The underlying receiver traces and material properties remain reduced-model inputs.

Thermal equilibrium pressure describes the tendency of a particular surface to emit or absorb carbon at its local temperature. It is not the chamber's instantaneous pressure. For example, the saved equilibrium pressure is 5.22×10⁻²⁵ Pa at 1,000 K and reaches 3.66×10⁻⁹ Pa during the short C13 prompt maximum. Material already in flight can persist after the wall cools. A cool graphitic surface has a strong thermodynamic incentive to retain carbon, but the effective capture probability includes reflection, reemission, ports and actual surface chemistry; it requires qualification.

The 400-fold sensitivity has a real experimental motivation: individual graphite/graphene nanoparticles initially sublimed about 400 times faster than bulk graphite, with rates changing as their surfaces evolved. Applying that number to this receiver is an explicit sensitivity, rather than a measured receiver enhancement or a bound on irradiation erosion. [Long et al. (2020)](https://arxiv.org/abs/2002.07835)

## Gas inventory and the column encountered by a packet

For diffuse cosine-law emission from a uniform spherical wall and flux-Maxwell speeds, the mean first return is

\[
\tau_i=\frac{4R}{3\bar v_i},\qquad
\Sigma_i=\frac{3\dot M_i\tau_i}{4\pi R^2}.
\]

At the inherited 1,700 K gas-speed assumption, C, C₂ and C₃ return in **9.81, 13.87 and 16.98 days** on average. Unit capture at the next wall encounter gives a time-averaged inventory of **5,570 kg**: approximately **40 gas molecules/m³**, an equivalent kinetic pressure of **9.35×10⁻¹⁹ Pa**, and a radial mass column of **1.10×10⁻¹⁵ kg/m²**. These are spatially averaged molecular-flow quantities; an individual fresh plume has a nonuniform velocity and density distribution.

The independent calculation also sums a conservative bound over all earlier scheduled pulses. Its radial mass-column bound is **2.573×10⁻¹⁵ kg/m²**, or **1.290×10¹¹ carbon atoms/m²**. The instantaneous single-pulse geometrical column factor peaks at **3.24997**, so the old factor-three illustration was slightly low. The new calculation uses 3.26 and the exact large-time inequality Γ(w)≤8/w⁴ to bound the old-pulse tail. Its numerical peak guard is 0.3%; it is a checked numerical bound within the specified diffuse model, not a theorem covering arbitrary directed plumes.

Let G multiply the saved thermal source and s be effective capture probability per wall encounter. A molecular-flow renewal screen scales the column approximately as G/s. A colder, rethermalized gas at 1,000 K increases residence and column by another factor **1.304** relative to 1,700 K. Neither 1,000 launchers nor subdivision into many packets divides the column seen along an individual incoming path.

| Scenario | Radial mass-column bound, kg/m² | Encounter probability at σ=10⁻¹⁸ m² |
|---|---:|---:|
| Saved source, s=1 | 2.57×10⁻¹⁵ | 1.29×10⁻⁷ |
| G=1, s=0.5 | 5.15×10⁻¹⁵ | 2.58×10⁻⁷ |
| G=400, s=1 | 1.03×10⁻¹² | 5.16×10⁻⁵ |
| G=400, s=0.1 | 1.03×10⁻¹¹ | 5.16×10⁻⁴ |

The encounter calculation is 1−exp(−Nσ), conservatively treating each carbon atom as a separate target. With an allocated **10⁻³ encounter fraction**, even the last stress case passes for σ≤**1.94×10⁻¹⁸ m²**; gas cooling to 1,000 K reduces this cross-section threshold to 1.49×10⁻¹⁸ m². Charge-changing and scattering cross sections depend on projectile species, charge, energy and target state. The 10⁻²⁰–10⁻¹⁸ m² range used here is a screening parameter range, not a fitted beam result. [OPEN-ADAS](https://open.adas.ac.uk/adf01) and the [NIFS carbon collision compilation](https://www.nifs.ac.jp/report/nifs-data091.html) illustrate the channel-resolved data needed once that beam state is specified.

## Dense packets: acceptance in terms of radius and column

A spherical condensed packet of radius a and density ρ has projected mass per area 4ρa/3. If fraction η of the encountered relative kinetic energy remains in the packet,

\[
\Delta u\leq\frac{3\eta\Sigma v^2}{8\rho a},\qquad
\frac{m_{\rm intercepted}}{m_{\rm packet}}=\frac{3\Sigma}{4\rho a}.
\]

The inequality assumes complete interception of the geometrical column; a penetrating projectile can transfer less energy. It does not resolve the stopping layer, surface evaporation, recoils or a charged halo.

For **a=1 mm, ρ=1,000 kg/m³ and v=10⁷ m/s**, complete retention gives **0.096 J/kg** under the saved unit-capture source and **386 J/kg** under the G=400, s=0.1 stress case. At an explicitly allocated **10 kJ/kg** extra-heat budget, the latter case requires a≥**38.6 μm**. The minimum radius scales as ηv²/ρ: at 2×10⁷ m/s it is 154 μm, or 1.54 mm if density is only 100 kg/m³. At 1,000 K gas speed these requirements increase by 30.4%. The 10 kJ/kg allocation must ultimately be compared with the actual cold-packet temperature, phase, allowable ablation and implosion-error budget.

These small bulk energy increments provide a reason to retain dense-packet delivery as a live option. They do not require every exposed atom to remain unperturbed. An already dissociated atomic or plasma beam instead needs the collision/charge-state test above and a field calculation.

## Central values and subjective 80% credible intervals

The following are **conditional analyst judgments**, given the saved thermal trace and an engineered condensing surface. They are not confidence intervals from the hydrodynamic simulation. No interval for total sputter/dust release is inferred.

| Parameter | Central value | Subjective 80% interval | Interpretation |
|---|---:|---:|---|
| Thermal source multiplier G | 1 | 0.1–400 | Broad allowance for accommodation and surface enhancement; current temperature history remains representative |
| Effective capture probability s | 0.5 | 0.1–0.95 | Condensing surface performance, including reemission and access geometry |
| Effective atomic encounter cross section | 10⁻¹⁹ m² | 10⁻²⁰–10⁻¹⁸ m² | Order-of-magnitude collision screen pending species/charge/energy selection |
| Retained fraction η of gas-collision kinetic energy | 0.5 | 0.1–0.9 | Heat partition in a condensed packet; η=1 also tested above |

[gas_assessment.py](gas_assessment.py) propagates explicit independent distributions having these median/10th/90th percentiles. The derived median phase-column bound is **9.40×10⁻¹⁵ kg/m²**, with an 80% interval **5.23×10⁻¹⁶–3.02×10⁻¹² kg/m²**. The corresponding encounter optical depth is **6.33×10⁻⁸**, with an 80% interval **1.14×10⁻⁹–2.33×10⁻⁵**, below the 10⁻³ allocation throughout that interval. For ρ=1,000 kg/m³ and v=10⁷ m/s, the formal minimum radius under the 10 kJ/kg allocation has median **14.5 nm** and an 80% interval **0.57 nm–4.75 μm**. The nanometre end is outside a reliable bulk-packet continuum description and supplies no manufacturing recommendation; the useful result is that micron-to-millimetre packets have substantial thermal-gas margin under these assumptions.

The derived medians differ from merely inserting all input medians because the source uncertainty is asymmetric and capture enters reciprocally. Independence and distribution shape are assumptions. A local hot spot, coating chemistry change, or impact-driven release can invalidate this conditional interval.

## Dust, local gas sources and charged material

The structural-carbon budget cannot establish a clean flight path. For a uniform steady neutral-C source at the specified gas speed and unit capture, the atomic encounter allocation corresponds to about **119 kg/s** at σ=10⁻¹⁸ m²; the 1 mm packet bulk-heat allocation corresponds to about **1,595 kg/s** with complete retention. Both scale down with capture probability. A localized port plume needs its own path-integrated column rather than a whole-sphere average.

Dust is more restrictive because energy arrives in discrete impacts. For the 1 mm packet above, with dust density 2,000 kg/m³ and η=0.5, a grain of radius **0.585 μm** carries enough coupled kinetic energy to consume the entire 10 kJ/kg packet allowance. Requiring fewer than **10⁻³ such impacts per packet** gives a critical-grain mass-column allowance of **5.33×10⁻¹³ kg/m²**. Under an illustrative isotropic 100 m/s dust speed and unit capture, this corresponds to **0.184 kg/s** of that grain population. This is an acceptance scale, not an estimated dust source or a crater calculation. A 1 cm packet raises the critical grain radius and source allowance tenfold. Smaller grains require an aggregate heat calculation; material larger than the critical grain size requires a size-resolved impact-risk calculation.

If any one damaged packet spoils an entire assembly of N packets, a 10⁻³ assembly-failure allocation requires approximately 10⁻³/N per packet. At N=1,000, the illustrative 1 mm dust-source allowance becomes **1.84×10⁻⁴ kg/s**. The launcher count is not N: each of the 1,000 launchers may emit many packets. Whether one affected packet spoils an assembly, can be rejected, or merely changes a small mass fraction remains a design tolerance. No dust-release rate is supplied by the thermal evaporation calculation.

Neutral thermal gas, surface-bound carbon and charged material also require different handling. Neutral molecules follow ballistic trajectories until a collision or ionization event; a magnetic field alone provides no neutral pumping. Ionized gas can form sheaths and collective fields even at low mass density, while electron recirculation and charging alter capture. The neutral-column result supplies a mass source to that calculation, without assuming an ionization fraction or a plasma clearing time.

## Bounded decision and reproduction

The current evidence closes the question **“does the modeled pure thermal graphite plume itself rule out incoming dense packets?”** in favor of continued development under the stated column and heat allocations. It leaves a concrete qualification requirement: measure or model the escaping sputter/dust size distribution and local port columns, then apply the acceptance equations to the selected packet population and acceptable assembly error. A global membrane should earn its added heat, damage and contamination burden through that result; local collectors or protected launch corridors remain candidates.

Run from this directory with Python, NumPy and SciPy:

```sh
python gas_assessment.py
```

The script reads existing `results.json`, records its hash and the `plume.py` hash, and writes [gas-assessment.json](gas-assessment.json). Four checks pass: source-rate reconstruction, guarded scalar maximum, an independent logarithmic-grid maximum check, and the exact tail coefficient of eight. It performs no hydrodynamic rerun, changes no earlier packet, and assumes no manufacturing energy-per-mass coefficient.
