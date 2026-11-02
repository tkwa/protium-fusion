# ✨ Finite-cap collection geometry: conditional mobile-particle routing

A continuous front cap and buried side channels must be solved together. This calculation combines them in a two-dimensional anisotropic diffusion domain and weights the result by the current H/He implantation spectra. It provides a quantitative material requirement for the candidate laminate.

## Model

The solid half-cell extends from its symmetry plane to a side collection boundary, with full channel pitch 20 µm and active depth 150 µm. The cavity-facing boundary absorbs an escaping mobile particle. The rear reflects it. The side collects it below the cap; the side above the cap reflects it. The two diffusivities are constant in space, with η=D_parallel/D_normal. The source is uniform laterally and follows each energy bin's NIST-derived mean projected depth. All source depths in the nominal current H/He spectrum lie inside this domain.

The backward first-passage equation is

\[
D_\perp\partial_z^2 u+D_\parallel\partial_y^2u=0,
\]

with u=0 at the cavity, u=1 at the collection boundary and zero normal derivative on reflecting boundaries. A conservative cell-centered finite-volume solve gives the eventual collection probability u. A second solve, −∇·D∇t=1 with zero first-exit time at absorbing boundaries, gives the mean exit-time scale.

This is an **already-mobile population** calculation. Trapping and detrapping determine when a particle joins that population. Reflection before stopping, range straggling, hydrogen chemistry, pressure-driven backflow and the damaged material's transport tensor must enter a complete model separately. The collection boundary is a perfect sink in this solve. Finite downstream pump speed is treated in the [barrier](../barrier-review/) and [gas](../gas/) accounts.

The side collection boundary has zero mathematical width. The physical proposal uses 100 nm-wide slots; particles passing over those slots can enter the sealed channel and stop on its floor. The [geometry](../geometry/) provides a stopping floor and dogleg gas route. Its finite roofs, supports, seals and directly irradiated floor need their own transport qualification.

## Results

Steady source-rate-weighted front-exit fractions of the already-mobile population are:

| Cap thickness | Diffusivity ratio η | H front exit | He front exit |
|---|---:|---:|---:|
| 100 nm | 10,000 | 0.22128% | 0.18190% |
| 100 nm | 100,000 | **0.10307%** | **0.08755%** |
| 10 nm | 10,000 | 0.13430% | 0.10994% |
| 10 nm | 100,000 | 0.02146% | 0.01734% |

These fractions are simulated responses of the specified mathematical geometry. No measured irradiated graphite has been assigned either diffusivity ratio. A thinner cap improves the shallow-population route while changing the cap's damage, permeability and repair requirements. A thicker uninterrupted cap can prevent a shallow particle from reaching a side channel before it returns to the cavity, even when lateral mobility is very large.

The nominal 100 nm/20 µm/η=100,000 row is a useful **qualification target**. Its direct front release already consumes a roughly 0.1% gas-return budget. Leakage from the pumped plenum, solid permeation, seams and collector regeneration must be added rather than included free in the same number. Stage-resolved and packet-resolved probabilities are exported in [fine-candidates.json](fine-candidates.json) for the chronological chamber calculation.

An unchannelled isotropic slab with absorbing front and rear sends a particle starting at depth z to the rear with probability z/L. Increasing its diffusivity shortens the waiting time without changing that branching fraction. Local directional extraction changes this geometry, allowing a substantial collection probability without requiring gas to diffuse through 17.5 cm of graphite.

## Numerical resolution and reproduction

[calculate.py](calculate.py) evaluates 144 pitch/cap/anisotropy combinations using a graded depth mesh and 64 half-pitch cells. It checks the uncapped limit against the independent infinite-strip Fourier solution; the maximum absolute discrepancy is 0.000159. Its broad grid-refinement test is below 0.000797 absolute probability, with the largest discrepancy in the thick 1 µm cap cases.

[fine_candidate.py](fine_candidate.py) resolves the four table rows using depth parameters 320 and 640 and 128 and 256 half-pitch cells. The largest source-weighted difference is 8.46×10⁻⁶ absolute probability. These are deterministic numerical discrepancies, not credible intervals for physical collection. Projected-range straggling and damaged-material properties can dominate them.

Run `calculate.py`, then `fine_candidate.py`, with Python, NumPy and SciPy. Both files resolve their inputs relative to the repository and save source/code hashes. The model and its boundary conditions received an independent review by the helium-evidence agent.
