# Pulse mixture comparison figure ✨

[PNG](mixture-gain.png) · [PDF](mixture-gain.pdf) · [SVG](mixture-gain.svg) · [CSV data](mixture-gain.csv)

## Caption

**Pulse gain versus hydrogen-to-catalyst composition for four prescribed spherical hot inflows.** Each stage retains its pulse mass and enclosed-mass density/temperature profiles. The pale curves show the completed 64-zone frozen-rate sweep and local refinements at fixed total prepared energy, with radial velocity adjusted as composition changes. Dashed orange curves retain the reference radial velocity profile and therefore change total prepared energy. Blue points show the tested 128-zone cases using median updates to the C12(p,γ)N13 and C13(p,γ)N14 forward/reverse rate pairs; their displacement from the pale curves combines changes in spatial resolution and these rates. Physical end times are fixed within each stage, at the values printed on the panels. All runs use grey opacity 0.001 m²/kg, compression-cap artificial viscosity, and the full 27-isotope, 194-rate network. Lines join tested cases as a visual guide. The four dotted reference recipes have H/catalyst = 1, 1, 2 and 1 for A, B1, B2 and C13, respectively.

Gain is the solver's `deposited_nuclear_over_initial_energy`: `(nuclear_energy_J − weak_energy_escaped_upper_J) / initial_energy_J`. The denominator includes initial internal and bulk kinetic energy. This measures nuclear deposition during the hot pulse; subsequent photon escape and whole-plant recovery, preparation efficiency and auxiliary loads have separate accounts. The plot contains completed nominal-duration cases. Source-formation profiles, opacity/temperature/energy sensitivities, longer-duration cases and rate uncertainty bands have their own comparisons.

## Data and reuse

The figure has **65 plotted rows from 61 unique completed records**. Four fixed-flow reference rows reuse the corresponding equal-energy baseline: at the reference composition the velocity multiplier is one and the initial states coincide. `shared_reference=True` and the original recorded constraint identify these rows explicitly in the CSV.

| Stage | Frozen-rate64 equal energy | Frozen-rate64 fixed flow, including shared reference | Median-carbon128 equal energy |
|---|---:|---:|---:|
| A | 11 | 4 | 2 |
| B1 | 9 | 4 | 3 |
| B2 | 9 | 4 | 2 |
| C13 | 10 | 4 | 3 |

The equal-energy64 data come from the `control`, `grid` and `refine` run directories. The 128-zone data come exclusively from `mixture-rate-runs/nominal128`. [Provenance](mixture-gain-provenance.json) records absolute original paths, hashes, profile/solver/runner hashes, nuclear-update metadata, per-stage durations and filtering conditions. The compact CSV is sufficient to replot without heavy simulation traces. No error bars or probability interval are inferred from these distinct model families.

## Reproduction

`mixture_gain.py` produces the PNG, PDF and SVG from the compact data. It requires Python and Matplotlib. The pinned environment used for these exports is recorded in [mixture-gain-exports.json](mixture-gain-exports.json).

```sh
PYTHONPATH=/tmp/protium-python PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/protium-mpl \
  /tmp/protium-python-interpreters/cpython-3.10.20-linux-x86_64-gnu/bin/python3.10 \
  docs/protium-fusion/mixture-design/plots/mixture_gain.py
```

Add `--rebuild-data` to re-read the original simulation records. `--run-root` and `--rate-root` allow equivalent source directories elsewhere. Rebuilding checks completion, physical duration, opacity, composition, resolution, nuclear-rate family, mass, exact gain definition and the declared energy/velocity constraints. Duplicate stage/series/composition records are rejected. The current case counts are asserted to prevent an accidental silent change of plotted family.

The PNG was visually reviewed at export resolution; labels, legend, plotted points and scope notes are separated. Export hashes are recorded. The original simulation files and archived research remain unchanged.
