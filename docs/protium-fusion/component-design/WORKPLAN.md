# ✨ Current component and mixture research plan

Started 2026-10-31. User authorized up to 24 hours. Work in this isolated checkout and preserve all archived research snapshots. Current comparison point is the corrected 64-zone four-stage pulse set plus the cold-sink redesign: 60% thermal conversion, 60% driver, 0.007 fresh carbon/hydrogen mass budget, 5.395 mW/kg conditional station output. This is a comparison point, not an engineered reactor qualification.

## Priorities and ownership

- Root: full-network mixture comparisons, numerical checks, integrated plant comparisons and final research/writer handoff.
- opaque_driver: physical source, delivery and merger architecture; quantify yield versus entropy instead of imposing an arbitrary near-perfect yield requirement.
- ion_dynamics: actual-spectrum material stopping, nuclear loss and recovery options; direct material collection, partial braking and batch extraction remain competing branches.
- pulse_transport: generalized actual-composition material balance, isotope/H recovery, decay and processing, executable scheduling and mixed-isotope tolerances.
- plant_redteam: receiver sector, support, thermal distribution, manufacturing and installed-mass consequences.
- gas_review: packet radiation environment, catcher vapor, opacity and radiation-transport sensitivities, consequential uncertainty review.

Prioritize uncertainties with a subjective chance above 10% of forcing substantial design change, or deserving explanation in the blog. This is a probability-of-consequence threshold, not a requirement to investigate only 10% changes in parameters. Do not invent purity, alignment, rigid-radius or immediate-per-pulse processing requirements.

## Mixture comparison contract

The one-proton A/B1/C13 and two-proton B2 cases are reference compositions. Start with fixed pulse mass, density profile and initial temperature. Compare both unchanged inflow velocity (prepared energy varies) and matched total prepared energy (velocity adjusts). Hold physical simulated duration fixed against each corrected baseline. Use the full 27-isotope network and compression-capped artificial viscosity at 64 zones. Promising cases receive longer-duration, tighter-tolerance and/or refined-grid checks chosen for actual numerical risks.

Report initial composition explicitly as abundances per baryon. Initial-rate formulas motivate candidate ratios but do not predict pulse gain. Re-solve actual isotope/H flows including excess hydrogen losses, chemical processing, decay inventories, pulse cadence and receiver/driver limits before selecting a station point. Fresh C:H denotes externally consumed/lost hydrogen, with no bypass dilution. The 0.0056 fuel / 0.0014 structural split is an allocation to revisit, not a physical constant. Constant-grey-opacity limitations must remain explicit.

## Publication and evidence

The writing thread was queued with the user's instruction to remove historical-version references and present one internally consistent current design. Writer owns narrative assets; research team owns new supporting results and a stable handoff. Distinguish published evidence, analytic screens, simulated outcomes, apparatus requirements and engineered operating points. Use central estimates and subjective 80% credible intervals when useful, with an explicit acceptance criterion and basis. Never present a subjective interval as a measured uncertainty.

## Closed-reservoir correction requested by Thomas

Ordinary matter intercepted by the massive closed receiver remains plant inventory. Incomplete first-pass extraction, vapor circulation, retired substrate and slow trajectories do not imply terminal mass loss. Separate actual escape and nuclear transformation from hold-up, repeated processing and deliberate export. Use eventual-return=1 as the closed-receiver comparison, then add only named terminal mechanisms. The original 0.0056 fuel allocation is a tolerated-loss boundary, not a prediction of carbon use. Gross recycling still carries energy, equipment, scheduling and inventory burdens.
