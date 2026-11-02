# Pulse radiation and plasma transport ✨

This packet postprocesses spherical pulse trajectories to estimate Compton scattering, true free–free absorption, thermalization, electron conduction and plasma equation-of-state indicators. It does not rerun hydrodynamics or determine an integrated bremsstrahlung escape fraction.

The sampled energy-releasing interiors showed large local thermalization margins using the diagnostic classical absorption coefficient and temperature-dependent Compton estimate. A 100-fold simultaneous reduction of both coefficients retained a large effective cell depth in the sampled cases. These are model sensitivities, not experimental error bounds. Frequency-dependent radiation transport, relativistic and dense-matter opacity refinements, and the evolving boundary remain important.

[Source references](SOURCES.md), [audit](audit.py), [microphysics](microphysics.py), [validation](validate.py), and the saved JSON files describe the assumptions and calculations. The retained [v0.1 snapshot](https://github.com/tkwa/protium-fusion/tree/v0.1/docs/protium-fusion/phase-two/pulse-transport) records the original paths and full historical report. The original results' source-path labels and hashes describe that snapshot; moving the packet does not recalculate them.
