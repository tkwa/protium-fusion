# Pulsed CNO reactor research ✨

This research develops an ordinary-hydrogen reactor using staged, optically thick inertial pulses. The [mixture design](mixture-design/README.md) describes the current pulse family and isotope flow. The [warm-buffer plant account](component-design/catcher-followup/plant-account/README.md) supplies the website's reference power and mass.

## Models and engineering

- [Mixtures and pulse methods](mixture-design/methods.md): fuel compositions, nuclear reaction network and spherical hydrodynamics.
- [Pulse transport](transport/pulse/README.md): scattering, free–free absorption, thermalization and equation-of-state diagnostics.
- [Launcher reference](launcher-reference/README.md): delivery geometry, timing and source requirements.
- [Component design](component-design/README.md): assembly, radiation transport, recovery, service equipment and the integrated plant.
- [Tamper study](mixture-design/tamper-study/TAMPER-STUDY.md): an untamped radius sensitivity and candidate liners for future simulations.
- [Sources](SOURCES.md): nuclear data and model references. Individual packets document additional sources and their scope.

The pulse calculations start from specified converging flows. A physically demonstrated driver, three-dimensional assembly and instability calculations, coupled spectral radiation transport, and a fully qualified component and recovery system remain open. Nuclear gain, captured heat and net electrical output have separate energy boundaries.

The reference reaction network, EOS and pulse scripts remain at this directory's top level. Prescribed profile inputs are in [inputs/pulse-profiles](inputs/pulse-profiles/); their recorded original paths and source hashes remain part of the results' provenance. Current transport diagnostics are in [transport/pulse](transport/pulse/). The [v0.1 release](https://github.com/tkwa/protium-fusion/releases/tag/v0.1) preserves the website-backed snapshot and earlier research history.
