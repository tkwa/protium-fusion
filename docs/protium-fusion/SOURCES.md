# Sources for the pulsed CNO design

The references below support individual physical mechanisms or numerical methods. The associated assumptions and evidence limits are recorded in the linked study packets.

## Nuclear rates and equation of state

- [Solar fusion III](https://arxiv.org/abs/2405.06470): hydrogen-burning reactions and evaluations.
- [JINA REACLIB](https://doi.org/10.1088/0067-0049/189/1/240) and [pynucastro](https://pynucastro.github.io/pynucastro/): reaction rates and numerical network tools.
- [TUNL evaluated nuclear half-lives](https://nucldata.tunl.duke.edu/nucldata/HalfLife.shtml): decay times for pulse products.
- [Nitrogen-13 proton-capture constraint](https://arxiv.org/abs/0903.3094): experimental input for one hot-CNO branch.
- [pynucastro stellar EOS](https://pynucastro.github.io/pynucastro/stellar-eos.html): electron, positron, radiation and ion models and their limitations.
- [Nuclear review](mixture-design/nuclear-review.md): evaluations, rate alternatives, uncertainties and the sources used for the current mixture calculations.

## Radiation, transport and hydrodynamics

- [Compton-scattering opacity](https://arxiv.org/abs/1606.09466) and [thermal bremsstrahlung](https://arxiv.org/abs/astro-ph/9906342): radiation physics for the opaque plasma.
- [Pulsed-plasma transport sources](transport/pulse/SOURCES.md): absorption, scattering, thermalization, conductivity and equation-of-state diagnostics.
- [Formation of a spherical plasma liner](https://laro.lanl.gov/esploro/outputs/journalArticle/Formation-of-a-spherical-plasma-liner/9916500393103761) and [multidimensional plasma-liner modeling](https://arxiv.org/abs/2508.09895): experiments and simulations relevant to converging flows.
- [Geant4 electromagnetic physics](https://geant4.web.cern.ch/documentation/dev/plg_html/PhysicsListGuide/electromagnetic/Opt4.html): transport models used in receiver calculations.
- [NIST graphite attenuation](https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z06.html): cold-material interaction coefficients; thick-shield transport requires a separate calculation.

## Engineering and integration

- [Launcher reference](launcher-reference/README.md) and [component design](component-design/README.md) link source-specific calculations for fuel delivery, receiver, recovery, thermal transport and equipment.
- [Proposed tamper study](mixture-design/tamper-study/TAMPER-STUDY.md) cites liner and fuel-mixing studies and records the untested candidate concepts.
- The [v0.1 release](https://github.com/tkwa/protium-fusion/releases/tag/v0.1) preserves the source index and research history accompanying the published website-backed calculation.
