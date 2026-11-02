# Sources and coefficient provenance

Accessed 2026-10-29. These references supply microphysical methods or their validity limits. The pulse numbers are this checkpoint's calculations, using the archived trajectories identified in the output JSON.

## Compton scattering

**Poutanen, J. (2017), “Rosseland and Flux Mean Opacities for Compton Scattering,” ApJ 835, 119.** [DOI](https://doi.org/10.3847/1538-4357/835/2/119), [author PDF](https://users.utu.fi/jurpou/wp-content/uploads/sites/1275/2021/10/2017_Poutanen_ApJ_835_119_rosseland.pdf), [arXiv](https://arxiv.org/abs/1606.09466).

Read the degeneracy convention, transport formulation, equation 36, Table 1, nondegenerate fit and absorption discussion. `poutanen()` implements the wide fit: 2–300 keV and −10≤η≤7, with published maximum fit error 6.5%. Its comparison with the separate nondegenerate approximation allows both published errors. Total lepton density approximates the pair contribution; this implementation does not reproduce the full frequency redistribution kernel. True absorption is calculated separately.

## Free-free emission and absorption

**Itoh et al. (2000), “Relativistic Thermal Bremsstrahlung Gaunt Factor for the Intracluster Plasma. II. Analytic Fitting Formulae,” ApJS 128, 125.** [DOI](https://doi.org/10.1086/313375), [primary paper](https://arxiv.org/pdf/astro-ph/9906342).

Equations 1 and 7 supply the classical emissivity normalization 1.426e−27 in cgs, with the Gaunt factor explicitly separated. The paper's relativistic fit covers log10(T/K)=6–8.5 and specified frequency/charge ranges; its fit is not extrapolated to the hotter pulse peaks. This checkpoint instead uses constant Gaunt factor 1.2 as a diagnostic. Kirchhoff's law yields the frequency-dependent absorption, independently checked against the integrated emissivity. Degeneracy, collective dielectric effects and relativistic corrections remain explicit limitations.

**NRL Plasma Formulary (2023).** [University-hosted copy](https://homepage.physics.uiowa.edu/~ghowes/teach/phys4731/resources/NRL_Plasma_Formulary_2023.pdf).

The previously retained local copy `/tmp/protium-sources/NRL_Plasma_Formulary_2023.pdf` and its extracted text were consulted for radiation coefficients, plasma parameters and electron–ion thermal equilibration. The equal-temperature exchange coefficient provides an independent normalization check. Coulomb logarithm one is a sensitivity parameter outside the large-logarithm asymptotic regime. An earlier official-host 403 was respected; no blocked endpoint was retried or access control bypassed.

## Conductivity

**Spitzer, L. Jr. & Härm, R. (1953), “Transport Phenomena in a Completely Ionized Gas,” Physical Review 89, 977.** [DOI](https://doi.org/10.1103/PhysRev.89.977), [paper copy](https://ayuba.fr/pdf/spitzer1953.pdf).

Read equations 36–38 and Table III. The classical comparison uses the Lorentz limit, epsilon=0.4 and delta_T=1. Substitution of the tabulated hydrogen values, epsilon=0.419 and delta_T=0.225, recovers the conventional hydrogen conductivity coefficient in a separate check. The calculation does not extend this classical operator into a validated partially degenerate transport theory.

**Cassisi et al. (2007), “Updated Electron-Conduction Opacities: The Impact on Low-Mass Stellar Models,” ApJ 661, 1094.** [DOI](https://doi.org/10.1086/516819), [primary paper](https://arxiv.org/pdf/astro-ph/0703011).

Read the conductivity/collision-frequency relations and discussion of partial degeneracy and mixed compositions. Equations 2 and 8 supply the degenerate electron–ion comparison. The implementation omits electron–electron resistance and takes the larger classical/degenerate conductivity as a sensitivity envelope. Neither that maximum nor Coulomb logarithm one is asserted to give a rigorous bound. The independent algebra check eliminates collision frequency to compare an equivalent expression.

## Coulomb free energy and screening

**Chugunov, A. I. & DeWitt, H. E. (2009), “Nuclear Fusion Reaction Rates for Strongly Coupled Ionic Mixtures,” Physical Review C 80, 014611.** [DOI](https://doi.org/10.1103/PhysRevC.80.014611), [primary paper](https://arxiv.org/pdf/0905.3844).

Read the coupling definitions, mean-field screening formulation and equation 24. The ion–ion energy diagnostic analytically differentiates that one-component-plasma free-energy fit and linearly mixes species. The derivative is checked against finite differences of the independent implementation. Screening uses `pynucastro.screening.chugunov_2009` from installed version 3.1.0; its implementation was inspected. The pressure/energy diagnostic omits electron polarization and deviations from linear mixing. The cold-contact calculation is a failure indicator for the fully stripped model, not a self-consistent chemical EOS.

**Potekhin, A. Y. & Chabrier, G. (2010), “Thermodynamic Functions of Dense Plasmas: Analytic Approximations for Astrophysical Applications,” Contributions to Plasma Physics 50, 82.** [primary paper](https://arxiv.org/pdf/1001.0690).

Consulted the free-energy decomposition and linear-mixing discussion to identify terms omitted by the diagnostic. No complete EOS implementation was imported. A related 2013 A&A paper was read through its indexed primary page; a terminal download returned 403 and was not retried. Ioffe conductivity-code pages timed out or returned a certificate error; no code was obtained, and TLS checks were not disabled.

## Partial ionization and dense contact

**Potekhin, A. Y., Massacrier, G. & Chabrier, G. (2005), “Equation of State for Partially Ionized Carbon at High Temperatures,” Physical Review E 72, 046402.** [primary paper](https://arxiv.org/pdf/physics/0510006).

Read the chemical free-energy treatment, bound-state contributions and pressure-ionization discussion. It establishes the type of thermodynamically consistent calculation needed for cold contact. Its carbon tables were not substituted for the actual evolving C/N/H mixture.

**Massacrier, G., Potekhin, A. Y. & Chabrier, G. (2011), “Equation of State for Partially Ionized Carbon and Oxygen Mixtures at High Temperatures,” Physical Review E 84, 056406.** [primary paper](https://arxiv.org/pdf/1111.0532).

Read the model domain and comparison of ionization/pressure-ionization structure, especially the shell transitions discussed with figures 7 and 8. This supports retaining bound electrons in the contact-state problem; it does not supply ionization fractions for the particular pulse mixture. Hydrogenic binding and orbital-radius values in `contact_screen.py` are explicitly isolated-atom warning scales.

## Reproducibility boundary

Installed numerical dependencies: NumPy 2.2.6, SciPy 1.15.3, pynucastro 3.1.0, Numba and Matplotlib from `/tmp/protium-python`. Full trace files remain in `/tmp/protium-research`; input hashes are recorded in `pulse-transport-audit.json`, checked against the existing trace archive, and linked to the compact output files by `manifest.json`. This directory distributes original code and calculated outputs, without copying the source papers.

## Prompt capture photons: small-pulse follow-up

**NIST X-Ray Mass Attenuation Coefficients, nitrogen.** [Element table](https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z07.html). Accessed 2026-10-29. At 6 and 8 MeV the tabulated total attenuation coefficients are 0.02511 and 0.02209 cm²/g, with mass energy-absorption coefficients 0.01639 and 0.01512 cm²/g. These atomic-material values provide a scale check on the independently evaluated stationary-electron Klein–Nishina expression; they are not a dense-plasma transport table. The small-pulse JSON records the actual nuclear mass-defect Q used in its comparison.
