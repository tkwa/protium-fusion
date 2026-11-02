# Charge state and reionization energy ✨

The accepted 60% thermal-conversion working case can collect neutral debris. Charge state mainly changes stopping, plasma transport and the optional direct-conversion apparatus. The existing hydrodynamics uses a fully ionized EOS; the stopping packet also uses fully stripped projectiles as a diagnostic. Neither calculation predicts the atomic charge distribution at the 1.1 × 10⁹ m receiver.

This checkpoint uses the archived per-cell isotope composition and velocity distribution, including its finite-flight beta-decay daughters. It calculates ballistic density, charge-conditional recombination timescales and atomic ionization-energy inventories. It does not integrate a collisional-radiative network. Its accompanying [temperature analysis](debris-temperature.md) treats thermal expansion and heating separately.

## What the saved state establishes

At the four-sound-time C13 endpoint, the ion thermal-energy component gives local temperatures 27.4–112.9 million K and mass densities 549–59,525 kg/m³. The analogous ranges are 30.5–124.1 MK for A, 30.5–125.9 MK for B2 and 13.0–39.8 MK for B1. These are hot plasma states. Their fully ionized EOS is a reasonable starting approximation for this screen, with atomic populations still unevolved.

The archived internal energy separately includes ions, electrons and trapped radiation. The radial velocity distribution describes coherent outward motion; it is not a local Maxwellian temperature. Residual pressure can accelerate individual shells beyond the four-time endpoint. Beta-decay recoil, recombination heating, photoelectron heating and electric fields are absent from the later ballistic projection. The projected daughter inventory retains the parent's bulk velocity.

For each original cell we propagate its two boundaries as $r_j(t)=r_{j0}+v_jt$, retain its mass, and calculate

\[
\rho(t)={m\over(4\pi/3)(r_o^3-r_i^3)},\quad
 t_{\rm exp}={V\over \dot V},\quad
 n_{e,\max}={\rho\over m_u}\sum_s Z_sY_s.
\]

The last expression is a fully stripped free-electron upper inventory. At C13 arrival, its energy-weighted 1st/50th/99th percentiles are **3.57 × 10⁴ / 1.56 × 10⁵ / 1.93 × 10⁵ cm⁻³**. The corresponding separate expansion-time percentiles are **19.6 / 31.3 / 89.4 s**. Across all four pulses, the maximum sampled upper electron density is 4.14 × 10⁵ cm⁻³. These conditional densities retain the frozen radial trajectory; pressure, recoil and thermal spreading can modify them.

## Recombination and charge freeze-out

A useful local test is $n_e\alpha(T_e)t_{\rm exp}$. The code evaluates the Verner–Ferland bare-to-hydrogen-like radiative recombination fits for H, He, C, N and O. At **$T_e=197$ K**, the maximum bare-O test over sampled packets is 0.00478 for C13, 0.00920 for A, 0.00644 for B2 and 0.000107 for B1. Values decrease on the warmer trial branches. This makes continued bare-ion radiative capture slow over a local expansion time on those branches. It does not establish which ions survived the earlier, far denser cooling phase. Partially bound ions also require dielectronic, metastable and charge-exchange channels. [Verner & Ferland's data](https://www.pa.uky.edu/~verner/rec.html), [Badnell's level-resolved RR and complementary DR framework](https://arxiv.org/abs/astro-ph/0604144).

For a separate singly ionized diagnostic, the NRL formula is $\alpha_3=8.75\times10^{-27}T_e^{-4.5}$ cm⁶/s, with $T_e$ in eV. Using one electron per nucleus, its largest $n_e^2\alpha_3t_{\rm exp}$ at 197 K is 2.68 × 10⁻⁷. At **10 K**, the C13 radiative bare-O test rises to 0.0245 and its singly ionized three-body test to 0.0169; A's three-body test reaches 0.179. Further extrapolation becomes sensitive to Rydberg populations, radiation, recombination reheating, correlations and the formula's plasma assumptions. At 10 K the hypothetical co-temperature, bare-N coupling parameter reaches 0.41 for C13 and 0.59 for A. The code uses 10 K as a diagnostic, not a predicted temperature. [NRL Plasma Formulary, pp. 54–56](https://homepage.physics.uiowa.edu/~ghowes/teach/phys4731/resources/NRL_Plasma_Formulary_2023.pdf).

**The 197 K radiator does not impose a 197 K electron-temperature floor on a transparent, non-LTE plume.** Consequently the warm-branch rate results cannot settle the cold branch. Saha equilibrium cannot replace the missing time-dependent calculation. Coulomb momentum equilibration is also a different question: for fixed-charge homologous $n\propto t^{-3}$, $T\propto t^{-2}$, the usual weak-coupling collision-frequency factor $nT^{-3/2}$ stays approximately constant. Dilution alone does not imply collisionless flow.

The 1,000 K receiver's blackbody spectrum contains an energy fraction only $6.3\times10^{-52}$ above neutral carbon's 11.26 eV threshold. Its ground-state H/C/N/O photoionization is negligible. The same radiation can affect weakly bound Rydberg levels. Prompt high-energy photons, beta particles and decay photons have different spectra and timing: the grey pulse luminosity does not supply their photoionization rates. A charge network needs $\Gamma_q=\int 4\pi J_\nu\sigma_q(\nu)/(h\nu)\,d\nu$, electron-impact ionization, radiative and three-body recombination, dielectronic channels, and energy feedback.

These constraints support neither “all neutral” nor “still fully stripped” as a calculated wall state. A numerical charge-fraction credible interval would require assumptions about the missing temperature and radiation histories. The thermal working case does not depend on selecting one.

## Atomic threshold inventory and the 50% production case

Ground-state first-ionization thresholds are approximately H 13.598, He 24.587, C 11.260, N 14.534 and O 13.618 eV. Complete stripping totals are approximately 13.598, 79.005, 1,030.1, 1,486.1 and 2,043.8 eV respectively. These are **atomic threshold-energy inventories**, rather than the thermodynamic free energy of a dilute plasma. H/He/C use the accessible [NIST CODATA table](https://physics.nist.gov/cuu/pdf/CODATA_RMP2000.pdf); N/O use NIST values reproduced in the author's [laser-plasma thesis, Table 5.1](https://ft.nephy.chalmers.se/~vojtech/Theses/Doctoral_state_exam_Thesis_Vojtech_Horny.pdf). Their rounding is immaterial to the engineering scale below.

The code sums these thresholds over each velocity packet's actual arrival mixture. H/He/C/N/O exceed 99.99999997% of its nuclei. Other elements receive a separate conservative $14Z^3$ eV binding bound; their station-average contribution is below 8 × 10⁻¹² W/m².

| Pulse | Neutral → singly charged threshold, J/m² | Neutral → fully stripped threshold, J/m² | Full stripping at 50% production efficiency / outgoing kinetic energy |
|---|---:|---:|---:|
| C13 | 80.3 | 7,311 | 0.0428% |
| A | 98.7 | 7,618 | 0.0762% |
| B2 | 92.2 | 7,879 | 0.0534% |
| B1 | 1.47 | 56.7 | 0.0823% |

At the selected complete-cycle cadence, threshold power averages **0.000456 W/m² for one charge per atom** and **0.02753 W/m² for full stripping**. A declared 50% production efficiency charges twice these energies: **0.000911 and 0.05506 W/m²**, respectively. The latter is 0.0694% of outgoing bulk kinetic power. Already ionized material requires less threshold energy.

This establishes that atomic thresholds alone are small. Delivering absorbed ionizing photons or electrons to a thin moving plume, charge-state control, heating, neutral bypass, electron removal, transport losses and hardware are additional burdens. A photon source's wall-plug efficiency does not establish its absorption efficiency. The 50% case is an explicit sensitivity assumption for threshold delivery.

Fully stripped positive-ion current averaged over all shots is **1.448 × 10⁻⁴ A/m²**; one charge per nucleus gives **3.024 × 10⁻⁵ A/m²**. A steady apparatus requires corresponding electron collection and return paths. At potential difference $V$, uncompensated current handling carries power $jV$, irrespective of the small atomic threshold. Electron kinetic/thermal energy and sheath work must remain in the conversion ledger. Recombination energy can be recovered once; it cannot be added as new nuclear output. Atomic binding references also need consistent treatment when refining the fully ionized EOS and atomic-mass nuclear ledger.

## Reproduction

Run `python charge_state.py --output /path/to/charge-state.json` with NumPy. The archived script uses only existing repository JSONs and embedded source coefficients; **21/21 checks pass**, including packet mass closure and an independent NRL hydrogenic-rate comparison. [Results](charge-state.json) retain source hashes, complete temperature-conditional diagnostics and the access record. The NIST ASD CGI request returned 403; the static NIST paper and author-hosted table were used instead. No old calculation files or article text were changed.
