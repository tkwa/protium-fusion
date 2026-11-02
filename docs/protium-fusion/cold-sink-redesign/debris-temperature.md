# C13 debris: local temperature and coherent kinetic energy ✨

This bounded calculation uses the archived C13 endpoint and expands its shell trajectories to the receiver at **1.1×10⁹ m**. The main result is a **recoil-only ion second-moment scale of about 3,000 K**, or **about 400 K if that energy is shared with fully stripped electrons**. The ordinary radial velocity profile remains overwhelmingly coherent motion. The electron temperature and actual degree of ion–electron sharing require a charge/kinetic calculation; the saved fluid state does not select one unconditional wall temperature.

The conditional 3,000 K estimate has a subjective 80% interval of **1,000–10,000 K** for the **pulse-number-averaged, ion-retained recoil-equivalent second moment**. This judgment assumes monotone expansion, insignificant additional beta energy transferred to ions, and negligible net recoil-energy transfer from ions to electrons or escaping radiation. It is neither a spatial percentile interval nor a claim that all parcels are Maxwellian. The equilibrated 400 K branch lies outside it because it changes the energy-sharing assumption. No unconditional electron-temperature interval is assigned.

## What is being measured

Temperature describes peculiar velocity relative to the **local** mean flow. At one radius, fast shells arrive first and slow shells arrive later. Their difference in mean speed is useful directed energy; treating the entire arrival history as a single Maxwellian would incorrectly count it as heat.

The endpoint is at 70.004 microseconds, with an outer radius of 812 m. Its particle temperature is 27–113 million K. Its residual internal-energy fluence, expressed per receiver area, is:

| Endpoint reservoir | J/m² |
|---|---:|
| Ions | 55,887 |
| Electron EOS contribution | 375,024 |
| Trapped radiation | 1,567,091 |
| Total internal energy | 1,998,002 |
| Smooth coherent kinetic energy | 34,108,878 |

Thus 78.4% of the residual internal energy is radiation. The full internal-energy number cannot be converted directly to an ion temperature. The archived finite-volume kinetic-energy total is 34,155,055 J/m²; this calculation uses the packet-consistent smooth value throughout.

The energy-weighted 1st, 50th and 99th arrival percentiles are 59.10, 94.34 and 268.58 seconds. The corresponding ion-density quantiles are 8.33×10⁹, 2.43×10¹⁰ and 2.99×10¹⁰ m⁻³. Fully stripped electron densities are 3.57×10¹⁰, 1.56×10¹¹ and 1.93×10¹¹ m⁻³. These are dilute plasmas, even though the full pulse has substantial mass. The local volume-expansion times are 19.6, 31.3 and 89.4 seconds.

## Cold limits and the surviving recoil source

A monatomic, fixed-composition gas-only adiabat gives `T ∝ rho^(2/3)`. Applying that diagnostic to the saved shells yields about **3–16 microkelvin** at the receiver. Radiation decoupling, recombination, charge-dependent coupling and nuclear-decay heating invalidate its use as an actual prediction. It establishes how strongly ideal expansion erases the original thermal temperature.

For strictly force-free motion beginning at the archived endpoint, `v − r/t = −r₀/t`. A finite initial radius therefore bounds the local velocity dispersion even when the global velocity distribution is broad. Using the entire 812 m endpoint radius gives an equivalent-temperature bound of roughly **0.005–0.052 K** across the 1–99% energy quantiles. Residual pressure and radiation still accelerate the endpoint flow, so this is a force-free limiting calculation, not a rigorous bound on the actual plume.

C13-stage nuclei before flight are approximately 72.90% N14, 13.55% C13, 10.73% protons and **2.82% O15**. O15 is the relevant delayed recoil source. N13 and O14 are each below 3×10⁻⁸ of nuclei here. The B2 stage has a very different radioactive inventory, so these temperatures must not be copied to every pulse type.

The O15 positron endpoint is 1.732 MeV. The accessible N15 daughter is its ground state: the first excited state is above 5.27 MeV. The calculation uses an allowed beta spectrum with a repulsive point-charge Fermi factor and samples angular-correlation coefficients −1/3, 0 and 1. These are sensitivity cases; an evaluated O15 mixing-ratio calculation has not been substituted. Mean daughter recoil is **84–111 eV per decay**, with the middle case 90.90 eV and a kinematic maximum of 170.69 eV. The mean positron kinetic energy in this approximation is 0.73585 MeV. [CNO endpoint review](https://lss.fnal.gov/archive/2024/pub/fermilab-pub-24-0489-sqms-v.pdf), [evaluated N15 levels](https://nucldata.tunl.duke.edu/nucldata/HTML/A=15/15_04_1991.pdf), [NNDC beta-Q conventions](https://www.nndc.bnl.gov/nudat3/help/betadecay.jsp).

A kick at time τ in initially compact homologous expansion changes both position and velocity. Conditioning on the arrival position and time leaves a peculiar velocity approximately `(τ/t) δv`. Early kicks largely become position–velocity correlation. For decay exponent `x = lambda*t/gamma`, the retained local random-energy weight is

`I₂(x) = integral₀¹ x exp(−xs) s² ds = [2 − exp(−x)(x² + 2x + 2)]/x²`.

The independent quadrature check includes the small-x limit. Applying this kernel to each archived packet gives:

| Quantity | Result |
|---|---:|
| O15 decays before arrival, per total nucleus | 0.01598 |
| Recoil energy before arrival conditioning | 8.339 J/m² |
| Recoil energy after conditioning | **2.106 J/m²** |
| Number-mean ion-only equivalent temperature | **2,838 K** |
| Common temperature if shared with all fully stripped electrons | **392 K** |
| Recoil-only random fraction of bulk plus random KE | **0.00000618%** |
| Recoil-only coherent fraction | **99.99999382%** |

Composition varies by shell. The recoil-only local equivalent temperature has energy-weighted 10th/50th/90th percentiles of **0.00024 / 2,942 / 3,892 K**. Some energetic outer parcels contain almost no O15 and remain much colder in this calculation. A cold core plus minority fast daughters also has a second moment without having a Maxwellian temperature: the filtered recoil per decayed daughter is about 23 eV, far above its average per total nucleus.

## Collisions and charge select the sharing branch

Expansion does not automatically freeze Coulomb collisions. With fixed charge, `n ∝ t⁻³` and `T ∝ t⁻²`, the usual `nu_ii ∝ n/T^(3/2)` stays approximately constant while the expansion time grows. Recombination, neutral fractions and electron heating determine whether that fixed-charge extrapolation remains applicable.

Using the [NRL Plasma Formulary](https://www.nrl.navy.mil/Portals/38/PDF%20Files/NRL_Plasma_Formulary_2023.pdf), with the previously archived 2023 text, the like-ion estimate is

`nu_ii = 4.8×10⁻⁸ Z⁴ n_i lnLambda / (sqrt(A) T_i[eV]^(3/2)) s⁻¹`,

with density in cm⁻³. At trial Ti=3,000 K and lnLambda=10, the median `nu_ii * expansion_time` is 0.835, 67.6 and 2,004 for nitrogen charge 1, 3 and 7. These are core relaxation diagnostics; a minority recoil tail needs velocity-dependent slowing-down. Evaluating the collision rate only at the population-average temperature does not demonstrate rapid tail thermalization.

The same source gives, near equal electron/ion temperatures,

`nu_ei/n_i = nu_ie/n_e = 3.2×10⁻⁹ Z² lnLambda / (A T[eV]^(3/2)) cm³/s`.

A nitrogen proxy with `n_e=Z*n_i` gives these median temperature-difference relaxation times:

| Trial common temperature | Z=1 | Z=3 | Z=7 |
|---|---:|---:|---:|
| 400 K | 57.5 s | 3.19 s | 0.293 s |
| 3,000 K | 1,181 s | 65.6 s | 6.02 s |
| 10,000 K | 7,186 s | 399 s | 36.7 s |

High ionization can therefore share recoil energy with electrons during flight, making the 400 K branch physically relevant. Low charge or beta-heated hot electrons can share much more slowly. These local, near-equal-temperature estimates do not replace integration from a cold core plus a fast recoil tail. The companion [charge-state audit](charge-state.md) deliberately does not infer a charge fraction from Saha equilibrium.

The receiver's roughly 200 K exterior is not a temperature floor for this dilute gas. Radiative absorption/emission and coupling would have to establish that floor.

## Positron energy is the principal unresolved heating branch

The direct receiver-normal column is only 7.44×10⁻⁷ kg/m². A trial positron stopping power of 1–10 MeV cm²/g loses about 10⁻⁷–10⁻⁶ of a typical positron's kinetic energy across one such column. This is a crossing diagnostic, not an angular or time-dependent transport calculation. The PDG gives the relevant positron stopping framework; its neutral-material examples do not establish stopping in this plasma. [PDG passage of particles through matter](https://pdg.lbl.gov/2024/reviews/rpp2024-rev-passage-particles-matter.pdf).

Low direct collisional stopping does not prove escape. Charge separation, ambipolar fields, instabilities and self-generated magnetic fields can alter paths and transfer lepton energy. Primary ejecta-transport work shows that magnetic geometry matters, and recent kinetic simulations demonstrate self-confinement mechanisms. Their quantitative efficiencies cannot be transferred to this C13 plume. [Barnes et al. 2016](https://arxiv.org/abs/1605.07218), [2025 kinetic-instability study](https://academic.oup.com/mnras/article/541/3/2197/8203299).

The following scenarios deposit a fraction f of the O15 **positron kinetic energy instantaneously and locally**, then apply nonrelativistic adiabatic dilution. They are not assigned confinement probabilities. Relativistic particles stored and thermalized late have a different kernel. Neutrinos and annihilation photons are excluded; ionization costs and recombination radiation are not resolved.

| Instantaneous deposited fraction f | Electron-only increment, full-strip denominator | Common particle temperature including recoil | Coherent KE fraction, ideal closed increment ledger |
|---|---:|---:|---:|
| 0 | 0 K | 392 K | 99.99999382% |
| 0.000001 | 3.69 K | 396 K | 99.99999377% |
| 0.001 | 3,686 K | 3,568 K | 99.99994383% |
| 0.01 | 36,857 K | 32,154 K | 99.99949393% |
| 0.1 | 368,567 K | 318,004 K | 99.99499574% |
| 1 | 3,685,672 K | 3,176,507 K | 99.95010180% |

The remaining energy after the time-squared dilution is random energy. In the ideal closed increment ledger, the removed part becomes additional coherent expansion energy. For ballistic recoil, the corresponding change is the correlation of recoil velocity with the revised trajectory. The code also reports fixed-archived-KE diagnostic ratios for comparison. Neither ledger recalculates the changed trajectory or the remaining endpoint radiation reservoir. They prevent the cooling calculation from simply deleting deposited energy.

The total mass-difference Q from all in-flight decays is only 0.741% of the archived coherent KE, even including neutrinos. This gives a loose scale for delayed-nuclear heating, while shocks, imposed braking fields or unresolved conversion of the endpoint internal reservoir remain separate energy sources. A common temperature of 10⁷ K would still leave roughly 99.84% of particle kinetic energy coherent under the full-strip number denominator; it would require substantially more late random energy than the recoil-only branch.

## What this checkpoint settles

The million-kelvin endpoint temperature does not survive free expansion unchanged, and the wide arrival-speed range is predominantly coherent motion. Late daughter recoil supplies a calculable few-hundred-to-few-thousand-kelvin equivalent scale depending on sharing. Beta transport and charge determine whether electrons and ions instead reach thousands or millions of kelvin. A single 10,000 K guess conceals those branches; the present evidence supports conditional estimates and the very high coherent-energy fractions above.

The next discriminator is a coupled charge/electron-energy and velocity-dependent recoil-tail calculation, including self-generated fields or an independently bounded lepton escape fraction. A complete kinetic simulation is outside this bounded checkpoint.

Reproduce with `PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 python3 docs/protium-fusion/cold-sink-redesign/debris_temperature.py` from the repository root; NumPy is required. [The JSON](debris-temperature.json) records input/code hashes and all numeric scenarios. Seventeen algebra, normalization, quadrature and energy-increment checks pass. The recoil kernel and beta angular term also received an independent agent formula review; that review identified the sharing and recoil-tail qualifications retained here.
