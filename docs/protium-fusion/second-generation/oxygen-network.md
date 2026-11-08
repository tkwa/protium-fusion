# Extended catalyst pool and oxygen return ✨

**Oxygen belongs in the recoverable catalyst account when its actual return route, inventory and operating cost are included.** The selected quartet's oxygen sidestream is a process-boundary choice. This investigation follows CNO-II/III, their hot branches and later returns before assigning a terminal catalyst debit. The results are conditional network calculations on saved pulse histories, with a separately stated rate-update gate.

## Network, reactions and evidence class

[rates.json](rates.json) exports 66 isotopes and 625 linked REACLIB reactions, including reverse reactions, neutron/proton/alpha reactions, light species and weak decays through Si28. [network.py](network.py) has two calculations:

1. A linear tagged-heavy-nucleus calculation in each saved n/p/alpha bath, with actual temperature and density histories. It follows changing isotope, branching, cold decay, proton/alpha/neutron transfers, nuclear charge and reaction Q. The tagged isotope supplies no feedback to its host bath.
2. A nonlinear 66-species abundance integration for pure controls, 1% O16 admixtures and mixed C12/C13/N14/N15/O pulses, holding the source's temperature and density histories fixed. All abundances change together. New mixture heating, pressure, opacity and source preparation still need coupled radiation hydrodynamics.

The histories have 128 radial shells and 400 saved times. Equal-mass stratified quadrature is compared at 32/64 bins for tracers and 8/16 for nonlinear mixtures. The 16-bin controls reproduce the archived hot nuclear Q to within about 1%, separately from source mesh error. Adaptive integration is capped at 1/400 of pulse duration to resolve the time-dependent forcing. Earlier unrestricted-step outputs missed some short burns; they are retained as rejected records.

Current carbon median corrections are applied in both directions, consistently with the selected source. Oxygen and fluorine retain the pinned REACLIB evaluations. The newer LUNA 2025 O17 destruction rates, October 2026 O16 capture measurement and JUNA F19 measurements were read. Complete replacement tables were not obtained in this run: the LUNA zip retrieval timed out and the alternative provider's TLS verification failed. Numerical values below are **frozen-rate comparisons**, with new-rate integration required before a precise loss specification. [Primary-source details](sources.md).

The modeled routes include:

```text
CNO-I:   C12(p,g)N13(beta+)C13(p,g)N14(p,g)O15(beta+)N15(p,a)C12
Hot I:             N13(p,g)O14(beta+)N14
CNO-II:  N15(p,g)O16(p,g)F17(beta+)O17(p,a)N14
CNO-III: O17(p,g)F18(beta+)O18(p,a)N15
Further: O18(p,g)F19(p,a)O16
Hot O:   F17(p,g)Ne18(beta+)F18(p,a)O15
                    Ne18(beta+)F18(beta+)O18
         F18(p,g)Ne19(beta+)F19
Breakout competition: F19(p,g)Ne20;
                      Ne19(p,g)Na20(beta+)Ne20;
                      O15(a,g)Ne19 with subsequent competing capture/decay;
                      Ne18(a,p)Na21 and other alpha/neutron paths.
```

The full coefficient list, rather than this diagram, defines the calculation. F17 and F18 have library half-lives 64.50 and 6590 s; Ne18 1.67 s; O15 122 s. Beta waiting occurs after a tens-to-hundreds-of-microseconds burn. For F18, 99% decay takes roughly 12.2 hours. A same-hot-pulse oxygen return cannot assume the beta-decayed O17 is immediately available.

Several alpha captures remain inside the extended pool: C12(α,γ)O16, C13(α,n)O16, N14(α,γ)F18, N15(α,γ)F19 and N13(α,p)O16. O14(α,p)F17 also has an ordinary return route. Temporary Ne18/Ne19 production can return by decay. Conversely F19(p,γ)Ne20 is a proton-capture escape route. Classification follows subsequent flows and return cost, rather than the projectile or the instantaneous element name.

## Per-pass oxygen conversion

| Tagged input | N14 host: cold result / input | C12 host: cold result / input | N15 host: cold result / input |
| --- | --- | --- | --- |
| O16 | 0.00555 leaves O16 | 0.000722 leaves O16 | 0.00000254 leaves O16 |
| O17 | 0.9687 reaches the quartet | 0.9682 reaches the quartet | 0.00308 reaches the quartet |
| O18 | 0.9758 reaches the quartet | 0.9813 reaches the quartet | 0.2139 reaches the quartet |

The 64-bin O16 advancement values are 0.005534, 0.000724 and 0.000002545 for these hosts; the C13 value is about 0.000620. These are numerical quadrature comparisons, not whole-design uncertainty intervals. The hosts have different actual hot histories: selected peaks span approximately 0.54–1.67 GK, densities and residence times differ, and abundances evolve. Substituting a universal 50% conversion erases these differences.

Forward capture counts also differ from retained conversion. In the return policy below, approximately 23.5 O16(p,γ) events occur per introduced O16 nucleus over repeated passes, while photodisintegration repeatedly sends F17 back to O16. The final retained F17/oxygen branch controls advancement. Counting every forward capture as permanent processing would overstate useful conversion.

## Closed return policy, charge and energy

The illustrative policy uses N14-host exposures for O16 and C12-host exposures for O17, O18 and F19. Cold separation/recovery occurs between passes. Let M contain the probabilities of returning to these four held states. The exact repeated-pass account is

`visits = (I − M)^−1 × [1,0,0,0]`.

Per introduced O16, the 32-bin calculation requires **180.205 O16 passes, 0.96864 O17 passes, 0.06339 O18 passes and 0.00000673 cold F19 passes**. The extended route returns 0.999999873 nuclei to the original quartet; about 1.26×10⁻⁷ enters beyond-F quarantine. Intermediate F19 can react within a hot pass without appearing in the cold F19 queue.

The reaction-integral account per introduced O16 consumes **2.59736 protons**, produces **1.06246 alpha particles**, has negligible net neutron uptake, and releases **9.80275 MeV** of total nuclear mass-defect energy. Nuclear charge changes by −1.48636 elementary-charge units through weak processes. Reaction branching and final isotope mixtures explain the noninteger mean values. Baryon, charge and mass-energy identities are checked independently.

Positive weak Q totals 4.24914 MeV in the mean history. Charging anywhere from zero to all of that energy to neutrinos brackets the nuclear heat without assigning invented beta spectra. A precise economic point should replace the bracket with measured decay spectra. The hot and cold weak terms are separate; neither is added again as fresh energy.

Allocating each repeat pass the corresponding host pulse's preparation per heavy nucleus costs **98.606 MeV per introduced O16**, and 164.34 MeV of driver electricity at 60% preparation efficiency. At 54.57% heat conversion, and crediting returned preparation energy once, the incremental net is **−107.50 to −105.18 MeV/O16** across that neutrino bracket. This is a deliberately expensive dedicated-work comparison. An oxygen-only pulse has not been shown to reproduce the host bath at that allocated work, so it is not a claim of a working pure-oxygen reactor.

A 1.131×10¹¹ kg H2/O equivalent shot contains 3.784×10³⁶ O nuclei. Fully stripped O plus its two protons carry about 6.06×10¹⁸ C of positive charge; the local electron path must close the same charge. At the earlier grain-delivery q/m = 0.002656 C/kg, accelerator charge is about 3.00×10⁸ C. These are distinct circuits and material states. No unneutralized shell charge is credited as confinement or recovered energy.

The integrated radioactive residence is about **757 nucleus-seconds per introduced O16**. At the actual source this is about 3.3×10³¹ radioactive nuclei distributed over the cold-return chain, before chemical residence and batch stocks. Long-lived C14/Na22-like products retain their actual decay/store requirements; the cold absorbing map alone supplies no rapid processing claim.

## Three meanings of “1% oxygen”

The reference circulates 3.514×10⁷ kg/s of returned catalyst and about 1.547×10³³ heavy-nucleus attempts/s. Its half-year CNO inventory is 5.544×10¹⁴ kg. Its actual O16 production is **1158 kg/s**, or 4.360×10²⁸ nuclei/s.

| Definition | O quantity | Consequence |
| --- | --- | --- |
| 1% of catalyst **inventory by mass** | 5.544×10¹² kg held O | At the actual production rate, net oxygen residence is about 152 years. Inventory share does not fix processing throughput. |
| 1% of catalyst **inventory by number** | 6.487×10¹² kg O | Net residence is about 177 years at that production rate. The difference follows oxygen's larger mass. |
| New O feed equal to 1% of catalyst **mass throughput** | 351,360 kg/s | Approximately 303× actual O production; repeated passes handle 6.37×10⁷ kg/s of O16-equivalent material. |
| New O feed equal to 1% of catalyst **number throughput** | 411,108 kg/s | Approximately 355× actual production; repeated passes handle 7.45×10⁷ kg/s equivalent. |
| **Gross repeated oxygen passes** equal to 1% of catalyst mass throughput | 351,360 kg/s of passes | Supports about 1939 kg/s of net new O return in this policy. |

For the 1%-new-feed-by-mass case, the allocated driver draw is 3.48×10²⁰ W and the incremental export penalty about 2.23–2.28×10²⁰ W at the baseline helium rate. At 100 eV of useful sorting work per passed heavy nucleus and 50% efficiency, sorting consumes 7.68×10¹⁶ W. At the inherited equipment coefficient of 10⁸ kg per kg/s processed, fuel-only handling equipment weighs 6.37×10¹⁵ kg. These are specified process/equipment allowances; treating receiver substrate can require much more mass throughput.

For the plant's **actual** oxygen production, the same allocated-work comparison costs only about 7.35–7.51×10¹⁷ W of net output, about 0.11% of baseline export. Sorting costs 2.53×10¹⁴ W, and repeated O-equivalent processing handles about 210,000 kg/s. For the simulated large-N15 plant at 16× helium rate, oxygen production is about **19,963 kg/s** and the same return policy's net penalty is about **1.27–1.29×10¹⁹ W**, roughly 0.13% of its conditional export. The reaction exposure in its unchanged N14/C12 hosts is the basis of that extrapolation.

This small plant fraction does not establish the cheapest resource policy. Allocated dedicated recovery spends approximately **0.85–0.86 PJ of net electricity per kilogram of replacement carbon avoided**. Fresh carbon acquisition would need a comparable shadow cost, or the return policy would need substantially cheaper shared work, before dedicated recovery wins that comparison. The selected fresh-C requirement already fits comfortably inside the available 0.007 C:H allowance. Keeping oxygen as useful recoverable inventory and choosing to process it immediately are therefore separate economic decisions. No carbon price or transport tariff is assumed.

[summary.json](summary.json) combines the separate allocated oxygen penalty with the primary thermal account at fixed main-quartet throughput: export becomes about 9.86×10²¹ W, or 0.2773 W/kg before extra hardware. This is a first-order work comparison. The oxygen ledger gives its added proton consumption and alpha production; full eight-state stationary dispatch and mixed-source feedback remain required before a final integrated specification.

Thus a 1% standing oxygen inventory can be inexpensive, while burning a new oxygen stream equal to 1% of catalyst throughput can take a large fraction of export. Those cases must have separate labels.

## Shared shots, mixed fuel and the 7% heuristic

The preliminary 7%-of-calendar estimate assumed new oxygen throughput at 1% of baseline catalyst, two independent N15-sized shots, 50% conversion in each stage and serialized windows. It is a heuristic comparison. Our modeled O16 first pass is much slower.

With the actual per-pass probabilities but the illustrative host windows and small 1.131×10¹¹ kg shots, the 1%-new-feed mass/number cases require **315% / 368%** of a single calendar. Scaling those equivalent oxygen shots to 1.131×10¹³ kg reduces the arithmetic loads to **3.15% / 3.68%**. At the reactor's actual oxygen production, the same comparisons are **1.04%** for small shots and **0.0104%** for large shots at baseline rate. Multiply continuous loads by throughput multiplier when keeping shot masses fixed. Pure-oxygen source formation, gain, radiation and window shapes remain untested, so these percentages are schedule-equivalence comparisons.

A pure large oxygen shot also needs enough standing oxygen stock. At H/O=3 it contains about 9.52×10¹² kg O, exceeding the baseline 1%-by-mass inventory of 5.54×10¹² kg. That strict inventory allows a pure shot of at most about 6.58×10¹² kg total at this ratio. Trace admixture in a host pulse uses much less oxygen per event. A 16×-throughput plant with the same half-year catalyst residence has about sixteen times the working inventory, so its 1% oxygen-stock case has a different batch constraint.

Piggybacking has a more direct opportunity. At 1% O16 substitution, the existing host attempts advance approximately 7.53×10⁻⁵ O16 nuclei per net helium, about 72% of the original 1.05×10⁻⁴ source. Roughly **1.4% uniform admixture** would balance that source in a linear exposure screen. The larger N15 family produces somewhat more O16 and needs a correspondingly larger return allocation. O17/O18/F19 can use small shared return fractions. This removes a separate-shot calendar debit while displacing some primary fuel and changing the host's entropy, heat release and material preparation.

The 66-species nonlinear tests explicitly compare p:(C12,C13,N14,N15,O16) = 2:(0.25,0.25,0.25,0.24,0.01). At fixed histories, their hot Q values are approximately 7.34, 0.989, 7.34 and 7.16 MeV per initial heavy nucleus in N14, N15, C12 and C13 hosts. The N15-history mixture performs poorly compared with its pure control. More significantly, cold O16 grows from an initial 1% to about **2.09%, 1.00%, 3.37% and 3.79%**, respectively. Captures on mixed nitrogen create oxygen during the same hot shot. Wholesale mixing is therefore not an established cheap-return solution.

A 1% trace O admixture has a smaller response: N14/C12/C13 histories return cold O16 around 0.009944/0.009993/0.009994 per initial heavy nucleus; the N15 case increases it to 0.010084. These screens support testing modest oxygen in hotter return/capture shots and rejecting a universal “all isotopes share one optimal pulse” assumption. Actual mixed-source and pressure feedback remain the next calculation.

## Effective loss and routes beyond the extended pool

The original 27-species quartet account has approximately 4.21×10⁻¹¹ beyond-F nuclei per helium. The new 66-species fixed-history controls also expose a B11-bearing sidestream, particularly from the C13 history. Including these products gives approximately **1.50×10⁻¹⁰ quarantined nuclei/He** before oxygen return and **1.63×10⁻¹⁰/He** with the illustrated oxygen return. These are operational quarantine estimates in the declared finite network and rate set. F16 is an unbound transient and is excluded from terminal loss; the frozen library's low-temperature treatment is not used to claim stable fluorine-16 stock.

The reaction-integrated oxygen return has about **6.94×10⁻⁸ F19(p,γ)Ne20 events per introduced O16**, a substantial contribution to its tiny quarantine fraction. A uniform 10× increase of that rate is a sensitivity of order 6×10⁻⁷ additional losses/O16 to first order, rather than a claimed JUNA correction at all temperatures. Alpha routes contribute separately. The actual JUNA enhancement is temperature dependent and cannot be pasted onto every pulse as a fixed multiplier.

Ne20(γ,α)O16 has a 4.730 MeV threshold; Ne20(p,α)F17 and Ne22(p,α)F19 require 4.130 and 1.673 MeV, respectively. Return is physically possible. The saved C13 baths give one-pass return probabilities of order 10⁻⁸–10⁻⁷, with competing onward captures. [quarantine.py](quarantine.py) therefore follows Ne20/21/22, Na23, Mg24/25/26 and Al26/27 until return or a declared Si28/other quarantine boundary. In this tested policy, only about **5×10⁻⁶** of an introduced neon nucleus returns to the extended CNO pool before that boundary; almost all continues outward. The mean is about 136–165 processing passes, not the inverse of an isolated favorable reaction rate. Its 8-bin tracer/finite-network result is a screen, with Si28 and incomplete later decays explicitly bounded.

This makes ordinary repeated hot-pulse processing an unattractive neon-return policy: it mostly creates heavier quarantine material. Selective gamma irradiation, energetic projectiles, neutron routes or a different confinement/temperature cycle would need their cross sections, competing yields, driver efficiency, separation and return stock calculated before credit. Their low total production rate permits expensive treatment in absolute energy terms, but no economically demonstrated selective return route has been supplied here. Quarantine and small catalyst makeup remain the costed boundary.

Eventual material loss usually dominates this nuclear term. A 10⁻⁴ irreversible loss of baseline circulating CNO is about 3514 kg/s; a 10⁻³ loss is about 35,136 kg/s. The full 0.007 fresh-C/fresh-H allowance corresponds to an approximate mass-loss fraction 5.54×10⁻⁴, before allocating other fresh needs. Held gas, sorting residue, damaged receiver stock and redeposited material remain inventory until a specified route actually removes them.

The useful conclusion is an extended **recoverable catalyst pool with explicit return policies**, small operational nuclear quarantine in these short pulses, and separate material retention accounting. Precise irreversible loss awaits updated oxygen/F19 rates, isotope-boundary extensions, complete mixed-source burns and qualified material recovery.
