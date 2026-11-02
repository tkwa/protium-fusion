# Carbon interception channel screen ✨

`stopping_channel_thresholds.py` calculates Q values and laboratory thresholds from the same atomic-mass table used by the pulse energy ledger. It covers illustrative capture and neutron/proton/alpha emission channels for proton, alpha, C12, C13, N14 and N15 projectiles on C12/C13. Elastic channels are labeled separately. The target is at rest; the exact invariant-mass threshold formula is recorded in `stopping-channel-thresholds.json`.

Slowing the outgoing ions before contact with graphite cannot be described as taking them below **every** carbon-transmutation threshold. Exothermic examples include:

| Reaction | Q, MeV |
|---|---:|
| C12(p,gamma)N13 | +1.94349 |
| C13(p,gamma)N14 | +7.55056 |
| C12(alpha,gamma)O16 | +7.16192 |
| C13(alpha,n)O16 | +2.21561 |
| C12(C12,p)Na23 | +2.24088 |
| C12(C12,alpha)Ne20 | +4.61702 |

For comparison, specified endothermic channels do have useful laboratory thresholds: C13(p,n)N13 at 3.23592 MeV; C12(p,alpha)B9 at 8.18929 MeV; C12(alpha,p)N15 at 6.62284 MeV; C12(alpha,n)O15 at 11.34108 MeV; and C12(C12,n)Mg23 at 5.19560 MeV. Closing these channels does not close the exothermic alternatives.

The appropriate claim is that electromagnetic slowing can reduce integrated charged-ion nuclear reaction probabilities through Coulomb suppression. A quantitative loss limit requires cross sections along the actual slowing trajectory, the target isotope mixture, recoil products, the emitted particle spectrum and stopping powers. The threshold calculation supplies no cross section, ion-neutral fraction, magnetic topology, direct-conversion efficiency or equipment mass.

The low-energy reactions are experimentally established. LUNA measured C12/C13 proton capture down to 60 keV centre-of-mass energy, with systematic uncertainties of 7–8%. [Skowronski et al., Physical Review Letters 131, 162701 (2023)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.131.162701).

LUNA also measured C13(alpha,n)O16 at centre-of-mass energies around 230–300 keV. The collaboration-associated data entry distinguishes its 230–305 keV centre-of-mass interval from alpha laboratory energies 305–400 keV. These measurements establish that a few-hundred-keV alpha cutoff is not a zero-reaction threshold. [Ciani et al., Physical Review Letters 127, 152701 (2021)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.127.152701), [Ciani rate data](https://chanureps.chetec-infra.eu/?p=195).

The process .0056 fuel-carbon case leaves .0014 kg of fresh structural carbon per kg fresh hydrogen, about .000471 carbon atoms per net helium atom. Ordinary material reuse can preserve displaced carbon atoms, but nuclear transmutation must be routed and counted by species. In particular, C12(p,gamma)N13 is followed by beta decay to C13, so it changes the isotope and activates the wall without permanently removing elemental carbon. Conversion into another recoverable CNO isotope need not be permanent loss from a more extensive catalyst loop; conversion into discarded heavier ash is. No such additional recycling loop is credited in the present four-stage process.
