# ✨ Current catcher incidence, stopping and elastic-transfer screen

The exact delivered-B1/C13H2 family reproduces the requested H/He fluences. An intact **350 kg/m² graphite backing at R=1.1×10⁹ m** has a large stopping margin: the largest modeled endpoint-speed CSDA column is **0.1793 kg/m²**, about **1952 times smaller** than the backing column. Halving electronic stopping reduces that margin to977. The material-handling question is therefore concentrated in the shallow implanted layer and its connected collection surfaces.

These are calculated collisionless arrival spectra and continuous-slowing estimates. Implantation depth, prompt backscatter, trapping, subsequent gas release, collection and permanent escape are separate quantities. This packet supplies the first two spectrum/depth inputs; the [hydrogen](../hydrogen-evidence/README.md), [helium](../helium-evidence/README.md) and geometry packets address the subsequent steps.

## Exact family and load

The source is [delivered-B1-C13H2-stage-map.json](../../final-family/delivered-B1-C13H2-stage-map.json); the cadence is the fixed [350 kg/m² dispatch input](../../final-family-dispatch/sigma350/input.json). All four profiles have128 cells. B1 uses the actual delivered κprep=.001 endpoint. C13 uses the H2 source profile. No earlier H4, nominal hot256 or alternate B1 profile is substituted.

Mass in this table follows the hydro convention A mᵤ. Neutral atomic masses in the gas/plant accounting give slightly different kilograms for the same number of nuclei; protium alone differs by0.7825%.

|Stage|H, µg/m²|He, µg/m²|H+He, µg/m²|All material, µg/m²|Mean period, s|
|---|---:|---:|---:|---:|---:|
|A|96.13243|0.00001144|96.13244|743.80165|810809|
|B1|0.525969|1.395869|1.921837|7.438017|7698.20|
|B2|20.74323|0.00000000559|20.74323|743.80165|988256|
|C13|30.33186|0.00007255|30.33194|743.80165|11207481|

The user's rounded96.1/1.92/20.7/30.3 µg/m² values and0.744 mg/m² large-pulse mass agree. Here H=p and He=⁴He; the27-species source network does not separately evolve deuterium, tritium or³He. Endpoint, finite-flight and fully decayed isotope inventories are saved separately. Weak decay during flight is included, with daughter nuclei inheriting their parent's bulk velocity. The H/He values are effectively unchanged by that decay; β particles and recoil heating need their separate transport treatment.

|Stage|H nuclei/m²|He nuclei/m²|H kinetic energy, MJ/m²|He kinetic energy, MJ/m²|
|---|---:|---:|---:|---:|
|A|5.78923×10¹⁹|1.72275×10¹²|3.64382|0.0000002394|
|B1|3.16746×10¹⁷|2.10153×10¹⁷|0.0136454|0.0233496|
|B2|1.24919×10¹⁹|8.41043×10⁸|1.50063|8.456×10⁻¹¹|
|C13|1.82663×10¹⁹|1.09228×10¹³|3.02217|0.0000015007|

## Range calculation and primary data

[The NIST ICRU90 update](https://www.nist.gov/system/files/documents/2017/04/26/newstar.pdf), Tables A.8/A.11, supplies graphite proton/alpha electronic and elastic stopping, ranges and detour factors. The PDF, extracted text, exact tables and hashes are retained in[sources/](sources/). The proton electronic stopping incorporates measured low-energy graphite data. C/N/O retain the separately archived [IAEA measured electronic stopping](https://www-nds.iaea.org/stopping/) and screened elastic model; other trace ions retain the earlier analytic treatment.

For mass stopping powers Sₑ and Sₙ, the code integrates

`RΣ(E) = ∫₀ᴱ dE′/[Sₑ(E′)+Sₙ(E′)]`,

`Eelastic(E) = ∫₀ᴱ Sₙ(E′)/[Sₑ(E′)+Sₙ(E′)] dE′`.

Units are converted explicitly:1 mg/cm²=0.01 kg/m². Dividing column by the assumed2000 kg/m³ gives path length. NIST distinguishes CSDA path length from the mean projected penetration along the incident direction; its detour factor is their ratio. [NIST definitions](https://physics.nist.gov/PhysRefData/Star/Text/appendix.html).

An internal consistency diagnostic is preserved: integrating the tabulated proton stopping powers gives ranges above the printed range column by13.5% near30 keV,6.68% at0.3 MeV,1.21% at1 MeV and0.244% at3 MeV. We do not identify this as a proven NIST error. The nominal model consistently uses integrated stopping above1 keV, anchored to the published1 keV range; both range choices are exported. Below1 keV, explicit matched electronic/range completions and screened elastic stopping are used. The unresolved elastic energy below that threshold is bounded by1 keV per incident H/He nucleus in the saved results.

|Stage|H CSDA median /99.9%, µm|H mean-projected proxy median /99.9%, µm|He CSDA median /99.9%, µm|CNO CSDA median /99.9%, µm|
|---|---:|---:|---:|---:|
|A|2.880 /32.989|2.783 /32.814|2.296 /10.178|3.382 /17.104|
|B1|2.104 /9.813|2.008 /9.696|2.016 /5.348|2.587 /7.991|
|B2|9.159 /29.501|9.044 /29.335|1.829 /5.449|3.395 /16.922|
|C13|11.216 /87.667|11.095 /87.380|2.069 /17.686|5.400 /42.285|

Each percentile is over incident nucleus number. The H/He projected proxy is the integrated CSDA path times the evaluated detour factor at each fixed packet energy. It replaces that energy's stochastic depth distribution with its mean. It is consequently unsuitable as an exact straggling tail or a measured retention profile. CNO depth uses straight CSDA as a screening approximation. The current C13 H2 tail reaches~88 µm; the earlier~40 µm reference does not cover it.

The normal-incidence intact-slab model transmits zero primary charged particles through17.5 cm. This numerical result does not assign a statistical zero to straggling, cracks, pores, joints or apertures, and does not bound neutral secondary transport. Surface roughness, entry charge equilibration, temperature, crystallographic channeling and evolving microstructure remain relevant to the implantation profile.

For an illustrative isotropic diffusion slab with perfectly absorbing front/rear boundaries, the rear-exit probability is `min(x/L,1)`. The B1 helium mean-depth distribution gives8.26%,4.13%,2.06% for L25/50/100 µm; C13 hydrogen gives52.09%,33.64%,18.42%. These ideal geometry calculations demonstrate the importance of local collection paths. They do not predict the actual transport tensor or collection efficiency.

## Prompt reflection and material transfers

No retrieved primary study establishes an absolute prompt MeV H/He backscatter fraction for the proposed hot, mixed-ion, irradiated graphite. [Reichart et al.](https://arxiv.org/abs/cond-mat/0512240) measured a localized implanted-H distribution following2.25 MeV proton irradiation of HOPG. Their missing local microspot inventory supports lateral migration as a possibility; it is not a measured escape yield. The separate evidence packets retain conditions and access limitations for release experiments.

A prompt fraction f can be applied directly to the saved nucleus/energy packets in a gas-transport scenario; it should remain an input until an appropriate scattering/retention calculation or measurement exists. Elastic recoil energy alone does not determine sputter yield, ejected size distribution, chemically released H or open-hole production. Reflected particles and sputtered material inside the closed sphere remain recoverable stock; terminal escape requires a connected outward path or a separately identified transformation outside the useful inventory.

## Elastic-transfer and refurbishment exposure

The NRT-style screening count is `0.8 Eelastic/(2 Ed)`, with Ed25 or40 eV. It uses total recoil transfer before recoil electronic partition, cascade recombination or hot annealing. It is an exposure proxy rather than surviving defects or a requirement to replace that number of atoms. Electronic-excitation damage remains separate. The original framework is [Norgett, Robinson and Torrens1975](https://doi.org/10.1016/0029-5493(75)90035-7).

|Stage|Elastic transfer, J/m²/pulse|NRT proxy events/m²/pulse, Ed25 eV|
|---|---:|---:|
|A|342772|3.42307×10²²|
|B1|2769.94|2.76618×10²⁰|
|B2|361554|3.61063×10²²|
|C13|396549|3.96010×10²²|

At the actual cadence, elastic transfer totals**37.358 MJ/m²/year**, giving**3.73072×10²⁴ proxy events/m²/year** at25 eV and2.33170×10²⁴ at40 eV. Averaging the proxy over the first50 µm gives0.7441 dpa/year at25 eV; the most exposed1–2 µm bin gives8.308 dpa/year. These distinct averages must not be interchanged.

Electronic-stopping multipliers0.5/2 give65.007/20.901 MJ/m²/year and first50 µm averages1.293/0.4163 dpa/year. They are sensitivity scenarios, not an80% credible interval. Depth bins,50/100 µm skins and per-species recoil energies are available for selective restoration accounting without assuming wholesale layer atomization.

## Thin cover and incoming-gas comparators

[cover-and-gas.json](cover-and-gas.json) screens100 nm and1 µm intact covers. For B1, carbon stops0.107%/20.44% of H and0.0897%/13.25% of He by number, respectively. All-species primary energy removal is4.40%/40.41%. Thus an ion-transmitting cap can greatly reduce direct implantation into its own material, while still experiencing substantial ionization and recoil exposure. Removal from primary kinetic energy is not automatically local retained heat: secondary electrons and photons have their own transport.

At100 nm the carbon NRT25 proxy per pulse is0.0295,0.000319,0.0316 and0.0197 for A/B1/B2/C13. Neither track overlap nor NRT count specifies a gas-permeable hole. A proposed sealed plenum still needs cap continuity, thermal response, support, permeability, joint sealing and repair performance. Aluminium rows use older analytic H/He stopping and are a lower-confidence comparator; they do not supply an independently evaluated thin-cover acceptance. Their assumed2300 kg/m³ density is inherited from the earlier collection-layer screen, not a demonstrated solid foil state at1000 K.

For gas hitting incoming packets at2.4/4.6 Mm/s, H kinetic energies are30.1/110.5 keV and He119.5/438.9 keV. Graphite CSDA depths are0.424/1.011 µm for H and0.630/1.500 µm for He. Heating is therefore concentrated near the front of the proposed2.7–27 µm carbon capsule skins. Full energy retention is a separate upper deposition assumption. The0.3 Mm/s H row lies below the primary table's1 keV lower limit and uses the stated continuation.

## Reproduction and interfaces

Run `spectra.py`, `parse_nist.py`, `calculate.py`, `auxiliary.py`, then `verify.py` with the project's pinned Python/science runtime. The primary PDF/text are already archived. No hydro rerun is involved. `surface-spectra.json` contains32 quadrature packets per cell with incidence times and exact isotope counts; the integrated report uses64 and records their agreement. `surface-stopping.json` matches the32-point packet order exactly for H/He release and geometry calculations. `stopping-results.json` contains all charged species, depth/column quantiles, bins, sensitivity rows and current-cadence accumulation.

[validation.json](validation.json) passes297 source-hash, mass, energy, unit, cadence, packet-order, depth-integration and NRT identities. This verifies the calculation's provenance and bookkeeping. It does not select a hot-wall retention fraction, damage lifetime or statistical interval for the future collection apparatus. Earlier model and article files remain unchanged.

`charged_flux.py` additionally saves[surface-charged-flux.json](surface-charged-flux.json), the summed all-charged-ion packet counts and arrival flux. Its64-point sampled peaks are8.05585×10¹⁷,6.61324×10¹⁵,6.37869×10¹⁷ and5.75863×10¹⁷ nuclei/m²/s for A/B1/B2/C13;32→64 changes are below7.1×10⁻⁵. These are inputs for a conditional damage/healing model, with no inferred hole yield. An independent read-only review by the sector/geometry agent confirmed the column/depth, elastic-energy and NRT unit conversions.


## Secondary-neutron scope check

[secondary-screen.json](secondary-screen.json) checks the specific¹³C(α,n)¹⁶O channel along the actual He stopping paths. The [2023 direct measurements](https://doi.org/10.1103/PhysRevC.108.L061601) at2.9–8 MeV lie below1 barn; the dominant B1 alpha distribution is mostly below2 MeV, so a precise yield needs the low-energy cross-section curve. Instead, constant1/10 barn values are explicit generous conditional ceilings, not fitted rates or rigorously established bounds at every energy. At10 barn and1.07%¹³C, the B1 probability is2.29×10⁻⁶ per incident He and total reaction kinetic energy averages3.11×10⁻⁵ W/m². Assigning the entire target to¹³C gives0.00268 W/m². Every reaction's full incident energy plusQ2.215 MeV is assigned to the escape allowance, although the neutron receives only part. This channel is unlikely to change plant energetics. The calculation preserves possible neutron escape; other channels, secondary transport and local dose remain separate. The daughter oxygen remains material inventory rather than automatically becoming a terminal catalyst loss.
