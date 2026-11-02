# External-carbon process balance ✨

The open-loop design supplies carbon from outside the station and recycles its circulating C12, C13, N14 and N15. It removes the local triple-alpha carbon factory. The chosen supply constraint is **0.007 kg of fresh carbon per kg of fresh protium**, including hydrogen actually lost by the process. This is a design budget supplied by the user, not a measurement of universal composition.

`process.py` exposes `load_stages(root, stage_paths)`, `case(stages, C13fraction=.0107, carbon_loss=..., hydrogen_loss=None, eta_driver=.6, eta_heat=.4, capture=.99)` and `solve_supply_target(stages, target=.007, C13fraction=.0107, hydrogen_loss=None)`. The default hydrogen loss equals the CNO loss. Explicit stage paths map `c12`, `c13`, `n14`, `n15` to saved pulse summaries. There is no factory term.

## Provenance and present status

`process.json` uses four completed corrected 64-cell stages, each integrated for four sound times with the `compression_cap` viscosity. The explicit source map is `corrected-stage-map.json`; all source summaries are archived in `b2-completion/`. `process-legacy-diagnostic.json` preserves the historical comparison, including the incomplete B2 endpoint and the incorrect expanding-cell scalar-viscosity activation. Source acceptance checks the viscosity mode, effective termination, solver success and elapsed sound times. This is a completed **provisional 64-cell family**, not a spatially converged solution; the subsequent B2 energy-matched grid comparison is recorded in `process-scenarios.json` and `process-scenarios.md` without replacing these source results. The algebra has also been checked by an independently formulated six-variable material solve: all 21 operating points passed 84 checks (`b2-completion/open-process-comparison-corrected.json`). No apparatus mass is inferred here.

## Material equations

Let `M[i,j]` be the number of cold catalyst nuclei of species `i` leaving one attempt of stage `j`. Stage order is C12, C13, N14, N15. Each attempt consumes one nucleus of its named isotope; C12 consumes two protons and the other stages one proton. Cold decay maps the short-lived products to their daughter isotopes, with the inherited conservative neutrino-energy treatment.

For one incoming external carbon atom, `f=(1-x13,x13,0,0)`, CNO return `rC`, and the attempt vector `a`,

```
(I - rC M) a = f.
```

The fresh hydrogen requirement is `sum(a*p_input) - rH*sum(a*p_output)`. Thus it includes both nuclear consumption and actual unrecovered hydrogen. There is no bypass hydrogen added to dilute the C:H ratio. CNO return applies to all four catalyst species after every stage, and hydrogen return applies to cold protons after every stage. Nuclear transmutation into other species is already present in `M`; perfect process recovery cannot undo that leakage.

Every nonreturned output enters the ash ledger. He4 is product; other output species, unrecovered protons and discarded catalysts remain explicit. The code checks baryon conservation and independently reconstructs the summed nuclear Q from external input mass minus product/ash mass. This avoids treating replacement carbon as an energy-free catalyst if part of it is consumed in net transmutation.

The pulse simulation uses baryonic mass `A*m_u` to normalize attempts per pulse. External supply masses use atomic masses: C12 exactly 12 u, C13 13.00335483507 u, and protium 1.00782503223 u. A C:H **mass** ratio of .007 corresponds to an atom ratio about .0005874 for the representative carbon mixture. [CIAAW carbon](https://www.ciaaw.org/carbon.htm), [CIAAW hydrogen](https://www.ciaaw.org/hydrogen.htm), [NIST carbon](https://pml.nist.gov/cgi-bin/Compositions/stand_alone.pl?ascii=ascii&ele=C), [NIST hydrogen](https://pml.nist.gov/cgi-bin/Compositions/stand_alone.pl?ascii=ascii&ele=H).

## Energy and process boundary

`H_nonnu_MeV_per_He` is the total mass-defect nuclear energy minus neutrino energy. It includes photons which escaped the plasma; the closed receiving structure can collect those photons. `U_prepared_MeV_per_He` is the injected thermal-plus-kinetic energy of all stage attempts. The net electric energy before other auxiliaries is

```
W = eta_heat * energy_capture * (H_nonnu + U_prepared)
    - U_prepared / eta_driver.
```

Prepared energy is credited to collected heat once, then its electrical preparation cost is subtracted once. Bulk kinetic energy, hot internal energy, escaped prompt photons and delayed decay heat are reported separately and checked to sum to `H_nonnu + U_prepared`. The diagnostic default is 60% wall-to-prepared efficiency, 40% thermal conversion and 99% energy capture. The process module does not establish these equipment efficiencies.

Energy capture and material return are independent parameters: a stopped ion can deposit its heat and then leave as discarded ash. If material actually escapes the first receiver, its carried energy must also be removed by the plant model. Separation, isotope assay, recycling chemistry, launch hardware, source conditioning, radiators, confinement, supplies, replacement and transport require their own power and mass budgets. A chosen return fraction is an equipment performance requirement, not an experimental result.

The inherited summary-to-ledger conversion clips minute negative numerical abundances and renormalizes total abundance; it rejects total negative mass fraction above 1e-7. The actual correction must be inspected in the stage data and compared with the numerical error budget. It does not justify clipping or injecting energy inside the hydrodynamic integration.

## Corrected process result

At the .007 supply target, coupled CNO/H return is **.9995280479 (99.952805%)**. Fresh H is 4.005476 atoms per He: 4.004599 consumed by nuclear reactions and .0008775 unrecovered. Fresh carbon is .002352706 atoms per He. The completed corrected pulse process gives H=24.733287 MeV/He, U=2.295113 MeV/He and W=6.878058 MeV/He at the illustrative conversion efficiencies. Nuclear Q is 26.762439 MeV/He, including the small net external-carbon transmutation contribution; neutrinos take 2.029152 MeV/He.

At one exawatt (10^18 W) of preauxiliary electrical output, this is 42.580 kg/s fresh carbon and 6082.921 kg/s fresh hydrogen. Gross pulse-fuel circulation is about 107,551 kg/s and returned fuel about 101,504 kg/s. These are scale factors, not a selected plant rating. The prepared/hot/prompt/delayed energy partition is 19.914238 MeV/He bulk kinetic, 1.088454 hot internal, .345531 prompt escaped photons and 5.680177 delayed decay heat, summing to 27.028400 MeV/He including recovered prepared energy.

The returned streams require a real inventory and beta-decay holdup before the assumed cold daughter routing. This steady balance does not assume instantaneous physical decay or zero startup catalyst stock. Extracting and delivering fresh carbon, constructing the receiver and filling the circulating stock are separate burdens.

The corrected recovery scan shows why an external carbon source does not eliminate recycling under this supply constraint:

| CNO and H return | Fresh C:H mass ratio |
|---:|---:|
| 0% | 5.99 |
| 90% | 1.2651 |
| 99% | .13959 |
| 99.9% | .014420 |
| 99.99% | .0017680 |
| 99.999% | .00050145 |
| 100% process return | .00036071 |

The zero-return case produces negligible helium and mainly converts fresh C/H into heavier ash; its very large energy per He is a normalization artifact. The supply constraint excludes it. At perfect process recovery, the remaining carbon demand is nuclear leakage.

The `H_recovery_sensitivity` scan separately tests lower hydrogen recovery while solving the same mass-ratio boundary. Lower H return permits larger C losses mathematically because more fresh H is genuinely consumed by replacement; it increases absolute supply burdens. It should not be called a reactor improvement without counting that extra supply and the changed separation hardware.

The representative fresh carbon is 1.07% C13 by atom. The CIAAW normal terrestrial interval is .0096–.0116; the code scans that interval and pure C12/pure C13 endpoints. Terrestrial isotope composition is a feed assumption, not a universal interstellar value. Under the corrected process, the terrestrial interval changes preauxiliary net energy by less than one part per million and changes the required CNO return by less than 9e-8; enriched C13 provides no compelling process advantage. A decision to use isotopic enrichment would need to pay for its equipment and energy.

## Reproduction

Run from the worktree with the preinstalled scientific dependencies:

```
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 NUMBA_DISABLE_JIT=1 \
MPLCONFIGDIR=/var/tmp/tkwa-chatgpt/protium-mpl \
PYTHONPATH=/tmp/protium-python:docs/protium-fusion \
python3 docs/protium-fusion/open-loop-redesign/process.py \
  --stage-map <explicit-corrected-stage-map.json> \
  --output docs/protium-fusion/open-loop-redesign/process.json
```

Omitting `--stage-map` deliberately selects and labels the invalid historical diagnostic. Sources and hashes are recorded in every report so future pulse revisions cannot silently change the calculation.
