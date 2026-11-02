# Ordinary hydrogen fusion research

This repository contains the calculations, input data and engineering studies behind [How to fuse ordinary hydrogen](https://tkwa.me/blog/ordinary-hydrogen-fusion/). The current reference is a staged CNO-cycle reactor with optically thick inertial pulses. Its performance depends on specified engineering assumptions; no complete reactor has been demonstrated.

The [v0.1 release](https://github.com/tkwa/protium-fusion/releases/tag/v0.1) preserves the research snapshot corresponding to the figures on the website, including the earlier research history. New work on the repository's main branch can change those figures.

## Start here

- [Design and research map](docs/protium-fusion/README.md): the physics models, interfaces and open questions.
- [Pulse mixtures and reactor design](docs/protium-fusion/mixture-design/README.md): stage calculations, isotope balance and performance comparisons.
- [Proposed tamper study](docs/protium-fusion/mixture-design/tamper-study/TAMPER-STUDY.md): candidate liners and an untamped radius sensitivity. The tamper candidates have not been simulated.
- [Launcher](docs/protium-fusion/launcher-reference/README.md) and [driver-physics studies](docs/protium-fusion/launcher-reference/driver-physics/README.md): packet delivery, cold compression, timing and geometry.
- [Receiver and service studies](docs/protium-fusion/component-design/README.md): radiation, fuel recovery, materials and supporting components.
- [Warm-buffer plant account](docs/protium-fusion/component-design/catcher-followup/plant-account/README.md): the complete balance behind the article's headline output, with [results](docs/protium-fusion/component-design/catcher-followup/plant-account/results.json) and [independent review](docs/protium-fusion/component-design/catcher-followup/plant-account/INDEPENDENT-REVIEW.md).

## Website reference

The [warm-buffer interface data](public/blog/ordinary-hydrogen-fusion/catcher-interface.json) contain the article's reference: 6.897 × 10²⁰ W net output, 3.339 × 10²² kg, and 20.6576 mW/kg. The [selected plant and isotope data](public/blog/ordinary-hydrogen-fusion/selected-plant.json) describe an earlier allocation of the same pulse family; use the warm-buffer case for the article's 20.66 mW/kg figure.

## Repository layout

```text
docs/protium-fusion/
  mixture-design/          Pulse mixtures, reaction rates and proposed tamper study
  inputs/pulse-profiles/   Prescribed hot-inflow inputs and associated controls
  transport/pulse/         Radiation and plasma transport diagnostics
  launcher-reference/      Launcher geometry and delivery
  component-design/        Receiver, recovery, assembly and plant integration
  b2-completion/           Corrected spherical pulse solver and supporting results
  cold-sink-redesign/      Radiator and heat-engine assumptions
  open-loop-redesign/      Fuel-cycle and material recovery calculations
  *.py, *.json             Core pulse, nuclear network, EOS and numerical records
public/blog/ordinary-hydrogen-fusion/
  selected-plant.json      Pulse and isotope reference data
  catcher-interface.json   Warm-buffer interface and plant balance
tools/verify_article_balance.py
```

Research outputs retain their original input hashes and any original path labels. Source-specific READMEs describe the numerical environments and data dependencies. Large full simulation traces and some local environments are not included.

## Verification

Run the standard-library check from the repository root:

```sh
python3 tools/verify_article_balance.py
```

It checks the archived plant-account hash and the warm-buffer energy balance. Other calculations generally require the scientific packages recorded in their individual READMEs. Numerical consistency does not establish the untested engineering assumptions.

External data are identified in the [sources](docs/protium-fusion/SOURCES.md) and individual studies. Consult original licenses before redistributing them.
