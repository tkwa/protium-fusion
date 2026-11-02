# Ordinary hydrogen fusion research

This repository contains the calculations, source data and engineering studies behind [How to fuse ordinary hydrogen](https://tkwa.me/blog/ordinary-hydrogen-fusion/). It investigates a proposed, very large pulsed CNO-cycle reactor and earlier alternatives. The calculations are exploratory: the plant performance depends on specified engineering assumptions, and the complete reactor has not been demonstrated.

## Start here

- [Research overview](docs/protium-fusion/README.md) and [sources](docs/protium-fusion/SOURCES.md): the investigation, evidence and historical branches.
- [Current mixture and pulse work](docs/protium-fusion/mixture-design/): the nuclear network, selected stage calculations and diagnostics.
- [Launcher reference](docs/protium-fusion/launcher-reference/README.md): packet delivery, timing and geometry.
- [Receiver and service studies](docs/protium-fusion/component-design/): radiation, materials, collection, pumps, power and supporting components.
- [Warm-buffer plant account](docs/protium-fusion/component-design/catcher-followup/plant-account/README.md): the complete energy and mass balance used for the article's headline output, with [results](docs/protium-fusion/component-design/catcher-followup/plant-account/results.json), [calculation](docs/protium-fusion/component-design/catcher-followup/plant-account/calculate.py), and [independent review](docs/protium-fusion/component-design/catcher-followup/plant-account/INDEPENDENT-REVIEW.md).

The investigation evolved through many designs. Earlier estimates in the research log, magnetic-confinement studies, and historical or fictional-worldbuilding documents are retained as context; check each document's date, inputs and supersession notes before comparing its figures with the article.

## Following the article's numbers

The [warm-buffer interface data](public/blog/ordinary-hydrogen-fusion/catcher-interface.json) summarize the current plant account. The first `case_rows` entry, `warm_compression_buffer_reference`, gives 6.897 × 10²⁰ W net output, 3.339 × 10²² kg and 20.6576 mW/kg. Its `energy_balance` gives the components used in the article's energy table. Its `source_sha256` identifies the underlying plant-account results above.

The [selected plant and isotope data](public/blog/ordinary-hydrogen-fusion/selected-plant.json) describe the earlier allocated reference and its pulse/isotope ledger. Its plant output is about 20.78 mW/kg; the warm-buffer account adds the catcher-interface requirements. These are different cases, so use the warm-buffer case for the article's 20.66 mW/kg and energy balance.

## Repository layout

```text
docs/protium-fusion/
  README.md and SOURCES.md       Research overview and references
  mixture-design/                Pulse mixtures, reaction rates and diagnostics
  launcher-reference/            Launcher geometry and incoming timing
  component-design/              Component calculations and plant integration
  cold-sink-redesign/            Radiator and heat-engine assumptions
  b2-completion/                 Carbon-stage numerical studies
  open-loop-redesign/            Material recovery and earlier plant design
  phase-two/, reassessment-2049/ Earlier investigations and historical scenarios
  *.py, *.json, *.md             Additional analyses, results and research notes
public/blog/ordinary-hydrogen-fusion/
  selected-plant.json            Selected pulse and isotope reference data
  catcher-interface.json         Warm-buffer interface and plant balance
tools/verify_article_balance.py    Check the article case against its source
```

The original `docs/protium-fusion/` path is retained because scripts and result manifests refer to it. Individual study READMEs document inputs, assumptions, commands and checks. Results preserve source hashes and often include precision for numerical reconciliation. Large full simulation traces and some local research environments are not included; compact traces and their manifests are retained where available.

## Running checks and calculations

Many analyses require Python 3.10 and the scientific packages recorded in their individual READMEs. The research overview lists versions used for the original calculations. Some specialized studies also require external data or tools and the full local traces; their READMEs describe these dependencies. Commands may regenerate results, so use a separate checkout when reproducing a calculation.

The article balance can be checked against its source data using the Python standard library:

```sh
python3 tools/verify_article_balance.py
```

One independent accounting check for an earlier branch also runs with the standard library from the repository root:

```sh
python3 docs/protium-fusion/reassessment-2049/cara_independent_audit.py
```

It checks an earlier ion-design branch rather than validating the current warm-buffer reactor. For the current plant, see the warm-buffer account and its independent review above. The numerical checks demonstrate consistency within their stated models; they do not establish the untested physical assumptions.

External nuclear datasets and other source materials are identified in [SOURCES.md](docs/protium-fusion/SOURCES.md) and individual study records. Consult their original licenses and terms before redistributing or reusing them.
