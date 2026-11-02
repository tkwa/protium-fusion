# ✨ Selected constructed-source plant account

This packet selects the 350 kg/m² receiver point on the near-tied 250–350 kg/m² plateau. It exports **6.9394 × 10²⁰ W** from **3.3389 × 10²² kg**, or **20.7837 mW/kg**, under the stated source, radiation, conversion and equipment assumptions.

- [Point, physical stages and energy/material ledger](point.json)
- [Receiver columns and heat-routing sensitivities](receiver-tradeoff.json)
- [Service-allocation sensitivity](services-sensitivity.json)
- [Manufacturing and economics](economics/README.md)
- [Scientific figure](plant-tradeoffs.png), [PDF](plant-tradeoffs.pdf), [SVG](plant-tradeoffs.svg), [data](plant-tradeoffs.csv)
- [Complete mixture interpretation](../README.md)
- [Exact-family thermal/support assessment](../../component-design/final-receiver.md)
- [Chronological dispatch](../../component-design/final-family-dispatch/README.md)

The main account uses the delivered nitrogen-15 precursor and carbon-13 polymer source, full eventual return of modeled reusable material, and the explicit Geant4 hard-source scenario. All energy missed by the hot receiver is charged to cold equipment. The primary row adds no arbitrary extra missed-energy fraction; separate0.5% and1% sensitivities are preserved. A separate1W/m² thermal bypass remains included in every row.

The350 choice has lower calculated receiver temperature and panel bow than250 while retaining essentially identical specific power. The propagated80% Monte Carlo sampling interval is20.754–20.813mW/kg, with all physical and equipment assumptions fixed. It measures sampling precision, while source construction, spectrum, opacity, mass and lifetime remain design uncertainties.

`build_selection.py` selects and checks the point from the canonical `receiver-tradeoff.json`. The parent `coupled-receiver-tradeoff/integrate.py` rebuilds that table from `component-design/final-family/delivered-B1-C13H2/integration.json` and `component-design/coupled-photon-extension/combined-basis.json`, with columns200,250,350,500,750,1000,2000. `services_sensitivity.py` re-solves the complete mass/power account for each service allocation. `plot.py` generates the standalone figure and CSV. Input and implementation hashes are embedded in each output. All calculations use the recorded Python3.10 research runtime with the existing NumPy/SciPy stack.
