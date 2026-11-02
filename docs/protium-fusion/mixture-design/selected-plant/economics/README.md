# Selected-plant manufacturing and economics ✨

Source: [point.json](/home/tkwa/.codex/worktrees/protium-fusion/website/docs/protium-fusion/mixture-design/selected-plant/point.json), selection `/ (document root)`. Source scope: See the exact selected-input record.

This explicitly selected point exports **6.93942e+20 W** from **3.33888e+22 kg**, or **20.78 mW/kg**. It requires **48.11 installed kg/W**. The largest named mass is `processing_structure_thermal_service` at **45.5%** of the total.

## Manufacturing energy

At the declared availability 1, every **1 MJ/kg** of actual whole-plant embodied energy takes **1.525 years** of export to repay. Manufacturing efficiency applies to a specified useful free-energy change, E = ΔG/η, with η ≥ 0.5. It does not set a universal fabrication cost per kilogram.

| Actual total embodied energy | Equivalent export years |
|---:|---:|
| 1 MJ/kg | 1.525 |
| 2 MJ/kg | 3.049 |
| 10 MJ/kg | 15.25 |
| 100 MJ/kg | 152.5 |
| 1000 MJ/kg | 1525 |

The receiver contains **5.32186e+21 kg C** under the stated carbon-fraction assumption. For CH₄ → graphite + 2 H₂ at 298.15 K and 1 bar, the standard useful chemical free-energy increase is **4.227 MJ/kg C**. [NIST-JANAF methane](https://janaf.nist.gov/tables/C-067.html).

At 50% efficiency, this specified chemical step requires 8.454 MJ/kg C and equals 2.054 years of this plant's export.
At 99% efficiency, this specified chemical step requires 4.27 MJ/kg C and equals 1.038 years of this plant's export.

Existing condensed carbon changes the starting chemical endpoint. Acquisition, selected material structure, construction machinery and delivery have separate accounts. Recycled material returning to the same bulk state can have a small net chemical free-energy change.

## Material access and operating supply

At resource C:H = 0.007, the receiver carbon corresponds to **7.60265e+23 kg of associated resource hydrogen**, at complete carbon recovery. That is a resource-access scale; only the selected products need be delivered. Services and reserve stocks need an elemental bill before total construction carbon can be stated.

Operation uses **2.77895e+06 kg/s fresh H** and **869.585 kg/s fresh C**, giving **C:H = 0.000312918**. The unused part of the 0.007 feed allowance is **18583.1 kg/s C** before separately identified terminal structural makeup. The former 0.0056 fuel / 0.0014 structural split is an allocation, while the complete material balance determines actual makeup.

Each kilogram of fresh hydrogen supplies **69.36 GWh** of net electricity. Retained skins, sorting tails and material in processing remain inventory; refurbishment throughput alone does not imply fresh carbon consumption.

## Factories and recurring work

| Construction duration | Average installed kg/s | Multiple of installed repair throughput |
|---:|---:|---:|
| 5 years | 2.11605e+14 | 2000 |
| 20 years | 5.29013e+13 | 500 |
| 100 years | 1.05803e+13 | 100 |

These are construction schedules supplied as scenarios. Energy efficiency does not establish a factory growth rate or deployment date. Delivery depends on site and trajectory, momentum exchange and braking recovery; a kinetic-energy scale alone is not a dissipative transport charge.

Routine repair handles **1.05803e+11 kg/s**, equivalent to **0.01%** of installed mass annually. The selected plant already charges its included refurbishment load. Fixed-throughput alternatives in the data report only the increment above that load; thermal resizing would require a new plant solve.

An additional 1 GJ/kg applied to all recovered hydrogen would draw **0.000465%** of reference export. This sensitivity does not replace geometric capture or processing-residence requirements.

## Economic comparison

With an illuminated collector at 30% efficiency, 1 kg/m² complete areal mass and 100% delivery efficiency, collector mass-specific output equals this plant at **140.2 AU**. The irradiance normalization is the IAU nominal 1,361 W/m² at 1 AU. [IAU 2015 Resolution B3](https://www.iau.org/common/Uploaded%20files/IAUGA2015-Resolution-B3-recommended-nominal-conversion.pdf). Equal-service comparison still needs siting, availability, storage and delivery apparatus.

Useful roles depend on the value of independent long-lived power and local resource access. Thermodynamic bounds establish energy scales; a money price requires the chosen supply chain, useful life, industrial throughput and opportunity cost. Numerical scenarios here are conditional inputs, not subjective probability intervals.

## Provenance

The companion JSON retains the exact plant selection, source hash, economic assumptions and arithmetic checks. Receiver capture, cold rejection, isotope flows and source preparation come unchanged from the selected upstream plant record. This calculation adds no operating debit or mass to that record.
