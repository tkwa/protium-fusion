# Final-family chronological dispatch and retained inventories ✨

The representative 350 kg/m² receiver operates at 4.15089401438×10³² net helium nuclei per second in the selected plant ledger. Using the actual delivered-B1/C13H2 source family, the unchanged chronological scheduler produces 4,228 pulses in a 370-day request horizon. Their main delivery/capture windows occupy 68.3917% of that horizon; the corresponding continuous-rate calculation is 68.3806%. The largest queueing delay is 83.51 minutes. The last capture finishes 298.66 seconds after the request horizon, and that completion is included in the inventory calculation.

The 250 kg/m² comparison is preserved separately. Its slightly larger calculated specific output, 20.7926 versus 20.7837 mW/kg, is a 0.043% difference below the transport sampling uncertainty. The main project selects 350 for its more favorable thermal response. This dispatch packet fixes the stated plant points and phase choices; it does not optimize receiver thickness or assign a statistical probability to those choices.

| Stage | Main incoming/outgoing window | 350 kg/m² nominal cadence | Pulses in 370 days | 250 kg/m² nominal cadence | Pulses in 370 days |
|---|---:|---:|---:|---:|---:|
| A: H3/N14 | 3,475 / 1,233 s | 9.384 days | 40 | 9.615 days | 39 |
| B1: H2/N15 | 3,652 / 1,529 s | 128.303 minutes | 4,153 | 131.458 minutes | 4,053 |
| B2: H2/C12 | 3,187 / 1,254 s | 11.438 days | 32 | 11.719 days | 32 |
| C13: H2/C13 | 3,545 / 958 s | 129.716 days | 3 | 132.905 days | 3 |

The radius is 1.1×10⁹ m. B1's incoming window uses the actual outer precursor profile, while its outgoing window uses the completed burn. The source-preparation radiation and driver energy remain in the final-family energy ledger; changing dispatch does not credit them again. Exact phase choices, every launch/focus/capture time, original input hashes, composition, material outputs and plant accounts are embedded in `sigma250/` and `sigma350/`. The integer pulse counts change total helium produced over this finite period by −3.85 ppm and +81.91 ppm respectively relative to the stationary-rate targets.

Each pulse requests a nominal focus time. In that order, its main launch begins after the preceding main capture window ends; an early request waits in the queue. This serializes the mass99.9% windows. It permits occasional closely spaced large pulses: the minimum focus separation between successive non-B1 shots is 9,911.15 seconds (2.75 hours), and their main-window gap is 5,181.68 seconds. A one-day large-pulse separation is not imposed. The finite-backing thermal assessment must use these actual gaps.

The outer 99.9% incoming mass and first 99.9% captured outgoing mass define these timing windows. The innermost incoming material therefore needs the separately proposed seed/protected delivery arrangement, whose hydrodynamic equivalence remains an engineering qualification. Slow outgoing material remains inside the closed material account, pending capture or re-entrainment. Window serialization does not establish that later packets can pass through those tails without unwanted interactions.

The inventory calculation conservatively withholds, for every recycled isotope separately, up to the entire omitted 0.1% pulse mass plus a 1% atomic/baryonic mass margin throughout the 370-day horizon. Adding those separate isotope bounds overallocates a physically common tail; it is a conservative stock screen rather than a predicted composition of escaped debris. No permanent material loss or fresh-carbon debit is assigned to that retained stock. Eventual hydrogen and CNO return fractions remain one. Nonquartet nuclear products stay in the held side-stream account of the material cycle.

At 350 kg/m², the calculated initial ready stock is 3.775×10¹³ kg for a one-day solid-processing hold, 1.238×10¹⁴ kg for 30 days, and 1.125×10¹⁵ kg for 365.25 days. Hydrogen has a separate selected one-hour local gas hold in all three cases. These are finite chronological ready-stock bounds under the specified delayed returns, not a qualification of gas capture or isotope-processing hardware. A steady annual rotating-skin plant also requires its resident age distribution. At the selected rate the gross fuel flow is 4.109×10⁷ kg/s, returned unburned hydrogen is 3.227×10⁶ kg/s, fresh hydrogen is 2.779×10⁶ kg/s, and fresh carbon is about 870 kg/s; the supplied carbon allowance remains mostly unused.

`mixture_schedule.py` and `material_cycle.py` are unchanged snapshots of the existing executable scheduler and isotope ledger. Re-run with the project's scientific Python environment:

```sh
python mixture_schedule.py --input sigma350/input.json --output rerun350.json --events-csv rerun350.csv --days 370
```

`python3 verify_dispatch.py` independently reconstructs the chronological events and cadence/accounting identities using only the archived inputs and CSVs. The scheduler's archived `events_csv` field is redirected to its adjacent `events.csv`; its original filename is preserved in the archive note. The independent review and source hashes are recorded beside the outputs.
