# ✨ Primary benchmarks for warm mechanical gas cleanup

Current molecular pumps establish useful capture and compression scales. A 1,000 K inlet rotor and its bearings, clearances, radiation enclosure and downstream stages require a new design. [calculate.py](calculate.py) and [results.json](results.json) preserve the source values and the conditional calculations below.

## Official hardware data

| Model | Nominal inlet | H₂ / He pumping speed | H₂ / He maximum compression | Rotation | Pump mass |
|---|---:|---:|---:|---:|---:|
| Pfeiffer HiPace 700 H | DN160 | 555 / 655 L/s | 2×10⁷ / 10⁹ | 49,200 rpm | about11.8kg |
| Edwards nEXT85H ISO63 | DN63 | 54 / 78 L/s | 5×10⁵ / 2×10⁷ | 90,000 rpm | about3.0kg |

Sources: [Pfeiffer official leaflet](https://www.pfeiffervacuum.com/global/assets/134975/1/Product%20Leaflet%20HiPace%20700%20H%20Global%20EN.pdf), [Pfeiffer operating manual](https://www.pfeiffervacuum.com/de/assets/123614/1/Instruction%20Manual%20HiPace%20700%20H%20with%20TC%20400_DE_de%2CEN_en.pdf), [Edwards official datasheet](https://www.edwardsvacuum.cn/content/dam/brands/edwards-vacuum/general-vacuum/downloads/turbomolecular-pumps/3601-0243-01-nEXT85-datasheet.pdf), and [Edwards speed specifications](https://digitalbrochure.edwardsvacuum.com/nEXT-TURBOMOLECULAR-PUMPS/en/technical-specifications.html). These are ordinary operating specifications; the table gives no hot-rotor qualification. DN is a nominal diameter, used below as a transparent aperture approximation.

## Compression at finite flow

The free-molecular linear balance is

\[
Q=S_0\left(p_{\rm inlet}-p_{\rm outlet}/K_0\right).
\]

For each species separately, a 10⁻³ Pa outlet and 10⁻¹⁰ Pa inlet require a pressure ratio10⁷. Setting \(K_0=10^7\) leaves zero forward throughput. Retaining50%,90%,99% of the negligible-backpressure speed requires \(K_0\ge2\times10^7,10^8,10^9\), respectively. Consequently the room-temperature HiPace example retains50% of its H₂ speed and99% of its He speed at those chosen partial pressures. Total chamber pressure must not be substituted independently for both partial pressures.

[Pfeiffer's manufacturer theory](https://www.pfeiffervacuum.com/no/en/knowledge/vacuum-technology/knowledge-book/4-vacuum-generation/4.9-turbomolecular-pumps.html) gives the exponential dependence on blade speed divided by molecular thermal speed and warns that interchannel leakage reduces ideal Holweck compression. Its effective-speed approximation for heavy gases should not be applied to H₂; the aperture comparison here uses the actual light-gas specifications.

## Capture area and temperature

Using \(\alpha=4S/(A\bar v)\), the HiPace nominal opening implies H₂/He effective captures0.0622/0.1034 at a 300 K reference. Holding its specified speeds fixed while admitting 1,000 K gas gives0.0341/0.0567, requiring nominal mouth fractions0.294/0.176 to achieve global \(g=0.01\). This is a comparison, not a hot performance prediction or a lower bound on future equipment. Net capture includes return through the inlet, hub obstruction and internal transmission. A 1% mouth supplies \(g=0.01\) only if its net capture is unity.

For an unchanged geometry with *every relevant stage* heated300→1000K, the similarity estimate \(\ln K\propto u/\sqrt T\) reduces the HiPace compression scales to about10⁴(H₂) and8.5×10⁴(He). Maintaining the same ratio instead requires about1.826 times the speed, or the effective logarithmic compression length at fixed speed. A hot first stage followed by shielded cooler stages must be evaluated stage by stage; the entire train need not suffer the all-hot factor.

A [primary Pfeiffer R&D presentation at CERN](https://indico.cern.ch/event/161398/contributions/235875/attachments/187008/262586/Turbopumps_from_theory_to_practice.pdf), slide12, gives a300–400m/s blade-speed scale. Its indexed slide text was accessible; direct PDF retrieval returned403. Scaling that reference gives548–730m/s for preserved \(u/\sqrt T\). These speeds specify a material/rotordynamic requirement; they are not demonstrated1,000K pump ratings.

## Useful next calculation

A unit-cell test-particle calculation should launch1000K flux-Maxwellian H₂ and He into the actual rotating blade geometry, retain inlet returns, and tally forward/reverse transmission at several rotor speeds. Couple hot inlet cells to the proposed cooler or hot drag stages and compute finite-flow compression. Include clearances and the same optical baffling used by the thermal design.

A relevant [primary DSMC study](https://arxiv.org/abs/physics/0401028) treats a Holweck pump with moving geometry and compares collisionless and collisional calculations. Its helium example uses400rps and13.3→133Pa, finding about20% less flow when collisions are included. Those pressures differ greatly from the present first stage; the20% is evidence for checking the regime, not a loss factor to copy. No subjective performance interval is assigned before a blade geometry and accommodation model exist.
