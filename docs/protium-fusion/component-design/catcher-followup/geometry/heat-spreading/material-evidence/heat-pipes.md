# ✨ Primary sodium heat-pipe evidence and property inputs

Sodium phase-change transport near 1,000 K has experimental precedent. The evidence supports developing a distributed heat spreader; it does not by itself establish the receiver's transient conductance, mass or capillary margin. Two full primary reports were inspected. Their measured quantities need different interpretations.

## Tests and their useful scope

**Rosenfeld et al., NASA/TM-2004-212959.** The one-tenth Stirling heater-head segment accumulated 42,700 operating hours by September 2003; the report describes the life operation near 700°C (973 K). Its Table 1 lists a **design** operating point of 1023 K, 4500 W throughput, 200 kW/m² condenser flux and 100 g sodium charge. These design numbers must not be combined into a claim of 42,700 hours at that full point. A separate 316L/nickel-wick bench unit operated more than 115,000 hours near 650°C; its heaters were rated at 4 kW. Liquid return there used gravity. Another otherwise similar unit failed after 200 hours because a defective heater braze caused local overheating and a leak. This supplies both positive compatibility evidence and a concrete interface failure mode. The report does not give a transferable measured whole-device conductance or receiver areal mass. See printed pp.3–6, Tables 1–2. [Full NASA report](https://ntrs.nasa.gov/api/citations/20130013063/downloads/20130013063.pdf).

**Boo, Kim and Kang, Energy Procedia 69 (2015), 608–617.** This sodium loop used roughly 0.5 m source-to-condenser separation and gravity-assisted return, usually at 70° inclination. The useful loaded point is estimated condenser output **671 W**, thermal resistance **0.13 K/W**, and heated-face flux **57.4 kW/m²**, with 1100 W heater input. Thus device conductance is about **7.69 W/K**. Output was inferred from surface temperatures and convection/radiation correlations. The quoted 47,000 W/(m·K) effective value is normalized to vapor-line cross section, not to the full receiver area. The striking 902,000–973,000 W/(m·K) and 0.006 K/W results used an insulated condenser and heater input despite an unresolved condenser/other-loss split. They cannot certify an 800 W useful load at a 4.8 K drop. The tests reached heater limits rather than measuring a capillary limit. Frozen startup took roughly 90–140 minutes. No long-duration or orientation-independent performance is established here. [Paper DOI](https://doi.org/10.1016/j.egypro.2015.03.070), [full primary-paper text](https://www.researchgate.net/publication/279215308_An_Experimental_Study_on_a_Sodium_Loop-type_Heat_Pipe_for_Thermal_Transport_from_a_High-temperature_Solar_Receiver), §§1–4 and Table 3.

## Sodium inputs at 1,000 K

The following are evaluations of published correlations, rather than new measurements:

| Property | Value at 1000 K | Source |
|---|---:|---|
| Saturation pressure | 19,218.5 Pa | SAS4A/SASSYS equations 12.13-2 |
| Vaporization enthalpy | 4.01517 MJ/kg | Equation 12.13-1 |
| Saturated vapor density | 0.0582399 kg/m³ | Equation 12.13-6 |
| Liquid density | 779.284 kg/m³ | Equation 12.13-5 |
| Liquid viscosity | 1.858273×10⁻⁴ Pa·s | Equation 12.13-12 |
| Liquid surface tension | 0.135455 N/m | Fink–Leibowitz §1.6 |
| Saturated vapor viscosity | 1.86965×10⁻⁵ Pa·s | Golden–Tokar equation 5.20 |

The first five evaluations use the [official Argonne SAS4A/SASSYS sodium property documentation](https://sas-doc.nse.anl.gov/5.4/Part04/Ch12/12.13.html). Its polynomial approximations cover 590–2270 K. Quoted approximation errors relative to the underlying correlations are approximately 0.5–1.5%; these are not complete experimental uncertainties. The older safety-analysis correlations and the later Fink–Leibowitz assessment should retain separate provenance when more precision is needed.

For surface tension, [Fink and Leibowitz, ANL/RE-95/2](https://digital.library.unt.edu/ark:/67531/metadc793620/m2/1/high_res_d/94649.pdf), §1.6, use

    sigma[N/m] = 0.2405 (1 − T/2503.7)^1.126.

Their table gives 135.4 mN/m at 1000 K. The reported uncertainty is ±11% at two standard deviations across the measured 371–1600 K range. Wetting/contact angle remains an independent interface property; using 2 sigma/r as capillary pressure assumes complete wetting.

For vapor viscosity, [Golden and Tokar, ANL-7323](https://digital.library.unt.edu/ark:/67531/metadc1033805/m2/1/high_res_d/4511962.pdf), printed pp.51–52, fit smoothed Stefanov data as

    mu[lb/(h ft)] = 0.03427 + 8.176e−6 temperature[°F].

At 1000 K, temperature is 1340.33°F; multiply by 0.45359237/(3600×0.3048) to obtain Pa·s. The reported experiments span 686–870°C, which includes 1000 K, but the assessment notes missing apparatus/calculation details. This is a useful preliminary friction input with weaker provenance than a modern precision transport measurement. The report distinguishes saturated-mixture viscosity from its monatomic-vapor calculation.

## Application to the receiver

The current receiver requirement is a transient local heat-spreading problem. A candidate must carry the actual pulse through its evaporator contact, vapor space, condenser contact and liquid-return path while retaining the backing's storage capacity. A high effective conductivity assigned to the whole slab is a requirement proxy until those interfaces close.

A heat-pipe network adds envelope, wick, fluid and contact mass; replacing graphite with void changes shielding and heat capacity. Its capillary test is available meniscus pressure versus the complete liquid/vapor pressure loss and any acceleration head. Vapor friction, entrainment, sonic transport, local boiling and dryout must be evaluated at pulse loading. None of the cited tests measures that receiver-specific limit. Maintaining a hot standby state can avoid a frozen start on every pulse, while cold startup and local repair still need a procedure.

No probability distribution over whole-device performance is inferred from these two demonstrations. Their value is a physically realized transport mechanism, documented long-lived material combinations, and specific interface/return-path requirements for the next design calculation.
