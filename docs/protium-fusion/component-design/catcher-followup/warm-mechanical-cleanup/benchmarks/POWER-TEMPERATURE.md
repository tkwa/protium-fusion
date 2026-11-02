# ✨ Measured pump power and temperature limits

This supplement preserves the earlier [benchmark results](results.json). [power-temperature.json](power-temperature.json) stores the specifications, exact source links and nominal-area arithmetic. The useful measured power anchor is the Edwards pump's **20 W at ultimate pressure**; the Pfeiffer source supplies a maximum rating only.

| Pump | Ultimate-pressure power | Rated or maximum power | Relevant thermal limits |
|---|---:|---:|---|
| Edwards nEXT85H ISO63, B8G410B01 | 20 W | 80 W total rated | Ambient 5–40°C; inlet flange 120°C; increase cooling when operating pump surface exceeds 60°C |
| Pfeiffer HiPace 700 H | Unspecified in sources inspected | 420 W maximum | Rotor 90°C; HV-flange bakeout 120°C; external incident heat 4.2 W maximum |

The [Edwards official product page](https://us.my.edwardsvacuum.com/en_US/USD/Catalog/Pumps/Turbomolecular-Pumps-Mechanical/nEXT-Pumps/nEXT85/p/B8G410B01) directly lists the 20 W and 80 W values, 24/48 V supply options, and 5–40°C ambient range. Its specification does not partition the ultimate-pressure draw among bearings, motor and controller. Cooling accessories and the backing pump require their own account.

The [official Edwards B8G000880_H manual](https://us.my.edwardsvacuum.com/medias/B8G000880-Instructions.pdf?attachment=true&context=bWFzdGVyfHJvb3R8Mjc0NzAxNnxhcHBsaWNhdGlvbi9wZGZ8YUdOakwyZzFZUzh4TVRneU56ZzROakEzTXpRME9UUXZRamhITURBd09EZ3dJQzBnU1c1emRISjFZM1JwYjI1ekxuQmtaZ3w2M2NkMjBlNGQ0NjRhNGEzMzAxYWE5MTU2ZjY2NGE1MDNmMGY2YzkwYjc3MjNmMTdlODk0MDM5YThkMmVkMWVi), printed pages 20, 56 and 64, supplies the 120°C flange limit and the 60°C surface cooling instruction. It permits an adjustable **maximum** power setting of 50–120 W, with 80 W default. That setting limits the controller; it does not impose 50 W continuous consumption. The analogue rotor-temperature display spans 0–100°C, while the serial interface has a different reporting range. Neither display range establishes a permissible rotor temperature. An explicit continuous rotor maximum was not identified. The flange limit likewise does not qualify the whole housing or rotor at 120°C.

The [Pfeiffer official manual](https://www.pfeiffervacuum.com/de/assets/123614/1/Instruction%20Manual%20HiPace%20700%20H%20with%20TC%20400_DE_de%2CEN_en.pdf), English printed pages 79 and 121, gives the 420 W consumption ceiling, 48 V supply and the thermal limits tabulated above. Ambient operation is 5–35°C with air cooling or 5–40°C with water cooling. Controller gas-mode vertices and backing-pump switching thresholds are not idle-power measurements.

Using nominal circular openings, DN63 gives 0.003117245 m² and the Edwards ultimate draw is **6416W per nominal inlet m²**; its rated value gives 25664 W/m². DN160 gives 0.020106193 m² and the Pfeiffer maximum corresponds to 20889 W/m². These are complete small-device comparisons. Hub blockage, motor sizing, module scale and controller overhead prevent treating them as invariant losses per future aperture area.

For the proposed warm first stage, the material and heat-flow design must support 1,000 K gas-facing structures while assigning actual rotor, bearings, motor and electronics temperatures. These specifications neither measure such a stage nor establish its minimum parasitic power. The blade-transport calculation can constrain capture and compression; a separate rotor-loss model must close the electrical budget.
