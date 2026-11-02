# ✨ Sealed fuel capsules: heat, pressure and mixture burden

Checkpoint, 2026-10-31. Allowing controlled melting removes the apparent need for every pellet to retain its initial solid phase. The nitrogen recipes have useful bulk enthalpy margins. Full-flight methane is possible only with a substantial pressure vessel or a changed preparation path; its added carbon naturally changes the delivered H:C ratio. Isotope-enriched polyethylene is a particularly useful simpler competitor for the C13 stage. None of these bulk screens qualifies the coating, neutralizer or three-dimensional assembly.

`screen.py` / `screen.json` preserve a reproducible calculation with 48 algebra, energy, domain and geometry checks. An independent read-only review confirmed the thick-sphere boundary stresses and the outer-area/field-energy scalings. They add no hydro trajectory and leave all existing source families unchanged. CoolProp 8.0.0 was installed in a task-local directory; its version/revision and input/code hashes are saved. The instructions in the script reproduce the result with that isolated runtime.

## Heat is density- and geometry-dependent

The reference is a **1 cm inner cavity**, a 1.1 billion metre radial flight at 2.4 million m/s, and the previously computed smooth-Al absorptivity 0.02182 for an isotropic 1000 K cavity. Flight time is 458 s. This is one representative speed; slower packets and the separately timed seed need their own exposure. For outer/inner radius ratio β and condensed-volume filling f,

\[
q_{IR}=3\alpha\sigma T_w^4R\beta^2/(\rho f a v),\quad
q_E=3p_E\beta^3/(\rho f),\quad p_E=\epsilon_0E^2/2.
\]

At 1 GV/m, p_E=4.427 MPa. For a fully filled 1000 kg/m³ core, the two doses are 170.1 and 13.28 kJ/kg. The electric term therefore rises to 28.47 kJ/kg for the 466.44 kg/m³ methane reference, before headspace. All stored field energy is conservatively deposited in the core for the bulk calculation. Electron escape, radiation, electrical recovery and reradiation would reduce that term. The core receives no credit for the shell's heat capacity.

The structural carbon shell needs a bonded reflective exterior; the Al optical result cannot be assigned to bare carbon. A 100 nm Al layer adds about 0.01% of core mass in the nitrogen/polymer examples. Its aluminum nuclei and recovery route belong in the actual source composition. They are absent from the present 27-species pulse network.

The absorptivity judgment already in [radiation-interface.md](../radiation-interface.md), central 0.022 and subjective 80% interval 0.007–0.06, remains conditional on a manufactured reflective surface. The present tables use its reference value. That optical uncertainty changes pressure and shell requirements materially; it is separate from the deterministic stress calculations below.

## Nitrogen fuels and solid polymer can use phase change

The table uses 20% cavity headspace. The selected cold densities are the stored source conventions, including isotope mass scaling; they are not new measurements of capsules. The shell pressure for NH3 and hydrazine is their small equilibrium pressure near melting. For polyethylene, 1 MPa is an explicit allowable internal-pressure **design scenario** requiring sufficient expansion space, not a measured EOS prediction.

| Core | Initial density / temperature | Full-flight IR + field dose | Bulk phase screen | Carbon shell/core at 1 / 10 GPa allowable |
|---|---|---:|---|---:|
| N14H3 | 850 kg/m³ / 150 K | 270 kJ/kg | About 46% melting near 195 K in the equilibrium enthalpy model | 1.97% / 0.196% |
| N15₂H4 | 1062.5 kg/m³ / 250 K | 216 kJ/kg | Below complete fusion energy; at most 58% melted if sensible heat is ignored | 1.57% / 0.156% |
| C12H2 polymer | 944.2 kg/m³ / 300 K | 243 kJ/kg | Around the melting interval for a sufficiently crystalline core | 2.17% / 0.216% |
| C13H2 polymer | 1011.64 kg/m³ / 300 K | 227 kJ/kg | Corresponding isotope-scaled polymer case | 2.03% / 0.201% |

[Overstreet and Giauque's primary ammonia calorimetry](https://doi.org/10.1021/ja01281a008) supplies the solid heat capacities used from 150 K to melting and approximately 5.65 kJ/mol fusion energy. The [indexed original table](https://electronicsandbooks.com/edt/manual/Magazine/J/Journal%20of%20the%20American%20Chemical%20Society%20US/1937%20%20%28vol%20059%29/02%20%20%28219-428%29/254-259.pdf) is transcribed in the script. [NIST](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7664417&Mask=4) gives a melting reference of 194.95 K and triple pressure 6060 Pa. The natural-ammonia liquid EOS near 195.5 K gives 734 kg/m³: even complete melting would fill only 92.7% of the chosen cavity. This supports a finite headspace solution without claiming a resolved transient temperature field.

Hydrazine's 12.657 kJ/mol fusion value and thermal-stability discussion are recorded in [the ACS supporting data](https://pstorage-acs-6854636.s3.amazonaws.com/3977605/op400120g_si_001.pdf), also traced to AIAA/NBS data. The isotope scaling gives 372 kJ/kg for N15₂H4. The dose is below that latent heat alone. Hydrazine's chemical decomposition can release several MJ/kg; compatibility with the actual liner and neutralization hot spots is consequently a real acceptance test. A low mean temperature does not qualify a catalytic surface. No synthesis or handling procedure is proposed here.

For polymer, the bounded phase model uses c_p=1.6 kJ/kg/K, a chosen 400 K melting interval and crystallinity 25%, 50% or 100%. A [direct fast-calorimetry measurement](https://doi.org/10.1016/j.polymer.2017.08.042) finds 281±6 kJ/kg for fully crystalline PE. At 50% crystallinity the reference dose reaches partial melting at 400 K; at 25%, the same piecewise model reaches about 408 K. The constant heat capacity and sharp melting point are screening assumptions. The actual chosen polymer needs calorimetry, thermal expansion and permeability data. Earlier 20→300 K calorimetry cannot be reused as the allowance of a core that starts at 300 K.

A self-supporting polymer core may carry some charging stress and reduce its separate structural skin. The table conservatively assigns the entire charge/pressure load to carbon. For same-isotope skins, the resulting polymer H:C ratios are **1.955 at 1 GPa and 1.995 at 10 GPa** for C13. This small composition shift is much easier to carry into a new pulse input than a large methane vessel. The new C13-polymer contact/hydro comparison is a separate task; no polymer capsule result has been substituted into the frozen methane source family.

A useful way to retain the mapped 300 K cold reference is to **pre-cool the polymer before launch**. The measured 20→300 K enthalpy range is 269–302 kJ/kg for the natural polymer samples, or about 251–282 kJ/kg with the same 14/15 isotope scaling used here. The 80%-filled C13 polymer reference absorbs 227 kJ/kg, below even the lower bound. It can therefore reach assembly below 300 K in the bulk equilibrium model. The actual optical coefficient and surface gradients still decide the margin. Cooling energy is small beside the approximately3 TJ/kg prepared pulse energy; the refrigerator and capsule-throughput apparatus remain source components. The 300 K-start melting option and the precooled option are distinct preparation paths.

## Methane and ethane pressure vessels

The natural-methane [Setzmann–Wagner EOS](https://doi.org/10.1063/1.555898), as implemented in [CoolProp](https://coolprop.org/fluid_properties/fluids/Methane.html), is evaluated at fixed cavity density and specific internal energy. C13 mass is represented by the same 17/16 volume-preserving convention used in the stored source material. This is not an isotope-specific molecular EOS. The original 102.4 K/439 kg/m³ liquid reference is already pressurized: the EOS gives 4.50 MPa. [Younglove's calorimetry](https://pmc.ncbi.nlm.nih.gov/articles/PMC6742817/) independently shows the strong pressure rise along this dense-liquid path. Headspace turns it into a different, initially two-phase trajectory.

For a thick elastic spherical shell with inner pressure p, outward electric traction p_E and allowable hoop stress S,

\[
X=(r_o/r_i)^3=(S+p)/(S-p/2-3p_E/2),\qquad
M_s/M_c=\rho_s(X-1)/(f\rho_c).
\]

The calculations use carbon density 2000 kg/m³. **1 and 10 GPa are required macroscopic allowables**, including defects, seams and temperature; neither is supplied by a qualified capsule material. Adding final pressure and charging traction is a conservative simultaneous load. Actual sequencing usually requires the larger of the charged cold-state and uncharged hot-state loads. Thermoelastic strain and membrane prestress require a subsequent material model.

The 10 GPa rows below iterate the increased outer absorbing area and stored electric energy consistently. Shell heat capacity and radiation cooling remain omitted. The 1 GPa entries are mechanical mass estimates at the uncorrected core heat dose; their large shell thicknesses require a new coupled thermal design.

| Methane fill fraction | Full-flight final T / p, 10 GPa geometry | Carbon skin/core | Actual H:C with C13 skin | 1 GPa mechanical mass/core |
|---|---:|---:|---:|---:|
| 1.0 | 307 K / 295 MPa | 19.5% | 3.19 | 221% |
| 0.8 | 331 K / 133 MPa | 11.1% | 3.49 | 117% |
| 0.5 | 428 K / 63.8 MPa | 8.80% | 3.59 | 90.3% |

These temperatures are compatible with molecular fluid states in the reference EOS. The design burden is containment, thermal contact and changed composition, rather than a universal requirement to remain below the normal boiling point. Full filling is especially expensive. A C12 skin instead of C13 changes the isotope mixture as well as total H:C; its exact added atom count is saved in the results.

A useful comparator is **C13-rich ethane**, initially a 100 K liquid. Its natural saturated-liquid density is 641 kg/m³, or 684 kg/m³ under 32/30 isotope scaling. The [Bücker–Wagner ethane EOS](https://doi.org/10.1063/1.1859286) gives, with the same 10 GPa geometry iteration:

| Ethane fill | Final T / p | Carbon skin/core | Actual H:C |
|---|---:|---:|---:|
| 1.0 | 279 K / 382 MPa | 17.3% | 2.47 |
| 0.8 | 276 K / 65.8 MPa | 3.87% | 2.86 |
| 0.5 | 320 K / 9.40 MPa | 1.21% | 2.96 |

Half-filled ethane preserves an approximately H:C=3 payload with modest vessel mass. Its lower average density increases capsule count/volume and charge per payload mass. A new cold-contact map is still required. Polymer avoids this volatile-fluid pressure path altogether and deserves priority when its whole-plant performance is close.

The existing **nominal128** carbon-median C13 runs have Q=24.5463 at H:C=3 and 24.6986 at H:C=4, at the same prepared energy. They establish nearby high-gain hot solutions. They do not prove the new shell-bearing cold trajectories. The root's integrated comparison finds the corresponding whole-plant difference much smaller than the pulse-gain difference, because C13 is only part of the cycle. Source simplicity can therefore outweigh the nominal H:C=4 optimum.

## Carbon inventory and local heat must be carried explicitly

The nitrogen skins add 0.00277 C12 atoms per N14 and 0.00222 per N15 at the 10 GPa reference. The 1 GPa counts are approximately ten times larger. The old fixed-history passive tracer sends most added C12 in A into the ordinary CNO routing and leaves almost all C12 in B1 unburned. For the 10 GPa A skin, its old tracer proton uptake would be about 0.0045 extra H per N14. This is evidence that the contamination has a tractable nuclear route; it is **not** a current source128 impurity tolerance or a nonlinear gain result. The tracer uses older mixtures/rates and has a finite native reconstruction error.

Against the unchanged source128 flow, A/B1/B2 skins would circulate 0.0266 kg carbon per kg net fresh H at 10 GPa, or 0.268 at 1 GPa. Those are gross circulating flows. They are not fresh consumption. Allocating the entire remaining 0.007 C:H allowance to this subtotal would allow terminal loss fractions of 25% and 2.5%, respectively; real receiver and methane/polymer-stage losses also use that allowance. A re-solved isotope ledger with actual shell nuclei is necessary before choosing a fresh-feed value. Coatings retained in a receiver or purification stock remain recoverable inventory.

The reference 0.33 ms neutralization time is a flux estimate, not a heat-transfer calculation. In a charge-limited thin carbon shell, assigning all field energy to the shell gives approximately **2S/ρ_s**, or 1 MJ/kg shell at 1 GPa and 10 MJ/kg at 10 GPa. With an illustrative 1 MJ/kg local sensible-heat allowance, the latter can retain only about 10% in the shell before exceeding that allowance. Primary-electron penetration, escaping electrons and contact with the core decide the partition. Slower discharge, a thermally bonded liner, lower surface field or more shell heat capacity are available design changes. The field energy must be counted once regardless of the route.

For an illustrative core diffusivity 10⁻⁷ m²/s, diffusion reaches about 6 micrometres during 0.33 ms and 6.8 mm during the full flight. Headspace can thermally isolate an exterior film and make it much hotter than the core. Conductive ribs, a wetted/attached core geometry or distributed small cells must supply the required contact. Gross melting margin does not establish surface integrity. Hollow-space collapse, shell fragmentation and the changed mean density also belong in the eventual three-dimensional merger calculation; they are not represented by silently keeping the old condensed-material entropy map.

## Protected exits have finite hot loads

A protected final path is an alternative, with hardware charged to the source. At only 20% of the full unshielded flight, 80%-filled methane reaches about 140 K/0.65 MPa in the 10 GPa model, needs a 0.408% carbon skin and delivers H:C=3.979. The last unshielded distance is 220 million metres. A hot exterior in equilibrium with the 1000 K cavity, multilayer insulation and a cold inner radiation shield is the physically coherent tube concept. An entirely cold exterior would continuously absorb a large cavity heat load.

For 100,000 aligned outlets of 250 m bore radius ending at r_tip=0.2R, the total cylindrical area is 1.38×10¹⁷ m², or 0.00909 m² per square metre of receiver. Illustrative 10–100 kg/m² tube/insulation areal allocations add 0.091–0.909 kg/m² at the receiver scale; 1–10 W/m² insulation leakage adds 0.0091–0.091 W/m². These exclude a resolved support, accelerator, cooling and terminal collector. Such tubes require distributed guidance/support rather than a rigid billion-metre cantilever.

The tip is exposed to outgoing pulses. Fluence grows as (R/r_tip)², and a homologous ballistic peak flux grows as (R/r_tip)³ because its duration also shortens. Using the actual coherent source128 C13 endpoint, the r_tip=0.2R loads are **1.26 GJ/m² kinetic, 63.6 MJ/m² advected-internal upper reference, 1.70 MJ/m² prompt photons, and 121 MW/m² peak kinetic flux**. The projected intercepted pulse fraction is N b²/(4r_tip²)=3.23×10⁻⁸, or 2.48×10¹⁹ J kinetic across all tips per C13 pulse. It is globally small but locally severe.

A separately armored hot shutter with volumetric/renewable stopping, or a charged-particle deflector and local catcher, is therefore an explicit terminal subsystem. An ordinary cold foil or bare graphite face does not inherit qualification at that concentrated flux. Incoming/outgoing tail overlap and shutter timing still need a real operating sequence. The saved half-radius example reduces these loads and shortens the protected tube. Protected flight trades capsule difficulty for a finite source-hardware problem; it is not free shielding.

## Next useful decision

Prefer a source-aware comparison of isotope-rich CH2 polymer, headspace ethane near H:C=3, and shell-bearing methane near H:C=3.5. Keep exact H:C=4 only if its small plant benefit pays for the corresponding capsule or protected-exit system. The immediate material discriminator is a coupled shell/core heat-deposition, phase/pressure and stress calculation with the manufactured reflector and neutralizer spectrum. That calculation should accept melting and finite pressure, while testing leakage, thermal contact, chemical compatibility and payload fragmentation. A successful average enthalpy calculation does not replace it.
