# Launcher ports and shielding ✨

The provisional **1,000 launchers** can each emit many packets. This number defines an organizational and hardware partition; the present spherical pulse calculation does not validate it as the number of independently controllable angular sources. A useful reference separates chamber-facing emitters, supplied through protected service passages, from emitters outside the receiver that send full-speed fuel through open receiver ports. These arrangements have different aperture and shutter requirements.

All numbers below use `R = 1.1e9 m`. They are deterministic calculations or explicit scenarios, not statistical estimates. The 100/1,000/10,000 count scan is a sensitivity study. There is presently no defensible central estimate and subjective 80% interval for the optimum count, muzzle radius, emittance or shutter reliability: the input design that would support such judgments is missing.

`geometry.py` reads the frozen `open-loop-redesign/station-inputs.json` and `collector-inputs.json`, writes `geometry.json`, and records input and code hashes. Reproduction from this directory:

```sh
python3 geometry.py --root /absolute/path/to/website
```

## Territory, emitting area and aperture area

For `N` equally allocated stations on the spherical receiver,

```
A_receiver = 4 pi R² = 1.520530844e19 m²
A_cell = A_receiver/N
Omega_cell = 4 pi/N
r_equal_area = sqrt(A_cell/pi) = 2R/sqrt(N)
theta_cap = acos(1 - 2/N)
s_characteristic = sqrt(A_cell)
```

| Launcher count | Area per allocation, m² | Equal-area disk radius, km | Characteristic spacing, km |
|---:|---:|---:|---:|
| 100 | 1.521e17 | 220,000 | 389,940 |
| 1,000 | 1.521e16 | 69,570 | 123,310 |
| 10,000 | 1.521e15 | 22,000 | 38,994 |

The spacing is a scale derived from area, rather than a solved spherical packing distance. The cell radius is the radius of an equal-area flat disk; the corresponding spherical-cap radius for 1,000 allocations is 69,582 km. These are service territories, not apertures or emitter footprints.

A compact point emitter aimed into a small central volume supplies a narrow set of approach directions. Firing additional packets from it supplies additional material and timing freedom; this alone does not fill its entire angular service cell. Broad directional coverage can instead come from distributed emitting sites within a launcher installation, relocated launch positions, or a separately accounted steering/relay system. For example, the cap half-angle for 1,000 equal territories is 3.624 degrees, while a 170 m central assembly radius subtends only `170/R = 1.55e-7 rad` at the receiver. Confusing these angles would hide the angular assembly problem.

Let total active emitting area be `A_emit`. If it is fixed by required throughput and qualified areal delivery capability, partitioning it between `N` launchers gives `A_emit/N` per launcher. It does **not** multiply `A_emit` by `N`. The same rule applies to a fixed total clear aperture area `A_port`. It does not establish that `A_emit = A_port`: beam transport, focusing, incidence angle, supports, shutters and stand-off geometry determine that mapping.

For example, hold `A_port = 1000 pi (1000 m)² = 3.14159e9 m²`. Per-port radii become 3.162 km, 1 km and 0.3162 km for 100, 1,000 and 10,000 ports. Total projected opening and equivalent plug mass remain unchanged. Total perimeter grows as `sqrt(N)`, from 1,987 km through 6,283 km to 19,869 km, so edge seals, mechanisms and service connections can still penalize a larger count.

## Permanently open, full-speed injection ports

For small circular ports normal to radial rays, isotropic emission from a central point gives

```
f_open = sum(A_projected)/(4 pi R²) = N a²/(4R²).
```

For a planar disk at distance `R`, the exact solid-angle fraction is `N[1 - 1/sqrt(1+(a/R)²)]/2`; the numerical implementation uses an algebraically equivalent expression without subtractive cancellation. Oblique openings use projected area, and extended/asymmetric sources require a view-factor calculation. The following fractions assume every ray entering a port escapes the station and is not recovered outside it.

| Radius of each of 1,000 ports | Fraction of isotropic output entering ports | Share of entire 1% energy-loss budget | Share of entire model fuel-loss budget |
|---:|---:|---:|---:|
| 1 km | 2.066e−10 | 2.066e−8 | 5.548e−7 |
| 10 km | 2.066e−8 | 2.066e−6 | 5.548e−5 |
| 100 km | 2.066e−6 | 2.066e−4 | 0.005548 |
| 1,000 km | 2.066e−4 | 0.02066 | 0.5548 |

The selected fresh-carbon accounting uses 0.007 total C:H by mass, including 0.0056 for fuel processing. Its common CNO/hydrogen nonreturn fraction is **0.0003723931 per cold-stage output**. For isotropic, unselected matter, letting all this nonreturn occur through the openings gives

```
a_full_budget = 2R sqrt(f_budget/N).
```

| N | Radius using all 1% energy-loss budget, km | Radius using all 0.0003723931 fuel-loss budget, km |
|---:|---:|---:|
| 100 | 22,000 | 4,245.45 |
| 1,000 | 6,957.01 | 1,342.53 |
| 10,000 | 2,200 | 424.545 |

These are **maximum full-budget diagnostics**. The existing loss allowances already cover energy capture and process nonreturn; they are not spare port allocations. A real port assessment needs the unallocated residual budget, time-dependent open fraction, angular mass/energy distribution, and any external catcher's recovery efficiency, mass and power. Temporarily closing a port restores capture only if the closure itself functions as part of that receiving system.

Small global escape fractions do not make the exposed accelerator hardware lightly loaded. The C13 reference pulse supplies 48.850 kJ/m² of prompt radiation and 34.155 MJ/m² of kinetic fluence at this radius, plus separate residual internal and decay energy. A single 1 km-radius opening admits **1.535e11 J prompt** and **1.073e14 J kinetic** on the isotropic assumption. Delayed and scattered radiation, secondary particles and vapor can also reach the opening after the initial pulse. The saved pulse spectrum does not yet determine a qualified shield.

## Chamber-facing launchers with protected refill

A launcher mounted on the chamber-facing side can accelerate and emit fuel without sending the fast packet through the main receiver thickness. Its feed, power and maintenance passages can follow shielded, non-line-of-sight routes behind the emitting head. Slow feed can be transferred into local buffers using paired isolation valves or a staged transfer lock; the emitting head then operates while the refill path is closed. This is a concrete architectural option, rather than a demonstrated valve design or a claim that all front-side apparatus is protected.

The chamber-facing nozzle, steering surfaces, supports and exposed shield still receive the local pulse. They need renewable or cooled receiving surfaces and a closed muzzle or equivalent protection during the outgoing phase. Their projected shadow and intercepted fuel need a return route. A launcher installation can contain many muzzles, gates and packet channels without changing the organizational count of 1,000.

A passage carrying slow feed can have bends and radiation traps. A straight ballistic path for a full-speed neutral packet cannot simultaneously be blocked by a fixed baffle. Charged-beam steering can remove line of sight only with a specified field/charge-state transport section; neutral pellets do not automatically follow that path. Multiple shutter stages can separate prompt-radiation shielding, later debris interception and vapor isolation. Fusion research has studied these functions separately: Kondo et al. propose multiple shutters to control chamber-gas backflow in heavy-ion systems, while LLNL reports continuing debris-shield and gap-management work in its actual final optics. These examples support the engineering functions, not this reactor’s performance. [Kondo et al.](https://arxiv.org/abs/1608.01063), [LLNL optics and photonics](https://engineering.llnl.gov/research/optics).

The 1,000 K receiver is also a persistent thermal environment. A black surface viewing an equilibrium 1,000 K cavity receives `sigma_SB T⁴ = 56.7 kW/m²` before accounting for its own emission, reflection or restricted view. Cold feed protection and local heat rejection therefore remain necessary between pulses.

## Closing before the return pulse

The central prompt signal reaches the receiver after `R/c = 3.6692 s`. In the prescribed incoming profiles, the fastest outer material implies a last ballistic emission roughly 237–268 s before central assembly. Under that specific radial schedule:

| Stage | Last emission before focus, s | Time from last emission to prompt arrival, s | First endpoint-ballistic outgoing ions after focus, s |
|---|---:|---:|---:|
| B2 | 237.19 | 240.86 | 61.81 |
| C13 | 262.46 | 266.13 | 58.50 |
| A | 256.19 | 259.86 | 80.79 |
| B1 | 267.80 | 271.46 | 118.05 |

This schedule allows preplanned closure after the final incoming packet. It does not establish a mechanism or permit waiting for the flash at the receiver before blocking that same flash. Earlier/slower material may already occupy the chamber; later source corrections, different acceleration locations and packet train lengths change the available interval. The actual packet tail and clearance sensor must confirm that a closing gate does not intercept incoming fuel. Slow outgoing material and local vapor require a longer protection/cleanup schedule than the first-ion arrival time.

As a scale check, a monolithic shield with radius `a`, areal mass `Sigma = 5000 kg/m²`, and stroke `2a`, closing from rest to rest in `T = 250 s` with a triangular velocity profile, has

```
m = pi a² Sigma
v_peak = 4a/T
acceleration = 8a/T²
K_peak = m v_peak²/2 ∝ a⁴.
```

At `a = 1 km`, this means 1.571e10 kg, 16 m/s peak speed and 2.011e12 J peak translational energy per plug. At `a = 1000 km`, it means 1.571e16 kg, 16 km/s and 2.011e24 J. Kinetic energy can be recovered; peak actuator power, precision, momentum transfer and stored energy still matter. This monolithic example is not a lower bound. Many short-stroke tiles, rotating closures or a shutter at a narrower beam section can change the scaling, with their own gaps, redundancy, heat paths and failure requirements.

Closing with a 10 nm membrane does not substitute for the absorber: the separate ion screen removes only about 0.2% of C13 kinetic energy in intact silica of that thickness. A plug intended to restore full receiving coverage needs the relevant stopping, heat transport and shielding functions. At 10 and 20 MeV, NIST's carbon attenuation coefficients give uncollided transmissions `exp[-(mu/rho) Sigma]` of 5.57e−5 and 3.80e−4 for 5000 kg/m²; achieving 1e−6 uncollided transmission would require 7,052 and 8,772 kg/m². These are diagnostic monoenergetic ray attenuations, not dose, full shower containment or proof of the prompt spectrum. Scattered radiation and neutron transport require their own analysis. [NIST carbon coefficients](https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z06.html).

## Inputs that could favor another launcher count

A count decision needs the following inputs together: packet mass/size and train length; emitter footprint and number of independently directed channels; the covariance or transverse emittance of each transported beam; allowable halo/interception and focus error; steering range and energy spread; aperture dimensions and stand-off; closure time and gap leakage; qualified areal launch power and heat rejection; repair access, redundancy and failure correlations; and the acceptable angular perturbations of the converging flow.

For a drift, the transverse variance obeys the directly reproducible relation

```
sigma_x²(L) = sigma_x0² + 2L Cov(x0,theta0) + L² sigma_theta0².
```

A clear aperture must contain the chosen fraction of this distribution plus centroid error, pointing jitter and mechanical margins. A bound-pellet trajectory uses its size, attitude and tracking errors rather than assigning it a thermal-gas emittance. Neither aperture follows from the hydro pulse's 30 m reference radius alone.

More independently positioned emitting sites may improve angular control and failure tolerance; fewer may reduce repeated mechanisms and service overhead. At fixed total active area, smaller apertures may ease closure while increasing total edge length. Only a three-dimensional assembly calculation, a transport/aperture design and a mass/power/reliability comparison can distinguish those effects. The present evidence supports **retaining 1,000 as a provisional reference**, with multiple packets and potentially multiple emitting channels per launcher; it supplies no evidence that 100, 1,000 or 10,000 is optimal.
