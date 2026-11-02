# Incoming timing and discrete source directions ✨

The reference is **1,000 launcher complexes at a 1.1 × 10⁹ m radius**. Each may emit many packets. The frozen corrected64 calculations constrain a spherical radial mass, velocity, temperature and composition profile. They provide no direct count of launchers, physical packets or independent angular channels. This note reconstructs their incoming timing and gives conditional geometry comparisons. No new hydrodynamics was run.

## Flight time and sequence

For material at radius (r), inward speed (u), and hot-snapshot time zero, ballistic launch from radius (R) occurs at

\[
t_{\rm launch}=-(R-r)/u.
\]

Slow inner material leaves first. Faster outer material follows. Within each saved cell, this extraction preserves its mass, assumes uniform spatial density, and interpolates nodal velocity linearly in radius. The flight-time quantiles below describe the incoming **mass**, with the longest flight times corresponding to the earliest launches.

| Stage | RMS speed, km/s | Median flight, s | Launch window containing the middle 90% of mass, seconds before snapshot | 99.9th-percentile flight, s |
|---|---:|---:|---:|---:|
| A | 2,397 | 484 | 320–933 | 3,437 |
| B2 | 2,588 | 448 | 297–864 | 3,183 |
| C13 | 2,339 | 496 | 328–956 | 3,521 |
| B1 | 2,273 | 512 | 336–1,012 | 3,771 |

C13 therefore requires roughly ten minutes of scheduled delivery for its middle 90% of mass. Its innermost resolved cell centroid takes 8.35 hours to arrive. The other stages' corresponding times span 7.55–8.95 hours. A continuous central profile with (u\to0) has an arbitrarily long flight-time tail; the centroid is a discretization diagnostic.

The hot initial outer radii are 170–191 m. Straight-line boundary focus times after this snapshot span 26.6–46.5 μs across the stages. The supplied 30 m reference radius is a scaling parameter; pressure and burning determine the subsequent actual compression.

A backward extrapolation reaches maximum averaged density 2,000 kg/m³ at C13 outer radius 1,755 m, 0.378 ms before the hot snapshot. A and B2 give nearly the same radius; B1 gives 406 m and 0.0523 ms. These are **packing-density backmaps**. They do not locate the first collision: crossing streams can collide at much lower volume filling fraction, and pressure would alter a sufficiently collisional extrapolation. For equal hard grains of radius (a), solid density (\rho_s), and local bulk density (\rho_b), the geometrical mean free path is (\lambda=\rho_s a/(3\rho_b)), with collision frequency (n\sigma v_{\rm rel}). Stream occupancy, relative velocities, fragmentation, ionization and collective stopping are needed to turn this into a merger calculation.

## Angular channels and the count comparison

Repeated packets from one compact source refine its time and mass profile. Their trajectories toward a region only a few kilometres across differ by only a few microradians when launched from (R=1.1\times10^9) m. Additional approach directions require additional source positions, spatially distributed outlets within a complex, or trajectory changes between launch and merger.

For (N) uniformly distributed axes, an equal-area spherical cap has half-angle (\alpha=\cos^{-1}(1-2/N)\), while the characteristic axis spacing is approximately (\sqrt{4\pi/N}\). The following comparison holds total C13 pulse mass and prepared energy fixed and divides the source load equally. Peak electrical allocations retain the frozen source-power reconstruction and its assumed 60% driver efficiency.

| Angular source count | Cap half-angle | Typical axis spacing | Patch radius at the 1,755 m backmap | Mass per source per shot, kg | Peak prepared electrical power per source, W |
|---|---:|---:|---:|---:|---:|
| 100 | 11.48° | 20.31° | 349 m | 1.13 × 10¹¹ | 2.10 × 10²¹ |
| 1,000 | 3.62° | 6.42° | 111 m | 1.13 × 10¹⁰ | 2.10 × 10²⁰ |
| 10,000 | 1.15° | 2.03° | 35.1 m | 1.13 × 10⁹ | 2.10 × 10¹⁹ |

At 1,000 sources the equal-area sector corresponds to a transverse extent of about 69,500 km on the launch sphere. This describes a distributed sector's geometry; a launcher complex's physical outlet arrangement still needs to be specified.

One conditional energy diagnostic treats each sector as a parallel incoming flow which must lose its transverse motion relative to the locally radial direction. For mass uniform per solid angle within a cap,

\[
\langle\sin^2\theta\rangle={2\over N}-{4\over3N^2},\qquad
 e_\perp={v_{\rm rms}^2\over N}\left(1-{2\over3N}\right).
\]

For C13, (N=100,1{,}000,10{,}000) gives 54.4, 5.47 and 0.547 GJ/kg respectively. These are 1.99%, 0.200% and 0.0200% of initial kinetic energy, or **26.3%, 2.65% and 0.265% of the prescribed initial internal energy**. The saved C13 internal-to-kinetic ratio is 0.07542. The angular energy may contribute to the prescribed preheat budget; adding it again would double-count that contribution. Its deposition location and timing matter because heat generated early can increase entropy and resist compression. Shaped flows, staged mergers and transverse-energy recovery change the diagnostic.

This scaling motivates comparing 100, 1,000 and 10,000 genuinely different incoming directions. It neither supplies an optimized count nor establishes a yield penalty. A credible interval for a required minimum count would require a specified channel layout and a calibrated relation between merger errors and yield. The [proposed 95% yield-retention criterion](literature.md#proposed-comparison-criterion), at fixed prepared energy, supplies a concrete next comparison rather than an invented count posterior. These deterministic geometry results have no subjective credible interval attached.

## Packets, grains and tolerances

The C13 shot contains 1.131 × 10¹³ kg. At 1,000 launchers, each supplies 1.131 × 10¹⁰ kg. One packet per launcher assigns that whole mass to one packet; 1,000 packets per launcher gives one million packets per shot, each 1.131 × 10⁷ kg. A packet can itself contain many grains. Illustratively, 1 cm radius grains at 1,000 kg/m³ weigh 4.19 g, requiring 2.7 × 10¹⁵ grains per shot and about 3.87 × 10⁹ grains/s per launcher during the middle-90% delivery window. Neither this grain size nor this launcher throughput has been qualified. The 64 radial numerical cells are unrelated to physical packet count.

An uncompensated positional error budget (b) gives angular error (b/R), timing error (b/u), and fractional speed error approximately (b/R) over a complete free flight. With C13 RMS speed, (b=10) m corresponds to 9.1 nrad, 4.27 μs and (9.1\times10^{-9}) respectively. The choice of 10 m is an illustrative criterion. In-flight measurement and correction can trade these launch tolerances against later controls. The spherical calculation supplies no verified acceptable (b), correlation pattern or emittance.

## Collision overlap and cadence

For a preceding shot's outward material of speed (w), a following shot's inward material of speed (u), and focus separation (\Delta), the ballistic meeting point is

\[
t_{\rm cross}={u\Delta\over u+w},\qquad
r_{\rm cross}={uw\Delta\over u+w}.
\]

They meet inside the receiver if (\Delta<R/u+R/w). For illustrative C13 speeds (u=2.339\times10^6) and (w=1.166\times10^7) m/s, this threshold is 565 s. A 60 s focus separation gives a meeting radius (1.17\times10^8) m. Slower tails impose longer separations. Increasing launcher count alone leaves this radial overlap condition unchanged.

Within one incoming sequence, a later, faster packet launched at (t_2) behind a packet launched at (t_1) meets it at

\[
r_{\rm catch}=R-{u_1u_2(t_2-t_1)\over u_2-u_1}.
\]

A positive crossing radius identifies a kinematic encounter; collision probability and dissipated heat require the streams' actual transverse overlap. Desired convergence itself eventually produces such encounters.

The selected frozen `high_source_rating_surface` case has mean all-stage interval **6,662 s (111.0 minutes)** and B1 interval 6,774 s. Its mass-99.9% incoming-plus-outgoing occupancy screens are 72.6 minutes for B2, 79.9 for A, 80.0 for C13 and 89.0 for B1. These are compatible average allocations at the case's selected 80% occupancy; an executable mixed-stage schedule must allocate individual gaps. They retain the receiver's older uniform-mass-per-speed reconstruction, whose incoming 99.9% quantiles are about 0.6% longer than the spatial reconstruction above.

Full resolved centroid occupancies are 10.3–12.6 hours. The remaining tails therefore need explicit recovery, protected routing or altered preparation. Material outside a 99.9% scheduling window remains part of the C:H return and station-power accounting. Neither average pulse rates nor an angular count establishes a collision-free recurring schedule.

## Reproduction and remaining calculation

[`incoming_timing.py`](incoming_timing.py) reads only the three frozen files listed with SHA-256 hashes in [`incoming-timing.json`](incoming-timing.json). Run `python incoming_timing.py --output /path/to/result.json` with NumPy and SciPy installed. Its 36 checks cover ordered profiles, density-root residuals, cap areas, collision identities and an analytic homogeneous-flow quantile limit.

The next discriminating model needs explicit angular source positions and apertures, packet/grain phase-space distributions, kinetic-to-fluid transition, shock/entropy production and source-to-source errors. Compare recovered nuclear yield at the same prepared energy, then cycle return, source hardware, correction power and scheduling. The existing radial calculations supply a useful target profile and its timing; they do not select the spatial arrangement that can produce it.
