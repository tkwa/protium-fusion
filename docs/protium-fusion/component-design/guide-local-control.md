# ✨ Local guide dynamics checkpoint

A passive flexible rail is insufficient in the present moving-stream model. A representative unstable bend has wavelength 45.7 km and grows by a factor e in 38.9 seconds. The earlier global affine-mode check does not include this mode. This result makes distributed stream/rail feedback a component requirement; it does not establish that the shell requires a much heavier rigid structure.

The saved calculation is [guide_local_control.py](guide_local_control.py), with [results](guide-local-control.json). It uses the heaviest of the three support families, its paired 10 km/s streams, the full receiver tributary mass, transverse confinement frequency 1.403 rad/s, and a selected weak mesh pretension of 10,000 N/m. The local approximation treats a rail as straight and continuously loaded. Finite carriers, three-family intersections, the complete shell's elasticity and the detailed magnetic field remain outside this model.

For stream centroids y+ and y− and guide displacement Y, the equations are

```
(∂t ± v ∂s)² y± = −κ²(y± − Y)
Ytt + ωg²Y = ε κ²[(y+ + y−)/2 − Y].
```

Here ε is combined stream mass divided by receiver tributary mass. With a Fourier mode exp(st + ikx), b=(kv)², c=κ² and x=s², the characteristic polynomial is

```
[x² + 2(c+b)x + (c−b)²] [x + ωg² + εc]
    − εc²[x + c−b] = 0.
```

When guide stiffness is zero, its constant term is negative for 0 < kv < κ, proving a growing mode. The saved matrix calculation and independent polynomial agree to 2.5×10⁻¹⁵ after scaling. A separate agent independently verified the polynomial. Similar moving-material/structure instabilities are established in fluid-conveying pipes, including experiments; the reactor's particular numerical result is our own reduced-model calculation. Relevant primary examples are [Païdoussis's theory and experiments](https://journals.sagepub.com/doi/10.1243/JMES_JOUR_1970_012_017_02), [electromagnetic feedback calculations](https://doi.org/10.1016/j.jsv.2018.06.045), and an [experimental semi-active stabilization study](https://www.ippt.pan.pl/repository/open/o6415.pdf).

An ideal magnetic axis fixed to the planned trajectory removes the rail-following coupling. Displacement cancellation alone is too brittle to use as the design: finite response and small gain errors can leave a growing mode. The proposed control must observe the two stream positions and transverse velocities as well as guide motion, then command independent small dipole corrections in the two bores. Each actuator's equal-and-opposite reaction on the guide is included; no external stationary foundation supplies a free force.

A full-state linear controller stabilizes the examined modes. The script solves the complex continuous algebraic Riccati equation separately at each wavelength, then adds an explicit first-order actuator response. In the relevant 20–1000 km band, all 2,025 tested combinations remain stable with a 0.1-second actuator time, including ±20% confinement-frequency variation, a factor-two stream-mass range, and ±20% actuator gain. A one-second response fails some cases. These are deterministic sensitivity cases rather than a statistical confidence interval. The real controller must realize the spatial filters with finite sensing and communication; the mode-dependent calculation is an existence and scale check, not that implementation.

For selected 40–1000 km modes, initial guide or stream displacements up to 1 mm require at most 0.0093 m/s² corrective acceleration in the saved transient samples. Equivalent dipole trim stays below 0.038 T, and maximum stream-to-guide displacement stays below 4 mm. Both fit separately within the reference field and 29.6 mm geometric margin. Nonlinear superposition, noise, saturation and coil losses still require analysis. Applying the same aggressive controller to a 1 km wave creates unnecessarily large trim demands; those short waves are already stable in the passive screen. That arbitrary control choice must not become a new submillimetre alignment requirement for the whole station.

The useful next calculation is a real-space, finite-carrier guide segment with shared cross-bracing, sampled sensors, finite latency, trim saturation and localized field faults. It should determine the controller's required spatial bandwidth and ordinary AC loss. These belong inside the proposed guide package and auxiliary budget until a component model demonstrates that they exceed those allocations. Predictable kilometre-scale common breathing of the whole sphere follows the planned trajectory and does not consume the local guide clearance.
