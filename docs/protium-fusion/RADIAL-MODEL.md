# Spatial pulse model ✨

`radial_pulse.py` tests the initially compressed state with spatially varying burning and expansion. It does not supply that state or demonstrate a driver.

## Equations and numerical method

Fixed-mass spherical cells carry composition and temperature. Their boundaries carry radius and velocity. Density is cell mass divided by its spherical volume. A node's mass is half the adjacent cell mass on each side, with half a cell at the free outer surface. The center is fixed by symmetry.

The force on node j is `4*pi*r_j^2*(P_left-P_right)`. Each cell loses internal energy at `(P+q)*dV/dt`, where q is a standard compression-only artificial viscosity. The same P+q is used in both expressions. Summing over nodes and cells cancels their mechanical energy exchange in the continuous-time discretization. The BDF time integrator contributes finite tolerance error, which is reported against initial internal energy. Nuclear mass defects, escaped weak-decay energy and radiative escape enter the same ledger.

The viscosity is `rho*[2*(negative velocity jump)^2 + 0.5*c_s*(negative velocity jump)]`. Its effect on burning requires resolution and coefficient checks. Spatial symmetry suppresses multidimensional instabilities and mixing, so these runs cannot establish an achievable implosion.

The EOS includes ideal ions, radiation, and an interpolated arbitrary-degeneracy/relativity electron–positron contribution from pynucastro. Its composition derivatives are included when evolving temperature. The reaction network has 27 isotopes and 194 rates. Screening, Coulomb corrections and electron heat conduction are currently omitted.

Optional grey radiation diffusion uses the rational Levermore–Pomraning limiter `(2+R)/(6+3R+R^2)`. The outer boundary approaches blackbody escape when optically thin and diffusion when thick. Opacity is a fixed input, presently 0.001 m²/kg. This is a transport sensitivity parameter, not a calculated opacity table. LTE radiation is already included in the EOS, so its energy and pressure are not added a second time.

References for the method, not validation of this implementation:

- Compatible Lagrangian force/work construction: <https://doi.org/10.1006/jcph.1998.6029>.
- A published implementation of the same flux-limiter form: <https://doi.org/10.1051/0004-6361/201424954>.

## Verification

For a gamma=5/3 ideal gas, uniform initial density and initial pressure `P(r)=P_c*(1-r²/R_0²)`, homologous expansion has the exact scale factor `a(t)=sqrt(1+K*t²)`, with `K=2*P_c/(rho_0*R_0²)`. In initial central sound-crossing units, `a=sqrt(1+1.2*t²)`. Density falls as a^-3, temperature as a^-2, and velocity is proportional to initial radius. The test initializes cell-volume-averaged pressure.

At one sound-crossing time, 16, 32 and 64 cells give mass-weighted temperature errors of approximately 0.00207, 0.000547 and 0.000143. Their total-energy residuals are approximately 1.7e-6, 4.9e-7 and 2.8e-7 of initial energy. This establishes useful convergence for a smooth nonreacting solution. Reacting fronts require separate convergence tests.

An initially uniform nonreacting sphere also shows expansion starting at the surface while the center initially remains undisturbed. The 32-cell result at 0.3 sound-crossing times conserves energy to about 2.6e-8 of its initial value.

The default SciPy adaptive finite-difference Jacobian overflowed its perturbation factor for absent species in the stiff network. The code now supplies fixed-scale grouped finite differences. This is a numerical repair; it does not modify the nuclear rates.

## Preliminary reacting results

The 100 m radius N14+p case at rho=1e8 kg/m3 and 300 MK initially gives desired-product fractions 0.733, 0.786 and 0.790 at 16, 32 and 64 cells. The corresponding homogeneous calculation gave 0.897. Surface expansion therefore materially reduces burning at this scale.

At 300 m, radiation-enabled full-network results give conversion fractions 0.9351, 0.9038, 0.9169 and 0.9189 at 32, 64, 128 and 256 cells. A reduced N14+p capture/reverse network reproduces the full-network conversions at 64, 128 and 256 cells within 2e-6 absolute. It gives 0.91972 and 0.92030 at 512 and 1024 cells, a relative change of 0.063%. This is useful convergence of bulk conversion; the very small radiative boundary flux is still resolution dependent. At 1024 cells, escaped radiation is 5.2e-7 of the nuclear energy. The fuel mass is approximately 1.13e16 kg, while complete plant mass is still uncalculated.

The corresponding B2 carbon pulse at 300 m gives 0.9011 and 0.9435 conversion at 32 and 64 cells, so its spatial convergence remains open. The lower-density B1 case (rho=1e6 kg/m3, T=200 MK, R=1500 m, 64 cells) gives 0.9917 desired conversion and nuclear energy 10.52 times initial internal energy. Its lower peak temperature, 443 MK, suppresses unwanted reactions. Rare-isotope yields below integration tolerances need stricter verification before entering a loss ledger.

The EOS table was extended to electron-density coordinates 1e-6 through 1e7 g/cm3 and temperatures 1e5 through 3e9 K for expanded surface zones. A new 37-point comparison against the original EOS gives maximum relative errors 1.6e-6 in pressure, 9.9e-6 in energy, 1.6e-4 in heat capacity and 6.4e-5 in adiabatic index. These compare the interpolant with its fully ionized source EOS; they do not validate that source approximation at cold, partially ionized states.

These results support continued investigation of opaque pulses. They do not prove a complete reactor's net specific power, stability or economics.

## Convergent-flow experiment

`compression_probe.py` obtains a cool spatial profile by expanding a nonreacting reference sphere for two sound-crossing times, reverses its bulk velocity and follows the return with reactions enabled. Both initial internal energy and inward kinetic energy count as input. This is a controlled initial-condition experiment; it does not supply a way to prepare the flow. Its initial temperature is spatially variable, approximately 0.9–38 MK in the 256-cell example, and its outer radius is approximately 1.84 km.

At 64, 128 and 256 cells the first burning runs give nitrogen-stage conversion 0.786, 0.811 and 0.825. The 256-cell nuclear release is 12.53 times initial thermal-plus-kinetic energy. Decreasing the linear numerical-viscosity coefficient from 0.5 to 0.25 gives 0.827 conversion and ratio 12.56. The 512-cell result gives 0.8313 conversion and ratio 12.63; a full-network cross-check is underway. These results show that early burning changes the compressed-state screen without immediately extinguishing the candidate pulse.

The nonreacting return approximately recovers the original thermal state. Artificial viscosity creates entropy during compression, and the free-surface tail is imperfectly reversible at finite resolution. Its total energy is conserved to about 1e-8 relative in the sampled returns. Reacting energy residuals are a few parts in 1e6 of input energy.

The supplied Jacobian stencil was widened for radiation: a face's density depends on both bounding nodes, so the inner radiative flux reaches two cells left. The previous stencil was an incomplete approximate Jacobian, although the ODE itself already included that dependence. Successful earlier integrations still satisfy their reported conservation and tolerance checks; future integrations use the complete stencil.
