# Interpreting the requested 80% intervals ✨

The acceptance criterion is at least 90% of the corresponding smooth pulse's
integrated nuclear output. Timing error means a mass-weighted deviation from
the intended arrival schedule, after subtracting a common clock offset. A
packet is an independently aimed angular portion of the incoming pulse, with
a shaped outline/density/longitudinal profile and one bulk flight direction.
Its microscopic grains can be much more numerous.

For a fixed error pattern, define the timing threshold as the first downward
crossing of 90% retained output as RMS error increases. Small perturbations can
occasionally increase output. A fitted quadratic or a monotone interpolant must
not erase that behavior or manufacture a crossing beyond the tested range.
Unbracketed crossings remain bounds until another calculation brackets them.

For packet count, distinguish the meridional-sector result, its geometrically
equivalent number of spherical patches, the actual 3D Voronoi-packet calculations,
and the engineering estimate including earlier assembly. Angular direction
variance supplies the 2D geometric conversion. The 3D runs explicitly evolve
edges and junctions at their resolved scales. The thermodynamic history of cold
packet assembly remains the largest extrapolation.

The numerical calculations are deterministic. The requested credible intervals
express judgment about the threshold in the physical system, informed by these
calculations. The report will expose the numerical anchor separately from the
allowances used for unresolved physics. Sensitivity to those allowances matters
more than many digits in a Monte Carlo quantile.

The uncertainty assessment must include:

- Angular spectrum and correlation of timing errors. Coherent low-order errors
  and independent packet errors with the same RMS can have different effects.
- Resolution and hydrodynamic-flux sensitivity, using matched successful
  controls. An aborted calculation supplies no completed tolerance point.
- The construction and time of the perturbation. The first evolved-reference
  seeds omit angular interactions before 11 microseconds. The earlier-seeded
  branch begins at the prepared starting time, using a smooth prehistory checked
  by forward replay. Ballistic seeds preserve parcel energies and approximate
  their thermal history. All these branches begin after substantial preparation.
- Three-dimensional flow and turbulent mixing. Axisymmetric rings have fewer
  degrees of freedom than spherical packets, and the literature supplies no
  universal multiplicative correction.
- Early merger and compression. The simulation starts from a prepared hot
  inflow. Shaped packet edges and graded contact can reduce merger shocks, but
  this study does not yet predict the achieved entropy from a cold launcher.
  This uncertainty can raise the packet requirement by much more than the
  numerical resolution error.

For timing, the principal 3D error family uses Gaussian quantiles permuted among
128 angular regions, normalized by their solid-angle weights. The interval
concerns the tolerable standard deviation of residual timing errors, after the
intended longitudinal schedule and a common clock shift are removed. The
credible interval expresses uncertainty about that threshold; it is separate
from the distribution of individual packet arrival times. A fixed RMS can
describe very different outliers and correlation patterns, so the error family
is part of the result.

The engineering count assessment represents uncertainty about how late a
feasible shaped-packet design can defer substantial irreversible merger heat.
The effective linear compression between merger heating and the prepared initial
profile is a judgment input. Compression within the simulated pulse is already
included in the radial heat sensitivities and is not counted twice. The
radial heat sensitivities constrain the allowed heat at the prepared state;
adiabatic amplification supplies a conditional translation to earlier heating.
The distribution assigned to achievable compression and the additional model
factors are exposed in `subjective_intervals.py` and the final assessment JSON.
Their independence assumptions and prior sensitivity are part of the assessment.

A numerical threshold from the hot-flow model and a subjective engineering
interval should therefore both be visible. If the packet interval depends
strongly on a broad prior allowance for early merger, state that dependence
plainly. The interval is conditional on an architecture capable of producing
the intended smooth compression history. Failure of that architecture has no
finite timing or packet-count threshold and requires a separate feasibility
assessment.

Each isolated threshold uses the whole 10% output-loss budget. A combined
operating point needs smaller timing errors, more packets, or a directly tested
joint tolerance curve. Any additive loss-budget illustration will be labeled
as an approximation until joint simulations test it.

The proposed compression prior has median 4 and untruncated 10th/90th
percentiles approximately 1.07 and 15, truncated below 1. For an outer-radius
scale of 171 m this roughly spans major merger heating at 180–2,560 m; the
existing launcher note places a spherical packing-density backmap near
1,755 m. That backmap is context, not a measured collision location. The prior
expresses the judgment that shaped packets could delay substantial contact,
with considerable uncertainty. Changing its median to 2 or 8 exposes how
strongly this engineering assumption controls the inferred packet count.

The early-heat screen measures sensitivity of compression and burning to
entropy. It supplies no independent demonstration of instability growth.
Instability evidence comes from the multidimensional disturbed-flow histories.
The final output criterion includes both effects because either can consume
the allowed 10% loss.
