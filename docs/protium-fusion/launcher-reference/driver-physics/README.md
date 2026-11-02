# Driver physics and cold compression ✨

These studies investigate how a converging fuel flow might be delivered and formed from cold material. They retain useful physical requirements, numerical methods, input records, stopped calculations and failure cases. Their original pulse family is often substantially more massive than the current reference, so their reported operating points and gains must not be transferred directly.

- [Driver delivery report](REPORT.md): launch schedules, charging and discharge, source power, parcel survival, focusing, neutralization and first-contact constraints.
- [Cold-compression study](COLD-COMPRESSION.md): cold-flow surrogate, EOS limits, sensitivity to resolution and artificial viscosity, guard-stopped runs, failures, and a 30 m prescribed-hot comparison. The model is not a demonstrated physical driver. The [corrected full-network regression](../../b2-completion/regression.md) supplies the later comparison at this scale.
- [Delivery inputs and results](delivery-screen.json), [cold-driver summary](cold-driver-summary.json), source scripts, validation records and saved NPZ trajectories retain the original calculations. Manifests and hashes identify their original source; they do not by themselves certify the relocated package.
- [Small-pulse transport diagnostic](../../transport/pulse/small-pulse-check.json): the independent opacity analysis for the earlier reduced-network 30 m case. Its historical JSON labels a fictional scenario separately from the calculated physical diagnostics.

The earlier reports use explicitly labeled illustrative and fictional test-range scenarios. Their assumed performance, dates, and efficiencies are not results or adopted operating points. The [current launcher reference](../README.md), [mixture design](../../mixture-design/README.md), and [component design](../../component-design/README.md) contain the newer plant context.

The archived documents and input hashes are preserved for research continuity. Evaluate each physical model against its actual domain, compare it with later work, and rerun it with appropriate EOS, geometry and radiation transport before applying it to a new design.
