# Explicit source-profile integration ✨

The selected 64-zone source-profile family gives **6.568906 mW/kg**, compared with **6.572770 mW/kg** for the reference 003 recipes: a **0.0588% decrease** under the same conditional plant assumptions. The stage outputs and preparation energies are integrated through the actual isotope circulation, receiver reconstruction and plant mass/power ledger. No proportional pulse-gain derating is used.

| Quantity | Reference 003 | Selected source profiles |
|---|---:|---:|
| Prepared energy, MeV/net He | 1.960452 | 1.964126 |
| Net export, W | 8.01246×10²⁰ | 8.00776×10²⁰ |
| Plant mass, kg | 1.2190393×10²³ | 1.2190395×10²³ |
| Capacity net He/s | 4.473376×10³² | 4.472180×10³² |
| Main-window occupancy | 0.789876 | 0.789725 |
| Radiator area fraction | 0.900000 | 0.900000 |
| Fresh C:H mass ratio | 0.00035863621 | 0.00035863631 |

The unchanged comparison assumptions are radius 1.1×10⁹ m, 200 K radiator, 60% thermal conversion, 60% preparation efficiency, 99% energy capture and complete eventual ordinary-material return. The result retains the existing conditional source, heat-exchanger, receiver, services and support mass model. Receiver surface-temperature estimates remain material/spectrum-dependent diagnostics.

## Exact stage selection

The selected map is [selected-source-paths.json](../component-design/selected-source-paths.json), with profiles, trace identities and input hashes in [selected-source-stage-map.json](../component-design/selected-source-stage-map.json).

- **A-H3:** mapped local merger-entropy profile with controlled initial preheat to the reference temperature floor. The added initial internal energy is paid from radial kinetic energy, preserving total prepared energy. The floor is part of the prepared input profile.
- **C13-H4:** mapped local merger-entropy temperature profile, with total prepared energy preserved by the radial velocity adjustment.
- **B1-H2 and B2-H2.5:** original 003 profiles and completed full-network outputs.

These are prescribed spherical hot-flow sensitivities. Their equation of state during burning remains the archived FastEOS; FEOS supplies the preparation-profile mapping. The multidimensional preparation path is a separate engineering qualification.

The C13 profile changes its cold N14 yield from **0.222027 to 0.252914 per initial catalyst**, while N15 changes from **0.713202 to 0.670622**. Its retained C13 fraction also rises. The complete stationary balance consequently raises C13 attempts/net He only from **0.04002785 to 0.04053464**. A attempts rise from **1.266157 to 1.272424**; B1 and B2 attempts are essentially unchanged. The resulting prepared-energy increase is **0.1874%**, and gross fuel circulation per He rises **0.1902%**. These coupled routing changes explain why the plant response differs from a scalar pulse-gain correction.

At the original helium-production rate, net export decreases **0.0317%**. Respecting the 90% radiator boundary reduces the capacity rate by a further **0.0267%**. Peak installed source power remains controlled by B2 and is unchanged. All four receiver energy reconstructions pass, with no reported thermal-grid energy warnings.

## Utility and replay

[variant_integrate.py](variant_integrate.py) accepts exactly four named stage records, either as a simple map or as a `stage_map` member. Each map entry can specify `record`, `initial`, `inputs`, `trace`, `duration_s`, `initial_Y`, `eos` and expected file hashes. Relative paths are resolved against the supplied research root. An optional provenance document provides the explicitly named sidecars and trace paths. The utility never discovers alternative runs.

Original records are preserved. Temporary metadata adapters provide `trace_path` and the declared physical duration to the existing reducer. Actual initial composition is checked against the trace, and the original initial/internal/kinetic/radiation energy fields are independently reconstructed. Failed, guarded or short trajectories, stage-label mismatches, conflicting hashes and unsupported EOS declarations are rejected. Nonideal EOS variants require a matching energy-reconstruction adapter.

The result embeds the actual stage data and complete receiver reductions. It can therefore be replayed without the large traces:

```bash
python "$PROTIUM_RESEARCH/mixture-design/variant_integrate.py" \
  --research-root "$PROTIUM_RESEARCH" \
  --replay "$PROTIUM_RESEARCH/mixture-design/variant_source64.json" \
  --family selected64-source-response-N2000-replay \
  --work-dir /tmp/protium-variant-replay \
  --output /tmp/variant_source64_replay.json
```

For a new mesh or rate family, supply its four explicit records and matching initial profiles, input sidecars and traces. The current receiver projection uses a declared 10⁹ m reference radius before the plant scales it. Updated nuclear-rate families have their own input provenance and should be integrated as separate named variants.

## Verification and artifacts

The explicit reference map reproduces the 003 plant and cycle quantities exactly. The selected source integration passes **18 checks**; replay reproduces all physical cycle fields and the whole plant result. The atomic-mass subtraction residual changes only at approximately 10⁻¹³ MeV/net He between processes, remaining below 10⁻¹⁰. Eleven additional verification checks cover reference reproduction, replay and rejection of missing stages, wrong labels, mismatched hashes, unsupported EOS and unfinished physical durations.

- [Selected result](variant_source64.json) and [four receiver reductions](variant_source64_receivers/n14.receiver.json).
- [Reference reproduction](variant_reference003.json) and [explicit reference map](variant_reference003_paths.json).
- [Replay result](variant_source64_replay.json), [verification results](variant_validation.json) and [verification script](variant_verify.py).

The JSON files retain consumed implementation, parameter, profile and trace hashes. No archived material packet, original hydro record or plant comparison implementation was changed.
