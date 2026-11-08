"""Expose analyst judgments separately from deterministic simulation thresholds.

The draws express uncertainty about physical thresholds. They are not repeated
hydrodynamic experiments or a likelihood fitted to measured hardware.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.special import ndtri

p=argparse.ArgumentParser()
p.add_argument('--timing-anchor-us',type=float,required=True)
p.add_argument('--hot-packet-anchor',type=float,required=True)
p.add_argument('--out',type=Path,required=True)
a=p.parse_args();z90=float(ndtri(.9));rng=np.random.default_rng(20261031);n=400000

# These distributions are deliberately explicit analyst judgments. The effective
# compression after excess merger heat is the dominant packet-count uncertainty.
normal=rng.normal(size=(7,n))
exponent=rng.uniform(1.5,2,n)
heat_fraction=np.exp(rng.uniform(np.log(.026),np.log(.073),n))
hot_count=a.hot_packet_anchor*np.exp(.30*normal[0])
extra_count_factor=np.exp(np.log(2)/z90*normal[1])
timing_model_factor=np.exp(.25*normal[2])
timing_transfer_factor=.6*np.exp(.9*normal[3])

def summary(x):
    q=np.quantile(x,[.1,.5,.9])
    return dict(median=float(q[1]),interval80=[float(q[0]),float(q[2])])

counts={}
for median in [2.,4.,8.]:
    compression=np.maximum(1.,median*np.exp(np.log(15/4)/z90*normal[4]))
    heat_requirement=2*compression**exponent/heat_fraction
    required=np.maximum(hot_count,heat_requirement)*extra_count_factor
    counts[str(median)]=dict(**summary(required),probability_minimum_at_most_1000=float(np.mean(required<=1000)))

result=dict(
    meaning='Subjective, conditional engineering assessment; numerical anchors plus explicit judgment distributions',
    acceptance='At least 90% of matched smooth-pulse integrated nuclear output; each isolated tolerance uses the whole 10% loss budget',
    scope='Prepared 30 m-reference Stage A nitrogen/proton pulse; one bulk direction per shaped angular packet; adequate aiming; architecture can produce the intended hot compression history',
    timing_error_family='Light-tailed arrival errors with no common offset, at angular scales represented by the tested patterns; no adversarial or arbitrary-outlier RMS claim',
    anchors=dict(timing_us=a.timing_anchor_us,hot_packet_count=a.hot_packet_anchor,added_heat_fraction90=.026,credited_heat_fraction90=.073),
    judgments=dict(
        effective_compression_before_prepared_profile='max(1, lognormal with median 4 and untruncated 90th percentile 15)); sensitivity medians 2 and 8',
        thermal_amplification_exponent='uniform between 1.5 and 2',
        allowed_hot_heat_fraction='log-uniform between 0.026 and 0.073, spanning additional heating and local credit for planned preheat',
        hot_packet_model_factor='lognormal, median 1, log standard deviation 0.30',
        additional_packet_model_factor='lognormal, median 1, 80% interval [0.5,2], for geometry/optimization, unresolved mixing and transfer assumptions',
        timing_model_factor='lognormal, median 1, log standard deviation 0.25',
        timing_transfer_factor='lognormal, median 0.6, log standard deviation 0.9, for unmodeled cold assembly, error-pattern variation and physical/model uncertainty'),
    arrival_time_us=summary(a.timing_anchor_us*timing_model_factor*timing_transfer_factor),
    minimum_packet_count=counts['4.0'],
    packet_count_prior_sensitivity=counts,
    sample_count=n,seed=20261031,
    caveat='These are marginal intervals conditional on architecture feasibility. A failed preparation architecture has no finite acceptable threshold and is not assigned a probability by this calculation.')
a.out.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ['arrival_time_us','minimum_packet_count','packet_count_prior_sensitivity']},indent=2))
