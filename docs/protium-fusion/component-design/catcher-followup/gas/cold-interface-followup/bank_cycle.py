"""Pure integration API for an active and a regenerating He cryosorber bank.
All powers are per receiver area. Constant heat capacity and the specified
cooldown are apparatus scenarios; callers must charge other thermal paths.
"""
import math
from scipy.integrate import quad
SB=5.670374419e-8

def bank_cycle(sink_K=200.,fraction_Carnot=.3,active_area=4.8666935800053074e-5,
               heat_capacity_J_per_active_m2_K=1000.,period_s=3600.,
               warmup_s=300.,cooldown_s=100.,standby_K=None,cold_K=5.,
               electrical_warming=True):
    standby_K=sink_K if standby_K is None else standby_K
    assert cold_K<standby_K<=sink_K and 0<fraction_Carnot<=1
    assert warmup_s+cooldown_s<=period_s
    C=active_area*heat_capacity_J_per_active_m2_K
    hold=period_s-warmup_s-cooldown_s
    gamma=math.log(standby_K/cold_K)/cooldown_s
    def Qrad(T):return active_area*SB*(sink_K**4-T**4)
    def Wrad(T):return Qrad(T)*(sink_K/T-1)/fraction_Carnot
    Wsen=C*(sink_K*math.log(standby_K/cold_K)-(standby_K-cold_K))/(period_s*fraction_Carnot)
    Wramp=quad(lambda t:Wrad(standby_K*math.exp(-gamma*t)),0,cooldown_s,epsabs=1e-12)[0]/period_s
    Whold=Wrad(standby_K)*hold/period_s
    Wheater=C*(standby_K-cold_K)/period_s if electrical_warming else 0.
    old85=C*(sink_K*math.log(85./cold_K)-(85.-cold_K))/(period_s*fraction_Carnot)
    peakQ=Qrad(cold_K)+C*gamma*cold_K
    peakW=peakQ*(sink_K/cold_K-1)/fraction_Carnot
    return {'sink_K':sink_K,'standby_K':standby_K,'fraction_Carnot':fraction_Carnot,
      'active_bank_radiation_W_m2':Wrad(cold_K),'inactive_sensible_cooldown_W_m2':Wsen,
      'inactive_radiation_during_cooldown_W_m2':Wramp,'inactive_standby_radiation_W_m2':Whold,
      'inactive_electrical_warming_W_m2':Wheater,'inactive_total_W_m2':Wsen+Wramp+Whold+Wheater,
      'additional_to_old_85K_sensible_cycle_W_m2':Wsen+Wramp+Whold+Wheater-old85,
      'both_banks_peak_5K_cooling_W_m2':peakQ+Qrad(cold_K),
      'both_banks_peak_electric_W_m2':peakW+Wrad(cold_K),
      'source_hot_to_cold_bypass_increment_W_m2':0.}
