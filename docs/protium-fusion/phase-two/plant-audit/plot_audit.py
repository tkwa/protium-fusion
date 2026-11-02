"""Static figures from saved audit outputs; no recomputation of plant cases."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parent
c=json.loads((P/'capture-competition.json').read_text());h=json.loads((P/'heat-transport.json').read_text());w=json.loads((P/'warm-distribution.json').read_text())
fig,ax=plt.subplots(1,2,figsize=(12.8,5.5),layout='constrained')
colors={'C':'#347956','Zr':'#258cad','Ti':'#d77b29','Ti30Zr30V40':'#7a5ba7','Cu':'#aa4755'}
for coat in colors:
 rows=[x for x in c['cases'] if x['geometry']['coating']==coat and x['geometry']['thickness_m']>0]
 ax[0].plot([x['geometry']['thickness_m']*1e6 for x in rows],[x['plant']['net_W_kg']*1000 for x in rows],marker='o',linestyle='none',color=colors[coat],label=coat)
ax[0].set_xscale('log');ax[0].set_ylim(-.25,1.6);ax[0].axhline(1,color='#333',ls='--',lw=1,label='1 mW/kg target')
ax[0].set_xlabel('Coating thickness (μm)');ax[0].set_ylabel('Conditional net specific power (mW/kg)');ax[0].set_title('Capture competition at actual channel-wall area');ax[0].legend(fontsize=9,ncol=2);ax[0].grid(alpha=.15)
central=next(x for x in h['central_loops'] if x['relative_Carnot']==.85 and x['span_fraction']==.5)
hot=next(x for x in h['alternative_exchanger_loops'] if x['relative_Carnot']==.85 and x['exchanger']['engine_hot_K']==290)
series=[('Central, minimum UA',central['pipe_and_fluid_mass_kg'],central['pump_electric_W'],'#4c6a9a'),('Central, hot side 290 K',hot['pipe_and_fluid_mass_kg'],hot['pump_electric_W'],'#243f69')]
for gas,col in [('H2','#bf6e2c'),('He4','#7a5ba7')]:
 x=w['equal_quality_best']['0.85_'+gas];series.append(('Warm '+('H₂' if gas=='H2' else 'He'),x['total_transport_mass_kg'],x['total_pump_electric_W']+x['lost_reference_generation_from_warm_radiation_W'],col))
for name,m,p,col in series:
 ax[1].scatter(m/1e14,p/1e9,color=col,s=85);ax[1].annotate(name,(m/1e14,p/1e9),xytext=(7,-12 if 'minimum' in name else 8),textcoords='offset points',fontsize=9)
ax[1].axvline(h['reference']['radiator_mass_kg']/1e14,color='#777',ls=':',lw=1);ax[1].text(h['reference']['radiator_mass_kg']/1e14-.035,48,'Existing radiator mass allowance',rotation=90,va='bottom',ha='right',fontsize=9,color='#666')
ax[1].set_xlim(.75,2.25);ax[1].set_ylim(30,140);ax[1].set_xlabel('Transport + insulation mass (10¹⁴ kg)');ax[1].set_ylabel('Pumps + radiation generation loss (GW)');ax[1].set_title('Matched engine quality: 85% of local Carnot');ax[1].grid(alpha=.15)
fig.suptitle('Full-plant requirements at the saved 140 T operating point',fontsize=15)
fig.text(.5,-.025,'Assigned material/engine properties; common thermal-neutron flux; explicit transport layouts. Exchanger and machinery mass remain unallocated.',ha='center',fontsize=9)
fig.savefig(P/'plant-audit.png',dpi=180,bbox_inches='tight');fig.savefig(P/'plant-audit.svg',bbox_inches='tight')
