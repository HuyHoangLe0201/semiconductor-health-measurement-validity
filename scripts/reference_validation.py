"""Healthy-reference intervention under stable versus changed machine parameters."""
import json
import numpy as np
import matplotlib.pyplot as plt
import plant_validation as plant

def render_bias(result):
    # A scale-degenerate unrestricted fit has no meaningful plotted health value.
    references=result['healthy_reference']['records']
    shown=[item for item in result['conditions'] if item['name']!='extended_nuisance']
    shown.append({'name':'stationary_healthy_reference','estimate':references[0]['relative_estimate'],'nominal_SE':references[0]['nominal_relative_SE']})
    labels=['Matched','Noisy\nfeedback','1% model\nmismatch','Known L,\nunknown E','Stationary\nreference']
    fig,axes=plt.subplots(1,2,figsize=(7.2,3.1),layout='constrained')
    for j,ax in enumerate(axes):
        ax.errorbar(np.arange(len(shown)),[1000*item['estimate'][j] for item in shown],yerr=[1.96*1000*item['nominal_SE'][j] for item in shown],fmt='o',capsize=3)
        ax.axhline(50 if j==0 else .2,color='gray',ls='--',label='Injected contrast')
        ax.set_xticks(range(len(shown)),labels,fontsize=7)
        ax.set_ylabel('Threshold contrast (mV)' if j==0 else r'Resistance contrast (m$\Omega$)')
        if j==0:
            ax.set_yscale('symlog',linthresh=100)
            ax.set_ylim(30,5000)
            ax.set_yticks([50,100,1000,3000],['50','100','1000','3000'])
        ax.legend(fontsize=7);ax.grid(alpha=.2)
    plant.save(fig,'plant_bias')

def main():
    result=json.loads((plant.DATA/'plant_results.json').read_text())
    wrong=dict(plant.CFG,L=plant.CFG['L']*1.01,E=plant.CFG['E']*1.01)
    references=[]
    for label,cfg in [('stationary_machine',wrong),('changed_machine',plant.CFG)]:
        trace=plant.run(cfg=cfg)
        X,y,mask=plant.observations(trace)
        estimate,se,rank,defect,_,details=plant.fit(X,y,keep=mask)
        references.append({'name':label,'healthy_estimate':estimate.tolist(),'nominal_SE':se.tolist(),'rank':rank,'details':details})
        if label=='stationary_machine':
            np.savez_compressed(plant.DATA/'healthy_reference_trace.npz',current=trace['current'],states=trace['states'],X=X,y=y,retained=mask)
        print(label,estimate,flush=True)
    aged=next(item for item in result['conditions'] if item['name']=='L_and_EMF_mismatch')
    for ref in references:
        ref['relative_estimate']=(np.array(aged['estimate'])-ref['healthy_estimate']).tolist()
        ref['nominal_relative_SE']=np.sqrt(np.array(aged['nominal_SE'])**2+np.array(ref['nominal_SE'])**2).tolist()
    result['healthy_reference']={'scope':'Two independently generated healthy records; both contrasts are differences of fitted health representatives. Stationary mismatch cancels common offset; changed machine violates that stationarity assumption. Nominal SE assumes independent measurement records and excludes discrepancy.','records':references}
    original=np.load(plant.DATA/'plant_trace.npz')
    dv=np.zeros((3,8));dr=dv.copy();dv[0,0]=.05;dr[0,0]=.0002
    fine=plant.run(dv=dv,dr=dr,refine=4,replay=original['states'])
    XX,yy,valid=plant.observations(fine)
    finer,_,_,_,_,_=plant.fit(XX,yy,keep=valid)
    matched=next(item for item in result['conditions'] if item['name']=='matched')
    result['integration_target_difference']=np.abs(finer-np.array(matched['estimate'])).tolist()
    result['integration']={'method':'RK4 with bisection at stage sign changes','hold_subdivisions':2,'blanking_subdivisions':1,'maximum_sign_bisection_depth':10,'convergence_hold_subdivisions':4}
    render_bias(result)
    (plant.DATA/'plant_results.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result['healthy_reference'],indent=2))

if __name__=='__main__':main()
