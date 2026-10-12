"""Adapt the reviewed Npp53 offline viewer; embed only verified current artifacts."""
import json,re,hashlib,datetime,os
from pathlib import Path
import numpy as np,pandas as pd
from ete3 import Tree
P=Path(os.environ.get('RB1_PROJECT',str(Path(__file__).resolve().parents[1])))
R=P/'lineage_tracing/tree_viewer_20260923_seed42'
C=P/'lineage_tracing/cassiopeia_eight_20260923_seed42'
F=P/'analysis/tumor_states_20260922/exclude_cl5_cl8_seed42/exclude_current_c17_c19_seed42'
REF=Path(os.environ.get('RB1_VIEWER_TEMPLATE',str(P.parent/'Npp53_tracing/reports/Npp53_tree_viewer_interactive.html')))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def transform_tree(p,sample,cells):
    t=Tree(p.read_text(),format=1);leaves=t.get_leaves()
    assert len(leaves)==len(cells) and set(t.get_leaf_names())==set(cells)
    for leaf in leaves:leaf.name=sample+'|'+leaf.name
    return t.write(format=1)
def replace_function(s,name,following,body):
    a=s.index('function '+name+'(');b=s.index(following,a)
    return s[:a]+body+'\n\n'+s[b:]

def build():
    manifest=json.loads((R/'PREPARED.json').read_text());samples=manifest['samples']
    palettes=json.loads(Path(os.environ.get('RB1_VIEWER_PALETTE',str(F/'SPLIT_PLOTS_COMPLETE.json'))).read_text())['colour_mapping']
    pal={s:'#'+''.join(f'{round(v*255):02x}' for v in rgb) for s,rgb in palettes.items()}
    modes=['state','fitness','pseudotime','treatment','run']
    data=dict(samples=samples,treeTypes=['Greedy','NJ','Lineage Trajectory','Spectral','MaxCut','Percolation'],
        colorModes=modes,colorLabels=dict(state='Putative cell state (overlay only)',fitness='LBI fitness (Greedy)',pseudotime='PhyloVelo pseudotime (Greedy)',treatment='Treatment',run='Library / run'),
        stateOrder=list(pal),trees={},meta={},palettes=dict(state=pal,treatment={'Intact':'#2ca02c','Castration':'#ff7f0e','Castration+ENZA':'#d62728'},run={'MJZ005':'#9467bd','MJZ008':'#17becf','MJZ025':'#bcbd22'}),
        stateToMajor={},cellPlasticity={},plasticityFitness={},sampleInfo={},treatmentColors={'Intact':'#2ca02c','Castration':'#ff7f0e'},carta={},status={},counts={})
    scheduler=json.loads((R/'SCHEDULER_STATUS.json').read_text()) if (R/'SCHEDULER_STATUS.json').exists() else {}
    def checked_status(sample,method,status):
        state=scheduler.get(sample+'/'+method,{}).get('state','')
        if status.startswith(('Pending','Started')) and state in ['FAILED','TIMEOUT','OUT_OF_MEMORY','CANCELLED','NODE_FAIL']:
            return state+' — no verified result'
        return status
    sources={};summary=[]
    for sample in samples:
        cas=R/'subsample_cassiopeia' if sample.startswith('JZ135_n') else C
        meta=pd.read_csv(cas/'inputs'/sample/'cell_metadata.csv',index_col=0)
        data['sampleInfo'][sample]=meta.treatment.iloc[0];data['counts'][sample]=len(meta)
        data['trees'][sample]={};data['status'][sample]={}
        result=R/'results'/sample
        for method in ['Greedy','NJ','Spectral','MaxCut','Percolation']:
            out=cas/'results'/sample if method=='Greedy' else result/method
            qc=out/'COMPLETE.json';nwk=out/(f'{sample}_tree_greedy.nwk' if method=='Greedy' else 'tree.nwk')
            status='Pending / not yet available'
            if qc.exists():
                d=json.loads(qc.read_text());assert d['status']=='COMPLETE' and sha(nwk)==d['tree_sha256']
                data['trees'][sample][method]=transform_tree(nwk,sample,meta.index)
                sources[str(nwk.relative_to(P))]=sha(nwk);status='Verified'
            elif (out/'FAILED.json').exists():status='Failed — review required'
            elif sample.startswith('JZ135_n') and method in ['Spectral','MaxCut','Percolation']:status='Unavailable: ambiguous edits (not resolved)'
            elif (out/'STARTED.json').exists():status='Started; completion not yet verified'
            data['status'][sample][method]=checked_status(sample,method,status)
        fit={};pt={}
        for method,filename,column,destination in [('LBI','lbi_fitness.csv','fitness',fit),('PhyloVelo','pseudotime.csv','pseudotime',pt)]:
            out=result/method;qc=out/'COMPLETE.json';target=out/filename
            status='Pending / not yet available'
            if qc.exists():
                d=json.loads(qc.read_text())
                if method=='PhyloVelo' and not d['pseudotime_valid']:status='Withheld: fewer than 20 distinct rounded values'
                else:
                    df=pd.read_csv(target)
                    if method=='LBI':
                        assert sha(target)==d['fitness_sha256'];df=df[df.is_leaf];df=df.set_index('node')
                    else:df=df.set_index('cell')
                    assert set(df.index)==set(meta.index) and len(df)==len(meta) and np.isfinite(df[column]).all()
                    destination.update(df[column].to_dict());sources[str(target.relative_to(P))]=sha(target);status='Verified'
            elif (out/'FAILED.json').exists():status='Failed — review required'
            elif (out/'STARTED.json').exists():status='Started; completion not yet verified'
            data['status'][sample][method]=checked_status(sample,method,status)
        for cell,row in meta.iterrows():
            m=dict(cell=cell,s=row.putative_state,t=row.treatment,run=row.run_id,sample=sample)
            if cell in fit:m['f']=float(fit[cell])
            if cell in pt:m['pt']=float(pt[cell])
            data['meta'][sample+'|'+cell]=m
        if pt and 'Greedy' in data['trees'][sample]:data['trees'][sample]['Lineage Trajectory']=data['trees'][sample]['Greedy']
        data['status'][sample]['Lineage Trajectory']='Verified: Greedy topology + pseudotime reveal' if pt else data['status'][sample]['PhyloVelo']
        summary.append(dict(analysis=sample,cells=len(meta),**data['status'][sample]))
    html=REF.read_text()
    html=re.sub(r'const DATA = .*;\n','const DATA = '+json.dumps(data,separators=(',',':'),allow_nan=False)+';\n',html,count=1)
    html=html.replace('Npp53 Interactive Tree Viewer','Npp53Rb1 Interactive Tree Viewer').replace('All view (25 samples)','All view (10 datasets)')
    note=('8 biological samples · JZ135: nested 1,000 / 2,000 / 3,000-cell versions, not independent replicates. Full JZ135 held. '
          'Cell states are display overlays and sampling strata only; not used to fit trees, LBI or PhyloVelo. '
          'LBI is a sampling-sensitive relative proxy, not measured growth; Greedy-derived LBI/pseudotime overlays remain the same across alternative trees. Gray = unavailable. CARTA not run. '
          'Snapshot: '+datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%d %H:%M UTC')+'.')
    html=html.replace('<div id="sampleInfoBar"',f'<div class="sub" style="max-width:1500px;line-height:1.5">{note}</div><div id="sampleInfoBar"',1)
    html=html.replace('<h3>Hover info</h3>','<h3>Analysis status</h3><div id="analysisStatus" class="info"></div><h3>Hover info</h3>')
    html=replace_function(html,'dropCyclingSub','// --- Cell count','function dropCyclingSub(root) { return root; }')
    html=html.replace("if (!keys.length) { cellCountCache[sample] = 0; return 0; }","if (!keys.length) return DATA.counts[sample] || 0;")
    html=html.replace('c.length || 1','Number.isFinite(c.length) ? c.length : 1')
    # Parenthesize ternary inside sums; preserve genuine zero-length NJ edges.
    html=html.replace('x + (Number.isFinite(c.length) ? c.length : 1)','x + (Number.isFinite(c.length) ? c.length : 1)')
    html=html.replace("svg.addEventListener('click', () => clearCtx(ctx));","svg.onclick = () => clearCtx(ctx);")
    html=replace_function(html,'isHighlighted','function getColor',"""function isHighlighted(cell, mode) {
      const m=DATA.meta[cell]; if(!m) return false;
      const key=mode==='state'?m.s:mode==='treatment'?m.t:mode==='run'?m.run:null;
      return !highlightSet.size || key===null || highlightSet.has(key);
    }""")
    html=replace_function(html,'getColor','function applyLegendHighlight',"""function getColor(cell,mode) {
      mode=mode||document.getElementById('selColor').value;const m=DATA.meta[cell];if(!m)return '#ddd';
      if(!isHighlighted(cell,mode))return '#e6e6e6';
      if(mode==='fitness'){if(m.f===undefined)return '#ddd';return viridis(1-(m.f-FIT_USER_RANGE.min)/(FIT_USER_RANGE.max-FIT_USER_RANGE.min||1));}
      if(mode==='pseudotime')return m.pt===undefined?'#ddd':viridis(m.pt);
      const key=mode==='state'?m.s:mode==='treatment'?m.t:m.run;return DATA.palettes[mode]?.[key]||'#ddd';
    }""")
    html=html.replace("if (mode === 'fitness') {\n    const d", "if (mode === 'pseudotime') { leg.innerHTML='<b>PhyloVelo pseudotime:</b> 0 <span class=\"legend-gradient\" style=\"transform:scaleX(-1)\"></span> 1 · gray = unavailable; within-sample scale, not chronological time'; return; }\n  if (mode === 'fitness') {\n    const d")
    html=replace_function(html,'showTooltipLeaf','function showTooltipNode',"""function showTooltipLeaf(leaf,e) {
      const m=DATA.meta[leaf.name]||{};const f=m.f===undefined?'unavailable':m.f.toFixed(4);const pt=m.pt===undefined?'unavailable':m.pt.toFixed(4);
      const text=`<b>${m.sample||''}</b><br>${m.cell||leaf.name}<br>${m.s||''}<br>${m.t||''} · ${m.run||''}<br>LBI fitness: ${f}<br>PhyloVelo: ${pt}`;
      tooltip.innerHTML=text;tooltip.style.display='block';tooltip.style.left=Math.min(e.clientX+10,innerWidth-260)+'px';tooltip.style.top=Math.min(e.clientY+10,innerHeight-180)+'px';
      document.getElementById('hoverInfo').innerHTML=text;
    }""")
    html=replace_function(html,'showComposition','// ==================== CARTA PANEL',"""function showComposition(node,sampleFP) {
      const leaves=getLeaves(node),mode=getPrimaryMode();document.getElementById('selInfo').textContent=leaves.length+' cells in selected subclone';
      const comp=document.getElementById('composition');
      if(['fitness','pseudotime'].includes(mode)){const key=mode==='fitness'?'f':'pt',values=leaves.map(l=>DATA.meta[l.name]?.[key]).filter(Number.isFinite);comp.textContent=values.length?'Mean '+mode+': '+(values.reduce((a,b)=>a+b,0)/values.length).toFixed(4):'No verified values';return;}
      const counts={};leaves.forEach(l=>{const m=DATA.meta[l.name];const key=mode==='state'?m.s:mode==='treatment'?m.t:m.run;counts[key]=(counts[key]||0)+1;});
      comp.innerHTML=Object.entries(counts).sort((a,b)=>b[1]-a[1]).map(([k,v])=>`<div class="count-row"><span>${k}</span><span>${v} (${(100*v/leaves.length).toFixed(1)}%)</span></div>`).join('');
    }""")
    html=replace_function(html,'renderCarta','function renderSampleInfoBar',"""function renderCarta(){document.getElementById('cartaContent').textContent='Not run — reserved for CARTA results.';}
    // ==================== MODE MANAGEMENT ====================
    let currentMode='focus';const selectedSamples=[];const treatFilter=new Set();
    function renderAnalysisStatus(s){
      document.getElementById('analysisStatus').innerHTML=Object.entries(DATA.status[s]||{}).map(([k,v])=>`<div style="padding-bottom:5px"><b>${k}</b>: ${v}</div>`).join('');
    }""")
    html=html.replace("const nwk = (DATA.trees[s] || {})[tt];\n  if (!nwk)","renderCarta();renderAnalysisStatus(s);document.getElementById('ptWrap').style.display=tt==='Lineage Trajectory'?'inline-block':'none';document.getElementById('hoverInfo').textContent='Move mouse over a leaf.';document.getElementById('selInfo').textContent='Click a node to select a subclone.';document.getElementById('composition').innerHTML='';\n  const nwk = (DATA.trees[s] || {})[tt];\n  if (!nwk)",1)
    html=html.replace("const ptCutoff = parseFloat(document.getElementById('selPT').value);","const ptCutoff = tt==='Lineage Trajectory'?parseFloat(document.getElementById('selPT').value):100;")
    html=html.replace('function loadFocus() {','function loadFocus() {\n  tooltip.style.display="none"; renderSampleInfoBar();')
    html=html.replace("const availTT = Object.keys(DATA.trees[s] || {}); if (!availTT.length) return;","const availTT = Object.keys(DATA.trees[s] || {});")
    html=html.replace("No tree available","Tree pending / unavailable — see analysis status")
    html=html.replace("if ((mode === 'state' && (k === 'c14_Cycling' || k === 'Cycling')) || (mode === 'major' && k === 'Cycling')) return;",'')
    html=html.replace("syncPrimaryShim(); renderLegend(); rerenderAll();","highlightSet.clear(); syncPrimaryShim(); renderLegend(); rerenderAll();")
    html=html.replace('max-height:60px','max-height:100px')
    # Standalone/offline: reference template has no external script dependencies.
    assert '<script src=' not in html and 'JZ136' not in html
    out=P/'Npp53Rb1_tree_viewer_interactive.html';out.write_text(html)
    pd.DataFrame(summary).to_csv(R/'viewer_analysis_status.csv',index=False)
    receipt=dict(built_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),html_sha256=sha(out),
        reference_sha256=sha(REF),sources=sources,analyses=len(samples),biological_samples=8,
        full_JZ135='HELD by user; cancelled702686_0',carta='NOT RUN',state_independent_fits=True,status=summary)
    (R/'VIEWER_BUILD.json').write_text(json.dumps(receipt,indent=2))
    print(out);print(pd.DataFrame(summary).to_string(index=False))
if __name__=='__main__':build()
