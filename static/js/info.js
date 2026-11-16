"use strict";
// Structure information: a one-line summary under the file field and a card with cell, space group, composition and data quality.
// Everything comes from the CIF (model.items) and from what this tool computed (prep); missing values are left out or computed.
(function(){
const SUB='₀₁₂₃₄₅₆₇₈₉',sub=s=>String(s).replace(/\d/g,d=>SUB[+d]);
const ok=v=>v!==undefined&&v!==null&&String(v).trim()!==''&&!/^[?.]$/.test(String(v).trim());
const g=(...k)=>{for(const x of k){const v=model.items[x];if(ok(v))return String(v).trim().replace(/\s+/g,' ');}return null;};
const el=document.getElementById.bind(document);
// 'C20 H18 Fe N2 O' -> 'C₂₀H₁₈FeN₂O'
const prettyFormula=s=>s.split(' ').map(t=>{const m=/^([A-Za-z]+)([\d.]*)$/.exec(t);return m?m[1]+(m[2]&&m[2]!=='1'?sub(m[2].replace(/\.0+$/,'')):''):t;}).join('');
// 'P 21 21 21' -> 'P 2₁ 2₁ 2₁', '-3' -> '3̄'
const prettySG=s=>s.replace(/\b([2346])([1-5])\b/g,(m,a,b)=>a+SUB[+b]).replace(/-(\d)/g,'$1\u0304');
const SYSTEM=n=>n<=2?'triclinic':n<=15?'monoclinic':n<=74?'orthorhombic':n<=142?'tetragonal':n<=167?'trigonal':n<=194?'hexagonal':'cubic';
const pct=x=>isFinite(x)?(100*x).toFixed(2)+' %':null;

function rows(){
  const R=[],S=t=>R.push({sec:t}),A=(k,v)=>{if(v!==null&&v!==undefined&&v!=='')R.push({k,v});};
  const c=model.cell,raw=k=>g('_cell_'+k);
  S('Cell');
  A('a, b, c',['length_a','length_b','length_c'].map((k,i)=>raw(k)||[c.a,c.b,c.c][i].toFixed(4)).join('  ')+' Å');
  A('α, β, γ',['angle_alpha','angle_beta','angle_gamma'].map((k,i)=>raw(k)||[c.al,c.be,c.ga][i].toFixed(3)).join('  ')+' °');
  A('Volume',(g('_cell_volume')||model.V.toFixed(2))+' Å³');
  const sgName=g('_space_group_name_h-m_alt','_symmetry_space_group_name_h-m','_space_group_name_h-m'),sgNo=g('_space_group_it_number','_symmetry_int_tables_number');
  S('Space group');
  A('Symbol',sgName?prettySG(sgName)+(sgNo?'  (No. '+sgNo+')':''):(sgNo?'No. '+sgNo:null));
  A('Hall',g('_space_group_name_hall','_symmetry_space_group_name_hall'));
  A('Crystal system',g('_space_group_crystal_system','_symmetry_cell_setting')||(sgNo&&isFinite(+sgNo)?SYSTEM(+sgNo):null));
  const inv=model.ops.some(o=>o.R.every((r,i)=>r.every((x,j)=>Math.abs(x-(i===j?-1:0))<1e-6)));
  A('Operations',model.ops.length+' ('+(inv?'centrosymmetric':'non-centrosymmetric')+')');
  S('Composition');
  const sum=g('_chemical_formula_sum'),Z=g('_cell_formula_units_z');
  A('Formula',sum?prettyFormula(sum):null);A('Formula weight',g('_chemical_formula_weight'));A('Z',Z);
  A('F(000)',g('_exptl_crystal_f_000')||(model.f000cif?String(model.f000cif):null));
  const mw=+g('_chemical_formula_weight'),zz=+Z;
  const dens=g('_exptl_crystal_density_diffrn')||(mw&&zz?(zz*mw/(0.602214*model.V)).toFixed(3):null);A('Density (calc.)',dens?dens+' g/cm³':null);
  const cnt={};let nH=0;for(const a of model.atoms){cnt[a.el]=(cnt[a.el]||0)+a.occ;if(a.el==='H')nH++;}
  const els=Object.keys(cnt),order=els.includes('C')?['C','H',...els.filter(e=>e!=='C'&&e!=='H').sort()].filter(e=>els.includes(e)):els.sort(); // Hill order
  A('Asymmetric unit',order.map(e=>e+(Math.abs(cnt[e]-1)>1e-6?sub((Math.round(cnt[e]*100)/100).toString()):'')).join('')+'  ('+model.atoms.length+' atoms, '+nH+' H)');
  S('Measurement');
  const wl=g('_diffrn_radiation_wavelength'),rt=g('_diffrn_radiation_type');
  A('Temperature',(()=>{const t=g('_diffrn_ambient_temperature','_cell_measurement_temperature');return t?t+' K':null;})());
  A('Radiation',(rt?rt.replace(/\\a/g,'α').replace(/\\b/g,'β'):'')+(wl?(rt?'  ':'')+'λ = '+wl+' Å':'')||null);
  A('Reflections',[g('_diffrn_reflns_number')&&g('_diffrn_reflns_number')+' collected',prep.nRaw+' read, '+prep.U.length+' unique, '+prep.nSel+' with I>2σ'].filter(Boolean).join('; '));
  A('d(min)',Math.min(...prep.U.map(u=>u.d)).toFixed(2)+' Å'+(g('_diffrn_reflns_theta_max')?'  (θmax '+g('_diffrn_reflns_theta_max')+' °)':''));
  A('Rint',[prep.Rint!==undefined?'this tool '+pct(prep.Rint):null,isFinite(model.ref.Rint)?'CIF '+pct(model.ref.Rint):null].filter(Boolean).join('   '));
  S('Refinement');
  A('R1 (I>2σ)',[pct(prep.R1)?'this tool '+pct(prep.R1):null,isFinite(model.ref.R1)?'CIF '+pct(model.ref.R1):null].filter(Boolean).join('   '));
  A('wR2',[pct(prep.wR2)?'this tool '+pct(prep.wR2):null,isFinite(model.ref.wR2)?'CIF '+pct(model.ref.wR2):null].filter(Boolean).join('   '));
  A('GooF / parameters',[g('_refine_ls_goodness_of_fit_ref'),g('_refine_ls_number_parameters')?g('_refine_ls_number_parameters')+' parameters':null].filter(Boolean).join('  /  '));
  return R;
}
function shortLine(){
  const sum=g('_chemical_formula_sum'),sgName=g('_space_group_name_h-m_alt','_symmetry_space_group_name_h-m','_space_group_name_h-m'),Z=g('_cell_formula_units_z');
  return[sum?prettyFormula(sum):null,sgName?prettySG(sgName):null,Z?'Z = '+Z:null].filter(Boolean).join('  ·  ');
}
window.infoReset=()=>{el('strline').textContent='';el('strline').hidden=true;el('strBody').textContent='';el('strCopy').disabled=true;};
window.infoUpdate=()=>{
  if(!model)return;const R=rows(),b=el('strBody');b.textContent='';
  const dl=document.createElement('dl');dl.className='kv';
  for(const r of R){if(r.sec){const d=document.createElement('div');d.className='sec';d.textContent=r.sec;dl.appendChild(d);}else{const t=document.createElement('dt'),v=document.createElement('dd');t.textContent=r.k;v.textContent=r.v;dl.append(t,v);}}
  b.appendChild(dl);const sl=shortLine(),e=el('strline');e.textContent=sl;e.hidden=!sl;
  el('strCopy').disabled=false;el('strCopy')._text=R.map(r=>r.sec?'\n'+r.sec:r.k+': '+r.v).join('\n').trim();
};
el('strCopy').addEventListener('click',async e=>{
  const b=e.currentTarget,t=b._text;if(!t)return;const t0=b.dataset.label||(b.dataset.label=b.textContent),flash=m=>{b.textContent=m;clearTimeout(b._t);b._t=setTimeout(()=>{b.textContent=t0;},1600);};
  try{await navigator.clipboard.writeText(t);flash('Copied');}catch(_){flash('Failed');}
});
window.InfoCore={rows:()=>model?rows():[],shortLine:()=>model?shortLine():'',prettySG,prettyFormula};
})();
