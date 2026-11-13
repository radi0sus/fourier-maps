// Simulated refinement (teaching mode). Nothing is refined: positions and ADPs come from the CIF.
// Only the atoms ticked in the model list contribute to Fc. The atoms that are left out are still in Fo,
// so they show up as peaks in Fo-Fc. Adding atoms step by step makes the peaks disappear and R1 fall.
// Fo keeps the scale of the full model (K fixed), unless "Refit scale" is on. The phases of the Fo map come
// from the partial model, as in a real structure solution (model bias included).
(function(){
const SIM={on:false,sel:new Set(),hist:[],cache:null,ghost:true,refit:false,auto:true};
const ZS='H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi'.split(' ');
const Zof=el=>ZS.indexOf(el)+1;
const P1=v=>(100*v).toFixed(1)+' %';

// Fc, phases, scale and R values for the ticked atoms (cached per selection)
function simCalc(){
  const key=[...SIM.sel].sort().join('|')+(SIM.refit?'#r':'');
  if(SIM.cache&&SIM.cache.key===key)return SIM.cache;
  const sel=model.atoms.filter(a=>SIM.sel.has(a.lab)),m2=Object.assign({},model,{atoms:sel});
  const U=prep.U.map(u=>{const r=sel.length?calcF(m2,u.h):[0,0],f2=r[0]*r[0]+r[1]*r[1];
    return Object.assign({},u,{Fc2:f2,Fc:Math.sqrt(f2),phi:Math.atan2(r[1],r[0])});});
  const wt=(u,G)=>{const P=(Math.max(u.I*G,0)+2*u.Fc2)/3,s=u.s*G;return 1/(s*s+(model.wa*P)**2+model.wb*P);};
  let G=prep.K*prep.K;
  if(SIM.refit&&sel.length){let a=0,b=0;for(const u of U){a+=u.I*u.Fc2;b+=u.I*u.I;}G=a/b;
    for(let it=0;it<6;it++){a=b=0;for(const u of U){const w=wt(u,G);a+=w*u.I*u.Fc2;b+=w*u.I*u.I;}if(b>0)G=a/b;}}
  let nw=0,dw=0,rn=0,rd=0;
  for(const u of U){const f2=u.I*G;nw+=wt(u,G)*(f2-u.Fc2)**2;dw+=wt(u,G)*f2*f2;
    if(f2>2*u.s*G){const f=Math.sqrt(f2);rn+=Math.abs(f-u.Fc);rd+=f;}}
  return SIM.cache={key,U,K:Math.sqrt(G),n:sel.length,R1:rd>0?rn/rd:NaN,wR2:dw>0?Math.sqrt(nw/dw):NaN};
}
// Fo-Fc density at a fractional position (trilinear, from the 3D grid that is already on screen)
function densAt(f){
  const[N1,N2,N3]=vol.N,p=[f[0]*N1,f[1]*N2,f[2]*N3],b=p.map(Math.floor),md=(a,n)=>((a%n)+n)%n,NN=N1*N2;
  const x0=md(b[0],N1),x1=md(b[0]+1,N1),y0=md(b[1],N2)*N1,y1=md(b[1]+1,N2)*N1,z0=md(b[2],N3)*NN,z1=md(b[2]+1,N3)*NN,fx=p[0]-b[0],fy=p[1]-b[1],fz=p[2]-b[2];
  return tri(vol.fo,x0,x1,y0,y1,z0,z1,fx,fy,fz)-tri(vol.fc,x0,x1,y0,y1,z0,z1,fx,fy,fz);
}
// atoms that are not in the model, strongest Fo-Fc peak first
function missingPeaks(){
  if(!vol)return[];
  return model.atoms.filter(a=>!SIM.sel.has(a.lab)).map(a=>({lab:a.lab,el:a.el,h:densAt(a.x)})).sort((p,q)=>q.h-p.h);
}

// ---- hooks used by app.js / map.js ----
window.simPrep=p=>{if(!SIM.on)return p;const c=simCalc();return Object.assign({},p,{U:c.U,K:c.K});};
window.simMissing=lab=>SIM.on&&!SIM.sel.has(lab);
window.simGhosts=()=>SIM.on&&SIM.ghost;
window.simInfo=()=>{
  if(!SIM.on||!vol||!model)return'';
  const c=simCalc();let s=`\nPartial model: ${c.n} of ${model.atoms.length} atoms, R1 = ${P1(c.R1)}, wR2 = ${P1(c.wR2)} (full model: ${P1(prep.R1)}, ${P1(prep.wR2)})${SIM.refit?', scale refitted':', scale of the full model'}`;
  const mp=missingPeaks().slice(0,4);
  if(mp.length)s+='\nStrongest Fo−Fc at missing atoms: '+mp.map(q=>`${q.lab} ${q.h.toFixed(2)}`).join(', ')+' e/Å³';
  else s+='\nAll atoms are in the model.';
  const h=SIM.hist.slice(-8),o=SIM.hist.length-h.length;
  s+='\n'+h.map((e,i)=>`  ${o+i+1}. ${e.desc} → R1 ${P1(e.R1)}`).join('\n');
  return s;
};

// ---- panel ----
function ui(){
  const c=simCalc(),N=model.atoms.length;
  $('simStat').innerHTML=`<span class="big">R1 ${P1(c.R1)}</span><span class="big">wR2 ${P1(c.wR2)}</span><span class="hint">full model ${P1(prep.R1)} / ${P1(prep.wR2)} &middot; ${c.n} of ${N} atoms</span>`;
  $('simAtoms').querySelectorAll('button').forEach(b=>b.classList.toggle('sel',SIM.sel.has(b.dataset.lab)));
  $('simAdd').disabled=SIM.sel.size>=N;$('simUndo').disabled=SIM.hist.length<2;
}
// contour step that follows the model: about 12 lines up to the strongest missing peak, snapped to a round value
function autoStep(){
  if(!SIM.auto||$('type').value!=='diff'||!vol)return;
  let ref=missingPeaks()[0];ref=ref?ref.h:0;
  if(!(ref>0.02)){const A=vol.fo,B=vol.fc;ref=0;for(let i=0;i<A.length;i++){const d=Math.abs(A[i]-B[i]);if(d>ref)ref=d;}}
  const want=Math.max(ref,0.1)/12,nice=[0.01,0.02,0.03,0.05,0.1,0.2,0.3,0.5,1];
  let best=nice[0];for(const n of nice)if(Math.abs(Math.log(n/want))<Math.abs(Math.log(best/want)))best=n;
  $('cstep').value=best;syncOut();
}
function refresh(){computeVolume();autoStep();render();ui();}
function commit(desc){
  const c=simCalc();
  SIM.hist.push({desc,sel:[...SIM.sel],R1:c.R1,wR2:c.wR2});if(SIM.hist.length>60)SIM.hist.shift();
  refresh();
}
const setSel=(labels,desc,restart)=>{SIM.sel=new Set(labels);if(restart)SIM.hist=[];commit(desc);};
function toggle(lab,add){
  if(add===SIM.sel.has(lab))return;
  if(add)SIM.sel.add(lab);else SIM.sel.delete(lab);
  commit((add?'+ ':'− ')+lab);
}
function heaviest(){const z=Math.max(...model.atoms.map(a=>Zof(a.el))),A=model.atoms.filter(a=>Zof(a.el)===z);return{labs:A.map(a=>a.lab),el:A[0].el};}

window.simReset=()=>{
  SIM.on=false;SIM.sel=new Set();SIM.hist=[];SIM.cache=null;
  $('simOn').checked=false;$('simOn').disabled=true;$('simBody').hidden=true;$('simAtoms').textContent='';
};
window.simReady=()=>{
  const box=$('simAtoms');box.textContent='';
  for(const a of model.atoms){const b=document.createElement('button');b.className='chip';b.dataset.lab=a.lab;b.textContent=a.lab;b.dataset.tip=a.el+' – click to add to / remove from the model';box.appendChild(b);}
  $('simOn').disabled=false;
};

$('simOn').addEventListener('change',()=>{
  if(!model)return;
  SIM.on=$('simOn').checked;$('simBody').hidden=!SIM.on;
  if(SIM.on){
    if($('type').value!=='diff'){$('type').value='diff';$('cstep').value=0.05;syncOut();}
    if(!SIM.hist.length){const h=heaviest();setSel(h.labs,'start: '+h.el,true);}else refresh();
  }else refresh();
});
$('simAtoms').addEventListener('click',e=>{const b=e.target.closest('button');if(b&&model)toggle(b.dataset.lab,!SIM.sel.has(b.dataset.lab));});
$('simHeavy').addEventListener('click',()=>{if(!model)return;const h=heaviest();setSel(h.labs,'start: '+h.el,true);});
$('simNonH').addEventListener('click',()=>{if(model)setSel(model.atoms.filter(a=>a.el!=='H').map(a=>a.lab),'all non-H atoms');});
$('simAll').addEventListener('click',()=>{if(model)setSel(model.atoms.map(a=>a.lab),'all atoms');});
$('simNone').addEventListener('click',()=>{if(model)setSel([],'no atoms',true);});
$('simAdd').addEventListener('click',()=>{const q=missingPeaks()[0];if(q)toggle(q.lab,true);});
$('simUndo').addEventListener('click',()=>{
  if(SIM.hist.length<2)return;SIM.hist.pop();SIM.sel=new Set(SIM.hist[SIM.hist.length-1].sel);refresh();
});
$('simAuto').addEventListener('change',()=>{SIM.auto=$('simAuto').checked;if(model&&SIM.on)refresh();});
$('cstep').addEventListener('input',()=>{if(SIM.auto){SIM.auto=false;$('simAuto').checked=false;}}); // moving the slider by hand switches the automatic step off
$('simGhost').addEventListener('change',()=>{SIM.ghost=$('simGhost').checked;if(model)render();});
$('simRefit').addEventListener('change',()=>{SIM.refit=$('simRefit').checked;if(model&&SIM.on)refresh();});

// click a peak in the map: the nearest missing atom (also one lying up to ~1.2 A beside the plane) is added; Shift+click removes an atom of the model
let dn=null;
cv.addEventListener('pointerdown',e=>{dn={x:e.clientX,y:e.clientY,b:e.button};});
cv.addEventListener('pointerup',e=>{
  if(!SIM.on||!dn||!last||dn.b>0)return;
  const moved=Math.hypot(e.clientX-dn.x,e.clientY-dn.y);dn=null;if(moved>4)return;
  const r=cv.getBoundingClientRect();if(!r.width)return;
  const ext=last[2],u=-ext+(e.clientX-r.left)/r.width*2*ext,v=ext-(e.clientY-r.top)/r.height*2*ext,add=!e.shiftKey,sl=+$('slab').value||0.4;
  let best=null,bd=add?0.9:0.8; // score = distance within the plane + a small penalty for the height above/below it
  for(const a of last[3]){if(SIM.sel.has(a.lab)===add)continue;if(!add&&Math.abs(a.dz)>=sl)continue; // atoms of the model are removed only when drawn
    const d=Math.hypot(a.u-u,a.v-v)+(add?0.25*Math.abs(a.dz):0);if(d<bd){bd=d;best=a;}}
  if(best)toggle(best.lab,add);
});
})();
