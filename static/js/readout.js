"use strict";
// Live read-out under the map (Fourier tab): fractional coordinates and the density of the map that is shown, at the pointer.
(function(){
const el=$('freadout'),NAME={diff:'Δρ','2fofc':'2ρ(Fo)−ρ(Fc)',fo:'ρ(Fo)',fc:'ρ(Fc)',fo_abs:'ρ(Fo)',fc_abs:'ρ(Fc)'};
let smp=null,key='';
cv.addEventListener('pointermove',e=>{
  if(!model||!vol||!last||!last[5])return;
  const fr=last[5],r=cv.getBoundingClientRect(),ext=last[2];if(!r.width)return;
  const u=-ext+(e.clientX-r.left)/r.width*2*ext,v=ext-(e.clientY-r.top)/r.height*2*ext,P=[0,1,2].map(i=>fr.o[i]+u*fr.e1[i]+v*fr.e2[i]),f=mv(model.Mi,P);
  const t=$('type').value,abs=t.endsWith('_abs'),k=t+'|'+(vol&&vol.N.join());
  if(!smp||key!==k+'|'+(window.simSig?simSig():'')+'|'+$('dmin').value){smp=makeSampler(vol,abs?f000(model).v/model.V:0,1);key=k+'|'+(window.simSig?simSig():'')+'|'+$('dmin').value;}
  const d=smp(P),val=t==='diff'?d[0]-d[1]:t==='2fofc'?2*d[0]-d[1]:t.startsWith('fo')?d[0]:d[1];
  let near='',bd=0.6;for(const a of last[3]){const s=Math.hypot(a.u-u,a.v-v)+0.15*Math.abs(a.dz);if(s<bd){bd=s;near=`   |   near ${a.lab} (${Math.hypot(a.u-u,a.v-v,a.dz).toFixed(2)} Å)`;}}
  el.textContent=`x y z = ${f.map(x=>x.toFixed(3)).join('  ')}   |   ${NAME[t]||t} = ${val.toFixed(2)} e/Å³${near}`;
});
cv.addEventListener('pointerleave',()=>{el.textContent='';});
})();
