"use strict";
// Turning the section. Cmd / Ctrl / Alt + drag rotates the section plane around the midpoint of the three plane atoms;
// Cmd / Ctrl / Alt + click picks atoms (3 picks define a new plane and bring it into the picture plane).
// The molecule (the bonded cluster around the plane) is drawn as a pale skeleton with atom names while the key is held,
// or permanently with "Show molecule". The map itself is still cut out of the 3D grid, so any plane works.
(function(){
const I3=()=>[[1,0,0],[0,1,0],[0,0,1]];
const O={R:I3(),custom:null,picks:[],mod:false,hover:false,rot:false,base:null,mol:null,molKey:''};
const mulM=(A,B)=>A.map(r=>[0,1,2].map(j=>r[0]*B[0][j]+r[1]*B[1][j]+r[2]*B[2][j]));
const mulV=(A,v)=>A.map(r=>r[0]*v[0]+r[1]*v[1]+r[2]*v[2]);
const rodrigues=(a,t)=>{const c=Math.cos(t),s=Math.sin(t),x=a[0],y=a[1],z=a[2],k=1-c;
  return[[c+x*x*k,x*y*k-z*s,x*z*k+y*s],[y*x*k+z*s,c+y*y*k,y*z*k-x*s],[z*x*k-y*s,z*y*k+x*s,c+z*z*k]];};

// ---- the frame of the section: plane through three points (atoms of the list, or picked atoms), then turned by R ----
function basePoints(m,labels){
  if(O.custom)return O.custom.pts;
  return labels.map(l=>{const a=m.atoms.find(x=>x.lab===l);if(!a)throw new Error('Atom not found: '+l);return mv(m.M,a.x);});
}
function baseFrame(P){
  const o=[0,1,2].map(i=>(P[0][i]+P[1][i]+P[2][i])/3),e1=unit(P[1].map((v,i)=>v-P[0][i]));
  const n=unit(cross(e1,P[2].map((v,i)=>v-P[0][i])));return{o,e1,e2:cross(n,e1),n,P0:P[0]};
}
window.orientFrame=(m,labels,off)=>{ // same as planeFrame when nothing is turned
  const b=baseFrame(basePoints(m,labels)),e1=mulV(O.R,b.e1),e2=mulV(O.R,b.e2),n=mulV(O.R,b.n);O.base=b;
  return{o:b.o.map((v,i)=>v+off*n[i]),e1,e2,n}; // the centre of the plane atoms is the pivot; the offset moves the section along the turned normal
};

// ---- the bonded molecule around the pivot (all symmetry images and cell translations, bonds by covalent radii as in the map) ----
function molecule(b){
  const key=[model.atoms.length,model.V.toFixed(2),...b.o.map(x=>x.toFixed(2)),...b.P0.map(x=>x.toFixed(2))].join('|');
  if(O.mol&&O.molKey===key)return O.mol;
  const R=14,lo=[1e9,1e9,1e9],hi=[-1e9,-1e9,-1e9],cand=[],seen=new Set();
  for(const sx of[-1,1])for(const sy of[-1,1])for(const sz of[-1,1]){const f=mv(model.Mi,[b.o[0]+sx*R,b.o[1]+sy*R,b.o[2]+sz*R]);for(let i=0;i<3;i++){lo[i]=Math.min(lo[i],f[i]);hi[i]=Math.max(hi[i],f[i]);}}
  for(const at of model.atoms)for(const o of model.ops){
    const x0=[0,1,2].map(i=>dot(o.R[i],at.x)+o.t[i]);
    for(let tx=Math.ceil(lo[0]-x0[0]);tx<=Math.floor(hi[0]-x0[0]);tx++)for(let ty=Math.ceil(lo[1]-x0[1]);ty<=Math.floor(hi[1]-x0[1]);ty++)for(let tz=Math.ceil(lo[2]-x0[2]);tz<=Math.floor(hi[2]-x0[2]);tz++){
      const c=mv(model.M,[x0[0]+tx,x0[1]+ty,x0[2]+tz]);if(Math.hypot(c[0]-b.o[0],c[1]-b.o[1],c[2]-b.o[2])>R)continue;
      const k=at.lab+'|'+c.map(v=>Math.round(v*50)).join(',');if(seen.has(k))continue;seen.add(k);cand.push({lab:at.lab,el:at.el,c});}}
  let seed=0,sd=1e9;cand.forEach((a,i)=>{const d=Math.hypot(a.c[0]-b.P0[0],a.c[1]-b.P0[1],a.c[2]-b.P0[2]);if(d<sd){sd=d;seed=i;}});
  const rad=a=>RAD[a.el]||1.4,bonded=(p,q)=>Math.hypot(p.c[0]-q.c[0],p.c[1]-q.c[1],p.c[2]-q.c[2])<rad(p)+rad(q)+0.4;
  const inC=new Set([seed]),queue=[seed],bonds=[];
  while(queue.length){const i=queue.pop();for(let j=0;j<cand.length;j++){if(j===i||!bonded(cand[i],cand[j]))continue;
    if(i<j||!inC.has(j))bonds.push([i,j]);if(!inC.has(j)){inC.add(j);queue.push(j);}}}
  const idx=[...inC],map=new Map(idx.map((v,k)=>[v,k])),seenB=new Set(),bb=[];
  for(const[i,j]of bonds){const a=map.get(i),c=map.get(j);if(a===undefined||c===undefined)continue;const k=a<c?a+','+c:c+','+a;if(!seenB.has(k)){seenB.add(k);bb.push([a,c]);}}
  O.molKey=key;return O.mol={atoms:idx.map(i=>cand[i]),bonds:bb};
}

// ---- drawing: called from draw() in app.js, below the atoms of the section ----
const showMol=()=>$('orMol').checked||(O.mod&&O.hover)||O.rot;
window.orientDraw=(ux,vy,fg,bg,mut)=>{
  if(!model||!last||!last[5]||!O.base)return;
  const fr=last[5],sl=+$('slab').value||0.4,proj=c=>{const d=[c[0]-fr.o[0],c[1]-fr.o[1],c[2]-fr.o[2]];return{u:dot(d,fr.e1),v:dot(d,fr.e2),w:dot(d,fr.n)};};
  ctx.save();ctx.textAlign='left';ctx.font='11px system-ui';
  if(showMol()){
    const M=molecule(O.base),P=M.atoms.map(a=>Object.assign({a},proj(a.c)));
    ctx.lineWidth=1.6;ctx.strokeStyle=mut;
    for(const[i,j]of M.bonds){const p=P[i],q=P[j],w=(Math.abs(p.w)+Math.abs(q.w))/2;ctx.globalAlpha=Math.max(0.15,Math.min(0.7,0.8-0.12*w));
      ctx.beginPath();ctx.moveTo(ux(p.u),vy(p.v));ctx.lineTo(ux(q.u),vy(q.v));ctx.stroke();}
    for(const p of P){const inS=Math.abs(p.w)<sl,al=inS?1:Math.max(0.3,Math.min(0.85,0.95-0.12*Math.abs(p.w)));
      if(!inS){ctx.globalAlpha=al;ctx.fillStyle=ATOM_COLOR[p.a.el]||'#a050c0';ctx.beginPath();ctx.arc(ux(p.u),vy(p.v),3.2,0,7);ctx.fill();}
      if(p.a.el==='H'&&!inS)continue; // no names for hydrogen atoms outside the slab
      ctx.globalAlpha=al;ctx.fillStyle=fg;ctx.fillText(p.a.lab,ux(p.u)+(inS?9:5),vy(p.v)+(inS?4:-4));}
  }
  ctx.globalAlpha=1;
  O.picks.forEach((pk,i)=>{const p=proj(pk.c);ctx.strokeStyle='#e8a200';ctx.lineWidth=2.5;ctx.beginPath();ctx.arc(ux(p.u),vy(p.v),11,0,7);ctx.stroke();ctx.fillStyle='#e8a200';ctx.fillText(String(i+1),ux(p.u)+13,vy(p.v)-9);});
  if(O.rot){const p=proj(O.base.o),x=ux(p.u),y=vy(p.v);ctx.strokeStyle=fg;ctx.lineWidth=1.5;ctx.beginPath();ctx.moveTo(x-9,y);ctx.lineTo(x+9,y);ctx.moveTo(x,y-9);ctx.lineTo(x,y+9);ctx.stroke();} // pivot
  ctx.restore();
};
window.orientInfo=()=>{
  let s='';if(O.custom)s+=`\nSection plane through the picked atoms ${O.custom.labs.join(', ')}`;
  if(O.base&&last&&last[5]){const ang=Math.acos(Math.max(-1,Math.min(1,dot(last[5].n,O.base.n))))*180/Math.PI;
    if(ang>0.05)s+=`\nSection turned by ${ang.toFixed(1)}° against the plane through its three atoms (Reset tilt in the Plane card)`;}
  return s;
};
function ui(){
  const e=$('orInfo');
  e.textContent=O.picks.length?`Picked ${O.picks.length}/3: ${O.picks.map(p=>p.lab).join(', ')}. Cmd/Ctrl/Alt + click the next atom.`:
    'Cmd / Ctrl / Alt + drag turns the section, + click picks three atoms for a new plane.';
}
window.orientReset=()=>{O.R=I3();O.custom=null;O.picks=[];O.mol=null;O.molKey='';O.base=null;if($('orMol'))ui();};

// ---- gestures ----
const turn=(dx,dy,k)=>{ // dx, dy: pointer movement in px, k: radians per px. The molecule follows the pointer: the frame turns the other way
  const fr=last&&last[5];if(!fr)return;
  const w=[0,1,2].map(i=>fr.e2[i]*dx*k+fr.e1[i]*dy*k),ang=Math.hypot(...w);if(ang<1e-9)return;
  O.R=mulM(rodrigues(w.map(x=>x/ang),-ang),O.R);
};
function pick(e){
  const fr=last&&last[5];if(!fr||!O.base)return;
  const r=cv.getBoundingClientRect(),ext=last[2],u=-ext+(e.clientX-r.left)/r.width*2*ext,v=ext-(e.clientY-r.top)/r.height*2*ext;
  const cand=last[3].map(a=>({lab:a.lab,c:a.c,u:a.u,v:a.v,w:a.dz}));
  for(const a of molecule(O.base).atoms){const d=[0,1,2].map(i=>a.c[i]-fr.o[i]);cand.push({lab:a.lab,c:a.c,u:dot(d,fr.e1),v:dot(d,fr.e2),w:dot(d,fr.n)});}
  let best=null,bd=0.8;for(const a of cand){const s=Math.hypot(a.u-u,a.v-v)+0.15*Math.abs(a.w);if(s<bd){bd=s;best=a;}}
  if(!best)return;
  if(O.picks.length>=3)O.picks=[];
  if(O.picks.some(p=>Math.hypot(p.c[0]-best.c[0],p.c[1]-best.c[1],p.c[2]-best.c[2])<0.05))return;
  O.picks.push({lab:best.lab,c:best.c.slice()});
  if(O.picks.length===3){ // new plane through the three picked atoms (the images that were clicked), lying in the picture plane
    O.custom={pts:O.picks.map(p=>p.c),labs:O.picks.map(p=>p.lab)};O.R=I3();view.x=view.y=0;
    O.picks.forEach((p,i)=>{$('a'+(i+1)).value=p.lab;});O.picks=[];}
  ui();live();
}
let gs=null;
cv.addEventListener('pointerdown',e=>{
  if(!model||!vol||e.button>0||!(e.metaKey||e.ctrlKey||e.altKey))return;
  gs={x:e.clientX,y:e.clientY,moved:0};O.rot=false;try{cv.setPointerCapture(e.pointerId);}catch(_){}e.preventDefault();
});
cv.addEventListener('pointermove',e=>{
  if(!gs)return;const r=cv.getBoundingClientRect();if(!r.width)return;
  const dx=e.clientX-gs.x,dy=e.clientY-gs.y;gs.moved+=Math.hypot(dx,dy);gs.x=e.clientX;gs.y=e.clientY;
  if(gs.moved>4){O.rot=true;drag=true;turn(dx,dy,Math.PI/r.width);live();} // dragging over the whole width turns by 180 degrees
});
const endGesture=(e,click)=>{
  if(!gs)return;const moved=gs.moved;gs=null;
  if(O.rot){O.rot=false;drag=false;live();}else if(click&&moved<=4)pick(e);
};
cv.addEventListener('pointerup',e=>endGesture(e,true));
cv.addEventListener('pointercancel',e=>endGesture(e,false));
cv.addEventListener('contextmenu',e=>{if(e.ctrlKey)e.preventDefault();}); // Ctrl+click is a right click on macOS
cv.addEventListener('pointerenter',()=>{O.hover=true;if(O.mod&&model&&vol&&!$('orMol').checked)live();});
cv.addEventListener('pointerleave',()=>{O.hover=false;if(O.mod&&model&&vol&&!$('orMol').checked)live();});
const setMod=v=>{if(O.mod===v)return;O.mod=v;cv.classList.toggle('rotmode',v);if(O.hover&&model&&vol&&!$('orMol').checked)live();};
window.addEventListener('keydown',e=>{if(e.key==='Meta'||e.key==='Control'||e.key==='Alt')setMod(true);});
window.addEventListener('keyup',e=>{if(e.key==='Meta'||e.key==='Control'||e.key==='Alt')setMod(e.metaKey||e.ctrlKey||e.altKey);});
window.addEventListener('blur',()=>setMod(false));

// ---- panel ----
['a1','a2','a3'].forEach(id=>$(id).addEventListener('change',()=>{O.R=I3();O.custom=null;O.picks=[];ui();})); // atoms chosen from the lists: plane through those atoms again
$('orMol').addEventListener('change',()=>{if(model&&vol)live();});
$('orReset').addEventListener('click',()=>{O.R=I3();view.x=view.y=0;if(model&&vol)live();});
ui();
window.OrientCore={molecule:()=>molecule(O.base),state:O};
})();
