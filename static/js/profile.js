"use strict";
// Density profile along a line (Fo, Fc, Fo-Fc). Shift + drag on the map from atom to atom draws the line; its end points snap
// to atoms (also to atoms that lie beside the section). The profile is cut from the 3D grid, so a bond that leaves the plane works too.
// The strip below the map shows the curves; hovering it (or the line in the map) moves a marker in the other one.
(function(){
const PF={line:null,draft:null,hover:null,dkey:'',dat:null,ver:0};
const pc=$('pfc'),stage=$('main-fourier'),LINE='#0ea5a4',COL={fo:'#2f6df6',fc:'#8a909c',df:'#e03b3b'};
const css=n=>getComputedStyle(document.documentElement).getPropertyValue(n).trim();
const fmt=(x,d=1)=>(Math.abs(x)<10?x.toFixed(d+1):x.toFixed(d));

// ---- density at a Cartesian point: [Fo, Fc] of the maps on screen (relative to the cell mean, plus F000/V for absolute map types) ----
function sampler(){
  const [N1,N2,N3]=vol.N,Mi=model.Mi,NN=N1*N2,[w0,w1,w2]=wb,ix=new Int32Array(4),iy=new Int32Array(4),iz=new Int32Array(4),md=(a,b)=>{a%=b;return a<0?a+b:a;};
  const c0=$('type').value.endsWith('_abs')?f000(model).v/model.V:0;
  return c=>{
    const p0=(Mi[0][0]*c[0]+Mi[0][1]*c[1]+Mi[0][2]*c[2])*N1,p1=(Mi[1][0]*c[0]+Mi[1][1]*c[1]+Mi[1][2]*c[2])*N2,p2=(Mi[2][0]*c[0]+Mi[2][1]*c[1]+Mi[2][2]*c[2])*N3;
    const b0=Math.floor(p0),b1=Math.floor(p1),b2=Math.floor(p2);cubInto(p0-b0,w0);cubInto(p1-b1,w1);cubInto(p2-b2,w2);
    for(let a=0;a<4;a++){ix[a]=md(b0+a-1,N1);iy[a]=md(b1+a-1,N2)*N1;iz[a]=md(b2+a-1,N3)*NN;}
    return[cub(vol.fo,ix,iy,iz,w0,w1,w2)+c0,cub(vol.fc,ix,iy,iz,w0,w1,w2)+c0];
  };
}
// all symmetry images of the atoms within r of the segment A-B, with position s along the line and distance to it
function atomsNear(A,B,r){
  const d=B.map((v,i)=>v-A[i]),len=Math.hypot(...d),u=d.map(x=>x/len),lo=[0,0,0],hi=[0,0,0],out=[];
  // bounding box of the segment in fractional coordinates, padded by r
  const fa=mv(model.Mi,A),fb=mv(model.Mi,B),pad=[0,1,2].map(i=>r*Math.hypot(model.Mi[i][0],model.Mi[i][1],model.Mi[i][2])*1.2);
  for(let i=0;i<3;i++){lo[i]=Math.min(fa[i],fb[i])-pad[i];hi[i]=Math.max(fa[i],fb[i])+pad[i];}
  const seen=new Set();
  for(const at of model.atoms)for(const o of model.ops){
    const x0=[0,1,2].map(i=>dot(o.R[i],at.x)+o.t[i]);
    for(let tx=Math.ceil(lo[0]-x0[0]);tx<=Math.floor(hi[0]-x0[0]);tx++)for(let ty=Math.ceil(lo[1]-x0[1]);ty<=Math.floor(hi[1]-x0[1]);ty++)for(let tz=Math.ceil(lo[2]-x0[2]);tz<=Math.floor(hi[2]-x0[2]);tz++){
      const c=mv(model.M,[x0[0]+tx,x0[1]+ty,x0[2]+tz]),q=c.map((v,i)=>v-A[i]),s=dot(q,u),perp=Math.hypot(q[0]-s*u[0],q[1]-s*u[1],q[2]-s*u[2]);
      if(perp>r)continue;const k=at.lab+'|'+c.map(v=>Math.round(v*50)).join(',');if(seen.has(k))continue;seen.add(k);
      out.push({lab:at.lab,el:at.el,c,s,perp});}}
  return out.sort((a,b)=>a.s-b.s);
}

// ---- the profile of the current line (cached; the key changes with the map, the line and the extension) ----
function data(){
  const L=PF.line;if(!L||!vol||!model)return null;
  const ext=+$('pfext').value,key=[PF.ver,$('type').value,ext,...L.A,...L.B].join('|');
  if(PF.dat&&PF.dkey===key)return PF.dat;
  const d=L.B.map((v,i)=>v-L.A[i]),len=Math.hypot(...d),u=d.map(x=>x/len),P0=L.A.map((v,i)=>v-ext*u[i]),tot=len+2*ext,n=Math.max(2,Math.ceil(tot/0.02)+1),smp=sampler();
  const s=new Float64Array(n),fo=new Float64Array(n),fc=new Float64Array(n),df=new Float64Array(n);
  for(let i=0;i<n;i++){s[i]=tot*i/(n-1);const c=P0.map((v,k)=>v+u[k]*s[i]),r=smp(c);fo[i]=r[0];fc[i]=r[1];df[i]=r[0]-r[1];}
  // atoms near the line (positions measured from the start of the extended line); values at the atom centres
  const near=atomsNear(P0,P0.map((v,i)=>v+u[i]*tot),0.9).map(a=>{const r=smp(a.c);return Object.assign(a,{fo:r[0],fc:r[1],df:r[0]-r[1]});});
  PF.dkey=key;return PF.dat={s,fo,fc,df,n,len,tot,ext,u,P0,near,abs:$('type').value.endsWith('_abs')};
}
const at=(D,s)=>D.P0.map((v,i)=>v+D.u[i]*s); // point on the extended line
const valAt=(D,s)=>{const t=Math.max(0,Math.min(D.n-1,s/D.tot*(D.n-1))),i=Math.min(D.n-2,Math.floor(t)),f=t-i,L=a=>a[i]*(1-f)+a[i+1]*f;return{fo:L(D.fo),fc:L(D.fc),df:L(D.df)};};

// ---- the strip ----
function drawStrip(){
  const D=data();pc.hidden=!D;stage.classList.toggle('haspf',!!D);if(!D)return;
  const r=pc.getBoundingClientRect(),dpr=Math.min(2,window.devicePixelRatio||1),W=Math.round((r.width||760)*dpr),H=Math.round((r.height||170)*dpr);
  if(pc.width!==W||pc.height!==H){pc.width=W;pc.height=H;}
  const g=pc.getContext('2d'),fg=css('--fg'),mut=css('--mut'),bd=css('--bd'),bgc=css('--card');
  g.setTransform(dpr,0,0,dpr,0,0);const w=W/dpr,h=H/dpr;g.clearRect(0,0,w,h);g.fillStyle=bgc;g.fillRect(0,0,w,h);
  const show=[['fo','Fo',$('pfFo').checked],['fc','Fc',$('pfFc').checked],['df','Fo−Fc',$('pfDif').checked]].filter(x=>x[2]);
  const mL=50,mR=14,mT=30,mB=30,pw=w-mL-mR,ph=h-mT-mB;
  let lo=0,hi=0;for(const[k]of show)for(let i=0;i<D.n;i++){const v=D[k][i];if(v<lo)lo=v;if(v>hi)hi=v;}
  if(hi-lo<1e-6){hi=1;lo=-1;}const pad=(hi-lo)*0.08;lo-=pad;hi+=pad;
  const X=s=>mL+s/D.tot*pw,Y=v=>mT+(hi-v)/(hi-lo)*ph;
  g.font='11px system-ui';g.textBaseline='middle';
  // grid and axes
  const step=(()=>{const raw=(hi-lo)/4,p=Math.pow(10,Math.floor(Math.log10(raw))),m=raw/p;return(m<1.5?1:m<3.5?2:m<7.5?5:10)*p;})();
  g.textAlign='right';g.lineWidth=1;
  for(let v=Math.ceil(lo/step)*step;v<=hi+1e-9;v+=step){g.strokeStyle=Math.abs(v)<step/1e3?mut:bd;g.setLineDash(Math.abs(v)<step/1e3?[4,3]:[]);g.beginPath();g.moveTo(mL,Y(v));g.lineTo(mL+pw,Y(v));g.stroke();g.fillStyle=mut;g.fillText(fmt(v,Math.abs(v)<1?2:1).replace(/\.0+$/,''),mL-6,Y(v));}
  g.setLineDash([]);g.strokeStyle=mut;g.beginPath();g.moveTo(mL,mT);g.lineTo(mL,mT+ph);g.lineTo(mL+pw,mT+ph);g.stroke();
  g.textAlign='center';g.fillStyle=mut;const xs=D.tot>6?1:D.tot>3?0.5:0.25;for(let s=0;s<=D.tot+1e-9;s+=xs){g.fillText(s.toFixed(xs<1?(xs<0.5?2:1):0),X(s),mT+ph+11);}
  g.textAlign='right';g.fillText('Å',w-4,mT+ph+22);g.save();g.translate(11,mT+ph/2);g.rotate(-Math.PI/2);g.textAlign='center';g.fillText('e/Å³',0,0);g.restore();
  // end points of the real line (the rest is the extension), atom marks
  g.strokeStyle=LINE;g.lineWidth=1.5;g.setLineDash([2,3]);for(const s of[D.ext,D.ext+D.len]){g.beginPath();g.moveTo(X(s),mT);g.lineTo(X(s),mT+ph);g.stroke();}g.setLineDash([]);
  let rowEnd=[-1e9,-1e9];g.textAlign='center';g.font='11px system-ui';
  for(const a of D.near){if(a.perp>0.6)continue;const x=X(a.s),row=x-rowEnd[0]>34?0:x-rowEnd[1]>34?1:0;rowEnd[row]=x;
    g.strokeStyle=bd;g.lineWidth=1;g.setLineDash([1,3]);g.beginPath();g.moveTo(x,mT+11+row*11);g.lineTo(x,mT+ph);g.stroke();g.setLineDash([]);g.fillStyle=fg;g.fillText(a.lab,x,mT+5+row*11);}
  // curves
  g.lineJoin='round';
  for(const[k]of show){g.strokeStyle=COL[k];g.lineWidth=k==='fc'?1.6:2;g.setLineDash(k==='fc'?[5,3]:[]);g.beginPath();
    for(let i=0;i<D.n;i+=Math.max(1,Math.floor(D.n/(pw*1.5)))){const x=X(D.s[i]),y=Y(D[k][i]);i?g.lineTo(x,y):g.moveTo(x,y);}g.stroke();}
  g.setLineDash([]);
  // title, legend, close mark
  g.textAlign='left';g.fillStyle=fg;g.font='600 12px system-ui';const L=PF.line;
  g.fillText(`${L.labA||'•'} → ${L.labB||'•'}   ${D.len.toFixed(2)} Å`,mL,13);g.font='11px system-ui';
  let lx=w-mR-26;g.textAlign='right';for(const[k,name]of show.slice().reverse()){g.fillStyle=fg;g.fillText(name,lx,13);const tw=g.measureText(name).width;g.strokeStyle=COL[k];g.lineWidth=2;g.setLineDash(k==='fc'?[4,2]:[]);g.beginPath();g.moveTo(lx-tw-22,13);g.lineTo(lx-tw-5,13);g.stroke();g.setLineDash([]);lx-=tw+34;}
  g.fillStyle=mut;g.font='16px system-ui';g.textAlign='right';g.fillText('×',w-6,12);
  // hover cursor with values
  if(PF.hover!=null&&PF.hover>=0&&PF.hover<=D.tot){
    const s=PF.hover,v=valAt(D,s),x=X(s);g.strokeStyle=LINE;g.lineWidth=1.2;g.beginPath();g.moveTo(x,mT);g.lineTo(x,mT+ph);g.stroke();
    const lines=[`${s.toFixed(2)} Å`,...show.map(([k,n])=>`${n} ${v[k].toFixed(2)}`)];g.font='11px system-ui';
    const bw=Math.max(...lines.map(t=>g.measureText(t).width))+14,bh=lines.length*14+8,bx=x+bw+10>mL+pw?x-bw-8:x+8,by=mT+22;
    g.fillStyle=bgc;g.globalAlpha=0.92;g.fillRect(bx,by,bw,bh);g.globalAlpha=1;g.strokeStyle=bd;g.strokeRect(bx,by,bw,bh);
    g.textAlign='left';lines.forEach((t,i)=>{g.fillStyle=i?COL[show[i-1][0]]:fg;g.fillText(t,bx+7,by+11+i*14);});
  }
  pc._geom={mL,pw,w,h,tot:D.tot};
}

// ---- hooks used by app.js ----
window.profileInvalidate=()=>{PF.ver++;PF.dat=null;};
window.profileReset=()=>{PF.line=PF.draft=null;PF.hover=null;PF.dat=null;pc.hidden=true;stage.classList.remove('haspf');ui();};
window.profileUpdate=()=>{if(PF.line)drawStrip();};
window.profileInfo=()=>{
  const D=PF.line&&data();if(!D)return'';
  let s=`\nProfile ${PF.line.labA||'•'} → ${PF.line.labB||'•'}, ${D.len.toFixed(2)} Å (extended by ${D.ext} Å at both ends)${D.abs?', absolute density':''}`;
  const on=D.near.filter(a=>a.perp<0.4&&a.s>=D.ext-0.3&&a.s<=D.ext+D.len+0.3).slice(0,8);
  if(on.length)s+='\n'+on.map(a=>`  ${a.lab}: Fo ${a.fo.toFixed(1)}, Fc ${a.fc.toFixed(1)}, Fo−Fc ${a.df.toFixed(1)} e/Å³`).join('\n');
  return s;
};
const proj=(fr,c)=>{const d=[c[0]-fr.o[0],c[1]-fr.o[1],c[2]-fr.o[2]];return{u:dot(d,fr.e1),v:dot(d,fr.e2),w:dot(d,fr.n)};};
window.profileDraw=(ux,vy,fg,bg,mut)=>{ // the line (or the line being dragged) and the hover marker, in the map
  const L=PF.draft||PF.line,fr=last&&last[5];if(!L||!fr)return;
  const a=proj(fr,L.A),b=proj(fr,L.B),al=w=>Math.max(0.35,Math.min(1,1.1-Math.abs(w)/3));
  ctx.save();ctx.lineCap='round';
  const D=PF.line&&!PF.draft?data():null;
  if(D&&D.ext>0){const s0=proj(fr,at(D,0)),s1=proj(fr,at(D,D.tot));ctx.setLineDash([3,4]);ctx.lineWidth=1.8;ctx.strokeStyle=LINE;ctx.globalAlpha=0.6;ctx.beginPath();ctx.moveTo(ux(s0.u),vy(s0.v));ctx.lineTo(ux(s1.u),vy(s1.v));ctx.stroke();ctx.setLineDash([]);}
  ctx.globalAlpha=Math.min(al(a.w),al(b.w));ctx.strokeStyle=LINE;ctx.lineWidth=3;ctx.beginPath();ctx.moveTo(ux(a.u),vy(a.v));ctx.lineTo(ux(b.u),vy(b.v));ctx.stroke();
  for(const[p,lab]of[[a,L.labA],[b,L.labB]]){ctx.globalAlpha=al(p.w);ctx.lineWidth=2.5;ctx.beginPath();ctx.arc(ux(p.u),vy(p.v),8,0,7);ctx.stroke();}
  if(D&&PF.hover!=null&&PF.hover>=0&&PF.hover<=D.tot){const m=proj(fr,at(D,PF.hover));ctx.globalAlpha=1;ctx.fillStyle=LINE;ctx.strokeStyle=bg;ctx.lineWidth=2;ctx.beginPath();ctx.arc(ux(m.u),vy(m.v),5.5,0,7);ctx.fill();ctx.stroke();}
  ctx.restore();
};

// ---- gestures: Shift + drag draws the line ----
function pointInfo(e){ // in-plane point under the pointer and the atom it snaps to (atoms of the section, and the molecule around it)
  const fr=last&&last[5];if(!fr)return null;
  const r=cv.getBoundingClientRect(),ext=last[2],u=-ext+(e.clientX-r.left)/r.width*2*ext,v=ext-(e.clientY-r.top)/r.height*2*ext;
  const cand=last[3].map(a=>({lab:a.lab,c:a.c,u:a.u,v:a.v,w:a.dz}));
  if(window.OrientCore){try{for(const a of OrientCore.molecule().atoms){const p=proj(fr,a.c);cand.push({lab:a.lab,c:a.c,u:p.u,v:p.v,w:p.w});}}catch(_){}}
  let best=null,bd=0.7;for(const a of cand){const s=Math.hypot(a.u-u,a.v-v)+0.15*Math.abs(a.w);if(s<bd){bd=s;best=a;}}
  return best?{c:best.c.slice(),lab:best.lab}:{c:[0,1,2].map(i=>fr.o[i]+u*fr.e1[i]+v*fr.e2[i]),lab:null};
}
let gs=null;
cv.addEventListener('pointerdown',e=>{
  if(!e.shiftKey||e.metaKey||e.ctrlKey||e.altKey||e.button>0||!model||!vol)return;
  const p=pointInfo(e);if(!p)return;gs={A:p,moved:0,x:e.clientX,y:e.clientY};try{cv.setPointerCapture(e.pointerId);}catch(_){}e.preventDefault();
});
cv.addEventListener('pointermove',e=>{
  if(gs){gs.moved=Math.hypot(e.clientX-gs.x,e.clientY-gs.y);
    if(gs.moved>4){const q=pointInfo(e);if(q){PF.draft={A:gs.A.c,B:q.c,labA:gs.A.lab,labB:q.lab};drag=true;live();}}return;}
  if(!PF.line||!last||!last[5]||PF.draft)return; // hover on the line in the map -> cursor in the strip
  const D=data();if(!D)return;const fr=last[5],r=cv.getBoundingClientRect(),ext=last[2];if(!r.width)return;
  const u=-ext+(e.clientX-r.left)/r.width*2*ext,v=ext-(e.clientY-r.top)/r.height*2*ext,a=proj(fr,at(D,0)),b=proj(fr,at(D,D.tot)),dx=b.u-a.u,dy=b.v-a.v,l2=dx*dx+dy*dy;
  let h=null;if(l2>1e-9){const t=Math.max(0,Math.min(1,((u-a.u)*dx+(v-a.v)*dy)/l2)),d=Math.hypot(u-(a.u+t*dx),v-(a.v+t*dy));if(d<0.35)h=t*D.tot;}
  if(h!==PF.hover){PF.hover=h;drawStrip();}
});
const endDrag=e=>{
  if(!gs)return;const g=gs;gs=null;const d=PF.draft;PF.draft=null;drag=false;
  if(d&&g.moved>4){const len=Math.hypot(d.B[0]-d.A[0],d.B[1]-d.A[1],d.B[2]-d.A[2]);if(len>0.3){PF.line=d;PF.hover=null;PF.dat=null;ui();}}
  live();
};
cv.addEventListener('pointerup',endDrag);cv.addEventListener('pointercancel',endDrag);
cv.addEventListener('pointerleave',()=>{if(PF.hover!=null&&!gs){PF.hover=null;if(PF.line){drawStrip();live();}}});

// ---- strip: hover shows the marker in the map, click on the x removes the line ----
let raf=0;
pc.addEventListener('pointermove',e=>{
  if(!pc._geom||!PF.line)return;const r=pc.getBoundingClientRect(),G=pc._geom,x=(e.clientX-r.left)/r.width*G.w;
  const s=(x-G.mL)/G.pw*G.tot;PF.hover=s>=0&&s<=G.tot?s:null;drawStrip();
  if(!raf)raf=requestAnimationFrame(()=>{raf=0;live();});
});
pc.addEventListener('pointerleave',()=>{PF.hover=null;if(PF.line){drawStrip();live();}});
pc.addEventListener('click',e=>{const r=pc.getBoundingClientRect(),G=pc._geom;if(!G)return;if((e.clientX-r.left)/r.width*G.w>G.w-24&&(e.clientY-r.top)/r.height*G.h<24)clear();});
function clear(){PF.line=PF.draft=null;PF.hover=null;PF.dat=null;pc.hidden=true;stage.classList.remove('haspf');ui();if(model&&vol)live();}

// ---- panel ----
function ui(){const D=PF.line&&model&&vol?data():null;$('pfInfo').textContent=D?`${PF.line.labA||'•'} → ${PF.line.labB||'•'}, ${D.len.toFixed(2)} Å`:'Shift + drag on the map from atom to atom.';}
['pfFo','pfFc','pfDif'].forEach(i=>$(i).addEventListener('change',()=>{if(PF.line)drawStrip();}));
$('pfext').addEventListener('input',()=>{if(PF.line){PF.dat=null;drawStrip();live();}});
$('pfClear').addEventListener('click',clear);
window.addEventListener('resize',()=>{if(PF.line)drawStrip();});
ui();
window.ProfileCore={state:PF,data};
})();
