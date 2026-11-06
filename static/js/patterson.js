// Patterson tab: P(u,v,w) = 1/V * sum_hkl |F|^2 cos(2 pi (hu+kv+lw)), computed once as a 3D grid (FFT, via density3D from map.js).
// Sections are axis-parallel (u, v or w = const); Harker lines / planes are derived from the symmetry operators and overlaid.
(function(){
'use strict';
const $=id=>document.getElementById(id),cvp=$('cvp'),cx=cvp.getContext('2d');
let pvol=null,pscale=1,hk=[],last=null,active=false,pq=false,pdr=false,pt=null,pqt=null,pdT=null,ppend=false,poff=null;
const pv={x:0,y:0},DEF={plain:{cs:3,sat:25,u:' %'},sharp:{cs:1,sat:6,u:' σ'}};

// ---------- Harker elements from the symmetry operators ----------
const gcd=(a,b)=>b?gcd(b,a%b):a;
function nullspace(A){ // basis of {n : A n = 0}
  const m=A.map(r=>r.slice()),piv=[];let r=0;
  for(let c=0;c<3&&r<3;c++){
    let p=r;for(let i=r+1;i<3;i++)if(Math.abs(m[i][c])>Math.abs(m[p][c]))p=i;
    if(Math.abs(m[p][c])<1e-9)continue;
    [m[r],m[p]]=[m[p],m[r]];const d=m[r][c];m[r]=m[r].map(v=>v/d);
    for(let i=0;i<3;i++)if(i!==r){const f=m[i][c];m[i]=m[i].map((v,j)=>v-f*m[r][j]);}
    piv.push(c);r++;
  }
  const out=[];for(const f of[0,1,2].filter(c=>!piv.includes(c))){const v=[0,0,0];v[f]=1;piv.forEach((pc,i)=>{v[pc]=-m[i][f];});out.push(v);}
  return{rank:piv.length,basis:out};
}
function prim(v){ // smallest integer vector with first non-zero entry > 0
  for(let k=1;k<=12;k++){const w=v.map(x=>x*k);
    if(w.every(x=>Math.abs(x-Math.round(x))<1e-6)){let r=w.map(Math.round),g=0;for(const x of r)g=gcd(g,Math.abs(x));
      if(g>1)r=r.map(x=>x/g);if(r.find(x=>x!==0)<0)r=r.map(x=>-x);return r.map(x=>x+0);}}
  return v.map(x=>Math.round(x*100)/100);
}
const fm1=x=>{let r=x-Math.floor(x);return(r<1e-6||r>1-1e-6)?0:r;};
function fs(c){c=fm1(c);if(c===0)return'0';for(const d of[2,3,4,6,8,12]){const k=c*d;if(Math.abs(k-Math.round(k))<1e-6)return Math.round(k)+'/'+d;}return c.toFixed(3);}
function opStr(o){
  return[0,1,2].map(i=>{let s='';for(let j=0;j<3;j++){const r=o.R[i][j];if(Math.abs(r)<1e-9)continue;
    s+=(r<0?'-':(s?'+':''))+(Math.abs(Math.abs(r)-1)<1e-9?'':Math.abs(r))+'xyz'[j];}
    const t=o.t[i];if(Math.abs(t)>1e-9)s+=(t<0?'-':'+')+fs(Math.abs(t));return s||'0';}).join(', ');
}
const nTerms=n=>n.map((x,i)=>x?(x===1?'+':x===-1?'-':(x>0?'+':'')+x)+'uvw'[i]:'').join('').replace(/^\+/,'');
function harker(m){
  const els=[],seen=new Set();
  for(const o of m.ops){
    const M=[0,1,2].map(i=>[0,1,2].map(j=>(i===j?1:0)-o.R[i][j])),MT=[0,1,2].map(i=>[0,1,2].map(j=>M[j][i]));
    const ns=nullspace(MT);if(ns.rank===0||ns.rank===3)continue; // identity / inversion-like: no Harker element
    const N=ns.basis.map(prim),c1=N.map(n=>fm1(-(n[0]*o.t[0]+n[1]*o.t[1]+n[2]*o.t[2]))),c2=c1.map(c=>fm1(-c));
    const kind=ns.rank===2?'plane':'line',vars=[c1];if(c2.some((c,i)=>Math.abs(c-c1[i])>1e-6))vars.push(c2); // Patterson is centrosymmetric: u and -u
    for(const cs of vars){
      const cons=N.map((n,i)=>({n,c:cs[i]})),key=kind+cons.map(c=>c.n.join()+':'+c.c.toFixed(4)).sort().join('|');
      if(seen.has(key))continue;seen.add(key);
      const al=cons.map(c=>{const nz=c.n.map((x,i)=>x?i:-1).filter(i=>i>=0);return nz.length===1&&c.n[nz[0]]===1?{axis:nz[0],pos:c.c}:null;});
      let label;
      if(kind==='plane')label=al[0]?`${'uvw'[al[0].axis]} = ${fs(cons[0].c)}`:`${nTerms(cons[0].n)} = ${fs(cons[0].c)}`;
      else if(al.every(Boolean)){const co=['u','v','w'];al.forEach((a,i)=>{co[a.axis]=fs(cons[i].c);});label=`[${co.join(', ')}]`;}
      else label=cons.map(c=>`${nTerms(c.n)} = ${fs(c.c)}`).join(' & ');
      els.push({kind,cons,label,align:al.find(Boolean)||null,
        tip:`from  ${opStr(o)}  (${kind==='plane'?'rotation / screw axis':'mirror / glide plane'}): Harker ${kind}`});
    }
  }
  els.sort((a,b)=>(a.kind===b.kind?0:a.kind==='plane'?-1:1)||a.label.localeCompare(b.label));
  return els;
}
// elements as they appear in the current section: whole section, lines, points
function hInSec(sec){
  const out={whole:[],lines:[],pts:[]};
  for(const el of hk){
    const L=[];let ok=true;
    for(const c of el.cons){const np=c.n[sec.pi],nq=c.n[sec.qi],rhs=c.c-c.n[sec.a]*sec.s;
      if(!np&&!nq){if(Math.abs(rhs-Math.round(rhs))>0.01)ok=false;}else L.push({np,nq,rhs,el});}
    if(!ok)continue;
    if(!L.length)out.whole.push(el);
    else if(L.length===1)out.lines.push(L[0]);
    else{const det=L[0].np*L[1].nq-L[0].nq*L[1].np;
      if(Math.abs(det)<1e-9){out.lines.push(L[0],L[1]);continue;}
      const dup=new Set();
      for(let k1=-6;k1<=6;k1++)for(let k2=-6;k2<=6;k2++){const r1=L[0].rhs+k1,r2=L[1].rhs+k2,p=(r1*L[1].nq-r2*L[0].nq)/det,q=(L[0].np*r2-L[1].np*r1)/det;
        if(p<-1.6||p>2.6||q<-1.6||q>2.6)continue;const key=Math.round(p*500)+','+Math.round(q*500);if(dup.has(key))continue;dup.add(key);out.pts.push({p,q,el});}}
  }
  return out;
}

// ---------- coefficients ----------
// Wilson-type normalisation: <|F|^2>(s) = exp(a + slope*s^2) * sum_j f_j^2 ; E^2 - 1 gives a sharpened, origin-removed Patterson
function wilson(){
  const U=prep.U,K2=prep.K*prep.K,els=[...new Set(model.atoms.map(a=>a.el))];
  for(const u of U){const s2=1/(4*u.d*u.d),fe={};for(const e of els)fe[e]=f0(e,s2);let S=0;
    for(const at of model.atoms){const f=fe[at.el]+at.fp;S+=model.ops.length*at.occ/at.stab*(f*f+at.fpp*at.fpp);}u.S=S;u.s2=s2;}
  const idx=U.map((_,i)=>i).sort((a,b)=>U[a].s2-U[b].s2),nb=Math.max(1,Math.min(20,Math.floor(U.length/60))),X=[],Y=[];
  for(let b=0;b<nb;b++){const lo=Math.floor(b*U.length/nb),hi=Math.floor((b+1)*U.length/nb);let si=0,ss=0,s2=0,n=0;
    for(let q=lo;q<hi;q++){const u=U[idx[q]];si+=u.I*K2;ss+=u.S;s2+=u.s2;n++;}
    if(n&&si>0){X.push(s2/n);Y.push(Math.log(si/ss));}}
  let a=0,sl=0;
  if(X.length>=3){const mx=X.reduce((p,c)=>p+c,0)/X.length,my=Y.reduce((p,c)=>p+c,0)/Y.length;let sxy=0,sxx=0;
    for(let i=0;i<X.length;i++){sxy+=(X[i]-mx)*(Y[i]-my);sxx+=(X[i]-mx)**2;}sl=sxy/sxx;a=my-sl*mx;}
  else if(Y.length)a=Y.reduce((p,c)=>p+c,0)/Y.length;
  for(const u of U)u.nrm=Math.exp(a+sl*u.s2)*u.S;
  prep.wil={a,B:-sl/2};
}
function patTerms(dmin,sharp){ // |F|^2 (or E^2-1) as real coefficients, phase 0, F000 left out
  const T=new Map(),K2=prep.K*prep.K;
  for(const u of prep.U){if(u.d<dmin)continue;let co=u.I*K2,cc=u.Fc2;if(sharp){co=co/u.nrm-1;cc=cc/u.nrm-1;}
    for(const o of model.ops){let g=hRt(o.R,u.h).map(Math.round);const f=g.find(x=>x!==0);if(f===undefined)continue;if(f<0)g=g.map(x=>-x+0);
      const k=g.join();if(!T.has(k))T.set(k,{g,fo:co,fc:cc,p:0});}}
  return[...T.values()];
}
function computeP(){
  const sharp=$('pcoef').value==='sharp';if(sharp&&!prep.wil)wilson();
  pvol=density3D(model,patTerms(+$('pdmin').value||0,sharp));
  const g=pvol.fo;
  if(sharp){let s=0,s2=0;for(let i=0;i<g.length;i++){s+=g[i];s2+=g[i]*g[i];}const mu=s/g.length;pscale=1/(Math.sqrt(s2/g.length-mu*mu)||1);}
  else pscale=g[0]>0?100/g[0]:1; // origin peak = 100
}

// ---------- section, rendering ----------
function secFrame(){
  const a=+document.querySelector('input[name=paxis]:checked').value,[pi,qi]=[0,1,2].filter(i=>i!==a),s=+$('ppos').value;
  const col=i=>[model.M[0][i],model.M[1][i],model.M[2][i]],Ap=col(pi),Aq=col(qi),e1=unit(Ap),n=unit(cross(Ap,Aq)),e2=cross(n,e1);
  const f=[0,0,0];f[a]=s;f[pi]=.5;f[qi]=.5;
  return{a,pi,qi,s,fr:{o:mv(model.M,f),e1,e2,n},ext0:(Math.hypot(...Ap)+Math.hypot(...Aq))/2};
}
function prender(){
  if(!pvol||!active)return;
  try{
    const fast=pdr||pq,sec=secFrame(),ext=sec.ext0/(+$('pzoom').value||1),step=Math.max(+$('pgstep').value,2*ext/((fast?121:241)-1)),cs=+$('pcstep').value;
    const fr={o:sec.fr.o.map((v,i)=>v+pv.x*sec.fr.e1[i]+pv.y*sec.fr.e2[i]),e1:sec.fr.e1,e2:sec.fr.e2,n:sec.fr.n};
    const mp=samplePlane(model,pvol,fr,ext,step,$('pmap').value,0,fast);
    for(let i=0;i<mp.g.length;i++)mp.g[i]*=pscale;
    sec.fr=fr;last={mp,ext,sec};$('ppng').disabled=false;
    const txt=pdraw(mp,cs,ext,sec),H=$('pshowH').checked?hInSec(sec):null;
    let s=`section ${'uvw'[sec.a]} = ${sec.s.toFixed(3)}   ·   ${txt}`;
    if(H){const w=H.whole.map(e=>e.label),l=[...new Set([...H.lines.map(x=>x.el.label),...H.pts.map(x=>x.el.label)])];
      if(w.length)s+=`\nthis section is a Harker section: ${w.join(', ')}`;
      if(l.length)s+=`\nHarker elements in this section: ${l.join(', ')}`;
      if(!w.length&&!l.length)s+='\nno Harker element lies in this section';}
    if(!fast){const P=peaks(mp,ext,sec);s+=P.length?'\npeaks (u v w, height):\n'+P.map(p=>`  ${p.f.map(x=>x.toFixed(3)).join('  ')}   ${p.v.toFixed(1)}${DEF[$('pcoef').value==='sharp'?'sharp':'plain'].u}`).join('\n'):'\nno peaks in this section';}
    $('pinfo').textContent=s;
  }catch(err){$('pinfo').textContent='Error: '+err.message;}
}
function peaks(mp,ext,sec){
  const{g,n,gs}=mp,fr=sec.fr,out=[];
  for(let j=1;j<n-1;j++)for(let i=1;i<n-1;i++){const v=g[j*n+i];if(v<=0)continue;
    if(v>g[j*n+i-1]&&v>=g[j*n+i+1]&&v>g[(j-1)*n+i]&&v>=g[(j+1)*n+i]&&v>g[(j-1)*n+i-1]&&v>=g[(j+1)*n+i+1]&&v>g[(j-1)*n+i+1]&&v>=g[(j+1)*n+i-1]){
      const X=-ext+i*gs,Y=-ext+j*gs,f=mv(model.Mi,[0,1,2].map(k=>fr.o[k]+X*fr.e1[k]+Y*fr.e2[k])).map(x=>{const r=x-Math.floor(x);return r>1-1e-6?0:r;});
      const dd=mv(model.M,f.map(x=>x-Math.round(x)));if(Math.hypot(...dd)<0.8)continue; // origin peak
      if(out.some(o=>o.f.every((x,k)=>{const d=Math.abs(x-f[k]);return Math.min(d,1-d)<0.02;})))continue;out.push({f,v});}}
  return out.sort((a,b)=>b.v-a.v).slice(0,8);
}
function pdraw(mp,cs,ext,sec){
  const css=getComputedStyle(document.documentElement),fg=css.getPropertyValue('--fg'),bg=css.getPropertyValue('--bg'),mut=css.getPropertyValue('--mut');
  const W=cvp.width,n=mp.n,sx=i=>i/(n-1)*W,sy=j=>W-j/(n-1)*W,fr=sec.fr,unitS=$('pcoef').value==='sharp'?' σ':' % of origin';
  const px=f=>{const c=mv(model.M,f),d=[c[0]-fr.o[0],c[1]-fr.o[1],c[2]-fr.o[2]];return[(dot(d,fr.e1)+ext)/(2*ext)*W,W-(dot(d,fr.e2)+ext)/(2*ext)*W];};
  const F=(p,q)=>{const f=[0,0,0];f[sec.a]=sec.s;f[sec.pi]=p;f[sec.qi]=q;return px(f);};
  cx.setLineDash([]);cx.globalAlpha=1;cx.fillStyle=bg;cx.fillRect(0,0,W,W);
  let mn=1e9,mx=-1e9;for(const x of mp.g){if(x<mn)mn=x;if(x>mx)mx=x;}
  const col=$('pcm1').checked,vm=+$('psat').value||1;
  if(col){const b=hexRGB(bg);if(!poff||poff.width!==n){poff=document.createElement('canvas');poff.width=poff.height=n;}
    const o=poff.getContext('2d'),im=o.createImageData(n,n);
    for(let j=0;j<n;j++)for(let i=0;i<n;i++){const t=Math.max(-1,Math.min(1,mp.g[j*n+i]/vm)),c=t>=0?[47,109,246]:[224,59,59],a=Math.abs(t)**0.7,p=((n-1-j)*n+i)*4;
      for(let k=0;k<3;k++)im.data[p+k]=b[k]+(c[k]-b[k])*a;im.data[p+3]=255;}
    o.putImageData(im,0,0);cx.imageSmoothingEnabled=true;const q=W/(2*(n-1));cx.drawImage(poff,0,0,n,n,-q,-q,W*n/(n-1),W*n/(n-1));}
  const k0=Math.max(Math.ceil(mn/cs),-80),k1=Math.min(Math.floor(mx/cs),80);
  cx.strokeStyle=fg;cx.lineWidth=1.3;cx.globalAlpha=col?0.6:1;
  const segs=contourAll(mp.g,mp.n,cs,k0,k1);
  [[0,[]],[1,[6,4]],[2,[1.5,3.5]]].forEach(([t,dash])=>{const S=segs[t];if(!S.n)return;cx.setLineDash(dash);cx.beginPath();
    for(let q=0;q<S.n;q+=4){cx.moveTo(sx(S.a[q]),sy(S.a[q+1]));cx.lineTo(sx(S.a[q+2]),sy(S.a[q+3]));}cx.stroke();});
  cx.globalAlpha=1;cx.setLineDash([]);
  // Harker overlay
  const H=$('pshowH').checked?hInSec(sec):null,ORG='#e8890c',GRN='#12a150';
  if(H){
    cx.lineWidth=2;
    for(const L of H.lines){const pl=L.el.kind==='plane',nn=L.np*L.np+L.nq*L.nq,dn=Math.sqrt(nn);cx.strokeStyle=pl?ORG:GRN;cx.setLineDash(pl?[9,5]:[]);
      for(let k=-12;k<=12;k++){const r=L.rhs+k,bp=L.np*r/nn,bq=L.nq*r/nn,dp=-L.nq/dn*10,dq=L.np/dn*10,A=F(bp-dp,bq-dq),B=F(bp+dp,bq+dq);
        cx.beginPath();cx.moveTo(A[0],A[1]);cx.lineTo(B[0],B[1]);cx.stroke();}}
    cx.setLineDash([]);cx.strokeStyle=GRN;cx.lineWidth=2.5;
    for(const P of H.pts){const[x,y]=F(P.p,P.q);cx.beginPath();cx.arc(x,y,9,0,7);cx.stroke();}
    if(H.whole.length){cx.strokeStyle=ORG;cx.lineWidth=6;cx.strokeRect(3,3,W-6,W-6);}
  }
  // unit-cell outline of the section with ticks
  cx.strokeStyle=mut;cx.lineWidth=1;cx.setLineDash([3,3]);cx.beginPath();
  const c=[F(0,0),F(1,0),F(1,1),F(0,1)];cx.moveTo(c[0][0],c[0][1]);for(let i=1;i<4;i++)cx.lineTo(c[i][0],c[i][1]);cx.closePath();cx.stroke();cx.setLineDash([]);
  cx.fillStyle=fg;cx.font='12px system-ui';cx.textAlign='center';
  for(const t of[.25,.5,.75]){const a=F(t,0),b=F(0,t);cx.fillText(t,a[0],a[1]+15);cx.textAlign='right';cx.fillText(t,b[0]-6,b[1]+4);cx.textAlign='center';}
  cx.font='italic 14px system-ui';cx.fillText('uvw'[sec.pi],c[1][0]+4,c[1][1]+17);cx.textAlign='right';cx.fillText('uvw'[sec.qi],c[3][0]-6,c[3][1]-4);
  cx.textAlign='right';cx.font='12px system-ui';cx.fillText('0',c[0][0]-6,c[0][1]+15);
  // scale bar
  const raw=2*ext/6,nice=[0.5,1,2,5,10,20,50].filter(v=>v<=raw).pop()||0.5,L=nice*W/(2*ext);
  cx.strokeStyle=fg;cx.fillStyle=fg;cx.lineWidth=2;cx.beginPath();cx.moveTo(20,W-24);cx.lineTo(20+L,W-24);cx.stroke();cx.textAlign='left';cx.font='13px system-ui';cx.fillText(nice+' Å',20,W-32);
  if($('pshowLegend').checked){
    const bw=300,bh=(col?64:44)+(H?18:0),bx=W-bw-16,by=W-bh-16;cx.globalAlpha=.9;cx.fillStyle=bg;cx.fillRect(bx,by,bw,bh);cx.globalAlpha=1;
    cx.strokeStyle=mut;cx.lineWidth=1;cx.strokeRect(bx,by,bw,bh);cx.font='12px system-ui';cx.textAlign='left';cx.fillStyle=fg;cx.lineWidth=1.3;
    let x=bx+10;for(const[d,t,w]of[[[],'positive',85],[[6,4],'negative',90],[[1.5,3.5],'zero',0]]){cx.strokeStyle=fg;cx.setLineDash(d);cx.beginPath();cx.moveTo(x,by+16);cx.lineTo(x+22,by+16);cx.stroke();cx.setLineDash([]);cx.fillText(t,x+27,by+20);x+=w;}
    let y=by+36;
    if(col){const g=cx.createLinearGradient(bx+10,0,bx+150,0),b=hexRGB(bg);
      for(let k=0;k<=8;k++){const t=k/4-1,cc=t>=0?[47,109,246]:[224,59,59],a=Math.abs(t)**0.7;g.addColorStop(k/8,`rgb(${[0,1,2].map(i=>Math.round(b[i]+(cc[i]-b[i])*a)).join()})`);}
      cx.fillStyle=g;cx.fillRect(bx+10,y-8,140,10);cx.strokeStyle=mut;cx.strokeRect(bx+10,y-8,140,10);cx.fillStyle=fg;
      cx.textAlign='left';cx.fillText('−'+vm,bx+10,y+14);cx.textAlign='center';cx.fillText('0',bx+80,y+14);cx.textAlign='right';cx.fillText('+'+vm,bx+150,y+14);cx.textAlign='left';
      cx.fillText(`step ${cs}${unitS}`,bx+160,y+1);y+=36;}
    else{cx.fillStyle=fg;cx.fillText(`step ${cs}${unitS}`,bx+10,y);y+=18;}
    if(H){cx.setLineDash([9,5]);cx.strokeStyle=ORG;cx.lineWidth=2;cx.beginPath();cx.moveTo(bx+10,y-4);cx.lineTo(bx+32,y-4);cx.stroke();cx.setLineDash([]);cx.fillStyle=fg;cx.fillText('Harker plane',bx+37,y);
      cx.strokeStyle=GRN;cx.beginPath();cx.moveTo(bx+135,y-4);cx.lineTo(bx+157,y-4);cx.stroke();cx.fillText('Harker line / point',bx+162,y);}
  }
  return `range ${mn.toFixed(1)} … ${mx.toFixed(1)}${unitS}, step ${cs} (solid +, dashed −, dotted 0)`;
}

// ---------- UI wiring ----------
const live=()=>{if(!ppend){ppend=true;requestAnimationFrame(()=>{ppend=false;prender();});}};
const quick=()=>{pq=true;clearTimeout(pqt);pqt=setTimeout(()=>{pq=false;live();},150);};
function recompute(){clearTimeout(pdT);pdT=setTimeout(()=>{if(!model||!prep||!active)return;computeP();prender();},150);}
function ensure(){
  if(!model||!prep){$('pinfo').textContent='Load a CIF first (drop zone above).';return;}
  if(!pvol){$('pinfo').textContent='Computing Patterson map ...';setTimeout(()=>{computeP();prender();},20);}else prender();
}
function setMode(){ // sensible defaults per coefficient type
  const d=DEF[$('pcoef').value==='sharp'?'sharp':'plain'];
  $('pcstep').value=d.cs;$('pcstep').dataset.u=d.u;$('psat').value=d.sat;$('psat').dataset.u=d.u;syncOut();
}
document.querySelectorAll('.tab').forEach(b=>b.addEventListener('click',()=>{
  document.querySelectorAll('.tab').forEach(x=>x.classList.toggle('on',x===b));
  $('tab-fourier').hidden=b.dataset.tab!=='fourier';$('tab-patterson').hidden=b.dataset.tab!=='patterson';
  active=b.dataset.tab==='patterson';if(active)ensure();
}));
['pcoef'].forEach(id=>$(id).addEventListener('change',()=>{setMode();pvol=null;pv.x=pv.y=0;recompute();}));
$('pdmin').addEventListener('input',recompute);
['pmap','ppos','pzoom','pgstep','pcstep','psat','pshowH','pshowLegend','pcm0','pcm1','pax0','pax1','pax2'].forEach(id=>['input','change'].forEach(e=>$(id).addEventListener(e,live)));
['ppos','pzoom','pgstep'].forEach(id=>$(id).addEventListener('input',quick));
document.querySelectorAll('input[name=paxis]').forEach(r=>r.addEventListener('change',()=>{pv.x=pv.y=0;}));
document.querySelectorAll('[data-pos]').forEach(b=>b.addEventListener('click',()=>{$('ppos').value=b.dataset.pos;syncOut();live();}));
$('hkList').addEventListener('click',e=>{const b=e.target.closest('button[data-i]');if(!b)return;const el=hk[+b.dataset.i];if(!el.align)return;
  $('pax'+el.align.axis).checked=true;$('ppos').value=el.align.pos;pv.x=pv.y=0;syncOut();live();});
// pan / zoom / reset / cursor read-out
cvp.addEventListener('pointerdown',e=>{if(!pvol||e.button>0)return;pt={x:e.clientX,y:e.clientY};pdr=true;cvp.classList.add('drag');try{cvp.setPointerCapture(e.pointerId);}catch(_){}});
cvp.addEventListener('pointermove',e=>{
  const r=cvp.getBoundingClientRect();if(!r.width||!last)return;
  if(pt){const k=2*last.ext/r.width;pv.x-=(e.clientX-pt.x)*k;pv.y+=(e.clientY-pt.y)*k;pt={x:e.clientX,y:e.clientY};live();return;}
  const{mp,ext,sec}=last,X=-ext+(e.clientX-r.left)/r.width*2*ext,Y=ext-(e.clientY-r.top)/r.height*2*ext,fr=sec.fr;
  const f=mv(model.Mi,[0,1,2].map(k=>fr.o[k]+X*fr.e1[k]+Y*fr.e2[k])).map(x=>x-Math.floor(x)),i=Math.round((X+ext)/mp.gs),j=Math.round((Y+ext)/mp.gs);
  const v=(i>=0&&j>=0&&i<mp.n&&j<mp.n)?mp.g[j*mp.n+i]:NaN;
  $('preadout').textContent=`u v w = ${f.map(x=>x.toFixed(3)).join('  ')}   |   P = ${isNaN(v)?'–':v.toFixed(1)}${$('pcoef').value==='sharp'?' σ':' %'}`;});
const endDrag=()=>{if(!pt)return;pt=null;pdr=false;cvp.classList.remove('drag');live();};
cvp.addEventListener('pointerup',endDrag);cvp.addEventListener('pointercancel',endDrag);
cvp.addEventListener('wheel',e=>{if(!pvol||!last)return;e.preventDefault();const r=cvp.getBoundingClientRect();if(!r.width)return;
  const z=$('pzoom'),z0=+z.value,z1=Math.min(+z.max,Math.max(+z.min,z0*Math.exp(-e.deltaY*0.0015))),e0=last.sec.ext0/z0,e1=last.sec.ext0/z1;
  const fx=(e.clientX-r.left)/r.width*2-1,fy=1-(e.clientY-r.top)/r.height*2;pv.x+=fx*(e0-e1);pv.y+=fy*(e0-e1);z.value=z1;syncOut();quick();live();},{passive:false});
cvp.addEventListener('dblclick',()=>{pv.x=pv.y=0;$('pzoom').value=1;syncOut();live();});
$('ppng').addEventListener('click',()=>{if(!last)return;cvp.toBlob(b=>{const a=document.createElement('a');a.href=URL.createObjectURL(b);
  a.download=`${fileBase}_patterson_${'uvw'[last.sec.a]}${last.sec.s.toFixed(3)}.png`;document.body.appendChild(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(a.href),1000);},'image/png');});

// hooks called from app.js
window.patReset=()=>{pvol=null;last=null;hk=[];pv.x=pv.y=0;$('ppng').disabled=true;$('hkList').textContent='';$('preadout').textContent='';};
window.patReady=()=>{
  hk=harker(model);const L=$('hkList');L.textContent='';
  hk.forEach((el,i)=>{const b=document.createElement('button');b.className='chip'+(el.align?'':' off');b.dataset.i=i;b.dataset.tip=el.tip+(el.align?'':' (not axis-parallel: drawn as overlay only)');
    b.textContent=el.label+(el.kind==='plane'?' (plane)':' (line)');L.appendChild(b);});
  if(!hk.length)L.textContent='none (no rotation, screw or mirror/glide operators)';
  if(active)ensure();
};
setMode();
window.PattersonCore={harker,nullspace,patTerms,wilson,computeP:()=>{computeP();return{pvol,pscale};}};
})();
