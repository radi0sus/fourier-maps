// CPK-like colours (readable on light and dark backgrounds) and covalent radii for bond drawing
const ATOM_COLOR={H:'#9aa0a6',B:'#e8a0a0',C:'#6b7280',N:'#3b6fe0',O:'#e03b3b',F:'#4cc3d9',S:'#d4a017',Cl:'#2eb82e',P:'#e07b2b',Fe:'#e0762b',Cu:'#b87333',Ni:'#3fae3f',Co:'#d96ba0',Zn:'#7d80b0',Mn:'#9c7ac7'};
const RAD={H:.31,B:.84,C:.76,N:.71,O:.66,F:.57,S:1.05,Cl:1.02,P:1.07,Fe:1.32,Cu:1.32,Ni:1.24,Co:1.26,Zn:1.22,Mn:1.39};
let model=null,prep=null,last=null,baseInfo='',fileBase='map';
function hexRGB(s){s=s.trim().replace('#','');if(s.length===3)s=s.replace(/./g,'$&$&');return[0,2,4].map(i=>parseInt(s.slice(i,i+2),16));}
const $=id=>document.getElementById(id),cv=$('cv'),ctx=cv.getContext('2d');
async function loadFile(f){
  if(!f)return;$('dzname').textContent=f.name;fileBase=f.name.replace(/\.[^.]+$/,'');last=null;vol=null;view.x=view.y=0;$('png').disabled=true;
  try{
    $('info').textContent='Reading '+f.name+' ...';
    model=buildModel(parseCIF(await f.text()));
    $('info').textContent='Merging and computing structure factors ...';
    await new Promise(r=>setTimeout(r));
    prep=prepare(model);
    $('info').textContent=`${f.name}: ${prep.nRaw} measured -> ${prep.U.length} unique (${prep.nAbs} systematic absences dropped, redundancy ${prep.red.toFixed(1)}), ${prep.nSel} with I>2σ, dmin = ${Math.min(...prep.U.map(u=>u.d)).toFixed(2)} Å
this tool: Rint = ${(100*prep.Rint).toFixed(2)} %, R1 = ${(100*prep.R1).toFixed(2)} %, wR2 = ${(100*prep.wR2).toFixed(2)} %  |  CIF: Rint = ${(100*model.ref.Rint).toFixed(2)} %, R1 = ${(100*model.ref.R1).toFixed(2)} %, wR2 = ${(100*model.ref.wR2).toFixed(2)} %, ${model.ref.n} unique`;
    for(const id of['a1','a2','a3'])$(id).innerHTML=model.atoms.map(a=>`<option>${a.lab}</option>`).join('');
    baseInfo=$('info').textContent;
    model.atoms.filter(a=>a.el!=='H').slice(0,3).forEach((a,i)=>{$('a'+(i+1)).value=a.lab;}); // default plane: first three non-H atoms
    $('info').textContent+='\nComputing 3D density ...';await new Promise(r=>setTimeout(r));
    computeVolume();render();
  }catch(err){$('info').textContent='Error: '+err.message;}
}
$('file').addEventListener('change',e=>loadFile(e.target.files[0]));
const dz=$('drop');
['dragenter','dragover'].forEach(t=>dz.addEventListener(t,e=>{e.preventDefault();dz.classList.add('over');}));
['dragleave','drop'].forEach(t=>dz.addEventListener(t,e=>{e.preventDefault();dz.classList.remove('over');}));
dz.addEventListener('drop',e=>loadFile(e.dataTransfer.files[0]));
// slider read-outs
function syncOut(){document.querySelectorAll('input[type=range]').forEach(r=>{const o=$(r.id+'Out');if(o)o.textContent=(+r.value).toFixed((r.step.split('.')[1]||'').length)+(r.dataset.u||'');});}
document.querySelectorAll('input[type=range]').forEach(r=>r.addEventListener('input',syncOut));syncOut();
$('type').addEventListener('change',()=>{$('cstep').value=$('type').value==='diff'?0.05:0.5;syncOut();});
// ---- live view: 3D density once per file / resolution limit, sections are cut out instantly ----
let offc=null,vol=null,pending=false,dTimer=null,drag=false,quick=false,qt=null;const view={x:0,y:0};
function computeVolume(){vol=density3D(model,volTerms(model,prep,+$('dmin').value||0));}
function render(){
  if(!vol)return;
  try{
    const fast=drag||quick,ext=4/(+$('zoom').value||1),step=Math.max(+$('gstep').value,2*(4/(+$('zoom').value||1))/((fast?121:241)-1)),cs=+$('cstep').value,type=$('type').value,F=f000(model);
    const fr=planeFrame(model,['a1','a2','a3'].map(i=>$(i).value),+$('off').value||0);
    fr.o=fr.o.map((v,i)=>v+view.x*fr.e1[i]+view.y*fr.e2[i]); // pan: shifts the centre of the section within its plane
    const mp=samplePlane(model,vol,fr,ext,step,type,type.endsWith('_abs')?F.v/model.V:0,fast);
    last=[mp,cs,ext,overlay(model,fr,ext,+$('slab').value||0.4),type];$('png').disabled=false;
    $('info').textContent=baseInfo+'\n'+draw(...last)+(type.endsWith('_abs')?`\nF000 = ${F.v.toFixed(1)} e (${F.src}) -> +${(F.v/model.V).toFixed(3)} e/Å³ added`:'');
  }catch(err){$('info').textContent='Error: '+err.message;}
}
const live=()=>{if(!pending){pending=true;requestAnimationFrame(()=>{pending=false;render();});}};
['a1','a2','a3','off','slab','zoom','gstep','cstep','type','showAtoms','showLegend','cm0','cm1'].forEach(id=>['input','change'].forEach(e=>$(id).addEventListener(e,live)));
$('dmin').addEventListener('input',()=>{clearTimeout(dTimer);dTimer=setTimeout(()=>{if(!model)return;computeVolume();render();},150);});
function draw(mp,cs,ext,atoms,type){
  const css=getComputedStyle(document.documentElement),fg=css.getPropertyValue('--fg'),bg=css.getPropertyValue('--bg'),mut=css.getPropertyValue('--mut');
  const W=cv.width,sx=i=>i/(mp.n-1)*W,sy=j=>W-j/(mp.n-1)*W,ux=u=>(u+ext)/(2*ext)*W,vy=v=>W-(v+ext)/(2*ext)*W;
  ctx.fillStyle=bg;ctx.fillRect(0,0,W,W);
  let mn=1e9,mx=-1e9;for(const x of mp.g){if(x<mn)mn=x;if(x>mx)mx=x;}
  const col=document.querySelector('input[name=cmode]:checked').value==='col',vm=Math.max(Math.abs(mn),Math.abs(mx))||1;
  if(col){ // shaded density: blue = positive, red = negative, background colour = zero
    const n=mp.n;if(!offc||offc.width!==n){offc=document.createElement('canvas');offc.width=offc.height=n;}const off=offc,o=off.getContext('2d'),im=o.createImageData(n,n),b=hexRGB(bg);
    for(let j=0;j<n;j++)for(let i=0;i<n;i++){const t=mp.g[j*n+i]/vm,c=t>=0?[47,109,246]:[224,59,59],a=Math.min(1,Math.abs(t))**0.7,p=((n-1-j)*n+i)*4;
      for(let k=0;k<3;k++)im.data[p+k]=b[k]+(c[k]-b[k])*a;im.data[p+3]=255;}
    o.putImageData(im,0,0);ctx.imageSmoothingEnabled=true;const q=W/(2*(n-1));ctx.drawImage(off,0,0,n,n,-q,-q,W*n/(n-1),W*n/(n-1));
  }
  const k0=Math.max(Math.ceil(mn/cs),-40),k1=Math.min(Math.floor(mx/cs),40);
  ctx.strokeStyle=fg;ctx.lineWidth=1.3;ctx.globalAlpha=col?0.6:1;
  const segs=contourAll(mp.g,mp.n,cs,k0,k1);
  [[0,[]],[1,[6,4]],[2,[1.5,3.5]]].forEach(([t,dash])=>{const S=segs[t];if(!S.n)return;ctx.setLineDash(dash);ctx.beginPath();
    for(let q=0;q<S.n;q+=4){ctx.moveTo(sx(S.a[q]),sy(S.a[q+1]));ctx.lineTo(sx(S.a[q+2]),sy(S.a[q+3]));}ctx.stroke();});
  ctx.globalAlpha=1;ctx.setLineDash([]);ctx.font='13px system-ui';ctx.textAlign='center';ctx.strokeStyle=mut;ctx.lineWidth=2.5;
  const at2=$('showAtoms').checked?atoms:[];
  for(let i=0;i<at2.length;i++)for(let j=i+1;j<at2.length;j++){const p=at2[i],q=at2[j];
    if(Math.hypot(p.c[0]-q.c[0],p.c[1]-q.c[1],p.c[2]-q.c[2])<(RAD[p.el]||1.4)+(RAD[q.el]||1.4)+0.4){ctx.beginPath();ctx.moveTo(ux(p.u),vy(p.v));ctx.lineTo(ux(q.u),vy(q.v));ctx.stroke();}}
  for(const a of at2){ctx.fillStyle=ATOM_COLOR[a.el]||'#a050c0';ctx.beginPath();ctx.arc(ux(a.u),vy(a.v),6,0,7);ctx.fill();
    ctx.strokeStyle=bg;ctx.lineWidth=1.5;ctx.stroke();ctx.fillStyle=fg;ctx.fillText(a.el,ux(a.u),vy(a.v)-10);}
  ctx.strokeStyle=fg;
  const L=W/(2*ext);ctx.strokeStyle=fg;ctx.fillStyle=fg;ctx.lineWidth=2;ctx.beginPath();ctx.moveTo(20,W-24);ctx.lineTo(20+L,W-24);ctx.stroke();
  ctx.textAlign='left';ctx.fillText('1 Å',20,W-32);
  if($('showLegend').checked){ // legend box, bottom right
    const bx=W-292,by=W-80,bw=276,bh=64;ctx.globalAlpha=.9;ctx.fillStyle=bg;ctx.fillRect(bx,by,bw,bh);ctx.globalAlpha=1;
    ctx.strokeStyle=mut;ctx.lineWidth=1;ctx.strokeRect(bx,by,bw,bh);ctx.font='12px system-ui';ctx.textAlign='left';ctx.fillStyle=fg;ctx.lineWidth=1.3;
    let x=bx+10;for(const[d,t,w]of[[[],'positive',85],[[6,4],'negative',90],[[1.5,3.5],'zero',0]]){
      ctx.strokeStyle=fg;ctx.setLineDash(d);ctx.beginPath();ctx.moveTo(x,by+16);ctx.lineTo(x+22,by+16);ctx.stroke();ctx.setLineDash([]);ctx.fillText(t,x+27,by+20);x+=w;}
    if(col){const g=ctx.createLinearGradient(bx+10,0,bx+150,0),b=hexRGB(bg);
      for(let k=0;k<=8;k++){const t=k/4-1,c=t>=0?[47,109,246]:[224,59,59],a=Math.abs(t)**0.7;g.addColorStop(k/8,`rgb(${[0,1,2].map(i=>Math.round(b[i]+(c[i]-b[i])*a)).join()})`);}
      ctx.fillStyle=g;ctx.fillRect(bx+10,by+28,140,10);ctx.strokeStyle=mut;ctx.strokeRect(bx+10,by+28,140,10);
      ctx.fillStyle=fg;ctx.textAlign='left';ctx.fillText('−'+vm.toFixed(2),bx+10,by+54);ctx.textAlign='center';ctx.fillText('0',bx+80,by+54);ctx.textAlign='right';ctx.fillText('+'+vm.toFixed(2),bx+150,by+54);ctx.textAlign='left';}
    ctx.fillStyle=fg;ctx.fillText(`step ${cs} e/Å³`,col?bx+160:bx+10,col?by+37:by+42);
  }
  return `map range ${mn.toFixed(2)} … ${mx.toFixed(2)} e/Å³, contour step ${cs} (solid +, dashed −, dotted 0)`;
}

// PNG export of the canvas as displayed
$('png').addEventListener('click',()=>{if(!last)return;
  cv.toBlob(b=>{const a=document.createElement('a');a.href=URL.createObjectURL(b);a.download=`${fileBase}_${last[4]}.png`;
    document.body.appendChild(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(a.href),1000);},'image/png');});

// ---- pan (drag), zoom (wheel, around the pointer), reset (double-click) ----
['a1','a2','a3'].forEach(id=>$(id).addEventListener('change',()=>{view.x=view.y=0;}));
let pt=null;
cv.addEventListener('pointerdown',e=>{if(!vol||e.button>0)return;pt={x:e.clientX,y:e.clientY};drag=true;cv.classList.add('drag');try{cv.setPointerCapture(e.pointerId);}catch(_){}});
cv.addEventListener('pointermove',e=>{if(!pt)return;const r=cv.getBoundingClientRect();if(!r.width)return;const k=2*(4/(+$('zoom').value||1))/r.width;
  view.x-=(e.clientX-pt.x)*k;view.y+=(e.clientY-pt.y)*k;pt={x:e.clientX,y:e.clientY};live();});
const endDrag=()=>{if(!pt)return;pt=null;drag=false;cv.classList.remove('drag');live();};
cv.addEventListener('pointerup',endDrag);cv.addEventListener('pointercancel',endDrag);
cv.addEventListener('wheel',e=>{if(!vol)return;e.preventDefault();const r=cv.getBoundingClientRect();if(!r.width)return;
  const z=$('zoom'),z0=+z.value,z1=Math.min(+z.max,Math.max(+z.min,z0*Math.exp(-e.deltaY*0.0015))),e0=4/z0,e1=4/z1;
  const fx=(e.clientX-r.left)/r.width*2-1,fy=1-(e.clientY-r.top)/r.height*2; // pointer position in the view, -1..1
  view.x+=fx*(e0-e1);view.y+=fy*(e0-e1);z.value=z1;syncOut();live();},{passive:false});
cv.addEventListener('dblclick',()=>{view.x=view.y=0;$('zoom').value=1;syncOut();live();});

// while a slider or the wheel is moving: cheap preview (coarser grid, trilinear), full quality 150 ms after the last change
const quickMode=()=>{quick=true;clearTimeout(qt);qt=setTimeout(()=>{quick=false;live();},150);};
['off','slab','zoom','gstep'].forEach(id=>$(id).addEventListener('input',quickMode));
cv.addEventListener('wheel',quickMode,{passive:true});
