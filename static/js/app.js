// CPK-like colours (readable on light and dark backgrounds) and covalent radii for bond drawing
const ATOM_COLOR={H:'#9aa0a6',B:'#e8a0a0',C:'#6b7280',N:'#3b6fe0',O:'#e03b3b',F:'#4cc3d9',S:'#d4a017',Cl:'#2eb82e',P:'#e07b2b',Fe:'#e0762b',Cu:'#b87333',Ni:'#3fae3f',Co:'#d96ba0',Zn:'#7d80b0',Mn:'#9c7ac7'};
const RAD={H:.31,B:.84,C:.76,N:.71,O:.66,F:.57,S:1.05,Cl:1.02,P:1.07,Fe:1.32,Cu:1.32,Ni:1.24,Co:1.26,Zn:1.22,Mn:1.39};
let model=null,prep=null,last=null,baseInfo='',fileBase='map';
function hexRGB(s){s=s.trim().replace('#','');if(s.length===3)s=s.replace(/./g,'$&$&');return[0,2,4].map(i=>parseInt(s.slice(i,i+2),16));}
const $=id=>document.getElementById(id),cv=$('cv'),ctx=cv.getContext('2d');
async function loadFile(f){
  if(!f)return;$('dzname').textContent=f.name;fileBase=f.name.replace(/\.[^.]+$/,'');last=null;$('png').disabled=true;
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
    $('go').disabled=false;
  }catch(err){$('info').textContent='Error: '+err.message;$('go').disabled=true;}
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
$('go').addEventListener('click',async()=>{
  try{
    $('go').disabled=true;
    const ext=+$('ext').value,step=+$('gstep').value,cs=+$('cstep').value,type=$('type').value;
    const fr=planeFrame(model,['a1','a2','a3'].map(i=>$(i).value),+$('off').value||0);
    const terms=mapTerms(model,prep,type,+$('dmin').value||0);
    const mp=await computeMap(model,terms,fr,ext,step,p=>{$('go').textContent=Math.round(p*100)+' %';},type.endsWith('_abs')?f000(model).v/model.V:0);
    last=[mp,cs,ext,overlay(model,fr,ext,+$('slab').value||0.4),type];
    const F=f000(model);$('png').disabled=false;
    $('info').textContent=baseInfo+'\n'+draw(...last)+(type.endsWith('_abs')?`\nF000 = ${F.v.toFixed(1)} e (${F.src}) -> +${(F.v/model.V).toFixed(3)} e/Å³ added`:'');
  }catch(err){$('info').textContent='Error: '+err.message;}
  $('go').textContent='Calculate';$('go').disabled=false;
});
function draw(mp,cs,ext,atoms,type){
  const css=getComputedStyle(document.documentElement),fg=css.getPropertyValue('--fg'),bg=css.getPropertyValue('--bg'),mut=css.getPropertyValue('--mut');
  const W=cv.width,sx=i=>i/(mp.n-1)*W,sy=j=>W-j/(mp.n-1)*W,ux=u=>(u+ext)/(2*ext)*W,vy=v=>W-(v+ext)/(2*ext)*W;
  ctx.fillStyle=bg;ctx.fillRect(0,0,W,W);
  let mn=1e9,mx=-1e9;for(const x of mp.g){if(x<mn)mn=x;if(x>mx)mx=x;}
  const col=document.querySelector('input[name=cmode]:checked').value==='col',vm=Math.max(Math.abs(mn),Math.abs(mx))||1;
  if(col){ // shaded density: blue = positive, red = negative, background colour = zero
    const n=mp.n,off=document.createElement('canvas');off.width=off.height=n;const o=off.getContext('2d'),im=o.createImageData(n,n),b=hexRGB(bg);
    for(let j=0;j<n;j++)for(let i=0;i<n;i++){const t=mp.g[j*n+i]/vm,c=t>=0?[47,109,246]:[224,59,59],a=Math.min(1,Math.abs(t))**0.7,p=((n-1-j)*n+i)*4;
      for(let k=0;k<3;k++)im.data[p+k]=b[k]+(c[k]-b[k])*a;im.data[p+3]=255;}
    o.putImageData(im,0,0);ctx.imageSmoothingEnabled=true;const q=W/(2*(n-1));ctx.drawImage(off,0,0,n,n,-q,-q,W*n/(n-1),W*n/(n-1));
  }
  const k0=Math.max(Math.ceil(mn/cs),-40),k1=Math.min(Math.floor(mx/cs),40);
  ctx.strokeStyle=fg;ctx.lineWidth=1.3;ctx.globalAlpha=col?0.6:1;
  for(let k=k0;k<=k1;k++){
    ctx.setLineDash(k>0?[]:k<0?[6,4]:[1.5,3.5]);ctx.beginPath();
    for(const s of contour(mp.g,mp.n,k*cs)){ctx.moveTo(sx(s[0]),sy(s[1]));ctx.lineTo(sx(s[2]),sy(s[3]));}
    ctx.stroke();
  }
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

// options re-draw the last map without recalculating
['showAtoms','showLegend','cm0','cm1'].forEach(id=>$(id).addEventListener('change',()=>{if(last)draw(...last);}));

// PNG export of the canvas as displayed
$('png').addEventListener('click',()=>{if(!last)return;
  cv.toBlob(b=>{const a=document.createElement('a');a.href=URL.createObjectURL(b);a.download=`${fileBase}_${last[4]}.png`;
    document.body.appendChild(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(a.href),1000);},'image/png');});
