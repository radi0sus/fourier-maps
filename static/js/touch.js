"use strict";
// Touch / tablet mode. Additive: mouse input goes through the old code untouched.
// Touch and pen pointers get (1) a mode bar (Pan / Turn / Profile) that stands in for the modifier keys,
// (2) two-finger pinch zoom + pan, (3) tap read-out, (4) tap tooltips. Add ?touch=1 to the URL to force the bar, ?touch=0 to hide it.
(function(){
const A=$('cv'),B=$('cvp'),T=e=>e.pointerType==='touch'||e.pointerType==='pen',isPat=()=>!$('main-patterson').hidden;
const Q=location.search,OFF=/[?&]touch=0/.test(Q),MODS={pan:{},turn:{altKey:true},line:{shiftKey:true}};
let mode='pan',synth=false,pinch=null,tipT=0;const fing=new Map();
const show=()=>{if(!OFF)document.body.classList.add('touch');};
if(matchMedia('(pointer:coarse)').matches||/[?&]touch=1/.test(Q))show();
const eff=()=>isPat()&&mode==='turn'?'pan':mode;

// ---- mode bar ----
const bar=document.createElement('div');bar.id='touchbar';
bar.innerHTML='<div><button data-m="pan">Pan</button><button data-m="turn">Turn</button><button data-m="line">Profile</button><button data-r="1">Reset view</button></div><div class="tb-hint"></div>';
document.body.appendChild(bar);
const HINT={pan:['Drag pans · pinch zooms · tap adds a model atom (Model card)','Drag pans · pinch zooms · tap reads out a point'],
  turn:['Drag turns the section · tap picks atoms for a new plane',''],
  line:['Drag from atom to atom draws a profile line · tap removes a model atom','Drag draws a profile line']};
function setMode(m){
  if(m)mode=m;const e=eff(),pat=isPat();
  bar.querySelectorAll('[data-m]').forEach(b=>{b.classList.toggle('on',b.dataset.m===e);if(b.dataset.m==='turn')b.disabled=pat;});
  bar.querySelector('.tb-hint').textContent=HINT[e][pat?1:0];
  try{const O=window.OrientCore&&OrientCore.state;if(O){O.mod=!pat&&e!=='pan';O.hover=true;if(O.mod&&model&&vol)live();}}catch(_){} // like a held modifier key: shows the molecule skeleton
}
bar.addEventListener('click',e=>{const b=e.target.closest('button');if(!b)return;
  if(b.dataset.r)(isPat()?B:A).dispatchEvent(new MouseEvent('dblclick',{bubbles:true}));else setMode(b.dataset.m);});
document.querySelectorAll('.tab').forEach(t=>t.addEventListener('click',()=>setTimeout(()=>setMode(),0)));
setMode();

// ---- pointers ----
const fire=(el,ev)=>{synth=true;try{el.dispatchEvent(ev);}finally{synth=false;}};
const ptr=(el,type,x,y,id=9001)=>fire(el,new PointerEvent(type,{pointerId:id,pointerType:'mouse',button:0,buttons:type==='pointerup'?0:1,clientX:x,clientY:y,bubbles:true,cancelable:true}));
const shadow=(e,m)=>{for(const k of['shiftKey','metaKey','ctrlKey','altKey'])Object.defineProperty(e,k,{value:e[k]||!!m[k],configurable:true});}; // the mode acts like a held key
const mine=e=>!synth&&T(e)&&(e.target===A||e.target===B);
function down(e){
  if(!mine(e))return;show();
  fing.set(e.pointerId,{x:e.clientX,y:e.clientY,x0:e.clientX,y0:e.clientY,t:e.pointerType});
  if(pinch||fing.size>2){e.stopImmediatePropagation();return;}
  if(fing.size===2){ // second finger: cancel the first finger's gesture, then pinch/pan with both
    const[[id,f],[,g]]=[...fing.entries()],c=e.target,cx=(f.x+g.x)/2,cy=(f.y+g.y)/2;
    fire(c,new PointerEvent('pointercancel',{pointerId:id,pointerType:f.t,bubbles:true}));
    pinch={el:c,d:Math.hypot(f.x-g.x,f.y-g.y)||1,over:false};ptr(c,'pointerdown',cx,cy);
    e.stopImmediatePropagation();e.preventDefault();return;
  }
  shadow(e,MODS[eff()]);
}
function move(e){
  if(synth||!T(e)||!fing.has(e.pointerId))return;
  const f=fing.get(e.pointerId);f.x=e.clientX;f.y=e.clientY;
  if(!pinch)return;e.stopImmediatePropagation();
  if(pinch.over||fing.size<2)return;
  const[a,b]=[...fing.values()],cx=(a.x+b.x)/2,cy=(a.y+b.y)/2,d=Math.hypot(a.x-b.x,a.y-b.y)||1,s=d/pinch.d;
  if(Math.abs(Math.log(s))>0.002){fire(pinch.el,new WheelEvent('wheel',{deltaY:-Math.log(s)/0.0015,clientX:cx,clientY:cy,bubbles:true,cancelable:true}));pinch.d=d;} // same zoom code as the mouse wheel
  ptr(pinch.el,'pointermove',cx,cy);
}
function up(e){
  if(synth||!T(e)||!fing.has(e.pointerId))return;
  const f=fing.get(e.pointerId);fing.delete(e.pointerId);
  if(pinch){e.stopImmediatePropagation();
    if(!pinch.over){pinch.over=true;ptr(pinch.el,'pointerup',-9e3,-9e3);} // far away, so that no atom is added
    if(!fing.size)pinch=null;return;}
  shadow(e,MODS[eff()]);
  if(e.type==='pointerup'&&Math.hypot(e.clientX-f.x0,e.clientY-f.y0)<=6){const c=e.target,x=e.clientX,y=e.clientY;setTimeout(()=>ptr(c,'pointermove',x,y),0);} // tap: read-out at that point
}
document.addEventListener('pointerdown',down,true);document.addEventListener('pointermove',move,true);
document.addEventListener('pointerup',up,true);document.addEventListener('pointercancel',up,true);
document.addEventListener('pointerleave',e=>{if(mine(e))e.stopImmediatePropagation();},true); // keep the read-out after the finger is lifted

// ---- tooltips on tap ----
document.addEventListener('click',e=>{
  if(!T(e))return;
  document.querySelectorAll('.tip-open').forEach(x=>x.classList.remove('tip-open'));
  const t=e.target.closest&&e.target.closest('[data-tip]');if(!t)return;
  t.classList.add('tip-open');clearTimeout(tipT);tipT=setTimeout(()=>t.classList.remove('tip-open'),5000);
});
})();
