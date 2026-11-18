"use strict";
// Full screen for the maps. "Full screen" here means the whole browser window (no Fullscreen API, so it behaves the same
// on desktop and tablet): the stage of the active tab is pinned over the window and everything else is hidden.
// The canvases, their state (section, zoom, tilt, profile line) and all handlers stay untouched; only the CSS box changes.
// Mouse: button in the top right corner of the map (Esc leaves). Touch: button in the touch bar, exit button top right.
(function(){
const $=id=>document.getElementById(id);
const IN='<svg viewBox="0 0 24 24" width="18" height="18" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><path d="M4 9V4h5M20 9V4h-5M4 15v5h5M20 15v5h-5"/></svg>';
const OUT='<svg viewBox="0 0 24 24" width="18" height="18" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><path d="M9 4v5H4M15 4v5h5M9 20v-5H4M15 20v-5h5"/></svg>';
const TIP_IN='Full screen: the map fills the browser window. The controls are hidden; Esc or the button at the top right goes back. Section, zoom and profile line stay as they are.';
const TIP_OUT='Back to the normal view (Esc).';
const stages=[$('main-fourier'),$('main-patterson')];
const active=()=>$('main-patterson').hidden?$('main-fourier'):$('main-patterson');
let cur=null;

function mk(st){
  const d=document.createElement('div');d.className='fsctl';
  d.innerHTML='<button type="button" class="fspf" data-tip="Shows or hides the profile strip over the map.">Hide profile</button><button type="button" class="fsbtn" aria-label="Full screen" data-tip="'+TIP_IN+'">'+IN+'</button>';
  st.appendChild(d);
  d.querySelector('.fsbtn').addEventListener('click',()=>cur?leave():enter(st));
  d.querySelector('.fspf').addEventListener('click',e=>{const off=st.classList.toggle('pfoff');e.currentTarget.textContent=off?'Show profile':'Hide profile';after();});
}
stages.forEach(mk);

// touch bar (made by touch.js): an extra button for entering and leaving
let tb=null;
const bar=$('touchbar');
if(bar){
  tb=document.createElement('button');tb.type='button';tb.dataset.fs='1';tb.textContent='Full screen';
  const row=bar.querySelector('div');if(row)row.appendChild(tb);
  tb.addEventListener('click',()=>cur?leave():enter(active()));
}

function sync(){
  const on=!!cur;
  stages.forEach(st=>{
    const b=st.querySelector('.fsbtn'),isCur=st===cur;
    b.innerHTML=isCur?OUT:IN;b.setAttribute('aria-label',isCur?'Exit full screen':'Full screen');b.dataset.tip=isCur?TIP_OUT:TIP_IN;
  });
  if(tb)tb.textContent=on?'Exit full screen':'Full screen';
  if(tb)tb.classList.toggle('on',on);
}
// the maps are CSS-scaled canvases and the profile strips redraw on 'resize', so one resize event after the layout change is all they need
function after(){requestAnimationFrame(()=>window.dispatchEvent(new Event('resize')));}

function enter(st){
  if(cur||!st)return;
  cur=st;
  st.classList.remove('pfoff');const pf=st.querySelector('.fspf');if(pf)pf.textContent='Hide profile'; // the profile starts as in the normal view
  st.classList.add('fs');document.body.classList.add('fsmode');
  sync();after();
}
function leave(){
  if(!cur)return;
  cur.classList.remove('fs','pfoff');document.body.classList.remove('fsmode');
  cur=null;sync();after();
}
document.addEventListener('keydown',e=>{if(e.key==='Escape'&&cur){e.preventDefault();leave();}});
sync();
})();
