// Collapsible cards: remember which ones are open (localStorage, if available)
(function(){
let st={};try{st=JSON.parse(localStorage.getItem('fm-cards')||'{}');}catch(_){}
document.querySelectorAll('details[data-key]').forEach(d=>{
  const k=d.dataset.key;if(k in st)d.open=!!st[k];
  d.addEventListener('toggle',()=>{st[k]=d.open;try{localStorage.setItem('fm-cards',JSON.stringify(st));}catch(_){}});
});
})();

// Redraw the maps when the colour scheme changes: the canvases read the theme colours only when they are drawn
(function(){
const mq=matchMedia('(prefers-color-scheme: dark)');
const redraw=()=>{try{
  if(typeof render==='function')render();
  if(window.profileUpdate)profileUpdate();
  const p=document.getElementById('ppos');if(p)p.dispatchEvent(new Event('input')); // Patterson tab redraws through its own handler
}catch(_){}};
if(mq.addEventListener)mq.addEventListener('change',redraw);else if(mq.addListener)mq.addListener(redraw);
})();
