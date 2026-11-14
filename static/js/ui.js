// Collapsible cards: remember which ones are open (localStorage, if available)
(function(){
let st={};try{st=JSON.parse(localStorage.getItem('fm-cards')||'{}');}catch(_){}
document.querySelectorAll('details[data-key]').forEach(d=>{
  const k=d.dataset.key;if(k in st)d.open=!!st[k];
  d.addEventListener('toggle',()=>{st[k]=d.open;try{localStorage.setItem('fm-cards',JSON.stringify(st));}catch(_){}});
});
})();
