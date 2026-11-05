// Minimal CIF reader (no modules, file:// safe). Tags are lower-cased.
function tokenize(text){
  const out=[],L=text.split(/\r?\n/);
  for(let i=0;i<L.length;i++){
    const s=L[i];
    if(s[0]===';'){const b=[s.slice(1)];while(++i<L.length&&L[i][0]!==';')b.push(L[i]);out.push({v:b.join('\n'),q:1});continue;}
    let p=0;
    while(p<s.length){
      const c=s[p];
      if(c===' '||c==='\t'){p++;continue;}
      if(c==='#')break;
      if(c==="'"||c==='"'){let e=p+1;while(e<s.length&&!(s[e]===c&&(e+1===s.length||/\s/.test(s[e+1]))))e++;out.push({v:s.slice(p+1,e),q:1});p=e+1;continue;}
      let e=p;while(e<s.length&&s[e]!==' '&&s[e]!=='\t')e++;
      out.push({v:s.slice(p,e)});p=e;
    }
  }
  return out;
}
function parseCIF(text){
  const T=tokenize(text),items={},loops=[];let i=0;
  const isTag=t=>!t.q&&t.v[0]==='_',isKey=t=>!t.q&&(isTag(t)||/^(loop_|data_|save_|global_)/i.test(t.v));
  while(i<T.length){
    const t=T[i];
    if(!t.q&&t.v.toLowerCase()==='loop_'){
      i++;const tags=[];while(i<T.length&&isTag(T[i]))tags.push(T[i++].v.toLowerCase());
      const vals=[];while(i<T.length&&!isKey(T[i]))vals.push(T[i++].v);
      const n=tags.length,rows=[];for(let r=0;r+n<=vals.length;r+=n)rows.push(vals.slice(r,r+n));
      loops.push({tags,rows});
    }else if(isTag(t)){items[t.v.toLowerCase()]=i+1<T.length?T[i+1].v:'';i+=2;}
    else i++;
  }
  return{items,loops};
}
const num=s=>parseFloat(String(s).replace(/\(.*?\)/,''));
function parseSym(s){
  const R=[],t=[];
  s.toLowerCase().replace(/\s/g,'').split(',').forEach(c=>{
    const r=[0,0,0];let tr=0,m;const re=/([+-]?)(\d+\/\d+|\d*\.?\d+)?\*?([xyz])?/g;
    while((m=re.exec(c))&&m[0]!==''){
      const sg=m[1]==='-'?-1:1;let v=1;
      if(m[2]){if(m[2].includes('/')){const[a,b]=m[2].split('/');v=a/b;}else v=parseFloat(m[2]);}
      if(m[3])r['xyz'.indexOf(m[3])]+=sg*v;else tr+=sg*v;
    }
    R.push(r);t.push(tr);
  });
  return{R,t};
}
