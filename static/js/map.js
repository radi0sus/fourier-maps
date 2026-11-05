const TP=2*Math.PI;
// f0 via XFilePrep's IT92 table (s = sin(theta)/lambda, here given as s^2)
function f0(el,s2){const f=scatteringFactor(el,Math.sqrt(s2));if(f===null)throw new Error('No scattering factors for '+el);return f;}
function inv3(m){const[a,b,c]=m[0],[d,e,f]=m[1],[g,h,i]=m[2];const A=e*i-f*h,B=-(d*i-f*g),C=d*h-e*g,det=a*A+b*B+c*C;
 return[[A/det,-(b*i-c*h)/det,(b*f-c*e)/det],[B/det,(a*i-c*g)/det,-(a*f-c*d)/det],[C/det,-(a*h-b*g)/det,(a*e-b*d)/det]];}
const mv=(M,v)=>[0,1,2].map(i=>M[i][0]*v[0]+M[i][1]*v[1]+M[i][2]*v[2]);
const dot=(a,b)=>a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
const cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
const unit=v=>{const n=Math.hypot(...v);return v.map(x=>x/n);};

function buildModel(cif){
  const it=cif.items,L=t=>cif.loops.find(l=>l.tags.includes(t));
  const r=Math.PI/180,P=k=>num(it['_cell_'+k]);
  const a=P('length_a'),b=P('length_b'),c=P('length_c'),al=P('angle_alpha'),be=P('angle_beta'),ga=P('angle_gamma');
  const ca=Math.cos(al*r),cb=Math.cos(be*r),cg=Math.cos(ga*r),sg=Math.sin(ga*r);
  const V=a*b*c*Math.sqrt(1-ca*ca-cb*cb-cg*cg+2*ca*cb*cg);
  const M=[[a,b*cg,c*cb],[0,b*sg,c*(ca-cb*cg)/sg],[0,0,V/(a*b*sg)]],Mi=inv3(M);
  const rs=[b*c*Math.sin(al*r)/V,a*c*Math.sin(be*r)/V,a*b*sg/V];
  const sl=L('_space_group_symop_operation_xyz')||L('_symmetry_equiv_pos_as_xyz');
  const st=sl.tags.find(t=>/xyz$/.test(t));
  const ops=sl.rows.map(row=>parseSym(row[sl.tags.indexOf(st)]));
  const disp={},dispI={},dl=L('_atom_type_symbol');
  if(dl&&dl.tags.includes('_atom_type_scat_dispersion_real'))for(const w of dl.rows){disp[w[0]]=num(w[dl.tags.indexOf('_atom_type_scat_dispersion_real')]);const k=dl.tags.indexOf('_atom_type_scat_dispersion_imag');if(k>=0)dispI[w[0]]=num(w[k]);}
  const an={},al2=L('_atom_site_aniso_label');
  if(al2)for(const w of al2.rows){const g=t=>num(w[al2.tags.indexOf('_atom_site_aniso_u_'+t)]);an[w[0]]=[g('11'),g('22'),g('33'),g('23'),g('13'),g('12')];}
  const sx=L('_atom_site_fract_x'),ix=t=>sx.tags.indexOf(t),atoms=[];
  for(const w of sx.rows){
    const lab=w[ix('_atom_site_label')];let ts=ix('_atom_site_type_symbol')>=0?w[ix('_atom_site_type_symbol')]:lab;
    ts=ts.replace(/[^A-Za-z]/g,'');const el=ts[0].toUpperCase()+ts.slice(1).toLowerCase();
    const oc=ix('_atom_site_occupancy')>=0?num(w[ix('_atom_site_occupancy')]):1;
    atoms.push({lab,el,x:['x','y','z'].map(q=>num(w[ix('_atom_site_fract_'+q)])),occ:isNaN(oc)?1:oc,
      fp:disp[ts]||0,fpp:dispI[ts]||0,Ui:num(w[ix('_atom_site_u_iso_or_equiv')]),Ua:an[lab]||null});
  }
  const wm=/\(([\d.]+)P\)\^2\^\s*\+\s*([\d.]+)P/.exec(it['_refine_ls_weighting_details']||'');
  const ref={R1:num(it._refine_ls_r_factor_gt),wR2:num(it._refine_ls_wr_factor_ref),n:num(it._reflns_number_total),Rint:num(it._diffrn_reflns_av_r_equivalents),dmax:num(it._diffrn_reflns_theta_max)};
  for(const at of atoms){let k=0;for(const o of ops){let ok=true;for(let i=0;i<3;i++){const d=o.R[i][0]*at.x[0]+o.R[i][1]*at.x[1]+o.R[i][2]*at.x[2]+o.t[i]-at.x[i];if(Math.abs(d-Math.round(d))>0.002)ok=false;}if(ok)k++;}at.stab=k||1;}
  return{cell:{a,b,c,al,be,ga},V,M,Mi,rs,ops,atoms,wa:wm?+wm[1]:0,wb:wm?+wm[2]:0,ref,f000cif:num(it._exptl_crystal_f_000),refl:readRefl(cif)};
}
// F000 = electrons per cell (incl. f'): value from the CIF if given, else summed from the model
function f000(m){
  if(m.f000cif>0)return{v:m.f000cif,src:'CIF'};
  let s=0;for(const at of m.atoms)s+=m.ops.length*at.occ/at.stab*(f0(at.el,0)+at.fp);
  return{v:s,src:'model'};
}
function readRefl(cif){
  const L=t=>cif.loops.find(l=>l.tags.includes(t)),out=[];
  let l=L('_refln_index_h'),pre='_refln_',fi='f_squared_meas',si='f_squared_sigma';
  if(!l){l=L('_diffrn_refln_index_h');pre='_diffrn_refln_';fi='intensity_net';si='intensity_u';}
  if(l){
    const q=t=>l.tags.indexOf(pre+t),h=q('index_h'),k=q('index_k'),m=q('index_l'),I=q(fi),S=q(si)>=0?q(si):q('intensity_sigma');
    if(I<0||S<0)throw new Error('Reflection loop has no intensity/sigma columns');
    for(const w of l.rows)out.push([+w[h],+w[k],+w[m],num(w[I]),num(w[S])]);
  }else if(cif.items._shelx_hkl_file){
    for(const r of parseHKL(cif.items._shelx_hkl_file).reflections)out.push([r.h,r.k,r.l,r.f2,r.sigma]);
  }
  if(!out.length)throw new Error('No reflection data found in this CIF');
  return out;
}
const dsq=(g,h)=>{let s=0;for(let i=0;i<3;i++){const v=g.Mi[0][i]*h[0]+g.Mi[1][i]*h[1]+g.Mi[2][i]*h[2];s+=v*v;}return s;}; // 1/d^2
const hRt=(R,h)=>[0,1,2].map(j=>h[0]*R[0][j]+h[1]*R[1][j]+h[2]*R[2][j]); // R^T h
const cmp=(a,b)=>a[0]-b[0]||a[1]-b[1]||a[2]-b[2];

function mergeRefl(m){
  const mp=new Map();
  for(const[h,k,l,I,s]of m.refl){
    if(!(s>0)||isNaN(I))continue;let best=null;
    for(const o of m.ops){const g=hRt(o.R,[h,k,l]).map(Math.round);
      for(const sg of[1,-1]){const c=g.map(x=>sg*x+0);if(!best||cmp(c,best)>0)best=c;}}
    const key=best.join(),w=1/(s*s);let e=mp.get(key);
    if(!e){e={h:best,sw:0,swi:0,n:0,v:[]};mp.set(key,e);}e.sw+=w;e.swi+=w*I;e.n++;e.v.push(I);e.ws=(e.ws||[]);e.ws.push(w);
  }
  return[...mp.values()].map(e=>{const I=e.swi/e.sw;let s=1/Math.sqrt(e.sw);
    if(e.n>1){let q=0;for(let i=0;i<e.n;i++)q+=e.ws[i]*(e.v[i]-I)**2;s=Math.max(s,Math.sqrt(q/((e.n-1)*e.sw)));} // larger of internal/external esd, as in SHELXL
    return{h:e.h,I,s,n:e.n,v:e.v};});
}
function calcF(m,h){
  const s2=dsq(m,h)/4,X=2*Math.PI*Math.PI;let A=0,B=0;
  const ops=m.ops.map(o=>({hr:hRt(o.R,h),ht:dot(h,o.t)}));
  for(const at of m.atoms){
    const f=f0(at.el,s2)+at.fp;
    for(const o of ops){
      const g=o.hr;let q;
      if(at.Ua){const u=at.Ua,x=g[0]*m.rs[0],y=g[1]*m.rs[1],z=g[2]*m.rs[2];
        q=x*x*u[0]+y*y*u[1]+z*z*u[2]+2*(x*y*u[5]+x*z*u[4]+y*z*u[3]);}
      else q=at.Ui*4*s2;
      const w=at.occ/at.stab*f*Math.exp(-X*q),p=TP*(dot(g,at.x)+o.ht);
      const c=Math.cos(p),s=Math.sin(p),wi=at.occ/at.stab*at.fpp*Math.exp(-X*q);A+=w*c-wi*s;B+=w*s+wi*c;
    }
  }
  return[A,B];
}
function prepare(m){
  let U=mergeRefl(m);let sa=0,sb=0,nm=0;
  for(const u of U){const r=calcF(m,u.h);u.Fc2=r[0]*r[0]+r[1]*r[1];}
  const fmax=Math.max(...U.map(u=>u.Fc2)),nAll=U.length;
  U=U.filter(u=>u.Fc2>1e-9*fmax); // systematic absences (Fc = 0 by symmetry) are dropped, as SHELXL does
  for(const u of U){const r=calcF(m,u.h);u.Fc2=r[0]*r[0]+r[1]*r[1];u.Fc=Math.sqrt(u.Fc2);u.phi=Math.atan2(r[1],r[0]);
    u.d=1/Math.sqrt(dsq(m,u.h));u.Fo=Math.sqrt(Math.max(u.I,0));nm+=u.v.length;
    if(u.v.length>1){const mean=u.v.reduce((x,y)=>x+y,0)/u.v.length;for(const x of u.v){sa+=Math.abs(x-mean);sb+=Math.abs(x);}}}
  let G=0;{let a=0,b=0;for(const u of U){a+=u.I*u.Fc2;b+=u.I*u.I;}G=a/b;}
  const wt=u=>{const P=(Math.max(u.I*G,0)+2*u.Fc2)/3,s=u.s*G;return 1/(s*s+(m.wa*P)**2+m.wb*P);};
  for(let it=0;it<6;it++){let a=0,b=0;for(const u of U){const w=wt(u);a+=w*u.I*u.Fc2;b+=w*u.I*u.I;}G=a/b;} // SHELXL-like OSF on F^2
  let nw=0,dw=0,rn=0,rd=0,ng=0;
  for(const u of U){const w=wt(u),f2=u.I*G;nw+=w*(f2-u.Fc2)**2;dw+=w*f2*f2;
    if(f2>2*u.s*G){const f=Math.sqrt(f2);rn+=Math.abs(f-u.Fc);rd+=f;ng++;}}
  return{U,K:Math.sqrt(G),R1:rn/rd,wR2:Math.sqrt(nw/dw),Rint:sa/sb,nSel:ng,nRaw:m.refl.length,nAbs:nAll-U.length,red:nm/U.length};
}
function planeFrame(m,labels,off=0){
  const P=labels.map(l=>{const a=m.atoms.find(x=>x.lab===l);if(!a)throw new Error('Atom not found: '+l);return mv(m.M,a.x);});
  const o=[0,1,2].map(i=>(P[0][i]+P[1][i]+P[2][i])/3),e1=unit(P[1].map((v,i)=>v-P[0][i]));
  const n=unit(cross(e1,P[2].map((v,i)=>v-P[0][i]))),e2=cross(n,e1);return{o:o.map((v,i)=>v+off*n[i]),e1,e2,n}; // off = parallel shift of the section along its normal (A)
}
function mapTerms(m,pr,type,dmin){
  const T=new Map();
  for(const u of pr.U){
    if(u.d<dmin)continue;const Fo=pr.K*u.Fo;
    const amp=type==='diff'?Fo-u.Fc:type.startsWith('fo')?Fo:u.Fc;
    for(const o of m.ops){
      let g=hRt(o.R,u.h).map(Math.round),p=u.phi-TP*dot(u.h,o.t);
      const f=g.find(x=>x!==0);if(f===undefined)continue;if(f<0){g=g.map(x=>-x+0);p=-p;}
      const key=g.join();if(!T.has(key))T.set(key,{g,amp,p});
    }
  }
  return[...T.values()];
}
async function computeMap(m,terms,fr,ext,step,onProg,c0=0){
  const n=Math.round(2*ext/step)+1,gs=2*ext/(n-1),N=terms.length;
  const fo=mv(m.Mi,fr.o),f1=mv(m.Mi,fr.e1),f2=mv(m.Mi,fr.e2);
  const A=new Float64Array(N),Ph=new Float64Array(N),a=new Float64Array(N),p=new Float64Array(N),q=new Float64Array(N);
  terms.forEach((t,i)=>{A[i]=t.amp;Ph[i]=t.p;a[i]=TP*dot(t.g,fo);p[i]=TP*dot(t.g,f1);q[i]=TP*dot(t.g,f2);});
  const g=new Float32Array(n*n);
  for(let j=0;j<n;j++){
    const v=-ext+j*gs;
    for(let i=0;i<n;i++){const u=-ext+i*gs;let s=0;
      for(let k=0;k<N;k++)s+=A[k]*Math.cos(Ph[k]-a[k]-u*p[k]-v*q[k]);
      g[j*n+i]=2*s/m.V+c0;} // c0 = F000/V for absolute maps
    if(j%8===0){onProg&&onProg(j/n);await new Promise(r=>setTimeout(r));}
  }
  return{g,n,gs};
}
function contour(g,n,lv){
  const S=[];
  for(let j=0;j<n-1;j++)for(let i=0;i<n-1;i++){
    const a=g[j*n+i],b=g[j*n+i+1],c=g[(j+1)*n+i+1],d=g[(j+1)*n+i],P=[];
    const I=(p,q,x1,y1,x2,y2)=>{const t=(lv-p)/(q-p);P.push([x1+(x2-x1)*t,y1+(y2-y1)*t]);};
    if((a>lv)!==(b>lv))I(a,b,i,j,i+1,j);
    if((b>lv)!==(c>lv))I(b,c,i+1,j,i+1,j+1);
    if((d>lv)!==(c>lv))I(d,c,i,j+1,i+1,j+1);
    if((a>lv)!==(d>lv))I(a,d,i,j,i,j+1);
    for(let k=0;k+1<P.length;k+=2)S.push([P[k][0],P[k][1],P[k+1][0],P[k+1][1]]);
  }
  return S;
}
function overlay(m,fr,ext,slab=0.4){ // atoms (incl. symmetry copies) within 0.4 A of the plane
  const out=[];
  for(const at of m.atoms)for(const o of m.ops)for(let tx=-1;tx<=1;tx++)for(let ty=-1;ty<=1;ty++)for(let tz=-1;tz<=1;tz++){
    const x=[0,1,2].map(i=>dot(o.R[i],at.x)+o.t[i]+[tx,ty,tz][i]),c=mv(m.M,x),d=c.map((v,i)=>v-fr.o[i]);
    const u=dot(d,fr.e1),v=dot(d,fr.e2);
    if(Math.abs(dot(d,fr.n))<slab&&Math.abs(u)<ext&&Math.abs(v)<ext&&!out.some(z=>Math.hypot(z.u-u,z.v-v)<.05))out.push({u,v,el:at.el,c});
  }
  return out;
}
