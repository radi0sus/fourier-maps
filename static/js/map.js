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
// ---- live view: the density is computed ONCE as a 3D grid (FFT); sections are then cut out by interpolation ----
function volTerms(m,pr,dmin){ // half-sphere reflections with phase; Fo and Fc amplitudes kept separately
  const T=new Map();
  for(const u of pr.U){if(u.d<dmin)continue;const fo=pr.K*u.Fo;
    for(const o of m.ops){let g=hRt(o.R,u.h).map(Math.round),p=u.phi-TP*dot(u.h,o.t);
      const f=g.find(x=>x!==0);if(f===undefined)continue;if(f<0){g=g.map(x=>-x+0);p=-p;}
      const k=g.join();if(!T.has(k))T.set(k,{g,fo,fc:u.Fc,p});}}
  return[...T.values()];
}
function fftLine(re,im,n,cs,sn,rev){ // in-place radix-2 forward DFT, exp(-i...)
  for(let i=0;i<n;i++){const j=rev[i];if(j>i){let t=re[i];re[i]=re[j];re[j]=t;t=im[i];im[i]=im[j];im[j]=t;}}
  for(let len=2;len<=n;len<<=1){const h=len>>1,st=n/len;
    for(let i=0;i<n;i+=len)for(let k=0;k<h;k++){const wr=cs[k*st],wi=-sn[k*st],a=i+k,b=a+h,xr=re[b]*wr-im[b]*wi,xi=re[b]*wi+im[b]*wr;
      re[b]=re[a]-xr;im[b]=im[a]-xi;re[a]+=xr;im[a]+=xi;}}
}
// Fo + i*Fc are packed into one complex array; both densities are real, so one FFT gives rho_o (real) and rho_c (imaginary)
function density3D(m,terms){
  const hm=[0,0,0];for(const t of terms)for(let i=0;i<3;i++)hm[i]=Math.max(hm[i],Math.abs(t.g[i]));
  const N=hm.map(h=>{let n=32;while(n<4*h+2&&n<256)n*=2;return n;}),[N1,N2,N3]=N,T=N1*N2*N3,re=new Float64Array(T),im=new Float64Array(T);
  const ix=g=>((((g[2]%N3)+N3)%N3)*N2+(((g[1]%N2)+N2)%N2))*N1+(((g[0]%N1)+N1)%N1);
  for(const t of terms){const c=Math.cos(t.p),s=Math.sin(t.p);let k=ix(t.g);
    re[k]+=t.fo*c-t.fc*s;im[k]+=t.fo*s+t.fc*c;
    k=ix(t.g.map(x=>-x));re[k]+=t.fo*c+t.fc*s;im[k]+=-t.fo*s+t.fc*c;} // Friedel partner: conjugates
  for(const[ax,n,stride]of[[0,N1,1],[1,N2,N1],[2,N3,N1*N2]]){
    const cs=new Float64Array(n/2),sn=new Float64Array(n/2),rev=new Uint32Array(n),bits=Math.log2(n);
    for(let k=0;k<n/2;k++){cs[k]=Math.cos(TP*k/n);sn[k]=Math.sin(TP*k/n);}
    for(let i=0;i<n;i++){let r=0;for(let j=0;j<bits;j++)if((i>>j)&1)r|=1<<(bits-1-j);rev[i]=r;}
    const a=new Float64Array(n),b=new Float64Array(n);
    for(let s0=0;s0<T;s0++){const c=ax===0?s0%N1:ax===1?Math.floor(s0/N1)%N2:Math.floor(s0/(N1*N2));if(c)continue;
      for(let i=0;i<n;i++){a[i]=re[s0+i*stride];b[i]=im[s0+i*stride];}
      fftLine(a,b,n,cs,sn,rev);
      for(let i=0;i<n;i++){re[s0+i*stride]=a[i];im[s0+i*stride]=b[i];}}}
  const fo=new Float32Array(T),fc=new Float32Array(T);for(let i=0;i<T;i++){fo[i]=re[i]/m.V;fc[i]=im[i]/m.V;}
  return{N,fo,fc};
}
const wb=[new Float64Array(4),new Float64Array(4),new Float64Array(4)];
function cubInto(t,w){const t2=t*t,t3=t2*t;w[0]=(-t3+2*t2-t)/2;w[1]=(3*t3-5*t2+2)/2;w[2]=(-3*t3+4*t2+t)/2;w[3]=(t3-t2)/2;} // Catmull-Rom weights
function cub(A,ix,iy,iz,w0,w1,w2){let s=0;for(let c=0;c<4;c++)for(let d=0;d<4;d++){const r=iz[c]+iy[d],wy=w2[c]*w1[d];for(let a=0;a<4;a++)s+=wy*w0[a]*A[r+ix[a]];}return s;}
function tri(A,x0,x1,y0,y1,z0,z1,fx,fy,fz){
  const a=A[z0+y0+x0]*(1-fx)+A[z0+y0+x1]*fx,b=A[z0+y1+x0]*(1-fx)+A[z0+y1+x1]*fx,c=A[z1+y0+x0]*(1-fx)+A[z1+y0+x1]*fx,d=A[z1+y1+x0]*(1-fx)+A[z1+y1+x1]*fx;
  return(a*(1-fy)+b*fy)*(1-fz)+(c*(1-fy)+d*fy)*fz;}
// cut a section out of the 3D grid: tricubic when idle, trilinear (fast=true) while dragging; no per-pixel allocations
function samplePlane(m,vol,fr,ext,step,type,c0=0,fast=false){
  const n=Math.round(2*ext/step)+1,gs=2*ext/(n-1),[N1,N2,N3]=vol.N,g=new Float32Array(n*n),[w0,w1,w2]=wb,FO=vol.fo,FC=vol.fc;
  const o=mv(m.Mi,fr.o),e1=mv(m.Mi,fr.e1),e2=mv(m.Mi,fr.e2),ix=new Int32Array(4),iy=new Int32Array(4),iz=new Int32Array(4);
  const md=(a,b)=>{a%=b;return a<0?a+b:a;},dif=type==='diff',fo=type.startsWith('fo'),NN=N1*N2;
  for(let j=0;j<n;j++)for(let i=0;i<n;i++){
    const u=-ext+i*gs,v=-ext+j*gs,p0=(o[0]+u*e1[0]+v*e2[0])*N1,p1=(o[1]+u*e1[1]+v*e2[1])*N2,p2=(o[2]+u*e1[2]+v*e2[2])*N3;
    const b0=Math.floor(p0),b1=Math.floor(p1),b2=Math.floor(p2);let val;
    if(fast){const x0=md(b0,N1),x1=md(b0+1,N1),y0=md(b1,N2)*N1,y1=md(b1+1,N2)*N1,z0=md(b2,N3)*NN,z1=md(b2+1,N3)*NN,fx=p0-b0,fy=p1-b1,fz=p2-b2;
      val=dif?tri(FO,x0,x1,y0,y1,z0,z1,fx,fy,fz)-tri(FC,x0,x1,y0,y1,z0,z1,fx,fy,fz):tri(fo?FO:FC,x0,x1,y0,y1,z0,z1,fx,fy,fz);}
    else{cubInto(p0-b0,w0);cubInto(p1-b1,w1);cubInto(p2-b2,w2);
      for(let a=0;a<4;a++){ix[a]=md(b0+a-1,N1);iy[a]=md(b1+a-1,N2)*N1;iz[a]=md(b2+a-1,N3)*NN;}
      val=dif?cub(FO,ix,iy,iz,w0,w1,w2)-cub(FC,ix,iy,iz,w0,w1,w2):cub(fo?FO:FC,ix,iy,iz,w0,w1,w2);}
    g[j*n+i]=val+c0;}
  return{g,n,gs};
}
// all contour levels k0..k1 (step cs) in one pass; segments are collected in 3 flat buffers (positive, negative, zero)
// so that only three canvas strokes are needed: returns [{a:Float32Array,n}, ...]
function contourAll(g,n,cs,k0,k1){
  const B=[0,1,2].map(()=>({a:new Float32Array(4096),n:0})),P=new Float64Array(8);
  for(let j=0;j<n-1;j++)for(let i=0;i<n-1;i++){
    const a=g[j*n+i],b=g[j*n+i+1],c=g[(j+1)*n+i+1],d=g[(j+1)*n+i];
    let lo=a,hi=a;if(b<lo)lo=b;if(b>hi)hi=b;if(c<lo)lo=c;if(c>hi)hi=c;if(d<lo)lo=d;if(d>hi)hi=d;
    const ka=Math.max(k0,Math.ceil(lo/cs)),kb=Math.min(k1,Math.floor(hi/cs));
    for(let k=ka;k<=kb;k++){const lv=k*cs;let q=0;
      if((a>lv)!==(b>lv)){P[q++]=i+(lv-a)/(b-a);P[q++]=j;}
      if((b>lv)!==(c>lv)){P[q++]=i+1;P[q++]=j+(lv-b)/(c-b);}
      if((d>lv)!==(c>lv)){P[q++]=i+(lv-d)/(c-d);P[q++]=j+1;}
      if((a>lv)!==(d>lv)){P[q++]=i;P[q++]=j+(lv-a)/(d-a);}
      const t=B[k>0?0:k<0?1:2];
      for(let r=0;r+3<q;r+=4){if(t.n+4>t.a.length){const u=new Float32Array(t.a.length*2);u.set(t.a);t.a=u;}
        t.a[t.n++]=P[r];t.a[t.n++]=P[r+1];t.a[t.n++]=P[r+2];t.a[t.n++]=P[r+3];}}}
  return B;
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
function overlay(m,fr,ext,slab=0.4){ // atoms (symmetry copies + cell translations) near the section; follows the view, so panning stays populated
  const out=[],seen=new Set(),R=ext+2,lo=[1e9,1e9,1e9],hi=[-1e9,-1e9,-1e9];
  for(const su of[-R,R])for(const sv of[-R,R])for(const sn of[-slab,slab]){
    const f=mv(m.Mi,[0,1,2].map(i=>fr.o[i]+su*fr.e1[i]+sv*fr.e2[i]+sn*fr.n[i]));
    for(let i=0;i<3;i++){lo[i]=Math.min(lo[i],f[i]);hi[i]=Math.max(hi[i],f[i]);}}
  for(const at of m.atoms)for(const o of m.ops){
    const x0=[0,1,2].map(i=>dot(o.R[i],at.x)+o.t[i]);
    for(let tx=Math.ceil(lo[0]-x0[0]);tx<=Math.floor(hi[0]-x0[0]);tx++)for(let ty=Math.ceil(lo[1]-x0[1]);ty<=Math.floor(hi[1]-x0[1]);ty++)for(let tz=Math.ceil(lo[2]-x0[2]);tz<=Math.floor(hi[2]-x0[2]);tz++){
      const c=mv(m.M,[x0[0]+tx,x0[1]+ty,x0[2]+tz]),d=c.map((v,i)=>v-fr.o[i]),u=dot(d,fr.e1),v=dot(d,fr.e2);
      if(Math.abs(dot(d,fr.n))>=slab||Math.abs(u)>=R||Math.abs(v)>=R)continue;
      const key=Math.round(u*20)+','+Math.round(v*20);if(seen.has(key))continue;seen.add(key);out.push({u,v,el:at.el,lab:at.lab,c,dz:dot(d,fr.n)});}}
  return out;
}
