(() => {
"use strict";

const $ = id => document.getElementById(id);
const MASS_NAMES = ["M₁","M₂","M₃","M₄","M₅"];

const DEFAULTS = {
  F0:10, fExc:2,
  M1:1, M2:.8, M3:1.2, M4:1.2, M5:1,
  K1:1000, K2:700, K3:500, K4:1200,
  B1:5, B2:3, B3:4, B4:5, B5:5,
  speed:.35, visualScale:55
};

const PARAM_IDS = [
  "F0","fExc","M1","M2","M3","M4","M5",
  "K1","K2","K3","K4","B1","B2","B3","B4","B5"
];

const state = {
  params:{...DEFAULTS},
  M:null,C:null,K:null,Minv:null,
  q:new Float64Array(5),
  v:new Float64Array(5),
  a:new Float64Array(5),
  t:0,
  running:false,
  accumulator:0,
  fixedDt:1/1600,
  sampleDt:1/240,
  sampleAccumulator:0,
  history:[],
  historyLimit:4800,
  naturalFrequencies:[],
  force:0,
  fps:0,
  fpsCount:0,
  fpsClock:0,
  rafTime:performance.now(),
  applyQueued:false
};

// ------------------------------------------------------------
// Matemáticas básicas
// ------------------------------------------------------------
function zeros(n,m){return Array.from({length:n},()=>Array(m).fill(0))}
function diag(v){const A=zeros(v.length,v.length);for(let i=0;i<v.length;i++)A[i][i]=v[i];return A}
function clone(A){return A.map(r=>r.slice())}
function dot(a,b){let s=0;for(let i=0;i<a.length;i++)s+=a[i]*b[i];return s}
function maxAbs(v){let m=0;for(const x of v)m=Math.max(m,Math.abs(x));return m}
function clamp(x,a,b){return Math.max(a,Math.min(b,x))}
function matVec(A,x){const y=new Float64Array(A.length);for(let i=0;i<A.length;i++){let s=0;for(let j=0;j<x.length;j++)s+=A[i][j]*x[j];y[i]=s}return y}

function addPair(A,i,j,value){
  A[i][i]+=value; A[j][j]+=value;
  A[i][j]-=value; A[j][i]-=value;
}
function addGround(A,i,value){A[i][i]+=value}

function human(x,dec=4){
  if(!Number.isFinite(x)) return "—";
  if(Math.abs(x)<1e-12) return "0.0000";
  const ax=Math.abs(x);
  let digits=dec;
  if(ax>=1000) digits=2;
  else if(ax>=100) digits=3;
  else if(ax>=10) digits=3;
  else if(ax>=1) digits=4;
  else if(ax>=.1) digits=5;
  else if(ax>=.01) digits=6;
  else digits=7;
  let s=Number(x).toFixed(digits);
  s=s.replace(/\.?0+$/,"");
  return s===" -0" || s==="-0" ? "0" : s;
}

function humanUnit(x,unit){
  return `${human(x)} ${unit}`;
}

function matrixRange(A,r0,r1,c0,c1){
  const out=[];
  for(let i=r0;i<r1;i++)out.push(A[i].slice(c0,c1));
  return out;
}

// ------------------------------------------------------------
// Modelo físico
//
// q = [u1,u2,u3,u4,u5]
//
// PARED-F(t)-M1
// M1-B1-piso
// M1-K1-nodo
// nodo -> K2-M2-B2 -> nodo de reunión -> M4
// nodo -> K3-M3-B3 -> nodo de reunión -> M4
// M4-B4-M5
// M5-B5-piso
// M5-K4-pared
//
// El nodo entre K1,K2,K3 no tiene masa.
// ------------------------------------------------------------
function buildModel(p){
  const M=diag([p.M1,p.M2,p.M3,p.M4,p.M5]);

  const C=zeros(5,5);
  addGround(C,0,p.B1);
  addPair(C,1,3,p.B2);
  addPair(C,2,3,p.B3);
  addPair(C,3,4,p.B4);
  addGround(C,4,p.B5);

  // Kfull incluye el nodo sin masa en índice 5.
  const Kfull=zeros(6,6);
  addPair(Kfull,0,5,p.K1);
  addPair(Kfull,5,1,p.K2);
  addPair(Kfull,5,2,p.K3);
  addGround(Kfull,4,p.K4);

  // Condensación estática del nodo ideal:
  // K_eff = Kqq - Kqx*Kxx^-1*Kxq
  const Kqq=matrixRange(Kfull,0,5,0,5);
  const Kqx=Kfull.slice(0,5).map(r=>[r[5]]);
  const Kxq=[Kfull[5].slice(0,5)];
  const Kxx=Kfull[5][5];

  const K=zeros(5,5);
  if(Kxx>0){
    for(let i=0;i<5;i++){
      for(let j=0;j<5;j++){
        K[i][j]=Kqq[i][j]-Kqx[i][0]*Kxq[0][j]/Kxx;
      }
    }
  }else{
    for(let i=0;i<5;i++)for(let j=0;j<5;j++)K[i][j]=Kqq[i][j];
  }

  const Minv=[1/p.M1,1/p.M2,1/p.M3,1/p.M4,1/p.M5];
  return {M,C,K,Minv};
}

// ------------------------------------------------------------
// Eigenproblema Kφ = ω²Mφ
// ------------------------------------------------------------
function jacobi(Ain,maxIter=300,tol=1e-11){
  const A=clone(Ain), n=A.length, V=diag(Array(n).fill(1));
  for(let it=0;it<maxIter;it++){
    let p=0,q=1,largest=0;
    for(let i=0;i<n;i++)for(let j=i+1;j<n;j++){
      const x=Math.abs(A[i][j]);
      if(x>largest){largest=x;p=i;q=j}
    }
    if(largest<tol)break;

    const tau=(A[q][q]-A[p][p])/(2*A[p][q]);
    const sgn=tau>=0?1:-1;
    const t=sgn/(Math.abs(tau)+Math.sqrt(1+tau*tau));
    const c=1/Math.sqrt(1+t*t), s=t*c;

    const App=A[p][p], Aqq=A[q][q], Apq=A[p][q];
    A[p][p]=App-t*Apq;
    A[q][q]=Aqq+t*Apq;
    A[p][q]=0;A[q][p]=0;

    for(let k=0;k<n;k++){
      if(k===p||k===q)continue;
      const Akp=A[k][p], Akq=A[k][q];
      A[k][p]=c*Akp-s*Akq;
      A[p][k]=A[k][p];
      A[k][q]=s*Akp+c*Akq;
      A[q][k]=A[k][q];
    }
    for(let k=0;k<n;k++){
      const Vkp=V[k][p],Vkq=V[k][q];
      V[k][p]=c*Vkp-s*Vkq;
      V[k][q]=s*Vkp+c*Vkq;
    }
  }
  const vals=Array.from({length:n},(_,i)=>A[i][i]);
  const order=vals.map((_,i)=>i).sort((a,b)=>vals[a]-vals[b]);
  return {values:order.map(i=>vals[i]), vectors:order.map(i=>V.map(r=>r[i]))};
}

function modalAnalysis(M,K){
  const n=5;
  const s=Array.from({length:n},(_,i)=>1/Math.sqrt(M[i][i]));
  const A=zeros(n,n);
  for(let i=0;i<n;i++)for(let j=0;j<n;j++)A[i][j]=s[i]*K[i][j]*s[j];

  const eig=jacobi(A);
  const freqs=[],modes=[];
  for(let m=0;m<n;m++){
    const lambda=Math.max(0,eig.values[m]);
    const f=Math.sqrt(lambda)/(2*Math.PI);
    let phi=eig.vectors[m].map((z,i)=>s[i]*z);
    const scale=maxAbs(phi)||1;
    phi=phi.map(z=>z/scale);
    freqs.push(f);modes.push(phi);
  }
  return {freqs,modes};
}

// ------------------------------------------------------------
// Fuente + aceleración
// ------------------------------------------------------------
function forceAt(t,p){
  if(p.F0===0)return 0;
  return p.F0*Math.sin(2*Math.PI*p.fExc*t);
}

function nodeDisplacement(q,p){
  const s=p.K1+p.K2+p.K3;
  if(s<=0)return (q[0]+q[1]+q[2])/3;
  return (p.K1*q[0]+p.K2*q[1]+p.K3*q[2])/s;
}

function nodeVelocity(v,p){
  const s=p.K1+p.K2+p.K3;
  if(s<=0)return (v[0]+v[1]+v[2])/3;
  return (p.K1*v[0]+p.K2*v[1]+p.K3*v[2])/s;
}

function acceleration(q,v,t){
  const p=state.params;
  const f=new Float64Array(5);
  f[0]=forceAt(t,p);state.force=f[0];

  const cv=matVec(state.C,v);
  const kq=matVec(state.K,q);
  const a=new Float64Array(5);
  for(let i=0;i<5;i++)a[i]=state.Minv[i]*(f[i]-cv[i]-kq[i]);
  return a;
}

// ------------------------------------------------------------
// RK4
// ------------------------------------------------------------
function rk4(dt){
  const q0=state.q,v0=state.v,t0=state.t;
  const a1=acceleration(q0,v0,t0);

  const q2=new Float64Array(5),v2=new Float64Array(5);
  for(let i=0;i<5;i++){q2[i]=q0[i]+dt*.5*v0[i];v2[i]=v0[i]+dt*.5*a1[i]}
  const a2=acceleration(q2,v2,t0+.5*dt);

  const q3=new Float64Array(5),v3=new Float64Array(5);
  for(let i=0;i<5;i++){q3[i]=q0[i]+dt*.5*v2[i];v3[i]=v0[i]+dt*.5*a2[i]}
  const a3=acceleration(q3,v3,t0+.5*dt);

  const q4=new Float64Array(5),v4=new Float64Array(5);
  for(let i=0;i<5;i++){q4[i]=q0[i]+dt*v3[i];v4[i]=v0[i]+dt*a3[i]}
  const a4=acceleration(q4,v4,t0+dt);

  for(let i=0;i<5;i++){
    state.q[i]+=dt*(v0[i]+2*v2[i]+2*v3[i]+v4[i])/6;
    state.v[i]+=dt*(a1[i]+2*a2[i]+2*a3[i]+a4[i])/6;
  }
  state.t+=dt;
  state.a=acceleration(state.q,state.v,state.t);
}

// ------------------------------------------------------------
// Controles: slider <-> campo numérico realmente vinculados
// ------------------------------------------------------------
function niceMax(x){
  if(x<1)return 1;
  if(x<10)return Math.ceil(x*2);
  if(x<100)return Math.ceil(x*1.5/10)*10;
  if(x<1000)return Math.ceil(x*1.3/100)*100;
  return Math.ceil(x*1.2/1000)*1000;
}

function syncControl(id,value){
  const r=$(id),n=$(id+"num"),o=$(id+"Out");
  if(!r||!n)return;
  const v=Number(value);
  if(!Number.isFinite(v))return;

  if(v>Number(r.max))r.max=String(niceMax(v));
  if(v<Number(r.min))r.min=String(v);

  r.value=String(v);
  n.value=String(v);

  if(o){
    const digits=id.startsWith("K")?0:2;
    o.textContent=v.toFixed(digits);
  }
}

function scheduleApply(){
  if(state.applyQueued)return;
  state.applyQueued=true;
  requestAnimationFrame(()=>{
    state.applyQueued=false;
    applyParameters(false);
  });
}

function setupControls(){
  for(const id of PARAM_IDS){
    const r=$(id),n=$(id+"num");
    const from=(source)=>{
      const v=Number(source==="range"?r.value:n.value);
      if(!Number.isFinite(v))return;
      syncControl(id,v);
      state.params[id]=v;
      if(state.M)scheduleApply();
    };
    r.addEventListener("input",()=>from("range"));
    n.addEventListener("input",()=>from("number"));
    r.addEventListener("change",()=>applyParameters(false));
    n.addEventListener("change",()=>applyParameters(false));
    syncControl(id,DEFAULTS[id]);
  }

  $("speed").addEventListener("input",()=>{
    state.params.speed=Number($("speed").value);
    $("speedOut").textContent=state.params.speed.toFixed(2)+"×";
  });
  $("visualScale").addEventListener("input",()=>{
    state.params.visualScale=Number($("visualScale").value);
    $("visualScaleOut").textContent=state.params.visualScale.toFixed(0)+"×";
    drawMechanism();
  });
}

function readParameters(){
  const p={...state.params};
  for(const id of PARAM_IDS)p[id]=Number($(id+"num").value);
  p.speed=Number($("speed").value);
  p.visualScale=Number($("visualScale").value);
  return p;
}

function validate(p){
  const invalid=[];
  for(const id of PARAM_IDS)if(!Number.isFinite(p[id]))invalid.push(id);
  for(const id of ["M1","M2","M3","M4","M5"])if(p[id]<=0)invalid.push(id);
  if(p.K1+p.K2+p.K3<=0)invalid.push("K1+K2+K3");
  for(const id of ["F0","fExc","K1","K2","K3","K4","B1","B2","B3","B4","B5"])if(p[id]<0)invalid.push(id);
  return [...new Set(invalid)];
}

function resetState(){
  state.q.fill(0);state.v.fill(0);state.a.fill(0);
  state.t=0;state.accumulator=0;state.sampleAccumulator=0;
  state.history=[];state.force=0;
}

function chooseDt(){
  const maxF=Math.max(...state.naturalFrequencies,1);
  state.fixedDt=Math.min(1/1600,1/(Math.max(20,maxF*45)));
  state.sampleDt=Math.max(state.fixedDt*2,1/240);
}

function applyParameters(reset){
  const p=readParameters();
  const bad=validate(p);
  if(bad.length){
    $("errorBox").textContent="Parámetros inválidos: "+bad.join(", ");
    return false;
  }

  $("errorBox").textContent="";
  state.params=p;
  const model=buildModel(p);
  state.M=model.M;state.C=model.C;state.K=model.K;state.Minv=model.Minv;

  const modal=modalAnalysis(state.M,state.K);
  state.naturalFrequencies=modal.freqs;
  chooseDt();

  if(reset)resetState();
  else state.a=acceleration(state.q,state.v,state.t);

  updateTelemetry();
  drawAll();
  return true;
}

function resetDefaults(){
  for(const id of PARAM_IDS)syncControl(id,DEFAULTS[id]);
  $("speed").value=DEFAULTS.speed;
  $("visualScale").value=DEFAULTS.visualScale;
  $("speedOut").textContent=DEFAULTS.speed.toFixed(2)+"×";
  $("visualScaleOut").textContent=DEFAULTS.visualScale.toFixed(0)+"×";
  state.running=false;
  applyParameters(true);
  updateEngineUI();
}

function updateEngineUI(){
  const box=$("engineState");
  const text=$("engineStateText");
  if(state.running){
    box.classList.add("running");text.textContent="EJECUTANDO";$("runBtn").textContent="Ⅱ Pausar";
  }else{
    box.classList.remove("running");text.textContent="DETENIDO";$("runBtn").textContent="▶ Iniciar";
  }
}

// ------------------------------------------------------------
// Gráficos / Canvas
// ------------------------------------------------------------
function canvasSetup(canvas){
  const dpr=Math.min(window.devicePixelRatio||1,2);
  const rect=canvas.getBoundingClientRect();
  canvas.width=Math.max(10,Math.floor(rect.width*dpr));
  canvas.height=Math.max(10,Math.floor(rect.height*dpr));
  const ctx=canvas.getContext("2d");
  ctx.setTransform(dpr,0,0,dpr,0,0);
  return {ctx,w:rect.width,h:rect.height};
}

function line(ctx,x1,y1,x2,y2,color,width=2,dash=[]){
  ctx.save();ctx.strokeStyle=color;ctx.lineWidth=width;ctx.setLineDash(dash);
  ctx.beginPath();ctx.moveTo(x1,y1);ctx.lineTo(x2,y2);ctx.stroke();ctx.restore();
}

function wall(ctx,x,y0,y1){
  line(ctx,x,y0,x,y1,"#dce4ea",5);
  for(let y=y0;y<y1;y+=11)line(ctx,x-9,y,x,y+8,"#6e7e8c",1);
}

function floor(ctx,x0,x1,y){
  line(ctx,x0,y,x1,y,"#dce4ea",3);
  for(let x=x0;x<x1;x+=12)line(ctx,x,y,x-7,y+7,"#5f707e",1);
}

function spring(ctx,x1,x2,y,label){
  if(x2-x1<20)return;
  const lead=14, amp=6, coils=8;
  const length=Math.max(5,x2-x1-2*lead);
  const step=length/(coils*2);
  ctx.save();ctx.strokeStyle=state.params? "#9de393":"#9de393";ctx.lineWidth=2;ctx.beginPath();
  ctx.moveTo(x1,y);ctx.lineTo(x1+lead,y);
  let x=x1+lead,sgn=-1;
  for(let i=0;i<coils*2;i++){x+=step;ctx.lineTo(x,y+sgn*amp);sgn*=-1}
  ctx.lineTo(x2-lead,y);ctx.lineTo(x2,y);ctx.stroke();
  ctx.fillStyle="#7b8996";ctx.font="9px system-ui";ctx.textAlign="center";ctx.fillText(label,(x1+x2)/2,y-12);
  ctx.restore();
}

function damperH(ctx,x1,x2,y,label){
  const mid=(x1+x2)/2;
  ctx.save();ctx.strokeStyle="#c5a3ff";ctx.lineWidth=2;
  line(ctx,x1,y,mid-18,y,"#c5a3ff",2);
  ctx.strokeRect(mid-18,y-8,36,16);
  line(ctx,mid-8,y-8,mid-8,y+8,"#c5a3ff",2);
  line(ctx,mid+18,y,x2,y,"#c5a3ff",2);
  ctx.fillStyle="#7b8996";ctx.font="9px system-ui";ctx.textAlign="center";ctx.fillText(label,mid,y-12);
  ctx.restore();
}

function damperV(ctx,x,y1,y2,label){
  const mid=(y1+y2)/2;
  ctx.save();ctx.strokeStyle="#c5a3ff";ctx.lineWidth=2;
  line(ctx,x,y1,x,mid-16,"#c5a3ff",2);
  ctx.strokeRect(x-8,mid-16,16,32);
  line(ctx,x-8,mid-7,x+8,mid-7,"#c5a3ff",2);
  line(ctx,x,mid+16,x,y2,"#c5a3ff",2);
  ctx.fillStyle="#7b8996";ctx.font="9px system-ui";ctx.textAlign="left";ctx.fillText(label,x+11,mid+3);
  ctx.restore();
}

function mass(ctx,x,y,label,secondary=false){
  const w=62,h=44;
  ctx.save();
  ctx.fillStyle=secondary?"#182129":"#14283a";
  ctx.strokeStyle=secondary?"#8696a3":"#79baff";
  ctx.lineWidth=2;
  ctx.beginPath();ctx.roundRect(x-w/2,y-h/2,w,h,7);ctx.fill();ctx.stroke();
  ctx.fillStyle="#edf3f7";ctx.font="700 12px system-ui";ctx.textAlign="center";ctx.textBaseline="middle";
  ctx.fillText(label,x,y);ctx.restore();
}

function actuator(ctx,x1,x2,y){
  const mid=(x1+x2)/2;
  ctx.save();
  ctx.strokeStyle="#ffd071";ctx.lineWidth=2;
  line(ctx,x1,y,mid-21,y,"#ffd071",2);
  ctx.strokeRect(mid-21,y-10,42,20);
  line(ctx,mid-9,y-10,mid-9,y+10,"#ffd071",2);
  line(ctx,mid+21,y,x2,y,"#ffd071",2);
  ctx.fillStyle="#ffd071";ctx.font="700 10px system-ui";ctx.textAlign="center";ctx.fillText("F(t)",mid,y-17);
  ctx.restore();
}

function arrow(ctx,x,y,dx,color="#79baff"){
  if(Math.abs(dx)<2)return;
  const tip=x+dx,sgn=dx>=0?1:-1;
  line(ctx,x,y,tip,y,color,2);
  ctx.save();ctx.fillStyle=color;ctx.beginPath();ctx.moveTo(tip,y);ctx.lineTo(tip-sgn*7,y-4);ctx.lineTo(tip-sgn*7,y+4);ctx.closePath();ctx.fill();ctx.restore();
}

function drawMechanism(){
  const {ctx,w,h}=canvasSetup($("mechanismCanvas"));
  ctx.clearRect(0,0,w,h);

  const p=state.params, left=45,right=w-45,usable=right-left;
  const cy=h*.48;
  const scale=p.visualScale;

  const x1=left+usable*.16+state.q[0]*scale;
  const xSplit=left+usable*.34+nodeDisplacementSafe()*scale;
  const xBranch=left+usable*.53;
  const xJoin=left+usable*.70;
  const x4=left+usable*.84+state.q[3]*scale;
  const x5=left+usable*.84+state.q[4]*scale;

  const yTop=cy-84;
  const yBottom=cy+84;
  const y5=cy+76;
  const mw=62,mh=44;

  wall(ctx,left,cy-145,cy+158);
  wall(ctx,right,cy-145,cy+158);
  floor(ctx,left,right,cy+145);

  ctx.fillStyle="#4e5d69";ctx.font="9px system-ui";ctx.textAlign="left";
  ctx.fillText("desplazamiento positivo →",left,18);

  // F(t) - M1
  actuator(ctx,left+10,x1-mw/2,cy);
  mass(ctx,x1,cy,"M₁");
  damperV(ctx,x1,cy+mh/2,cy+145,"B₁");

  // K1
  spring(ctx,x1+mw/2,xSplit-4,cy,"K₁");

  // nodo separación
  line(ctx,xSplit,yTop,xSplit,yBottom,"#79baff",2);
  ctx.fillStyle="#79baff";ctx.beginPath();ctx.arc(xSplit,cy,5,0,Math.PI*2);ctx.fill();

  // rama superior K2 - M2 - B2
  line(ctx,xSplit,yTop,xSplit+12,yTop,"#79baff",2);
  spring(ctx,xSplit+12,xBranch-mw/2,yTop,"K₂");
  mass(ctx,xBranch,yTop,"M₂");
  damperH(ctx,xBranch+mw/2,xJoin,yTop,"B₂");

  // rama inferior K3 - M3 - B3
  line(ctx,xSplit,yBottom,xSplit+12,yBottom,"#79baff",2);
  spring(ctx,xSplit+12,xBranch-mw/2,yBottom,"K₃");
  mass(ctx,xBranch,yBottom,"M₃");
  damperH(ctx,xBranch+mw/2,xJoin,yBottom,"B₃");

  // colector de las dos ramas
  line(ctx,xJoin,yTop,xJoin,yBottom,"#79baff",2);
  ctx.fillStyle="#79baff";
  for(const y of [yTop,yBottom,cy]){ctx.beginPath();ctx.arc(xJoin,y,y===cy?5:4,0,Math.PI*2);ctx.fill()}
  line(ctx,xJoin,cy,x4-mw/2,cy,"#79baff",2);

  // M4 arriba y M5 debajo
  mass(ctx,x4,cy,"M₄");
  mass(ctx,x5,y5,"M₅",true);

  // B4 entre M4 y M5
  damperV(ctx,x4,cy+mh/2,y5-mh/2,"B₄");

  // B5 de M5 al piso
  damperV(ctx,x5,y5+mh/2,cy+145,"B₅");

  // K4 de M5 a pared derecha
  spring(ctx,x5+mw/2,right-4,y5,"K₄");

  const vmax=Math.max(maxAbs(state.v),1e-9);
  arrow(ctx,x1,cy-32,clamp(state.v[0]/vmax*42,-42,42));
  arrow(ctx,xBranch,yTop-32,clamp(state.v[1]/vmax*42,-42,42));
  arrow(ctx,xBranch,yBottom-32,clamp(state.v[2]/vmax*42,-42,42));
  arrow(ctx,x4,cy-32,clamp(state.v[3]/vmax*42,-42,42));
  arrow(ctx,x5,y5+32,clamp(state.v[4]/vmax*42,-42,42));

  ctx.fillStyle="#647584";ctx.font="9px system-ui";ctx.textAlign="center";
  ctx.fillText("ramas suspendidas · sin conexión al piso",xJoin,cy+118);

  ctx.textAlign="left";
  ctx.fillStyle="#7b8996";
  ctx.fillText("verde = resorte   violeta = fricción viscosa   ámbar = fuerza   azul = conexión",left,h-16);
}

function nodeDisplacementSafe(){
  return nodeDisplacement(state.q,state.params);
}

// ------------------------------------------------------------
// Gráfico temporal
// ------------------------------------------------------------
const COLORS=["#79baff","#9de393","#c5a3ff","#ffd071","#ff9494"];

function drawTimeChart(){
  const {ctx,w,h}=canvasSetup($("timeChart"));
  ctx.clearRect(0,0,w,h);
  const pad={l:54,r:16,t:22,b:34};
  ctx.fillStyle="#0a1118";ctx.fillRect(0,0,w,h);

  for(let i=0;i<5;i++){
    const y=pad.t+(h-pad.t-pad.b)*i/4;
    line(ctx,pad.l,y,w-pad.r,y,"#18232c",1);
  }
  for(let i=0;i<6;i++){
    const x=pad.l+(w-pad.l-pad.r)*i/5;
    line(ctx,x,pad.t,x,h-pad.b,"#18232c",1);
  }

  const qty=$("timeQuantity").value;
  const wnd=Number($("timeWindow").value);
  const title=qty==="velocity"?"Velocidad":qty==="displacement"?"Desplazamiento":"Aceleración";
  const unit=qty==="velocity"?"m/s":qty==="displacement"?"m":"m/s²";
  ctx.fillStyle="#7b8996";ctx.font="10px system-ui";ctx.textAlign="left";
  ctx.fillText(`${title} · ${unit}`,pad.l,13);

  if(state.history.length<2)return;
  const end=state.history[state.history.length-1].t;
  const start=Math.max(0,end-wnd);
  const data=state.history.filter(r=>r.t>=start);
  if(data.length<2)return;

  const vals=r=>qty==="velocity"?r.v:qty==="displacement"?r.q:r.a;
  let ymax=1e-12;
  for(const r of data)for(const z of vals(r))ymax=Math.max(ymax,Math.abs(z));
  ymax*=1.08;

  for(let s=0;s<5;s++){
    ctx.save();ctx.strokeStyle=COLORS[s];ctx.lineWidth=1.5;ctx.beginPath();
    for(let i=0;i<data.length;i++){
      const x=pad.l+(data[i].t-start)/(end-start||1)*(w-pad.l-pad.r);
      const y=h-pad.b-(vals(data[i])[s]+ymax)/(2*ymax)*(h-pad.t-pad.b);
      i===0?ctx.moveTo(x,y):ctx.lineTo(x,y);
    }
    ctx.stroke();ctx.restore();
  }

  ctx.fillStyle="#697886";ctx.font="8px system-ui";
  ctx.fillText(human(ymax,4),5,pad.t+3);
  ctx.fillText(human(-ymax,4),5,h-pad.b);
  ctx.fillText(start.toFixed(2)+" s",pad.l,h-9);
  ctx.textAlign="right";ctx.fillText(end.toFixed(2)+" s",w-pad.r,h-9);
  ctx.textAlign="left";

  for(let i=0;i<5;i++){
    ctx.fillStyle=COLORS[i];
    ctx.fillText(MASS_NAMES[i],pad.l+i*46,13);
  }
}

// ------------------------------------------------------------
// UI de estado, modos y telemetría
// ------------------------------------------------------------
function updateLive(){
  $("liveClock").textContent=`t = ${state.t.toFixed(3)} s`;
  $("forceNow").textContent=humanUnit(state.force,"N");

  const mf=Math.max(Math.abs(state.params.F0),1e-9);
  $("forceFill").style.width=(clamp(Math.abs(state.force)/mf,0,1)*100).toFixed(1)+"%";

  const vmax=Math.max(maxAbs(state.v),1e-12);
  $("velocityGrid").innerHTML="";
  for(let i=0;i<5;i++){
    const el=document.createElement("div");
    el.className="velocity-item";
    const pct=clamp(Math.abs(state.v[i])/vmax,0,1)*100;
    el.innerHTML=`
      <div class="v-label">${MASS_NAMES[i]}</div>
      <div class="v-value">${human(state.v[i],6)}</div>
      <div class="v-unit">m/s</div>
      <div class="v-bar"><div class="v-fill" style="width:${pct}%"></div></div>
    `;
    $("velocityGrid").appendChild(el);
  }

  $("dtReadout").textContent=(state.fixedDt*1000).toFixed(3)+" ms";
  $("fpsReadout").textContent=state.fps?state.fps.toFixed(0):"—";
  $("samplesReadout").textContent=String(state.history.length);

  let avgSpeed=0;
  let speedSamples=0;
  for(const r of state.history){
    for(let i=0;i<5;i++){
      avgSpeed+=Math.abs(r.v[i]);
      speedSamples++;
    }
  }
  avgSpeed=speedSamples?avgSpeed/speedSamples:0;
  $("avgSpeedReadout").textContent=human(avgSpeed,5)+" m/s";

  const excitation=state.params.fExc;
  let nearest=Infinity;
  for(const fn of state.naturalFrequencies)if(fn>0)nearest=Math.min(nearest,Math.abs(fn-excitation));
  $("regimeReadout").textContent=nearest<Math.max(.1,excitation*.02)?"cercana a modo":"forzada";

}

function updateTelemetry(){
  const peaks=Array(5).fill(0);
  if(state.history.length)for(const r of state.history)for(let i=0;i<5;i++)peaks[i]=Math.max(peaks[i],Math.abs(r.v[i]));

  $("telemetryBody").innerHTML="";
  for(let i=0;i<5;i++){
    const tr=document.createElement("tr");
    tr.innerHTML=`
      <td><b>${MASS_NAMES[i]}</b></td>
      <td>${human(state.q[i],6)}</td>
      <td>${human(state.v[i],6)}</td>
      <td>${human(state.a[i],6)}</td>
      <td>${human(peaks[i],6)}</td>`;
    $("telemetryBody").appendChild(tr);
  }
}

// ------------------------------------------------------------
// Exportación
// ------------------------------------------------------------
function downloadCSV(){
  if(!state.history.length){
    $("errorBox").textContent="No hay datos almacenados todavía.";
    return;
  }
  const lines=["t_s,u1_m,u2_m,u3_m,u4_m,u5_m,v1_mps,v2_mps,v3_mps,v4_mps,v5_mps"];
  for(const r of state.history){
    lines.push([r.t,...r.q,...r.v].map(x=>String(x).replace(",",".")).join(","));
  }
  const blob=new Blob([lines.join("\n")],{type:"text/csv;charset=utf-8"});
  const url=URL.createObjectURL(blob),a=document.createElement("a");
  a.href=url;a.download="simulacion_mecanica_5GDL.csv";
  document.body.appendChild(a);a.click();a.remove();URL.revokeObjectURL(url);
}

// ------------------------------------------------------------
// Loop
// ------------------------------------------------------------
function recordSample(){
  state.history.push({t:state.t,q:Array.from(state.q),v:Array.from(state.v)});
  if(state.history.length>state.historyLimit)state.history.shift();
}

function loop(now){
  const elapsed=Math.min((now-state.rafTime)/1000,.05);
  state.rafTime=now;

  state.fpsCount++;state.fpsClock+=elapsed;
  if(state.fpsClock>=.5){
    state.fps=state.fpsCount/state.fpsClock;
    state.fpsCount=0;state.fpsClock=0;
  }

  if(state.running){
    state.accumulator+=elapsed*state.params.speed;
    let guard=0;
    while(state.accumulator>=state.fixedDt && guard<320){
      rk4(state.fixedDt);

      state.accumulator-=state.fixedDt;
      state.sampleAccumulator+=state.fixedDt;

      if(state.sampleAccumulator>=state.sampleDt){
        recordSample();
        state.sampleAccumulator=0;
      }
      guard++;
    }
    if(guard>=320)state.accumulator=0;
  }

  updateLive();
  updateTelemetry();
  drawMechanism();
  drawTimeChart();

  requestAnimationFrame(loop);
}

function drawAll(){
  drawMechanism();
  drawTimeChart();
}

// ------------------------------------------------------------
// Eventos
// ------------------------------------------------------------
$("runBtn").addEventListener("click",()=>{
  state.running=!state.running;
  updateEngineUI();
});
$("resetDefaultsBtn").addEventListener("click",resetDefaults);
$("applyBtn").addEventListener("click",()=>applyParameters(false));
$("zeroBtn").addEventListener("click",()=>{
  resetState();
  updateLive();updateTelemetry();drawAll();
});
$("timeQuantity").addEventListener("change",drawTimeChart);
$("timeWindow").addEventListener("change",drawTimeChart);
$("exportBtn").addEventListener("click",downloadCSV);
window.addEventListener("resize",drawAll);

// ------------------------------------------------------------
// Arranque
// ------------------------------------------------------------
setupControls();
resetDefaults();
updateEngineUI();
requestAnimationFrame(loop);

})();
