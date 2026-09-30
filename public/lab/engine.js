export const DEFAULTS={model:'sopdt',scenario:'step',kp:2,ti:3,td:0.5,enableI:true,enableD:false,gain:0.8,tau:3,tau2:0.5,delay:0.35,target:60,disturbance:-8,duration:30,band:0.5,noise:0,filter:0.1,umin:0,umax:100,antiWindup:true};
export const SPECS={kp:[0,20,0.1],ti:[0.1,20,0.1],td:[0,3,0.05],gain:[0.02,3,0.02],tau:[0.2,15,0.1],tau2:[0.1,5,0.1],delay:[0,3,0.05],target:[40,75,0.5],disturbance:[-20,20,0.5],duration:[15,90,1],band:[0.1,3,0.1],noise:[0,1,0.02],filter:[0.02,1,0.02],umin:[0,35,1],umax:[45,100,1]};
export function clean(raw={}){const c={...DEFAULTS};for(const [k,[a,b]]of Object.entries(SPECS)){const n=Number(raw[k]);if(raw[k]!==undefined&&Number.isFinite(n))c[k]=Math.max(a,Math.min(b,n));}for(const k of ['enableI','enableD','antiWindup'])if(typeof raw[k]==='boolean')c[k]=raw[k];if(['sopdt','fopdt','integrating'].includes(raw.model))c.model=raw.model;if(['step','disturbance','both'].includes(raw.scenario))c.scenario=raw.scenario;return c;}
export function units(c){return c.model==='integrating'?'%':'℃';}
export function eventTime(c){return c.scenario==='step'?1:c.scenario==='disturbance'?5:10;}
export function simulate(input,stepOverride){
  const c=clean(input),dt=stepOverride||Math.min(0.02,c.tau/80,c.tau2/40,c.filter/5),steps=Math.ceil(c.duration/dt),step=c.duration/steps;
  const initial=c.scenario==='disturbance'?c.target:50,bias=40,delaySteps=c.delay/step;
  let y=initial,act=bias,integral=0,filteredRate=0,prevMeasured=initial,ia=0,ita=0,saturation=0,peak=-Infinity,low=Infinity,lastOutside=null,stopped=false;
  const commandHistory=[],rows=[];const event=eventTime(c);
  for(let n=0;n<=steps;n++){
    const t=n*step,sp=c.scenario==='disturbance'||t>=1?c.target:50;
    const disturbed=c.scenario!=='step'&&t>=(c.scenario==='both'?10:5);
    // Fixed 0.02 min sample-and-hold noise, independent of the integration step.
    let noiseKey=(Math.floor((t+1e-9)/0.02)+29371)|0;
    noiseKey=Math.imul(noiseKey^(noiseKey>>>16),0x45d9f3b);
    noiseKey=Math.imul(noiseKey^(noiseKey>>>16),0x45d9f3b);
    const measured=y+c.noise*(((noiseKey^(noiseKey>>>16))>>>0)/4294967296*2-1);
    const e=sp-measured,rawRate=(measured-prevMeasured)/step;
    filteredRate+=(1-Math.exp(-step/c.filter))*((n?rawRate:0)-filteredRate);
    // Apply the same first-order filter at the first and subsequent samples.
    if(n===0)filteredRate=0;
    const p=c.kp*e,d=c.enableD?-c.kp*c.td*filteredRate:0;
    let raw=bias+p+integral+d;const deltaI=c.enableI?c.kp/c.ti*e*step:0;
    if(c.enableI&&(!c.antiWindup||!((raw>=c.umax&&deltaI>0)||(raw<=c.umin&&deltaI<0))))integral+=deltaI;
    if(!c.enableI)integral=0;
    raw=bias+p+integral+d;const u=Math.max(c.umin,Math.min(c.umax,raw));commandHistory.push(u);
    const q=n-delaySteps,j=Math.floor(q),fraction=q-j;const h=idx=>idx<0?bias:commandHistory[Math.min(idx,n)];const delayed=h(j)*(1-fraction)+h(j+1)*fraction;
    if(c.model!=='sopdt')act=delayed;
    const trueError=sp-y;
    if(t>=event){peak=Math.max(peak,y);low=Math.min(low,y);if(Math.abs(trueError)>c.band)lastOutside=t;if(n<steps){ia+=Math.abs(trueError)*step;ita+=(t-event)*Math.abs(trueError)*step;if(raw>c.umax||raw<c.umin)saturation+=step;}}
    rows.push({t,sp,pv:y,mv:measured,u,actual:act,p,i:integral,d,raw,load:disturbed?c.disturbance:0});
    if(!Number.isFinite(y)||!Number.isFinite(raw)||Math.abs(y)>10000||(c.model==='integrating'&&(y<0||y>100))){stopped=true;break;}
    prevMeasured=measured;
    if(c.model==='sopdt')act+=(delayed-act)*(1-Math.exp(-step/c.tau2));
    if(c.model==='integrating')y+=(c.gain*(act-bias)+(disturbed?c.disturbance:0))*step;
    else y+=(initial+c.gain*(act-bias)+(disturbed?c.disturbance:0)-y)*(1-Math.exp(-step/c.tau));
  }
  const last=rows.at(-1),recovery=lastOutside===null?event:lastOutside+step;
  const settled=!stopped&&Math.abs(last.sp-last.pv)<=c.band&&last.t-recovery>=3?Math.max(0,recovery-event):null;
  const maxDev=Math.max(Math.abs(peak-c.target),Math.abs(low-c.target));
  const over=c.target>=initial?Math.max(0,peak-c.target):Math.max(0,c.target-low);
  return {config:c,dt:step,rows,metrics:{peak,low,error:last.sp-last.pv,settled,iae:ia,itae:ita,saturation,overshoot:over,maxDeviation:maxDev,stopped,lastTime:last.t,final:last.pv},event};
}
export const PRESETS={
 p:{label:'仅比例',note:'先看窗口末偏差；增强 P 可以减小偏差，但可能增加过冲。',values:{kp:2,enableI:false,enableD:false}},
 pi:{label:'PI 消偏',note:'保持 P 不变，启用积分；比较最终偏差和恢复过程。Ti 越小，积分累积越快。',values:{kp:2,ti:3,enableI:true,enableD:false}},
 aggressive:{label:'过强积分',note:'积分累积快于对象响应，可能出现过冲、振荡及输出限幅。',values:{kp:8,ti:0.5,enableI:true,enableD:false}},
 pid:{label:'微分修正',note:'D 根据测量变化趋势修正输出；对比最高偏离与调节时间。',values:{kp:4,ti:4,td:0.5,enableI:true,enableD:true}},
 noise:{label:'测量噪声',note:'观察测量值与真实值的差异，以及 D 对输出的影响；可改变滤波时间。',values:{kp:4,ti:4,td:0.5,enableI:true,enableD:true,noise:0.25,filter:0.04}},
 disturbance:{label:'抗扰响应',note:'设定值不变，第 5 分钟施加持续冷负荷；观察最低点、恢复时间与输出。',values:{scenario:'disturbance',kp:2,ti:3,enableI:true,enableD:false}}
};
