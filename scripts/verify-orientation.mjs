import assert from 'node:assert/strict';
import { readFileSync, mkdtempSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { Buffer } from 'node:buffer';
import { URL } from 'node:url';
import { BufferAttribute, Vector3 } from 'three/webgpu';
import { compileKernel } from './compile-kernel.mjs';
import { loadModel } from './load-model.mjs';
import { SoftBody } from '../src/physics/soft-body.js';
import { PHYS } from '../src/physics/constants.js';

const directory=mkdtempSync(join(tmpdir(),'jelly-orientation-'));
try {
  const wrapper=readFileSync('src/physics/soft-body-kernel.js','utf8');
  async function makeBody(reference) {
    const wasm=join(directory,reference?'reference.wasm':'indexed.wasm');
    compileKernel(wasm,['VERIFY_ORIENTATION',...(reference?['REFERENCE_ORIENTATION_SCAN']:[])]);
    // Use the real memory layout and wrapper; expose test exports only in this
    // in-memory copy. Neither test variant replaces the checked-in payload.
    const pattern=/const KERNEL_BASE64='[^']*';/;
    assert(pattern.test(wrapper)&&wrapper.includes('setCenter(center)'));
    const source=wrapper.replace(pattern,`const KERNEL_BASE64='${readFileSync(wasm).toString('base64')}';`)
      // A data: module has no directory from which to resolve wrapper imports.
      .replace("'./collision-kernel.js'",JSON.stringify(new URL('../src/physics/collision-kernel.js',import.meta.url).href))
      .replace('setCenter(center)','verification:ex, setCenter(center)');
    const {createSoftBodyKernel}=await import(`data:text/javascript;base64,${Buffer.from(source).toString('base64')}`);
    const body=new SoftBody(loadModel());
    body.kernel=createSoftBodyKernel(body);assert(body.kernel);
    for(const key of ['x','previous','candidate','velocity','contact','nodalF'])body[key]=body.kernel[key];
    body.surface.positions=body.kernel.surfacePositions;
    body.surface.geometry.setAttribute('position',new BufferAttribute(body.kernel.surfacePositions,3));
    body.surface.geometry.setAttribute('normal',new BufferAttribute(body.kernel.surfaceNormals,3));
    body.canSleep=false;
    return body;
  }
  const indexed=await makeBody(false),reference=await makeBody(true),shipped=new SoftBody(loadModel());
  assert(shipped.kernel);shipped.canSleep=false;
  assert.equal(shipped.x.buffer.byteLength,16*1024*1024,'fixed WASM memory budget is unchanged');
  const repairBodies=[indexed,reference],bodies=[...repairBodies,shipped];
  const bytes=array=>Buffer.from(array.buffer,array.byteOffset,array.byteLength);
  function compare(label,surface=false,targets=bodies) {
    for(const body of targets) {
      for(const key of ['x','previous','velocity','contact'])assert(bytes(body[key]).equals(bytes(reference[key])),`${label}: exact ${key}`);
      assert(bytes(body.kernel.meta.subarray(0,7)).equals(bytes(reference.kernel.meta.subarray(0,7))),`${label}: exact solver metadata`);
      assert.equal(body.stepFraction,reference.stepFraction);
    }
    if(surface) {
      for(const body of targets)body.updateSurface();
      for(const body of targets)for(const key of ['nodalF','surfacePositions','surfaceNormals','meta'])assert(bytes(body.kernel[key]).equals(bytes(reference.kernel[key])),`${label}: exact ${key}`);
    }
  }

  // Surface stencils at the left, right and crown: ordinary and violent motion
  // use the same compliant constraints as touch input, including shared nodes.
  const surface=indexed.surface,anchors=[];
  for(const direction of [new Vector3(-1,.3,0),new Vector3(1,.3,0),new Vector3(0,1,0)]) {
    let best=0,score=-Infinity;
    for(let i=0;i<surface.positions.length/3;i++) {
      const value=new Vector3().fromArray(surface.positions,i*3).dot(direction);
      if(value>score){score=value;best=i;}
    }
    anchors.push(best);
  }
  let indexedWork=0,referenceWork=0,repairSteps=0,maximumStepWork=0;
  const timings=[[],[]];
  for(const body of bodies)body.grabs=anchors.map(vertex=>{
    const weights=Array.from({length:4},(_,k)=>[surface.bindingIds[vertex*4+k],surface.bindingWeights[vertex*4+k]]);
    const point=new Vector3();for(const [id,w] of weights)point.addScaledVector(new Vector3().fromArray(body.x,id*3),w);
    return {weights,point,target:point.clone(),lambda:new Float64Array(3)};
  });
  const origins=indexed.grabs.map(grab=>grab.target.clone());
  for(let step=0;step<360;step++) {
    if(step===240)for(const body of bodies)body.grabs=[];
    for(let b=0;b<bodies.length;b++) {
      const body=bodies[b];
      for(let g=0;g<body.grabs.length;g++) {
        const phase=step*.39+g*2.1,amplitude=step<24?.004:.10;
        body.grabs[g].target.copy(origins[g]).add(new Vector3(Math.sin(phase)*amplitude,.065+Math.cos(phase*.83)*amplitude,Math.sin(phase*1.31)*amplitude));
      }
      const before=body.kernel.verification?.verify_evaluations()??0,start=performance.now();
      body.step(PHYS.step);if(b<2)timings[b].push(performance.now()-start);
      // WASM i32 exports are signed JS numbers; differences wrap modulo 2^32.
      const work=((body.kernel.verification?.verify_evaluations()??0)-before)>>>0;
      if(b===0){indexedWork+=work;maximumStepWork=Math.max(maximumStepWork,work);}else if(b===1)referenceWork+=work;
      assert(body.isFinite());assert(body.lastMinJacobian>=.12,'hard grabs/release preserve orientation');
      assert.equal(body.stepFraction,1,'hard grabs/release consume the full timestep');
    }
    repairSteps+=indexed.kernel.meta[2];compare(`step ${step}`,step%30===0);
  }
  assert(repairSteps>20,'trajectory exercises sustained orientation repair');
  assert(indexedWork<referenceWork*.5,'local indexing still avoids redundant full scans');
  assert(maximumStepWork<indexed.elements.length*32,'every stressed substep has bounded work, including fallback');
  assert(indexed.guardedSteps>20,'trajectory exercises admissible-motion recovery');

  // Force stalled projection and equal minima. This exercises local blending,
  // deduplication, tie order and a fresh cache after external/reset position edits.
  for(const scenario of ['collapsed cluster','flat body','collapsed body','reset']) {
    for(const body of repairBodies) {
      body.reset();body.previous.set(body.x);
      if(scenario==='collapsed cluster') {
        const ids=body.elements[body.elements.length>>1].ids,point=body.x.slice(ids[0]*3,ids[0]*3+3);
        for(const id of ids)body.x.set(point,id*3);
      } else if(scenario==='flat body')for(let i=1;i<body.x.length;i+=3)body.x[i]=.04;
      else if(scenario==='collapsed body')body.x.fill(.04);
      body.kernel.verification.verify_orientation(PHYS.floor);
      assert(body.minimumJacobian()>=.12,`${scenario}: even a collapsed candidate must produce an admissible state`);
      assert(body.isFinite(),`${scenario}: finite repair output`);
    }
    compare(scenario,true,repairBodies);
    assert.equal(indexed.minimumJacobian(),reference.minimumJacobian(),`${scenario}: exact final minimum`);
    if(scenario==='collapsed body')assert(indexed.kernel.verification.verify_fallbacks()>0,`${scenario}: stalled projection falls back locally`);
    if(scenario==='reset')assert.equal(indexed.kernel.verification.verify_evaluations(),indexed.elements.length,'valid state only scans once');
  }
  const summary=values=>({mean:values.reduce((sum,v)=>sum+v,0)/values.length,p95:[...values].sort((a,b)=>a-b)[Math.floor(values.length*.95)],max:Math.max(...values)});
  console.log('PASS — byte-identical native state, normals and bounds; hard grabs, release, stalled repair, ties and reset.',
    {repairSteps,guardedSteps:indexed.guardedSteps,maximumStepWork,indexedWork,referenceWork,indexedMs:summary(timings[0]),referenceMs:summary(timings[1])});
} finally {rmSync(directory,{recursive:true,force:true});}
