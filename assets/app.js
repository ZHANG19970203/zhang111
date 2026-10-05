/* Original mathematics and animation binding implementation. MIT. */
'use strict';
window.SM={
 clamp:(x,a,b)=>Math.max(a,Math.min(b,x)),
 add:(a,b)=>a.map((v,i)=>v+b[i]),sub:(a,b)=>a.map((v,i)=>v-b[i]),mul:(a,s)=>a.map(v=>v*s),
 dot:(a,b)=>a[0]*b[0]+a[1]*b[1]+a[2]*b[2],cross:(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]],
 norm:a=>{const l=Math.hypot(...a);return l>1e-12?a.map(v=>v/l):[1,0,0]},
 ident:()=>new Float32Array([1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1]),
 matmul:(a,b)=>{const c=new Float32Array(16);for(let j=0;j<4;j++)for(let i=0;i<4;i++)for(let k=0;k<4;k++)c[4*j+i]+=a[4*k+i]*b[4*j+k];return c},
 point:(m,p)=>[m[0]*p[0]+m[4]*p[1]+m[8]*p[2]+m[12],m[1]*p[0]+m[5]*p[1]+m[9]*p[2]+m[13],m[2]*p[0]+m[6]*p[1]+m[10]*p[2]+m[14]],
 vec:(m,p)=>[m[0]*p[0]+m[4]*p[1]+m[8]*p[2],m[1]*p[0]+m[5]*p[1]+m[9]*p[2],m[2]*p[0]+m[6]*p[1]+m[10]*p[2]],
 rotate(axis,radians,center=[0,0,0]){const [x,y,z]=this.norm(axis),c=Math.cos(radians),s=Math.sin(radians),t=1-c;const m=new Float32Array([t*x*x+c,t*x*y+s*z,t*x*z-s*y,0,t*x*y-s*z,t*y*y+c,t*y*z+s*x,0,t*x*z+s*y,t*y*z-s*x,t*z*z+c,0,0,0,0,1]);const v=this.vec(m,center);for(let i=0;i<3;i++)m[12+i]=center[i]-v[i];return m},
 fromTo(a,b){a=this.norm(a);b=this.norm(b);let dot=this.clamp(this.dot(a,b),-1,1),cross=this.cross(a,b);if(dot>.9999999)return this.ident();if(dot<-.9999999)cross=this.cross(a,Math.abs(a[0])<.8?[1,0,0]:[0,1,0]);return this.rotate(cross,Math.acos(dot))},
 mat3(m){return new Float32Array([m[0],m[1],m[2],m[4],m[5],m[6],m[8],m[9],m[10]])},
 perspective(fov,aspect,near,far){const f=1/Math.tan(fov/2);return new Float32Array([f/aspect,0,0,0,0,f,0,0,0,0,(far+near)/(near-far),-1,0,0,2*far*near/(near-far),0])},
 lookAt(eye,target){const z=this.norm(this.sub(eye,target)),x=this.norm(this.cross([0,1,0],z)),y=this.cross(z,x);return new Float32Array([x[0],y[0],z[0],0,x[1],y[1],z[1],0,x[2],y[2],z[2],0,-this.dot(x,eye),-this.dot(y,eye),-this.dot(z,eye),1])},
};

/* Source-MJCF prescribed kinematics. Not a physics/activation solver.
   Joint ranges are model constraints, NOT population normal ROM. */
'use strict';
window.NativeKinematics=class {
 constructor(model){this.model=model;this.jmap=Object.fromEntries(model.joints.map((j,i)=>[j.name,i]));this.bmap=Object.fromEntries(model.bodies.map((b,i)=>[b.name,i]));this.smap=Object.fromEntries(model.sites.map((s,i)=>[s.name,i]));this.tmap=Object.fromEntries(model.tendons.map((t,i)=>[t.name,i]));this.q=model.joints.map(j=>j.ref);this.world=SM.rotate([1,0,0],-Math.PI/2);this.pose({});}
 values(controls={}){const q=this.model.joints.map(j=>j.ref);for(const [k,v] of Object.entries(controls)){const i=this.jmap[k];if(i===undefined||!Number.isFinite(v))continue;const j=this.model.joints[i];q[i]=SM.clamp(v,...j.range);}const eq=this.model.equalities;for(let round=0;round<4;round++){for(const e of eq){const s=e.source<0?0:q[e.source]-this.model.joints[e.source].ref;let val=0;for(let p=e.poly.length-1;p>=0;p--)val=val*s+e.poly[p];q[e.target]=this.model.joints[e.target].ref+val;}}return q;}
 pose(controls={},rootTilt=0){this.controls={...controls};this.q=this.values(controls);this.transforms=[];this.jointFrames=[];const world=SM.matmul(this.world,SM.rotate([0,1,0],rootTilt));
  this.model.bodies.forEach((b,i)=>{let T=SM.matmul(b.parent<0?world:this.transforms[b.parent],b.base);for(const ji of b.joints){const j=this.model.joints[ji];this.jointFrames[ji]={pos:SM.point(T,j.pos),axis:SM.norm(SM.vec(T,j.axis))};let J=SM.ident();const v=this.q[ji]-j.ref;if(j.type==='slide'){let ax=SM.mul(SM.norm(j.axis),v);J[12]=ax[0];J[13]=ax[1];J[14]=ax[2];}else J=SM.rotate(j.axis,v,j.pos);T=SM.matmul(T,J);}this.transforms[i]=T;});
  this.sitePositions=this.model.sites.map(s=>SM.point(this.transforms[s.body],s.pos));return this;
 }
 site(name){const i=this.smap[name];return i===undefined?null:this.sitePositions[i];}
 joint(name){return this.jointFrames[this.jmap[name]];}
 anchorPath(name){const t=this.model.tendons[this.tmap[name]];if(!t)return [];return t.steps.filter(x=>'site'in x).map(x=>this.sitePositions[x.site]);}
 pathLength(name){const t=this.model.tendons[this.tmap[name]];if(!t||t.hasWrap)return null;const p=this.anchorPath(name);let length=0;for(let i=1;i<p.length;i++)length+=Math.hypot(...SM.sub(p[i],p[i-1]));return length;}
 pathDerivative(name,joint){const t=this.model.tendons[this.tmap[name]];if(!t||t.hasWrap||!(joint in this.jmap))return null;const copy={...this.controls},q=copy[joint]??this.model.joints[this.jmap[joint]].ref;const eps=.0008;this.pose({...copy,[joint]:q+eps});const hi=this.pathLength(name),qh=this.q[this.jmap[joint]];this.pose({...copy,[joint]:q-eps});const lo=this.pathLength(name),ql=this.q[this.jmap[joint]];this.pose(copy);return qh>ql? (hi-lo)/(qh-ql):null;}
};

/* Offline WebGL mesh viewer. The renderer does not invent muscle forces. */
'use strict';
window.LabViewer=class {
 constructor(canvas,callbacks={}){this.canvas=canvas;this.callbacks=callbacks;this.gl=canvas.getContext('webgl',{antialias:true,preserveDrawingBuffer:true,alpha:false});if(!this.gl||!this.gl.getExtension('OES_element_index_uint'))throw Error('未能启用三维显示。请在支持WebGL的浏览器中打开，并启用硬件加速。');this.datasets={};this.selected=-1;this.target=[-.17,1.2,-.18];this.yaw=-.72;this.pitch=.08;this.radius=1.25;this.fov=.61;this.paths=[];this.axes=[];this.planes=[];this.pointers=new Map();this.setup();this.bind();this.observer=new ResizeObserver(()=>this.invalidate());this.observer.observe(canvas);}
 setup(){const g=this.gl;const vs=`precision highp float;attribute vec3 aP;attribute vec3 aN;uniform mat4 uVP,uM;varying vec3 N,W;void main(){vec4 p=uM*vec4(aP,1.);W=p.xyz;N=mat3(uM)*aN;gl_Position=uVP*p;gl_PointSize=8.;}`;const fs=`precision highp float;varying vec3 N,W;uniform vec3 uColor,uEye,uPick;uniform float uSelected,uAlpha;uniform bool uPicking,uUnlit;void main(){if(uPicking){gl_FragColor=vec4(uPick,1.);return;}if(uUnlit){gl_FragColor=vec4(uColor,uAlpha);return;}vec3 n=normalize(N);if(!gl_FrontFacing)n=-n;vec3 v=normalize(uEye-W),l=normalize(vec3(-.35,.8,1.1));float d=max(dot(n,l),0.);float f=max(dot(n,normalize(vec3(.8,.25,-.4))),0.);float s=pow(max(dot(n,normalize(l+v)),0.),45.)*.16;float head=max(dot(n,v),0.);vec3 c=uColor*(.57+.26*d+.24*head+.10*f)+s*.65;float rim=pow(1.-abs(dot(v,n)),4.);c+=vec3(.045,.065,.07)*rim;c=mix(c,vec3(.07,.55,.45),uSelected*.53);gl_FragColor=vec4(c,uAlpha);}`;const shader=(t,s)=>{let a=g.createShader(t);g.shaderSource(a,s);g.compileShader(a);if(!g.getShaderParameter(a,g.COMPILE_STATUS))throw Error(g.getShaderInfoLog(a));return a;};this.prog=g.createProgram();g.attachShader(this.prog,shader(g.VERTEX_SHADER,vs));g.attachShader(this.prog,shader(g.FRAGMENT_SHADER,fs));g.linkProgram(this.prog);if(!g.getProgramParameter(this.prog,g.LINK_STATUS))throw Error('三维着色器启动失败');this.at={p:g.getAttribLocation(this.prog,'aP'),n:g.getAttribLocation(this.prog,'aN')};this.u={};['VP','M','Color','Eye','Pick','Selected','Alpha','Picking','Unlit'].forEach(n=>this.u[n]=g.getUniformLocation(this.prog,'u'+n));this.lineBuffer=g.createBuffer();this.normalBuffer=g.createBuffer();g.enable(g.DEPTH_TEST);g.disable(g.CULL_FACE);}
 addDataset(name,model,buffer){const g=this.gl;const upload=(target,a)=>{let b=g.createBuffer();g.bindBuffer(target,b);g.bufferData(target,a,g.STATIC_DRAW);return b;};const meshes=model.parts.map(p=>({p:upload(g.ARRAY_BUFFER,new Float32Array(buffer,p.positions,p.vertexCount*3)),n:upload(g.ARRAY_BUFFER,new Float32Array(buffer,p.normals,p.vertexCount*3)),i:upload(g.ELEMENT_ARRAY_BUFFER,new Uint32Array(buffer,p.indices,p.indexCount)),count:p.indexCount}));this.datasets[name]={model,meshes};}
 use(name,transforms){const d=this.datasets[name];if(!d)throw Error('模型尚未载入');this.dataset=name;this.model=d.model;this.meshes=d.meshes;this.transforms=transforms;this.visible=d.model.parts.map(()=>true);this.selected=-1;this.paths=[];this.axes=[];this.planes=[];this.invalidate();}
 color(p){if(p.color)return p.color.slice(0,3);return p.type==='muscle'?[.70,.29,.25]:p.type==='tendon'?[.88,.83,.64]:p.type==='ligament'?[.83,.70,.44]:[.87,.84,.73];}
 camera(){const M=SM;this.eye=M.add(this.target,[this.radius*Math.cos(this.pitch)*Math.sin(this.yaw),this.radius*Math.sin(this.pitch),this.radius*Math.cos(this.pitch)*Math.cos(this.yaw)]);this.vp=M.matmul(M.perspective(this.fov,this.canvas.width/this.canvas.height,Math.max(.0001,this.radius/2000),20),M.lookAt(this.eye,this.target));}
 resize(){const c=this.canvas,d=Math.min(devicePixelRatio||1,1.8);let w=Math.max(1,Math.round(c.clientWidth*d)),h=Math.max(1,Math.round(c.clientHeight*d));if(c.width!==w||c.height!==h){c.width=w;c.height=h;this.dropFbo();}this.gl.viewport(0,0,w,h);}
 project(p){if(!this.vp)return [0,0,false];let m=this.vp,x=p[0],y=p[1],z=p[2],w=m[3]*x+m[7]*y+m[11]*z+m[15];return [(m[0]*x+m[4]*y+m[8]*z+m[12])/w*this.canvas.clientWidth/2+this.canvas.clientWidth/2,-(m[1]*x+m[5]*y+m[9]*z+m[13])/w*this.canvas.clientHeight/2+this.canvas.clientHeight/2,w>0];}
 fit(){const mn=[Infinity,Infinity,Infinity],mx=[-Infinity,-Infinity,-Infinity];this.model.parts.forEach((p,i)=>{if(!this.visible[i])return;const T=this.transforms[i]||SM.ident();for(let a=0;a<8;a++){let pt=p.bounds[0].map((v,j)=>p.bounds[(a>>j)&1][j]);pt=SM.point(T,pt);pt.forEach((v,j)=>{mn[j]=Math.min(mn[j],v);mx[j]=Math.max(mx[j],v);});}});if(!Number.isFinite(mn[0]))return;this.target=mn.map((v,j)=>(v+mx[j])/2);const size=mx.map((v,j)=>v-mn[j]);const ar=Math.max(.35,this.canvas.clientWidth/Math.max(1,this.canvas.clientHeight));this.radius=Math.max(.12,Math.max(size[1],size[0]/ar,size[2]/ar)/(2*Math.tan(this.fov/2))*1.23+size[2]*.25);this.saveView();this.invalidate();}
 saveView(){this.home={target:this.target.slice(),radius:this.radius,yaw:this.yaw,pitch:this.pitch};}
 reset(){if(this.home){Object.assign(this,this.home);this.target=this.home.target.slice();}this.invalidate();}
 orient(v){this.yaw=v==='front'?0:v==='back'?Math.PI:-Math.PI/2;this.pitch=.06;this.invalidate();}
 invalidate(){if(this.lost||this.pending)return;this.pending=requestAnimationFrame(()=>{this.pending=0;this.draw();});}
 draw(pick=false){if(!this.model||this.lost)return;let g=this.gl;this.resize();this.camera();g.useProgram(this.prog);g.uniformMatrix4fv(this.u.VP,false,this.vp);g.uniform3fv(this.u.Eye,this.eye);g.enableVertexAttribArray(this.at.p);g.enableVertexAttribArray(this.at.n);g.clearColor(pick?0:.944,pick?0:.954,pick?0:.945,1);g.depthMask(true);g.enable(g.DEPTH_TEST);g.disable(g.BLEND);g.clear(g.COLOR_BUFFER_BIT|g.DEPTH_BUFFER_BIT);g.uniform1i(this.u.Picking,pick?1:0);g.uniform1i(this.u.Unlit,0);g.uniform1f(this.u.Alpha,1);pick?g.disable(g.DITHER):g.enable(g.DITHER);
  this.meshes.forEach((m,i)=>{if(!this.visible[i])return;g.bindBuffer(g.ARRAY_BUFFER,m.p);g.vertexAttribPointer(this.at.p,3,g.FLOAT,false,0,0);g.bindBuffer(g.ARRAY_BUFFER,m.n);g.vertexAttribPointer(this.at.n,3,g.FLOAT,false,0,0);g.bindBuffer(g.ELEMENT_ARRAY_BUFFER,m.i);g.uniformMatrix4fv(this.u.M,false,this.transforms[i]||SM.ident());g.uniform3fv(this.u.Color,this.color(this.model.parts[i]));g.uniform1f(this.u.Selected,i===this.selected?1:0);g.uniform3fv(this.u.Pick,[(i+1)%256/255,Math.floor((i+1)/256)/255,0]);g.drawElements(g.TRIANGLES,m.count,g.UNSIGNED_INT,0);});
  if(!pick){g.uniform1i(this.u.Unlit,1);g.uniformMatrix4fv(this.u.M,false,SM.ident());g.uniform1f(this.u.Selected,0);g.disableVertexAttribArray(this.at.n);g.vertexAttrib3f(this.at.n,0,1,0);g.lineWidth(2);for(const plane of (this.planes||[])){const ps=plane.points;if(ps.length!==4)continue;g.enable(g.BLEND);g.blendFunc(g.SRC_ALPHA,g.ONE_MINUS_SRC_ALPHA);g.enable(g.DEPTH_TEST);g.depthMask(false);g.uniform1f(this.u.Alpha,plane.alpha||.12);g.uniform3fv(this.u.Color,plane.color||[.15,.5,.7]);g.bindBuffer(g.ARRAY_BUFFER,this.lineBuffer);g.bufferData(g.ARRAY_BUFFER,new Float32Array([ps[0],ps[1],ps[2],ps[0],ps[2],ps[3]].flat()),g.DYNAMIC_DRAW);g.vertexAttribPointer(this.at.p,3,g.FLOAT,false,0,0);g.drawArrays(g.TRIANGLES,0,6);}g.depthMask(true);g.disable(g.BLEND);g.uniform1f(this.u.Alpha,1);for(const p of [...this.paths,...this.axes]){if(p.points.length<2)continue;g.bindBuffer(g.ARRAY_BUFFER,this.lineBuffer);g.bufferData(g.ARRAY_BUFFER,new Float32Array(p.points.flat()),g.DYNAMIC_DRAW);g.vertexAttribPointer(this.at.p,3,g.FLOAT,false,0,0);g.uniform3fv(this.u.Color,p.color||[.73,.26,.18]);if(p.overlay)g.disable(g.DEPTH_TEST);else g.enable(g.DEPTH_TEST);g.drawArrays(g.LINE_STRIP,0,p.points.length);if(p.dots)g.drawArrays(g.POINTS,0,p.points.length);}g.enable(g.DEPTH_TEST);this.callbacks.draw?.();}
 }
 dropFbo(){const g=this.gl;if(this.fbo){g.deleteFramebuffer(this.fbo);g.deleteTexture(this.tex);g.deleteRenderbuffer(this.depth);this.fbo=null;}}
 pick(x,y){this.resize();const g=this.gl;if(!this.fbo){this.fbo=g.createFramebuffer();this.tex=g.createTexture();this.depth=g.createRenderbuffer();g.bindTexture(g.TEXTURE_2D,this.tex);g.texParameteri(g.TEXTURE_2D,g.TEXTURE_MIN_FILTER,g.NEAREST);g.texParameteri(g.TEXTURE_2D,g.TEXTURE_MAG_FILTER,g.NEAREST);g.texImage2D(g.TEXTURE_2D,0,g.RGBA,this.canvas.width,this.canvas.height,0,g.RGBA,g.UNSIGNED_BYTE,null);g.bindRenderbuffer(g.RENDERBUFFER,this.depth);g.renderbufferStorage(g.RENDERBUFFER,g.DEPTH_COMPONENT16,this.canvas.width,this.canvas.height);g.bindFramebuffer(g.FRAMEBUFFER,this.fbo);g.framebufferTexture2D(g.FRAMEBUFFER,g.COLOR_ATTACHMENT0,g.TEXTURE_2D,this.tex,0);g.framebufferRenderbuffer(g.FRAMEBUFFER,g.DEPTH_ATTACHMENT,g.RENDERBUFFER,this.depth);if(g.checkFramebufferStatus(g.FRAMEBUFFER)!==g.FRAMEBUFFER_COMPLETE)throw Error('结构选择缓冲区不可用');}g.bindFramebuffer(g.FRAMEBUFFER,this.fbo);this.draw(true);const r=this.canvas.getBoundingClientRect(),v=new Uint8Array(4);g.readPixels(SM.clamp(Math.round((x-r.left)*this.canvas.width/r.width),0,this.canvas.width-1),SM.clamp(Math.round((r.bottom-y)*this.canvas.height/r.height),0,this.canvas.height-1),1,1,g.RGBA,g.UNSIGNED_BYTE,v);g.bindFramebuffer(g.FRAMEBUFFER,null);let id=v[0]+v[1]*256-1;this.selected=id>=0&&id<this.model.parts.length?id:-1;this.invalidate();this.callbacks.pick?.(this.selected);}
 bind(){const c=this.canvas;let origin=null,moved=false,pinch=null;c.addEventListener('contextmenu',e=>e.preventDefault());c.addEventListener('pointerdown',e=>{c.setPointerCapture(e.pointerId);this.pointers.set(e.pointerId,[e.clientX,e.clientY]);origin=[e.clientX,e.clientY];moved=false;pinch=null;});c.addEventListener('pointermove',e=>{let old=this.pointers.get(e.pointerId);if(!old)return;let dx=e.clientX-old[0],dy=e.clientY-old[1];this.pointers.set(e.pointerId,[e.clientX,e.clientY]);moved=moved||Math.hypot(e.clientX-origin[0],e.clientY-origin[1])>5;if(this.pointers.size===2){const a=[...this.pointers.values()];let d=Math.hypot(a[0][0]-a[1][0],a[0][1]-a[1][1]);if(pinch)this.radius=SM.clamp(this.radius*pinch/d,.08,9);pinch=d;}else if(e.shiftKey||e.buttons===2){const s=this.radius*.0014;this.target[0]-=dx*Math.cos(this.yaw)*s;this.target[2]+=dx*Math.sin(this.yaw)*s;this.target[1]+=dy*s;}else{this.yaw-=dx*.006;this.pitch=SM.clamp(this.pitch+dy*.006,-1.4,1.4);}this.invalidate();});c.addEventListener('pointerup',e=>{let one=this.pointers.size===1;this.pointers.delete(e.pointerId);if(!moved&&one)this.pick(e.clientX,e.clientY);pinch=null;});c.addEventListener('pointercancel',e=>this.pointers.delete(e.pointerId));c.addEventListener('wheel',e=>{e.preventDefault();this.radius=SM.clamp(this.radius*Math.exp(e.deltaY*.001),.08,9);this.invalidate();},{passive:false});c.addEventListener('webglcontextlost',e=>{e.preventDefault();this.lost=true;this.callbacks.error?.('三维上下文丢失。请导出学习记录后刷新。');});c.addEventListener('keydown',e=>{if(e.key==='ArrowLeft')this.yaw-=.12;else if(e.key==='ArrowRight')this.yaw+=.12;else if(e.key==='ArrowUp')this.pitch-=.1;else if(e.key==='ArrowDown')this.pitch+=.1;else return;e.preventDefault();this.invalidate();});}
};

/* CPU z-buffer compatibility renderer; no mesh decimation or fake anatomy.
 * Used only when a browser cannot create WebGL. Exact source triangles are
 * rasterized at a capped display resolution; playback is correspondingly slower.
 */
'use strict';
window.SoftwareViewer=class {
 constructor(canvas,callbacks={}){this.canvas=canvas;this.callbacks=callbacks;this.ctx=canvas.getContext('2d',{alpha:false});if(!this.ctx)throw Error('浏览器无法创建显示画布。');this.datasets={};this.selected=-1;this.target=[-.17,1.2,-.18];this.yaw=-.72;this.pitch=.08;this.radius=1.25;this.fov=.61;this.paths=[];this.axes=[];this.pointers=new Map();this.software=true;this.bind();this.observer=new ResizeObserver(()=>this.invalidate());this.observer.observe(canvas);}
 addDataset(name,model,buffer){this.datasets[name]={model,meshes:model.parts.map(p=>({p:new Float32Array(buffer,p.positions,p.vertexCount*3),n:new Float32Array(buffer,p.normals,p.vertexCount*3),i:new Uint32Array(buffer,p.indices,p.indexCount),count:p.indexCount}))};}
 resize(){const c=this.canvas;const scale=Math.min(1,940/Math.max(1,c.clientWidth));const w=Math.max(1,Math.round(c.clientWidth*scale)),h=Math.max(1,Math.round(c.clientHeight*scale));if(c.width!==w||c.height!==h){c.width=w;c.height=h;this.image=this.ctx.createImageData(w,h);this.depth=new Float32Array(w*h);this.owners=new Int32Array(w*h);}this.scale=scale;}
 draw(){if(!this.model)return;this.resize();this.camera();const w=this.canvas.width,h=this.canvas.height,pix=this.image.data,dep=this.depth,owner=this.owners;dep.fill(Infinity);owner.fill(-1);for(let i=0;i<pix.length;i+=4){pix[i]=241;pix[i+1]=244;pix[i+2]=240;pix[i+3]=255;}const vp=this.vp;
  for(let part=0;part<this.meshes.length;part++){
   if(!this.visible[part])continue;const mesh=this.meshes[part],T=this.transforms[part]||SM.ident(),P=SM.matmul(vp,T),n=mesh.p.length/3;let sx=new Float32Array(n),sy=new Float32Array(n),sz=new Float32Array(n),lit=new Float32Array(n),valid=new Uint8Array(n);
   for(let i=0;i<n;i++){let a=i*3,x=mesh.p[a],y=mesh.p[a+1],z=mesh.p[a+2],pw=P[3]*x+P[7]*y+P[11]*z+P[15];if(pw<=0)continue;valid[i]=1;sx[i]=(P[0]*x+P[4]*y+P[8]*z+P[12])/pw*w*.5+w*.5;sy[i]=-(P[1]*x+P[5]*y+P[9]*z+P[13])/pw*h*.5+h*.5;sz[i]=(P[2]*x+P[6]*y+P[10]*z+P[14])/pw;
    x=mesh.n[a];y=mesh.n[a+1];z=mesh.n[a+2];let nx=T[0]*x+T[4]*y+T[8]*z,ny=T[1]*x+T[5]*y+T[9]*z,nz=T[2]*x+T[6]*y+T[10]*z;const ln=Math.hypot(nx,ny,nz)||1;nx/=ln;ny/=ln;nz/=ln;lit[i]=.45+.43*Math.max(0,-nx*.25+ny*.55+nz*.80)+.14*Math.max(0,nx*.8+ny*.25-nz*.3);
   }
   const color=this.color(this.model.parts[part]).slice();if(part===this.selected)for(let k=0;k<3;k++)color[k]=color[k]*.55+[.07,.55,.45][k]*.45;
   for(let k=0;k<mesh.i.length;k+=3){const a=mesh.i[k],b=mesh.i[k+1],c=mesh.i[k+2];if(!valid[a]||!valid[b]||!valid[c])continue;let x0=sx[a],y0=sy[a],x1=sx[b],y1=sy[b],x2=sx[c],y2=sy[c];const denom=(y1-y2)*(x0-x2)+(x2-x1)*(y0-y2);if(Math.abs(denom)<.001)continue;const loX=Math.max(0,Math.ceil(Math.min(x0,x1,x2)-.5)),hiX=Math.min(w-1,Math.floor(Math.max(x0,x1,x2)-.5)),loY=Math.max(0,Math.ceil(Math.min(y0,y1,y2)-.5)),hiY=Math.min(h-1,Math.floor(Math.max(y0,y1,y2)-.5));if(loX>hiX||loY>hiY)continue;const inv=1/denom,dax=(y1-y2)*inv,dbx=(y2-y0)*inv;for(let y=loY;y<=hiY;y++){let wa=((y1-y2)*(loX+.5-x2)+(x2-x1)*(y+.5-y2))*inv,wb=((y2-y0)*(loX+.5-x2)+(x0-x2)*(y+.5-y2))*inv;let offset=y*w+loX;for(let x=loX;x<=hiX;x++,offset++,wa+=dax,wb+=dbx){const wc=1-wa-wb;if(wa<-.001||wb<-.001||wc<-.001)continue;const z=wa*sz[a]+wb*sz[b]+wc*sz[c];if(z>=dep[offset])continue;dep[offset]=z;owner[offset]=part;const l=wa*lit[a]+wb*lit[b]+wc*lit[c];let pi=offset*4;pix[pi]=Math.min(255,color[0]*l*255);pix[pi+1]=Math.min(255,color[1]*l*255);pix[pi+2]=Math.min(255,color[2]*l*255);}}}
  }
  this.ctx.putImageData(this.image,0,0);const ctx=this.ctx;ctx.lineWidth=2;ctx.lineCap='round';for(const p of (this.planes||[])){const ps=p.points.map(pt=>this.project(pt));if(ps.some(pt=>!pt[2]))continue;ctx.save();ctx.globalAlpha=p.alpha||.12;ctx.fillStyle='rgb(40,125,180)';ctx.beginPath();ps.forEach((pt,i)=>i?ctx.lineTo(pt[0]*this.scale,pt[1]*this.scale):ctx.moveTo(pt[0]*this.scale,pt[1]*this.scale));ctx.closePath();ctx.fill();ctx.restore();}for(const p of [...this.paths,...this.axes]){if(p.points.length<2)continue;const color=p.color||[.7,.3,.2];ctx.strokeStyle=ctx.fillStyle=`rgb(${color.map(x=>Math.round(x*255)).join(',')})`;for(let i=1;i<p.points.length;i++){const p0=p.points[i-1],p1=p.points[i];const a=this.project(p0),b=this.project(p1);if(!a[2]||!b[2])continue;ctx.beginPath();let pen=false;const steps=Math.max(2,Math.ceil(Math.hypot(a[0]-b[0],a[1]-b[1])*this.scale));for(let j=0;j<=steps;j++){let t=j/steps;const pt=p0.map((v,k)=>v+(p1[k]-v)*t),x=(a[0]+(b[0]-a[0])*t)*this.scale,y=(a[1]+(b[1]-a[1])*t)*this.scale;const zz=(vp[2]*pt[0]+vp[6]*pt[1]+vp[10]*pt[2]+vp[14])/(vp[3]*pt[0]+vp[7]*pt[1]+vp[11]*pt[2]+vp[15]),px=Math.round(x),py=Math.round(y);const inside=px>=0&&px<w&&py>=0&&py<h;if(inside&&(p.overlay||zz<=dep[py*w+px]+.00008)){if(pen)ctx.lineTo(x,y);else ctx.moveTo(x,y);pen=true;}else pen=false;}ctx.stroke();}if(p.dots)for(const pt of p.points){const a=this.project(pt);if(!a[2])continue;const x=a[0]*this.scale,y=a[1]*this.scale;const zz=(vp[2]*pt[0]+vp[6]*pt[1]+vp[10]*pt[2]+vp[14])/(vp[3]*pt[0]+vp[7]*pt[1]+vp[11]*pt[2]+vp[15]),off=Math.round(y)*w+Math.round(x);if(p.overlay||zz<=(dep[off]??Infinity)+.00008){ctx.beginPath();ctx.arc(x,y,3.2,0,Math.PI*2);ctx.fill();}}}
  this.callbacks.draw?.();
 }
 pick(x,y){this.draw();const r=this.canvas.getBoundingClientRect(),px=Math.floor((x-r.left)*this.canvas.width/r.width),py=Math.floor((y-r.top)*this.canvas.height/r.height);this.selected=px>=0&&px<this.canvas.width&&py>=0&&py<this.canvas.height?this.owners[py*this.canvas.width+px]:-1;this.invalidate();this.callbacks.pick?.(this.selected);}
};
for(const name of ['use','color','camera','project','fit','saveView','reset','orient','invalidate','bind'])SoftwareViewer.prototype[name]=LabViewer.prototype[name];
window.createLabViewer=(canvas,callbacks)=>{try{return new LabViewer(canvas,callbacks);}catch(e){console.warn('WebGL不可用，改用同网格软件光栅显示。',e.message);return new SoftwareViewer(canvas,callbacks);}};

/* V6 double-role teaching scene. BodyParts3D exposed surface + original garments.
   Original approximate 16-bone skinning: qualitative role movement, NOT muscle
   simulation or a validated clinical measure. Anatomical analysis remains in V5.
   All role commands are controlled by ScenarioEngine; this renderer awards no scores. */
'use strict';
window.ClinicGeometry={
 mesh(name,v,i,color){const p=new Float32Array(v.flat()),ix=new Uint32Array(i.flat()),n=new Float32Array(p.length);for(let k=0;k<ix.length;k+=3){const a=ix[k]*3,b=ix[k+1]*3,c=ix[k+2]*3,u=[p[b]-p[a],p[b+1]-p[a+1],p[b+2]-p[a+2]],w=[p[c]-p[a],p[c+1]-p[a+1],p[c+2]-p[a+2]],q=SM.cross(u,w);for(const x of [a,b,c])for(let j=0;j<3;j++)n[x+j]+=q[j];}for(let k=0;k<n.length;k+=3){const q=Math.hypot(n[k],n[k+1],n[k+2])||1;for(let j=0;j<3;j++)n[k+j]/=q;}return {name,p,n,i:ix,color};},
 box(name,center,size,color){let v=[],i=[];for(const [normal,axes] of [[[1,0,0],[1,2]], [[-1,0,0],[2,1]],[[0,1,0],[2,0]],[[0,-1,0],[0,2]],[[0,0,1],[0,1]],[[0,0,-1],[1,0]]]){let base=v.length,axis=normal.findIndex(x=>x!==0);for(const [u,w] of [[-1,-1],[1,-1],[1,1],[-1,1]]){let p=center.slice();p[axis]+=normal[axis]*size[axis]/2;p[axes[0]]+=u*size[axes[0]]/2;p[axes[1]]+=w*size[axes[1]]/2;v.push(p);}i.push([base,base+1,base+2],[base,base+2,base+3]);}return this.mesh(name,v,i,color);},
 sphere(name,c,s,color,st=12,sl=20){const v=[],i=[];for(let j=0;j<=st;j++){const a=Math.PI*j/st;for(let k=0;k<=sl;k++){const b=2*Math.PI*k/sl;v.push([c[0]+s[0]*Math.sin(a)*Math.cos(b),c[1]+s[1]*Math.cos(a),c[2]+s[2]*Math.sin(a)*Math.sin(b)]);}}for(let j=0;j<st;j++)for(let k=0;k<sl;k++){const a=j*(sl+1)+k;i.push([a,a+1,a+sl+1],[a+1,a+sl+2,a+sl+1]);}return this.mesh(name,v,i,color);},
 rod(name,a,b,r,color,steps=14){let v=[],i=[],d=SM.sub(b,a),u=SM.norm(SM.cross(SM.norm(d),Math.abs(d[1])<.9*Math.hypot(...d)?[0,1,0]:[1,0,0])),w=SM.norm(SM.cross(d,u));for(const c of [a,b])for(let j=0;j<steps;j++){const x=2*Math.PI*j/steps;v.push(c.map((z,k)=>z+r*(Math.cos(x)*u[k]+Math.sin(x)*w[k])));}for(let j=0;j<steps;j++){const k=(j+1)%steps;i.push([j,k,j+steps],[k,k+steps,j+steps]);}v.push(a,b);for(let j=0;j<steps;j++){i.push([2*steps,j,(j+1)%steps]);i.push([2*steps+1,(j+1)%steps+steps,j+steps]);}return this.mesh(name,v,i,color);}
};
window.ClinicalScene=class {
 constructor(canvas,model,raw,cb={}){this.canvas=canvas;this.cb=cb;this.actorModel=model;this.raw=raw;this.bones=model.bones;this.bmap=Object.fromEntries(model.bones.map((b,i)=>[b.name,i]));this.roles={};this.parts=[];this.meshData=[];this.transforms=[];this.animation=null;this.case=null;this.support=false;this.seated=false;this.active=false;this.poseVersion=0;this.lastTick=0;this.materials={therapist:{skin:[.77,.60,.47],hair:[.095,.09,.077],lip:[.58,.37,.33],shirt:[.11,.44,.42],trim:[.08,.32,.31],pants:[.12,.24,.27],shoe:[.20,.27,.29],sole:[.79,.81,.76]},patient:{skin:[.76,.58,.44],hair:[.19,.17,.15],lip:[.58,.36,.31],shirt:[.73,.49,.27],trim:[.49,.31,.20],pants:[.22,.29,.35],shoe:[.26,.28,.31],sole:[.83,.83,.79]}};
 this.viewer=createLabViewer(canvas,{draw:()=>this.drawLabels(),pick:i=>{if(i>=0){const role=this.parts[i].role;this.cb.pick?.(role||'room',this.parts[i].name);}},error:m=>cb.error?.(m)});this.makeRoom();this.makeActor('therapist');this.makeActor('patient');this.finalize();this.setPose('therapist','rest',0);this.setPose('patient','rest',0);this.frame('both');this.onResize=new ResizeObserver(()=>{if(this.active)this.frame(this.viewMode||'both',true);});this.onResize.observe(canvas);this.tick=this.tick.bind(this);requestAnimationFrame(this.tick);}
 push(mesh,extra={}){const index=this.parts.length;const p={name:mesh.name,label:mesh.name,color:mesh.color,...extra,vertexCount:mesh.p.length/3,indexCount:mesh.i.length,bounds:[[0,0,0],[0,0,0]]};for(let j=0;j<3;j++){let mn=Infinity,mx=-Infinity;for(let k=j;k<mesh.p.length;k+=3){mn=Math.min(mn,mesh.p[k]);mx=Math.max(mx,mesh.p[k]);}p.bounds[0][j]=mn;p.bounds[1][j]=mx;}this.parts.push(p);this.meshData.push(mesh);this.transforms.push(SM.ident());return index;}
 makeRoom(){const G=ClinicGeometry;this.roomIDs={chair:[],rails:[],plate:[],table:[],props:[]};const add=(m,group)=>{let i=this.push(m,{group});if(group)this.roomIDs[group].push(i);return i;};
 add(G.box('治疗室地面',[0,-.055,0],[5.4,.10,4.8],[.86,.87,.82]));
 for(let x=-2.5;x<2.6;x+=.45)add(G.box('浅木地板',[x,-.002,0],[.435,.012,4.8],[.90+(Math.round(x*10)%2)*.006,.92,.90]));
 add(G.box('后墙',[0,1.20,-1.35],[5.4,2.45,.10],[.965,.975,.965]));add(G.box('墙面下护板',[0,.38,-1.28],[5.4,.73,.045],[.83,.89,.86]));
 add(G.box('通光窗框',[-1.45,1.60,-1.26],[1.4,1.04,.08],[.97,.97,.94]));add(G.box('日光窗',[-1.45,1.6,-1.20],[1.26,.92,.024],[.76,.86,.86]));
 add(G.box('窗中框',[-1.45,1.60,-1.17],[.035,.96,.04],[.96,.97,.93]));add(G.box('窗横框',[-1.45,1.64,-1.17],[1.30,.035,.04],[.96,.97,.93]));
 add(G.box('教学墙卡',[.55,1.64,-1.20],[1.13,.53,.033],[.96,.97,.94]));for(let j=0;j<3;j++)add(G.box('教学墙卡线',[.45,1.76-j*.11,-1.17],[.72-j*.12,.024,.014],[.40,.59,.57]));
 // Background examination couch: no unvalidated contact maneuvers are shown.
 add(G.box('评估床垫',[1.63,.63,-.65],[.72,.12,1.20],[.17,.43,.42]));add(G.box('床架',[1.63,.52,-.65],[.62,.12,1.07],[.81,.86,.83]));for(const x of [1.35,1.91])for(const z of [-1.10,-.20])add(G.rod('床脚',[x,.02,z],[x,.52,z],.031,[.61,.68,.66]));
 // Patient's stable chair; excluded in standing cases.
 const cx=.58,cz=.15;add(G.box('稳定椅坐垫',[cx,.455,cz],[.46,.075,.43],[.26,.43,.43]),'chair');add(G.box('稳定椅靠背',[cx,.81,cz-.208],[.46,.59,.055],[.26,.43,.43]),'chair');for(const x of [cx-.18,cx+.18])for(const z of [cz-.16,cz+.16])add(G.rod('椅脚',[x,.01,z],[x,.43,z],.022,[.51,.59,.57]),'chair');
 // Therapist sits opposite in seated observation; stool is not a treatment contact.
 add(G.box('治疗师凳面',[-.63,.45,.03],[.42,.065,.40],[.22,.40,.39]),'chair');for(const x of [-.79,-.47])for(const z of [-.12,.18])add(G.rod('治疗师凳脚',[x,.01,z],[x,.43,z],.02,[.50,.61,.58]),'chair');
 // A rail is available beside the participant, not through them.
 for(const z of [-.13,.68])add(G.rod('支撑立柱',[1.03,.02,z],[1.03,.94,z],.022,[.52,.64,.62]),'rails');add(G.rod('支撑扶手',[1.03,.94,-.15],[1.03,.94,.71],.028,[.24,.41,.40]),'rails');
 add(G.box('教学测力板',[.58,.013,.20],[.83,.024,.72],[.28,.47,.47]),'plate');add(G.box('板面',[.58,.031,.20],[.76,.012,.65],[.71,.81,.77]),'plate');
 // Small desk stays behind therapist, avoids patient movement envelope.
 add(G.box('教学台面',[-1.56,.79,.10],[.6,.075,.82],[.77,.64,.47]),'table');for(const x of [-1.8,-1.3])for(const z of [-.22,.42])add(G.rod('台脚',[x,.02,z],[x,.75,z],.024,[.40,.52,.49]),'table');
 add(G.box('平板显示器',[-1.57,.98,-.08],[.31,.25,.045],[.18,.30,.31]));add(G.box('平板屏幕',[-1.57,.99,-.05],[.27,.19,.01],[.70,.82,.78]));
 add(G.rod('花盆',[-2.13,0,-.72],[-2.13,.24,-.72],.15,[.71,.61,.47]));for(let j=0;j<7;j++){let a=j*2.4;add(G.rod('绿植枝',[-2.13,.23,-.72],[-2.13+Math.sin(a)*.17,.52+j*.04,-.72+Math.cos(a)*.12],.008,[.25,.41,.30]));add(G.sphere('绿植叶',[-2.13+Math.sin(a)*.20,.55+j*.04,-.72+Math.cos(a)*.14],[.105,.032,.07],[.23,.44,.32],6,10));}
 // Contact shadows are illustrations, not simulated forces.
 for(const [x,z] of [[-.63,.1],[.58,.3]])add(G.sphere('接触阴影',[x,.005,z],[.35,.005,.30],[.66,.69,.63],3,28));
 this.cupID=add(G.rod('无质量教学杯',[0,0,0],[0,.09,0],.033,[.93,.90,.77]),'props');this.badgeID=add(G.box('治疗师胸牌',[-.085,1.25,.143],[.060,.076,.012],[.96,.97,.93]),'props');
 }
 makeActor(role){let indices=[];for(const p of this.actorModel.parts){const P=new Float32Array(this.raw,p.positions,p.vertexCount*3),N=new Float32Array(this.raw,p.normals,p.vertexCount*3),I=new Uint32Array(this.raw,p.indices,p.indexCount),W=new Float32Array(this.raw,p.weights,p.vertexCount*4),J=new Uint16Array(this.raw,p.joints,p.vertexCount*4);let i=this.push({name:p.label||p.name,p:P.slice(),n:N.slice(),i:I.slice(),color:this.materials[role][p.material]}, {role,material:p.material});this.meshData[i].restP=P;this.meshData[i].restN=N;this.meshData[i].weights=W;this.meshData[i].joints=J;indices.push(i);}this.roles[role]={indices,pose:'rest',p:0,palette:[],world:[],base:[role==='patient'?.58:-.63,0,role==='patient'?.20:.03],yaw:role==='patient'?.03:.32};}
 finalize(){let bytes=0;for(let i=0;i<this.parts.length;i++){const p=this.parts[i],m=this.meshData[i];p.positions=bytes;bytes+=m.p.byteLength;p.normals=bytes;bytes+=m.n.byteLength;p.indices=bytes;bytes+=m.i.byteLength;}const raw=new ArrayBuffer(bytes);for(let i=0;i<this.parts.length;i++){const p=this.parts[i],m=this.meshData[i];new Float32Array(raw,p.positions,p.vertexCount*3).set(m.p);new Float32Array(raw,p.normals,p.vertexCount*3).set(m.n);new Uint32Array(raw,p.indices,p.indexCount).set(m.i);m.p=new Float32Array(raw,p.positions,p.vertexCount*3);m.n=new Float32Array(raw,p.normals,p.vertexCount*3);}this.viewer.addDataset('clinic',{parts:this.parts},raw);this.viewer.use('clinic',this.transforms);this.buffer=raw;}
 setCase(c){this.case=c;this.stop();this.support=false;this.seated=!!c.seated;this.roles.patient.base=[.58,this.seated?-.415:0,.20];this.roles.therapist.base=[-.63,this.seated?-.415:0,.03];this.roles.therapist.yaw=.32;this.setPose('patient',c.motion,0);this.setPose('therapist','rest',0);for(let i=0;i<this.parts.length;i++){let g=this.parts[i].group;this.viewer.visible[i]=g==='chair'?this.seated:g==='rails'?['balance','hip','gait'].includes(c.motion):g==='plate'?['balance','foot'].includes(c.motion):g==='props'?false:true;}this.viewer.visible[this.badgeID]=true;this.updateBadge();this.frame('both');}
 static translation(p){const m=SM.ident();m[12]=p[0];m[13]=p[1];m[14]=p[2];return m;}
 poseAngles(kind,t,role){let q={},r=[0,0,0];const isPatient=role==='patient',seated=this.seated; if(seated){q.thigh_r=[-Math.PI/2,0,0];q.thigh_l=[-Math.PI/2,0,0];q.shin_r=[Math.PI/2,0,0];q.shin_l=[Math.PI/2,0,0];}
 const wave=Math.sin(Math.PI*t);switch(kind){
 case 'greet':q.upper_r=[-.22,0,-.20];q.fore_r=[-.55-.28*wave,0,0];q.head=[.035*wave,-.08,0];break;
 case 'shoulder':q.upper_r=[0,0,-1.13*t];break;
 case 'shoulderFlex':q.upper_r=[-1.13*t,0,0];break; // qualitative, not source scapulohumeral rhythm
 case 'elbow':case 'contraction':q.fore_r=[-1.70*t,0,0];break;
 case 'wrist':q.fore_r=[-1.50,0,0];q.hand_r=[-.52*t,0,0];break;
 case 'hip':q.thigh_r=[-.82*t,0,0];q.shin_r=[.92*t,0,0];break;
 case 'knee':if(!seated){q.thigh_r=[-.72,0,0];q.shin_r=[1.20-.72*t,0,0];}else q.shin_r=[Math.PI/2-1.24*t,0,0];break;
 case 'ankle':q.foot_r=[-.30*t,0,0];break;
 case 'spine':q.chest=[.25*t,0,0];q.neck=[-.04*t,0,0];break;
 case 'breathing':q.chest=[-.012*wave,0,0];break;
 case 'cardio':q.fore_r=[-.9,0,0];q.hand_r=[0,0,0];break;
 case 'foot':q.chest=[0,0,0];break;
 case 'balance':q.chest=[0,0,.013*Math.sin(t*8)];q.neck=[0,0,-.013*Math.sin(t*8)];break;
 case 'gait': {const cycle=t*2*Math.PI;let a=.27*Math.sin(cycle),b=.27*Math.sin(cycle+Math.PI);q.thigh_r=[-a,0,0];q.thigh_l=[-b,0,0];q.shin_r=[Math.max(0,-Math.sin(cycle))*.50,0,0];q.shin_l=[Math.max(0,Math.sin(cycle))*.50,0,0];q.upper_r=[a*.5,0,0];q.upper_l=[b*.5,0,0];r[2]=t*.6;break;}
 }
 if(isPatient&&this.support&&['hip','balance','gait'].includes(kind)){q.upper_l=[-.15,0,.28];q.fore_l=[-.62,0,0];}
 return {q,root:r};}
 setPose(role,kind,t){const a=this.roles[role];if(!a)return;t=SM.clamp(t,0,1);a.pose=kind;a.p=t;const {q,root}=this.poseAngles(kind,t,role);a.rootOffset=root;const world=[],pal=[];for(let j=0;j<this.bones.length;j++){const b=this.bones[j],parent=b.parent>=0?this.bones[b.parent].pos:[0,0,0];let local=ClinicalScene.translation(SM.sub(b.pos,parent)),angles=q[b.name]||[0,0,0];for(let k=0;k<3;k++){const axis=[0,0,0];axis[k]=1;if(angles[k])local=SM.matmul(local,SM.rotate(axis,angles[k]));}world[j]=b.parent<0?local:SM.matmul(world[b.parent],local);pal[j]=SM.matmul(world[j],ClinicalScene.translation(b.pos.map(x=>-x)));}a.world=world;a.palette=pal;
 const tr=SM.matmul(ClinicalScene.translation(SM.add(a.base,root)),SM.rotate([0,1,0],a.yaw));a.transform=tr;
 for(const idx of a.indices){const m=this.meshData[idx],P=m.restP,N=m.restN,W=m.weights,J=m.joints;for(let vi=0;vi<P.length/3;vi++){let x=0,y=0,z=0,nx=0,ny=0,nz=0,p=vi*3,w=vi*4;for(let k=0;k<4;k++){const weight=W[w+k];if(weight<=0)continue;const T=pal[J[w+k]],xx=P[p],yy=P[p+1],zz=P[p+2];x+=weight*(T[0]*xx+T[4]*yy+T[8]*zz+T[12]);y+=weight*(T[1]*xx+T[5]*yy+T[9]*zz+T[13]);z+=weight*(T[2]*xx+T[6]*yy+T[10]*zz+T[14]);const ux=N[p],uy=N[p+1],uz=N[p+2];nx+=weight*(T[0]*ux+T[4]*uy+T[8]*uz);ny+=weight*(T[1]*ux+T[5]*uy+T[9]*uz);nz+=weight*(T[2]*ux+T[6]*uy+T[10]*uz);}m.p[p]=x;m.p[p+1]=y;m.p[p+2]=z;m.n[p]=nx;m.n[p+1]=ny;m.n[p+2]=nz;}this.transforms[idx]=tr;if(this.viewer.gl){const gl=this.viewer.gl,buf=this.viewer.meshes[idx];gl.bindBuffer(gl.ARRAY_BUFFER,buf.p);gl.bufferSubData(gl.ARRAY_BUFFER,0,m.p);gl.bindBuffer(gl.ARRAY_BUFFER,buf.n);gl.bufferSubData(gl.ARRAY_BUFFER,0,m.n);}}
 this.updateBadge();this.updateCup();this.poseVersion++;this.viewer.invalidate();}
 updateBadge(){const a=this.roles.therapist;if(a?.palette?.[1])this.transforms[this.badgeID]=SM.matmul(a.transform,a.palette[1]);}
 updateCup(){const a=this.roles.patient;if(!a?.world[6])return;this.viewer.visible[this.cupID]=this.case?.motion==='contraction';const pt=this.bonePoint('patient',6,[-.285,.80,.052]);const t=ClinicalScene.translation(SM.add(pt,[0,-.05,.025]));this.transforms[this.cupID]=t;}
 bonePoint(role,i,rest){const a=this.roles[role];return SM.point(a.transform,SM.point(a.palette[i],rest||this.bones[i].pos));}
 focusPoint(){const m=this.case?.motion;let i=5,p=[-.237,1.1,0];if(m==='shoulder'){i=4;p=[-.20,1.34,0];}else if(m==='wrist'){i=6;p=[-.285,.83,.02];}else if(m==='hip'){i=10;p=[-.08,.84,0];}else if(m==='knee'){i=11;p=[-.08,.48,.02];}else if(m==='ankle'||m==='foot'){i=12;p=[-.08,.08,.07];}else if(m==='spine'||m==='breathing'||m==='cardio'){i=1;p=[0,1.20,.10];}else if(m==='gait'||m==='balance'){i=0;p=[0,.95,0];}return this.bonePoint('patient',i,p);}
 frame(mode='both',keepAngle=false){this.viewMode=mode;const v=this.viewer,w=Math.max(1,this.canvas.clientWidth),h=Math.max(1,this.canvas.clientHeight),ar=w/h;let target,size;if(mode==='detail'){target=this.focusPoint();size=['wrist','ankle','foot'].includes(this.case?.motion)?[.58,.68,.6]:[1.10,1.15,.85];}else if(mode==='patient'){target=[.58,this.seated?.79:.94,.32];size=[1.05,this.seated?1.45:1.88,1.2];}else{target=[-.04,this.seated?.73:.96,.20];size=[2.35,this.seated?1.65:2.07,1.8];}v.target=target;if(!keepAngle){v.yaw=mode==='detail'?.55:.12;v.pitch=.065;}v.radius=Math.max(size[1],size[0]/ar)/(2*Math.tan(v.fov/2))*1.09+size[2]*.30;v.saveView();v.invalidate();}
 orient(side){this.viewer.yaw=side==='front'?0:side==='side'?-.92:Math.PI;this.viewer.pitch=.065;this.viewer.invalidate();}
 positionTherapist(which){this.stop();this.roles.therapist.base=this.seated?[-.63,-.415,.03]:(which==='side'?[-.65,0,-.48]:[-.63,0,.03]);this.roles.therapist.yaw=which==='side'?.95:.32;this.setPose('therapist','rest',0);}
 setSupport(value){this.support=value;for(const i of this.roomIDs.rails)this.viewer.visible[i]=value||['balance','hip','gait'].includes(this.case?.motion);this.setPose('patient',this.case.motion,this.roles.patient.p);}
 play(role,kind,{from=0,to=1,duration=4000,cycle=false,speed=1,onEnd=()=>{},onTick=()=>{}}={}){this.stop();this.animation={role,kind,from,to,duration,cycle,speed,elapsed:0,onEnd,onTick};this.setPose(role,kind,from);this.cb.play?.(role,true);}
 stop(){const a=this.animation;this.animation=null;if(a)this.cb.play?.(a.role,false);}
 tick(time){const dt=Math.min(80,Math.max(0,time-(this.lastTick||time)));this.lastTick=time;const a=this.animation;if(a&&this.active&&!document.hidden){a.elapsed+=dt*a.speed;const f=SM.clamp(a.elapsed/a.duration,0,1),e=f*f*(3-2*f),phase=a.cycle?(f<.5?2*f:2*(1-f)):e;let t=a.from+(a.to-a.from)*phase;if(a.onTick(t,f)!==false&&this.animation===a)this.setPose(a.role,a.kind,t);if(f>=1&&this.animation===a){this.animation=null;this.cb.play?.(a.role,false);a.onEnd();}}if(this.active)requestAnimationFrame(this.tick);}
 drawLabels(){this.cb.frame?.({therapist:this.viewer.project(this.bonePoint('therapist',3,[0,1.80,0])),patient:this.viewer.project(this.bonePoint('patient',3,[0,1.80,0])),target:this.viewer.project(this.focusPoint()),renderer:this.viewer.software?'software':'webgl'});}
 photo(){try{const c=document.createElement('canvas');c.width=360;c.height=Math.round(this.canvas.height/this.canvas.width*360);const ctx=c.getContext('2d');ctx.drawImage(this.canvas,0,0,c.width,c.height);return c.toDataURL('image/jpeg',.65);}catch{return null;}}
};

/* Real user-supplied mesh stage. Prescribed source kinematics, not muscle simulation. */
(function(g){'use strict';
const M=()=>g.SM;
const corners=(p,T)=>Array.from({length:8},(_,i)=>SM.point(T,p.bounds[0].map((x,j)=>p.bounds[(i>>j)&1][j])));
class ModelStage{
 constructor(canvas,assets,onPick){this.canvas=canvas;this.assets=assets;this.viewer=createLabViewer(canvas,{pick:i=>{this.viewer.selected=i;this.viewer.invalidate();onPick?.(i>=0?this.viewer.model.parts[i]:null);},error:e=>console.warn(e)});this.kin=new NativeKinematics(assets.native.model);this.native=assets.native.model;this.anatomy=assets.anatomy.model;this.motions=Object.fromEntries(assets.motions.map(m=>[m.id,m]));this.viewer.addDataset('native',this.native,assets.native.raw);this.viewer.addDataset('anatomy',this.anatomy,assets.anatomy.raw);this.action=null;this.mode='native';this.layer='both';this.showGuides=true;this.t=.0;this.basePose=this.kin.pose({}).transforms.map(x=>x.slice());this.fullMask=[];}
 supportsAnatomy(a){const r=a.region||this.motions[a.motion]?.region;return ['upper','shoulder','hand','knee','lower','foot'].includes(r);}
 nativeMask(a){let r=a.region||this.motions[a.motion]?.region;return this.native.parts.map(p=>{let b=p.body;if(a.kind==='fsu')return [4,5].includes(b);if(r==='head')return b===10||b===11; if(r==='trunk')return b>=1&&b<=11||b===88;if(r==='shoulder')return b>=16&&b<=51||b===8;if(r==='upper')return b>=22&&b<=51;if(r==='hand')return b>=(a.motion.startsWith('finger')||a.motion==='thumb'?25:23)&&b<=51;if(r==='lower')return [1,88,89,90,91,92,93,94].includes(b);if(r==='knee')return [89,90,94].includes(b);if(r==='foot')return [90,91,92,93].includes(b);return false;});}
 anatomyMask(a){const r=a.region||this.motions[a.motion]?.region;return this.anatomy.parts.map(p=>{const s=p.regions||[],n=(p.label||p.name||'').toLowerCase();let ok=r==='shoulder'?(s.includes('shoulder')||(p.type==='bone'&&(s.includes('elbow-forearm')||s.includes('wrist-hand')))):r==='upper'?(s.includes('elbow-forearm')||(p.type==='bone'&&s.includes('wrist-hand'))):r==='hand'?(s.includes('wrist-hand')||(p.type==='bone'&&s.includes('elbow-forearm')&&(p.english||'').toLowerCase().includes('ulna'))):r==='lower'?s.includes('hip-pelvis')||s.includes('knee'):r==='knee'?s.includes('knee'):r==='foot'?s.includes('ankle-foot'):false;return ok&&(this.layer==='both'||(this.layer==='muscle'?p.type==='muscle'||p.type==='tendon':p.type==='bone'));});}
 controls(a,t){const m=this.motions[a.motion];if(!m)return{};const q={...(m.base||{})},r=a.range||[m.min,m.max];q[m.joint]=(r[0]+(r[1]-r[0])*t)*Math.PI/180;
 if(a.motion==='fingerFlex'){for(let d=2;d<=5;d++){q[`mcp${d}_flexion_r`]=1.2*t;q[`pm${d}_flexion_r`]=1.2*t;q[`md${d}_flexion_r`]=.8*t;}}
 if(a.motion==='fingerSpread'){for(let d=2;d<=5;d++)q[`mcp${d}_abduction_r`]=d===3?0:(d===2?.30:-.30)*t;}
 if(a.motion==='thumb'){q.cmc_abduction_r=(t-.5)*.4;q.mp_flexion_r=.55*t;q.ip_flexion_r=.35*t;}
 return q;}
 fsuTransforms(a,t){const lower=this.native.parts.find(p=>p.body===4),upper=this.native.parts.find(p=>p.body===5);const pts=[...corners(lower,this.basePose[4]),...corners(upper,this.basePose[5])];const min=[0,1,2].map(j=>Math.min(...pts.map(p=>p[j]))),max=[0,1,2].map(j=>Math.max(...pts.map(p=>p[j])));const c=min.map((v,j)=>(v+max[j])/2);let T=SM.ident();const u=2*t-1;
 if(a.motion==='flexion')T=SM.rotate([1,0,0],u*.087266,c);
 if(a.motion==='lateral')T=SM.rotate([0,0,1],u*.087266,c);
 if(a.motion==='axial')T=SM.rotate([0,1,0],u*.087266,c);
 if(a.motion==='ap')T[14]=u*.002;
 if(a.motion==='ml')T[12]=u*.002;
 if(a.motion==='axialT')T[13]=u*.002;
 this.fsuCenter=c;this.fsuTransform=T;
 return this.native.parts.map(p=>p.body===5?SM.matmul(T,this.basePose[p.body]):this.basePose[p.body]);}
 set(a,t=0,mode='native'){this.action=a;this.t=t;this.mode=mode;this.viewer.use(mode==='anatomy'?'anatomy':'native',[]);this.viewer.visible=mode==='anatomy'?this.anatomyMask(a):this.nativeMask(a);this.fullMask=this.viewer.visible.slice();this.viewer.yaw= mode==='anatomy'?-.45: a.kind==='fsu'?-.80:(a.region||this.motions[a.motion]?.region)==='trunk'?-.85: -.62;this.viewer.pitch=.1;this.update(t);this.fitAll();}
 update(t){this.t=t;if(!this.action)return;const a=this.action;if(this.mode==='anatomy'){this.viewer.transforms=this.anatomy.parts.map(()=>SM.ident());this.viewer.paths=[];this.viewer.axes=[];this.viewer.planes=[];}
 else{if(a.kind==='fsu'){this.viewer.transforms=this.fsuTransforms(a,t);this.viewer.planes=[];const c=this.fsuCenter,l=.027;this.viewer.axes=[{points:[SM.add(c,[-l,0,0]),SM.add(c,[l,0,0])],color:[.74,.30,.25],overlay:true},{points:[SM.add(c,[0,-l,0]),SM.add(c,[0,l,0])],color:[.15,.48,.40],overlay:true},{points:[SM.add(c,[0,0,-l]),SM.add(c,[0,0,l])],color:[.26,.43,.65],overlay:true}];}
 else{this.kin.pose(this.controls(a,t));this.viewer.transforms=this.native.parts.map(p=>this.kin.transforms[p.body]);this.viewer.axes=[];this.viewer.planes=[];if(this.showGuides)this.motionGuides(a,t);}
 this.viewer.paths=[];}
 this.viewer.invalidate();}

 motionGuides(a,t){
 const m=this.motions[a.motion];if(!m)return;const frame=this.kin.joint(m.joint);if(!frame)return;
 const n=SM.norm(frame.axis),c=frame.pos.slice();let seed=Math.abs(n[1])<.85?[0,1,0]:[0,0,1];
 const u=SM.norm(SM.cross(n,seed)),v=SM.norm(SM.cross(n,u));
 const r=['hand','foot'].includes(a.region||m.region)?.075:(a.region||m.region)==='trunk'?.18:.14;
 const at=(x,y)=>SM.add(c,SM.add(SM.mul(u,x),SM.mul(v,y)));
 this.guide={joint:m.joint,center:c,normal:n,plane:[at(-r,-r),at(r,-r),at(r,r),at(-r,r)]};
 this.viewer.planes=[{points:this.guide.plane,color:[.12,.51,.72],alpha:.11}];
 const axis=[SM.add(c,SM.mul(n,-r*.8)),SM.add(c,SM.mul(n,r*.8))];
 this.viewer.axes=[{points:axis,color:[.08,.43,.66],overlay:true},{points:[...this.guide.plane,this.guide.plane[0]],color:[.44,.65,.76],overlay:false}];
 // Arc follows increasing action-slider coordinate. It is a direction guide, not muscle force.
 const range=a.range||[m.min,m.max],sign=range[1]>=range[0]?1:-1;
 const theta=sign*(.45+t*1.8),rr=r*.63,arc=Array.from({length:24},(_,i)=>at(Math.cos(theta-sign*.7+sign*.7*i/23)*rr,Math.sin(theta-sign*.7+sign*.7*i/23)*rr));
 const tip=arc[arc.length-1],prev=arc[arc.length-3],dir=SM.norm(SM.sub(tip,prev)),side=SM.norm(SM.cross(n,dir));
 const back=SM.sub(tip,SM.mul(dir,r*.12));
 this.viewer.axes.push({points:arc,color:[.72,.29,.16],overlay:true},{points:[SM.add(back,SM.mul(side,r*.055)),tip,SM.sub(back,SM.mul(side,r*.055))],color:[.72,.29,.16],overlay:true});
 }
 fitAll(){const v=this.viewer,a=this.action;if(!a)return;let points=[];const t=this.t;for(const u of this.mode==='anatomy'?[t]:[0,.2,.4,.6,.8,1]){this.update(u);this.viewer.model.parts.forEach((p,i)=>{if(v.visible[i])points.push(...corners(p,v.transforms[i]||SM.ident()));});}this.update(t);if(!points.length)return;const min=[0,1,2].map(j=>Math.min(...points.map(p=>p[j]))),max=[0,1,2].map(j=>Math.max(...points.map(p=>p[j])));v.target=min.map((n,j)=>(n+max[j])/2);const z=SM.norm([Math.cos(v.pitch)*Math.sin(v.yaw),Math.sin(v.pitch),Math.cos(v.pitch)*Math.cos(v.yaw)]),x=SM.norm(SM.cross([0,1,0],z)),y=SM.cross(z,x),ar=Math.max(.35,this.canvas.clientWidth/Math.max(1,this.canvas.clientHeight)),tv=Math.tan(v.fov/2),th=tv*ar;let d=.08;for(const p of points){const q=SM.sub(p,v.target),dep=SM.dot(q,z);d=Math.max(d,Math.abs(SM.dot(q,x))/th+dep,Math.abs(SM.dot(q,y))/tv+dep);}v.radius=d*1.17;v.saveView();v.invalidate();}
 orient(side){this.viewer.orient(side);this.fitAll();}
 setLayer(layer){this.layer=layer;if(this.mode==='anatomy'){this.viewer.visible=this.anatomyMask(this.action);this.fullMask=this.viewer.visible.slice();this.fitAll();}}
 isolate(){let i=this.viewer.selected;if(i<0)return false;this.viewer.visible=this.fullMask.map((v,k)=>k===i&&v);this.fitAll();return true;}
 restore(){this.viewer.visible=this.fullMask.slice();this.viewer.selected=-1;this.fitAll();}
}
g.ModelStage=ModelStage;
})(window);

/* Original teaching diagrams and explicitly synthetic examples. No patient estimation. */


/* V8.2 adapters. No mirroring of anatomy; left/right select existing source parts. */
(function(g){'use strict';
const P=ModelStage.prototype,oldMask=P.nativeMask,oldControls=P.controls,oldGuide=P.motionGuides;
P.nativeMask=function(a){
 const direct=oldMask.call(this,a);if(a.side!=='left')return direct;
 const byName=Object.fromEntries(this.native.bodies.map((b,i)=>[b.name,i]));
 const visibleBodies=new Set(this.native.parts.filter((p,i)=>direct[i]).map(p=>p.body));
 const swapped=new Set([...visibleBodies].map(i=>{const n=this.native.bodies[i].name;return n.endsWith('_r')?(byName[n.slice(0,-2)+'_l']??i):i;}));
 return this.native.parts.map(p=>swapped.has(p.body));
};
P.controls=function(a,t){let q=oldControls.call(this,a,t);if(a.side==='left')q=Object.fromEntries(Object.entries(q).map(([k,v])=>[k.endsWith('_r')?k.slice(0,-2)+'_l':k,v]));return q;};
P.motionGuides=function(a,t){const m=this.motions[a.motion];if(!m)return;const orig=m.joint;try{if(a.side==='left'&&orig.endsWith('_r'))m.joint=orig.slice(0,-2)+'_l';oldGuide.call(this,a,t);}finally{m.joint=orig;}};
P.supportsSide=function(a){const mo=this.motions[a.motion];return !!mo?.joint?.endsWith('_r');};

class TeachingClinic extends ClinicalScene{
 makeRoom(){super.makeRoom();this.coatIDs=[];const G=ClinicGeometry,white=[.975,.978,.972];
  // Continuous open-front coat surface instead of floating rectangular panels.
  const v=[],ix=[],n=28;for(let j=0;j<=n;j++){const q=.17+(Math.PI*2-.34)*j/n;for(let k=0;k<2;k++){const rx=k?.228:.208,rz=k?.173:.160;v.push([rx*Math.sin(q),k?.715:1.040,-.023+rz*Math.cos(q)]);}}
  for(let j=0;j<n;j++){const i=j*2;ix.push([i,i+1,i+2],[i+2,i+1,i+3]);}
  this.coatIDs.push(this.push(G.mesh('白大褂连续下摆',v,ix,white),{group:'coat',follow:'pelvis'}));
  for(const side of ['r','l']){const sign=side==='r'?-1:1;
   this.coatIDs.push(this.push(G.rod('白大褂上袖延伸',[sign*.214,1.230,-.016],[sign*.237,1.099,-.015],.046,white,18),{group:'coat',follow:'upper_'+side}));
   this.coatIDs.push(this.push(G.rod('白大褂前臂袖',[sign*.237,1.098,-.015],[sign*.274,.917,.006],.041,white,18),{group:'coat',follow:'fore_'+side}));
  }
 }
 makeActor(role){if(role==='therapist'){this.materials.therapist.shirt=[.97,.975,.97];this.materials.therapist.trim=[.83,.87,.88];this.materials.therapist.pants=[.17,.25,.31];}super.makeActor(role);}
 setCase(c){super.setCase(c);this.side=c.side||'right';this.active=true;
  this.roles.therapist.base=[-.76,0,.04];this.roles.patient.base=[.65,this.seated?-.415:0,.08];
  this.roles.therapist.yaw=Math.PI/2;this.roles.patient.yaw=-Math.PI/2;
  const center=[.58,0,.15];const turn=SM.rotate([0,1,0],-Math.PI/2,center);turn[12]+=.07;turn[14]-=.07;
  this.parts.forEach((p,i)=>{
   if(p.name.includes('治疗师凳')||p.name.includes('评估床')||p.name.includes('床架')||p.name.includes('床脚')||p.name.includes('绿植')||p.name==='花盆'||p.name==='教学墙卡线')this.viewer.visible[i]=false;
   else if(p.group==='chair'){this.viewer.visible[i]=this.seated;this.transforms[i]=turn;}
   else if(p.group==='plate'||p.group==='rails')this.viewer.visible[i]=false;
  });
  this.setPose('therapist','rest',0);this.setPose('patient','rest',0);this.frame('both');
 }
 poseAngles(kind,t,role){let prior=this.seated;this.seated=prior&&role==='patient';let out;try{out=super.poseAngles(kind,t,role);}finally{this.seated=prior;}
  if(kind==='knee'){out.q.thigh_r=[0,0,0];out.q.shin_r=[1.20*t,0,0];}
  if(this.side==='left'){
   const q={};for(const [name,vec]of Object.entries(out.q)){
    if(name.endsWith('_r'))q[name.slice(0,-2)+'_l']=[vec[0],-vec[1],-vec[2]];
    else if(name.endsWith('_l'))q[name.slice(0,-2)+'_r']=[vec[0],-vec[1],-vec[2]];
    else q[name]=vec;
   }out.q=q;
  }return out;
 }
 setPose(role,kind,t){super.setPose(role,kind,t);if(role!=='therapist'||!this.coatIDs)return;const a=this.roles.therapist;
  for(const i of this.coatIDs){const follow=this.parts[i].follow,bi=this.bmap[follow];if(a.palette[bi])this.transforms[i]=SM.matmul(a.transform,a.palette[bi]);}
  this.viewer.invalidate();
 }
 frame(mode='both',keep=false){super.frame(mode,keep);if(mode==='both'){this.viewer.target=[0,.96,.03];this.viewer.radius=Math.max(3.40,2.42/(Math.max(.35,this.canvas.clientWidth/Math.max(1,this.canvas.clientHeight))*2*Math.tan(this.viewer.fov/2)))+.15;this.viewer.yaw=.15;this.viewer.pitch=.06;this.viewer.saveView();this.viewer.invalidate();}}
}
g.TeachingClinic=TeachingClinic;
})(window);

/* Optional source tendon anchor reference: not wrapping or computed muscle force. */
(function(){const oldUpdate=ModelStage.prototype.update;
ModelStage.prototype.update=function(t){
 oldUpdate.call(this,t);
 if(this.mode!=='anatomy'&&this.action?.kind==='native'&&this.action.showMusclePaths){
  const names=this.motions[this.action.motion]?.paths||[];
  const paths=[];for(let name of names){
   // The right source names use no side suffix; the left copy is named *_l.
   if(this.action.side==='left')name=name+'_left';
   try{const points=this.kin.anchorPath(name);if(points&&points.length>1&&points.flat().every(Number.isFinite))paths.push({points,color:[.71,.26,.18],dots:true,overlay:true});}catch{}
  }
  this.viewer.paths=paths;this.viewer.invalidate();
 }
};})();


/* Source-mesh pelvis orientation experiment. No fabricated pelvis outline.
 * A prescribed rigid pelvis rotates around source hip centres. */
(function(){'use strict';
const P=ModelStage.prototype,oldMask=P.nativeMask,oldUpdate=P.update,oldSet=P.set;
P.nativeMask=function(a){if(a.kind!=='pelvis3d')return oldMask.call(this,a);const support=a.motion==='pelvisSupport';return this.native.parts.map(p=>[1,2,88,89,95].includes(p.body)||(support&&p.body>=90&&p.body<=100));};
P.pelvisPose=function(a,t){const support=a.motion==='pelvisSupport',side=a.supportSide||'right',other=side==='right'?'l':'r';
 const q={};if(support){q['hip_flexion_'+other]=.28;q['knee_angle_'+other]=.78;}
 this.kin.pose(q);const T=this.kin.transforms[88],H0=this.kin.joint('hip_flexion_r').pos,H1=this.kin.joint('hip_flexion_l').pos;
 const centre=support?(side==='right'?H0:H1):SM.mul(SM.add(H0,H1),.5),u=2*SM.clamp(t,0,1)-1,deg=u*(support?8:12),axis=support?SM.norm(SM.vec(T,[1,0,0])):SM.norm(SM.sub(H1,H0));
 const R=SM.rotate(axis,deg*Math.PI/180*(support?(side==='right'?-1:1):1),centre),held=side==='right'?new Set([89,90,91,92,93,94]):new Set([95,96,97,98,99,100]);
 const transform=p=>[1,2,88].includes(p.body)||support&&!held.has(p.body)?SM.matmul(R,this.kin.transforms[p.body]):this.kin.transforms[p.body];
 const transforms=this.native.parts.map(transform);
 // Landmarks are source-mesh vertices selected for clear visual orientation,
 // not palpated ASIS/PSIS measurements and not patient-specific calibration.
 if(!this.pelvisSourceMarks){const raw=this.assets.native.raw;const find=(right,posterior=false)=>{const p=this.native.parts.find(p=>p.body===88&&p.name===(right?'r_pelvis':'l_pelvis')),v=new Float32Array(raw,p.positions,p.vertexCount*3);let best=null,score=-Infinity;
  for(let i=0;i<v.length;i+=3){const pt=[v[i],v[i+1],v[i+2]];if(pt[1]<.0||Math.abs(pt[2])<.045)continue;const z=(posterior?-pt[0]:pt[0])+.22*pt[1];if(z>score){score=z;best=pt;}}
  return best;};this.pelvisSourceMarks={right:find(true),left:find(false),posterior:find(true,true)};}
 const toWorld=p=>SM.point(T,p),toMoved=p=>SM.point(R,toWorld(p)),marks=this.pelvisSourceMarks,A=toMoved(marks.right),B=toMoved(marks.left),post=toMoved(marks.posterior);
 const mid=support?SM.mul(SM.add(toWorld(marks.right),toWorld(marks.left)),.5):SM.mul(SM.add(toWorld(marks.right),toWorld(marks.posterior)),.5),hdir=support?[1,0,0]:SM.norm([SM.vec(T,[1,0,0])[0],0,SM.vec(T,[1,0,0])[2]]),horizontal=[SM.add(mid,SM.mul(hdir,-.24)),SM.add(mid,SM.mul(hdir,.24))];
 this.viewer.transforms=transforms;this.viewer.paths=[];this.viewer.planes=[];
 this.viewer.axes=[{points:horizontal,color:[.20,.43,.78],overlay:true},{points:support?[A,B]:[post,A],color:[.77,.42,.24],dots:true,overlay:true}];
 this.viewer.paths=[{points:[A,A],color:[.77,.42,.24],dots:true,overlay:true},{points:[B,B],color:[.77,.42,.24],dots:true,overlay:true}];
 const fixedHip=side==='right'?H0:H1;this.pelvisState={angle:deg,support,side,centre:Array.from(centre),markers:{right:A,left:B,posterior:post},transformation:Array.from(R),hipResidual:Math.hypot(...SM.sub(SM.point(R,fixedHip),fixedHip)),bilateralHipResidual:support?null:Math.max(...[H0,H1].map(h=>Math.hypot(...SM.sub(SM.point(R,h),h)))),label:support?(Math.abs(deg)<.2?'接近水平':deg>0?'对侧骨盆下降':'对侧骨盆上提'):(Math.abs(deg)<.2?'模型中立位':deg>0?'骨盆前倾':'骨盆后倾')};
 this.viewer.invalidate();};
P.update=function(t){if(this.action?.kind==='pelvis3d'&&this.mode!=='anatomy'){this.t=t;this.pelvisPose(this.action,t);return;}oldUpdate.call(this,t);};
P.set=function(a,t,mode){oldSet.call(this,a,t,mode);if(a.kind==='pelvis3d'){this.viewer.yaw=a.motion==='pelvisTilt'?Math.PI/2:0;this.viewer.pitch=0;this.fitAll();}};
})();


(function(g){'use strict';
const C={ink:'#203c3a',muted:'#58736d',green:'#167669',red:'#b8504b',blue:'#4276a3',bone:'#e6d8bf',line:'#c9dcd2',paper:'#f5f9f5'};
const clamp=(v,a,b)=>Math.max(a,Math.min(b,v)),esc=s=>String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const fmt=(v,n=1)=>Number.isFinite(v)?v.toFixed(n):'—';
const text=(x,y,s,size=20,color=C.ink,anchor='start')=>`<text x="${x}" y="${y}" fill="${color}" font-size="${size}" text-anchor="${anchor}">${esc(s)}</text>`;
const line=(a,b,color=C.green,width=4,dash='')=>`<line x1="${a[0]}" y1="${a[1]}" x2="${b[0]}" y2="${b[1]}" stroke="${color}" stroke-width="${width}" stroke-linecap="round" ${dash?'stroke-dasharray="'+dash+'"':''}/>`;
const arrow=(a,b,color=C.green,width=4)=>{const d=Math.atan2(b[1]-a[1],b[0]-a[0]),len=12;return line(a,b,color,width)+`<path d="M${b} L${b[0]-len*Math.cos(d-.5)},${b[1]-len*Math.sin(d-.5)} M${b} L${b[0]-len*Math.cos(d+.5)},${b[1]-len*Math.sin(d+.5)}" fill="none" stroke="${color}" stroke-width="${width}" stroke-linecap="round"/>`;};
const dot=(p,r=6,c=C.green)=>`<circle cx="${p[0]}" cy="${p[1]}" r="${r}" fill="${c}"/>`;
const svg=(label,body,h=480)=>`<svg viewBox="0 0 820 ${h}" role="img" aria-label="${esc(label)}" xmlns="http://www.w3.org/2000/svg"><rect width="820" height="${h}" rx="20" fill="${C.paper}"/><g font-family="system-ui,Microsoft YaHei,sans-serif">${body}</g></svg>`;
const card=(x,y,w,h)=>`<rect x="${x}" y="${y}" width="${w}" height="${h}" rx="14" fill="#fff" stroke="${C.line}"/>`;
const poly=(pts,c=C.green,w=4)=>`<polyline points="${pts.map(p=>p.join(',')).join(' ')}" fill="none" stroke="${c}" stroke-width="${w}" stroke-linecap="round" stroke-linejoin="round"/>`;
const curve=(fn,t,x,y,w,h,color,partial=false)=>poly(Array.from({length:81},(_,i)=>{const q=i/80*(partial?t:1);return[x+w*q,y+h/2-h*.44*fn(q)];}),color,3);
const bodySeg=(a,b,w=16,c=C.ink)=>line(a,b,c,w);
const Data={
 lever(kind,t,s={}){const load=Number.isFinite(+s.load)?clamp(+s.load,5,150):30;let de,dr;if(kind==='head'){de=.038+.008*Math.sin(Math.PI*t);dr=.028+.052*t;}else if(kind==='heel'){de=.16-.02*Math.sin(Math.PI*t);dr=.065+.03*t;}else{de=.035+.01*Math.sin(Math.PI*t);dr=.12+.18*t;}return{load,effort:load*dr/de,effortArm:de,resistanceArm:dr,torque:load*dr,advantage:de/dr};},
 com(t){const vertical=.95*Math.cos(4*Math.PI*t-0.45);const lateral=.75*Math.sin(2*Math.PI*t-0.2);return{vertical,lateral};},
 arch(t,s={}){const h0=({low:26,typical:43,high:61})[s.morphology]||43,k=clamp(+s.stiffness||1,.6,2),load=Math.sin(Math.PI*t)**2;return{h0,load,deflection:14*load/k,height:h0-14*load/k,stiffness:k,energy:.5*k*(14*load/k)**2};},
 cardio(kind,t){if(kind==='rest')return{hr:74+Math.round(Math.sin(4*Math.PI*t)),rr:14,bp:'118 / 76',speed:0,time:Math.round(60*t),label:'安静观察',gait:0,safety:'检查安全夹、鞋带与急停键'};if(kind==='walk')return{hr:Math.round(74+50*(1-Math.exp(-3*t))),rr:Math.round(14+12*(1-Math.exp(-3*t))),bp:'142 / 78',speed:3.5,time:Math.round(180*t),label:'跑步机行走',gait:(t*5)%1,safety:'安全夹已连接 · 旁站保护中'};return{hr:Math.round(74+48*Math.exp(-3*t)),rr:Math.round(14+9*Math.exp(-3*t)),bp:'124 / 76',speed:0,time:Math.round(180*t),label:'恢复与观察',gait:0,safety:'运动已停止，继续观察恢复反应'};},
 hr(age,rest,intensity,maxValue=null){if(!Number.isFinite(age)||age<18||age>90||!Number.isFinite(rest)||rest<35||rest>130||!Number.isFinite(intensity)||intensity<.2||intensity>.85)throw Error('请使用例题范围内的数值；这不是个人运动处方。');const maximum=maxValue===null?220-age:maxValue;if(!Number.isFinite(maximum)||maximum<=rest||maximum>240)throw Error('最大心率必须大于静息心率，且应说明来源。');return{maximum,reserve:maximum-rest,target:rest+(maximum-rest)*intensity,source:maxValue===null?'220－年龄的粗略预测':'用户提供的最大心率'};},
 index(rear,mid,fore){if(![rear,mid,fore].every(x=>Number.isFinite(x)&&x>=0)||rear+mid+fore===0)throw Error('三分区面积须非负，总面积须大于0。');return mid/(rear+mid+fore);},
 parseCOP(csv,label='导入数据'){if(csv.length>1000000)throw Error('CSV须小于1 MB。');const rows=csv.replace(/^\uFEFF/,'').trim().split(/\r?\n/);if(rows.length<4||rows.length>5001)throw Error('需要3—5000个采样点。');const keys=rows[0].split(',').map(x=>x.trim());const req=['time_s','x_mm','y_mm','fz_N'];if(req.some(x=>!keys.includes(x)))throw Error('列名应为 time_s,x_mm,y_mm,fz_N。');const out=rows.slice(1).filter(x=>x.trim()).map((r,i)=>{const v=r.split(',');const p=Object.fromEntries(req.map(k=>{const raw=v[keys.indexOf(k)];return[k,raw!==undefined&&raw.trim()!==''?Number(raw):NaN]}));if(!req.every(k=>Number.isFinite(p[k]))||p.fz_N<0||p.time_s<0||Math.abs(p.x_mm)>2000||Math.abs(p.y_mm)>2000)throw Error(`第${i+2}行：含无效数值或超出读取范围。`);return p;});for(let i=1;i<out.length;i++)if(out[i].time_s<=out[i-1].time_s)throw Error('时间必须严格递增。');if(!out.some(p=>p.fz_N>0))throw Error('全程无正向接触力，不能作为COP轨迹。');return{label:String(label).slice(0,100),status:'用户导入；来源和坐标尚待教师复核',rows:out};}
};
function lever(a,t,s){const d=Data.lever(a.motion,t,s),head=a.motion==='head',heel=a.motion==='heel';const title=head?'第一类杠杆：头颈前屈 / 后伸':heel?'第二类杠杆：提踵':'第三类杠杆：屈肘';let body=text(30,40,title,25)+text(30,72,'请跟随滑条，看动作过程中动力臂 dE 与阻力臂 dR 的变化。',18,C.muted);const F=head?[400,304]:heel?[198,332]:[200,295];const E=head?[303,255]:heel?[505-38*t,194]:[280,268];const R=head?[504,208+10*t]:heel?[335+118*t,307]:[535+95*t,287];
if(head){body+=`<path d="M312 246 Q290 198 308 126 Q340 65 428 67 Q509 70 529 128 L552 158 L529 170 L531 213 Q481 244 408 232 L408 275 L354 273Z" fill="${C.bone}" stroke="#9c8b73" stroke-width="3"/><path d="M404 242 Q435 288 418 354 Q383 386 382 406" fill="none" stroke="#b9ab92" stroke-width="20"/>`+poly([[302,314],[307,274],[340,220]],C.green,10);}
else if(heel){body+=`<g transform="translate(90,10)"><path d="M130 291 Q158 260 214 282 L288 307 Q360 300 395 273 Q432 302 486 312 L537 321 Q573 330 553 347 L166 347 Q126 341 130 291Z" fill="${C.bone}" stroke="#9c8b73" stroke-width="3"/></g>`+poly([[465-30*t,179],[457-15*t,214],[452-8*t,250]],C.green,10)+line([126,356],[604,356],C.line,3);}
else{body+=`<path d="M182 126 Q205 107 230 126 L238 252 Q258 276 234 300 L194 302 Q172 279 185 256Z" fill="${C.bone}" stroke="#9c8b73" stroke-width="3"/><path d="M226 282 Q258 265 294 280 L666 280 Q696 280 700 295 Q699 308 668 308 L244 311 Q218 306 226 282Z" fill="${C.bone}" stroke="#9c8b73" stroke-width="3"/>`+poly([[220,125],[250,168],[252,225],[281,277]],C.green,12);}
body+=arrow([E[0],E[1]+(head?-92:90)],[E[0],E[1]],C.green)+arrow([R[0],R[1]-84],[R[0],R[1]],C.red)+dot(F,8,C.blue)+`<path d="M${F[0]},${F[1]+8} l-15,26 h30Z" fill="${C.blue}"/>`;
body+=text(E[0],head?E[1]-104:E[1]+122,'E 动力',17,C.green,'middle')+text(R[0],R[1]-94,'R 阻力 / 负荷',17,C.red,'middle')+text(F[0]+6,F[1]+52,'F 支点',17,C.blue);
body+=line(F,E,'#7b9d93',2,'5 5')+line(F,R,'#c6a29c',2,'5 5');
body+=text((F[0]+E[0])/2-10,(F[1]+E[1])/2-10,'dE',18,C.green)+text((F[0]+R[0])/2+10,(F[1]+R[1])/2-10,'dR',18,C.red);
body+=card(32,430,756,118)+text(52,462,`动力臂 dE = ${(d.effortArm*1000).toFixed(1)} mm`,20,C.green)+text(286,462,`阻力臂 dR = ${(d.resistanceArm*1000).toFixed(1)} mm`,20,C.red)+text(542,462,`机械优势 dE/dR = ${d.advantage.toFixed(2)}`,20,C.ink)+text(52,500,`平衡所需动力 ≈ ${d.effort.toFixed(1)} N（教学示意）`,19)+text(52,532,head?'头颈前屈 / 后伸时，重点看头部重力线与颈后肌作用线怎样改变。':heel?'提踵时，重点看跟腱牵拉与前足接触区形成的力臂关系。':'屈肘时，重点看肱二头肌动力臂与手中负荷阻力臂的对比。',17,C.muted);
return svg(a.title,body,570);}

function bone(a,t){const map={compression:'压缩',tension:'拉伸',shear:'剪切',bending:'弯曲',torsion:'扭转'};const label=map[a.motion]||a.title;let body=text(30,40,`骨受载：${label}`,25)+text(30,72,'左侧看受力与形变，右侧看应力—应变或相应关系曲线。',18,C.muted);body+=card(28,98,360,390)+card(426,98,360,390);const top=154,bot=410,h=bot-top,m=36*t;let specimen='';if(a.motion==='compression'){specimen=`<rect x="170" y="${top+m/2}" width="78" height="${h-m}" rx="14" fill="${C.bone}" stroke="#a68e6e" stroke-width="3"/>`+arrow([209,112],[209,top+m/2-5],C.red)+arrow([209,452],[209,bot-m/2+5],C.red)+text(84,470,'长度缩短，横向膨胀示意',18,C.muted);}else if(a.motion==='tension'){specimen=`<rect x="178" y="${top-m/2}" width="62" height="${h+m}" rx="14" fill="${C.bone}" stroke="#a68e6e" stroke-width="3"/>`+arrow([209,top-m/2-6],[209,110],C.red)+arrow([209,bot+m/2+6],[209,454],C.red)+text(84,470,'长度增加，横向变窄示意',18,C.muted);}else if(a.motion==='shear'){specimen=`<polygon points="170,${top} 248,${top} 228,${bot} 150,${bot}" fill="${C.bone}" stroke="#a68e6e" stroke-width="3"/>`+arrow([104,160],[170+m,160],C.red)+arrow([312,404],[240-m,404],C.red)+text(88,470,'上下层相对错动',18,C.muted);}else if(a.motion==='bending'){specimen=`<path d="M154 146 Q${134-26*t} 285 154 422 L250 422 Q${270-26*t} 285 250 146Z" fill="${C.bone}" stroke="#a68e6e" stroke-width="3"/>`+arrow([320,282],[250-10*t,282],C.red)+text(84,470,'一侧受压，另一侧受拉',18,C.muted);}else{specimen=`<rect x="171" y="${top}" width="76" height="${h}" rx="14" fill="${C.bone}" stroke="#a68e6e" stroke-width="3"/><ellipse cx="209" cy="${top}" rx="39" ry="14" fill="none" stroke="${C.red}" stroke-width="3" transform="rotate(${26*t},209,${top})"/>`+Array.from({length:6},(_,i)=>poly([[180+i*12,top],[184+i*10,280],[180+(5-i)*12,bot]],'#b98650',2)).join('')+arrow([292,165],[264,136],C.red)+text(92,470,'端面相对旋转，出现扭角',18,C.muted);}body+=specimen+text(62,126,'受力与形变',21,C.ink);body+=text(456,126,'应力—应变 / 力学关系',21,C.ink)+line([485,430],[750,430],C.ink,2)+line([485,430],[485,158],C.ink,2)+text(758,438,'应变 ε',17)+text(470,154,'应力 σ',17);const curvePath='M495 422 C540 390 568 320 602 280 S676 208 736 178';body+=`<path d="${curvePath}" fill="none" stroke="${C.green}" stroke-width="4"/>`;const px=495+220*t,py=422-205*Math.pow(t,.8);body+=dot([px,py],8,C.red)+text(px+12,py-8,'当前阶段',16,C.red)+text(506,402,'弹性区',16,C.muted)+text(626,242,'更大载荷下形变增加',16,C.muted)+text(506,470,'应力 = 力 / 面积；应变 = 形变量 / 原长',18)+text(506,500,'刚度反映抵抗变形，强度反映承载而不失效。',18)+text(506,530,'课堂用于观察本质，不输出骨折阈值。',17,C.muted);return svg(a.title,body,560);}
function chain(a,t,s){const fixed=s.chain==='fixed';const title=fixed?'闭链：远端固定':'开链：远端游离';let body=text(30,40,title,24)+text(30,72,'同样的外观看起来都像“在动”，但分析对象并不相同。',18,C.muted);body+=line([160,398],[660,398],C.line,4);body+=`<path d="M258 368 Q280 343 326 353 L390 376 Q429 368 460 344 Q490 378 546 385 L576 390 Q600 397 586 410 L300 410 Q260 406 258 368Z" fill="${C.bone}" stroke="#9c8b73" stroke-width="3"/>`;if(fixed){body+=poly([[352,186],[410,272],[454,396]],C.green,16)+poly([[352,186],[332,132]],'#8fa7a0',12)+text(42,140,'闭链：足固定，小腿 / 身体相对足运动',20);}else{body+=poly([[352,186],[352,276],[418-85*t,354]],C.green,16)+text(42,140,'开链：远端游离，足相对小腿运动',20);}body+=card(34,436,744,90)+text(54,468,fixed?'课堂表达：小腿相对足背伸 / 跖屈，属于闭链任务。':'课堂表达：足相对小腿背伸 / 跖屈，属于开链任务。',18)+text(54,502,'关键是先说明“哪一段相对哪一段运动”。',18,C.muted);return svg(a.title,body,540);}
function breath(a,t){const f=(1-Math.cos(2*Math.PI*t))/2,forced=a.motion==='forced';const chestMode=a.motion==='thoracic',abdMode=a.motion==='abdominal';const chest=chestMode?26:14;const chestNow=chest*f;let body=text(30,40,forced?'用力呼气：主动呼气肌参与':chestMode?'胸式呼吸观察':abdMode?'腹式呼吸观察':'安静呼吸周期',24)+text(30,72,'把胸廓变化与腹壁 / 膈肌变化分开观察。',18,C.muted);body+=card(30,98,360,348)+card(428,98,360,348);body+=text(52,128,'侧面示意',20)+`<path d="M180 140 C${120-chestNow} 150 ${118-chestNow} 300 160 378 Q220 405 ${280+chestNow} 378 C${322+chestNow} 300 ${320+chestNow} 150 180 140Z" fill="#e7f1ec" stroke="#7ca295" stroke-width="3"/>`+Array.from({length:4},(_,i)=>`<path d="M${135-chestNow} ${188+i*38} Q180 ${202+i*34} ${225+chestNow} ${188+i*38}" fill="none" stroke="#c9bea9" stroke-width="9"/>`).join('')+`<path d="M${134-chestNow} 371 Q180 ${332-24*f} ${226+chestNow} 371" fill="none" stroke="${C.red}" stroke-width="9"/>`+line([254,336],[226+chestNow,365],C.red,2)+text(258,336,'膈肌',18,C.red)+arrow([85,206],[122-chestNow,206],C.green,3)+arrow([275,206],[236+chestNow,206],C.green,3)+text(76,188,chestMode?'胸廓扩张更明显':'胸廓变化',17,C.green);body+=text(450,128,'胸式 / 腹式对比',20)+card(448,158,150,118)+card(616,158,150,118)+text(523,188,'胸式',20,C.ink,'middle')+text(691,188,'腹式',20,C.ink,'middle')+arrow([478,228],[558,228],C.green,3)+text(523,250,'胸廓位移较明显',16,C.muted,'middle')+arrow([646,228],[726,228],C.green,3)+text(691,250,'腹壁位移较明显',16,C.muted,'middle');body+=card(448,300,318,116)+text(468,332,forced?'用力呼气时，腹肌等主动参与，帮助胸廓回缩。':'安静呼气主要依靠弹性回缩，不等于用力呼气。',18)+text(468,364,chestMode?'当前重点：胸廓提升与下降。':abdMode?'当前重点：腹壁隆起与回落。':'当前重点：一个完整呼吸周期。',18,C.green)+text(468,394,'说明卡：图示为功能示意，不作为肌纤维测量。',16,C.muted);return svg(a.title,body,470);}
function airflow(a,t){const dist=90*t;const hold=8*t;let body=text(30,40,'吹纸距离与持续时间',24)+text(30,72,'距离侧重输出效果，持续时间侧重持续输出表现。',18,C.muted);body+=card(30,108,360,320)+card(430,108,360,320);body+=text(52,140,'任务一：距离',20)+`<path d="M86 231 q34 -30 58 0 q-28 26 -58 0" fill="#d4b08e"/>`+Array.from({length:3},(_,i)=>arrow([160,198+i*28],[255+dist,198+i*28],C.green,3)).join('')+`<rect x="266" y="162" width="70" height="142" rx="12" fill="#fff" stroke="#91ada2" stroke-width="3"/><rect x="350" y="192" width="26" height="112" rx="4" fill="#dfccb1"/><path d="M363 192 Q${350+14*Math.sin(Math.PI*t)} 162 ${364+12*Math.sin(Math.PI*t)} 133 Q386 172 368 192" fill="#dfa24c"/>`+line([250,332],[368,332],C.ink,2)+text(248,356,'距离',16)+text(352,356,`${Math.round(25+35*t)} cm`,16)+text(52,390,'吹纸 / 虚拟蜡烛的距离任务：主要观察输出效果。',17,C.muted);body+=text(452,140,'任务二：持续时间',20)+`<path d="M486 231 q34 -30 58 0 q-28 26 -58 0" fill="#d4b08e"/>`+Array.from({length:3},(_,i)=>arrow([560,198+i*28],[650,198+i*28],C.green,3)).join('')+`<rect x="660" y="162" width="84" height="142" rx="12" fill="#fff" stroke="#91ada2" stroke-width="3" transform="rotate(${-24*Math.sin(Math.PI*t)},702,233)"/>`+card(500,320,224,76)+text(520,350,`纸张持续飘动时间：${hold.toFixed(1)} s`,18,C.red)+text(520,382,'持续时间任务：主要观察持续输出表现。',17,C.muted);body+=text(32,455,'课堂任务指标不直接等于正式呼吸肌力或耐力测量，也不安排真人明火操作。',17,C.muted);return svg(a.title,body,485);}
function gait(a,t){const p=t,stage=p<.1?'初始双支撑':p<.5?'右单支撑 / 左摆动':p<.6?'末期双支撑':'左单支撑 / 右摆动';const com=Data.com(t);let body=text(30,40,'步态周期、支撑状态与全身质心',24)+text(30,72,'以同侧初始接触到下一次同侧初始接触为一个周期。',18,C.muted);body+=card(28,100,360,392)+card(410,100,382,392);const floor=410,cx=206,cy=250;const rightSwing=p>.6,leftSwing=(p>.1&&p<.5);const rk=[cx+12+(rightSwing?45*Math.sin((p-.6)/.4*Math.PI):0), floor-90+(rightSwing?-20*Math.sin((p-.6)/.4*Math.PI):0)],ra=[cx+20+(rightSwing?88*Math.sin((p-.6)/.4*Math.PI):0), floor-(rightSwing?55*Math.sin((p-.6)/.4*Math.PI):0)],lk=[cx-22+(leftSwing?-44*Math.sin((p-.1)/.4*Math.PI):0), floor-88+(leftSwing?-18*Math.sin((p-.1)/.4*Math.PI):0)],la=[cx-30+(leftSwing?-84*Math.sin((p-.1)/.4*Math.PI):0), floor-(leftSwing?55*Math.sin((p-.1)/.4*Math.PI):0)];body+=line([70,floor],[340,floor],C.line,4)+line([cx,cy],[cx,170],C.green,18)+dot([cx-2,140],24,'#c69a78')+poly([[cx,cy],[rk[0],rk[1]],[ra[0],ra[1]]],C.ink,16)+poly([[cx,cy],[lk[0],lk[1]],[la[0],la[1]]],'#9fb6ad',14)+line([ra[0],ra[1]],[ra[0]+36,ra[1]+6],C.ink,12)+line([la[0],la[1]],[la[0]+36,la[1]+6],'#9fb6ad',10)+poly([[cx,186],[cx+28,220],[cx+48,260]],'#9fb6ad',10)+poly([[cx,186],[cx-34,222],[cx-56,260]],C.ink,10)+dot([cx+24,cy-12*com.vertical],8,C.red)+text(240,157,'● 全身质心（教学信号）',17,C.red)+text(52,126,'侧面行走示意',20)+text(52,458,stage,18,C.green);body+=text(430,130,'参考腿：右侧',20)+text(430,164,`周期进度 ${Math.round(100*t)}%`,18)+text(430,195,'足接触条',18,C.ink)+line([456,232],[750,232],C.line,12)+line([456,266],[750,266],C.line,12);const mapx=q=>456+294*q;const seg=(y,a,b,color)=>line([mapx(a),y],[mapx(b),y],color,12);seg(232,0,.6,C.green);seg(266,0,.1,C.blue);seg(266,.5,1,C.blue);dot([mapx(t),232],8,C.red);dot([mapx(t),266],8,C.red);body+=text(402,238,'右足',17)+text(402,272,'左足',17)+text(456,300,'0%',15)+text(mapx(.1),300,'10%',15,'#58736d','middle')+text(mapx(.5),300,'50%',15,'#58736d','middle')+text(mapx(.6),300,'60%',15,'#58736d','middle')+text(750,300,'100%',15,C.ink,'end');body+=card(430,324,324,132)+text(450,352,'全身质心变化',18)+text(450,380,`上下位移：${com.vertical.toFixed(2)}（一周期通常出现两次起伏）`,17,C.red)+text(450,408,`左右位移：${com.lateral.toFixed(2)}（向支撑侧转移）`,17,C.blue)+text(450,436,'异常步态线索请在“肌群时序与异常步态线索”页对照。',16,C.muted);return svg(a.title,body,520);}
function cardio(a,t){const d=Data.cardio(a.motion,t); const moving=a.motion==='walk'; const belt=110+180*(moving?d.gait:0); let body=text(30,40,'跑步机：运动反应观察',24)+text(30,72,'保留跑步机，并加入安全保护、急停与风险应对提示。',18,C.muted); body+=card(24,104,370,360)+card(418,104,372,360); body+=text(46,136,'设备与人物',20)+`<rect x="98" y="378" width="252" height="42" rx="10" fill="#5e7570"/>`+Array.from({length:7},(_,i)=>line([118+((i*34+belt)%228),390],[132+((i*34+belt)%228),408],'#a8b8b1',3)).join('')+line([324,376],[324,184],C.green,6)+line([332,184],[358,184],C.green,6)+`<circle cx="358" cy="184" r="8" fill="${C.red}"/>`+text(326,172,'急停',16,C.red)+line([130,376],[130,140],C.line,4)+line([130,140],[282,140],C.line,4)+line([282,140],[282,376],C.line,4)+line([206,172],[206,206],C.line,4)+dot([206,196],18,'#c69a78')+line([206,214],[206,290],C.green,28)+poly([[206,230],[176,274],[160,332]],'#8fa9a1',11)+poly([[206,230],[236,274],[252,332]],'#8fa9a1',11)+poly([[206,290],[182,344],[168,392]],C.ink,13)+poly([[206,290],[230,340],[244,392]],'#a8bbb4',11)+line([160,332],[190,336],'#8fa9a1',10)+line([252,332],[282,334],'#8fa9a1',10)+line([190,160],[190,270],C.red,2,'5 5')+text(142,160,'安全夹',16,C.red); body+=text(46,438,d.safety,17,C.muted);
body+=text(440,136,'监测与应对',20)+card(446,160,316,204)+text(468,196,d.label,24)+text(468,236,`心率 ${d.hr} 次/分`,21,C.red)+text(468,272,`呼吸 ${d.rr} 次/分`,21)+text(468,308,`血压 ${d.bp} mmHg`,20)+text(468,338,`速度 ${d.speed} km/h · 时间 ${d.time} s`,19)+card(446,378,316,64)+text(466,405,moving?'若出现头晕、胸闷、步态不稳：立即急停、扶持并反馈。':'无反应且无正常呼吸时，启动求助、CPR / AED的教学流程提示。',17);return svg(a.title,body,485);}
function arch(a,t,s){const d=Data.arch(t,s),base=345;let b=text(30,40,'足弓：形态、负重与等效刚度分开观察',23)+text(30,70,'外形低 / 常见 / 高不自动决定软硬，也不自动决定症状',17,C.muted);const hs=a.motion==='compare'?[26,43,61]:[d.h0],xs=a.motion==='compare'?[145,410,675]:[410];hs.forEach((h,i)=>{const x=xs[i],w=a.motion==='compare'?92:165,k=d.stiffness,def=14*d.load/k,hcur=h-def;b+=`<path d="M${x-w} ${base} Q${x} ${base-2*hcur} ${x+w} ${base}" fill="none" stroke="${C.bone}" stroke-width="18" stroke-linecap="round"/><path d="M${x-w} ${base} Q${x} ${base-2*h} ${x+w} ${base}" fill="none" stroke="#aebdb3" stroke-width="2" stroke-dasharray="6 6"/>`+line([x-w,base+15],[x+w,base+15],C.green,4)+arrow([x,base-hcur-65-35*d.load],[x,base-hcur-15],C.red)+text(x,base+56,a.motion==='compare'?['低弓形态','常见弓高示意','高弓形态'][i]:'实线：加载　虚线：起始',18,C.ink,'middle');});b+=text(36,135,`相对加载 ${fmt(d.load*100,0)}%　等效刚度 ${fmt(d.stiffness,1)}（模型单位）`,20)+text(36,169,t<.5?'加载：结构形变、储能':'卸载：回复、释放储能',20,C.green)+text(30,455,'未模拟组织黏弹性、主动肌肉调节、真实接触压力或临床足弓分类。',17,C.muted);return svg(a.title,b);}
function syntheticCOP(){return {label:'合成读图样例（不是某足型标准）',status:'synthetic',rows:Array.from({length:101},(_,i)=>{const t=i/100;return{time_s:t,x_mm:24*Math.sin(Math.PI*t)-18*t,y_mm:25+220*t,fz_N:i===0||i===100?0:600*Math.sin(Math.PI*t)**.45};})};}
function cop(a,t,s){
 const ds=s.copData||syntheticCOP(),rows=ds.rows,j=Math.min(rows.length-1,Math.floor(t*(rows.length-1))),p=rows[j],valid=rows.filter(r=>r.fz_N>0),imported=ds.status!=='synthetic';
 let minx=Math.min(...valid.map(r=>r.x_mm)),maxx=Math.max(...valid.map(r=>r.x_mm)),miny=Math.min(...valid.map(r=>r.y_mm)),maxy=Math.max(...valid.map(r=>r.y_mm));
 const span=Math.max(maxx-minx,maxy-miny,10)*1.18,cx=(minx+maxx)/2,cy=(miny+maxy)/2;
 minx=cx-span/2;maxx=cx+span/2;miny=cy-span/2;maxy=cy+span/2;
 const map=r=>[205+270*(r.x_mm-minx)/span,365-270*(r.y_mm-miny)/span];
 let b=text(30,40,'足底压力中心 COP：接触合力作用位置',23)+text(30,70,'与全身质心 COM 不同；足离地时不继续描画该足的COP',17,C.muted);
 b+=card(200,90,280,280);
 for(let i=0;i<=2;i++){const x=205+135*i,y=365-135*i;
 b+=line([x,95],[x,365],C.line,1)+line([205,y],[475,y],C.line,1);
 b+=text(x,386,fmt(minx+span*i/2,1),13,C.muted,'middle')+text(192,y+4,fmt(miny+span*i/2,1),13,C.muted,'end');}
 b+=text(487,371,'X / mm',14,C.muted)+text(200,86,'Y / mm',14,C.muted);
 let seg=[];for(let k=0;k<=j;k++){if(rows[k].fz_N>0)seg.push(map(rows[k]));else{if(seg.length>1)b+=poly(seg,C.green,4);seg=[];}}
 if(seg.length>1)b+=poly(seg,C.green,4);if(p.fz_N>0)b+=dot(map(p),8,C.red);
 b+=text(555,137,`时间 ${fmt(p.time_s,2)} s`,20)+text(555,177,p.fz_N>0?`X ${fmt(p.x_mm)} mm`:'离地 / 无有效COP',20)+text(555,215,p.fz_N>0?`Y ${fmt(p.y_mm)} mm`:'不延长上一接触轨迹',18)+text(555,252,`Fz ${fmt(p.fz_N,0)} N`,19);
 b+=text(555,296,'X、Y等比例显示',17,C.green)+text(555,326,'不自动配准到足形轮廓',16,C.muted);
 b+=text(30,420,ds.label,17,C.green)+text(30,456,imported?'导入数据：坐标方向、侧别、仪器与阈值需复核；未过滤噪声。':'合成坐标读图样例：不代表正常、低弓或高弓足的标准轨迹。',16,C.muted);
 return svg(a.title,b);
}
function scapula(a,t,s){const mode=s.scap||'elevation',tx=mode==='protraction'?60*(t-.5):0,ty=mode==='elevation'?-70*(t-.5):0,rot=mode==='rotation'?-45*(t-.5):0;return svg(a.title,`${text(30,40,'肩胛胸壁运动 · 参考胸廓',24)}${text(30,70,'功能性连接的方向示意；并非精确贴合或肩肱节律数据',17,C.muted)}<ellipse cx="410" cy="256" rx="112" ry="150" fill="#e4eee7" stroke="#b1c5ba" stroke-width="3"/>${line([410,109],[410,408],C.line,8)}<g transform="translate(${tx},${ty}) rotate(${rot},457,202)"><path d="M431 157 L532 187 L450 318 Q426 262 431 157" fill="${C.bone}" stroke="#93836c" stroke-width="3"/>${line([433,183],[521,191],'#ac9777',5)}</g>${text(50,453,mode==='elevation'?'观察上下方向：上提 / 下降':mode==='protraction'?'前伸 / 后缩是绕胸廓的空间运动；图中仅投影说明':'观察上回旋 / 下回旋方向；不设固定比例',19)}`);}
function pelvis(a,t){const support=a.motion==='support',ang=(t-.5)*(support?18:28);return svg(a.title,`${text(30,40,support?'单侧支持：观察骨盆相对水平线':'骨盆前后倾：相对固定参考',24)}${text(30,70,'方向示意；外观变化不独立确定病因或肌肉无力',17,C.muted)}${line([225,260],[596,260],'#a1b7aa',3,'7 7')}<g transform="rotate(${ang},410,260)"><path d="M302 231 Q320 191 370 228 L410 247 L450 228 Q501 191 518 231 L485 294 Q453 328 410 284 Q367 328 335 294Z" fill="${C.bone}" stroke="#9f8c72" stroke-width="4"/>${dot([341,263],7,C.red)}${dot([479,263],7,C.red)}</g>${bodySeg([479,290],[484,401],22,C.bone)}${support?bodySeg([342,290],[323,365],22,C.bone):bodySeg([342,290],[338,401],22,C.bone)}${text(31,454,support?'保持支撑条件记录；不把全部骨盆参与自动叫作代偿。':'前后倾是骨盆整体定向；与髋角和腰椎曲度分开描述。',18)}`);}
function ssc(a,t){const crouch=t<.5?1-t*1.2:(t-.5)*1.8; const y=310+42*Math.max(0,crouch), knee=390-44*Math.max(0,crouch), jump=t>.6?Math.sin((t-.6)/.4*Math.PI)*48:0; let body=text(30,40,'牵拉—缩短周期：起跳前下蹲到蹬伸离地',24)+text(30,72,'观察股四头肌与小腿三头肌如何先预牵拉，再迅速缩短发力。',18,C.muted); body+=line([92,420],[346,420],C.line,4)+dot([210,138-jump],22,'#c69a78')+line([210,160-jump],[210,240+y-jump],C.green,20)+poly([[210,184-jump],[176,222-jump],[150,274-jump]],'#9eb4ab',10)+poly([[210,184-jump],[244,222-jump],[268,274-jump]],'#9eb4ab',10)+poly([[210,240+y-jump],[196,300+knee-jump],[178,402-jump]],C.ink,14)+poly([[210,240+y-jump],[228,298+knee-jump],[246,402-jump]],C.ink,14)+line([178,402-jump],[210,405-jump],C.ink,10)+line([246,402-jump],[278,405-jump],C.ink,10);
body+=`<ellipse cx="196" cy="${322+knee/2-jump}" rx="18" ry="42" fill="#d66a5f" opacity="0.75"/><ellipse cx="248" cy="${348-jump}" rx="16" ry="38" fill="#d66a5f" opacity="0.75"/>`+text(82,146,'股四头肌',18,C.red)+line([152,150],[186,282+knee/2-jump],C.red,2)+text(264,146,'小腿三头肌',18,C.red)+line([314,150],[252,326-jump],C.red,2);
body+=card(390,126,394,262)+text(412,160,t<.45?'阶段1：起跳前下蹲（预牵拉）':t<.65?'阶段2：转换（时间短）':'阶段3：蹬伸离地（缩短发力）',22,C.green)+text(412,202,t<.45?'股四头肌与小腿三头肌被主动控制地拉长，储备弹性势能。':t<.65?'从离心到向心的转换越短，牵拉—缩短周期利用越充分。':'股四头肌伸膝、小腿三头肌跖屈，完成蹬伸与离地。',18)+text(412,252,'教学重点：不要只看“上下起跳”，要看肌肉工作方式的连续变化。',17,C.muted)+text(412,294,'对应课堂能力表现：爆发力可用摸高与滞空时间辅助观察。',17,C.muted); return svg(a.title,body,470);}
const Concept={capacity:[['肌力：1RM举重','以举起1RM重物的情境示意最大输出能力，重点看完成一次最大负荷任务所需的主要肌群。'],['爆发力：摸高与滞空','用摸高高度与滞空时间的课堂表现辅助观察快速发力能力；不是单纯看动作快慢。'],['肌耐力：平板支撑','以平板支撑保持时间示意持续输出能力，重点观察体位维持与疲劳后的质量下降。'],['肌耐力：深蹲次数','在规定节律下完成的深蹲次数，体现反复完成给定肌肉工作的能力。']],measure:[['测量前','明确体位、安静时段、仪器与计时方式。'],['测量中','区分心率、呼吸频率和血压，不互相替代。'],['测量后','保存时间点、单位、条件与主观感受。'],['实际技能','本页不验证真实袖带、听诊或传感器操作。']],curves:[['四个生理弯曲','颈曲、胸曲、腰曲、骶曲：先明确观察方向。'],['描述与解释分开','写清哪一方向、哪一部位、哪些证据。'],['相关肌群','伸肌、腹肌等按任务参与，不凭姿势判肌力。'],['侧弯观察','外观线索不等同完整诊断；不编造异常模型。']],bandage:[['先明确目的','支持、限制方向或本校规程所定的矫正目的。'],['解释与许可','说明接触范围与操作；不适时停止并反馈。'],['观察要点','松紧、皮肤与末梢情况，按本校规程核定。'],['尚未开放','无经核定路线，不将通用8字图当作矫正标准。']],gaitReviewPlus:[['先确定周期','以同侧初始接触到下一次同侧初始接触为一个周期。'],['支撑状态','明确双支撑、单支撑与摆动，不把瞬时事件和持续阶段混淆。'],['肌群时序','负重接受、支撑稳定和摆动清足阶段，主要工作肌群不同。'],['异常步态线索','可从抗痛步态、垂足步态、臀中肌无力相关步态等角度观察，但不直接给出诊断。']]};
function concept(a){const arr=Concept[a.motion]||[];return `<div class="concept-grid">${arr.map(([title,body])=>`<section><span class="micro">观察卡</span><h3>${esc(title)}</h3><p>${esc(body)}</p></section>`).join('')}</div>`;}
const Sim={Data,syntheticCOP,esc,render(a,t,s={}){switch(a.kind){case'lever':return lever(a,t,s);case'bone':return bone(a,t);case'chain':return chain(a,t,s);case'breath':return breath(a,t);case'airflow':return airflow(a,t);case'gait':return gait(a,t);case'cardio':return cardio(a,t);case'arch':return arch(a,t,s);case'cop':return cop(a,t,s);case'scapula':return scapula(a,t,s);case'pelvis':return pelvis(a,t);case'ssc':return ssc(a,t);case'concept':return concept(a);default:return'';}},hr:Data.hr};
g.Sim=Sim;if(typeof module!=='undefined')module.exports=Sim;
})(typeof window!=='undefined'?window:globalThis);



/* Original transparent classroom models. Every numeric display uses these functions.
   None of the constants below is a patient-specific measurement or injury threshold. */
(function(g){'use strict';
const clamp=(v,a=0,b=1)=>Math.max(a,Math.min(b,v));
const finite=(v,name='数值')=>{v=Number(v);if(!Number.isFinite(v))throw Error(name+'必须是有限数值');return v;};
const add=(a,b)=>a.map((v,i)=>v+b[i]),sub=(a,b)=>a.map((v,i)=>v-b[i]),mul=(a,s)=>a.map(v=>v*s),dot=(a,b)=>a.reduce((v,x,i)=>v+x*b[i],0),cross=(a,b)=>a[0]*b[1]-a[1]*b[0],rot=(p,q)=>[p[0]*Math.cos(q)-p[1]*Math.sin(q),p[0]*Math.sin(q)+p[1]*Math.cos(q)];
function lineFoot(F,P,u){const den=dot(u,u);if(den<1e-12)throw Error('作用线方向不能为零');return add(P,mul(u,dot(sub(F,P),u)/den));}
function lever(kind,t,s={}){
 t=clamp(finite(t));const load=clamp(finite(s.load??30),5,300),L=clamp(finite(s.distance??.32),.18,.42);let theta,E,R,O;
 if(kind==='head'){theta=(-(-15+45*t))*Math.PI/180;E=rot([-.046,.033],theta);R=rot([.055,.105],theta);O=[-.080,-.155];}
 else if(kind==='heel'){theta=-24*t*Math.PI/180;E=rot([-.205,.028],theta);R=rot([-.085,.075],theta);O=[-.15,.43];}
 else{theta=(18+102*t)*Math.PI/180;E=rot([0,-.044],theta);R=rot([0,-L],theta);O=[0,.29];}
 const F=[0,0],forceDir=sub(O,E),len=Math.hypot(...forceDir),u=mul(forceDir,1/len),r=[0,-1];
 const fe=lineFoot(F,E,u),fr=lineFoot(F,R,r),de=Math.hypot(...fe),dr=Math.hypot(...fr),signedE=cross(E,u),signedR=cross(R,r);
 const effort=Math.abs(signedE)>1e-7?-load*signedR/signedE:null;
 return {kind,t,theta,angle:kind==='head'?-theta*180/Math.PI:kind==='heel'?-theta*180/Math.PI:theta*180/Math.PI,F,E,R,O,u,r,fe,fr,effortArm:de,resistanceArm:dr,load,effort,torque:load*dr,advantage:dr>1e-7?de/dr:null,distance:L,balanceResidual:effort===null?null:effort*signedE+load*signedR};
}
function bone(kind,t,s={}){
 const u=clamp(finite(t)),v=s.loadCycle?(u<.5?u*2:(1-u)*2):u;
 const E=s.material==='B'?12000:6000,nu=.30,G=E/(2*(1+nu)),L=100,A=100,I=10000/12,J=10000/6;
 let load,x,y,xMax,yMax,xLabel,yLabel,displacement,info;
 if(kind==='tension'||kind==='compression'){
  load=600*v;const stress=load/A,eps=stress/E;displacement=eps*L;x=eps;y=stress;xMax=(600/A)/E;yMax=600/A;xLabel='轴向应变 |ε|';yLabel='轴向应力 |σ| / MPa';info='均匀轴向试件：σ=F/A，ε=ΔL/L。压缩时原始符号为负，本图画绝对值。';
 }else if(kind==='shear'){
  load=120*v;y=load/A;x=y/G;displacement=x*L;xMax=120/A/G;yMax=120/A;xLabel='剪应变 γ / rad';yLabel='平均剪应力 τ / MPa';info='理想均匀剪切：γ≈横向错动/原长；不等同于梁内真实剪应力分布。';
 }else if(kind==='bending'){
  load=5*v;y=load*L;x=y/(E*I);displacement=load*L**3/(3*E*I);xMax=500/(E*I);yMax=500;xLabel='固定端曲率 κ / mm⁻¹';yLabel='固定端弯矩 M / N·mm';info='小挠度悬臂梁，端部横向受力：κ=M/EI；内侧受压、外侧受拉。';
 }else{
  load=600*v;y=load;x=load*L/(G*J);displacement=x;xMax=600*L/(G*J);yMax=600;xLabel='相对扭角 θ / rad';yLabel='扭矩 T / N·mm';info='线弹性扭转：θ=TL/GJ。纵向网格的旋转来自同一扭角。';
 }
 return {kind,t:u,level:v,E,G,L,A,I,J,load,x,y,xMax,yMax,xLabel,yLabel,displacement,info,elastic:true,parameterSource:'原创线弹性例题，非人体骨参数'};
}
function capacity(mode,t,s={}){
 t=clamp(finite(t));
 if(mode==='rm'){
  const mass=clamp(finite(s.mass??12),8,14),successful=mass<=12,progress=successful?t:Math.min(t,.42);
  return {mode,phase:progress,angle:20+80*progress,mass,successful,validRep:successful&&t>=.99,caseMaximum:12,label:successful?(t>=.99?'完成一次规定动作':'举起中'):'本例无法完成全程',source:'虚构同一动作案例：8、10、12 kg成功；14 kg未完成'};
 }
 if(mode==='jump'||mode==='ssc'){
  const v=clamp(finite(s.takeoff??2.6),1.5,3.5),h=v*v/(2*9.81),flight=2*v/9.81;
  let squat=0,height=0,phase='准备',phi=0;
  if(t<.28){squat=t/.28;phase='下蹲制动';}
  else if(t<.40){squat=1-(t-.28)/.12;phase='转换 → 蹬伸';}
  else if(t<.82){phi=(t-.40)/.42;height=4*h*phi*(1-phi);phase='腾空';}
  else {squat=Math.sin((t-.82)/.18*Math.PI)*.55;phase='落地缓冲';}
  return {mode,t,squat,height,phase,flight,heightMax:h,standingReach:2.15,jumpReach:2.15+h,quad:phase==='下蹲制动'?'承受负荷 / 离心制动':phase==='转换 → 蹬伸'?'伸膝发力':phase==='落地缓冲'?'吸收冲击':'空中调整',calf:phase==='下蹲制动'?'控制背屈 / 肌腱参与':phase==='转换 → 蹬伸'?'跖屈推进':phase==='落地缓冲'?'控制着地':'空中调整',note:'阶段是功能解释；不将其当作每束肌纤维长度实测'};
 }
 const maxTime=mode==='plank'?45:30,elapsed=t*maxTime,quality=s.quality||'valid';
 if(mode==='plank'){const validUntil=quality==='valid'?45:22;return{mode,t,elapsed,validSeconds:Math.min(elapsed,validUntil),broken:elapsed>validUntil,label:elapsed>validUntil?'姿势要求已偏离：有效计时停止':'维持要求',phase:0};}
 const completed=Math.floor((elapsed+1e-9)/3),phase=(elapsed%3)/3;
 return{mode,t,elapsed,phase,squat:Math.sin(Math.PI*phase),total:completed,valid:quality==='valid'?completed:0,label:quality==='valid'?'每次完整下蹲并回到起始才计一次':'深度未达本题要求：不计有效重复'};
}
function kneeFromEnds(hip,ankle,L1=.42,L2=.43,bend=1){
 const H=[hip[0],hip[1],hip[2]||0],A=[ankle[0],ankle[1],ankle[2]||0],v=sub(A,H),dist=Math.hypot(...v);
 if(dist<1e-8||dist>L1+L2+1e-8||dist<Math.abs(L1-L2)-1e-8)throw Error('骨段目标超出两段运动链可达范围');
 const u=mul(v,1/dist),preferred=[1,0,0],b=sub(preferred,mul(u,dot(preferred,u))),bl=Math.hypot(...b);
 const n=bl>1e-9?mul(b,1/bl):[0,1,0],along=(L1*L1-L2*L2+dist*dist)/(2*dist),height=Math.sqrt(Math.max(0,L1*L1-along*along));
 return add(H,add(mul(u,along),mul(n,height*bend)));
}
function gait(t,pattern='normal'){
 t=clamp(finite(t));const periodic=t===1?0:t;const S=.9,shift=pattern==='short'?.42:.5,stanceR=pattern==='short'?.50:.6,stanceL=pattern==='short'?.73:.6;
 const hip=[S*t,.765-.012*Math.cos(4*Math.PI*periodic),.018*Math.sin(2*Math.PI*periodic)];
 const lean=pattern==='lean'?.065*Math.max(0,Math.sin(2*Math.PI*periodic)):0;
 const shoulder=[hip[0]-.012,hip[1]+.50,hip[2]+lean],head=[shoulder[0],shoulder[1]+.18,shoulder[2]];
 const legs={};
 for(const side of ['right','left']){
  const offset=side==='right'?0:shift,st=side==='right'?stanceR:stanceL,raw=t-offset,cycle=Math.floor(raw+1e-10),p=raw-cycle,contact=p<st||Math.abs(p)<1e-8;
  const stancePos=(cycle+offset)*S+.21,z=side==='right'?.095:-.095;
  let x=stancePos,y=.070,pitch=0;
  if(contact){
   if(p<.10){pitch=(1-p/.10)*.14;const heel=[-.055,-.070];const rr=rot(heel,pitch);x=stancePos+heel[0]-rr[0];y=-rr[1];}
   else if(p>st-.13){pitch=-.33*((p-(st-.13))/.13);const toe=[.17,-.070],rr=rot(toe,pitch);x=stancePos+toe[0]-rr[0];y=-rr[1];}
  }else{
   const u=(p-st)/(1-st),endPitch=.14,startPitch=-.33,pStart=rot([.17,-.07],startPitch),pEnd=rot([-.055,-.07],endPitch);
   const start=stancePos+.17-pStart[0],end=stancePos+S-.055-pEnd[0];
   const e=u*u*(3-2*u);x=start+(end-start)*e;
   y=(1-e)*(-pStart[1])+e*(-pEnd[1])+(pattern==='clearance'&&side==='right'?.012:.09)*Math.sin(Math.PI*u);
   pitch=startPitch+(endPitch-startPitch)*u;
  }
  const ankle=[x,y,z],h=[hip[0],hip[1],hip[2]+(side==='right'?.07:-.07)],k=kneeFromEnds(h,ankle,.42,.43,1);
  const heel=rot([-.055,-.070],pitch),toe=rot([.17,-.070],pitch);
  legs[side]={hip:h,knee:k,ankle,heel:[x+heel[0],y+heel[1],z],toe:[x+toe[0],y+toe[1],z],contact,phase:p,stance:st,pitch};
 }
 const arms={};for(const side of ['right','left']){const sg=side==='right'?1:-1;const a=.25*Math.sin(2*Math.PI*periodic)*sg;const s=[shoulder[0],shoulder[1]-.05,shoulder[2]+sg*.16];const elbow=[s[0]+Math.sin(a)*.27,s[1]-Math.cos(a)*.27,s[2]],wrist=[elbow[0]+Math.sin(a+.20)*.23,elbow[1]-Math.cos(a+.20)*.23,s[2]];arms[side]={s,elbow,wrist};}
 // Explicit teaching mass fractions sum to one. Model COM, not measured human COM.
 const weighted=[[head,.08],[mul(add(hip,shoulder),.5),.50]];
 for(const s of ['right','left']){const l=legs[s],a=arms[s];weighted.push([mul(add(l.hip,l.knee),.5),.10],[mul(add(l.knee,l.ankle),.5),.045],[mul(add(l.heel,l.toe),.5),.015],[mul(add(a.s,a.elbow),.5),.026],[mul(add(a.elbow,a.wrist),.5),.024]);}
 const com=weighted.reduce((acc,[v,m])=>add(acc,mul(v,m)),[0,0,0]);
 const R=legs.right.contact,L=legs.left.contact,support=R&&L?'双支撑':R?'右单支撑 · 左摆动':L?'左单支撑 · 右摆动':'双足离地';
 return {t,hip,shoulder,head,legs,arms,com,support,contact:[R,L],shift,stanceR,stanceL,massSum:weighted.reduce((s,v)=>s+v[1],0)};
}
function breath(t,kind='quiet'){
 t=clamp(finite(t));const expansion=(1-Math.cos(2*Math.PI*t))/2;return{phase:t<.5?'吸气':kind==='forced'?'主动用力呼气':'安静呼气',expansion,chest:expansion*(kind==='thoracic'?1:.55),abdomen:expansion*(kind==='abdominal'?1:.48),diaphragm:expansion*(kind==='thoracic'?.58:1),forced:kind==='forced'&&t>=.5};
}
function airflow(t,s={}){
 const distance=clamp(finite(s.airDistance??25),10,70),level=clamp(finite(s.airLevel??2),1,3),maxSeconds=clamp(finite(s.airSeconds??8),3,15),elapsed=t*15;
 const active=elapsed<=maxSeconds,potential=clamp(28*level*(25/distance)**1.4,0,65),deflection=active?potential:0;
 return{distance,level,maxSeconds,elapsed,deflection,active,validSeconds:potential>=15?Math.min(elapsed,maxSeconds):0,note:'任意教学传递关系，不校准为口压、肌力或肺功能'};
}
const M={clamp,finite,add,sub,mul,dot,cross,rot,lineFoot,lever,bone,capacity,gait,breath,airflow,kneeFromEnds};g.LabMath=M;if(typeof module!=='undefined')module.exports=M;
})(typeof window!=='undefined'?window:globalThis);


/* Auditable classroom reductions. This file uses one state for graphics and numbers.
 * No computed force or path length is a patient measurement. */
(function(g){'use strict';
const M=g.LabMath,{add,sub,mul,dot,cross,rot,clamp,finite,lineFoot}=M;
const fallback={heel:{F:[0,0],ankle:[-.130072087,.043846436],E:[-.174431251,.032861107],O:[-.127909871,.352248261]}};
M.lever=function(kind,t,s={}){
 t=clamp(finite(t));const load=clamp(finite(s.load??30),5,300),L=clamp(finite(s.distance??.32),.18,.42);let theta,E,R,O,ankle,alpha;
 if(kind==='head'){
  alpha=-10+30*t;theta=-alpha*Math.PI/180;
  E=rot([-.035,.027],theta);R=rot([.026,.090],theta);O=[-.049,-.100];
 }else if(kind==='heel'){
  const k=g.LabRegistration?.landmarks?.heel||fallback.heel;
  alpha=24*t;theta=-alpha*Math.PI/180;E=rot(k.E,theta);R=rot(k.ankle,theta);ankle=R;
  O=add(k.O,sub(ankle,k.ankle)); // shank and muscle origin translate with ankle.
 }else{
  alpha=20+85*t;theta=alpha*Math.PI/180;E=rot([.004,-.047],theta);R=rot([0,-L],theta);O=[.017,.285];
 }
 const F=[0,0],u=mul(sub(O,E),1/Math.hypot(...sub(O,E))),r=[0,-1],fe=lineFoot(F,E,u),fr=lineFoot(F,R,r),de=Math.hypot(...fe),dr=Math.hypot(...fr),se=cross(E,u),sr=cross(R,r);
 const effort=Math.abs(se)>1e-10?-load*sr/se:null;
 return{kind,t,theta,angle:alpha,F,E,R,O,u,r,fe,fr,ankle,effortArm:de,resistanceArm:dr,load,effort,torque:load*dr,advantage:dr>1e-10?de/dr:null,distance:L,balanceResidual:effort===null?null:effort*se+load*sr,
   note:kind==='heel'?'跖趾区为等效支承轴；趾骨保持，后足—中足转动，小腿随踝中心平移。':'骨网格与几何共用参考点；肌肉作用线为二维理想化。'};
};
M.contraction=function(mode,t,s={}){
 const u=clamp(finite(t)),delta=85,omega=clamp(finite(s.isoSpeed??60),30,180),reverse=mode==='eccentric'||mode==='isokinetic'&&s.isoDirection==='eccentric';
 const iso=mode==='isokinetic',fixed=mode==='isometric',duration=fixed?6:iso?delta/omega:3;
 const phase=fixed?.5:iso?u:(1-Math.cos(Math.PI*u))/2;
 const angle=fixed?75:reverse?105-delta*phase:20+delta*phase;
 const velocity=fixed?0:(reverse?-1:1)*(iso?omega:delta*Math.PI*Math.sin(Math.PI*u)/(2*duration));
 const q=angle*Math.PI/180,elbow=[0,0],shoulder=[0,.30],wrist=rot([0,-.265],q),hand=rot([0,-.39],q);
 const distal=rot([.004,-.047],q),proximal=[.019,.275],path=Math.hypot(...sub(distal,proximal));
 const type=fixed?'等长收缩':iso?'等速收缩':reverse?'离心收缩':'向心收缩';
 return{mode,u,t:u,time:u*duration,duration,angle,velocity,speed:omega,type,phase:fixed?'保持':reverse?'主动控制放下':'主动提起',
     lengthChange:fixed?'基本保持':reverse?'受力时拉长':'发力时缩短',work:fixed?'等长':reverse?'离心':'向心',elbow,shoulder,wrist,hand,distal,proximal,q,path,
     device:iso,relativeEffort:clamp(finite(s.effort??.6),.2,1),
     deviceTorque:iso?clamp(finite(s.effort??.6),.2,1)*(0.8+0.2*Math.sin(q)):null,
     note:'肌腹与纤维为几何示意。等速仅表示仪器约束的恒定角速度工作段，不等于恒力或恒张力。'};
};
M.gaitEvents=function(pattern='normal',side='right'){
 const d=M.gait(0,pattern),shift=d.shift,rs=d.stanceR,ls=d.stanceL;
 const events=[[0,'右初始接触'],[(shift+ls)%1,'左足离地'],[shift,'左初始接触'],[rs,'右足离地'],[1,'下次右初始接触']].sort((a,b)=>a[0]-b[0]);
 return side==='right'?events:events;
};
M.gaitPhase=function(t,pattern='normal',side='right'){
 const d=M.gait(t,pattern),l=d.legs[side],p=l.phase,st=l.stance;
 let phase=p<.10?'承重反应':p<st*.5?'支撑中期':p<st-.10?'支撑末期':p<st?'预摆动':p<st+(1-st)/3?'摆动初期':p<st+2*(1-st)/3?'摆动中期':'摆动末期';
 if(Math.abs(p)<1e-8)phase='初始接触（事件）';
 return{...d,referenceSide:side,phase,events:M.gaitEvents(pattern,side)};
};
})(typeof window!=='undefined'?window:globalThis);
/* Contact-consistent prescribed gait. Segment lengths remain constant;
 * posture is not calibrated as a population normal-value model. */
(function(g){'use strict';const M=g.LabMath,{add,sub,mul,rot,clamp}=M;
M.gait=function(t,pattern='normal'){
 t=clamp(M.finite(t));const periodic=t===1?0:t,S=1.0,shift=pattern==='short'?.42:.5,stR=pattern==='short'?.50:.6,stL=pattern==='short'?.73:.6;
 const hip=[S*t,.900-.014*Math.cos(4*Math.PI*periodic),.016*Math.sin(2*Math.PI*periodic)];const legs={};
 for(const side of ['right','left']){
  const offset=side==='right'?0:shift,st=side==='right'?stR:stL,raw=t-offset,cycle=Math.floor(raw+1e-10),p=raw-cycle,contact=p<st||Math.abs(p)<1e-8;
  const stancePos=(cycle+offset)*S+.23,z=side==='right'?.095:-.095;let x=stancePos,y=.07,pitch=0;
  if(contact){if(p<.10){pitch=(1-p/.10)*.14;const he=rot([-.055,-.070],pitch);x=stancePos-.055-he[0];y=-he[1];}
   else if(p>st-.13){pitch=-.33*(p-(st-.13))/.13;const toe=rot([.17,-.070],pitch);x=stancePos+.17-toe[0];y=-toe[1];}}
  else{const u=(p-st)/(1-st),ep=.14,sp=-.33,ps=rot([.17,-.07],sp),pe=rot([-.055,-.07],ep),e=u*u*(3-2*u);x=(1-e)*(stancePos+.17-ps[0])+e*(stancePos+S-.055-pe[0]);pitch=sp+(ep-sp)*u;
   const he=rot([-.055,-.07],pitch),to=rot([.17,-.07],pitch),clear=(pattern==='clearance'&&side==='right'?.006:.064)*Math.sin(Math.PI*u);
   // Whole foot stays above ground during swing; reduced clearance is not negative penetration.
   y=Math.max(-he[1],-to[1])+.0001+clear;
  }
  const ankle=[x,y,z];legs[side]={ankle,phase:p,stance:st,contact,pitch};
  const hx=hip[0],hz=hip[2]+(side==='right'?.07:-.07),horizontal=(x-hx)**2+(z-hz)**2;
  hip[1]=Math.min(hip[1],y+Math.sqrt(Math.max(.001,.847**2-horizontal)));
 }
 const lean=pattern==='lean'?.055*Math.max(0,Math.sin(2*Math.PI*periodic)):0,shoulder=[hip[0]-.01,hip[1]+.50,hip[2]+lean],head=[shoulder[0],shoulder[1]+.18,shoulder[2]],arms={};
 for(const side of ['right','left']){
  const l=legs[side],h=[hip[0],hip[1],hip[2]+(side==='right'?.07:-.07)];l.hip=h;l.knee=M.kneeFromEnds(h,l.ankle,.42,.43,1);
  const he=rot([-.055,-.07],l.pitch),to=rot([.17,-.07],l.pitch);l.heel=[l.ankle[0]+he[0],l.ankle[1]+he[1],l.ankle[2]];l.toe=[l.ankle[0]+to[0],l.ankle[1]+to[1],l.ankle[2]];
  const sg=side==='right'?1:-1,q=-.31*Math.sin(2*Math.PI*(periodic+.12))*sg,s=[shoulder[0],shoulder[1]-.05,shoulder[2]+sg*.16],el=[s[0]+Math.sin(q)*.27,s[1]-Math.cos(q)*.27,s[2]],w=[el[0]+Math.sin(q+.25)*.23,el[1]-Math.cos(q+.25)*.23,s[2]];arms[side]={s,elbow:el,wrist:w};
 }
 const weights=[[head,.08],[mul(add(hip,shoulder),.5),.50]];for(const side of ['right','left']){const l=legs[side],a=arms[side];weights.push([mul(add(l.hip,l.knee),.5),.10],[mul(add(l.knee,l.ankle),.5),.045],[mul(add(l.heel,l.toe),.5),.015],[mul(add(a.s,a.elbow),.5),.026],[mul(add(a.elbow,a.wrist),.5),.024]);}
 const com=weights.reduce((z,[v,w])=>add(z,mul(v,w)),[0,0,0]),R=legs.right.contact,L=legs.left.contact;
 return{t,hip,shoulder,head,arms,legs,com,support:R&&L?'双支撑':R?'右单支撑 · 左摆动':L?'左单支撑 · 右摆动':'双足离地',contact:[R,L],shift,stanceR:stR,stanceL:stL,massSum:weights.reduce((v,[p,w])=>v+w,0),stride:S,duration:1.1};
};
})(typeof window!=='undefined'?window:globalThis);


/* V8.2 dynamic teaching panels. Each diagram is synchronized with LabMath. */
(function(g){'use strict';const M=LabMath,old=Sim.render.bind(Sim),C={ink:'#183b3b',muted:'#536a69',green:'#087f76',blue:'#346a9e',red:'#c05a43',bone:'#e6dac5',line:'#d6e4df',bg:'#f7faf8'};
const esc=Sim.esc,txt=(x,y,s,size=20,col=C.ink,anchor='start')=>`<text x="${x}" y="${y}" font-size="${size}" fill="${col}" text-anchor="${anchor}">${esc(s)}</text>`;
const line=(a,b,c=C.green,w=3,dash='')=>`<line x1="${a[0]}" y1="${a[1]}" x2="${b[0]}" y2="${b[1]}" stroke="${c}" stroke-width="${w}" stroke-linecap="round" ${dash?`stroke-dasharray="${dash}"`:''}/>`;
const poly=(p,c=C.green,w=3)=>`<polyline points="${p.map(a=>a.join(',')).join(' ')}" fill="none" stroke="${c}" stroke-width="${w}" stroke-linejoin="round" stroke-linecap="round"/>`;
const dot=(p,r=5,c=C.red)=>`<circle cx="${p[0]}" cy="${p[1]}" r="${r}" fill="${c}"/>`;
const box=(x,y,w,h,fill='white')=>`<rect x="${x}" y="${y}" width="${w}" height="${h}" rx="15" fill="${fill}" stroke="${C.line}"/>`;
const arrow=(a,b,c=C.green,w=3)=>{const q=Math.atan2(b[1]-a[1],b[0]-a[0]),s=11;return line(a,b,c,w)+poly([[b[0]-s*Math.cos(q-.4),b[1]-s*Math.sin(q-.4)],b,[b[0]-s*Math.cos(q+.4),b[1]-s*Math.sin(q+.4)]],c,w);};
const wrap=(x,y,s,max=26,size=18,col=C.muted)=>String(s).match(new RegExp('.{1,'+max+'}','gu'))?.map((str,i)=>txt(x,y+i*(size+9),str,size,col)).join('')||'';
const svg=(a,body,h=590)=>`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1120 ${h}" role="img" aria-label="${esc(a.title)}"><rect width="1120" height="${h}" fill="${C.bg}"/><g font-family="system-ui,Microsoft YaHei,sans-serif">${body}</g></svg>`;
const metric=(x,y,label,value,col=C.ink)=>txt(x,y,label,17,C.muted)+txt(x,y+33,value,26,col);
const f=(v,n=1)=>Number.isFinite(v)?v.toFixed(n):'—';
function anatomyImage(name,x,y,w,h){const data=g.LabImages?.[name];return data?`<image href="${data}" x="${x}" y="${y}" width="${w}" height="${h}" preserveAspectRatio="none"/>`:'';}
function lever(a,t,s){const d=M.lever(a.motion,t,s),head=a.motion==='head',heel=a.motion==='heel';const origin=head?[300,370]:heel?[520,465]:[285,315],scale=head?1000:heel?840:770,to=p=>[origin[0]+p[0]*scale,origin[1]-p[1]*scale];
 let b=txt(30,42,a.title,27)+txt(30,72,'拖动动作 → 同步比较垂直力臂；绿色为动力，橙色为阻力。',18,C.muted);
 if(head){b+=anatomyImage('neck',235,346,96,191)+`<g transform="translate(${origin}) rotate(${-d.theta*180/Math.PI})">${anatomyImage('head',-103,-244,240,272)}</g>`;}
 else if(heel){b+=line([160,469],[558,469],C.line,4)+line(to(d.R),to([-.15,.39]),C.bone,24)+`<g transform="translate(${origin}) rotate(${-d.theta*180/Math.PI})">${anatomyImage('foot',-227,-73,248,73)}</g>`;}
 else{b+=anatomyImage('humerus',254,87,64,238)+`<g transform="translate(${origin}) rotate(${-d.theta*180/Math.PI})">${anatomyImage('forearm',-33,-7,66,290)}</g>`;}
 const F=to(d.F),E=to(d.E),R=to(d.R),O=to(d.O),FE=to(d.fe),FR=to(d.fr),sE=to(M.add(d.E,M.mul(d.u,.065)));
 b+=line(E,O,'#d79b89',8)+line(to(M.add(d.E,M.mul(d.u,-.08))),to(M.add(d.E,M.mul(d.u,.17))),C.green,1.5,'6 5');
 b+=line([R[0],Math.max(96,R[1]-65)],[R[0],Math.min(540,R[1]+85)],C.red,1.5,'6 5');
 b+=arrow(E,sE,C.green,4)+arrow(R,[R[0],R[1]+60],C.red,4)+line(F,FE,C.green,4,'7 5')+line(F,FR,C.red,4,'7 5')+dot(F,7,C.blue)+dot(E,5,C.green)+dot(R,5,C.red)+dot(FE,4,C.green)+dot(FR,4,C.red);
 b+=txt(F[0]-8,F[1]+25,'F',18,C.blue)+txt(E[0]-25,E[1]-12,'E',18,C.green)+txt(R[0]+12,R[1]-10,'R',18,C.red);
 b+=box(636,104,449,371)+metric(662,137,'动力臂 dE',f(d.effortArm*1000)+' mm',C.green)+metric(880,137,'阻力臂 dR',f(d.resistanceArm*1000)+' mm',C.red);
 b+=metric(662,221,'理想机械优势 dE / dR',f(d.advantage,2),C.blue)+metric(880,221,'平衡所需动力',f(d.effort)+' N');
 const samples=Array.from({length:61},(_,i)=>M.lever(a.motion,i/60,s)),max=Math.max(...samples.flatMap(q=>[q.effortArm,q.resistanceArm]))*1.1;
 const pt=(x,v)=>[670+370*x,431-102*v/max];b+=line([670,432],[1042,432],C.line,2)+poly(samples.map((q,i)=>pt(i/60,q.effortArm)),C.green,3)+poly(samples.map((q,i)=>pt(i/60,q.resistanceArm)),C.red,3)+dot(pt(t,d.effortArm),6,C.green)+dot(pt(t,d.resistanceArm),6,C.red)+line([670+370*t,319],[670+370*t,435],C.blue,1,'4 5')+txt(667,457,'动作起始',15,C.muted)+txt(1041,457,'动作终末',15,C.muted,'end');
 b+=txt(637,515,head?'支点取寰枕区教学参考，不是整段颈椎固定单轴。':heel?'R是传入足的阻力合力；不把它直接等同体重。':'前臂随屈肘旋转；手中阻力保持竖直向下。',17,C.muted);
 b+=txt(30,553,`动作示意角 ${f(d.angle)}°　·　骨网格投影；测力点与肌肉路径是理想化假设。`,17,C.muted);
 return svg(a,b,580);
}
function bone(a,t,s){const d=M.bone(a.motion,t,s),kind=a.motion,fmtx=x=>Math.abs(x)<.005?x.toExponential(2):x.toFixed(3);
 let b=txt(30,42,a.title,27)+txt(30,72,'理想试件：从未加载虚线，看到加载后的网格变化。不是人体骨材料测量。',18,C.muted);
 const cx=280,top=154,bottom=423,H=bottom-top,W=90;
 b+=`<rect x="${cx-W/2}" y="${top}" width="${W}" height="${H}" fill="none" stroke="#9caeaa" stroke-dasharray="6 5"/>`;
 let map=(x,y)=>[x,y],amp=1;
 if(kind==='compression'||kind==='tension'){amp=400;const sign=kind==='compression'?-1:1,delta=d.displacement*amp;map=(x,y)=>[cx+(x-cx)*(1-sign*.3*delta/H),top+H/2+(y-top-H/2)*(1+sign*delta/H)];}
 if(kind==='shear'){amp=650;map=(x,y)=>[x+d.displacement*amp*(bottom-y)/H,y];}
 if(kind==='bending'){amp=155;map=(x,y)=>{const u=(bottom-y)/H;return[x+d.displacement*amp*u*u*(3-u)/2,y];};}
 b+=`<polygon points="${[[cx-W/2,top],[cx+W/2,top],[cx+W/2,bottom],[cx-W/2,bottom]].map(p=>map(...p).join(',')).join(' ')}" fill="#eee3d1" stroke="none"/>`;
 for(let i=0;i<=6;i++)b+=poly(Array.from({length:21},(_,j)=>map(cx-W/2+W*i/6,top+H*j/20)),i===3?C.blue:'#a08765',i===3?2:1.5);
 for(let j=0;j<=10;j++)b+=poly(Array.from({length:9},(_,i)=>map(cx-W/2+W*i/8,top+H*j/10)),'#bfae92',1.3);
 if(kind==='torsion'){
  amp=20;b+=`<rect x="${cx-W/2}" y="${top}" width="${W}" height="${H}" rx="20" fill="#eee3d1" stroke="#a08765"/>`;
  for(let i=0;i<10;i++){const p=[];for(let j=0;j<=30;j++){const u=j/30,q=i*Math.PI/5+d.x*amp*(1-u);p.push([cx+W/2*Math.cos(q),top+H*u]);}b+=poly(p,'#a08765',1.5);}
  b+=`<ellipse cx="${cx}" cy="${top}" rx="45" ry="13" fill="#e1d0b5" stroke="#a08765"/>`+line([cx,top],[cx+44*Math.cos(d.x*amp),top-13*Math.sin(d.x*amp)],C.red,4)+arrow([cx+67,top+10],[cx+50,top-16],C.red,4);
 }else if(kind==='tension')b+=arrow(map(cx,top),[cx,106],C.red)+arrow(map(cx,bottom),[cx,473],C.red);
 else if(kind==='compression')b+=arrow([cx,104],map(cx,top),C.red)+arrow([cx,471],map(cx,bottom),C.red);
 else if(kind==='shear')b+=arrow([190,129],[330+40*d.level,129],C.red)+arrow([330,449],[199,449],C.red);
 else b+=arrow([153,top],[map(cx,top)[0]-10,top],C.red)+line([207,bottom+8],[358,bottom+8],C.ink,7)+txt(380,227,'外纤维拉压相反',17,C.muted);
 b+=txt(68,514,`显示放大 ×${amp}　｜　蓝线为参考中线`,18,C.muted);
 b+=box(605,101,481,408)+txt(633,133,d.yLabel,20)+line([660,412],[1037,412],C.ink,2)+line([660,412],[660,172],C.ink,2);
 const baseline=M.bone(kind,1,{...s,material:'A'}),xmax=baseline.xMax*1.15,ymax=baseline.yMax*1.15,xy=(x,y)=>[660+355*x/xmax,412-225*y/ymax];
 for(let j=1;j<5;j++)b+=line([660,412-45*j],[1030,412-45*j],C.line,1);
 for(const mat of ['A','B']){const arr=Array.from({length:31},(_,i)=>M.bone(kind,i/30,{...s,material:mat,loadCycle:false}));b+=poly(arr.map(q=>xy(q.x,q.y)),mat===s.material?C.green:mat==='A'?C.green:'#a4b6b2',mat===(s.material||'A')?4:2);}
 b+=dot(xy(d.x,d.y),7,C.red)+txt(1018,447,d.xLabel,18,C.ink,'end')+txt(660,439,'0',15,C.muted)+txt(662,474,`当前：${fmtx(d.x)}  →  ${f(d.y,2)}`,20,C.red);
 b+=wrap(631,540,d.info,30,17,C.muted)+txt(30,575,'深色线：所选材料；浅色线：另一材料。只比较弹性刚度，不据斜率判断强度或骨折。',17,C.muted);
 return svg(a,b,612);
}
function personSquat(depth=0,up=0,center=272,ground=489,scale=230,highlight=false,reach=false){
 const A=[0,.07],shin=(6+23*depth)*Math.PI/180,thigh=(5+67*depth)*Math.PI/180;
 const K=[A[0]+.43*Math.sin(shin),A[1]+.43*Math.cos(shin)],H=[K[0]-.42*Math.sin(thigh),K[1]+.42*Math.cos(thigh)],lean=(5+32*depth)*Math.PI/180,S=[H[0]+.5*Math.sin(lean),H[1]+.5*Math.cos(lean)],N=[S[0],S[1]+.18];
 const xy=p=>[center+p[0]*scale,ground-(p[1]+up)*scale];
 let b=line([center-133,ground+3],[center+153,ground+3],C.line,3);
 b+=poly([H,K,A].map(xy),'#b4c4be',23)+poly([H,K,A].map(p=>[xy(p)[0]+12,xy(p)[1]]),C.ink,17)+line(xy(A),xy([.18,.015]),C.ink,15)+line(xy(H),xy(S),C.green,29)+dot(xy(N),22,'#c89975');
 const EL=reach?[S[0]+.09,1.79-depth*.25]:[S[0]+.22,S[1]-.20],WR=reach?[S[0]+.03,2.15-depth*.30]:[EL[0]+.18,EL[1]-.15];b+=poly([S,EL,WR].map(xy),'#b4c4be',12);if(reach)b+=dot(xy(WR),6,C.red);
 if(highlight){b+=line([xy(H)[0]+12,xy(H)[1]+18],[xy(K)[0]+12,xy(K)[1]-15],C.red,10)+line([xy(K)[0]-12,xy(K)[1]+20],[xy(A)[0]-12,xy(A)[1]-20],C.red,9);}
 return {svg:b,H,K,A,S,N,xy};
}
function capacity(a,t,s){const d=M.capacity(a.motion,t,s);let b=txt(30,42,a.title,27)+txt(30,72,'给定任务、动作质量和指标一起观察。所有结果是预设案例或数学示例。',18,C.muted);
 if(a.motion==='rm'){
  const F=[278,277],ang=d.angle*Math.PI/180,R=[F[0]+137*Math.sin(ang),F[1]+137*Math.cos(ang)];
  b+=dot([247,139],22,'#c89975')+line([248,174],[248,312],C.green,36)+poly([[248,311],[309,344],[309,464]],C.ink,18)+line([309,464],[353,464],C.ink,15)+line([265,190],F,C.bone,22)+line(F,R,C.bone,18)+dot(F,7,C.blue)+box(R[0]-23,R[1]-12,46,39,'#d9dfdc')+line([R[0]-46,R[1]+9],[R[0]+46,R[1]+9],C.ink,8)+txt(R[0],R[1]-22,`${d.mass} kg`,22,C.ink,'middle');
  b+=box(593,106,491,399)+metric(621,139,'本例举起负荷',`${d.mass} kg`)+txt(621,220,d.label,22,d.successful?C.green:C.red)+txt(621,265,'同一规定动作的试次',19,C.muted);
  ['8 kg：完成一次','10 kg：完成一次','12 kg：完成一次','14 kg：未完成全程'].forEach((v,i)=>b+=txt(627,303+i*35,v,21));
  b+=txt(621,474,'本例已知1RM：12 kg（不是学生实测）',19,C.green);
 }else if(a.motion==='jump'){
  const actor=personSquat(d.squat,d.height,269,506,143,false,true);b+=actor.svg;
  // Visual reach levels use the same constant arm reach; not a claim of force/power measurement.
  const ry=506-(d.standingReach+d.heightMax)*143;b+=line([483,105],[483,506],C.line,4)+line([456,ry],[505,ry],C.red,4)+txt(466,ry+5,'最高触及',16,C.red,'end');
  b+=box(593,106,491,399)+metric(620,140,'当前阶段',d.phase)+metric(620,229,'腾空时间（等高起落假设）',f(d.flight,3)+' s')+metric(620,317,'由起跳速度计算的质心上升',f(d.heightMax*100,1)+' cm');
  b+=txt(620,406,`站立摸高 ${f(d.standingReach*100,0)} cm`,20)+txt(620,439,`跳起摸高 ${f(d.jumpReach*100,1)} cm（同伸手假设）`,20);
 }else if(a.motion==='plank'){
  const hip=[341,302+(d.broken?45:0)],shoulder=[150,265],toe=[533,371],elbow=[153,371];b+=line([84,396],[555,396],C.line,6)+dot([91,248],21,'#c89975')+line(shoulder,hip,C.green,28)+poly([hip,[441,336],toe],C.ink,19)+line(shoulder,elbow,'#95aea4',15)+line(elbow,[219,375],C.ink,13);
  b+=box(608,121,463,354)+metric(634,155,'案例经过时间',f(d.elapsed,1)+' s')+metric(634,249,'符合姿势要求的累计时间',f(d.validSeconds,1)+' s',d.broken?C.red:C.green)+wrap(634,353,d.label,23,22)+wrap(634,407,'用姿势标准决定计时停止，不把忍耐不适当作成功。',24,17);
 }else{
  b+=personSquat(d.squat,0,270,504,235).svg;b+=box(601,121,472,361)+metric(630,157,'完整重复数',String(d.total))+metric(846,157,'有效次数',String(d.valid),C.green)+metric(630,257,'规定节律', '3秒 / 次')+wrap(630,349,d.label,23,20)+wrap(630,414,'演示案例只展示规则；不按鼠标拖动计算使用者肌耐力。',24,17);
 }
 b+=txt(30,561,a.motion==='jump'?'摸高受身高、伸手方式与协调影响；腾空估高也受起落姿态影响，不能直接称为功率。':a.teach,17,C.muted);return svg(a,b,595);
}
function ssc(a,t,s){const d=M.capacity('ssc',t,s);let b=txt(30,42,'牵拉—缩短周期 · 从下蹲到蹬伸',27)+txt(30,73,'下蹲制动、转换、蹬伸和落地分开；股四头肌与小腿三头肌按功能阶段提示。',18,C.muted);b+=personSquat(d.squat,d.height,277,508,177,true).svg;
 b+=box(586,108,498,368)+metric(615,142,'当前阶段',d.phase,C.green)+txt(615,233,'股四头肌',22)+txt(615,267,d.quad,20,C.red)+txt(615,323,'小腿三头肌',22)+txt(615,357,d.calf,20,C.red)+wrap(615,408,'颜色只标记本阶段观察对象，不表示纤维长度、肌电幅值或肌力大小。',26,17);
 b+=txt(31,557,'转换效率还与预激活、力矩形成和动作协调有关；不把跳跃增益全部归为弹性储能。',17,C.muted);return svg(a,b);
}
function breathing(a,t){const d=M.breath(t,a.motion);let b=txt(30,42,a.title,27)+txt(30,72,'同一个呼吸周期中，对照胸廓扩张、膈肌下降和腹壁外移。',18,C.muted);
 function body(cx,cy,ch,ab,dia,label){let z=txt(cx,114,label,22,C.ink,'middle');const w=85+ch*20,front=32+ab*30;
  z+=`<path d="M${cx-40} 147 Q${cx-w} 161 ${cx-w} 304 Q${cx-w+8} 341 ${cx-56} 385 L${cx-52} 442 Q${cx+front} 465 ${cx+68} 425 L${cx+72} 347 Q${cx+w} 322 ${cx+w} 251 Q${cx+w} 174 ${cx+40} 147Z" fill="#e8f2ee" stroke="#84a397" stroke-width="2.5"/>`;
  for(let i=0;i<5;i++){const yy=187+i*30;z+=`<path d="M${cx} ${yy+21} Q${cx-w} ${yy+20-ch*10} ${cx-w+9} ${yy-4} M${cx} ${yy+21} Q${cx+w} ${yy+20-ch*10} ${cx+w-9} ${yy-4}" fill="none" stroke="#c7b69a" stroke-width="9"/>`;}
  z+=line([cx,171],[cx,308],C.bone,9);
  // More inferior/flat diaphragm during inspiration (screen y increases downward).
  const dome=309+32*dia;z+=`<path d="M${cx-w+15} 356 Q${cx} ${dome-26} ${cx+w-15} 356" fill="#edc9bb" stroke="${C.red}" stroke-width="5"/>`+arrow([cx,296],[cx,308+30*dia],C.red,3);
  z+=`<path d="M${cx+70} 359 Q${cx+80+front} 400 ${cx+62} 442" fill="none" stroke="${d.forced?C.red:C.green}" stroke-width="6"/>`;
  z+=arrow([cx-w+6,238],[cx-w-23*ch,238],C.green)+arrow([cx+w-6,238],[cx+w+23*ch,238],C.green)+txt(cx,486,`${d.phase} · 胸廓与腹壁都可参与`,17,C.muted,'middle');return z;
 }
 if(a.motion==='quiet'||a.motion==='forced'){
  b+=body(282,0,d.chest,d.abdomen,d.diaphragm,'胸腹联合示意')+box(617,123,462,330)+metric(646,158,'当前相位',d.phase,C.green)+wrap(646,248,d.forced?'腹肌等主动参与呼气；不是安静呼气时都主动压缩。':t<.5?'膈肌收缩下降，胸廓扩张。':'吸气肌放松后，以弹性回缩为主。',23,20)+wrap(646,369,'此图为功能结构示意，不显示真实肌纤维或精确器官形变。',24,17);
 }else{
  const e=d.expansion;b+=body(275,0,e,.48*e,.58*e,'胸廓活动较明显')+body(786,0,.55*e,e,e,'腹壁活动较明显');
 }
 b+=txt(30,556,'“胸式 / 腹式”表示贡献相对突出，不是把胸部与腹部画成完全互斥的两套呼吸。',17,C.muted);return svg(a,b);
}
function airflow(a,t,s){const d=M.airflow(t,s),mode=a.motion;let b=txt(30,42,a.title,27)+txt(30,72,mode==='distance'?'固定输出条件，只改变目标距离，比较纸张偏转。':'固定距离与偏转要求，观察持续达到要求的时间。',18,C.muted);
 const x=338+d.distance*3.8,y=220,ang=-d.deflection;
 b+=`<path d="M111 201 Q132 182 157 200 L183 220 L157 239 Q126 251 111 233" fill="#dab28f"/>`;
 if(d.active)for(let i=0;i<3;i++)b+=arrow([198,202+i*22],[Math.min(x-15,480),202+i*22],C.green,3);
 b+=line([x,129],[x,174],C.ink,3)+`<g transform="rotate(${ang},${x},174)"><rect x="${x}" y="174" width="77" height="172" rx="3" fill="white" stroke="#91a99e" stroke-width="2.5"/></g>`;
 b+=line([189,390],[x,390],C.line,3)+txt((189+x)/2,420,`${d.distance} cm`,23,C.blue,'middle')+txt(137,463,d.active?'正在输出':'本次输出结束，纸张回落',21,d.active?C.green:C.red);
 b+=box(750,129,337,360)+metric(774,163,'纸张偏转（示意）',f(d.deflection,0)+'°',C.green)+metric(774,255,'案例经过时间',f(d.elapsed,1)+' s')+metric(774,346,'达到本题阈值的持续时间',f(d.validSeconds,1)+' s');
 b+=txt(30,551,'距离和时间是课堂任务指标；不直接换算呼吸肌力、MIP / MEP或标准肌耐力。仅使用虚拟火焰。',17,C.muted);return svg(a,b);
}
function gait(a,t,s){const pattern=['clearance','lean','short'].includes(a.motion)?a.motion:'normal',d=M.gait(t,pattern),scale=235,to=p=>[250+(p[0]-(d.stride||.9)*t)*scale,497-p[1]*scale],sideColor={right:C.green,left:C.blue};
 let b=txt(30,42,a.title,27)+txt(30,72,'参考周期：右侧初始接触 → 下一次右侧初始接触。左、右支撑同时显示。',18,C.muted);
 for(let i=-2;i<=4;i++){const x=250+(i*.3-(d.stride||.9)*t)*scale;if(x>25&&x<537)b+=line([x,503],[x+27,503],C.line,3);}
 for(const side of ['left','right']){const l=d.legs[side],col=sideColor[side];b+=poly([l.hip,l.knee,l.ankle].map(to),col,side==='right'?18:13)+line(to(l.heel),to(l.toe),col,12);if(l.contact)b+=line([to(l.heel)[0],506],[to(l.toe)[0],506],col,5);}
 b+=line(to(d.hip),to(d.shoulder),C.ink,30)+dot(to(d.head),22,'#c89975');
 for(const side of ['left','right']){const a=d.arms[side];b+=poly([a.s,a.elbow,a.wrist].map(to),sideColor[side],11);}
 b+=dot(to(d.com),9,C.red)+txt(31,545,d.support,24,C.green);
 b+=box(590,106,501,170)+txt(617,140,'左右足接触',21)+txt(617,181,'右足',18,C.green)+txt(617,229,'左足',18,C.blue);
 for(const [side,y]of[['right',176],['left',224]]){
  const pts=[];for(let i=0;i<=100;i++){const on=M.gait(i/100,pattern).legs[side].contact;if(on)b+=line([674+i*3.6,y],[674+(i+1)*3.6,y],sideColor[side],12);}b+=dot([674+t*360,y],7,C.red);
 }
 b+=txt(674,260,'0%',15,C.muted)+txt(1034,260,'100%',15,C.muted,'end');
 const vals=Array.from({length:81},(_,i)=>M.gait(i/80,pattern));
 if(a.motion==='timing'){
  b+=box(590,301,501,228)+txt(617,332,'相关肌群的功能任务（不是肌电）',21);
  const bands=[['胫骨前肌',[[0,.10],[.60,1]]],['小腿三头肌',[[.12,.60]]],['股四头肌',[[0,.20]]],['腘绳肌',[[0,.12],[.85,1]]]];
  bands.forEach(([label,segs],i)=>{const y=364+i*40;b+=txt(617,y+6,label,17);for(const [l,r]of segs)b+=line([757+270*l,y],[757+270*r,y],C.green,12);});b+=line([757+270*t,345],[757+270*t,511],C.red,2);
 }else{
  for(const [idx,y,name,col]of[[1,318,'模型质心 · 上下',C.red],[2,428,'模型质心 · 左右',C.blue]]){
   const arr=vals.map(x=>x.com[idx]),min=Math.min(...arr),max=Math.max(...arr),xy=(i,v)=>[620+i*5.4,y+70-(v-min)/(max-min||1)*58];
   b+=txt(615,y-7,`${name}  幅差 ${f((max-min)*1000,1)} mm`,18,col)+poly(arr.map((v,i)=>xy(i,v)),col,3)+dot([620+t*432,y+70-(d.com[idx]-min)/(max-min||1)*58],6,col);
  }
 }
 b+=txt(30,580,'合成步态：足接触与几何同步，COM使用本模型明确假定的分段质量；不作为正常值或病例诊断。',16,C.muted);
 return svg(a,b,610);
}
function cardio(a,t,s){const safety=s.safety||{},moving=a.motion==='walk'&&safety.state==='running',d=M.gait(moving?(t*4)%1:0),baseline=a.motion==='rest',rest=a.motion==='recovery';const hr=baseline?74:rest?Math.round(74+48*Math.exp(-3*t)):Math.round(74+50*(1-Math.exp(-3*t))),rr=baseline?14:rest?Math.round(14+9*Math.exp(-3*t)):Math.round(14+12*(1-Math.exp(-3*t)));
 let b=txt(30,42,a.title,27)+txt(30,72,'运动反应样例 + 风险应对。安全状态会控制播放、皮带和采集。',18,C.muted);
 const sc=204,xy=p=>[257+(p[0]-(d.stride||.9)*d.t)*sc,469-p[1]*sc];
 b+=`<rect x="94" y="477" width="409" height="32" rx="12" fill="#627c76"/>`;
 for(let i=0;i<10;i++){const x=103+(i*39+(moving?t*470:0))%383;b+=line([x,483],[x+14,499],'#bfd3c9',3);}
 b+=poly([[464,476],[464,227],[494,211]],C.green,7)+box(464,178,91,47)+dot([487,202],9,C.red)+txt(482,168,'急停',17,C.red);
 for(const side of ['left','right']){const l=d.legs[side];b+=poly([l.hip,l.knee,l.ankle].map(xy),side==='right'?C.ink:'#a2b9ae',15)+line(xy(l.heel),xy(l.toe),C.ink,12);}
 b+=line(xy(d.hip),xy(d.shoulder),C.green,30)+dot(xy(d.head),23,'#c89975');for(const sd of ['right','left']){const ar=d.arms[sd];b+=poly([ar.s,ar.elbow,ar.wrist].map(xy),'#9cb3a9',11);}
 b+=line(xy(d.shoulder),[477,211],C.red,2,'5 5')+txt(366,277,'安全绳示意',16,C.red);
 b+=box(628,112,458,405)+txt(651,151,safety.state==='emergency'?'运动已终止 · 风险应对':moving?'皮带运行':'皮带停止',24,C.green)+metric(653,202,'合成心率',`${hr} 次/分`,C.red)+metric(875,202,'合成呼吸',`${rr} 次/分`)+txt(653,294,`给定血压样例：${baseline?'118/76':rest?'124/76':'142/78'} mmHg`,21)+txt(653,330,'血压不是连续从动画测得。',17,C.muted);
 b+=wrap(653,383,safety.message||'开始前核对设备、安全装置与旁站协助。',23,20)+txt(30,555,'急停/不适会锁定本次运动；恢复需新开试次。此处不评估实际CPR技能或给药。',17,C.muted);
 return svg(a,b);
}
function concept(a){return `<div class="readcard"><h2>${esc(a.title)}</h2><p>${esc(a.teach)}</p><div>${a.observe.map(x=>`<p>— ${esc(x)}</p>`).join('')}</div><p class="subtle">${esc(a.pitfall)}</p></div>`;}
const New={lever,bone,capacity,ssc,breathing,airflow,gait,cardio,concept};
Sim.render=function(a,t,s={}){if(a.kind==='lever')return lever(a,t,s);if(a.kind==='bone')return bone(a,t,s);if(a.kind==='capacity')return capacity(a,t,s);if(a.kind==='ssc')return ssc(a,t,s);if(a.kind==='breath')return breathing(a,t,s);if(a.kind==='airflow')return airflow(a,t,s);if(a.kind==='gait')return gait(a,t,s);if(a.kind==='cardio')return cardio(a,t,s);if(a.kind==='concept')return concept(a);return old(a,t,s);};
g.LabDiagrams=New;
})(window);


/* V8.4 unified graphic registration and muscle-focused observation.
 * Anatomy pixels are projections of the supplied meshes. Muscle shape is an
 * explicitly schematic overlay, never presented as measured fibres or force. */
(function(g){'use strict';
const M=LabMath,esc=Sim.esc,previous=Sim.render.bind(Sim);
const C={ink:'#193e3d',muted:'#607773',green:'#087f76',red:'#ba5b45',blue:'#426f99',line:'#d5e4de',bg:'#f7faf8',bone:'#ddd3be'};
const txt=(x,y,s,size=20,col=C.ink,anchor='start')=>`<text x="${x}" y="${y}" font-size="${size}" fill="${col}" text-anchor="${anchor}">${esc(s)}</text>`;
const line=(a,b,col=C.green,w=3,dash='')=>`<line x1="${a[0]}" y1="${a[1]}" x2="${b[0]}" y2="${b[1]}" stroke="${col}" stroke-width="${w}" stroke-linecap="round" ${dash?`stroke-dasharray="${dash}"`:''}/>`;
const poly=(ps,col=C.green,w=3)=>`<polyline points="${ps.map(p=>p.join(',')).join(' ')}" fill="none" stroke="${col}" stroke-width="${w}" stroke-linecap="round" stroke-linejoin="round"/>`;
const dot=(p,col=C.green,r=5)=>`<circle cx="${p[0]}" cy="${p[1]}" r="${r}" fill="${col}"/>`;
const arrow=(p,q,c=C.green,w=3)=>{const a=Math.atan2(q[1]-p[1],q[0]-p[0]);return line(p,q,c,w)+poly([[q[0]-11*Math.cos(a-.4),q[1]-11*Math.sin(a-.4)],q,[q[0]-11*Math.cos(a+.4),q[1]-11*Math.sin(a+.4)]],c,w);};
const box=(x,y,w,h)=>`<rect x="${x}" y="${y}" width="${w}" height="${h}" rx="14" fill="white" stroke="${C.line}"/>`;
const wrap=(x,y,s,max=23,size=18,col=C.muted)=>(String(s).match(new RegExp('.{1,'+max+'}','gu'))||[]).map((v,i)=>txt(x,y+i*(size+10),v,size,col)).join('');
function svg(title,body,w=1120,h=610){return `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${w} ${h}" role="img" aria-label="${esc(title)}"><rect width="${w}" height="${h}" fill="${C.bg}"/><g font-family="system-ui,Microsoft YaHei,sans-serif">${body}</g></svg>`;}
function registered(name,scale,alpha=1){const r=g.LabRegistration?.[name],src=g.LabImages?.[name];if(!r||!src)return'';return `<image href="${src}" x="${r.x*scale}" y="${-r.y*scale}" width="${r.width*scale}" height="${r.height*scale}" opacity="${alpha}" preserveAspectRatio="none" data-registered="${name}"/>`;}
function rightAngle(F,H,U,color){const n=M.sub(F,H),nlen=Math.hypot(...n);if(nlen<1e-8)return'';const u=M.mul(U,7),v=M.mul(n,7/nlen);return poly([M.add(H,u),M.add(M.add(H,u),v),M.add(H,v)],color,1.6);}
function lever(a,t,s){
 const d=M.lever(a.motion,t,s),heel=a.motion==='heel',head=a.motion==='head';
 const o=head?[310,336]:heel?[478,504]:[240,309],sc=head?1450:heel?1370:690,P=p=>[o[0]+sc*p[0],o[1]-sc*p[1]];
 let b=txt(28,42,a.title,27)+txt(28,74,'先动关节，再比较垂直力臂；骨骼、作用点和数值使用同一坐标。',18,C.muted);
 let art='';
 if(heel){
  art+=line([140,529],[589,529],C.line,4)+txt(35,557,'足底接触参考 · 未显示软组织',16,C.muted);
  art+=`<g transform="translate(${o})" data-part="toes-fixed">${registered('heel_toes',sc)}</g>`;
  art+=`<g transform="translate(${o}) rotate(${-d.theta*180/Math.PI})" data-part="foot-moving">${registered('heel_tarsals',sc)}</g>`;
  art+=`<g transform="translate(${P(d.ankle)})" data-part="shank-following">${registered('heel_shank',sc)}</g>`;
  art+=txt(36,124,'趾骨保持参考',18,C.blue)+txt(36,153,'跟部抬起，踝中心随骨段移动',18,C.muted);
  art+=line(P([-.19,-.022]),P([.072,-.022]),C.line,1);
 }else if(head){
  art+=`<g transform="translate(${o})">${registered('neck_aligned',sc)}</g>`;
  art+=`<g transform="translate(${o}) rotate(${-d.theta*180/Math.PI})">${registered('head_aligned',sc)}</g>`;
 }else{
  art+=`<g transform="translate(${o})">${registered('elbow_humerus',sc)}</g>`;
  art+=`<g transform="translate(${o}) rotate(${-d.theta*180/Math.PI})">${registered('elbow_forearm',sc)}</g>`;
 }
 const F=P(d.F),E=P(d.E),R=P(d.R),O=P(d.O),FE=P(d.fe),FR=P(d.fr),u=[d.u[0],-d.u[1]];
 art+=line(E,O,'#c99381',5)+line(P(M.add(d.E,M.mul(d.u,-.04))),P(M.add(d.E,M.mul(d.u,.13))),C.green,1.4,'5 5');
 art+=line([R[0],R[1]-90],[R[0],R[1]+80],C.red,1.4,'5 5');
 art+=arrow(E,M.add(E,M.mul(u,67)),C.green,4)+arrow(R,[R[0],R[1]+65],C.red,4);
 art+=line(F,FE,C.green,3,'6 4')+line(F,FR,C.red,3,'6 4')+rightAngle(F,FE,u,C.green)+rightAngle(F,FR,[0,1],C.red);
 art+=dot(F,C.blue,7)+dot(E,C.green,5)+dot(R,C.red,5)+txt(F[0]+11,F[1]+25,'F',19,C.blue)+txt(E[0]-22,E[1]-12,'E',19,C.green)+txt(R[0]+14,R[1]-12,'R',19,C.red);
 b+=`<defs><clipPath id="leverSceneClip"><rect x="16" y="92" width="605" height="468"/></clipPath></defs><g data-scene="lever" clip-path="url(#leverSceneClip)">${art}</g>`;
 b+=box(650,111,440,407)+txt(676,146,'同一姿态的力臂与平衡',23);
 b+=txt(677,191,'动力臂 dE',18,C.green)+txt(901,191,'阻力臂 dR',18,C.red)+txt(677,225,(d.effortArm*1000).toFixed(1)+' mm',28,C.green)+txt(901,225,(d.resistanceArm*1000).toFixed(1)+' mm',28,C.red);
 b+=txt(677,271,'机械优势 dE / dR',17,C.muted)+txt(905,271,'理想平衡动力',17,C.muted)+txt(677,306,d.advantage.toFixed(2),28,C.blue)+txt(905,306,d.effort.toFixed(1)+' N',28);
 const data=Array.from({length:81},(_,i)=>M.lever(a.motion,i/80,s)),max=Math.max(...data.flatMap(q=>[q.effortArm,q.resistanceArm]))*1.1,at=(i,v)=>[680+i*4.75,470-122*v/max];
 b+=line([680,472],[1060,472],C.line,2)+poly(data.map((q,i)=>at(i,q.effortArm)),C.green,3)+poly(data.map((q,i)=>at(i,q.resistanceArm)),C.red,3)+dot(at(t*80,d.effortArm),C.green,5)+dot(at(t*80,d.resistanceArm),C.red,5)+line([680+380*t,337],[680+380*t,476],C.blue,1,'4 4');
 b+=txt(680,499,'动作开始',16,C.muted)+txt(1060,499,'动作终末',16,C.muted,'end');
 b+=wrap(656,550,heel?'R为小腿传给足的阻力；F为跖趾区等效支承轴，不把R直接当作体重。':'反向拖动仅回看姿态。虚线是垂直力臂，不是骨长。',26,17);
 b+=txt(30,590,`教学角 ${d.angle.toFixed(1)}° · 准静态二维模型；不是人体实测肌力臂。`,17,C.muted);
 return svg(a.title,b,1120,628);
}
function muscleRibbon(P,Q,width,col,id){
 const dx=Q[0]-P[0],dy=Q[1]-P[1],len=Math.hypot(dx,dy),nx=-dy/len,ny=dx/len;
 const pt=(t,n)=>[P[0]+dx*t+nx*n,P[1]+dy*t+ny*n];
 const l1=pt(.28,width),l2=pt(.72,width*.86),r1=pt(.72,-width*.86),r2=pt(.28,-width);
 let b=`<path d="M${P} C${l1} ${l2} ${Q} C${r1} ${r2} ${P} Z" fill="${col}" opacity=".90" stroke="#713b32" stroke-width="1" data-muscle="${id}"/>`;
 for(let k=-2;k<=2;k++){const off=k*width*.25;b+=`<path d="M${pt(.10,0)} Q${pt(.50,off)} ${pt(.91,0)}" fill="none" stroke="#fff4e6" opacity=".48" stroke-width="1.3"/>`;}
 return b;
}
function contraction(a,t,s){
 const d=M.contraction(a.contractionMode,t,s),o=[222,337],sc=590,P=p=>[o[0]+sc*p[0],o[1]-sc*p[1]],flip=s.muscleSide==='left';
 let body=`<defs><linearGradient id="muscleFade" x1="0" y1="0" x2="1" y2="0"><stop stop-color="#b65644"/><stop offset=".5" stop-color="#da8c73"/><stop offset="1" stop-color="#a95244"/></linearGradient></defs>`;
 body+=txt(23,37,flip?'左侧方向示意（镜像）':'右肘侧面 · 肌肉为观察主体',23)+txt(23,67,'骨为定位；肌腹/肌腱几何示意，不给纤维实测值。',17,C.muted);
 body+=`<g transform="${flip?'translate(555 0) scale(-1 1)':''}">`;
 body+=`<g transform="translate(${o})" opacity=".64">${registered('elbow_humerus',sc)}</g><g transform="translate(${o}) rotate(${-d.angle})" opacity=".62">${registered('elbow_forearm',sc)}</g>`;
 const B0=P([.021,.287]),B1=P(M.add(d.proximal,[0,-.040])),B2=P(M.add(d.distal,[0,.027])),D=P(d.distal);
 const T0=P([-.040,.286]),T1=P([-.036,.215]),TD=P(M.rot([-.020,.024],d.q));
 body+=line(T0,T1,'#9caea9',5)+muscleRibbon(T1,P([-.030,.060]),18,'#9faeaa','triceps')+line(P([-.030,.060]),TD,'#c9c0a7',6);
 body+=line(B0,B1,'#cfc09b',5)+muscleRibbon(B1,B2,21*Math.sqrt(.31/d.path),'url(#muscleFade)','biceps')+line(B2,D,'#c8b387',5);
 const brA=P([.020,.124]),brB=P(M.rot([.021,-.034],d.q));body+=muscleRibbon(brA,brB,11,'#b7604c','brachialis');
 const radA=P([.045,.050]),radB=P(M.rot([.026,-.23],d.q));body+=muscleRibbon(radA,radB,9,'#cf8967','brachioradialis');
 const H=P(d.hand);if(!d.device)body+=`<g transform="translate(${H})"><rect x="-19" y="-7" width="38" height="14" rx="5" fill="#647b77"/><rect x="-28" y="-18" width="12" height="36" rx="4" fill="#93a9a2"/><rect x="16" y="-18" width="12" height="36" rx="4" fill="#93a9a2"/></g>`;
 if(!d.device)body+=arrow([H[0]+33,H[1]],[H[0]+33,H[1]+48],C.blue,3);
 if(d.device){
  body+=`<circle cx="${o[0]}" cy="${o[1]}" r="35" fill="none" stroke="#9fb5b3" stroke-width="6"/><circle cx="${o[0]}" cy="${o[1]}" r="46" fill="none" stroke="#b6c8c3" stroke-width="2" stroke-dasharray="6 5"/>`;
  body+=line(P([0,0]),P(d.wrist),'#7998a3',4)+line(P(d.wrist),H,'#7998a3',4);
 }
 body+=`<circle cx="${o[0]}" cy="${o[1]}" r="10" fill="white" stroke="${C.blue}" stroke-width="2"/>`+dot(o,C.blue,3);
 if(d.mode!=='isometric'){
  const aa=d.q+(d.work==='离心'?-.24:.24),r=.115,A=P([Math.sin(d.q)*r,-Math.cos(d.q)*r]),B=P([Math.sin(aa)*r,-Math.cos(aa)*r]);body+=arrow(A,B,C.red,3);
 }
 body+='</g>';
 body+=txt(24,126,'上臂固定',18)+txt(24,155,'前臂旋后',18)+txt(369,150,'① 肱二头肌',19,C.red)+txt(369,183,'② 肱肌（深层）',18,C.red)+txt(369,216,'③ 肱桡肌',18,C.red)+txt(24,224,'肱三头肌',18,C.muted)+txt(24,251,'相反功能参照',16,C.muted);
 body+=txt(25,602,'⊙ 近似内外侧轴垂直于图示矢状面',18,C.blue);
 let panel=`<div class="contraction-panel"><div class="panel-kicker">${d.device?'虚拟等速装置 · 恒速工作段':'主动工作 · 不等于被动拉长'}</div><h2>${esc(d.type)}</h2><p class="takeaway">${esc(d.lengthChange)} · ${esc(d.phase)}</p><div class="metric-grid"><div><small>肘屈曲角</small><b>${d.angle.toFixed(1)}°</b></div><div><small>案例时间</small><b>${d.time.toFixed(2)} s</b></div></div>`;
 if(d.device){
  panel+=`<div class="iso-box"><b>角速度 ${Math.abs(d.velocity).toFixed(0)}°/s</b><p>当前长度分类：<strong>${d.work}收缩</strong></p><div class="constant-speed"><span></span></div><small>速度由理想装置限制。改变发力级别，不改变设定速度。</small></div><p>装置反力矩响应（相对示意）</p><div class="effort-bar"><i style="width:${(100*d.deviceTorque).toFixed(1)}%"></i></div><small>反力矩是定义的模型响应，不是实测肌力。</small>`;
 }else{
  panel+=`<div class="step-explain"><strong>${d.mode==='isometric'?'关节不动，仍在维持负荷':d.work==='离心'?'肘在伸展，屈肘肌仍在控制':'肘在屈曲，屈肘肌缩短发力'}</strong><p>${d.mode==='isometric'?'拖动的是保持时间；角度与附着位置保持。':d.work==='离心'?'重力方向不变，肌肉主动制动；不要当成放松后自由下落。':'肌腹与肌腱连接到同一条运动链，前臂和手不会单独飘走。'}</p></div>`;
 }
 panel+=`<div class="muscle-summary"><b>重点肌群</b><p>肱二头肌、肱肌、肱桡肌</p><small>肱三头肌作为相反功能参照；未推断其激活或完全放松。</small></div><p class="subtle">${d.device?'等速是速度条件，可与向心、离心组合；不是按口令匀速就等于仪器等速训练。':'彩色肌腹为机制示意，精细外形可另看“局部骨肌”。无真实肌纤维变形求解。'}</p></div>`;
 return `<div class="contraction-layout"><div class="contraction-figure">${svg(a.title,body,575,626)}</div>${panel}</div>`;
}
function failure(a,t,s){
 const E=s.material==='B'?12000:6000,ey=s.material==='B'?.0017:.004,ef=s.material==='B'?.0045:.012,sy=E*ey,sm=s.material==='B'?26:34,e=t*.014;
 const stress=x=>x>ef?0:x<=ey?E*x:sy+(sm-sy)*(x-ey)/(ef-ey);
 const sig=stress(e),broken=e>ef,phase=broken?'试件失效':e>ey?'非线性阶段（模型）':'线弹性阶段';
 let b=txt(28,42,a.title,26)+txt(28,73,'示例材料有明确的模型屈服/失效点；这些数值不代表人体骨折阈值。',18,C.muted);
 const y=158,de=82*t;
 b+=`<rect x="225" y="${y-de/2}" width="99" height="${292+de}" rx="18" fill="#e5d8bd" stroke="#9b8a72" stroke-width="2"/>`;
 if(!broken)b+=arrow([274,y-de/2],[274,y-de/2-42*sig/sm],C.red)+arrow([274,y+292+de/2],[274,y+292+de/2+42*sig/sm],C.red);
 if(broken)b+=`<path d="M220 284 L243 275 L262 295 L281 282 L308 304 L330 292 L330 315 L308 327 L282 306 L262 319 L243 299 L220 308Z" fill="${C.bg}"/>`;
 b+=txt(275,552,phase,24,broken?C.red:C.green,'middle')+box(596,110,490,417)+txt(629,150,'应力 σ（MPa）',20)+line([641,466],[1040,466],C.ink,2)+line([641,466],[641,176],C.ink,2);
 const xy=(x,y)=>[641+390*x/.014,466-265*y/40];
 const points=Array.from({length:141},(_,i)=>xy(.014*i/140,stress(.014*i/140)));
 b+=poly(points,C.green,4)+dot(xy(e,sig),C.red,7)+txt(1015,501,'应变 ε',18,C.ink,'end')+txt(648,207,`E=${E} MPa（定义值）`,19)+txt(648,244,`屈服 ${sy.toFixed(1)} / 峰值 ${sm.toFixed(0)} MPa`,18)+txt(32,593,'材料B更硬不等于更强。切换材料，分别比较斜率、峰值与失效应变。',18,C.muted);
 return svg(a.title,b,1120,622);
}
function archStructure(a,t,s){
 const view=s.archView||'medial',name={medial:'内侧纵弓',lateral:'外侧纵弓',transverse:'横弓'}[view]||'内侧纵弓';
 let b=txt(28,43,'足弓结构 · '+name,26)+txt(28,75,'骨性结构用于定位；弓线表示结构关系，不是实测韧带路径。',18,C.muted);
 if(view==='transverse'){
  b+=txt(76,137,'前足横断教学示意',22);
  for(let i=0;i<5;i++){const x=111+i*84,y=350-52*Math.sin(i/4*Math.PI);b+=`<ellipse cx="${x}" cy="${y}" rx="28" ry="36" fill="#ded5c0" stroke="#9f927c" stroke-width="2"/>`+txt(x,y+70,`${i+1}`,20,C.ink,'middle');}
  b+=`<path d="M111 350 Q280 245 447 350" stroke="${C.green}" stroke-width="4" fill="none" stroke-dasharray="8 6"/>`+txt(96,461,'由内侧至外侧排列；不是俯视压力图。',18,C.muted);
 }else{
  const o=[487,389],sc=1570;b+=`<g transform="translate(${o})">${registered('heel_tarsals',sc)}${registered('heel_toes',sc)}</g>`;
  b+=`<path d="M192 409 Q315 ${view==='medial'?205:331} 514 409" stroke="${C.green}" stroke-width="5" fill="none" stroke-dasharray="8 6"/>`;
  b+=txt(90,147,'真实足骨网格侧投影',22)+txt(90,179,'叠加弓形关系，避免把整足当成一块骨。',17,C.muted);
 }
 const content=view==='medial'?['骨性联系','跟骨、距骨、舟骨、楔骨及内侧跖骨。','功能理解','由骨形态、韧带/腱膜及主动肌肉共同支持。']:view==='lateral'?['骨性联系','跟骨、骰骨及外侧跖骨。','功能理解','与内侧纵弓协同承重；不能只凭弓高判定刚度。']:['骨性联系','跗骨及跖骨横向排列形成的弓形结构。','功能理解','是三维结构，不等于某一张足迹中段的面积比。'];
 b+=box(626,121,459,402)+txt(650,167,content[0],23)+wrap(650,205,content[1],20,20)+txt(650,300,content[2],23)+wrap(650,338,content[3],20,20)+wrap(650,449,'负重回弹请切回“足弓负重与回弹”；足迹面积与压力轨迹分别操作。',24,17);
 return svg(a.title,b,1120,600);
}

// Split existing wide illustrations into stacked, readable viewports on phones.
// One state renders both views; no separate/mobile physiology is calculated.
function responsive(a,raw){
 if(!raw.startsWith('<svg'))return raw;
 const wide=raw.replace('<svg ','<svg class="diagram-desktop" ');
 let scenes={lever:[18,88,594,535],bone:[20,92,560,450],capacity:[15,90,570,450],ssc:[15,90,570,455],breath:[15,88,565,442],airflow:[55,105,660,405],gait:[15,90,565,492],cardio:[22,97,584,438],failure:[16,93,559,485],archStructure:[16,94,591,418]};
 let stats={lever:[637,95,470,511],bone:[602,96,493,515],capacity:[585,98,503,466],ssc:[580,100,512,458],breath:[602,98,490,406],airflow:[742,113,350,415],gait:[583,98,522,475],cardio:[620,103,473,437],failure:[590,98,504,440],archStructure:[615,104,480,432]};
 if(a.kind==='breath'&&['thoracic','abdominal'].includes(a.motion)){stats.breath=[527,88,565,442];}
 const q=scenes[a.kind],r=stats[a.kind];if(!q||!r)return raw;
 const narrow=(box,cls)=>raw.replace(/viewBox="[^"]+"/,`viewBox="${box.join(' ')}"`).replace('<svg ',`<svg class="${cls}" `);
 return `<div class="responsive-diagram">${wide}<div class="diagram-mobile">${narrow(q,'mobile-scene')}${narrow(r,'mobile-data')}</div></div>`;
}
const render=function(a,t,s={}){
 if(a.kind==='contraction')return contraction(a,t,s);
 let raw=a.kind==='lever'?lever(a,t,s):a.kind==='failure'?failure(a,t,s):a.kind==='archStructure'?archStructure(a,t,s):previous(a,t,s);
 return responsive(a,raw);
};
Sim.render=render;g.LabRefined={lever,contraction,failure,archStructure,responsive,registered};
})(window);


/* Functional references; these are not calculated activation or muscle-force outputs. */
window.MuscleTeaching={
 shoulderFlex:{plane:'矢状面',axis:'近似内外侧轴',direction:'上臂向前抬起；反向为从前屈位回落',joint:'肱骨相对肩带/躯干',main:['三角肌前部','喙肱肌','胸大肌锁骨部'],opposite:['三角肌后部','背阔肌','大圆肌'],note:'此源模型展示前举与返回；不把源抬举轴负角度当成已验证肩后伸。'},
 shoulderAbd:{plane:'冠状面',axis:'近似前后轴',direction:'上臂离开躯干；反向回到体侧',joint:'肱骨相对肩带，伴肩胛联动',main:['三角肌中部','冈上肌'],opposite:['胸大肌','背阔肌','大圆肌'],note:'肩胛上回旋还需前锯肌与斜方肌等协作。源模型比例不是所有人的固定肩肱节律。'},
 shoulderRotate:{plane:'垂直肱骨长轴的局部平面',axis:'肱骨长轴',direction:'内旋与外旋，前臂仅用来帮助观察',joint:'肱骨相对肩胛骨',main:['内旋：肩胛下肌、胸大肌等'],opposite:['外旋：冈下肌、小圆肌'],note:'这里只列方向对应的功能肌；不能据高亮认定肌肉正在产生某个力。'},
 elbow:{plane:'矢状面',axis:'近似内外侧轴',direction:'前臂靠近上臂为屈曲；远离为伸展',joint:'前臂相对上臂',main:['肱二头肌','肱肌','肱桡肌'],opposite:['肱三头肌','肘肌'],note:'重物缓慢下降时，屈肘肌可做离心控制；不能仅因肘在伸展就认定伸肌是唯一主动工作肌。'},
 forearm:{plane:'前臂局部横断参考面',axis:'前臂纵向斜轴',direction:'旋前：桡骨跨向尺骨前方；旋后：返回相对平行',joint:'桡骨相对尺骨',main:['旋前圆肌','旋前方肌'],opposite:['旋后肌','肱二头肌'],note:'与肩内外旋、腕桡尺偏区分。'},
 wristFlex:{plane:'矢状面',axis:'近似内外侧轴',direction:'手向掌侧弯为掌屈，向背侧为背伸',joint:'手相对前臂',main:['桡侧腕屈肌','尺侧腕屈肌','掌长肌'],opposite:['桡侧腕长/短伸肌','尺侧腕伸肌'],note:'先固定前臂，单独看腕；手指动作不能代替腕运动。'},
 wristDev:{plane:'冠状面（解剖参考）',axis:'近似前后轴',direction:'向拇指侧为桡偏，向小指侧为尺偏',joint:'手相对前臂',main:['桡偏：桡侧腕屈肌、桡侧腕伸肌群'],opposite:['尺偏：尺侧腕屈肌、尺侧腕伸肌'],note:'方向组合来自不同肌肉，不是每个方向各只有一块肌肉。'},
 fingerFlex:{plane:'各指局部屈伸平面',axis:'各关节近似内外侧轴',direction:'掌指及指间屈曲 / 伸展分开观察',joint:'相邻指骨与掌骨',main:['指浅屈肌','指深屈肌'],opposite:['指伸肌','部分手内在肌经伸肌腱装置参与'],note:'蚓状肌和骨间肌可协同掌指屈曲、指间伸展，不能把它们只归入单一屈伸组。'},
 fingerSpread:{plane:'手掌平面',axis:'掌指关节近似掌背向轴',direction:'四指以中指为参考展开 / 并拢',joint:'近节指骨相对掌骨',main:['背侧骨间肌'],opposite:['掌侧骨间肌'],note:'小指展肌等也参与特定指的动作；源联动是教学预设。'},
 thumb:{plane:'多平面组合',axis:'第一腕掌及相关关节多轴',direction:'拇指向指腹对合方向移动 / 复位',joint:'第一腕掌及拇指相关关节',main:['拇对掌肌','拇短展肌','拇短屈肌等'],opposite:['拇长/短伸肌等协助复位'],note:'源模型多轴路径不作为精确指腹接触或对掌等级判定。'},
 trunkFlex:{plane:'矢状面',axis:'近似内外侧轴',direction:'躯干前屈 / 后伸；滑条正负采用源模型坐标',joint:'躯干相对骨盆',main:['腹直肌','腹斜肌群'],opposite:['竖脊肌','多裂肌等'],note:'站立前屈可由背肌离心制动，不能将此功能表当作实时肌肉激活。'},
 trunkSide:{plane:'冠状面',axis:'近似前后轴',direction:'躯干向左或右侧屈',joint:'躯干相对骨盆',main:['同侧腰方肌','同侧躯干肌群协作'],opposite:['对侧相应躯干肌群'],note:'本程序分出方向，不给出精确单肌力或左右正常阈值。'},
 trunkRotate:{plane:'横断面',axis:'纵轴',direction:'躯干相对骨盆左旋 / 右旋',joint:'躯干相对骨盆',main:['腹内、外斜肌按方向协作'],opposite:['反向旋转的斜肌组合'],note:'向右旋可由左腹外斜肌与右腹内斜肌协作；多块深层肌也参与稳定。'},
 hipFlex:{plane:'矢状面',axis:'近似内外侧轴',direction:'股骨向前为前屈，向后为后伸',joint:'股骨相对骨盆',main:['髂腰肌','股直肌等'],opposite:['臀大肌','腘绳肌群（除股二头肌短头）'],note:'双关节肌的作用同时受相邻关节位置影响。'},
 hipAbd:{plane:'冠状面',axis:'近似前后轴',direction:'大腿远离中线为外展，靠近为内收',joint:'股骨相对骨盆',main:['臀中肌','臀小肌'],opposite:['内收长肌','内收短肌','大收肌等'],note:'单腿支撑时这些肌群还可参与骨盆稳定；骨盆偏移不直接证明单肌无力。'},
 hipRotate:{plane:'垂直股骨长轴的局部平面',axis:'股骨长轴',direction:'股骨内旋 / 外旋',joint:'股骨相对骨盆',main:['内旋：臀中/小肌前部、阔筋膜张肌等'],opposite:['外旋：臀大肌及髋深层外旋肌群'],note:'肌肉功能受髋角度影响；不能只看脚尖就把变化归为踝关节。'},
 knee:{plane:'以矢状面为主',axis:'随姿态变化的膝关节轴',direction:'屈膝：小腿接近大腿后侧；伸膝：返回伸直',joint:'小腿相对大腿',main:['屈膝：腘绳肌群等'],opposite:['伸膝：股四头肌'],note:'下降控制和起身推进需结合负重。源模型含伴随变换，不等于个体真实接触求解。'},
 ankle:{plane:'以矢状面为主',axis:'距小腿关节斜轴',direction:'背伸：脚背靠近小腿；跖屈：脚尖远离',joint:'足相对小腿',main:['背伸：胫骨前肌','趾长伸肌等'],opposite:['跖屈：腓肠肌、比目鱼肌等'],note:'跖屈型动作不要只看足趾下压；小腿三头肌可在不同阶段推进或控制。'},
 footTilt:{plane:'多平面耦合',axis:'后足斜轴（源模型）',direction:'足底转向内侧 / 外侧',joint:'后足相关关节',main:['内翻：胫骨后肌、胫骨前肌等'],opposite:['外翻：腓骨长肌、腓骨短肌等'],note:'不把足内翻等同跖屈，不用单一倾斜程度判断扭伤级别。'},
 toeFlex:{plane:'局部屈伸面',axis:'跖趾段近似内外侧轴',direction:'合并足趾段屈曲 / 伸展',joint:'源模型合并足趾段相对前足',main:['趾长/短屈肌等'],opposite:['趾长/短伸肌等'],note:'趾间关节未逐一绑定；这个模型不能用于逐趾精细动作考核。'}
};

MuscleTeaching.mvc=MuscleTeaching.elbow;MuscleTeaching.isokinetic=MuscleTeaching.elbow;


/* V8.6. Original illustrative task models. Numeric outputs are synthetic, not patient data.
 * Each drawing and readout shares its state function. */
(function(g){'use strict';
const M=g.LabMath,esc=g.Sim.esc,prev=g.Sim.render.bind(g.Sim),C={ink:'#29435d',muted:'#71869a',blue:'#356fcf',teal:'#188f88',orange:'#ce825b',red:'#bd654d',bone:'#e5d7bb',line:'#dbe5ef'};
const tx=(x,y,s,size=16,c=C.ink,anchor='start')=>`<text x="${x}" y="${y}" font-size="${size}" fill="${c}" text-anchor="${anchor}">${esc(s)}</text>`;
const ln=(a,b,c=C.blue,w=3,dash='')=>`<line x1="${a[0]}" y1="${a[1]}" x2="${b[0]}" y2="${b[1]}" stroke="${c}" stroke-width="${w}" stroke-linecap="round" ${dash?`stroke-dasharray="${dash}"`:''}/>`;
const poly=(p,c=C.blue,w=3,dash='')=>`<polyline points="${p.map(x=>x.join(',')).join(' ')}" fill="none" stroke="${c}" stroke-width="${w}" stroke-linecap="round" stroke-linejoin="round" ${dash?`stroke-dasharray="${dash}"`:''}/>`;
const circle=(p,r,c)=>`<circle cx="${p[0]}" cy="${p[1]}" r="${r}" fill="${c}"/>`;
const box=(x,y,w,h,c='#fff',r=9)=>`<rect x="${x}" y="${y}" width="${w}" height="${h}" rx="${r}" fill="${c}" stroke="${C.line}"/>`;
const arrow=(a,b,c=C.blue,w=2.3)=>{let q=Math.atan2(b[1]-a[1],b[0]-a[0]);return ln(a,b,c,w)+poly([[b[0]-9*Math.cos(q-.4),b[1]-9*Math.sin(q-.4)],b,[b[0]-9*Math.cos(q+.4),b[1]-9*Math.sin(q+.4)]],c,w);};
const svg=(title,body,w=500,h=420)=>`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${w} ${h}" role="img" aria-label="${esc(title)}"><g font-family="system-ui,Microsoft YaHei,sans-serif">${body}</g></svg>`;
const metric=(label,value)=>`<div><small>${esc(label)}</small><b>${esc(value)}</b></div>`;
const tile=(title,sub,figure,footer,extra='')=>`<section class="scene-tile ${extra}"><header><div><h3>${esc(title)}</h3><small>${esc(sub)}</small></div></header><div class="figure">${figure}</div><footer>${footer}</footer></section>`;
const title=(h,p,chip='操作观察')=>`<div class="learn-title"><div><h2>${esc(h)}</h2><p>${esc(p)}</p></div><span class="panel-chip">${esc(chip)}</span></div>`;
const clamp=(n,a=0,b=1)=>M.clamp(Number(n),a,b),smooth=x=>{x=clamp(x);return x*x*(3-2*x);};
const registered=(name,scale)=>g.LabRefined.registered(name,scale);

M.chainCompare=function(t){const u=smooth(t),beta=(90-78*u)*Math.PI/180,L=.40;
 const H=[-.36,.68],K=[.04,.68],A=[K[0]+L*Math.cos(beta),K[1]-L*Math.sin(beta)],td=[Math.sin(beta),Math.cos(beta)];
 const open={hip:H,knee:K,ankle:A,heel:[A[0]-.055*td[0],A[1]-.055*td[1]],toe:[A[0]+.16*td[0],A[1]+.16*td[1]],shoulder:[H[0]-.025,1.12],head:[H[0]-.025,1.26],kneeAngle:90-78*u};
 const q=(5+26*u)*Math.PI/180,r=(5+48*u)*Math.PI/180,cA=[0,.07],cK=[L*Math.sin(q),.07+L*Math.cos(q)],cH=[cK[0]-L*Math.sin(r),cK[1]+L*Math.cos(r)],tr=(5+28*u)*Math.PI/180;
 const sh=[cH[0]+.43*Math.sin(tr),cH[1]+.43*Math.cos(tr)];
 const closed={hip:cH,knee:cK,ankle:cA,heel:[-.055,0],toe:[.16,0],shoulder:sh,head:[sh[0]+.01,sh[1]+.14],kneeAngle:(q+r)*180/Math.PI};
 return {t,u,open,closed};};
function chain(a,t){const d=M.chainCompare(t),pnt=(p,closed=false)=>[(closed?258:196)+p[0]*224,370-p[1]*224];
 function fig(leg,closed){const P=p=>pnt(p,closed),col=closed?C.teal:C.blue;let b=ln([40,377],[458,377],C.line,2);
 if(!closed){b+=box(64,215,128,14,'#dce5ee',4)+ln([72,230],[72,373],'#95aabd',7)+ln([180,230],[180,373],'#95aabd',7)+ln([69,216],[69,121],'#95aabd',7);}
 b+=ln(P(leg.hip),P(leg.shoulder),'#6689a7',28)+circle(P(leg.head),19,'#caa585');
 if(closed){const h=P(leg.hip),k=P(leg.knee),an=P(leg.ankle);b+=poly([[h[0]-16,h[1]], [k[0]-12,k[1]+1],[an[0]-15,an[1]]],'#b3c9d7',13);}
 b+=poly([leg.hip,leg.knee,leg.ankle].map(P),col,19)+ln(P(leg.heel),P(leg.toe),col,12);
 const S=P(leg.shoulder),armEnd=closed?[S[0]+80,S[1]+33]:[P(leg.hip)[0]+60,P(leg.hip)[1]-15];b+=poly([S,[S[0]+28,S[1]+44],armEnd],'#83a4bd',11);
 for(const p of [leg.hip,leg.knee,leg.ankle])b+=circle(P(p),6,'white')+circle(P(p),3.5,col);
 const thighMid=M.mul(M.add(leg.hip,leg.knee),.5),tm=P(thighMid);b+=ln([tm[0]-8,tm[1]-10],[tm[0]+13,tm[1]+17],C.orange,11);
 if(closed){const mid=P(M.mul(M.add(leg.knee,leg.ankle),.5));b+=ln([mid[0]-13,mid[1]-12],[mid[0]-18,mid[1]+15],C.orange,9);b+=circle(P(leg.hip),12,C.orange);b+=ln(P(leg.heel),P(leg.toe),C.teal,5)+tx(327,360,'足底固定',15,C.teal);}
 let hist=Array.from({length:41},(_,i)=>P(closed?M.chainCompare(i/40).closed.hip:M.chainCompare(i/40).open.toe));b+=poly(hist,col,2,'4 5');const curr=P(closed?leg.hip:leg.toe);b+=circle(curr,6,col);
 b+=tx(26,30,closed?'远端：足与地面接触位置不变':'远端：足在空间中移动',16,col)+tx(26,57,closed?'髋、膝、踝联动':'本例：大腿相对固定，膝伸展',14,C.muted);
 b+=tx(26,405,closed?`当前膝屈曲 ${leg.kneeAngle.toFixed(0)}° · 下蹲阶段`:`当前膝屈曲 ${leg.kneeAngle.toFixed(0)}° · 伸膝阶段`,15,col);
 return svg(closed?'闭链深蹲':'开链坐姿伸膝',b,500,426);}
 return `<div class="learn-panel">${title('开链与闭链，看固定的是哪一端','比较腿屈伸与深蹲的远端约束。虚线仅为当前点的运动轨迹。')}<div class="comparison-grid">${tile('开链 · 坐姿腿屈伸','Open kinetic chain',fig(d.open,false),'<strong>重点肌群：</strong>股四头肌。<br>足游离，小腿与足相对大腿运动；不等于仅一块肌肉参与。')}${tile('闭链 · 深蹲','Closed kinetic chain',fig(d.closed,true),'<strong>控制任务：</strong>髋膝伸肌、小腿肌群与躯干协作。<br>足底固定，骨盆位置与相邻关节角度一起变化。')}</div><p class="lesson-callout">先判断远端约束，再分析动作。闭链也包含关节旋转；开链也可有多肌协作。当前例图不把“旋转 / 直线”当作开闭链的定义。</p></div>`;
}

M.deviceTrial=function(mode,t,s={}){const u=clamp(t),mvc=mode==='mvc',speed=clamp(s.isoSpeed??60,30,180),rev=s.isoDirection==='eccentric',duration=mvc?7:85/speed,ready=s.deviceStarted===true;
 const effort=clamp(s.deviceEffort??.7,.2,1),shape=s.deviceProfile==='steady'?1:s.deviceProfile==='wave'?.72+.25*Math.sin(2*Math.PI*u):.35+.65*smooth(u);
 const angle=mvc?75:20+85*(rev?1-u:u),q=angle*Math.PI/180;
 const envelope=v=>v<.12?0:v<.42?smooth((v-.12)/.30):v<.72?1-.025*Math.sin(Math.PI*(v-.42)/.3)**2:v<.92?1-smooth((v-.72)/.2):0;
 const maximum=(s.deviceSide==='left'?42:45)*(s.mvcEffort==='submax'?.6:1),forceAt=v=>maximum*envelope(v)/.265;
 const force=ready&&mvc?forceAt(u):0,peak=ready&&mvc?Math.max(force,...Array.from({length:201},(_,i)=>forceAt(u*i/200))):0;
 const userTorque=ready?(mvc?force*.265:40*effort*shape*(.8+.2*Math.sin(q))):0;
 return {u,mvc,ready,speed,duration,time:ready?u*duration:0,angle:mvc?75:ready?angle:rev?105:20,velocity:ready&&!mvc?(rev?-speed:speed):0,userTorque,deviceTorque:-userTorque,force,peak,forceAt,effort:effort*shape,reverse:rev,phase:s.deviceStopped?'试次已停止 · 保留末次值':!ready?'等待固定与校零':mvc?(u<.12?'准备':u<.42?'渐增用力':u<.72?'峰值保持':u<.92?'逐渐放松':'试次结束'):'恒定角速度工作段',condition:s.mvcEffort==='submax'?'次最大努力对照':'最大努力指令（合成）'};};
function deviceFigure(d){const o=[292,255],sc=440,P=p=>[o[0]+p[0]*sc,o[1]-p[1]*sc],ang=d.angle*Math.PI/180,wrist=M.rot([0,-.265],ang),hand=M.rot([0,-.39],ang);let b='';
 b+=ln([34,434],[514,434],C.line,2)+box(78,304,166,20,'#b9cce0',5)+ln([99,324],[99,430],'#728ba5',10)+ln([224,324],[224,430],'#728ba5',10)+box(70,146,22,164,'#acc2d9',8);
 b+=poly([[202,318],[265,340],[272,419]],'#9bb2c8',19)+ln([260,423],[305,423],'#65819b',12);
 b+=ln([204,311],[230,142],'#809cb7',40)+circle([233,104],24,'#cda887')+ln([244,146],P([0,.29]),'#9bb2c8',13);
 b+=ln([193,189],[237,199],'#536d86',12)+ln([190,281],[215,287],'#536d86',13)+tx(36,258,'固定带',13,C.muted)+ln([89,254],[185,277],C.muted,1.5);
 b+=box(410,258,55,94,'#d6e2ef',8)+ln([438,352],[438,427],'#7992ac',15)+ln([390,429],[480,429],'#7992ac',9)+ln(o,[428,275],'#8faabd',16)+circle(o,24,'#dce9f4')+circle(o,15,'white');
 b+=`<g transform="translate(${o})" opacity=".8">${registered('elbow_humerus',sc)}</g><g transform="translate(${o}) rotate(${-d.angle})" opacity=".8">${registered('elbow_forearm',sc)}</g>`;
 const B=P([.025,.236]),D=P(M.rot([.008,-.047],ang));b+=`<path d="M${B} Q${B[0]+26} ${B[1]+48} ${D}" stroke="#c97d62" fill="none" stroke-width="14" stroke-linecap="round"/>`;
 b+=ln(o,P(hand),'#5f83a6',6)+circle(P(wrist),9,'#6a89a7')+ln([P(wrist)[0]-7,P(wrist)[1]-12],[P(wrist)[0]+7,P(wrist)[1]+12],'#45647f',8)+circle(o,7,'white')+circle(o,3,C.blue);
 b+=tx(33,33,'肘关节测力装置 · 侧面示意',18)+tx(34,57,'右侧骨网格作定位；左右数据为独立合成案例。',12,C.muted);
 b+=tx(334,89,'轴线对齐',14,C.blue)+poly([[355,98],[325,158],o],C.blue,1.5,'4 4');
 if(d.mvc){b+=tx(365,356,'力臂锁定',14,C.blue)+tx(365,379,'角度保持75°',13,C.muted);}else{b+=tx(356,363,'测力头',14)+tx(356,388,'顺应性响应',13,C.muted);}
 const H=P(wrist);if(d.ready&&d.userTorque>0){const tang=[Math.cos(ang),-Math.sin(ang)],length=18+38*d.userTorque/45;
 b+=arrow([H[0]-7,H[1]-14],[H[0]-7+length*tang[0],H[1]-14+length*tang[1]],C.blue,3)+arrow([H[0]+7,H[1]+14],[H[0]+7-length*tang[0],H[1]+14-length*tang[1]],C.orange,3);}
 b+=tx(34,463,'蓝：受试者作用方向　橙：装置反作用方向',13,C.muted);return svg('虚拟肌力测试装置',b,540,480);
}
function dynamometer(a,t,s){const mode=a.motion,d=M.deviceTrial(mode,t,s);let screen=`<div class="device-state ${d.ready?'':'idle'}">${esc(d.phase)}<small>${d.mvc?'固定角度，用力增加时前臂不继续移动。':d.reverse?'装置带动伸肘，屈肘肌主动离心制动。':'保持匀速工作段，装置抵抗屈肘方向。'}</small></div><div class="metric-grid">${metric('肘角度',d.angle.toFixed(1)+'°')}${metric('试次时间',d.time.toFixed(2)+' s')}</div>`;
 if(d.mvc){
  let chart=ln([48,154],[450,154],C.line,1.5)+ln([48,26],[48,154],C.line,1.5)+tx(48,17,'合成力 / N',12,C.muted)+tx(450,185,'试次时间 / s',12,C.muted,'end');
  let sample=Array.from({length:121},(_,i)=>[48+402*i/120,154-d.forceAt(i/120)/190*117]);chart+=poly(sample,'#d7e1ee',2,'4 4');
  if(d.ready){const n=Math.max(2,Math.ceil(121*t));const xs=Array.from({length:n},(_,i)=>{const u=t*i/(n-1);return [48+402*u,154-d.forceAt(u)/190*117];});chart+=poly(xs,C.blue,3)+circle([48+402*t,154-d.force/190*117],4,C.blue);}
  chart+=tx(48,175,'0',12,C.muted)+tx(449,175,'7',12,C.muted,'end');screen+=`<div class="metric-grid">${metric('当前力',d.force.toFixed(1)+' N')}${metric('已观察峰值',d.peak.toFixed(1)+' N')}</div><div class="data-tile">${svg('等长力—时间曲线',chart,500,200)}</div><p class="numeric-note">${esc(d.condition)}。这里具体演示最大随意等长收缩（MVIC）；次最大试次不判作真正MVC。</p>`;
 }else{
  screen+=`<div class="metric-grid">${metric('设定角速度',d.speed+'°/s')}${metric('当前工作段速度',Math.abs(d.velocity).toFixed(0)+'°/s')}</div><div class="data-tile"><h3>发力与装置响应</h3><div class="bar-row"><span>受试者力矩</span><div class="bar-track"><i style="width:${100*d.userTorque/45}%"></i></div><b>${d.userTorque.toFixed(1)}</b></div><div class="bar-row"><span>装置力矩大小</span><div class="bar-track"><i style="width:${100*Math.abs(d.deviceTorque)/45}%"></i></div><b>${Math.abs(d.deviceTorque).toFixed(1)}</b></div><p class="numeric-note">单位 N·m，数值来自理想平衡假设；两者方向相反。</p></div>`;
  const pts=Array.from({length:61},(_,i)=>{const q=M.deviceTrial(mode,i/60,{...s,deviceStarted:true});return [36+406*i/60,124-q.userTorque/45*94];});let chart=tx(36,17,'力矩大小（合成）',12,C.muted)+ln([36,125],[442,125],C.line,1.5)+poly(pts,'#d7e1ee',2,'4 4');if(d.ready){const n=Math.max(2,Math.ceil(61*t));const seen=Array.from({length:n},(_,i)=>{const u=t*i/(n-1),q=M.deviceTrial(mode,u,{...s,deviceStarted:true});return [36+406*u,124-q.userTorque/45*94];});chart+=poly(seen,C.blue,2.7)+circle([36+406*t,124-d.userTorque/45*94],4,C.orange);}screen+=`<div class="data-tile">${svg('匀速段的力矩变化',chart,480,145)}<div class="mini-legend"><span><i></i>蓝线为已观察段；虚线为本合成条件的预览。</span></div></div>`;
 }
 return `<div class="learn-panel">${title(d.mvc?'最大随意收缩 · 等长测力':'等速肌力测试 · 让阻力跟随发力',d.mvc?'最大随意收缩测试（MVC）；本页明确采用等长条件。':'装置外形用于理解测试系统，不仿制具体品牌的控制程序。','合成试次')}<div class="device-layout"><div class="scene-tile"><div class="figure">${deviceFigure(d)}</div></div><div class="device-screen">${screen}</div></div><p class="panel-foot">本页无真实传感器。匀速段忽略加减速、重力补偿和设备摩擦，不能据图中力值评价个人肌力。实际测试须按设备规程与专业评估执行。</p></div>`;
}

M.candleTrial=function(t,s={},durationOverride=null,distanceOverride=null){const time=15*clamp(t),distance=clamp(distanceOverride??s.airDistance??55,10,80),level=clamp(s.airLevel??2,1,3),duration=clamp(durationOverride??s.airSeconds??8,3,15),active=time>.15&&time<duration;
 const potential=Math.min(58,20*level*(25/distance)**1.4),factor=active?smooth(Math.min((time-.15)/.3,(duration-time)/.25)):0,deflection=potential*factor;
 const threshold=15,reaches=potential>=threshold;let valid=0;
 // Integrate the same display threshold over the observed part; preserves duration after release.
 const step=.025;for(let x=.0125;x<Math.min(time,duration);x+=step){const f=x>.15?smooth(Math.min((x-.15)/.3,(duration-x)/.25)):0;if(potential*f>=threshold)valid+=Math.min(step,time-(x-.0125));}
 return{time,distance,level,duration,active,deflection,potential,validSeconds:Math.max(0,valid),reaches,threshold};};
function candleFig(d){const x=145+d.distance*3.25,y=142,bend=d.deflection*.72,flicker=1+Math.sin(d.time*6)*.018;let b=ln([34,287],[455,287],C.line,2);
 b+=`<path d="M40 131 Q59 110 82 129 L107 147 L82 165 Q54 174 40 155Z" fill="#d7b08d"/>`;
 if(d.active){for(let i=0;i<3;i++){const end=Math.min(x-13,122+200*(25/d.distance)**.35);b+=arrow([119,132+i*16],[end,132+i*16],i===1?C.teal:'#a8cfd0',i===1?2.4:1.5);}}
 b+=box(x-13,172,26,114,'#e8e2d8',5)+ln([x,173],[x,162],'#607281',2.4);
 b+=`<path d="M${x} 165 C${x-18+bend*.2} ${143} ${x-8+bend*.7} ${122*flicker} ${x+bend} ${106*flicker} C${x+9+bend*.8} ${135} ${x+23+bend*.1} ${154} ${x} 165Z" fill="#f5ae43"/><path d="M${x} 165 Q${x-10+bend*.3} 146 ${x+bend*.65} 127 Q${x+15+bend*.3} 158 ${x} 165Z" fill="#ffe3a0"/>`;
 b+=ln([107,312],[x,312],C.blue,1.6)+ln([107,306],[107,318],C.blue,1.4)+ln([x,306],[x,318],C.blue,1.4)+tx((107+x)/2,341,d.distance+' cm',20,C.blue,'middle');
 b+=tx(25,29,d.active?'正在输出':'尚未输出 / 已放松',15,d.active?C.teal:C.muted)+tx(25,55,`相对输出 ${d.level}级 · 火苗偏转 ${d.deflection.toFixed(0)}°`,13,C.muted);
 b+=tx(25,382,d.active?(d.deflection>=d.threshold?'达到明显偏转阈值':'未达到明显偏转阈值'):'火苗保持燃烧，未作吹灭判定',14,d.active&&d.deflection>=d.threshold?C.teal:C.muted);
 return svg('虚拟蜡烛火苗与距离',b,470,399);}
function candles(a,t,s){const distance=a.motion==='distance',left=M.candleTrial(t,s,distance?8:4,distance?15:s.airDistance??35),right=M.candleTrial(t,s,distance?8:s.airSeconds??8,s.airDistance??55);
 return `<div class="learn-panel">${title(distance?'同样用力，近处与远处有什么不同？':'同样距离，比较持续输出时间',distance?'两栏是独立试次：近处固定15 cm，远处可调；没有前一支蜡烛挡住气流。':'两栏保持相同输出级别与目标距离，只比较持续时长。','虚拟火焰')}<div class="comparison-grid">${tile(distance?'近距离试次':'短时试次','独立条件 A',candleFig(left),`<strong>${distance?'目标15 cm':'输出4 s'}</strong> · 当前达到偏转要求 ${left.validSeconds.toFixed(1)} s`)}${tile(distance?'远距离试次':'较长试次','独立条件 B',candleFig(right),`<strong>${distance?'目标'+right.distance+' cm':'输出'+right.duration+' s'}</strong> · 当前达到偏转要求 ${right.validSeconds.toFixed(1)} s`)}</div><p class="lesson-callout">${distance?'调整远处距离，比较能否让火苗达到同一个偏转要求；不是把没吹到的火苗画成熄灭。':'输出结束后保留已经达到要求的时间，不把累计时间清零。'}</p><p class="panel-foot">气流衰减与15°观察阈值均为教学设定，不换算成真实呼吸肌力、口腔压或临床耐力。仅在虚拟环境中观察火焰。</p></div>`;
}

const pressureNodes=[[0,.86,1],[-.19,.54,.42],[.22,.52,.75],[-.23,.27,1.1],[.23,.26,.9],[-.25,.075,.9],[.12,.105,.38]];
M.footPressure=function(u,kind='typical'){if(u<=0||u>=1)return{force:0,cop:null,weights:pressureNodes.map(()=>0)};
 const shift=[.10,.40,.40,.73,.69,.95,.89],width=[.16,.19,.20,.19,.20,.10,.13],amp=kind==='low'?[1,1.9,.9,1.05,.85,.95,.4]:kind==='high'?[1.2,.12,.28,.95,1.28,.65,.35]:[1,.48,.8,1.05,1,.85,.4];
 const env=Math.sin(Math.PI*u);const w=pressureNodes.map((p,i)=>env*amp[i]*Math.exp(-.5*((u-shift[i])/width[i])**2));const f=w.reduce((a,b)=>a+b,0),cop=[0,1].map(k=>w.reduce((z,x,i)=>z+x*pressureNodes[i][k],0)/f);return{force:f,cop,weights:w};};
M.footprintWalk=function(t,kind='typical'){const time=4.2*clamp(t),feet=Array.from({length:4},(_,i)=>{const onset=i*.95,u=(time-onset)/1.22;return{index:i,side:i%2?'left':'right',onset,u,started:time>=onset,done:u>=1,...M.footPressure(u,kind)};});return{time,feet};};
function footOutline(x,y,mirror=1,kind='typical'){const medial=kind==='low'?-24:kind==='high'?4:-9;
 return `<g transform="translate(${x} ${y}) scale(${mirror} 1)"><path d="M-18 99 C-36 93 -31 75 -27 65 Q${medial} 54 -28 35 C-40 19 -30 4 -10 7 C4 0 30 7 32 26 Q38 43 21 65 C19 80 29 102 5 106 Q-9 110 -18 99Z" fill="#e9f0f7" stroke="#cbd8e6" stroke-width="1.5"/>${[[-22,-2,10],[-3,-5,8],[13,0,7],[27,7,6],[36,17,5]].map(([a,b,r])=>`<ellipse cx="${a}" cy="${b}" rx="${r*.75}" ry="${r}" fill="#e9f0f7" stroke="#cbd8e6"/>`).join('')}</g>`;}
function footprintFig(t,kind){const d=M.footprintWalk(t,kind);let b=arrow([157,425],[157,44],'#d4deeb',1.5)+tx(172,40,'前进',12,C.muted);
 for(const foot of d.feet){const x=foot.side==='right'?221:99,y=370-foot.index*88,mirror=foot.side==='right'?1:-1;let part=footOutline(x,y,mirror,kind);
  if(foot.started){const u=clamp(foot.u,.001,.999);const P=p=>[x+mirror*p[0]*85,y+p[1]*103-7];const history=Array.from({length:51},(_,i)=>M.footPressure(.005+(u-.005)*i/50,kind).cop).filter(Boolean).map(P);
   if(foot.force>0){pressureNodes.forEach((n,i)=>{const w=foot.weights[i];if(w>.04){const at=P(n);part+=`<ellipse cx="${at[0]}" cy="${at[1]}" rx="${9+7*w}" ry="${12+9*w}" fill="${w>.6?'#e9a463':'#6eb9c5'}" opacity="${Math.min(.78,.25+w*.35)}"/>`;}});}
   part+=poly(history,kind==='low'?C.orange:kind==='high'?C.blue:C.teal,2.2);if(foot.cop)part+=circle(P(foot.cop),4,C.red)+circle(P(foot.cop),1.6,'white');
  }
  b+=`<g opacity="${!foot.started?.16:foot.done?.52:1}">${part}</g>`+tx(x+mirror*44,y+55,foot.side==='right'?'右':'左',11,C.muted,'middle');
 }
 b+=tx(25,503,'红点：当前COP　彩线：已走过的轨迹',12,C.muted);return svg('动态足印压力中心',b,320,520);}
function footprints(a,t,s){const labels=[['low','低弓情景','假设中足接触较多'],['typical','常见弓高情景','假设后足至前足平稳转移'],['high','高弓情景','假设中足接触较少']];
 return `<div class="learn-panel">${title('足印向前走，压力中心怎样转移？','足印落地后保持原位；同一足只有接触时才出现当前COP。','合成压力场景')}<div class="comparison-grid three">${labels.map(([k,n,d])=>tile(n,d,footprintFig(t,k),'<strong>COP</strong>由本情景的压力权重计算。<br>本曲线不是此足型的唯一或标准轨迹。','foot-tile')).join('')}</div><p class="lesson-callout">这里显示的是足底压力中心 COP，不是全身质心 COM。足型、接触分布和步态共同影响观测结果，不能只看弓高就推定一条“重心线”。</p></div>`;
}
function cycleFig(t,fast){const cycle=t*(fast?10:6)*2*Math.PI,hip=[211,173],hub=[283,290];let b=ln([37,377],[461,377],C.line,2)+poly([[121,373],[hub[0],hub[1]],[365,373]],'#8ba8bd',9)+ln([hip[0],185],[177,318],'#9bb4c7',7)+ln([179,184],[233,184],'#55748f',10)+circle(hub,37,'#e0eaf4')+circle(hub,26,'#f7faff')+ln([344,255],[364,131],'#9bb4c7',7)+ln([350,135],[390,135],'#6686a1',9);
 const sh=[247,91],head=[266,48];b+=ln(hip,sh,fast?'#e0a07e':'#7aa3c8',28)+circle(head,19,'#c7a17f')+poly([sh,[301,115],[354,134]],'#a2bacd',10);
 for(let i=1;i>=0;i--){const pedal=[hub[0]+34*Math.cos(cycle+i*Math.PI),hub[1]+34*Math.sin(cycle+i*Math.PI)],H=[hip[0],hip[1],0],A=[pedal[0],pedal[1],0],K=M.kneeFromEnds(H,A,104,103,1);b+=ln(hub,pedal,'#6e8aa1',4)+poly([hip,K.slice(0,2),pedal],i?'#bacbda':fast?C.orange:C.blue,13)+ln([pedal[0]-10,pedal[1]],[pedal[0]+15,pedal[1]],'#496982',8);}
 b+=tx(30,417,fast?'高踏频短时试次（不是运动处方）':'稳定踏频持续试次（时间压缩展示）',14,C.muted);return svg('功率自行车运动条件示意',b,500,440);}
function metabolism(a,t,s){const u=clamp(t),sprintP=u<.4?0:1;let continuous=tile('中等强度持续运动','例：持续骑行；保持可持续稳定负荷',cycleFig(t,false),`<strong>示意场景时间：</strong>${Math.round(600*u)} s / 600 s<br>观察持续任务与有氧代谢为主的特征。<div class="energy-strip"><span>磷酸原参与</span><span>糖酵解参与</span><span class="primary-energy">有氧代谢主导</span></div>`);
 let sprint=tile('短时高强度运动','例：短时冲刺；强度大、时间短',cycleFig(t,true),`<strong>示意场景时间：</strong>${(10*u).toFixed(1)} s / 10 s<br>${sprintP?'观察糖酵解贡献逐步增加。':'观察磷酸原系统的快速供能。'}<div class="energy-strip"><span class="${!sprintP?'primary-energy':''}">磷酸原较突出</span><span class="${sprintP?'primary-energy':''}">糖酵解增加</span><span>有氧同样参与</span></div>`);
 return `<div class="learn-panel">${title('有氧与无氧：先看强度，再看时间','两栏共用回放进度，不是相同现实持续时间。','定性供能说明')}<div class="comparison-grid">${continuous}${sprint}</div><p class="lesson-callout">主导供能随任务变化，并不等于其他系统关闭。本页不把心率、气促或动作名称当作实测乳酸阈，也不编造氧耗、乳酸或供能百分比。</p><p class="panel-foot">两个持续时间是为便于比较设置的虚拟场景。动作只是示范，不要求真实患者进行冲刺或极量测试。</p></div>`;
}
Sim.render=function(a,t,s={}){if(a.kind==='chain')return chain(a,t,s);if(a.kind==='dynamometer')return dynamometer(a,t,s);if(a.kind==='candle')return candles(a,t,s);if(a.kind==='footprints')return footprints(a,t,s);if(a.kind==='metabolism')return metabolism(a,t,s);return prev(a,t,s);};
g.LearningPanels={chain,dynamometer,candles,footprints,metabolism};
})(window);


/* V8.7 — four scoped revisions. No measured participant data are generated.
   Units and baselines are explicit. Renderers and readouts call these same functions. */
(function(g){'use strict';
const M=g.LabMath,C=(x,a=0,b=1)=>Math.max(a,Math.min(b,Number.isFinite(+x)?+x:a));
M.breath87=function(t,mode='thoracic'){
 const u=C(t),q=(1-Math.cos(2*Math.PI*u))/2;
 // Screen units are deliberately exaggerated; they are NOT millimetres or percentages of a lung volume.
 const amplitudes={thoracic:{width:84,drop:14,abd:8},abdominal:{width:18,drop:68,abd:29},quiet:{width:30,drop:32,abd:13},forced:{width:38,drop:46,abd:18}};
 const a=amplitudes[mode]||amplitudes.quiet;
 const phase=u===0?'呼气末 · 起始':Math.abs(u-.5)<1e-6?'吸气末':u===1?'呼气末':u<.5?'吸气':'呼气';
 return {t:u,q,mode,phase,inhaling:u>0&&u<.5,exhaling:u>.5&&u<1,
  width0:180,width:180+a.width*q,roof:115,apex0:285,apex:285+a.drop*q,
  edge:365+8*q,abdBulge:a.abd*q,deltaWidth:a.width*q,deltaHeight:a.drop*q,
  lengthUnits:'illustration_units_not_measurement',forcedExpiration:mode==='forced'&&u>.5};
};
function material(name){const B=name==='B';const E=B?12000:6000,ey=B?.0017:.004,ef=B?.0045:.012,sy=E*ey,peak=B?26:34;
 const ratio=E*(ef-ey)/(peak-sy);let lo=.000001,hi=20;
 for(let i=0;i<80;i++){const m=(lo+hi)/2;if(m/(1-Math.exp(-m))<ratio)lo=m;else hi=m;}
 return {name:B?'B':'A',E,ey,ef,sy,peak,k:(lo+hi)/2};
}
M.material87=material;
M.materialStress87=function(e,name='A'){
 const m=typeof name==='object'?name:material(name),x=C(e,0,m.ef);
 if(x<=m.ey)return m.E*x;
 const r=(x-m.ey)/(m.ef-m.ey);
 return m.sy+(m.peak-m.sy)*(-Math.expm1(-m.k*r))/(-Math.expm1(-m.k));
};
M.failure87=function(t,s={}){
 const u=C(t),m=material(s.material),mode=['elasticUnload','plasticUnload'].includes(s.failureMode)?s.failureMode:'loading';
 const peakEps=mode==='elasticUnload'?.70*m.ey:mode==='plasticUnload'?m.ey+.68*(m.ef-m.ey):m.ef;
 const peakSigma=M.materialStress87(peakEps,m),residual=Math.max(0,peakEps-peakSigma/m.E);
 const unloading=mode!=='loading'&&u>.5,broken=mode==='loading'&&u>=1;
 let eps,sig;if(mode==='loading'){eps=u*m.ef;sig=M.materialStress87(eps,m);}else if(!unloading){eps=peakEps*2*u;sig=M.materialStress87(eps,m);}else{sig=peakSigma*(2-2*u);eps=peakEps-(peakSigma-sig)/m.E;}
 const atYield=Math.abs(eps-m.ey)<1e-10;
 const phase=broken?'试件失效 / 曲线在断裂点终止':unloading?(u===1?(residual>1e-9?'卸载完成 · 保留残余变形':'卸载完成 · 回到原长'):'卸载 · 弹性部分回复'):atYield?'屈服点 B':eps<m.ey?'线弹性区 A—B':'塑性区 B—C · 弯曲上升';
 return {t:u,material:m,mode,eps,stress:sig,broken,unloading,phase,peakEps,peakSigma,residual,
  currentResidual:unloading?Math.max(0,eps-sig/m.E):Math.max(0,eps-sig/m.E),
  engineeringLength_mm:100*(1+eps),referenceLength_mm:100,deformationMagnification:20,
  currentLoad_N:broken?0:sig*10,area_mm2:10,graphStress:sig,synthetic:true};
};
const gaitCache=new Map();
M.gaitBounds87=function(pattern='normal'){
 if(gaitCache.has(pattern))return gaitCache.get(pattern);
 const samples=Array.from({length:241},(_,i)=>M.gait(i/240,pattern));
 const bounds=idx=>{const xs=samples.map(s=>s.com[idx]);const min=Math.min(...xs),max=Math.max(...xs);return {min,max,mid:(min+max)/2,span:(max-min)*1000};};
 const b={samples,vertical:bounds(1),lateral:bounds(2)};gaitCache.set(pattern,b);return b;
};
M.gaitSync87=function(t,pattern='normal'){
 const d=M.gait(C(t),pattern),b=M.gaitBounds87(pattern),h=.00005;
 const prev=M.gait((d.t-h+1)%1,pattern),next=M.gait((d.t+h)%1,pattern);
 const dy=(next.com[1]-prev.com[1])/(2*h),dz=(next.com[2]-prev.com[2])/(2*h);
 const leftOff=(d.shift+d.stanceL)%1,p=d.t===1?0:d.t;
 let phase=p===0?'右初始接触（事件）':p<leftOff?'右负重反应':p<(leftOff+d.shift)/2?'右支撑中期':p<d.shift?'右支撑末期':p<d.stanceR?'右预摆动':p<d.stanceR+(1-d.stanceR)/3?'右摆动初期':p<d.stanceR+2*(1-d.stanceR)/3?'右摆动中期':'右摆动末期';
 return {d,b,pattern,t:d.t,time:d.t*d.duration,phase,support:d.support,
  vertical_mm:(d.com[1]-b.vertical.min)*1000,lateral_mm:(d.com[2]-b.lateral.mid)*1000,
  verticalTrend:Math.abs(dy)<.00015?'上下转折附近':dy>0?'正在上升':'正在下降',
  lateralTrend:Math.abs(dz)<.00015?'侧向转折附近':dz>0?'正在向右迁移':'正在向左迁移',
  frame:'x forward, y up, z participant right; rear view screen right = participant right',
  synthetic:true};
};
// Literature summary is never used to rescale the model. Ranges = max minus min over a cycle.
M.gaitReference87={citation:'Orendurff et al., JRRD 2004;41(6A):829–834',doi:'10.1682/JRRD.2003.10.0150',n:10,
 rows:[{speed:0.7,verticalMean:27.4,verticalSD:5.2,lateralMean:69.9,lateralSD:13.4},
       {speed:1.6,verticalMean:48.3,verticalSD:9.2,lateralMean:38.5,lateralSD:14.1}],unit:'mm',
 note:'Study means ± SD, not diagnostic limits; not the animation time series.'};
})(typeof window!=='undefined'?window:globalThis);


/* V8.7 focused visual revisions. Original diagrams; no external fonts/CDN/assets.
   All panels read one normalized time from App; numerical definitions in revision87_math.js. */
(function(g){'use strict';
const M=g.LabMath,esc=g.Sim.esc,oldRender=g.Sim.render.bind(g.Sim);
const C={blue:'#326ad4',orange:'#c6763f',red:'#bf654f',teal:'#188b86',ink:'#263e58',muted:'#687f94',line:'#dce6ef',ghost:'#adbfcf',bone:'#d8c8a8',bg:'#f6f9fd'};
const tx=(x,y,s,size=17,c=C.ink,anc='start')=>`<text x="${x}" y="${y}" font-size="${size}" fill="${c}" text-anchor="${anc}">${esc(s)}</text>`;
const line=(a,b,c=C.blue,w=2,dash='')=>`<line x1="${a[0]}" y1="${a[1]}" x2="${b[0]}" y2="${b[1]}" stroke="${c}" stroke-width="${w}" stroke-linecap="round" ${dash?`stroke-dasharray="${dash}"`:''}/>`;
const poly=(ps,c=C.blue,w=2,dash='')=>`<polyline points="${ps.map(p=>p.join(',')).join(' ')}" fill="none" stroke="${c}" stroke-width="${w}" stroke-linecap="round" stroke-linejoin="round" ${dash?`stroke-dasharray="${dash}"`:''}/>`;
const dot=(p,c=C.red,r=5)=>`<circle cx="${p[0]}" cy="${p[1]}" r="${r}" fill="${c}"/>`;
const arrow=(a,b,c=C.blue,w=2.6)=>{const q=Math.atan2(b[1]-a[1],b[0]-a[0]),h=8;return line(a,b,c,w)+poly([[b[0]-h*Math.cos(q-.46),b[1]-h*Math.sin(q-.46)],b,[b[0]-h*Math.cos(q+.46),b[1]-h*Math.sin(q+.46)]],c,w);};
const dbl=(a,b,c=C.blue,w=2.6)=>arrow(a,b,c,w)+arrow(b,a,c,w);
const box=(x,y,w,h,fill='#fff',r=10)=>`<rect x="${x}" y="${y}" width="${w}" height="${h}" rx="${r}" fill="${fill}" stroke="${C.line}"/>`;
const svg=(name,b,w=520,h=440,attrs='')=>`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${w} ${h}" role="img" aria-label="${esc(name)}" ${attrs}><g font-family="system-ui,Microsoft YaHei,sans-serif">${b}</g></svg>`;
const heading=(title,sub)=>`<div class="v87-heading"><h2>${esc(title)}</h2><p>${esc(sub)}</p></div>`;
const metric=(label,v,unit='')=>`<div class="v87-metric"><span>${esc(label)}</span><strong>${esc(v)}<small>${esc(unit)}</small></strong></div>`;
function thorax(d,focus){
 const cx=199,w=d.width/2,xL=cx-w,xR=cx+w,cy=240,main=focus==='thoracic'?C.blue:C.orange;
 let b='';
 // Baseline silhouette and dome are not changed by current phase.
 const contour=(l,r)=>`M${cx-42} 105 Q${l+3} 125 ${l} 198 L${l+5} 310 Q${l+9} 354 ${cx-51} 384 L${cx+51} 384 Q${r-9} 354 ${r-5} 310 L${r} 198 Q${r-3} 125 ${cx+42} 105 Z`;
 b+=`<path d="${contour(xL,xR)}" fill="#eaf1fa" stroke="#698fbe" stroke-width="2.5"/>`;
 b+=`<path d="${contour(cx-90,cx+90)}" fill="none" stroke="${C.ghost}" stroke-width="2" stroke-dasharray="5 5"/>`;
 b+=`<path d="M${cx-74} 365 Q${cx} 205 ${cx+74} 365" fill="none" stroke="${C.ghost}" stroke-width="2" stroke-dasharray="5 5"/>`;
 // Each symmetric rib arc moves outward and slightly upward; sternum is only a locator.
 for(let i=0;i<6;i++){
  const yy=145+i*30-9*d.q,edge=i<2?.83:1;
  b+=`<path d="M${cx} ${yy+21} Q${cx-w*edge} ${yy+21} ${cx-w*edge+8} ${yy-5} M${cx} ${yy+21} Q${cx+w*edge} ${yy+21} ${cx+w*edge-8} ${yy-5}" fill="none" stroke="${C.bone}" stroke-width="7.5" stroke-linecap="round"/>`;
 }
 b+=line([cx,126-8*d.q],[cx,318-8*d.q],C.bone,8);
 // Quadratic Bezier apex at t=.5 is exactly d.apex.
 b+=`<path d="M${xL+17} ${d.edge} Q${cx} ${2*d.apex-d.edge} ${xR-17} ${d.edge} L${xR-17} ${d.edge+9} Q${cx} ${2*d.apex-d.edge+9} ${xL+17} ${d.edge+9} Z" fill="#efbd98" stroke="${C.orange}" stroke-width="2" data-diaphragm-apex="${d.apex.toFixed(4)}"/>`;
 b+=line([cx-89,90],[cx+89,90],C.ghost,1.5,'4 4');
 b+=dbl([xL,76],[xR,76],focus==='thoracic'?C.blue:'#9bb5d4',focus==='thoracic'?4:2);
 b+=line([xL,68],[xL,101],C.blue,1.4)+line([xR,68],[xR,101],C.blue,1.4);
 b+=tx(cx,50,'左右径',21,focus==='thoracic'?C.blue:C.muted,'middle');
 const ax=focus==='thoracic'?365:348;
 b+=dbl([ax,d.roof],[ax,d.apex],focus==='abdominal'?C.orange:'#bac4d1',focus==='abdominal'?4:2);
 b+=line([cx+50,d.roof],[ax+8,d.roof],C.line,1.2)+line([cx+12,d.apex],[ax+8,d.apex],C.orange,1.2,'4 4');
 b+=tx(ax+11,154,'上',20,focus==='abdominal'?C.orange:C.muted)+tx(ax+11,179,'下',20,focus==='abdominal'?C.orange:C.muted)+tx(ax+11,204,'径',20,focus==='abdominal'?C.orange:C.muted);
 if(d.inhaling){b+=arrow([xL+6,250],[xL-30,250],C.blue,focus==='thoracic'?4:2)+arrow([xR-6,250],[xR+30,250],C.blue,focus==='thoracic'?4:2)+arrow([cx,d.apex0+4],[cx,d.apex+18],C.orange,focus==='abdominal'?4:2);}
 else if(d.exhaling){b+=arrow([xL-28,250],[xL+4,250],C.blue,2.5)+arrow([xR+28,250],[xR-4,250],C.blue,2.5)+arrow([cx,d.apex+20],[cx,d.apex-3],C.orange,3);}
 b+=dot([cx,d.apex],C.orange,4)+tx(47,414,'虚线：呼气末原位',16,C.muted)+tx(47,438,'蓝：胸廓宽度　橙：胸腔纵向空间',15,C.muted);
 // Side cutaway: abdominal bulge is anterior, NOT the width of an abdomen in frontal view.
 b+=box(411,280,99,166,'#fff',9)+tx(460,303,'侧面',15,C.muted,'middle');
 b+=`<path d="M435 320 L435 420 L472 420 Q${484+d.abdBulge*.56} 382 470 353 L472 320" fill="#eef4fc" stroke="#859db5" stroke-width="2"/>`;
 b+=`<path d="M435 340 Q453 ${327+22*d.q} 472 340" fill="none" stroke="${C.orange}" stroke-width="3"/>`;
 b+=`<path d="M470 353 Q484 382 472 420" fill="none" stroke="${C.ghost}" stroke-width="1.5" stroke-dasharray="4 4"/>`;
 if(d.inhaling)b+=arrow([478,382],[494+Math.min(9,d.abdBulge*.2),382],C.orange,2);
 b+=tx(460,435,'腹壁前移',12,C.muted,'middle');
 return svg(focus==='thoracic'?'胸式呼吸左右径扩张':'腹式呼吸膈肌下降与上下径扩大',b,520,463,
  `data-breath-frame="${d.t}" data-mode="${focus}" data-width="${d.width}" data-height="${d.apex-d.roof}"`);
}
function breathing(a,t){
 const pair=['thoracic','abdominal'].includes(a.motion),d=M.breath87(t,pair?a.motion:a.motion==='forced'?'forced':'quiet');
 const cards=pair?['thoracic','abdominal']:[a.motion==='forced'?'forced':'quiet'];
 const titles={thoracic:'胸式呼吸',abdominal:'腹式呼吸',quiet:'自然呼吸',forced:'用力呼气参照'};
 return `<div class="v87-panel" data-revision="breath87">${heading(pair?'胸式与腹式：看清主要增加的径线':titles[cards[0]],pair?'同一呼吸周期、同一基线。主变化放大显示，不给人体毫米测量值。':'吸气时胸腔扩大，呼气时回复；箭头随呼吸方向改变。')}
  <div class="v87-phase"><b>${esc(d.phase)}</b><span>${pair?'两图共用播放进度；暂停可对照虚线原位。':'暂停可对照呼气末的虚线基线。'}</span></div>
  <div class="v87-breath-grid ${pair?'':'single'}">${cards.map(k=>{const q=M.breath87(t,k),f=k==='abdominal'?'abdominal':'thoracic';return `<section class="v87-card ${a.motion===k?'chosen':''}"><header><span class="v87-kicker">${k==='abdominal'?'膈肌下降更突出':k==='thoracic'?'肋骨上提、胸廓扩张':k==='forced'?'主动呼气 · 腹肌等参与':'胸廓与膈肌协调运动'}</span><h3>${titles[k]}</h3><p>${k==='abdominal'?'重点：胸腔上下径增大':k==='thoracic'?'重点：胸廓左右径扩大':k==='forced'?'重点：比较准备吸气与主动呼气':'重点：自然吸气与被动呼气'}</p></header><div class="v87-figure">${thorax(q,f)}</div><footer>${k==='abdominal'?'膈肌顶部向下，胸腔纵向空间增加；侧面小窗显示腹壁向前隆起。':k==='forced'?'前半程为吸气准备，后半程为主动用力呼气；腹肌等参与呼气。':'肋骨向上、向外运动，胸廓增宽；胸式呼吸也伴前后径变化。'}<small>${k==='abdominal'?'腹壁起伏不等于腹肌主动吸气。':'膈肌并未停止参与。图示肋骨数量为示意，不用于数骨。'}</small></footer></section>`;}).join('')}</div>
  <p class="v87-note">${pair?'这里比较的是相对突出程度，不把胸式与腹式画成两套互斥机制。':'安静呼气以弹性回缩为主；主动用力呼气另有肌肉参与。'}几何幅度为便于观察而设定，不表示精确器官形变。</p></div>`;
}
function walkingFigure(sync,view){
 const {d,b}=sync,side=view==='side',scale=239;
 const P=p=>side?[202+(p[0]-d.stride*d.t)*scale,442-p[1]*scale]:[209+p[2]*scale,442-p[1]*scale];
 let z=line([26,444],[412,444],C.line,2);
 if(side){for(let i=-3;i<=5;i++){let x=202+(i*.3-d.stride*d.t)*scale;if(x>18&&x<419)z+=line([x,449],[x+15,449],'#c8d6e5',2);}}
 else z+=line([209,48],[209,444],C.line,1.4,'5 5');
 // Same complete pose, different orthographic projections; no second gait clock.
 const L=d.legs.left,R=d.legs.right;
 for(const [name,col,w]of[['left',C.blue,13],['right',C.teal,15]]){
  const leg=d.legs[name],a=d.arms[name];
  z+=poly([leg.hip,leg.knee,leg.ankle].map(P),col,w);
  if(side)z+=line(P(leg.heel),P(leg.toe),col,9);
  else {const h=P(leg.heel),toe=P(leg.toe);z+=line([h[0]-10,h[1]],[toe[0]+10,toe[1]],col,8);}
  if(leg.contact){const pp=P(leg.ankle);z+=line([pp[0]-(side?14:15),446],[pp[0]+(side?32:15),446],col,3);}
  for(const p of [leg.knee,leg.ankle])z+=dot(P(p),'#fff',2.7);
 }
 if(side)z+=line(P(d.hip),P(d.shoulder),'#52789b',25);
 else{const h=P(d.hip),s=P(d.shoulder);z+=`<path d="M${s[0]-36} ${s[1]-5} Q${s[0]} ${s[1]-18} ${s[0]+36} ${s[1]-5} L${h[0]+25} ${h[1]} Q${h[0]} ${h[1]+9} ${h[0]-25} ${h[1]} Z" fill="#7595b5"/>`;}
 z+=dot(P(d.head),'#cba581',17);
 for(const [name,col]of[['left',C.blue],['right',C.teal]]){const arm=d.arms[name];z+=poly([arm.s,arm.elbow,arm.wrist].map(P),col,9);}
 const com=P(d.com);z+=`<circle cx="${com[0]}" cy="${com[1]}" r="11" fill="white" stroke="${C.red}" stroke-width="1.2"/>`+dot(com,C.red,6);
 // Actual (not magnified) min/max reference levels on the body.
 if(side){const ya=442-b.vertical.max*scale,yb=442-b.vertical.min*scale;
  z+=line([com[0]+14,ya],[350,ya],C.red,1,'4 4')+line([com[0]+14,yb],[350,yb],C.red,1,'4 4');
  z+=line([350,ya],[350,yb],C.red,2.5)+line([344,ya],[356,ya],C.red,2)+line([344,yb],[356,yb],C.red,2);
  z+=tx(319,ya-28,'上下全幅',16,C.red)+tx(319,ya-8,b.vertical.span.toFixed(1)+' mm',16,C.red);
 }else{const xl=209+b.lateral.min*scale,xr=209+b.lateral.max*scale;
  z+=line([xl,com[1]+17],[xl,336],C.red,1,'4 4')+line([xr,com[1]+17],[xr,336],C.red,1,'4 4');
  z+=line([xl,336],[xr,336],C.red,2.5)+line([xl,331],[xl,341],C.red,2)+line([xr,331],[xr,341],C.red,2);
  z+=tx(286,329,'左右全幅',16,C.red)+tx(286,351,b.lateral.span.toFixed(1)+' mm',16,C.red)+line([xr+3,338],[278,338],C.red,1.2);
  z+=tx(86,424,'左',15,C.blue)+tx(322,424,'右',15,C.teal);
 }
 z+=tx(28,28,side?'前进方向 →':'后面观察：左右同受试者',15,C.muted);
 return svg(side?'同步侧面步行动画':'同步后面步行动画',z,440,469,
 `data-gait-view="${view}" data-phase="${d.t}" data-com-y="${d.com[1]}" data-com-z="${d.com[2]}"`);
}
function timeChart(sync,axis){
 const b=sync.b,idx=axis==='vertical'?1:2,bd=axis==='vertical'?b.vertical:b.lateral,col=axis==='vertical'?C.orange:C.blue;
 const base=axis==='vertical'?bd.min:bd.mid;
 const lo=axis==='vertical'?0:(bd.min-base)*1000,hi=(bd.max-base)*1000;
 const pad=Math.max((hi-lo)*.12,1),min=lo-pad,max=hi+pad;
 const X=u=>63+409*u,Y=mm=>191-(mm-min)/(max-min)*139;
 let z='';
 for(let i=0;i<3;i++){const v=lo+(hi-lo)*i/2,y=Y(v);z+=line([63,y],[472,y],C.line,1,'3 4')+tx(54,y+4,v.toFixed(1),13,C.muted,'end');}
 z+=line([63,203],[472,203],C.muted,1.5)+line([63,45],[63,203],C.muted,1.5);
 const pts=b.samples.map((x,i)=>[X(i/(b.samples.length-1)),Y((x.com[idx]-base)*1000)]);
 z+=poly(pts,'#c4d4e5',2);
 const seen=pts.filter((p,i)=>i/(pts.length-1)<=sync.t);seen.push([X(sync.t),Y((sync.d.com[idx]-base)*1000)]);z+=poly(seen,col,3);
 const cur=[X(sync.t),Y((sync.d.com[idx]-base)*1000)];z+=line([cur[0],43],[cur[0],203],C.red,1.5,'4 4')+dot(cur,C.red,5);
 for(const u of [0,.25,.5,.75,1])z+=tx(X(u),224,Math.round(100*u)+'%',13,C.muted,'middle');
 z+=tx(20,25,axis==='vertical'?'上下位移：距本周期最低点 / mm':'左右位移：相对本周期中线 / mm',17,col);
 return svg('COM'+axis+'同步时间曲线',z,520,243,`data-com-chart="${axis}" data-phase="${sync.t}"`);
}
function contactChart(sync){
 const d=sync.d,pattern=sync.pattern;
 let z='';const X=u=>110+354*u;
 for(const [side,y,c,n]of[['right',39,C.teal,'右足'],['left',80,C.blue,'左足']]){
  z+=tx(26,y+6,n,17,c)+line([110,y],[464,y],C.line,10);
  // contact set sampled from exactly the same model, with independent segments after gaps
  let seg=[];
  for(let i=0;i<=200;i++){let q=M.gait(i/200,pattern);if(q.legs[side].contact)seg.push([X(i/200),y]);else if(seg.length){z+=poly(seg,c,9);seg=[];}}
  if(seg.length)z+=poly(seg,c,9);z+=dot([X(sync.t),y],C.red,5);
 }
 z+=line([X(sync.t),20],[X(sync.t),96],C.red,1.5)+tx(110,117,'0%',13,C.muted)+tx(464,117,'100%',13,C.muted,'end');
 return svg('左右足接触同步条',z,500,133,`data-contact-phase="${sync.t}"`);
}
function comInset(sync){
 const b=sync.b,ratio=3.0,midY=b.vertical.mid,midZ=b.lateral.mid;
 const P=p=>[250+(p[2]-midZ)*1000*ratio,113-(p[1]-midY)*1000*ratio];
 let z=line([170,113],[342,113],C.line,1,'3 4')+line([250,35],[250,188],C.line,1,'3 4');
 z+=poly(b.samples.map(s=>P(s.com)),'#c5d4e4',2)+poly(b.samples.filter((s,i)=>i/(b.samples.length-1)<=sync.t).map(s=>P(s.com)),C.blue,2.5)+dot(P(sync.d.com),C.red,6);
 const left=250-b.lateral.span*ratio/2,right=250+b.lateral.span*ratio/2,top=113-b.vertical.span*ratio/2,bottom=113+b.vertical.span*ratio/2;
 z+=`<rect x="${left}" y="${top}" width="${right-left}" height="${bottom-top}" fill="none" stroke="${C.ghost}" stroke-dasharray="4 5"/>`;
 z+=dbl([left,200],[right,200],C.blue,1.8)+dbl([360,top],[360,bottom],C.orange,1.8);
 z+=tx(250,228,'左右全幅 '+b.lateral.span.toFixed(1)+' mm',16,C.blue,'middle')+tx(370,108,b.vertical.span.toFixed(1),16,C.orange)+tx(370,131,'mm',15,C.orange);
 z+=tx(29,28,'同步质心局部放大窗',19)+tx(29,58,'位移比例放大，非身体摆幅放大',13,C.muted);
 return svg('质心冠状面放大轨迹',z,500,241,`data-com-inset-phase="${sync.t}"`);
}
function functionalTiming(sync){
 const bands=[['胫骨前肌',[[0,.10],[.60,1]]],['小腿三头肌',[[.12,.60]]],['股四头肌',[[0,.20]]],['腘绳肌',[[0,.12],[.85,1]]]];
 let z=tx(22,30,'右侧肌群的功能任务时序（不是肌电）',17);
 for(let i=0;i<bands.length;i++){let [name,segs]=bands[i],y=63+i*34;z+=tx(22,y+5,name,15)+line([147,y],[466,y],C.line,8);for(const [a,b]of segs)z+=line([147+319*a,y],[147+319*b,y],C.teal,8);}
 z+=line([147+319*sync.t,43],[147+319*sync.t,178],C.red,1.8);
 return svg('右侧肌群功能时序',z,500,193);
}
function walking(a,t,s){
 const pattern=['clearance','lean','short'].includes(a.motion)?a.motion:'normal',q=M.gaitSync87(t,pattern),b=q.b;
 return `<div class="v87-panel v87-gait" data-revision="gait87" data-phase="${q.t}">${heading('同一时刻：步行姿态与重心迁移',pattern==='normal'?'侧面看上下，后面看左右。红点和两条曲线均来自本模型同一姿态。':'这是参数化特征对照，不是真实患者步态；不据此诊断病因。')}
 <div class="v87-phase"><b>${esc(q.support)}</b><span>${esc(q.phase)} · 周期 ${(100*t).toFixed(0)}% · 案例时间 ${q.time.toFixed(2)} s</span></div>
 <div class="v87-walk-grid">
  <section class="v87-card"><header><h3>侧面 · 上下迁移</h3></header><div class="v87-walk-figure">${walkingFigure(q,'side')}</div><footer><b>${q.verticalTrend}</b><span>距本周期最低点 ${q.vertical_mm.toFixed(1)} mm</span></footer></section>
  <section class="v87-card"><header><h3>后面 · 左右迁移</h3></header><div class="v87-walk-figure">${walkingFigure(q,'rear')}</div><footer><b>${q.lateralTrend}</b><span>${q.lateral_mm>=0?'向右':'向左'}偏离中线 ${Math.abs(q.lateral_mm).toFixed(1)} mm</span></footer></section>
 </div>
 <div class="v87-statrow">${metric('本模型上下全幅',b.vertical.span.toFixed(1),'mm')}${metric('本模型左右全幅',b.lateral.span.toFixed(1),'mm')}</div>
 <p class="v87-note compact">全幅＝本周期最大值－最小值，不是累计路程，也不是单侧振幅。这里为假定分段质量模型的计算值，不能标成“正常人标准值”。</p>
 <div class="v87-gait-chartgrid"><section class="v87-card">${timeChart(q,'vertical')}${timeChart(q,'lateral')}</section><section class="v87-card"><h3 class="v87-padtitle">左右足接触与质心位置</h3>${contactChart(q)}${comInset(q)}</section></div>
 ${a.motion==='timing'?`<section class="v87-card">${functionalTiming(q)}</section>`:''}
 <div class="v87-reference"><button type="button" data-cmd="gaitReference87">健康成人文献参考</button><p>文献参考与本动画数值分开。步速不同，上下与左右全幅也会变化；不能把一组厘米数作为统一诊断界限。</p></div>
 <p class="v87-note">左右脚接触来自同一几何模型。原型仅为合成步态；红点是按明确质量比例计算的模型COM，不是COP，也不是直接将骨盆标记改名。局部放大窗使用同一组数据。</p></div>`;
}
function specimen(q){
 const L=230,top=84,base=top+L,stretch=L*q.eps*q.deformationMagnification,yEnd=base+stretch;
 let z=tx(28,30,'轴向拉伸试件（教学定义材料）',19)+tx(28,55,'仅形变量 ×20；不是骨折风险测试',14,C.muted);
 z+=`<rect x="183" y="${top}" width="92" height="${L}" rx="10" fill="none" stroke="${C.ghost}" stroke-width="2" stroke-dasharray="5 5"/>`;
 z+=`<rect x="183" y="${top}" width="92" height="${L+stretch}" rx="10" fill="#e5d8bc" stroke="#9d8c72" stroke-width="2.4"/>`;
 for(let i=1;i<5;i++){const yy=top+(L+stretch)*i/5;z+=line([190,yy],[268,yy],'#cfbd9b',1.4);}
 if(q.broken){z+=`<path d="M181 208 L197 200 L213 218 L230 205 L247 222 L276 210 L276 225 L247 237 L230 220 L213 233 L197 215 L181 223 Z" fill="white" stroke="white" stroke-width="1"/>`+tx(297,216,'断裂',17,C.red);}
 else if(q.stress>1e-8){const a=18+30*q.stress/q.material.peak;z+=arrow([230,top-2],[230,top-a],C.red,2.5)+arrow([230,yEnd+2],[230,yEnd+a],C.red,2.5);}
 z+=line([157,base],[331,base],C.ghost,1.5,'4 4')+tx(29,base+6,'原长末端',16,C.muted);
 if(stretch>1){z+=line([277,yEnd],[325,yEnd],C.orange,1.5)+line([325,base],[325,yEnd],C.orange,2.8)+line([319,base],[331,base],C.orange,2)+line([319,yEnd],[331,yEnd],C.orange,2);}
 z+=tx(28,417,q.unloading?'卸载：弹性部分回复':q.broken?'已失效：不继续绘制承载曲线':'加载：应力与形变同步',19,q.broken?C.red:C.blue);
 z+=tx(28,448,q.unloading&&q.t===1?(q.residual>1e-9?'仍比原长长：有残余变形':'回到原长：没有残余变形'):'虚线是未受载的长度基线',16,C.muted);
 return svg('同步试件形变',z,460,466,`data-specimen-eps="${q.eps}" data-material-phase="${q.phase}"`);
}
function failureChart(q){
 const m=q.material,maxX=m.ef*1.12,maxY=m.peak*1.24;
 const X=e=>65+398*e/maxX,Y=s=>322-245*s/maxY;
 let z='';
 z+=`<rect x="${X(0)}" y="68" width="${X(m.ey)-X(0)}" height="254" fill="#eef4ff"/><rect x="${X(m.ey)}" y="68" width="${X(m.ef)-X(m.ey)}" height="254" fill="#fff5eb"/>`;
 for(let i=0;i<5;i++){let y=m.peak*i/4;z+=line([65,Y(y)],[471,Y(y)],C.line,1,'3 4')+tx(55,Y(y)+4,y.toFixed(1),13,C.muted,'end');}
 for(let i=0;i<4;i++){let x=m.ef*i/3;z+=tx(X(x),345,(x*100).toFixed(2),13,C.muted,'middle');}
 z+=line([65,322],[478,322],C.ink,1.7)+line([65,66],[65,322],C.ink,1.7)+tx(24,29,'应力 σ / MPa',18)+tx(474,373,'应变 ε / %',17,C.ink,'end');
 z+=line([X(m.ey),Y(m.sy)],[X(m.ey),322],C.ghost,1.4,'4 4')+line([65,Y(m.sy)],[X(m.ey),Y(m.sy)],C.ghost,1.4,'4 4');
 const whole=Array.from({length:151},(_,i)=>[X(m.ef*i/150),Y(M.materialStress87(m.ef*i/150,m))]);
 z+=poly(whole,'#b8c9da',2,'5 4');
 const loadMax=q.mode==='loading'?q.eps:q.unloading?q.peakEps:q.eps;
 const elasticEnd=Math.min(loadMax,m.ey);z+=line([X(0),Y(0)],[X(elasticEnd),Y(m.E*elasticEnd)],C.blue,3.2);
 if(loadMax>m.ey){let p=Array.from({length:101},(_,i)=>{const e=m.ey+(loadMax-m.ey)*i/100;return[X(e),Y(M.materialStress87(e,m))];});z+=poly(p,C.orange,3.2);}
 if(q.mode!=='loading'){
  z+=line([X(q.peakEps),Y(q.peakSigma)],[X(q.residual),Y(0)],'#ae9dcc',1.7,'4 4');
  if(q.unloading)z+=line([X(q.peakEps),Y(q.peakSigma)],[X(q.eps),Y(q.stress)],'#8768ad',3);
  if(q.residual>1e-9){z+=line([X(0),395],[X(q.residual),395],'#8768ad',3)+tx(X(q.residual)/2+32.5,418,'残余应变',16,'#8768ad','middle');}
 }
 z+=dot([X(q.eps),Y(q.graphStress)],C.red,5.5);
 z+=dot([X(m.ef),Y(m.peak)],C.red,3.5)+line([X(m.ef),Y(m.peak)+4],[X(m.ef),322],C.ghost,1.2,'4 5');
 z+=tx(76,305,'A',16,C.blue)+tx(X(m.ey)+7,Y(m.sy)-7,'B 屈服',14,C.orange)+tx(X(m.ef)-4,Y(m.peak)-12,'C 断裂界限',14,C.red,'end');
 z+=tx(93,52,'弹性区',16,C.blue)+tx(306,52,'塑性区（曲线）',16,C.orange);
 z+=tx(29,447,q.broken?'点停在最后完整承载界限；破坏后不再接零值平线。':'实线：已观察过程；浅虚线：本定义材料的预览。',14,C.muted);
 return svg('弹性—塑性应力应变及卸载曲线',z,510,466,`data-curve-eps="${q.eps}" data-curve-stress="${q.graphStress}"`);
}
function failure(a,t,s){
 const q=M.failure87(t,s),m=q.material;
 return `<div class="v87-panel" data-revision="failure87">${heading('直线弹性区 → 曲线塑性区 → 失效','按课件第19—20页的分区组织；用定义材料演示，不伪装成真实骨试验数据。')}
 <div class="v87-phase"><b>${esc(q.phase)}</b><span>材料 ${m.name} · 弹性模量 ${m.E} MPa（定义值）</span></div>
 <div class="v87-failure-grid"><section class="v87-card"><div class="v87-figure">${specimen(q)}</div></section><section class="v87-card"><div class="v87-figure">${failureChart(q)}</div></section></div>
 <div class="v87-statrow">${metric(q.broken?'断裂前应力':'当前应力',q.stress.toFixed(2),'MPa')}${metric('当前应变',(q.eps*100).toFixed(3),'%')}${metric(q.unloading?'卸载残余应变':'此状态卸载后的残余应变',(q.currentResidual*100).toFixed(3),'%')}</div>
 <p class="v87-note">${q.mode==='elasticUnload'?'本试次在弹性范围内卸载，结束时回到原长。':q.mode==='plasticUnload'?'本试次先进入塑性区再卸载，结束时应力回到零，但长度不能完全恢复。':'屈服后应力随应变呈弯曲上升。可切换“弹性内卸载”与“塑性后卸载”直接比较残余变形。'}</p>
 <p class="v87-note subtle">曲线形状与教学参数为本程序定义；真实骨受组织、方向、速率等条件影响。是否发生塑性变形看卸载后能否完全恢复，不能仅凭“线弯了”判断。曲线在失效点终止，不延长成零应力水平段。</p></div>`;
}
g.Sim.render=function(a,t,s={}){if(a.kind==='breath')return breathing(a,t,s);if(a.kind==='gait')return walking(a,t,s);if(a.kind==='failure')return failure(a,t,s);return oldRender(a,t,s);};
g.Revision87={breathing,walking,failure,timeChart};
})(window);


/* BodyParts3D frame-preserving anatomical viewer. Motion is visual kinematics,
   not a physiological muscle-force/finite-element simulation. */
(function(g){'use strict';
const {ident,matmul,point,vec,norm,add,sub,mul,dot}=SM;const rotate=SM.rotate.bind(SM);
const near=(a,b)=>Math.hypot(...sub(a,b));const clamp=x=>Math.max(0,Math.min(1,x));
const shoulder=[-.160,1.315,.076],elbow=[-.219,1.036,.079],wrist=[-.249,.806,.115],hip=[-.052,.815,.090],knee=[-.079,.377,.082],ankle=[-.071,-.004,.076];
const handBone=p=>(p.regions||[]).includes('wrist-hand')&&p.structure!=='radius';
const footBone=p=>(p.regions||[]).includes('ankle-foot');
const spinal=p=>/vertebra|rib|sacrum|sternum|manubrium|xiphoid/.test(p.structure);
const isUpper=p=>['humerus','ulna','radius'].includes(p.structure)||handBone(p);
const isLower=p=>['femur','tibia','fibula','patella'].includes(p.structure)||footBone(p);
const trans=(v)=>{const m=ident();m[12]=v[0];m[13]=v[1];m[14]=v[2];return m;};
function policy(e){const mo=e.motion;let pivot=elbow,axis=[1,0,0],deg=-55,mask=p=>isUpper(p)||p.structure==='scapula',move=p=>['ulna','radius'].includes(p.structure)||handBone(p),joint='肘关节',name='屈肘',originHigh=true;
 if(mo.startsWith('shoulder')){pivot=shoulder;mask=p=>['scapula','humerus','clavicle'].includes(p.structure)||(e.id==='pecMajor'&&/rib|thoracic|sternum|manubrium|xiphoid/.test(p.structure));move=isUpper;joint='盂肱关节';deg=mo==='shoulderAbd'?-48:mo==='shoulderAdd'?20:mo==='shoulderFlex'?-45:mo==='shoulderExt'?22:mo==='shoulderIR'?-35:35;axis=['shoulderAbd','shoulderAdd'].includes(mo)?[0,0,1]:['shoulderIR','shoulderER'].includes(mo)?norm(sub(elbow,shoulder)):[1,0,0];name={shoulderAbd:'肩外展',shoulderAdd:'肩内收方向',shoulderFlex:'肩前屈',shoulderExt:'肩后伸',shoulderIR:'肩内旋',shoulderER:'肩外旋'}[mo];}
 if(mo==='elbowExt'){deg=55;name='伸肘';} // starts from a flexed reference, see pose below
 if(['pronation','supination'].includes(mo)){axis=norm(sub([-.238,.809,.105],[-.229,1.025,.080]));pivot=[-.229,1.025,.080];move=p=>p.structure==='radius'||handBone(p);deg=mo==='pronation'?-50:45;joint='近、远侧桡尺关节';name=mo==='pronation'?'前臂旋前':'前臂旋后';}
 if(mo.startsWith('wrist')){pivot=wrist;deg=mo==='wristFlex'?-30:30;move=handBone;mask=p=>isUpper(p);joint='腕关节';name=mo==='wristFlex'?'腕掌屈':'腕背伸';}
 if(mo.startsWith('finger')||mo.startsWith('thumb')){pivot=[-.267,.742,.137];deg=mo.endsWith('Ext')?25:-30;mask=p=>handBone(p)||['radius','ulna'].includes(p.structure);move=p=>handBone(p)&&/phalanx/.test(p.structure);if(mo.startsWith('thumb')){pivot=[-.28,.790,.14];move=p=>handBone(p)&&(/thumb|first-metacarpal/.test(p.structure));deg={thumbFlex:-20,thumbExt:20,thumbAbd:-22,thumbAdd:15,thumbOpp:-15}[mo];axis=['thumbAbd','thumbAdd','thumbOpp'].includes(mo)?[0,0,1]:[1,0,0];}joint=mo.startsWith('thumb')?'拇指关节复合体':'掌指与指间关节';name={fingerFlex:'屈指方向',fingerExt:'伸指方向',thumbFlex:'屈拇方向',thumbExt:'伸拇方向',thumbAbd:'拇指外展方向',thumbAdd:'拇指内收方向',thumbOpp:'拇指对掌方向'}[mo];}
 if(mo.startsWith('scap')){pivot=[0,1.29,.075];axis=[0,1,0];deg=mo==='scapRetract'?-12:12;move=p=>['scapula','clavicle'].includes(p.structure)||isUpper(p);mask=p=>['scapula','clavicle','humerus'].includes(p.structure)||/rib|thoracic|cervical/.test(p.structure);joint='肩胛胸壁功能关系';name={scapPro:'肩胛前伸',scapRetract:'肩胛后缩',scapElev:'肩胛上提'}[mo];}
 if(mo.startsWith('trunk')){pivot=[0,.93,.065];axis=[1,0,0];deg=mo==='trunkFlex'?12:-12;move=p=>spinal(p)&&!['sacrum','hip-bone'].includes(p.structure);mask=p=>spinal(p)||p.structure==='hip-bone';joint='躯干相对骨盆';name=mo==='trunkFlex'?'躯干前屈':'躯干后伸';}
 if(mo.startsWith('hip')){pivot=hip;deg={hipAbd:-25,hipAdd:15,hipFlex:-40,hipExt:18,hipER:25}[mo];axis=['hipAbd','hipAdd'].includes(mo)?[0,0,1]:mo==='hipER'?norm(sub(knee,hip)):[1,0,0];mask=p=>['hip-bone','sacrum','femur'].includes(p.structure)||(e.id==='psoas'&&/lumbar|twelfth-thoracic/.test(p.structure));move=isLower;joint='髋关节';name={hipAbd:'髋外展',hipAdd:'髋内收方向',hipFlex:'髋前屈',hipExt:'髋后伸',hipER:'髋外旋'}[mo];}
 if(mo.startsWith('knee')){pivot=knee;deg=mo==='kneeFlex'?45:-45;mask=p=>['hip-bone','sacrum','femur','patella','tibia','fibula'].includes(p.structure);move=p=>['patella','tibia','fibula'].includes(p.structure)||footBone(p);joint='膝关节';name=mo==='kneeFlex'?'屈膝':'伸膝';}
 if(mo.startsWith('ankle')||mo.startsWith('foot')||mo.startsWith('toe')){pivot=ankle;axis=mo.startsWith('foot')?[0,0,1]:[1,0,0];deg={ankleDorsi:18,anklePlantar:-25,footInvert:-12,footEvert:10,toeFlex:-18,toeExt:20}[mo];mask=p=>['tibia','fibula'].includes(p.structure)||footBone(p)||(e.id==='gastrocnemius'&&p.structure==='femur');move=footBone;joint=mo.startsWith('toe')?'跖趾与趾间关节':mo.startsWith('foot')?'足部复合关节':'距小腿关节';name={ankleDorsi:'踝背伸',anklePlantar:'踝跖屈',footInvert:'足内翻方向',footEvert:'足外翻方向',toeFlex:'屈趾方向',toeExt:'伸趾方向'}[mo];if(mo.startsWith('toe')){pivot=[-.122,-.057,.185];move=p=>footBone(p)&&/phalanx/.test(p.structure);}}
 if(e.boneRegion){const r=e.boneRegion;mask=p=>{const n=p.structure||'',regs=p.regions||[];if(r==='neck')return regs.includes('head-neck')||/cervical|clavicle|first-rib|second-rib|sternum|manubrium/.test(n);if(r==='upperBack')return /vertebra|rib|scapula|clavicle|humerus|occipital/.test(n);if(r==='thorax')return /rib|vertebra|sternum|manubrium|xiphoid|clavicle/.test(n);if(r==='trunk')return /vertebra|rib|sacrum|hip-bone|occipital/.test(n);if(r==='elbow')return isUpper(p)||n==='scapula';if(r==='hand')return isUpper(p);if(r==='pelvis')return /hip-bone|sacrum|femur|tibia|patella/.test(n);if(r==='foot')return footBone(p)||['tibia','fibula'].includes(n);return false;};move=()=>false;joint='结构定位';name='静态观察';}
 return {pivot,axis,deg,mask,move,joint,name,originHigh};
}
function motionMatrix(e,u){if(e.dynamic===false)return ident();const p=policy(e);if(e.motion==='scapElev')return trans([0,.035*u,0]);let d=p.deg*u;if(e.motion==='elbowExt')d=-55+55*u;if(e.motion==='shoulderAdd')d=-20+20*u;if(e.motion==='hipAdd')d=-15+15*u;if(e.motion==='kneeExt')d=45-45*u;return rotate(p.axis,d*Math.PI/180,p.pivot);}
function endRegions(e,geometry){const pts=[];for(const ar of geometry)for(let j=0;j<ar.length;j+=3)pts.push([ar[j],ar[j+1],ar[j+2]]);
 let seeds={supra:[[-.074,1.31,.017],[-.185,1.32,.086]],infra:[[-.09,1.27,.012],[-.18,1.305,.052]],subscap:[[-.10,1.265,.065],[-.15,1.303,.1]],teresMinor:[[-.125,1.255,.008],[-.187,1.305,.054]],teresMajor:[[-.111,1.207,.044],[-.16,1.257,.108]],gluteMed:[[-.122,.941,.04],[-.13,.824,.061]],gluteMin:[[-.116,.902,.08],[-.128,.824,.086]],gluteMax:[[-.058,.869,-.004],[-.114,.725,.075]],iliacus:[[-.09,.897,.07],[-.061,.761,.108]],psoas:[[-.037,1.03,.085],[-.061,.761,.108]],serratus:[[-.113,1.24,.18],[-.072,1.259,.045]],rhomboid:[[-.014,1.298,.016],[-.076,1.275,.022]],trapMid:[[-.017,1.35,.015],[-.139,1.322,.037]],pecMajor:[[-.028,1.25,.18],[-.163,1.257,.095]]}[e.id];
 function average(x){return x.reduce((s,p)=>s.map((v,j)=>v+p[j]),[0,0,0]).map(v=>v/x.length);}
 function nearest(seed){return average(pts.map(p=>[p,near(p,seed)]).sort((a,b)=>a[1]-b[1]).slice(0,Math.min(16,pts.length)).map(x=>x[0]));}
 if(seeds)return seeds.map(nearest);
 const axis=e.motion.startsWith('toe')?2:1,ascending=e.motion.startsWith('toe');pts.sort((a,b)=>ascending?a[axis]-b[axis]:b[axis]-a[axis]);const n=Math.max(5,Math.floor(pts.length*.018));return [average(pts.slice(0,n)),average(pts.slice(-n))];
}
function skin(positions,normals,O,I,T){const V=new Float32Array(positions.length),N=new Float32Array(normals.length),D=sub(I,O),den=dot(D,D)||1,delta=sub(point(T,I),I);let max=0;
 for(let i=0;i<positions.length;i+=3){const v=[positions[i],positions[i+1],positions[i+2]],u=clamp((dot(sub(v,O),D)/den-.12)/.76),w=u*u*(3-2*u),vw=add(v,mul(delta,w)),n0=[normals[i],normals[i+1],normals[i+2]],grad=mul(D,6*u*(1-u)/(.76*den)),divider=1+dot(grad,delta);let nn=Math.abs(divider)>.12?sub(n0,mul(grad,dot(delta,n0)/divider)):n0;
 // A small terminal attachment zone follows the distal bone exactly; the long
 // muscle belly responds to anchor displacement, not to an ankle-sized rotation
 // applied to its entire length. This is a visual deformation, not tissue physics.
 const terminal=clamp((u-.88)/.12);const rigid=point(T,v);const vv=vw.map((x,k)=>x+(rigid[k]-x)*terminal),nr=vec(T,n0);nn=norm(nn.map((x,k)=>x+(nr[k]-x)*terminal));V.set(vv,i);N.set(nn,i);max=Math.max(max,near(v,vv));}
 return {V,N,max};
}
const proto=ModelStage.prototype,oldSet=proto.set,oldUpdate=proto.update,oldFit=proto.fitAll,oldOrient=proto.orient;
proto.set=function(a,t,mode){if(a.kind!=='muscleAtlas')return oldSet.call(this,a,t,mode);this.action=a;this.mode='atlas';this.t=t;this.entry=g.MuscleCatalogue.find(e=>e.id===a.muscleID);if(!this.entry)throw Error('未找到肌肉条目');const e=this.entry,v=this.viewer;if(!v.atlasLabelCallback){const oldDraw=v.callbacks.draw;v.callbacks.draw=()=>{oldDraw?.();if(v.dataset==='atlas'&&this.action?.kind==='muscleAtlas'&&this.atlasReadout)this.atlasLabels();};v.atlasLabelCallback=true;}
 if(!this.atlasModel){this.atlasModel=JSON.parse(JSON.stringify(this.anatomy));this.atlasModel.parts.forEach(p=>{p.color=p.type==='bone'?[.86,.82,.69]:[.76,.78,.76];});v.addDataset('atlas',this.atlasModel,this.assets.anatomy.raw);}
 // Restore previously deformed buffers from immutable source before selecting another muscle.
 for(const i of this.deformedIndices||[])this.atlasUpload(i,new Float32Array(this.assets.anatomy.raw,this.anatomy.parts[i].positions,this.anatomy.parts[i].vertexCount*3),new Float32Array(this.assets.anatomy.raw,this.anatomy.parts[i].normals,this.anatomy.parts[i].vertexCount*3));
 v.use('atlas',[]);this.deformedIndices=e.partIndices.slice();this.atlasPolicy=policy(e);this.atlasSource=e.partIndices.map(i=>({i,P:new Float32Array(this.assets.anatomy.raw,this.anatomy.parts[i].positions,this.anatomy.parts[i].vertexCount*3),N:new Float32Array(this.assets.anatomy.raw,this.anatomy.parts[i].normals,this.anatomy.parts[i].vertexCount*3)}));
 [this.atlasO,this.atlasI]=endRegions(e,this.atlasSource.map(x=>x.P));
 this.atlasModel.parts.forEach((p,i)=>{p.label=e.partIndices.includes(i)?e.name:this.anatomy.parts[i].label;p.color=e.partIndices.includes(i)?[.80,.22,.15]:p.type==='bone'?[.92,.88,.76]:[.68,.72,.73];});
 v.visible=this.anatomy.parts.map((p,i)=>e.partIndices.includes(i)||(p.type==='bone'&&this.atlasPolicy.mask(p)));this.fullMask=v.visible.slice();
 v.yaw=['supra','infra','teresMinor','deltoidPost','rhomboid','trapMid','gluteMax','gluteMed','gluteMin','ilioLumbar','longissimus'].includes(e.id)?Math.PI+.35:-.45;v.pitch=e.id==='supra'?.65:e.id==='diaphragm'?.52:.18;if(['upperBack','trunk'].includes(e.boneRegion))v.yaw=Math.PI+.35;
 this.update(t);this.fitAll();};
proto.atlasUpload=function(i,P,N){const v=this.viewer,d=v.datasets.atlas,m=d.meshes[i];if(v.gl){const gl=v.gl;gl.bindBuffer(gl.ARRAY_BUFFER,m.p);gl.bufferData(gl.ARRAY_BUFFER,P,gl.DYNAMIC_DRAW);gl.bindBuffer(gl.ARRAY_BUFFER,m.n);gl.bufferData(gl.ARRAY_BUFFER,N,gl.DYNAMIC_DRAW);}else{m.p=P;m.n=N;}};
proto.update=function(t){if(this.action?.kind!=='muscleAtlas')return oldUpdate.call(this,t);this.t=t;const e=this.entry,v=this.viewer,p=this.atlasPolicy,T=motionMatrix(e,t);v.transforms=this.anatomy.parts.map(x=>x.type==='bone'&&p.move(x)?T:ident());
 for(const s of (e.dynamic===false?[]:this.atlasSource)){const d=skin(s.P,s.N,this.atlasO,this.atlasI,T);this.atlasUpload(s.i,d.V,d.N);}
 const O=this.atlasO,I=point(T,this.atlasI),mid=mul(add(O,I),.5),dl=norm(sub(O,I)),scale=Math.min(.055,near(O,I)*.4),end=add(I,mul(dl,scale));
 v.axes=[];if(e.markerMode!=='none'&&this.action.showAtlasMarkers!==false){v.axes.push({points:[O,add(O,[0,.012,0])],color:[.16,.43,.82],dots:true,overlay:true},{points:[I,add(I,[0,.012,0])],color:[.88,.47,.13],dots:true,overlay:true});}
 if(e.markerMode!=='none'&&this.action.showAtlasPull!==false){const side=norm(SM.cross(dl,[0,0,1]));v.axes.push({points:[I,end,add(sub(end,mul(dl,scale*.25)),mul(side,scale*.18)),end,sub(sub(end,mul(dl,scale*.25)),mul(side,scale*.18))],color:[.86,.2,.14],overlay:true});}
 this.atlasReadout={id:e.id,origin:O,insertion:I,referenceAngle:e.dynamic===false?0:p.deg*t,joint:p.joint,motion:p.name,geometryType:'source mesh, visual attachment binding',sourceParts:e.sourceIDs};v.invalidate();this.atlasLabels();};
proto.atlasLabels=function(){let wrap=document.getElementById('atlasLabels');if(!wrap)return;wrap.hidden=this.action?.kind!=='muscleAtlas'||this.entry?.markerMode==='none'||this.action.showAtlasMarkers===false;if(wrap.hidden)return;const v=this.viewer;const ensure=(key,text,pos,col)=>{let n=wrap.querySelector('[data-lab="'+key+'"]');if(!n){n=document.createElement('span');n.dataset.lab=key;n.style.setProperty('--tag',col);wrap.appendChild(n);}n.textContent=text;const q=v.project(pos);n.hidden=!q[2];n.style.left=Math.max(5,Math.min(this.canvas.clientWidth-102,q[0]+9))+'px';n.style.top=Math.max(5,Math.min(this.canvas.clientHeight-32,q[1]-14))+'px';};ensure('origin','起点区域',this.atlasReadout.origin,'#2e66c5');ensure('insertion','止点区域',this.atlasReadout.insertion,'#b35b14');};
proto.fitAll=function(){if(this.action?.kind!=='muscleAtlas')return oldFit.call(this);const v=this.viewer,points=[],keep=this.t;for(const t of [0,1]){const T=motionMatrix(this.entry,t);this.anatomy.parts.forEach((p,i)=>{if(!v.visible[i])return;const M=p.type==='bone'&&this.atlasPolicy.move(p)?T:ident();for(const x of [p.bounds[0][0],p.bounds[1][0]])for(const y of [p.bounds[0][1],p.bounds[1][1]])for(const z of [p.bounds[0][2],p.bounds[1][2]])points.push(point(M,[x,y,z]));});}
 const lo=[0,1,2].map(j=>Math.min(...points.map(p=>p[j]))),hi=[0,1,2].map(j=>Math.max(...points.map(p=>p[j])));v.target=lo.map((x,j)=>(x+hi[j])/2);const z=norm([Math.cos(v.pitch)*Math.sin(v.yaw),Math.sin(v.pitch),Math.cos(v.pitch)*Math.cos(v.yaw)]),x=norm(SM.cross([0,1,0],z)),y=SM.cross(z,x),ar=Math.max(.3,this.canvas.clientWidth/Math.max(1,this.canvas.clientHeight)),tv=Math.tan(v.fov/2),th=tv*ar;let dist=.03;for(const p of points){const q=sub(p,v.target),dep=dot(q,z);dist=Math.max(dist,Math.abs(dot(q,x))/th+dep,Math.abs(dot(q,y))/tv+dep);}v.radius=dist*1.10;v.saveView();v.invalidate();this.atlasLabels();};
proto.orient=function(side){if(this.action?.kind!=='muscleAtlas')return oldOrient.call(this,side);this.viewer.orient(side);this.fitAll();};
proto.atlasFocus=function(){if(this.action?.kind!=='muscleAtlas')return;const pts=[];for(const x of this.atlasSource){const d=this.entry.dynamic?skin(x.P,x.N,this.atlasO,this.atlasI,motionMatrix(this.entry,this.t)).V:x.P;for(let i=0;i<d.length;i+=3)pts.push([d[i],d[i+1],d[i+2]]);}const lo=[0,1,2].map(j=>Math.min(...pts.map(p=>p[j]))),hi=[0,1,2].map(j=>Math.max(...pts.map(p=>p[j])));this.viewer.target=lo.map((v,j)=>(v+hi[j])/2);const size=Math.max(...hi.map((v,j)=>v-lo[j]));this.viewer.radius=Math.max(.17,size*2.8);this.viewer.invalidate();};
// labels follow orbit camera on actual draw as well as slider updates
const draw=LabViewer.prototype.draw;LabViewer.prototype.draw=function(...args){const r=draw.apply(this,args);if(this.dataset==='atlas')g.App?.state.model?.atlasLabels();return r;};
g.AtlasMath={policy,motionMatrix,endRegions,skin};
})(window);


(function(g){'use strict';const M=g.LabMath,oldRender=Sim.render,oldDevice=M.deviceTrial,E=Sim.esc;
const clamp=(x,a=0,b=1)=>Math.max(a,Math.min(b,Number(x))),smooth=x=>{x=clamp(x);return x*x*(3-2*x);};
const S=(b,w=560,h=470)=>`<svg viewBox="0 0 ${w} ${h}" xmlns="http://www.w3.org/2000/svg" role="img">${b}</svg>`;
const line=(a,b,c='#406994',w=5,dash='')=>`<line x1="${a[0]}" y1="${a[1]}" x2="${b[0]}" y2="${b[1]}" stroke="${c}" stroke-width="${w}" stroke-linecap="round" ${dash?'stroke-dasharray="'+dash+'"':''}/>`;
const text=(x,y,s,n=17,c='#2c4b6c')=>`<text x="${x}" y="${y}" font-size="${n}" fill="${c}" font-family="system-ui,sans-serif">${E(s)}</text>`;
const dot=(p,r=7,c='#d17857')=>`<circle cx="${p[0]}" cy="${p[1]}" r="${r}" fill="${c}"/>`;
const rect=(x,y,w,h,c='#e7eef6')=>`<rect x="${x}" y="${y}" width="${w}" height="${h}" rx="9" fill="${c}"/>`;
const arrow=(a,b,c='#2c65bd',w=3)=>{const v=[b[0]-a[0],b[1]-a[1]],len=Math.hypot(...v)||1,u=v.map(x=>x/len),h=9;return line(a,b,c,w)+line(b,[b[0]-h*u[0]+h*.5*u[1],b[1]-h*u[1]-h*.5*u[0]],c,w)+line(b,[b[0]-h*u[0]-h*.5*u[1],b[1]-h*u[1]+h*.5*u[0]],c,w);};
const poly=(p,c='#326bba',w=3)=>`<polyline points="${p.map(q=>q.join(',')).join(' ')}" fill="none" stroke="${c}" stroke-width="${w}" stroke-linejoin="round"/>`;
function values(items){return '<div class="values">'+items.map(([l,v])=>`<div><small>${E(l)}</small><b>${E(v)}</b></div>`).join('')+'</div>';}
M.deviceTrial=function(mode,t,s={}){if(mode==='isokinetic')return oldDevice(mode,t,s);if(mode!=='mvc'&&mode!=='mvcBack'&&mode!=='mvcKnee')return oldDevice(mode,t,s);const u=clamp(t),ready=!!s.deviceStarted,back=mode!=='mvcKnee',duration=7,scale=(s.mvcEffort==='submax'?.6:1)*(s.deviceSide==='left'?.95:1),max=(back?430:285)*scale;
 const forceAt=x=>max*(x<.12?0:x<.42?smooth((x-.12)/.30):x<.72?1-.025*Math.sin(Math.PI*(x-.42)/.3)**2:x<.92?1-smooth((x-.72)/.20):0);
 const force=ready?forceAt(u):0,peak=ready?Math.max(...Array.from({length:151},(_,i)=>forceAt(u*i/150))):0,arm=back?.32:.34;
 return {u,mvc:true,back,ready,speed:0,duration,time:ready?duration*u:0,angle:back?0:60,velocity:0,userTorque:force*arm,deviceTorque:-force*arm,force,peak,forceAt,maximum:max,reverse:false,phase:s.deviceStopped?'试次停止':!ready?'准备固定体位':u<.12?'准备':u<.42?'逐渐用力':u<.72?'维持峰值':u<.92?'缓慢放松':'试次结束',condition:s.mvcEffort==='submax'?'次最大努力对照':'最大努力指令',jointLabel:back?'躯干相对中立位':'膝屈曲角'};
};
function mvcFigure(d){const back=d.back;let b=text(28,35,back?'躯干伸肌 · 固定体位等长用力':'坐位膝伸肌 · 固定体位等长用力',21);b+=line([34,433],[529,433],'#dce5ef',2)+rect(166,292,158,18,'#a4bbd1')+line([181,314],[181,428],'#819bb4',10)+line([310,314],[310,428],'#819bb4',10)+dot([242,95],25,'#c9a384');b+=line([237,132],[219,284],'#7795b2',44)+line([246,143],[282,210],'#acc0d1',14)+line([282,210],[315,252],'#acc0d1',12)+line([219,286],[320,312],'#87a2be',23);
 if(back){b+=line([320,312],[321,419],'#87a2be',18)+line([321,419],[367,419],'#51718e',13)+rect(136,148,24,107,'#d3e2f1')+line([149,255],[149,397],'#7898b6',10)+line([151,198],[199,198],'#7999b9',16)+rect(185,160,21,76,'#bdd2e5');b+=line([202,270],[245,279],'#3b607e',13)+line([290,310],[342,325],'#3b607e',12);b+=text(40,118,'传感器 / 背垫',16)+line([96,128],[174,171],'#7e96af',1.5);b+=text(356,277,'骨盆固定',17)+line([349,266],[236,271],'#7e96af',1.5);b+=text(365,367,'下肢固定',17)+line([359,350],[332,328],'#7e96af',1.5);if(d.force>0){const len=18+d.force/d.maximum*54;b+=arrow([193,182],[193-len,182],'#386fc2',4)+arrow([139,219],[139+len,219],'#cb7958',4);}b+=text(39,465,'尝试向后伸展，但固定装置不允许躯干继续后移。',15);
 }else{const K=[320,312],q=Math.atan2(26,101)+Math.PI/3,A=[320+110*Math.cos(q),312+110*Math.sin(q)];b+=line(K,A,'#86a3c0',20)+line(A,[A[0]+42,A[1]+7],'#567b9d',13)+rect(414,333,52,81,'#d5e3f1')+line([442,414],[442,430],'#7695b2',11)+line([410,430],[499,430],'#7695b2',10)+line([340,384],[440,370],'#6f91b1',8)+line([330,380],[351,391],'#345e85',13)+line([279,300],[291,326],'#345e85',13)+rect(155,168,21,115,'#d4e2ef');b+=text(357,237,'尝试伸膝',19)+text(358,263,'膝角保持60°',16);b+=text(54,350,'大腿与骨盆固定',16)+line([191,347],[282,316],'#8199ae',1.5);b+=text(368,459,'小腿前方固定传感器',15);if(d.force>0){const len=18+d.force/d.maximum*50;b+=arrow([340,384],[340+len,363],'#386fc2',4)+arrow([370,400],[370-len,421],'#cb7958',4);}}
 return S(b,560,480);}
function forceChart(d,t){const pts=Array.from({length:101},(_,i)=>[45+422*i/100,175-d.forceAt(i/100)/(d.maximum*1.15)*140]);let b=text(28,25,'传感器力—时间',18)+line([45,175],[475,175],'#8fa5bc',1.5)+line([45,47],[45,175],'#8fa5bc',1.5)+poly(pts,'#d4dfeb',2)+text(45,202,'0 s',14)+text(430,202,'7 s',14)+text(53,52,'N',13);if(d.ready){const n=Math.max(2,Math.ceil(t*100));b+=poly(Array.from({length:n},(_,i)=>{const u=t*i/(n-1);return [45+422*u,175-d.forceAt(u)/(d.maximum*1.15)*140];}),'#3470c4',3)+dot([45+422*t,175-d.force/(d.maximum*1.15)*140],5,'#bf6c48');}return S(b,520,220);}
function isoFigure(d){const O=[294,266],rad=d.angle*Math.PI/180,H=[O[0]+140*Math.sin(rad),O[1]+140*Math.cos(rad)],T=[O[0]+172*Math.sin(rad),O[1]+172*Math.cos(rad)];let b=text(24,35,'等速测力系统 · 肘关节',21)+line([27,435],[520,435],'#dae5ef',2);b+=rect(156,314,176,16,'#b7cbe0')+line([177,331],[177,429],'#7c99b2',9)+line([319,331],[319,429],'#7c99b2',9)+rect(132,170,22,143,'#c9dbed')+dot([232,102],24,'#cba583')+line([227,138],[205,302],'#7e9fbd',43)+line([208,301],[317,337],'#95adc3',21)+line([317,337],[325,423],'#95adc3',19)+line([325,423],[371,423],'#557997',12);b+=line([243,151],[294,174],'#a3bbce',15)+line([294,174],O,'#a3bbce',18)+line(O,H,'#a3bbce',17)+line(H,T,'#9bb2c5',12)+line([196,222],[238,228],'#45647f',12)+line([192,290],[230,300],'#45647f',12);b+=rect(450,256,62,80,'#d7e4f1')+line([481,336],[481,428],'#7c9db9',13)+line([449,430],[514,430],'#7c9db9',9)+line([471,275],O,'#a8c0d9',12)+dot(O,22,'#dce9f5')+dot(O,14,'white')+line(O,H,'#346cad',6)+dot(O,6,'#3474c7');const dx=Math.cos(rad)*12,dy=-Math.sin(rad)*12;b+=line([H[0]-dx,H[1]-dy],[H[0]+dx,H[1]+dy],'#3f6386',10);b+=text(33,383,'转轴与肘部对齐',15)+line([176,375],O,'#7f9db4',1.5);if(d.ready){const tang=[Math.cos(rad),-Math.sin(rad)],l=18+Math.abs(d.userTorque)/45*46;b+=arrow([H[0]-8,H[1]-8],[H[0]-8+tang[0]*l,H[1]-8+tang[1]*l],'#3376c6',4)+arrow([H[0]+8,H[1]+8],[H[0]+8-tang[0]*l,H[1]+8-tang[1]*l],'#c47858',4);}b+=text(30,465,'蓝：受试者力矩　橙：装置反向力矩',15);return S(b,560,485);}
function device(a,t,s){const d=M.deviceTrial(a.motion,t,s);let right=`<h2>${d.mvc?(d.back?'躯干伸肌最大等长收缩':'坐位膝伸肌最大等长收缩'):'等速肌力测试'}</h2><p class="atlas-note">${E(d.phase)}</p>`;
 if(d.mvc){right+=values([[d.jointLabel,d.angle.toFixed(0)+'°'],['试次时间',d.time.toFixed(2)+' s'],['当前传感器力',d.force.toFixed(1)+' N'],['已观察峰值',d.peak.toFixed(1)+' N']])+forceChart(d,t)+`<p>体位固定 → 逐渐用力 → 保持 → 放松。</p><p class="subtle-limit">传感器外力不等于某一块肌肉的力。${E(d.condition)}；读数为模拟案例。</p>`;}else{right+=values([['设定角速度',d.speed+'°/s'],['关节角度',d.angle.toFixed(1)+'°'],['当前角速度',Math.abs(d.velocity).toFixed(0)+'°/s'],['案例时间',d.time.toFixed(3)+' s']])+`<div class="data-tile"><h3>阻力随发力顺应变化</h3><div class="bar-row"><span>受试者力矩</span><div class="bar-track"><i style="width:${100*d.userTorque/45}%"></i></div><b>${d.userTorque.toFixed(1)}</b></div><div class="bar-row"><span>装置反力矩</span><div class="bar-track"><i style="width:${100*Math.abs(d.deviceTorque)/45}%"></i></div><b>${Math.abs(d.deviceTorque).toFixed(1)}</b></div><p class="subtle-limit">单位N·m。改变发力，反力矩改变；恒速工作段角速度不变。</p></div><p>工作方向：${d.reverse?'离心控制伸肘':'向心屈肘'}。上方设置的是关节角速度；下方倍速只改变回放快慢。</p><p class="subtle-limit">虚拟装置仅显示恒速工作段；未求解真实设备的加减速与重力补偿。</p>`;}
 return `<div class="instrument88"><div>${d.mvc?mvcFigure(d):isoFigure(d)}</div><div>${right}</div></div>`;}
M.vitals88=function(mode,t,s={}){const age=clamp(s.caseAge??22,16,80),restHR=clamp(s.caseRestHR??72,40,110),maxHR=220-age,pct=clamp(s.targetPct??65,40,90),intensity=clamp(s.caseIntensity??.60,.25,.8),sec=mode==='rest'?t*60:t*180,gain=1-Math.exp(-sec/48);let hr=restHR,sbp=118,dbp=76,rr=14,stage='安静';
 if(mode==='walk'){hr=restHR+(maxHR-restHR)*intensity*gain;sbp=118+80*intensity*(1-Math.exp(-sec/52));dbp=76-4*intensity*(1-Math.exp(-sec/65));rr=14+19*intensity*gain;stage='运动';}
 if(mode==='recovery'){const start=M.vitals88('walk',1,s),h=Math.exp(-sec/55);hr=restHR+(start.hr-restHR)*h;sbp=118+(start.sbp-118)*Math.exp(-sec/66);dbp=76+(start.dbp-76)*Math.exp(-sec/85);rr=14+(start.rr-14)*h;stage='恢复';}
 return {time:sec,age,restHR,maxHR,target:maxHR*pct/100,pct,hr,sbp,dbp,rr,stage,relativeHR:100*hr/maxHR,synthetic:true};};
function trendChart(mode,t,s){const ps=Array.from({length:81},(_,i)=>M.vitals88(mode,i/80,s)),map=(v,k,i)=>[38+i/80*454,177-(v[k]-40)/160*135];let b=text(25,25,'同步趋势',18)+line([38,177],[494,177],'#c4d4e5',1.5)+line([38,43],[38,177],'#c4d4e5',1.5);[['sbp','#b66b47','收缩压'],['dbp','#36a193','舒张压'],['hr','#306dcc','心率']].forEach(([k,c,l],i)=>{b+=poly(ps.map((v,j)=>map(v,k,j)),c,2)+text(45+i*150,218,l+(k==='hr'?' 次/分':' mmHg'),13,c);const v=M.vitals88(mode,t,s);b+=dot([38+t*454,177-(v[k]-40)/160*135],4,c);});b+=line([38+t*454,39],[38+t*454,177],'#8aa2bd',1,'4 4')+text(38,196,'0',13)+text(428,196,(mode==='rest'?60:180)+' s',13);return S(b,540,235);}
function treadmillFigure(t,s){const move=s.safety?.state==='running';const d=M.gait(move?(t*6)%1:0);const offset=d.stride*d.t;const P=p=>[225+(p[0]-offset)*192,429-p[1]*192];let b=rect(36,438,440,32,'#627e83')+line([449,437],[449,157],'#587e99',8)+rect(401,127,114,44,'#e0eaf4')+dot([478,150],9,'#bc6a4e');for(let i=0;i<10;i++){const x=48+(i*42+(move?t*360:0))%415;b+=line([x,447],[x+9,459],'#b6c8cb',3);}for(const sd of ['left','right']){const l=d.legs[sd];b+=poly([P(l.hip),P(l.knee),P(l.ankle)],sd==='right'?'#2e6aaa':'#a6bbd1',16)+line(P(l.heel),P(l.toe),'#527c9b',11);}b+=line(P(d.hip),P(d.shoulder),'#578aaf',28)+dot(P(d.head),23,'#c69e7d');for(const side of ['right','left']){const l=d.arms[side];b+=poly([P(l.s),P(l.elbow),P(l.wrist)],'#90adc3',10);}b+=line(P(d.shoulder),[478,150],'#bf7956',2,'5 5')+text(380,104,'急停',17,'#aa593e');return S(b,550,488);}
function cardio(a,t,s){const ready=s.safety?.state!=='idle';const u=ready?t:0,d=M.vitals88(a.motion,u,s),st=s.safety?.state;return `<header class="vital88-head"><h2>心肺与运动 · ${d.stage}观察</h2><p>${st==='running'?'本次正在运行':st==='emergency'?'运动已终止，请按风险情景处理':'等待启动或本次已停止'}　<span class="subtle-limit">模拟读数</span></p></header><div class="instrument88"><div>${treadmillFigure(u,s)}</div><div>${values([['心率',Math.round(d.hr)+' 次/分'],['呼吸频率',Math.round(d.rr)+' 次/分'],['收缩压',Math.round(d.sbp)+' mmHg'],['舒张压',Math.round(d.dbp)+' mmHg']])}${trendChart(a.motion,u,s)}<p class="vital88-note">动态运动中，收缩压通常随负荷上升；舒张压多保持相近或略降。恢复期逐渐向基线回落，不把两者画成相同走势。</p></div></div><div class="calculation88">案例年龄 ${d.age}岁 · 估算最大心率 <b>${d.maxHR}次/分</b> · 本题${d.pct}%最大心率 = <b>${d.target.toFixed(1)}次/分</b><br><small>当前约${d.relativeHR.toFixed(1)}%估算最大心率。此处为计算练习，不直接生成个体运动处方。</small></div>`;}
M.metabolic88=function(t,s={}){const age=clamp(s.caseAge??22,16,80),rest=clamp(s.caseRestHR??72,40,110),max=220-age;const p=clamp(t),contT=600*p,sprintT=10*p;const a=rest+(max-rest)*.58*(1-Math.exp(-contT/52)),b=rest+(max-rest)*.88*(1-Math.exp(-sprintT/8));return{age,rest,max,t:p,a:{time:contT,hr:a,rr:14+12*(1-Math.exp(-contT/55)),rpe:13,phase:contT<45?'起始调整':'持续稳定负荷'},b:{time:sprintT,hr:b,rr:14+22*(1-Math.exp(-sprintT/10)),rpe:19,phase:sprintT<4?'快速输出':'短时维持'},target:max*clamp(s.targetPct??65,40,90)/100};};
function cyclist(p,color,fast){const q=p*2*Math.PI*(fast?14:6),hip=[159,216],pedals=[[257+31*Math.cos(q),338+31*Math.sin(q)],[257-31*Math.cos(q),338-31*Math.sin(q)]];let b=line([30,443],[510,443],'#dbe5ef',2)+line([167,229],[124,399],'#93adc3',7)+line([144,229],[192,229],'#617f9f',9)+line([298,150],[339,164],'#7595b2',8)+line([335,164],[309,340],'#94afc7',7)+line([255,338],[109,434],'#92b0c9',9)+line([255,338],[354,434],'#92b0c9',9)+dot([255,338],42,'#e0eaf4');for(let i=1;i>=0;i--){const K=(()=>{const D=[pedals[i][0]-hip[0],pedals[i][1]-hip[1]],r=Math.hypot(...D),L=105,along=r/2,h=Math.sqrt(Math.max(0,L*L-along*along));return [hip[0]+D[0]*.5+D[1]/r*h,hip[1]+D[1]*.5-D[0]/r*h];})();b+=poly([hip,K,pedals[i]],i?'#afc5d7':color,15)+line([pedals[i][0]-15,pedals[i][1]],[pedals[i][0]+15,pedals[i][1]],'#527c9a',8);}b+=line(hip,[199,118],color,28)+dot([216,70],23,'#c9a17e')+poly([[198,120],[265,147],[321,163]],'#9db8d0',13)+text(45,481,fast?'短时高阻力快速踩踏':'中等阻力持续踩踏',17);return S(b,540,505);}
function metabolic(a,t,s){const d=M.metabolic88(t,s);function col(v,fast){return `<section><h3>${fast?'短时高强度任务':'中等强度持续任务'}</h3><p>${fast?'强度高、时间短，快速供能需求突出。':'稳定负荷持续骑行，氧化供能贡献逐渐突出。'}</p>${cyclist(t,fast?'#cd8863':'#3c79bf',fast)}<div class="readings"><span>时间 ${v.time.toFixed(1)} s</span><span>心率 ${v.hr.toFixed(0)}次/分</span><span>${(v.hr/d.max*100).toFixed(1)}%HRmax</span></div><p>呼吸 ${v.rr.toFixed(0)}次/分 · 主观用力示例RPE ${v.rpe}/20</p><p><b>${fast?(v.time<4?'磷酸原系统快速参与':'糖酵解贡献增加，氧化供能仍参与'):'持续阶段以有氧代谢为主，其他供能系统并未关闭'}</b></p><p class="subtle-limit">${fast?'心率有上升滞后，短时任务未必达到最高心率。':'任务强度和时间是理解主导供能的条件。'}</p></section>`;}
 return `<header class="vital88-head"><h2>有氧与无氧 · 比较强度和持续时间</h2><p>两种任务共享回放进度，但案例时长不同。心率为同一模型生成的练习数值。</p></header><div class="metab88">${col(d.a,false)}${col(d.b,true)}</div><div class="calculation88">案例年龄 ${d.age}岁：估算最大心率 <b>${d.max}次/分</b>；设定比例对应目标 <b>${d.target.toFixed(1)}次/分</b>。<br><small>最大心率估算采用220－年龄；HRmax百分比与心率储备百分比不是同一计算。心率不能单独测出乳酸阈或判定无氧供能比例。</small></div>`;}
Sim.render=function(a,t,s){if(a.kind==='dynamometer')return device(a,t,s);if(a.kind==='cardio')return cardio(a,t,s);if(a.kind==='metabolism')return metabolic(a,t,s);return oldRender(a,t,s);};
g.Revision88Views={device,cardio,metabolic};
})(window);


/* Classroom entry + per-learner local records. A teacher-name check is not
   identity authentication. Optional collection uses the included same-origin server. */
(function(g){'use strict';const $=id=>document.getElementById(id),E=Sim.esc;let service=false,csrf='',access=false,sessionStart=0,activeView='login',popping=false;
const blank=()=>({format:'kinesiology-observation-8.8',version:'8.8.0',profile:null,observations:[],explanations:[],events:[]});
const norm=s=>s.normalize('NFKC').trim();
function validate(name,id,teacher){name=norm(name);id=norm(id);teacher=norm(teacher);if(!name||name.length>60||/[<>\u0000-\u001f]/.test(name))return '请填写有效姓名（1—60个字符）。';if(!/^[A-Za-z0-9_-]{1,32}$/.test(id))return '请填写有效学号（字母、数字、下划线或短横线，最多32位）。';if(teacher!=='张积众')return '任课教师姓名不匹配，请核对后重新填写。';return {name,studentId:id,teacher};}
function key(){const p=App.state.profile;return p?'kine_v88_student_'+encodeURIComponent(p.teacher+'|'+p.studentId):null;}
function push(view){if(popping)return;const hs={kine:true,view,module:App.state.module?.id,action:App.state.action?.id};const s=history.state;if(s?.view===view&&s?.module===hs.module)return;history.pushState(hs,'','#'+view);}
function hidePortal(){['loginPage','welcomePage','finishPage'].forEach(id=>$(id).hidden=true);}
function onView(view){if(!access){loginView();return false;}hidePortal();activeView=view;$('profileBar').hidden=false;document.querySelector('.top nav').hidden=false;if(view==='projects')$('personalHello').textContent=App.state.profile.name+'同学，选择今天的学习项目。';push(view);return true;}
function loginView(){App.pause();if(App.state.clinic)App.state.clinic.active=false;access=false;activeView='login';hidePortal();$('loginPage').hidden=false;$('home').hidden=true;$('lesson').hidden=true;$('profileBar').hidden=true;document.querySelector('.top nav').hidden=true;document.body.classList.remove('focus-mode');}
function welcome(){hidePortal();$('welcomePage').hidden=false;$('home').hidden=true;$('lesson').hidden=true;activeView='welcome';$('welcomeName').textContent=App.state.profile.name+'同学';$('welcomeTeacher').textContent='任课教师：'+App.state.profile.teacher;$('profileName').textContent=App.state.profile.name+' · '+App.state.profile.studentId;$('profileBar').hidden=false;document.querySelector('.top nav').hidden=false;push('welcome');window.scrollTo(0,0);}
async function login(e){e.preventDefault();$('loginError').textContent='';const p=validate($('studentName').value,$('studentId').value,$('teacherName').value);if(typeof p==='string'){$('loginError').textContent=p;return;}
 const b=$('enterBtn');b.disabled=true;
 try{if(service){const r=await fetch('/api/classroom/login',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(p)});const d=await r.json();if(!r.ok)throw Error(d.error||'登记失败');csrf=d.csrf;}
 $('continueBtn')&&($('continueBtn').hidden=true);App.state.profile=p;App.state.record=blank();App.state.record.profile=p;
 try{const s=JSON.parse(localStorage.getItem(key())||'null');if(s){App.validateRecord(s);if(s.profile?.studentId!==p.studentId||s.profile?.teacher!==p.teacher)throw Error('记录身份不匹配');App.state.record=s;App.state.record.profile=p;}}catch(err){console.info('未载入旧记录：',err.message);}
 access=true;sessionStart=Date.now();App.state.record.events.push({type:'student_entry',at:new Date().toISOString()});App.save();welcome();
 }catch(err){$('loginError').textContent='未完成登记：'+err.message;}finally{b.disabled=false;}}
function finish(){if(!access)return loginView();App.pause();if(App.state.clinic)App.state.clinic.active=false;App.state.record.events.push({type:'session_end',at:new Date().toISOString(),seconds:Math.max(0,Math.round((Date.now()-sessionStart)/1000))});App.save();hidePortal();$('home').hidden=true;$('lesson').hidden=true;$('finishPage').hidden=false;activeView='finish';$('finishName').textContent=App.state.profile.name+'同学，感谢你的使用。';const n=new Set(App.state.record.events.filter(e=>e.type==='open').map(e=>e.module)).size;$('finishStats').textContent=`已访问 ${n} 个项目 · 已保存 ${App.state.record.observations.length} 个观察时刻 · ${App.state.record.explanations.length} 条学习解释`;$('submitRecord').hidden=!service;$('collectionStatus').textContent=service?'可将本次记录提交给任课教师。':'当前为本地使用；请导出学习记录后交给老师，信息不会自动发送。';push('finish');window.scrollTo(0,0);}
async function submit(){const b=$('submitRecord');b.disabled=true;try{const r=await fetch('/api/classroom/report',{method:'POST',headers:{'Content-Type':'application/json','X-Classroom-CSRF':csrf},body:JSON.stringify(App.state.record)});const d=await r.json();if(!r.ok)throw Error(d.error||'提交失败');$('collectionStatus').textContent='已提交。本次接收编号：'+d.receipt;}catch(e){$('collectionStatus').textContent='未提交成功：'+e.message+'。仍可导出本地备份。';}finally{b.disabled=false;}}
function logout(){App.pause();App.save();App.state.record=blank();App.state.profile=null;App.state.module=null;App.state.action=null;csrf='';loginView();$('loginForm').reset();$('loginError').textContent='';if(service)fetch('/api/classroom/logout',{method:'POST'}).catch(()=>{});history.replaceState({kine:true,view:'login'},'','#login');}
function init(){document.querySelector('.top nav').hidden=true;$('loginForm').addEventListener('submit',login);$('welcomeStart').onclick=()=>App.home();$('welcomeMuscle').onclick=()=>App.openLesson(201);$('endLearning').onclick=finish;$('endFromLesson').onclick=finish;$('backToProjects').onclick=()=>App.home();$('finishExport').onclick=()=>App.actionCommand('exportHTML');$('finishBackup').onclick=()=>App.actionCommand('exportJSON');$('submitRecord').onclick=submit;$('switchStudent').onclick=logout;$('finishLogout').onclick=logout;loginView();history.replaceState({kine:true,view:'login'},'','#login');
 window.addEventListener('popstate',e=>{popping=true;try{if(!access)return loginView();const s=e.state;if(s?.view==='lesson'&&s.module){App.openLesson(s.module);if(s.action)App.selectAction(s.action);}else if(s?.view==='welcome')welcome();else if(s?.view==='finish')finish();else App.home();}finally{popping=false;}});
 if(window.DXX_ENABLE_CLASSROOM_SERVICE===true&&(location.protocol==='http:'||location.protocol==='https:'))fetch('/api/classroom/status',{cache:'no-store'}).then(r=>r.ok?r.json():null).then(d=>{service=!!d?.classroom;$('entryPrivacy').textContent=service?'姓名、学号与学习记录将用于本课程教学登记，由本校部署服务保存。教师姓名仅为课堂入口校验。':'姓名、学号与学习记录仅保存在此浏览器。结束后可导出交给老师。教师姓名仅为课堂入口校验。';}).catch(()=>{});
}
g.Portal={init,onView,key,validate,finish,logout,welcome,get access(){return access;},get service(){return service;}};
})(window);


/* Observation controller. No automatic clinical or patient competency grading. */
(function(){'use strict';
const payload=window.__DXX_DATA__,course=payload.curriculum,$=id=>document.getElementById(id),E=Sim.esc,MM=LabMath;
window.LabImages=payload.images||{};window.LabRegistration=payload.registration||{};
window.MuscleCatalogue=payload.muscleCatalogue;
const state={atlasQuery:'',atlasGroup:'all',atlasEmpty:false,profile:null,module:null,action:null,mode:'main',side:'right',t:0,playing:false,dir:1,lastTime:0,raf:0,model:null,clinic:null,role:'therapist',clinicPhase:'intro',consented:false,record:{format:'kinesiology-observation-8.8',version:'8.8.0',observations:[],explanations:[],events:[]},params:{deviceStarted:false,deviceStopped:false,deviceChecks:[false,false,false],deviceEffort:.7,deviceProfile:'rise',deviceSide:'right',mvcEffort:'max',supportSide:'right',caseAge:22,caseRestHR:72,targetPct:65,caseIntensity:.60,showAtlasMarkers:true,showAtlasPull:true,isoSpeed:60,isoDirection:'concentric',effort:.6,failureMode:'loading',muscleSide:'right',archView:'medial',showPaths:false,load:30,distance:.32,material:'A',loadCycle:false,mass:12,takeoff:2.6,quality:'valid',airDistance:55,airLevel:2,airSeconds:8,chain:'free',stiffness:1,morphology:'typical',scap:'elevation',safety:{checked:[false,false,false],state:'idle',message:'先核对安全准备，再启动运动。'}},saved:false};
let toastTimer,selectedName='';
function toast(msg){$('toast').textContent=msg;$('toast').hidden=false;clearTimeout(toastTimer);toastTimer=setTimeout(()=>$('toast').hidden=true,4500);}
function readRecord(){state.saved=false;}
function save(){const k=window.Portal?.key();if(!k)return;state.record.profile=state.profile;try{localStorage.setItem(k,JSON.stringify(state.record));state.saved=true;}catch{state.saved=false;}$('storageHint').textContent=state.saved?'本浏览器按学号分别保存；结束后可导出或提交。':'暂不能持久保存，请在结束前导出记录。';}
function event(type,data={}){if(state.record.events.length<3000)state.record.events.push({at:new Date().toISOString(),module:state.module?.id,action:state.action?.id,type,...data});save();}
function modal(title,content){$('modalTitle').textContent=title;$('modalBody').innerHTML=content;$('modal').showModal();}
function download(name,type,content){const u=URL.createObjectURL(new Blob([content],{type})),a=document.createElement('a');a.href=u;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(u),1500);}
function decode(name){const a=payload[name];if(!a.raw){const bytes=Uint8Array.from(atob(a.b64),c=>c.charCodeAt(0));a.raw=bytes.buffer;}return a;}
function getModel(){if(!state.model){state.model=new ModelStage($('modelCanvas'),{native:decode('native'),anatomy:decode('anatomy'),motions:payload.motions},p=>{selectedName=p?.label||p?.english||p?.name||'';if(selectedName){toast(selectedName);if(state.action?.kind==='muscleAtlas')$('atlasSelected').textContent='当前点选：'+selectedName;}});}return state.model;}
function getClinic(){if(!state.clinic){const a=decode('actor');state.clinic=new TeachingClinic($('clinicCanvas'),a.model,a.raw,{error:toast});}return state.clinic;}
function roleMotion(a){const mo={shoulderAbd:'shoulder',shoulderFlex:'shoulderFlex',elbow:'elbow',knee:'knee',ankle:'ankle',trunkFlex:'spine'};return a.patient?mo[a.motion]||null:null;}
function assetSupported(a){return ['native','contraction'].includes(a.kind)&&['upper','shoulder','hand','lower','knee','foot'].includes(a.region||payload.motions.find(m=>m.id===a.motion)?.region);}
function nativeActive(){return !state.atlasEmpty&&(state.mode==='anatomy'||state.mode==='main'&&['native','fsu','pelvis3d','muscleAtlas'].includes(state.action?.kind)&&!state.action?.referenceOnly);}
function liveAction(){return {...state.action,side:state.side,supportSide:state.params.supportSide,showMusclePaths:state.params.showPaths,showAtlasMarkers:state.params.showAtlasMarkers,showAtlasPull:state.params.showAtlasPull};}
function home(){if(!window.Portal?.onView('projects'))return;pause();document.body.classList.remove('focus-mode');$('focusBtn').textContent='放大观察';$('focusBtn').setAttribute('aria-pressed','false');state.clinic&&(state.clinic.active=false);state.module=null;state.action=null;$('lesson').hidden=true;$('home').hidden=false;cards();window.scrollTo(0,0);}
let courseFilter='all';
const categories={201:'muscles',101:'basic',1:'basic',2:'basic',3:'basic',4:'physio',5:'physio',6:'joints',7:'joints',8:'joints',9:'joints',10:'joints',11:'gait',13:'gait'};
const categoryLabels={muscles:'肌肉学习',basic:'运动基础',physio:'心肺与呼吸',joints:'关节运动',gait:'足弓与步态'};
function cards(){const q=$('courseSearch').value.trim().toLowerCase();const list=course.modules.filter(m=>(courseFilter==='all'||categories[m.id]===courseFilter)&&(!q||(m.title+' '+m.actions.map(a=>a.title).join(' ')).toLowerCase().includes(q)));
 $('courseGrid').innerHTML=list.map(m=>`<button class="course" data-module="${m.id}"><span class="n">${(m.displayCode||String(m.id).padStart(2,'0'))}</span><span class="course-tag">${categoryLabels[categories[m.id]]}</span><h2>${E(m.title)}</h2><p>${E(({201:MuscleCatalogue.filter(e=>!e.referenceOnly).length+'项源模型 · 全身肌肉目录 · 起止与功能',101:'远端约束 · 腿屈伸与深蹲对照',1:'部位 · 动作 · 轴面方向',2:'动态力臂 · 骨受载与形变',3:'收缩形式 · 测力任务 · 肌肉能力',4:'跑步机 · 主导供能 · 安全观察',5:'胸腹运动 · 虚拟火焰对照',6:'肩带 · 肩 · 肘 · 前臂',7:'腕部与手指的精细运动',8:'整体运动 · FSU六自由度',9:'骨盆定向 · 髋膝联动',10:'踝背伸与跖屈 · 足内外翻',11:'弹性承重 · 动态足印与COP',13:'左右支撑 · 质心 · 步态特征'})[m.id])}</p><span class="small">${m.actions.length}项观察${state.record.observations.some(o=>o.module===m.id)?' · 已有本地记录':''}</span><span class="arrow">↗</span></button>`).join('');$('emptySearch').hidden=!!list.length;
 $('moduleCount').textContent=course.modules.length;$('actionCount').textContent=course.modules.reduce((n,m)=>n+m.actions.length,0);
}
function continueLast(){try{const last=JSON.parse(localStorage.getItem(Portal.key()+'_last')||'null');if(last?.module){openLesson(last.module);if(state.module?.actions.some(a=>a.id===last.action))selectAction(last.action);}}catch{}}
function remember(){try{localStorage.setItem(Portal.key()+'_last',JSON.stringify({module:state.module?.id,action:state.action?.id}));$('continueBtn').hidden=false;$('continueBtn').onclick=continueLast;}catch{}}

function openLesson(id){if(!window.Portal?.access)return;pause();const m=course.modules.find(x=>x.id===Number(id));if(!m)return;state.module=m;state.atlasQuery='';state.atlasGroup='all';state.atlasEmpty=false;Portal.onView('lesson');$('home').hidden=true;$('lesson').hidden=false;$('lessonNumber').textContent=`${[101,201].includes(m.id)?'独立专题':'实验'} ${(m.displayCode||String(m.id).padStart(2,'0'))} · 操作观察`;$('lessonTitle').textContent=m.title;
 $('actionSelect').innerHTML=m.id===201?[...new Set(m.actions.map(a=>a.group))].map(group=>`<optgroup label="${E(group)}">${m.actions.filter(a=>a.group===group).map(a=>`<option value="${a.id}">${E(a.title)}</option>`).join('')}</optgroup>`).join(''):m.actions.map(a=>`<option value="${a.id}">${E(a.title)}</option>`).join('');selectAction(m.id===201?m.actions.find(a=>a.muscleID==='supra').id:m.actions[0].id);window.scrollTo(0,0);}
function selectAction(id){if(!window.Portal?.access)return;pause();if(state.clinic)state.clinic.active=false;const a=state.module?.actions.find(x=>x.id===id);if(!a)return;state.atlasEmpty=false;state.action=a;state.t=a.kind==='fsu'?.5:a.default??0;state.mode='main';state.role='therapist';state.clinicPhase='intro';state.consented=false;state.side='right';state.params.muscleSide='right';selectedName='';$('atlasDetails').hidden=a.kind!=='muscleAtlas';$('atlasLabels').hidden=a.kind!=='muscleAtlas'||atlasEntry()?.markerMode==='none';
 $('actionSelect').disabled=false;$('actionSelect').value=a.id;$('command').textContent=a.command;$('recallBody').hidden=true;$('muscleAnswer').hidden=true;$('answer').value='';$('recall').hidden=!MuscleTeaching[a.motion];
 if(a.kind==='cardio')state.params.safety={checked:[false,false,false],state:'idle',message:'先核对安全准备，再启动运动。'};
 if(a.kind==='failure')state.params.failureMode='loading';
 if(a.kind==='candle')state.params.airDistance=a.motion==='duration'?30:55;
 if(a.kind==='dynamometer'){state.params.deviceStarted=false;state.params.deviceStopped=false;state.params.deviceChecks=[false,false,false];}
 if(a.kind==='muscleAtlas')renderAtlasDetails();if(a.kind==='dynamometer')$('speed').value='.25';else $('speed').value='1';renderModes();configure();if(a.kind==='muscleAtlas')syncAtlasOptions();notes();event('open');remember();}
function renderModes(){const a=state.action;let bs=[['main',a.kind==='muscleAtlas'?'肌肉位置与动作':a.kind==='contraction'?'肌肉收缩':a.kind==='native'?'动作模型':a.kind==='pelvis3d'?'三维骨盆':a.kind==='dynamometer'?'虚拟测力':a.kind==='fsu'?'FSU六自由度':'操作观察']];if(assetSupported(a))bs.push(['anatomy','局部骨肌 · 静态']);if(roleMotion(a))bs.push(['clinic','检查与模仿']);
 $('modes').innerHTML=bs.map(([id,s])=>`<button data-mode="${id}" aria-pressed="${state.mode===id}">${s}</button>`).join('');}
function switchMode(mode){pause();state.mode=mode;state.clinic&&(state.clinic.active=false);renderModes();configure();}
function configure(){const a=state.action,na=nativeActive();$('diagram').innerHTML='';$('modelCanvas').hidden=!na;$('clinicCanvas').hidden=state.mode!=='clinic';$('diagram').hidden=na||state.mode==='clinic';$('actorLabels').hidden=state.mode!=='clinic';$('viewport').classList.toggle('diagramMode',!na&&state.mode!=='clinic');$('viewButtons').hidden=!na&&state.mode!=='clinic';$('poseLabel').hidden=state.mode!=='clinic'&&a.kind!=='pelvis3d';$('clinicBar').hidden=state.mode!=='clinic';
 $('transport').hidden=state.mode==='anatomy'||a.static===true||a.kind==='concept';$('roundtrip').disabled=['contraction','cardio','capacity','ssc','gait','airflow','dynamometer','candle','metabolism','footprints','failure'].includes(a.kind);if($('roundtrip').disabled)$('roundtrip').checked=false;
 $('stageBadge').textContent=state.mode==='anatomy'?'来源网格 · 静态定位':state.mode==='clinic'?'双角色 · 给定动作':a.kind==='native'||a.kind==='pelvis3d'?'源骨网格 · 给定运动':a.kind==='fsu'?'独立节段分量示意':a.kind==='muscleAtlas'?(atlasEntry().referenceOnly?'解剖资料 · 无替代模型':atlasSide(atlasEntry())+'源网格 · '+(atlasEntry().dynamic?'附着动作':'结构定位')):a.kind==='contraction'?'肌腹机制示意 · 非肌纤维实测':a.kind==='lever'?'同坐标骨骼投影 · 准静态力学':'教学模型 / 合成案例';
 try{if(na){getModel().set(liveAction(),state.t,state.mode==='anatomy'?'anatomy':'native');}
 else if(state.mode==='clinic'){
  const c=getClinic(),motion=roleMotion(a),seated=['elbow','ankle'].includes(motion);c.setCase({motion,seated,side:state.side});c.setPose('therapist','rest',0);c.setPose('patient','rest',0);state.role='therapist';
  $('poseLabel').textContent=`治疗师站立面对患者；患者${seated?'坐位':'站位'}。仅主动方向演示，无手法接触。`;clinicUI();
 }
 }catch(e){$('modelCanvas').hidden=true;$('clinicCanvas').hidden=true;$('diagram').hidden=false;$('diagram').innerHTML=`<div class="readcard"><h2>模型未能加载</h2><p>${E(e.message)}</p><p>请使用支持WebGL的浏览器。文字和原理示意仍保留。</p></div>`;console.error(e);}
 $('caption').textContent=state.mode==='anatomy'?'静态结构用于定位，未与当前运动逐点配准；颜色不表示肌肉张力。':state.mode==='clinic'?'复用已有角色蒙皮，姿态为定性示范；不是临床量角、接触或患者能力测量。':a.caption;
 if(a.kind==='pelvis3d')$('poseLabel').textContent='蓝：固定水平参照　橙：骨性教学标记连线。可旋转核对三维关系。';
 renderExtras();update();}
function settingInput(key,label,min,max,step,value){return `<label>${label}<input type="number" data-setting="${key}" min="${min}" max="${max}" step="${step}" value="${value}"></label>`;}
function renderExtras(){const a=state.action,p=state.params;let h='';
 if(a.kind==='muscleAtlas'){const e=atlasEntry();h=atlasBrowser()+(!e.referenceOnly?`<div class="atlas-toolbar atlas89-controls"><span class="atlas89-controltitle">观察 ${E(e.name)}</span>${e.markerMode!=='none'?`<label><input type="checkbox" data-setting="showAtlasPull" ${p.showAtlasPull?'checked':''}>牵拉方向示意</label><label><input type="checkbox" data-setting="showAtlasMarkers" ${p.showAtlasMarkers?'checked':''}>起止区域示意</label>`:''}<button data-stage="0">起始定位</button>${!a.static?'<button data-stage="1">动作末位</button>':''}<button data-cmd="atlasFocus">聚焦肌肉</button><button data-cmd="atlasAll">看全骨段</button><button data-cmd="atlasRecord">记录结构</button></div>`:'');}
 else if(state.mode==='anatomy'){h='<label>结构层<select data-setting="layer"><option value="both">骨骼＋相关肌群</option><option value="bone">骨骼</option><option value="muscle">肌肉 / 肌腱</option></select></label><button data-cmd="isolate">单独观察所选</button><button data-cmd="restore">恢复局部</button><span class="chapterstatus">此参照为源资源原侧别，不镜像充作另一侧。</span>';}
 else if(state.mode==='clinic'||a.kind==='native'){
  if(state.mode==='clinic'||a.kind==='native'&&getModel().supportsSide(a))h+=`<label>观察侧<select id="sideSelect"><option value="right" ${state.side==='right'?'selected':''}>右侧</option><option value="left" ${state.side==='left'?'selected':''}>左侧</option></select></label><span class="chapterstatus">模型分别有左右结构；相同条件对照，不预设患侧诊断。</span>`;
  if(state.mode!=='clinic'&&a.kind==='native')h+=`<label><input type="checkbox" data-setting="showPaths" ${p.showPaths?'checked':''}>看源肌腱路径</label>`;
 }else if(a.kind==='pelvis3d'){
  h=(a.motion==='pelvisSupport'?`<label>支撑侧<select data-setting="supportSide"><option value="right" ${p.supportSide==='right'?'selected':''}>右侧支撑</option><option value="left" ${p.supportSide==='left'?'selected':''}>左侧支撑</option></select></label>`:'')+[[0,a.motion==='pelvisSupport'?'对侧上提':'后倾'],[.5,'中立参照'],[1,a.motion==='pelvisSupport'?'对侧下降':'前倾']].map(([t,l])=>`<button data-stage="${t}">${l}</button>`).join('');
 }else if(a.kind==='dynamometer'){
  h=`<div class="safetyline">${(a.motion==='isokinetic'?['核对关节与测力轴线','确认体位与固定带','完成模拟校零']:['核对固定体位与传感器位置','确认固定带与停止方式','完成模拟校零']).map((v,i)=>`<label><input type="checkbox" data-device-check="${i}" ${p.deviceChecks[i]?'checked':''} ${p.deviceStarted||p.deviceStopped?'disabled':''}>${v}</label>`).join('')}</div>`;
  if(a.motion!=='mvcBack')h+=`<label>示例侧别<select data-setting="deviceSide"><option value="right" ${p.deviceSide==='right'?'selected':''}>右侧案例</option><option value="left" ${p.deviceSide==='left'?'selected':''}>左侧案例</option></select></label>`;
  if(a.motion==='isokinetic')h+=`<label>设定角速度<select data-setting="isoSpeed">${[30,60,90,120,180].map(v=>`<option value="${v}" ${p.isoSpeed===v?'selected':''}>${v}°/s</option>`).join('')}</select></label><label>工作方向<select data-setting="isoDirection"><option value="concentric" ${p.isoDirection==='concentric'?'selected':''}>向心屈肘</option><option value="eccentric" ${p.isoDirection==='eccentric'?'selected':''}>离心控制</option></select></label><label>发力过程<select data-setting="deviceProfile">${[['rise','渐增用力'],['steady','相对稳定'],['wave','增强再减弱']].map(([k,l])=>`<option value="${k}" ${p.deviceProfile===k?'selected':''}>${l}</option>`).join('')}</select></label>`+settingInput('deviceEffort','相对发力',.2,1,.1,p.deviceEffort);
  else h+=`<label>努力条件<select data-setting="mvcEffort"><option value="max" ${p.mvcEffort==='max'?'selected':''}>最大努力指令</option><option value="submax" ${p.mvcEffort==='submax'?'selected':''}>次最大努力对照</option></select></label>`;
  h+='<button data-device="start" class="primary">开始试次</button><button data-device="stop" class="danger">停止试次</button><button data-device="new">新试次</button>';
 }else if(a.kind==='candle'){
  h=settingInput('airDistance',a.motion==='distance'?'远处距离（cm）':'相同目标距离（cm）',a.motion==='distance'?20:10,80,5,p.airDistance)+settingInput('airLevel','相对输出',1,3,1,p.airLevel)+(a.motion==='duration'?settingInput('airSeconds','较长输出（s）',5,15,1,p.airSeconds):'');
 }else if(a.kind==='footprints'){h='<span class="chapterstatus">三种假设接触情景同步播放；红点与轨迹由相同的压力权重计算。</span>';
 }else if(a.kind==='metabolism'){h=vitalSettings();
 }else if(a.kind==='contraction'){
  h=`<label>方向示意<select data-setting="muscleSide"><option value="right" ${p.muscleSide==='right'?'selected':''}>右侧</option><option value="left" ${p.muscleSide==='left'?'selected':''}>左侧（镜像）</option></select></label>`;
  if(a.contractionMode==='isokinetic')h+=`<label>恒速工作段<select data-setting="isoSpeed">${[30,60,90,120,180].map(v=>`<option value="${v}" ${p.isoSpeed===v?'selected':''}>${v}°/s</option>`).join('')}</select></label><label>长度变化<select data-setting="isoDirection"><option value="concentric" ${p.isoDirection==='concentric'?'selected':''}>向心阶段</option><option value="eccentric" ${p.isoDirection==='eccentric'?'selected':''}>离心阶段</option></select></label>`+settingInput('effort','相对发力',.2,1,.1,p.effort);
  }else if(a.kind==='breath'){
  h='<span class="chapterstatus">同一呼吸周期：</span>'+[[0,'呼气末'],[.25,'吸气中'],[.5,'吸气末'],[.75,'呼气中'],[1,'回到呼气末']].map(([t,l])=>`<button data-stage="${t}">${l}</button>`).join('');
 }else if(a.kind==='failure'){
  h=`<label>定义材料<select data-setting="material"><option value="A" ${p.material==='A'?'selected':''}>材料A</option><option value="B" ${p.material==='B'?'selected':''}>材料B</option></select></label><label>观察过程<select data-setting="failureMode">${[['loading','持续加载至失效'],['elasticUnload','弹性范围内卸载'],['plasticUnload','塑性区加载后卸载']].map(([k,l])=>`<option value="${k}" ${p.failureMode===k?'selected':''}>${l}</option>`).join('')}</select></label>`+[[0,'起始'],[.5,'过程中央'],[1,'过程终末']].map(([t,l])=>`<button data-stage="${t}">${l}</button>`).join('');
 }else if(a.kind==='archStructure'){h=`<label>观察弓形<select data-setting="archView"><option value="medial">内侧纵弓</option><option value="lateral">外侧纵弓</option><option value="transverse">横弓</option></select></label>`;
 }else if(a.kind==='lever')h=settingInput('load','阻力（N）',5,300,5,p.load)+(a.motion==='elbow'?settingInput('distance','负荷距肘中心（m）',.18,.42,.02,p.distance):'');
 else if(a.kind==='bone')h=`<label>例题材料<select data-setting="material"><option value="A" ${p.material==='A'?'selected':''}>材料A · E=6000 MPa</option><option value="B" ${p.material==='B'?'selected':''}>材料B · E=12000 MPa</option></select></label><label><input type="checkbox" data-setting="loadCycle" ${p.loadCycle?'checked':''}>同一滑条加载—卸载</label>`;
 else if(a.kind==='capacity'){
  if(a.motion==='rm')h=`<label>案例负荷<select data-setting="mass">${[8,10,12,14].map(x=>`<option value="${x}" ${x===p.mass?'selected':''}>${x} kg</option>`).join('')}</select></label>`;
  if(a.motion==='jump')h=settingInput('takeoff','合成起跳速度（m/s）',1.5,3.5,.1,p.takeoff);
  if(['plank','squat'].includes(a.motion))h=`<label>任务条件<select data-setting="quality"><option value="valid" ${p.quality==='valid'?'selected':''}>符合姿势 / 深度要求</option><option value="changed" ${p.quality!=='valid'?'selected':''}>姿势偏离 / 深度不足</option></select></label>`;
 }else if(a.kind==='airflow')h=settingInput('airDistance','目标距离（cm）',10,70,5,p.airDistance)+settingInput('airLevel','相对输出级别',1,3,1,p.airLevel)+settingInput('airSeconds','案例持续输出（s）',3,15,1,p.airSeconds);
 else if(a.kind==='cardio')h=vitalSettings()+safetyHTML();
 else if(a.kind==='chain')h='<span class="chapterstatus">两栏同一进度：足游离与足固定分开比较。轨迹线不是肌力箭头。</span>';
 else if(a.kind==='arch')h=`<label>弓高<select data-setting="morphology">${[['low','低弓'],['typical','常见弓高'],['high','高弓']].map(([v,l])=>`<option value="${v}" ${p.morphology===v?'selected':''}>${l}</option>`).join('')}</select></label>`+settingInput('stiffness','等效刚度（相对）',.6,2,.1,p.stiffness);
 else if(a.kind==='cop')h='<label class="filebutton">导入COP CSV<input id="copFile" type="file" accept=".csv"></label><button data-cmd="copExample">数据模板</button><button data-cmd="copReset">恢复合成样例</button>';
 else if(a.kind==='scapula')h=`<label>肩胛运动<select data-setting="scap">${[['elevation','上提 / 下降'],['protraction','前伸 / 后缩'],['rotation','上回旋 / 下回旋']].map(([v,l])=>`<option value="${v}" ${p.scap===v?'selected':''}>${l}</option>`).join('')}</select></label>`;
 else if(a.kind==='gait')h=`<span class="chapterstatus">关键事件：</span>${[[0,'右初始接触'],[.10,'左足离地'],[.50,'左初始接触'],[.60,'右足离地'],[1,'下次右初始接触']].map(([t,x])=>`<button data-gait-event="${t}">${x}</button>`).join('')}`;
 $('extras').innerHTML=h;
 if(a.kind==='gait'&&a.motion==='short')$('extras').innerHTML='<span class="chapterstatus">本特征例：左足离地15% / 左初始接触42% / 右足离地50%。不是正常标准。</span>'+[[0,'右初始接触'],[.15,'左足离地'],[.42,'左初始接触'],[.5,'右足离地']].map(([t,l])=>`<button data-gait-event="${t}">${l}</button>`).join('');
}
function vitalSettings(){const p=state.params,locked=state.action?.kind==='cardio'&&p.safety.state!=='idle';let h=settingInput('caseAge','案例年龄',16,80,1,p.caseAge)+settingInput('caseRestHR','案例安静心率',40,110,1,p.caseRestHR)+settingInput('targetPct','计算练习：HRmax比例（%）',40,90,1,p.targetPct);return locked?h.replaceAll('<input ','<input disabled '):h;}
function atlasEntry(){return MuscleCatalogue.find(e=>e.id===state.action?.muscleID);}
function atlasStatus(e){return e.referenceOnly?'功能资料 · 暂无网格':e.dynamic?'模型＋附着动作':'模型＋结构定位';}
function atlasSide(e){return e.side==='right'?'右侧':e.side==='left'?'左侧':e.side==='bilateral'?'双侧集合':'完整结构';}
function atlasNorm(s){return String(s||'').normalize('NFKC').toLowerCase().replace(/[\s·•（）()、，,\-_/]+/g,'');}
function atlasResults(){const q=atlasNorm(state.atlasQuery),grp=state.atlasGroup||'all';return state.module.actions.filter(a=>{const e=MuscleCatalogue.find(x=>x.id===a.muscleID);return (grp==='all'||e.group===grp)&&(!q||[e.name,e.english,...(e.aliases||[])].some(s=>atlasNorm(s).includes(q)));});}
function atlasBrowser(){
 const q=state.atlasQuery||'',group=state.atlasGroup||'all',e=atlasEntry(),items=atlasResults(),groups=[...new Set(MuscleCatalogue.map(x=>x.group))],real=MuscleCatalogue.filter(x=>!x.referenceOnly).length;
 return `<div class="atlas89-browser"><div class="atlas89-searchline"><label class="atlas89-search">全身查找肌肉<input id="atlasSearch" type="search" value="${E(q)}" placeholder="中文、英文或常用别名，如臀大肌 / biceps" autocomplete="off"></label><label>按部位浏览<select id="atlasGroup"><option value="all">全部部位</option>${groups.map(g=>`<option value="${E(g)}" ${g===group?'selected':''}>${E(g)}（${MuscleCatalogue.filter(x=>x.group===g).length}）</option>`).join('')}</select></label><button type="button" data-atlas-reset>查看全部肌肉</button></div>
 <div class="atlas89-resultline"><span id="atlasSearchCount" aria-live="polite">${q?'全身搜索':'当前目录'} ${items.length} / ${MuscleCatalogue.length} 项 · ${real}项有源网格</span><span class="atlas89-current" id="atlasCurrent">当前：${E(e?.name||'')} · ${E(e?atlasStatus(e):'')}</span></div>
 <details class="atlas89-directory" open><summary>主要肌肉目录 <span>点击名称即可切换模型与说明</span></summary><div id="atlasResults" class="atlas89-resultgrid">${atlasCards(items)}</div></details>
 <p class="atlas89-querynote">${q?'名称搜索覆盖全部部位，不会被上一次的部位筛选挡住。':'按部位浏览会清空旧关键词。分部和肌群条目已注明；并非全身每块肌肉均有模型。'}</p>
 </div>`;
}
function atlasCards(items){return items.map(a=>{const e=MuscleCatalogue.find(x=>x.id===a.muscleID);return `<button type="button" class="atlas89-card ${state.action?.id===a.id?'is-active':''} ${e.referenceOnly?'is-reference':''}" data-atlas-pick="${a.id}" aria-pressed="${state.action?.id===a.id}"><b>${E(e.name)}</b><span>${E(e.group)}</span><small>${atlasStatus(e)}</small></button>`;}).join('')||`<div class="atlas89-noresults">没有匹配项。请检查名称，或点击“查看全部肌肉”。</div>`;}
function syncAtlasOptions(){
 if(state.module?.id!==201)return;const items=atlasResults(),select=$('actionSelect');select.disabled=!items.length;
 if(!items.length){select.innerHTML='<option value="" selected>未找到匹配肌肉</option>';return;}
 select.innerHTML=[...new Set(items.map(a=>a.group))].map(g=>`<optgroup label="${E(g)}">${items.filter(a=>a.group===g).map(a=>`<option value="${a.id}">${E(a.title)}</option>`).join('')}</optgroup>`).join('');select.value=state.action?.id||'';
}
function renderAtlasDetails(){const e=atlasEntry();if(!e)return;$('atlasDetails').innerHTML=`<div><h3>${E(e.name)}</h3><strong>${E(e.group)} · ${E(atlasSide(e))}</strong><p>${E(e.location||'先旋转模型，对照邻近骨骼辨认该肌肉。')}</p><p class="atlas89-function"><b>主要功能：</b>${E(e.function)}</p><p id="atlasSelected" class="atlas-selected">${e.referenceOnly?'本条目无模型；不会借用其他肌肉替代。':'点击模型可查看当前结构的名称。'}</p></div><div><h3>起止点与附着</h3><p><strong>起点／近端附着：</strong>${E(e.origin)}</p><p><strong>止点／远端附着：</strong>${E(e.insertion)}</p><small>肌肉也可附着于筋膜、肌腱或腱膜；不能一律理解成两端都直接连骨。</small></div><div><h3>${e.dynamic?'附着关系与动作':'本条目显示范围'}</h3><p>${e.referenceOnly?'当前模型包没有该肌肉的独立、可核对网格。此处提供位置、附着与功能资料，不显示替代性假模型。':e.dynamic?'播放或拖动进度，观察一个主要动作中的骨段与附着关系。':'已有独立源网格，支持旋转、点选和局部观察；复杂动作未完成绑定，保留静态定位。'}</p><small>${e.dynamic?'蓝、橙标记为教学起止区域；不是实测骨性标志。动画不求解肌力、组织应力或肌纤维长度。':'起止点优先以解剖文字说明，不把网格最高点、最低点冒充精确附着点。'}</small><p class="atlas89-source">${e.sourceIDs.length?'源结构：'+E(e.sourceIDs.join(' / ')):'网格状态：未提供'} </p></div>`;}
function atlasReferenceView(){const e=atlasEntry();return `<section class="atlas89-reference-panel"><div class="atlas89-reference-mark">功能资料</div><h2>${E(e.name)}</h2><p>${E(e.location)}</p><div class="atlas89-reference-summary"><b>主要功能</b><p>${E(e.function)}</p></div><p>此肌肉的源网格尚未包含在当前资源包中。本页只展示解剖资料，未用其他结构代替。起止点见下方说明。</p><button type="button" data-atlas-reset>返回全部肌肉目录</button></section>`;}
function atlasEmptyView(){
 pause();state.atlasEmpty=true;$('modelCanvas').hidden=true;$('clinicCanvas').hidden=true;$('atlasLabels').hidden=true;$('atlasDetails').hidden=true;$('viewButtons').hidden=true;$('transport').hidden=true;$('diagram').hidden=false;$('numbers').hidden=true;$('stageBadge').textContent='当前搜索无匹配项';
 $('command').textContent='未找到“'+(state.atlasQuery||'')+'”。请改用标准名称、英文名称，或查看全部肌肉。';
 $('diagram').innerHTML='<div class="atlas89-reference-panel"><h2>没有匹配的肌肉</h2><p>为了避免名称与模型不一致，此时不继续展示上一块肌肉。</p><button type="button" data-atlas-reset>清除条件，显示全部肌肉</button></div>';
 $('atlasCurrent').textContent='当前无匹配，已隐藏上一块肌肉';
}
function filterMuscles(reason='input'){
 if(state.module?.id!==201)return;const focused=document.activeElement?.id==='atlasSearch',pos=$('atlasSearch')?.selectionStart;
 if(reason==='group'){state.atlasGroup=$('atlasGroup').value;state.atlasQuery='';$('atlasSearch').value='';}
 else if(reason==='reset'){state.atlasGroup='all';state.atlasQuery='';$('atlasSearch').value='';$('atlasGroup').value='all';}
 else{state.atlasQuery=$('atlasSearch')?.value||'';state.atlasGroup='all';$('atlasGroup').value='all';}
 const items=atlasResults();syncAtlasOptions();$('atlasSearchCount').textContent=(state.atlasQuery?'全身搜索':'当前目录')+' '+items.length+' / '+MuscleCatalogue.length+' 项';$('atlasResults').innerHTML=atlasCards(items);
 if(!items.length){atlasEmptyView();return;}
 const exact=items.find(a=>{const e=MuscleCatalogue.find(x=>x.id===a.muscleID);return [e.name,e.english,...(e.aliases||[])].some(x=>atlasNorm(x)===atlasNorm(state.atlasQuery));});
 const target=exact||items.find(a=>a.id===state.action?.id)||items[0];const needs=state.atlasEmpty||target.id!==state.action?.id;state.atlasEmpty=false;
 if(needs)selectAction(target.id);else {$('atlasResults').innerHTML=atlasCards(items);syncAtlasOptions();}
 if(focused){const input=$('atlasSearch');input.focus({preventScroll:true});try{input.setSelectionRange(pos,pos);}catch{}}
}

function safetyHTML(){const z=state.params.safety;return `<div class="safetyline">${['设备和急停位置','安全装置与衣鞋','患者说明与旁站协助'].map((txt,i)=>`<label><input type="checkbox" data-check="${i}" ${z.checked[i]?'checked':''} ${z.state==='running'||z.state==='emergency'?'disabled':''}>${txt}</label>`).join('')}</div><div class="row"><button data-risk="start" class="primary">启动本试次</button><button data-risk="stop" class="danger">急停 / 结束运动</button><button data-risk="dizzy">演练：头晕但清醒</button><button data-risk="breathing">演练：无反应，有正常呼吸</button><button data-risk="arrest">演练：无反应，无正常呼吸</button><button data-risk="new">新试次</button></div>${z.state==='emergency'?`<p class="risknote">${E(z.message)}</p><div class="row"><button data-risk="call">启动求助与急救系统</button><button data-risk="monitor">监护呼吸并按指令协助</button><button data-risk="cpr">按受训流程CPR并取AED</button></div>`:''}`;}
function risk(cmd){pause();const z=state.params.safety;if(cmd==='new'){z.checked=[false,false,false];z.state='idle';z.message='新的试次：重新核对条件。';z.type=null;state.t=0;}
 else if(cmd==='start'){if(!z.checked.every(Boolean)){toast('三项准备尚未全部核对。');return;}if(['emergency','stopped'].includes(z.state)){toast('此试次已结束，请新开试次。');return;}z.state='running';z.message='已按教学清单确认准备；这不等于验证了真实设备安全。';event('treadmill_start');play();}
 else if(cmd==='stop'){z.state='stopped';z.message='皮带与动画已停止；最后显示值不是停止后的新测量。';event('treadmill_stop');}
 else if(['dizzy','breathing','arrest'].includes(cmd)){z.state='emergency';z.type=cmd;z.called=false;z.message=cmd==='dizzy'?'已停止：患者清醒但头晕。稳定与评估，按机构规程求助；不能直接做胸外按压。':cmd==='breathing'?'已停止：无反应但仍有正常呼吸。启动求助，持续监护；不能因此自动判为心脏骤停。':'已停止：无反应、无正常呼吸。本病例给定：受训人员在10秒内不能明确触及脉搏。启动急救系统并取AED，依成人BLS流程行动。';event('risk',{condition:cmd,stage:state.t});}
 else if(cmd==='call'){if(z.state!=='emergency')return;z.called=true;z.message+=' 已记录求助（模拟，不拨出电话）。';event('emergency_call');}
 else if(cmd==='cpr'){if(z.type!=='arrest'){toast('本病例不满足这个分支。先核对意识和正常呼吸。');return;}if(!z.called){toast('请先启动急救系统；真实情境可由同伴同时呼救和取AED。');return;}z.message='已进入成人BLS教学路径：按受训流程CPR，尽快使用AED并遵从语音指令。此页面不判定按压深度、效果或复苏成功。';event('bls_path');}
 else if(cmd==='monitor'){if(z.type==='arrest'){toast('无正常呼吸时不能只等待监护；应及时启动复苏路径。');return;}z.message='持续监护与协助，按机构规程和急救人员指令处理；症状未被程序自动清除。';event('monitor_path');}
 renderExtras();update();}
function calcHTML(a){if(a.kind==='hr')return `<div class="calc"><h2>心率储备法 · 给定例题</h2><label>安静心率（次/分）<input id="restHR" type="number" value="70"></label><label>题目给定最大心率（次/分）<input id="maxHR" type="number" value="185"></label><label>目标强度（%）<input id="intensity" type="number" value="50"></label><button data-cmd="calcHR">计算</button><p id="calcResult" class="result">靶心率 = 安静心率 +（最大心率 − 安静心率）×强度</p><p class="subtle">最大心率来源已在题干给定；不根据年龄自动生成个人运动处方。</p></div>`;
 return ''; }
function update(){const a=state.action;if(!a||state.atlasEmpty)return;$('stageRange').value=Math.round(state.t*1000);$('stageOutput').textContent=Math.round(state.t*100)+'%';
 if(nativeActive()){if(state.model?.action){state.model.action.showMusclePaths=state.params.showPaths;state.model.action.showAtlasMarkers=state.params.showAtlasMarkers;state.model.action.showAtlasPull=state.params.showAtlasPull;}state.model?.update(state.t);}
 else if(state.mode==='clinic'){
  if(state.clinic&&state.consented&&!['stopped','closed'].includes(state.clinicPhase)){let t=state.t;if(a.kind==='contraction'){const d=MM.contraction(a.contractionMode,state.t,state.params);t=(d.angle-20)/85;}else{if(a.id==='03-3')t=1-state.t;if(a.id==='03-2')t=.5;}state.clinic.setPose(state.role,roleMotion(a),t);}
 }else if(a.kind==='muscleAtlas'&&a.referenceOnly){$('diagram').innerHTML=atlasReferenceView();}
 else if(a.kind==='hr'){if(!$('calcResult'))$('diagram').innerHTML=calcHTML(a);}
 else {$('diagram').innerHTML=Sim.render(a,state.t,state.params);}
 $('playBtn').textContent=state.playing?'Ⅱ 暂停':'▶ 播放';guide();numbers();}
function guide(){const a=state.action,o=$('motionGuide'),m=MuscleTeaching[a.motion];o.hidden=!(a.kind==='native'&&state.mode==='main');if(o.hidden)return;
 o.innerHTML=`<b>${E(m?.joint||'当前骨段')}</b><span>面：${E(m?.plane||'源局部平面')}</span><span>轴：${E(m?.axis||'源模型关节轴')}</span><span>方向：${E(m?.direction||'沿箭头观察')}</span><button id="guideToggle">${state.model?.showGuides?'隐藏':'显示'}轴面</button>`;
 $('guideToggle').onclick=()=>{state.model.showGuides=!state.model.showGuides;state.model.update(state.t);guide();};}
function numbers(){const a=state.action,p=state.params;let ar=[];
 if(state.mode==='main'){
  if(a.kind==='gait'){const v=MM.gaitPhase(state.t,['clearance','lean','short'].includes(a.motion)?a.motion:'normal',p.gaitSide||'right');$('command').textContent=state.action.command+' 当前：'+v.support+'；右侧 '+MM.gaitPhase(state.t,['clearance','lean','short'].includes(a.motion)?a.motion:'normal','right').phase+'。';}
  if(a.kind==='pelvis3d'){const d=state.model?.pelvisState;if(d)ar=[d.label,`相对模型角 ${Math.abs(d.angle).toFixed(1)}°`,d.support?`${d.side==='right'?'右':'左'}侧支撑`:'骨盆与骶骨整体观察'];}
  if(a.kind==='dynamometer'){const d=MM.deviceTrial(a.motion,state.t,p);ar=[d.phase,`${d.jointLabel||'肘角'} ${d.angle.toFixed(1)}°`,d.mvc?`已观察峰值 ${d.peak.toFixed(1)} N`:`装置力矩 ${Math.abs(d.deviceTorque).toFixed(1)} N·m`];}
  if(a.kind==='candle'){const d=MM.candleTrial(state.t,p);ar=[`场景时间 ${d.time.toFixed(1)} s`,`目标距离 ${d.distance} cm`,'仅虚拟火焰观察'];}
  if(a.kind==='lever'){const d=MM.lever(a.motion,state.t,p);ar=[`动力臂 ${d.effortArm.toFixed(4)} m`,`阻力臂 ${d.resistanceArm.toFixed(4)} m`,`机械优势 ${d.advantage.toFixed(2)}`,`理想平衡动力 ${d.effort.toFixed(1)} N`];}
  if(a.kind==='bone'){const d=MM.bone(a.motion,state.t,p);ar=[`相对载荷 ${(100*d.level).toFixed(0)}%`,`${d.xLabel}：${d.x.toExponential(2)}`,`${d.yLabel}：${d.y.toFixed(2)}`];}
  if(a.kind==='capacity'){const d=MM.capacity(a.motion,state.t,p);ar=a.motion==='rm'?[`${d.mass} kg`,d.label]:a.motion==='jump'?[d.phase,`腾空 ${d.flight.toFixed(3)} s`,`质心上升 ${(d.heightMax*100).toFixed(1)} cm`]:a.motion==='plank'?[d.label,`有效时间 ${d.validSeconds.toFixed(1)} s`]:[d.label,`完整 ${d.total} 次 / 有效 ${d.valid} 次`];}
  if(a.kind==='contraction'){const d=MM.contraction(a.contractionMode,state.t,p);ar=[d.type,`肘角 ${d.angle.toFixed(1)}°`,`案例时间 ${d.time.toFixed(2)} s`,d.device?`匀速工作段 ${Math.abs(d.velocity).toFixed(0)}°/s`:d.lengthChange];}
  if(a.kind==='ssc'){const d=MM.capacity('ssc',state.t,p);ar=[d.phase,`股四头肌：${d.quad}`,`小腿三头肌：${d.calf}`];}
  if(a.kind==='gait'){const d=MM.gaitSync87(state.t,['clearance','lean','short'].includes(a.motion)?a.motion:'normal');ar=[d.support,`模型上下全幅 ${d.b.vertical.span.toFixed(1)} mm`,`模型左右全幅 ${d.b.lateral.span.toFixed(1)} mm`,`距最低点 ${d.vertical_mm.toFixed(1)} mm`,d.lateralTrend];}
  if(a.kind==='breath'){const d=MM.breath87(state.t,a.motion);ar=[d.phase,'胸式：左右径扩大更突出','腹式：膈肌下降、上下径增加更突出'];}
  if(a.kind==='failure'){const d=MM.failure87(state.t,p);ar=[d.phase,`应变 ${(d.eps*100).toFixed(3)}%`,`此状态卸载残余 ${(d.currentResidual*100).toFixed(3)}%`];}
  if(a.kind==='airflow'){const d=MM.airflow(state.t,p);ar=[`距离 ${d.distance} cm`,`偏转 ${d.deflection.toFixed(0)}°`,`达到示意阈值 ${d.validSeconds.toFixed(1)} s`];}
  if(a.kind==='cardio'){const d=MM.vitals88(a.motion,state.t,p);ar=[p.safety.state==='running'?'本次运行':'本次未运行',`心率 ${Math.round(d.hr)}次/分`,`血压 ${Math.round(d.sbp)}/${Math.round(d.dbp)} mmHg`];}
  if(a.kind==='muscleAtlas'&&!a.referenceOnly){const d=state.model?.atlasReadout;if(d)ar=[d.motion,atlasSide(atlasEntry())+' · '+d.joint,state.action.static?'结构定位 · 暂不播放动作':'附着区随骨段变化'];}
  if(a.kind==='arch'){const d=Sim.Data.arch(state.t,p);ar=[`加载 ${Math.round(d.load*100)}%`,`等效刚度 ${d.stiffness}`,state.t<.5?'加载储能':'卸载释能'];}
  if(a.kind==='cop'){const ds=p.copData||Sim.syntheticCOP(),r=ds.rows[Math.min(ds.rows.length-1,Math.floor(state.t*(ds.rows.length-1)))];ar=[ds.label,r.fz_N>0?`COP (${r.x_mm.toFixed(1)}, ${r.y_mm.toFixed(1)}) mm`:'无正接触力，不显示COP'];}
 }
 $('numbers').hidden=!ar.length;$('numbers').innerHTML=ar.map(x=>`<span>${E(x)}</span>`).join('');}
function eligible(){if(state.atlasEmpty)return false;if(state.mode==='anatomy'||state.action?.static||state.action?.kind==='concept')return false;if(state.mode==='clinic')return state.consented&&!['stopped','closed','intro'].includes(state.clinicPhase);if(state.action?.kind==='cardio')return state.params.safety.state==='running';if(state.action?.kind==='dynamometer')return state.params.deviceStarted&&!state.params.deviceStopped;return true;}
function setStage(t,user=true){if(user){pause();if(!eligible()){toast(state.mode==='clinic'?'请先确认愿意参与，再示范或模仿。':state.action?.kind==='cardio'?'请先完成安全准备，启动新试次。':state.action?.kind==='dynamometer'?'请先核对轴线、固定与校零，再开始新试次。':'静态参考不播放动作。');return false;}}if(!Number.isFinite(Number(t)))return false;state.t=MM.clamp(Number(t));update();return true;}
function pause(){state.playing=false;cancelAnimationFrame(state.raf);state.raf=0;state.lastTime=0;state.clinic?.stop();$('playBtn').textContent='▶ 播放';}
function play(){if(state.playing){pause();return;}if(!eligible()){toast(state.mode==='clinic'?'先取得同意，再演示。':state.action?.kind==='dynamometer'?'先完成轴线核对、体位固定和校零，再开始试次。':state.action?.kind==='cardio'?'先完成安全准备，再启动试次。':'当前页面只作静态观察。');return;}
 if(state.t>=1){state.t=0;state.dir=1;}state.playing=true;state.lastTime=0;
 const duration=state.action.kind==='dynamometer'?MM.deviceTrial(state.action.motion,state.t,state.params).duration:state.action.kind==='footprints'?9:state.action.kind==='candle'?10:state.action.kind==='metabolism'?10:state.action.kind==='contraction'?MM.contraction(state.action.contractionMode,state.t,state.params).duration:state.action.kind==='gait'?7:state.action.kind==='ssc'?7:state.action.kind==='capacity'?9:7;
 function tick(now){if(!state.playing)return;const dt=state.lastTime?Math.min(.08,(now-state.lastTime)/1000):0;state.lastTime=now;state.t+=dt/duration*Number($('speed').value)*state.dir;
  if(state.t>=1||state.t<=0&&state.dir<0){if($('roundtrip').checked&&!['contraction','cardio','capacity','ssc','gait','airflow','dynamometer','candle','metabolism','footprints','failure'].includes(state.action.kind)){state.t=MM.clamp(state.t);state.dir*=-1;}else{state.t=MM.clamp(state.t);update();pause();if(state.mode==='clinic'){state.clinicPhase=state.role==='therapist'?'demonstrated':'observed';clinicUI();}return;}}
  update();state.raf=requestAnimationFrame(tick);
 }state.raf=requestAnimationFrame(tick);}
function clinicUI(){const a=state.action,phase=state.clinicPhase;let text,buttons;
 if(phase==='intro'){text='物理治疗师：您好，我是林老师。现在做个简短活动观察，我先示范，请您在舒适范围内模仿；不舒服随时告诉我，可以吗？';buttons=[['agree','患者同意 · 看示范'],['decline','患者暂不参与']];}
 else if(phase==='stopped'||phase==='closed'){text=phase==='stopped'?'物理治疗师：好的，我们停止。谢谢您告诉我感受，没有观察到的内容不补记。':'物理治疗师：今天的观察到这里。谢谢您的配合，先回到舒服的姿势休息。';buttons=[['new','重新开启一次交流']];}
 else {text=state.role==='therapist'?'物理治疗师：请看我示范动作方向。稍后请您模仿，我们再看另一侧。':'物理治疗师：'+a.command+'　模拟患者：好的，我会慢慢做，不舒服就停。';buttons=[['demo','治疗师示范'],['perform','请患者模仿'],['stop','患者请求停止'],['finish','感谢并结束']];}
 $('clinicBar').innerHTML=`<p>${E(text)}</p><div class="row">${buttons.map(([cmd,label])=>`<button data-clinic="${cmd}" class="${cmd==='stop'?'danger':''}">${label}</button>`).join('')}</div>`;
}
function clinicCommand(cmd){pause();if(cmd==='new'){state.clinicPhase='intro';state.consented=false;state.t=0;configure();return;}
 if(cmd==='agree'){if(state.clinicPhase!=='intro'){toast('需先重新开启交流并确认参与。');return;}state.consented=true;state.clinicPhase='demonstrating';state.role='therapist';event('consent');}
 else if(cmd==='decline'){state.clinicPhase='closed';state.consented=false;event('declined');clinicUI();return;}
 else if(cmd==='stop'){state.clinicPhase='stopped';state.consented=false;event('patient_stop',{stage:state.t,side:state.side});clinicUI();return;}
 else if(cmd==='finish'){state.clinicPhase='closed';state.consented=false;state.clinic?.setPose('therapist','rest',0);state.clinic?.setPose('patient','rest',0);event('closing');clinicUI();return;}
 else if(!state.consented||['closed','stopped','intro'].includes(state.clinicPhase)){toast('当前不可继续患者动作。');return;}
 else if(cmd==='perform'){state.role='patient';state.clinicPhase='observing';event('patient_instruction',{side:state.side});}
 else if(cmd==='demo'){state.role='therapist';state.clinicPhase='demonstrating';event('demonstration',{side:state.side});}
 state.t=0;state.dir=1;const other=state.role==='therapist'?'patient':'therapist';state.clinic.setPose(other,'rest',0);clinicUI();play();}
function sourceHTML(ids){return [...new Set(ids)].map(id=>{const s=payload.sources[id];return s?`<p><b>${E(id)}</b> ${s.url?`<a target="_blank" rel="noopener" href="${E(s.url)}">${E(s.title)}</a>`:E(s.title)}<br><small>${E(s.supports||'')} ${E(s.limit||'')}</small></p>`:'';}).join('');}
function muscleHTML(){const m=MuscleTeaching[state.action.motion];if(!m)return'';return `<div class="musclecols"><div><b>功能方向 / 主要参与</b>${m.main.map(x=>`<div>${E(x)}</div>`).join('')}</div><div><b>相反方向的功能肌群</b>${m.opposite.map(x=>`<div>${E(x)}</div>`).join('')}</div></div><p class="scope">${E(m.note)} 此处是功能对照，不是实际肌电或力的计算。</p>`;}
function notes(){const a=state.action;$('notesBody').innerHTML=`<h3>操作条件</h3><p>${E(a.position)}</p><p>${a.observe.map(E).join('；')}</p><h3>理解这个动作</h3><p>${E(a.teach)}</p><p class="scope">${E(a.pitfall)}</p><details><summary>资料与模型说明</summary>${sourceHTML(a.refs||[])}</details>`;}
function revisionReadout87(){if(state.action?.kind==='muscleAtlas'){const e=MuscleCatalogue.find(x=>x.id===state.action.muscleID);return {muscle:e.name,sourceIDs:e.sourceIDs,originText:e.origin,insertionText:e.insertion,function:e.function,dynamic:e.dynamic,reference:state.model?.atlasReadout,scope:e.animationScope};}const a=state.action,p=state.params;
 if(a.kind==='gait'){const d=MM.gaitSync87(state.t,['clearance','lean','short'].includes(a.motion)?a.motion:'normal');return {synthetic:true,phase:d.phase,support:d.support,com:d.d.com,time:d.time,vertical_mm:d.vertical_mm,lateral_mm:d.lateral_mm,verticalFull_mm:d.b.vertical.span,lateralFull_mm:d.b.lateral.span};}
 if(a.kind==='failure'){const d=MM.failure87(state.t,p);return {synthetic:true,phase:d.phase,mode:d.mode,material:d.material.name,strain:d.eps,stress_MPa:d.stress,broken:d.broken,unloadResidualStrain:d.currentResidual};}
 if(a.kind==='breath'){const d=MM.breath87(state.t,a.motion);return {synthetic:true,phase:d.phase,mode:d.mode,scope:'illustration geometry, not patient measurements'};}return undefined;}
function snapshot(){pause();if(!state.action)return;if(state.record.observations.length>=1000){toast('达到1000条，请先导出。');return;}const{copData,safety,...p}=state.params;state.record.observations.push({at:new Date().toISOString(),module:state.module.id,action:state.action.id,title:state.action.title,t:state.t,side:state.side,mode:state.mode,actor:state.mode==='clinic'?state.role:'独立教学模型',condition:state.mode==='clinic'?(state.clinic?.seated?'患者坐位，治疗师站立面对患者':'患者站立，治疗师站立面对患者'):state.action.position,clinicState:state.mode==='clinic'?state.clinicPhase:undefined,parameters:p,safety:state.action.kind==='cardio'?JSON.parse(JSON.stringify(safety)):undefined,dataSource:copData?copData.label:state.action.caption,revisionReadout:revisionReadout87(),deviceReadout:state.action.kind==='dynamometer'?(()=>{const d=MM.deviceTrial(state.action.motion,state.t,state.params);return {phase:d.phase,synthetic:true,time:d.time,angle:d.angle,velocity:d.velocity,userTorque:d.userTorque,deviceTorque:d.deviceTorque,force:d.force,peak:d.peak};})():undefined,exampleTime:state.action.kind==='contraction'?MM.contraction(state.action.contractionMode,state.t,state.params).time:undefined});save();toast('已保存这个观察时刻，不计为临床技能得分。');}
function validateRecord(x){
 if(!x||!['kinesiology-observation-8.2','kinesiology-observation-8.4','kinesiology-observation-8.6','kinesiology-observation-8.7','kinesiology-observation-8.8'].includes(x.format)||!Array.isArray(x.observations)||x.observations.length>1000||!Array.isArray(x.explanations)||x.explanations.length>1000||!Array.isArray(x.events)||x.events.length>3000)throw Error('不是支持的观察记录备份，或条数超过上限。');
 for(const r of x.observations){if(!Number.isFinite(r.t)||r.t<0||r.t>1||typeof r.title!=='string')throw Error('记录包含无效阶段或标题。');
  if(r.module===1&&r.action==='01-15'){r.legacyModule=1;r.module=101;}
  const m=course.modules.find(m=>m.id===r.module);if(!m)throw Error('记录包含未知项目。');
  if(r.module===3&&x.format==='kinesiology-observation-8.2'){
   const original=r.action;
   if(/等速/.test(r.title))r.action='03-9';else if(/牵拉|缩短周期/.test(r.title))r.action='03-4';else if(/1RM/.test(r.title))r.action='03-5';else if(/摸高|滞空|腾空|快速发力/.test(r.title))r.action='03-6';else if(/平板/.test(r.title))r.action='03-7';else if(/深蹲/.test(r.title))r.action='03-8';
   if(r.action!==original)r.legacyAction=original;
  }
  if(r.module===11&&r.action==='11-4'){r.legacyArchived=true;}else if(!m.actions.some(a=>a.id===r.action))throw Error('备份包含此版本没有的观察项。');
 }
 for(const r of x.explanations){if(!r||typeof r.answer!=='string'||r.answer.length>1600||typeof r.action!=='string')throw Error('解释记录格式无效。');}
 for(const r of x.events){if(!r||typeof r.type!=='string'||r.type.length>100)throw Error('过程事件格式无效。');}
 if(x.format==='kinesiology-observation-8.2'){
  const known={};for(const r of x.observations)if(r.legacyAction){known[r.legacyAction]=r.action;}
  for(const r of [...x.explanations,...x.events]){r.legacyVersion=x.version||'8.2系列';if(known[r.action]){r.legacyAction=r.action;r.action=known[r.action];}else if(r.action?.startsWith('03-'))r.legacyReferenceUnresolved=true;}
 }
 x.importedFrom=x.importedFrom||x.version;x.format='kinesiology-observation-8.8';x.version='8.8.0';return true;
}
function recordTable(){return `<table><tr><th>项目与观察</th><th>阶段 / 侧别</th><th>来源与方式</th></tr>${state.record.observations.map(r=>`<tr><td>${E(r.action)} ${E(r.title)}</td><td>${Math.round(r.t*100)}% / ${E(r.side)}</td><td>${E(r.mode)}<br>${E(r.dataSource||'')}</td></tr>`).join('')}</table>`;}
function records(){if(!Portal.access)return;modal('我的观察记录',`<p><b>${E(state.profile.name)} · ${E(state.profile.studentId)} · ${E(state.profile.teacher)}</b></p><p>${state.record.observations.length} 个时刻，${state.record.explanations.length} 条解释。没有自动临床评分。</p><div class="row"><button data-cmd="exportHTML">导出阅读报告</button><button data-cmd="exportJSON">导出备份</button><label class="filebutton">恢复备份<input id="backupFile" type="file" accept=".json"></label></div>${recordTable()}<h3>我的解释</h3>${state.record.explanations.map(r=>`<p><b>${E(r.action)}</b> ${r.legacyVersion?'〔旧版 '+E(r.legacyVersion)+(r.legacyReferenceUnresolved?'：项目编号待核对':'')+'〕 ':''}${E(r.answer)}</p>`).join('')}`);}
function about(){modal('动小学 · 人体运动学虚拟实验室 · V8.9',`<p>先登记姓名、学号和任课教师，再选择学习项目。独立肌肉模块包含104项肌肉、分部或肌群资料，其中98项有来源网格、6项仅功能资料；并非全身每块肌肉的完整仿真。研发者：张积众。</p><p>教师姓名只是课堂入口校验，不验证学生真实身份。默认数据仅在此浏览器按学号保存；使用配套校内部署服务时，可在结束页主动提交。没有学习通成绩接口。</p><p>骨肌形态来自既有BodyParts3D资源。新肌肉动画采用同源坐标和可视化附着绑定，不是肌肉力学、精确肌纤维形变或临床评估。</p><p>详细实施范围与来源见单独教师说明。未经真实手机、校内学习通及学生试教验收。</p><details><summary>模型许可</summary>${Object.entries(payload.notices).map(([n,v])=>`<h4>${E(n)}</h4><pre>${E(v)}</pre>`).join('')}</details>`);}

function actionCommand(cmd){if(cmd==='atlasRecord'){snapshot();return;}if(cmd==='atlasFocus'){state.model?.atlasFocus();return;}if(cmd==='atlasAll'){state.model?.fitAll();return;}try{
 if(cmd==='gaitReference87'){const r=MM.gaitReference87;modal('健康成人研究参考 · 与本动画分开',`<p>Orendurff等，2004。样本10名健康成人，15节段方法计算全身COM。下列为周期全幅（最大减最小）的组均值±标准差；原文厘米换算为毫米。</p><table class="source-table"><tr><th>步速</th><th>上下全幅</th><th>左右全幅</th></tr>${r.rows.map(x=>`<tr><td>${x.speed} m/s</td><td>${x.verticalMean.toFixed(1)} ± ${x.verticalSD.toFixed(1)} mm</td><td>${x.lateralMean.toFixed(1)} ± ${x.lateralSD.toFixed(1)} mm</td></tr>`).join('')}</table><p>这不是正常值区间或诊断阈值；样本、步速与算法都会影响结果。本动画没有导入该研究原始轨迹，也没有用这些均值缩放人体或曲线。</p>${sourceHTML(['V87_GAIT_STUDY'])}`);}

 if(cmd==='isolate'&&!state.model?.isolate())toast('先点击模型中的一个结构。');if(cmd==='restore')state.model?.restore();
 if(cmd==='calcHR'){const v=['restHR','maxHR','intensity'].map(id=>{if(!$(id).value.trim())throw Error('请填写所有值。');return MM.finite($(id).value);});const[r,m,i]=v;if(r<=0||m<=r||i<0||i>100)throw Error('最大值须大于静息值，强度须在0—100%。');$('calcResult').textContent=`例题靶心率 ${(r+(m-r)*i/100).toFixed(1)} 次/分；仅为给定数据计算。`;}
 if(cmd==='calcArch'){$('calcResult').textContent='面积比 '+Sim.Data.index(...['rearArea','midArea','foreArea'].map(id=>{if(!$(id).value.trim())throw Error('请填写三分区面积。');return MM.finite($(id).value);})).toFixed(3)+'；不作自动诊断。';}
 if(cmd==='copExample')download('cop_synthetic_template.csv','text/csv;charset=utf-8','time_s,x_mm,y_mm,fz_N\n'+Sim.syntheticCOP().rows.map(r=>[r.time_s,r.x_mm,r.y_mm,r.fz_N].join(',')).join('\n'));
 if(cmd==='copReset'){delete state.params.copData;state.t=0;update();}
 if(cmd==='exportJSON')download('kinesiology_v88_backup.json','application/json',JSON.stringify(state.record,null,2));
 if(cmd==='exportHTML')download('dongxiaoxue_v89_report.html','text/html;charset=utf-8',`<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>人体运动学观察记录</title><style>body{font:16px/1.8 system-ui;max-width:1000px;margin:30px auto;padding:15px}td,th{border-bottom:1px solid #ccc;padding:10px;text-align:left}table{width:100%}</style><h1>动小学 · 人体运动学V8.9观察记录</h1><p><b>${E(state.profile.name)} · 学号 ${E(state.profile.studentId)} · 任课教师 ${E(state.profile.teacher)}</b></p><p>不包含自动技能评分，不替代真实患者测量。</p>${recordTable()}<h2>自己的解释</h2>${state.record.explanations.map(r=>`<p>${E(r.action)}：${E(r.answer)}</p>`).join('')}<h2>过程事件</h2>${state.record.events.map(e=>`<p>${E(e.at)} — ${E(e.type)} — ${E(e.action||'')}</p>`).join('')}</html>`);
 }catch(e){if($('calcResult'))$('calcResult').textContent=e.message;else toast(e.message);}}
async function importCOP(file){try{if(file.size>1000000)throw Error('文件超过1MB。');const d=Sim.Data.parseCOP(await file.text(),file.name);state.params.copData=d;pause();state.t=0;update();toast('已读取。来源、坐标和测试条件仍需核定。');}catch(e){toast(e.message);}}
async function importBackup(file){try{if(file.size>4000000)throw Error('文件超过4MB。');const x=JSON.parse(await file.text());validateRecord(x);if(x.profile&&(x.profile.studentId!==state.profile.studentId||x.profile.teacher!==state.profile.teacher))throw Error('只能恢复当前学号的备份。');x.profile=state.profile;if(!confirm('将替换当前观察记录，继续吗？'))return;state.record=x;save();records();}catch(e){toast(e.message);}}
$('homeLink').onclick=e=>{e.preventDefault();home();};$('backBtn').onclick=home;$('courseGrid').onclick=e=>{const b=e.target.closest('[data-module]');if(b)openLesson(b.dataset.module);};$('actionSelect').onchange=e=>selectAction(e.target.value);$('modes').onclick=e=>{const b=e.target.closest('[data-mode]');if(b)switchMode(b.dataset.mode);};$('recordsBtn').onclick=records;$('aboutBtn').onclick=about;$('closeModal').onclick=()=>$('modal').close();
$('playBtn').onclick=play;$('prevFrame').onclick=()=>setStage(state.t-1/60);$('nextFrame').onclick=()=>setStage(state.t+1/60);$('resetBtn').onclick=()=>setStage(state.action?.kind==='fsu'?.5:0);$('stageRange').oninput=e=>setStage(Number(e.target.value)/1000);$('snapshotBtn').onclick=snapshot;
$('viewButtons').onclick=e=>{const b=e.target.closest('[data-view]');if(!b)return;const v=b.dataset.view;if(state.mode==='clinic'){v==='fit'?state.clinic.frame('both'):state.clinic.orient(v);}else v==='fit'?state.model.fitAll():state.model.orient(v);};
$('nextLesson').onclick=()=>{const i=course.modules.indexOf(state.module);if(i===course.modules.length-1)Portal.finish();else openLesson(course.modules[i+1].id);};
$('askMuscles').onclick=()=>{$('recallBody').hidden=!$('recallBody').hidden;};$('showAnswer').onclick=()=>{$('muscleAnswer').innerHTML=muscleHTML();$('muscleAnswer').hidden=false;event('reference_revealed');};$('saveAnswer').onclick=()=>{const v=$('answer').value.trim();if(!v){toast('先写下你的解释。');return;}if(state.record.explanations.length<1000)state.record.explanations.push({at:new Date().toISOString(),action:state.action.id,answer:v,side:state.side,t:state.t});save();toast('解释已保存，未自动评分。');};
document.addEventListener('click',e=>{const dx=e.target.closest('[data-device]');if(dx)deviceCommand(dx.dataset.device);const ds=e.target.closest('[data-stage]');if(ds)setStage(Number(ds.dataset.stage));const cm=e.target.closest('[data-catalog-module]');if(cm){$('modal').close();openLesson(cm.dataset.catalogModule);}const cmd=e.target.closest('[data-cmd]');if(cmd)actionCommand(cmd.dataset.cmd);const r=e.target.closest('[data-risk]');if(r)risk(r.dataset.risk);const c=e.target.closest('[data-clinic]');if(c)clinicCommand(c.dataset.clinic);const g=e.target.closest('[data-gait-event]');if(g)setStage(Number(g.dataset.gaitEvent));});
document.addEventListener('click',e=>{const pick=e.target.closest('[data-atlas-pick]');if(pick){selectAction(pick.dataset.atlasPick);return;}if(e.target.closest('[data-atlas-reset]'))filterMuscles('reset');});
document.addEventListener('change',e=>{const x=e.target;if(x.id==='atlasGroup')filterMuscles('group');if(x.dataset.deviceCheck!==undefined)state.params.deviceChecks[Number(x.dataset.deviceCheck)]=x.checked;if(x.id==='sideSelect'){pause();if(state.mode==='clinic'&&['stopped','closed'].includes(state.clinicPhase)){x.value=state.side;toast('本次患者观察已结束，不能通过换侧继续。');return;}state.side=x.value;if(state.mode==='clinic'){state.clinic.side=state.side;state.clinic.setPose('patient','rest',0);state.clinic.setPose('therapist','rest',0);state.t=0;}else state.model.set(liveAction(),state.t,'native');event('side_change',{side:state.side});update();}
 if(x.dataset.check!==undefined){state.params.safety.checked[Number(x.dataset.check)]=x.checked;}
 if(x.id==='copFile'&&x.files[0])importCOP(x.files[0]);if(x.id==='backupFile'&&x.files[0])importBackup(x.files[0]);});
$('extras').addEventListener('input',e=>{const x=e.target,k=x.dataset.setting;if(x.id==='atlasSearch'){filterMuscles();return;}if(!k)return;pause();if(k==='layer'){state.model.setLayer(x.value);return;}let v=x.type==='checkbox'?x.checked:['load','distance','mass','takeoff','airDistance','airLevel','airSeconds','stiffness','isoSpeed','effort','deviceEffort','caseAge','caseRestHR','targetPct'].includes(k)?Number(x.value):x.value;
 if(typeof v==='number'&&(!x.value.trim()||!Number.isFinite(v)||x.min!==''&&v<Number(x.min)||x.max!==''&&v>Number(x.max))){toast('请使用输入框所示范围内的有效数值。');x.setAttribute('aria-invalid','true');return;}x.removeAttribute('aria-invalid');state.params[k]=v;if(k==='failureMode'||(k==='material'&&state.action?.kind==='failure'))state.t=0;if(k==='muscleSide')state.side=v;if(k==='supportSide'){state.model.action.supportSide=v;state.model.fitAll();}if(state.action?.kind==='dynamometer'&&['deviceSide','isoSpeed','isoDirection','mvcEffort'].includes(k)){state.params.deviceStarted=false;state.params.deviceStopped=false;state.t=0;renderExtras();}update();});
document.addEventListener('visibilitychange',()=>{if(document.hidden)pause();});window.addEventListener('pagehide',pause);let resizeTimer;window.addEventListener('resize',()=>{clearTimeout(resizeTimer);resizeTimer=setTimeout(()=>{if(nativeActive()&&state.action)state.model?.fitAll();if(state.mode==='clinic')state.clinic?.frame('both');},150);});
function deviceCommand(cmd){pause();const p=state.params;if(state.action?.kind!=='dynamometer')return;
 if(cmd==='new'){p.deviceStarted=false;p.deviceStopped=false;p.deviceChecks=[false,false,false];state.t=0;}
 else if(cmd==='start'){if(!p.deviceChecks.every(Boolean)){toast('请先核对轴线、体位固定与模拟校零。');return;}if(p.deviceStopped){toast('当前试次已结束，请先新建试次。');return;}p.deviceStarted=true;state.t=0;state.dir=1;event('device_start',{mode:state.action.motion});renderExtras();play();}
 else if(cmd==='stop'){p.deviceStopped=true;event('device_stop',{stage:state.t});}
 renderExtras();update();}
$('catalogBtn').onclick=()=>Portal.access&&modal('选择学习项目',`<div class="dialogmenu">${course.modules.map(m=>`<button data-catalog-module="${m.id}"><span>${(m.displayCode||String(m.id).padStart(2,'0'))}</span>${E(m.title)}</button>`).join('')}</div>`);
$('courseSearch').oninput=cards;$('filters').onclick=e=>{const b=e.target.closest('[data-filter]');if(!b)return;courseFilter=b.dataset.filter;for(const x of $('filters').querySelectorAll('button'))x.setAttribute('aria-pressed',x===b?'true':'false');cards();};
$('focusBtn').onclick=()=>{const on=document.body.classList.toggle('focus-mode');$('focusBtn').textContent=on?'退出放大':'放大观察';$('focusBtn').setAttribute('aria-pressed',on?'true':'false');setTimeout(()=>{if(nativeActive())state.model?.fitAll();if(state.mode==='clinic')state.clinic?.frame('both');},60);};
document.addEventListener('keydown',e=>{if(e.key==='Escape'&&document.body.classList.contains('focus-mode'))$('focusBtn').click();});
window.App={state,course,payload,openLesson,selectAction,switchMode,setStage,deviceCommand,play,pause,risk,clinicCommand,sourceHTML,actionCommand,validateRecord,importCOP,importBackup,snapshot,home,renderExtras,update,save,records};
cards();Portal.init();
})();
