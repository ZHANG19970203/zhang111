/* Dongxiaoxue V8.9 static-deployment loader. Packaging only; no anatomical changes. */
(function () {
  'use strict';
  const root = new URL('./', document.baseURI);
  const diag = window.DXX_DEPLOY_DIAGNOSTICS = {
    state: 'loading', loadedBytes: 0, modelHashes: {}, checks: [], backend: false
  };
  const el = id => document.getElementById(id);
  const canHash = !!(globalThis.crypto && crypto.subtle);
  const mib = value => (value / 1048576).toFixed(1);
  let total = 1;
  function status(message) { el('dxxBootStatus').textContent = message; }
  function progress(n) {
    diag.loadedBytes += n;
    const value = Math.min(100, diag.loadedBytes / total * 100);
    el('dxxBootProgress').value = value;
    el('dxxBootBytes').textContent = mib(diag.loadedBytes) + ' / ' + mib(total) + ' MiB';
  }
  async function hash(buffer) {
    if (!canHash) return null;
    const digest = await crypto.subtle.digest('SHA-256', buffer);
    return Array.from(new Uint8Array(digest), b => b.toString(16).padStart(2, '0')).join('');
  }
  function checkedURL(path) {
    const url = new URL(path, root);
    if (url.origin !== root.origin || !url.pathname.startsWith(root.pathname)) {
      throw new Error('资源路径不属于当前部署目录：' + path);
    }
    return url;
  }
  async function download(path, size, sha) {
    let last;
    for (let attempt = 0; attempt < 2; attempt++) {
      const controller = new AbortController();
      const timeout = setTimeout(() => controller.abort(), 120000);
      try {
        const response = await fetch(checkedURL(path), {
          signal: controller.signal, credentials: 'same-origin', cache: 'default'
        });
        if (!response.ok) {
          const e = new Error('资源加载失败（HTTP ' + response.status + '）：' + path);
          e.noRetry = response.status === 404;
          throw e;
        }
        const data = await response.arrayBuffer();
        if (size != null && data.byteLength !== size) throw new Error('资源大小不完整：' + path);
        const actual = sha && canHash ? await hash(data) : null;
        if (actual && actual !== sha) throw new Error('资源校验不一致，请重新上传：' + path);
        return data;
      } catch (e) {
        last = e;
        if (e.noRetry) break;
        if (!attempt) await new Promise(r => setTimeout(r, 700));
      } finally { clearTimeout(timeout); }
    }
    throw last;
  }
  async function runApp(spec) {
    status('资源已就绪，正在启动实验室…');
    return new Promise((resolve, reject) => {
      let runtimeError = null;
      function capture(event) { runtimeError = event.error || new Error(event.message); }
      window.addEventListener('error', capture);
      const script = document.createElement('script');
      script.src = checkedURL(spec.path).href;
      script.integrity = spec.integrity;
      script.crossOrigin = 'anonymous';
      script.onload = () => {
        window.removeEventListener('error', capture);
        if (runtimeError || !window.App || !window.Portal) {
          reject(runtimeError || new Error('主程序没有完成初始化。'));
        } else resolve();
      };
      script.onerror = () => {
        window.removeEventListener('error', capture);
        reject(new Error('主程序加载失败，请核对 assets/app.js 是否完整上传。'));
      };
      document.body.appendChild(script);
    });
  }
  async function start() {
    if (location.protocol === 'file:') {
      throw new Error('这是分文件部署版，请先上传整个文件夹，再打开网站的 HTTPS 地址。本地直接双击请使用原来的 dongxiaoxue_v89_lab.html。');
    }
    const manifestBytes = await download('assets/deploy-manifest.json');
    const manifest = JSON.parse(new TextDecoder('utf-8').decode(manifestBytes));
    if (manifest.format !== 'dxx-static-split-v1') throw new Error('部署清单格式不匹配。');
    diag.sourceSHA256 = manifest.sourceSHA256;
    diag.integrity = canHash ? 'SHA-256 + byte lengths' : 'byte lengths only';
    total = manifest.metadata.bytes + manifest.models.reduce((s, a) => s + a.bytes, 0);
    status('正在加载课程、目录与模型索引…');
    const metaBytes = await download(manifest.metadata.path, manifest.metadata.bytes, manifest.metadata.sha256);
    const payload = JSON.parse(new TextDecoder('utf-8').decode(metaBytes));
    progress(metaBytes.byteLength);
    const buffers = {};
    const jobs = [];
    for (const model of manifest.models) {
      buffers[model.key] = new Uint8Array(model.bytes);
      let cursor = 0;
      for (const chunk of model.chunks) {
        if (chunk.offset !== cursor || chunk.offset + chunk.bytes > model.bytes) {
          throw new Error('模型分片顺序不完整：' + model.key);
        }
        cursor += chunk.bytes;
        jobs.push({key: model.key, chunk});
      }
      if (cursor !== model.bytes) throw new Error('模型分片缺失：' + model.key);
    }
    let next = 0, completed = 0;
    async function worker() {
      while (next < jobs.length) {
        const {key, chunk} = jobs[next++];
        const data = await download(chunk.path, chunk.bytes, chunk.sha256);
        buffers[key].set(new Uint8Array(data), chunk.offset);
        progress(data.byteLength);
        completed++;
        status('正在准备完整模型：' + completed + ' / ' + jobs.length + ' 个资源分片');
      }
    }
    await Promise.all(Array.from({length: Math.min(3, jobs.length)}, worker));
    for (const model of manifest.models) {
      const data = buffers[model.key].buffer;
      const actual = await hash(data);
      if (actual && actual !== model.sha256) throw new Error('合并后的模型校验失败：' + model.key);
      payload[model.key].raw = data;
      diag.modelHashes[model.key] = actual || 'not available in this browser';
    }
    // Existing synchronous readers can use the original bytes without base64 decoding.
    window.DXX_ENABLE_CLASSROOM_SERVICE = false;
    window.__DXX_DATA__ = payload;
    await runApp(manifest.application);
    diag.state = 'ready';
    diag.checks.push('all assets loaded', 'model buffers reconstructed', 'App and Portal initialized');
    el('dxxBootProgress').value = 100;
    document.body.classList.remove('dxx-loading');
    el('dxxBoot').remove();
  }
  start().catch(error => {
    diag.state = 'error';
    diag.error = error && error.message ? error.message : String(error);
    status('实验室尚未载入完成');
    el('dxxBootError').textContent = diag.error;
    el('dxxBootError').hidden = false;
    el('dxxBootRetry').hidden = false;
    el('dxxBootRetry').onclick = () => location.reload();
    el('dxxBootHelp').hidden = false;
    console.error('[动小学部署加载]', diag.error);
  });
})();
