/* BoilFlux landing page: a live, real-time jet quench.
 *
 * Explicit finite-volume conduction in the length x thickness section of a
 * 120 mm x 5 mm AA6082 plate (constant properties), plate frame. A water jet
 * sweeps along the quenched face; its local heat flux follows the zone-blended
 * Specht boiling curve of boilflux.synthetic.JetQuenchScenario, evaluated at the
 * surface-cell temperature. Top: surface heat flux q(x). Middle: temperature
 * field (thickness exaggerated). Bottom: the IR-face temperature, the data an
 * inverse solver works from.
 */
(function () {
  "use strict";

  // ---- physics ---------------------------------------------------------
  const NX = 120, NZ = 12, LX = 0.12, LZ = 0.005;
  const DX = LX / NX, DZ = LZ / NZ;
  const K = 180, RHOC = 2.5e6, ALPHA = K / RHOC;
  const DT = 0.9 / (2 * ALPHA * (1 / (DX * DX) + 1 / (DZ * DZ)));
  const T0 = 500, TF = 20, QDOT = 100;
  const V_JET = 0.008, X_START = 0.025, T_END = 10.5;
  const R_IZ = 0.009, R_FFZ = 0.011, R_WET = 0.035, WET_EDGE = 0.004;
  const IZ = { qmax: 6.5e6, tdnb: 200, qleid: 2.0e6, tleid: 450 };
  const FFZ = { qmax: 3.5e6, tdnb: 200, qleid: 1.0e6, tleid: 350 };

  const fc = (Ts) => (14.6 * (Ts + 273) + 68.5 * (TF + 273) - 21500) * Math.pow(QDOT, 0.33);
  const Q100 = fc(100) * (100 - TF);

  function specht(c, Ts) {
    if (Ts <= 100) return fc(Ts) * (Ts - TF);
    if (Ts <= c.tdnb) { const u = (c.tdnb - Ts) / (c.tdnb - 100); return c.qmax - u * u * (c.qmax - Q100); }
    if (Ts <= c.tleid) { const u = (Ts - c.tdnb) / (c.tleid - c.tdnb); return c.qmax - u * u * (c.qmax - c.qleid); }
    return c.qleid;
  }

  function createSim() {
    const T = new Float64Array(NX * NZ).fill(T0);
    const Tn = new Float64Array(NX * NZ);
    const q = new Float64Array(NX);
    let t = 0;
    const ax = ALPHA * DT / (DX * DX), az = ALPHA * DT / (DZ * DZ), bc = DT / (RHOC * DZ);

    function step() {
      const xs = X_START + V_JET * t;
      for (let i = 0; i < NX; i++) {
        const r = Math.abs((i + 0.5) * DX - xs);
        const w = Math.min(1, Math.max(0, (R_FFZ - r) / (R_FFZ - R_IZ)));
        const chi = 0.5 * (1 - Math.tanh((r - R_WET) / WET_EDGE));
        const Ts = T[i];
        q[i] = chi * (w * specht(IZ, Ts) + (1 - w) * specht(FFZ, Ts));
      }
      for (let k = 0; k < NZ; k++) {
        for (let i = 0; i < NX; i++) {
          const p = k * NX + i, Tp = T[p];
          let d = 0;
          if (i > 0) d += ax * (T[p - 1] - Tp);
          if (i < NX - 1) d += ax * (T[p + 1] - Tp);
          if (k > 0) d += az * (T[p - NX] - Tp);
          if (k < NZ - 1) d += az * (T[p + NX] - Tp);
          if (k === 0) d -= bc * q[i];
          Tn[p] = Tp + d;
        }
      }
      T.set(Tn);
      t += DT;
    }
    return { T, q, step, get t() { return t; }, reset() { T.fill(T0); q.fill(0); t = 0; } };
  }

  // ---- colour ----------------------------------------------------------
  // inferno, sampled at 11 points (perceptually uniform, standard in thermography)
  const INFERNO = ["#000004", "#160b39", "#420a68", "#6a176e", "#932667", "#bc3754",
                   "#dd513a", "#f37819", "#fca50a", "#f6d746", "#fcffa4"].map((h) =>
    [1, 3, 5].map((j) => parseInt(h.slice(j, j + 2), 16)));
  function inferno(s) {
    s = Math.min(1, Math.max(0, s)) * 10;
    const i = Math.min(9, Math.floor(s)), f = s - i, a = INFERNO[i], b = INFERNO[i + 1];
    return [a[0] + f * (b[0] - a[0]), a[1] + f * (b[1] - a[1]), a[2] + f * (b[2] - a[2])];
  }
  const cssVar = (el, name, fb) => (getComputedStyle(el).getPropertyValue(name).trim() || fb);

  // ---- drawing ---------------------------------------------------------
  function init(canvas) {
    if (canvas.dataset.bfInit) return;
    canvas.dataset.bfInit = "1";
    const ctx = canvas.getContext("2d");
    const sim = createSim();
    const field = document.createElement("canvas");
    field.width = NX; field.height = NZ;
    const fctx = field.getContext("2d");
    const img = fctx.createImageData(NX, NZ);
    const reduce = window.matchMedia("(prefers-reduced-motion: reduce)").matches;
    const stepsPerFrame = Math.max(1, Math.round((1 / 60) / DT));
    let W = 0, H = 0, visible = true, raf = 0, hold = 0;

    function resize() {
      const dpr = Math.min(2, window.devicePixelRatio || 1);
      W = canvas.clientWidth; H = canvas.clientHeight;
      canvas.width = Math.round(W * dpr); canvas.height = Math.round(H * dpr);
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    }

    function draw() {
      const fg = cssVar(canvas, "--bf-ink", "#172331");
      const muted = cssVar(canvas, "--bf-muted", "#6b7685");
      const grid = cssVar(canvas, "--bf-rule", "#d9dee5");
      const heat = cssVar(canvas, "--bf-heat", "#b8462a");
      const blue = cssVar(canvas, "--bf-blue", "#1f5fa6");
      const jet = cssVar(canvas, "--bf-jet", "#6aa9d8");
      ctx.clearRect(0, 0, W, H);
      ctx.font = "500 11px Inter, system-ui, sans-serif";
      ctx.textBaseline = "alphabetic";

      const padL = 8, padR = 8, x0 = padL, pw = W - padL - padR;
      const xOf = (i) => x0 + ((i + 0.5) / NX) * pw;
      const qTop = 22, qH = H * 0.24;
      const plateY = qTop + qH + 34, plateH = Math.max(34, H * 0.16);
      const tTop = plateY + plateH + 44, tH = H - tTop - 14;

      // q(x) panel
      ctx.strokeStyle = grid; ctx.lineWidth = 1;
      ctx.beginPath(); ctx.moveTo(x0, qTop + qH + 0.5); ctx.lineTo(x0 + pw, qTop + qH + 0.5); ctx.stroke();
      const narrow = pw < 380;
      ctx.fillStyle = heat; ctx.fillText(narrow ? "surface flux q(x)" : "surface heat flux  q(x)", x0, qTop - 8);
      ctx.fillStyle = muted; ctx.textAlign = "right";
      ctx.fillText("0 – 7 MW/m²", x0 + pw, qTop - 8); ctx.textAlign = "left";
      ctx.beginPath();
      for (let i = 0; i < NX; i++) {
        const y = qTop + qH - Math.min(1, sim.q[i] / 7e6) * qH;
        i ? ctx.lineTo(xOf(i), y) : ctx.moveTo(xOf(i), y);
      }
      ctx.strokeStyle = heat; ctx.lineWidth = 1.8; ctx.stroke();

      // jet
      const xs = X_START + V_JET * sim.t;
      const xj = x0 + (xs / LX) * pw;
      ctx.fillStyle = muted; ctx.fillRect(xj - 9, plateY - 30, 18, 7);
      ctx.fillStyle = jet; ctx.globalAlpha = 0.75;
      ctx.beginPath(); ctx.moveTo(xj - 4, plateY - 23); ctx.lineTo(xj + 4, plateY - 23);
      ctx.lineTo(xj + 7, plateY); ctx.lineTo(xj - 7, plateY); ctx.closePath(); ctx.fill();
      ctx.globalAlpha = 1;
      ctx.fillStyle = muted; ctx.fillText("water jet", xj + 14, plateY - 18);

      // temperature field
      for (let k = 0; k < NZ; k++) for (let i = 0; i < NX; i++) {
        const c = inferno((sim.T[k * NX + i] - TF) / (T0 - TF)), o = 4 * (k * NX + i);
        img.data[o] = c[0]; img.data[o + 1] = c[1]; img.data[o + 2] = c[2]; img.data[o + 3] = 255;
      }
      fctx.putImageData(img, 0, 0);
      ctx.imageSmoothingEnabled = true;
      ctx.drawImage(field, x0, plateY, pw, plateH);
      ctx.strokeStyle = fg; ctx.globalAlpha = 0.25; ctx.strokeRect(x0 + 0.5, plateY + 0.5, pw - 1, plateH - 1);
      ctx.globalAlpha = 1;

      // colour bar
      const cbW = Math.min(110, pw * 0.28), cbX = x0 + pw - cbW - 46, cbY = plateY + plateH + 8;
      for (let j = 0; j < cbW; j++) {
        const c = inferno(j / cbW); ctx.fillStyle = `rgb(${c[0] | 0},${c[1] | 0},${c[2] | 0})`;
        ctx.fillRect(cbX + j, cbY, 1.2, 6);
      }
      ctx.fillStyle = muted; ctx.font = "500 10px Inter, system-ui, sans-serif";
      ctx.textAlign = "right"; ctx.fillText("20", cbX - 4, cbY + 7);
      ctx.textAlign = "left"; ctx.fillText(`${T0} °C`, cbX + cbW + 4, cbY + 7);
      ctx.font = "500 11px Inter, system-ui, sans-serif";
      ctx.fillText(`t = ${sim.t.toFixed(1)} s`, x0, cbY + 8);

      // IR-face temperature panel
      ctx.strokeStyle = grid; ctx.lineWidth = 1;
      ctx.beginPath(); ctx.moveTo(x0, tTop + tH + 0.5); ctx.lineTo(x0 + pw, tTop + tH + 0.5); ctx.stroke();
      ctx.fillStyle = blue; ctx.fillText(narrow ? "IR face T(x)" : "IR camera: back-face temperature  T(x)", x0, tTop + 2);
      ctx.fillStyle = muted; ctx.textAlign = "right"; ctx.fillText("0 – 500 °C", x0 + pw, tTop + 2);
      ctx.textAlign = "left";
      ctx.beginPath();
      for (let i = 0; i < NX; i++) {
        const v = sim.T[(NZ - 1) * NX + i];
        const y = tTop + tH - Math.min(1, Math.max(0, v / T0)) * (tH - 10);
        i ? ctx.lineTo(xOf(i), y) : ctx.moveTo(xOf(i), y);
      }
      ctx.strokeStyle = blue; ctx.lineWidth = 1.8; ctx.stroke();
    }

    function frame() {
      raf = 0;
      if (!visible || document.hidden) return;
      if (hold > 0) { hold--; } else {
        for (let s = 0; s < stepsPerFrame; s++) sim.step();
        if (sim.t > T_END) { sim.reset(); hold = 45; }
      }
      draw();
      raf = requestAnimationFrame(frame);
    }
    const start = () => { if (!raf && !reduce) raf = requestAnimationFrame(frame); };

    resize();
    if (reduce) {
      while (sim.t < 3.0) sim.step();
      draw();
    } else {
      start();
    }
    new ResizeObserver(() => { resize(); draw(); }).observe(canvas);
    new IntersectionObserver((es) => { visible = es[0].isIntersecting; if (visible) start(); }).observe(canvas);
    document.addEventListener("visibilitychange", start);
    new MutationObserver(draw).observe(document.body, { attributes: true, attributeFilter: ["data-md-color-scheme"] });
  }

  // ---- reveal-on-scroll for cards --------------------------------------
  function reveal() {
    if (window.matchMedia("(prefers-reduced-motion: reduce)").matches || !("IntersectionObserver" in window)) return;
    const els = document.querySelectorAll(".bf-card, .bf-step, .bf-shot, .bf-part, .bf-figure");
    if (!els.length) return;
    document.documentElement.classList.add("bf-anim");
    const io = new IntersectionObserver((es) => es.forEach((e) => {
      if (e.isIntersecting) { e.target.classList.add("is-visible"); io.unobserve(e.target); }
    }), { rootMargin: "0px 0px -8% 0px" });
    els.forEach((el, i) => { el.classList.add("bf-reveal"); el.style.transitionDelay = `${(i % 4) * 60}ms`; io.observe(el); });
  }

  function boot() {
    document.querySelectorAll("canvas.bf-sim-canvas").forEach(init);
    reveal();
  }
  if (window.document$ && typeof window.document$.subscribe === "function") window.document$.subscribe(boot);
  else if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", boot);
  else boot();
})();
