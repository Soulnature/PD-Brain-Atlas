/* PD Brain Atlas viewer. Data: window.PD_ATLAS (built by biomni.pd.atlas.site). All remote-derived text is escaped. */
(function () {
  "use strict";
  const D = window.PD_ATLAS;
  const $ = (id) => document.getElementById(id);
  const esc = (s) => String(s ?? "").replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
  const REG = Object.fromEntries(D.regions.map((r) => [r.id, r]));
  const LIT = D.literature;
  const molKey = (m) => String(m || "").trim().toUpperCase();

  const DIVERGING = ["#2a78d6", "#86b6ef", "#e8e7e3", "#f4a37c", "#d9481c"];
  const SEQ = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"];
  const DIR_COLOR = { increased: "#d9481c", decreased: "#2a78d6", "no change": "#a8a7a2", present: "#a8a7a2", altered: "#a8a7a2", mixed: "#5b5a56" };
  const NODATA_CTX = "#ebeae6", NODATA_MESH = "#d4d3ce";

  // ---------- indices ----------
  const litByRegion = {};
  const litByMol = {};
  LIT.forEach((f, i) => {
    f._i = i;
    f.regions.forEach((r) => (litByRegion[r] = litByRegion[r] || []).push(f));
    (litByMol[molKey(f.molecule)] = litByMol[molKey(f.molecule)] || []).push(f);
  });
  // PD map (LCSB, CC BY 4.0) mechanism index: molecule -> pathways / compartments / reactions / PMIDs / partners
  const PM = D.pdmap || {};
  const PMM = PM.molecules || {};
  const pmGet = (m) => PMM[molKey(m)] || PMM[String(m || "").trim()] || null;
  const PATHWAYS = Object.keys(PM.pathways || {}).sort();
  // ParkDB (Brain 2026): prefrontal bulk RNA-seq, contrasts "D1:PD" -> [log2FC, padj]
  const PK = D.parkdb || {};
  const pkGet = (m) => (PK.genes || {})[molKey(m)] || null;
  // PD genetics: GWAS Locus Browser (IPDGC/GP2) + GWAS Catalog
  const GN = D.genetics || {};
  const gnGet = (m) => (GN.genes || {})[molKey(m)] || null;
  const isRisk = (g) => g && (g.m || g.n || (g.c && g.c[0] > 0));
  const PFC = ["superiorfrontal", "rostralmiddlefrontal", "caudalmiddlefrontal"].flatMap((p) => [`ctx_L_${p}`, `ctx_R_${p}`]);
  const evByRegion = {};
  D.evidence.forEach((e) => (evByRegion[e.region_id] = evByRegion[e.region_id] || []).push(e));

  const state = { map: "lit_count", molecule: null, selection: null };

  // ---------- colour ----------
  const hexToRgb = (h) => [1, 3, 5].map((i) => parseInt(h.slice(i, i + 2), 16));
  const rgbToHex = (c) => "#" + c.map((v) => Math.round(v).toString(16).padStart(2, "0")).join("");
  function ramp(stops, t) {
    t = Math.max(0, Math.min(1, t));
    const x = t * (stops.length - 1), i = Math.min(Math.floor(x), stops.length - 2), f = x - i;
    const a = hexToRgb(stops[i]), b = hexToRgb(stops[i + 1]);
    return rgbToHex(a.map((v, k) => v + (b[k] - v) * f));
  }
  function pathwayValues(pw) {
    // distinct pathway molecules with a finding in each region (not finding counts, so SNCA cannot dominate)
    const sets = {};
    LIT.forEach((f) => { const m = pmGet(f.molecule); if (m && m.pw.includes(pw)) f.regions.forEach((r) => (sets[r] = sets[r] || new Set()).add(molKey(f.molecule))); });
    return Object.fromEntries(Object.entries(sets).map(([r, s]) => [r, s.size]));
  }
  function mapSpec() {
    if (state.map === "genetics") {
      const sets = {};
      LIT.forEach((f) => { if (isRisk(gnGet(f.molecule))) f.regions.forEach((r) => (sets[r] = sets[r] || new Set()).add(molKey(f.molecule))); });
      const values = Object.fromEntries(Object.entries(sets).map(([r, s]) => [r, s.size]));
      return { kind: "count", label: "PD genetic risk genes with published findings", values, max: Math.max(1, ...Object.values(values)) };
    }
    if (state.map === "pathway") {
      const values = pathwayValues(state.pathway);
      return { kind: "count", label: `PD map ${state.pathway}: pathway molecules with findings`, values, max: Math.max(1, ...Object.values(values)) };
    }
    if (state.map === "molecule") return { kind: "molecule", label: `${state.molecule}: reported direction in PD` };
    const m = D.maps[state.map];
    const vals = Object.values(m.values);
    return { ...m, max: Math.max(1e-9, ...vals.map((v) => Math.abs(v))) };
  }
  // exact: the paper names this structure; broad: a lobe / "cortex" / "brainstem" spread over many atlas regions
  function molDirections(rids) {
    const exact = new Set(), broad = new Set();
    (litByMol[molKey(state.molecule)] || []).forEach((f) => {
      if (f.regions.some((r) => rids.includes(r))) (f.approximate_region ? broad : exact).add(f.direction);
    });
    const pk = pkGet(state.molecule);  // replicated ParkDB PD change -> prefrontal parcels (broad: lobe-level tissue)
    if (pk && pk.rep && rids.some((r) => PFC.includes(r))) broad.add((pk.c["D1:PD"] || [0])[0] > 0 ? "increased" : "decreased");
    return { exact, broad };
  }
  const tint = (hex, w) => rgbToHex(hexToRgb(hex).map((v) => Math.max(0, w >= 0 ? v + (255 - v) * w : v * (1 + w))));
  function summarizeDirs(s) {
    if (!s.size) return null;
    if (s.has("increased") && s.has("decreased")) return "mixed";
    for (const d of ["increased", "decreased", "altered", "present", "no change"]) if (s.has(d)) return d;
    return null;
  }
  // value + colour for a set of atlas regions (a mesh may stand for several regions)
  function regionsValue(rids, spec) {
    if (spec.kind === "molecule") {
      const { exact, broad } = molDirections(rids);
      const d = summarizeDirs(exact);
      if (d) return { color: DIR_COLOR[d], label: d };
      const b = summarizeDirs(broad);
      return b ? { color: tint(DIR_COLOR[b], 0.6), label: `${b} (reported for a broader area)` } : { color: null, label: "not reported" };
    }
    const vs = rids.map((r) => spec.values[r]).filter((v) => v !== undefined);
    if (!vs.length) return { color: null, label: "no data" };
    if (spec.kind === "count") {
      const v = Math.max(...vs);
      return { color: ramp(SEQ, Math.sqrt(v / spec.max)), label: state.map === "pathway" || state.map === "genetics" ? `${v} ${state.map === "genetics" ? "risk gene" : "pathway molecule"}${v === 1 ? "" : "s"} reported` : `${v} finding${v === 1 ? "" : "s"}` };
    }
    const v = vs[0];
    if (spec.kind === "ordinal") return { color: SEQ[Math.round(v) - 1], label: `Braak R${v}` };
    return { color: ramp(DIVERGING, 0.5 + v / (2 * spec.max)), label: `d = ${v.toFixed(3)}` };
  }

  // ---------- traces ----------
  const traces = [];
  const T = { context: null, cortex: [], meshes: [], markers: null };
  T.context = traces.length;
  traces.push({ type: "mesh3d", ...D.context, color: "#b8b7b2", opacity: 0.07, hoverinfo: "skip", flatshading: false,
    lighting: { ambient: 0.9, diffuse: 0.2 }, name: "outline", showlegend: false });
  D.cortex.forEach((c) => {
    T.cortex.push(traces.length);
    traces.push({ type: "mesh3d", x: c.x, y: c.y, z: c.z, i: c.i, j: c.j, k: c.k, vertexcolor: c.x.map(() => NODATA_CTX),
      opacity: 0.18, text: c.vregion.map(() => ""), hovertemplate: "%{text}<extra></extra>", name: c.id,
      lighting: { ambient: 0.55, diffuse: 0.7, specular: 0.05, roughness: 0.9 }, lightposition: { x: 100, y: 200, z: 150 } });
  });
  D.meshes.forEach((m) => {
    T.meshes.push(traces.length);
    traces.push({ type: "mesh3d", x: m.x, y: m.y, z: m.z, i: m.i, j: m.j, k: m.k, color: NODATA_MESH, opacity: 1,
      hovertemplate: esc(m.name) + "<extra></extra>", name: m.id,
      lighting: { ambient: 0.5, diffuse: 0.8, specular: 0.15, roughness: 0.6 }, lightposition: { x: 100, y: 200, z: 150 } });
  });
  T.markers = traces.length;
  const mk = { x: [], y: [], z: [], rid: [] };
  D.markers.forEach((m) => m.x.forEach((_, i) => { mk.x.push(m.x[i]); mk.y.push(m.y[i]); mk.z.push(m.z[i]); mk.rid.push(m.region); }));
  traces.push({ type: "scatter3d", mode: "markers", x: mk.x, y: mk.y, z: mk.z, customdata: mk.rid,
    marker: { size: 9, symbol: "diamond", color: mk.rid.map(() => NODATA_MESH), line: { color: "#1d1d1b", width: 1 } },
    text: mk.rid.map((r) => REG[r].name), hovertemplate: "%{text}<br><i>approximate position</i><extra></extra>", name: "nuclei" });

  const axis = { visible: false, showspikes: false };
  const layout = {
    margin: { l: 0, r: 0, t: 0, b: 0 }, paper_bgcolor: "#fcfcfb", showlegend: false, uirevision: "keep",
    scene: { xaxis: axis, yaxis: axis, zaxis: axis, aspectmode: "data", camera: { eye: { x: -1.15, y: 0.7, z: 0.38 } } },
    hoverlabel: { bgcolor: "#ffffff", bordercolor: "#c9c8c3", font: { color: "#1d1d1b", size: 12 } },
  };
  Plotly.newPlot("brain", traces, layout, { responsive: true, displaylogo: false, modeBarButtonsToRemove: ["toImage"] });

  // ---------- recolour ----------
  function recolor() {
    const spec = mapSpec();
    const cache = {};
    const val = (rids) => cache[rids.join()] || (cache[rids.join()] = regionsValue(rids, spec));
    D.cortex.forEach((c, n) => {
      const per = c.regions.map((r) => val([r]));
      const vc = c.vregion.map((ix) => (ix < 0 ? NODATA_CTX : per[ix].color || NODATA_CTX));
      c.border.forEach((v) => { if (c.vregion[v] >= 0) vc[v] = tint(vc[v], -0.3); });
      const tx = c.vregion.map((ix) => (ix < 0 ? "medial wall" : `${REG[c.regions[ix]].name}<br>${per[ix].label}`));
      Plotly.restyle("brain", { vertexcolor: [vc], text: [tx] }, [T.cortex[n]]);
    });
    const colors = [], hts = [];
    D.meshes.forEach((m) => { const v = val(m.regions); colors.push(v.color || NODATA_MESH); hts.push(`${esc(m.name)}<br>${esc(v.label)}<extra></extra>`); });
    Plotly.restyle("brain", { color: colors, hovertemplate: hts }, T.meshes);
    const mv = mk.rid.map((r) => val([r]));
    Plotly.restyle("brain", { "marker.color": [mv.map((v) => v.color || NODATA_MESH)],
      text: [mk.rid.map((r, i) => `${REG[r].name}<br>${mv[i].label}`)] }, [T.markers]);
    $("outsideChips").innerHTML = D.outside.map((r) => {
      const v = val([r]);
      return `<span class="chip" data-r="${esc(r)}" style="border-color:${v.color || "#e4e3df"}">${esc(REG[r].name)} · ${esc(v.label)}</span>`;
    }).join("");
    renderLegend(spec);
  }

  function renderLegend(spec) {
    let h = `<div><b>${esc(spec.label)}</b></div>`;
    if (spec.kind === "molecule") {
      h += ["increased", "decreased", "mixed", "no change"].map((d) =>
        `<div class="row"><span class="sw" style="background:${DIR_COLOR[d]}"></span>${d === "mixed" ? "mixed (increased and decreased reported)" : d === "no change" ? "no change / present / altered" : d}</div>`).join("") +
        `<div class="row"><span class="sw" style="background:${tint(DIR_COLOR.increased, 0.6)}"></span><span class="sw" style="background:${tint(DIR_COLOR.decreased, 0.6)}"></span>light: reported only for a broader area (lobe, “cortex”)</div>`;
    } else if (spec.kind === "count") {
      h += `<div class="bar" style="background:linear-gradient(90deg,${SEQ.join(",")})"></div><div class="ticks"><span>1</span><span>${Math.round(spec.max)}</span></div>`;
    } else if (spec.kind === "ordinal") {
      h += `<div class="row">${SEQ.slice(0, 6).map((c, i) => `<span class="sw" style="background:${c}"></span>R${i + 1}`).join(" ")}</div><div class="muted">R1 medulla → R6 neocortex</div>`;
    } else {
      h += `<div class="bar" style="background:linear-gradient(90deg,${DIVERGING.join(",")})"></div><div class="ticks"><span>−${spec.max.toFixed(2)}</span><span>0</span><span>+${spec.max.toFixed(2)}</span></div><div class="muted">d &lt; 0: thinner / smaller in PD</div>`;
    }
    h += `<div class="row"><span class="sw" style="background:${NODATA_MESH}"></span><span class="muted">no data</span></div>`;
    $("legend").innerHTML = h;
  }

  // ---------- panel ----------
  const pubmed = (p) => `<a href="https://pubmed.ncbi.nlm.nih.gov/${encodeURIComponent(p)}/" target="_blank" rel="noopener noreferrer">PMID ${esc(p)}</a>`;
  const dirChip = (d) => `<span class="dir ${esc(String(d).replace(/\s+/g, "-"))}">${esc(d)}</span>`;
  const nStr = (f) => (f.n_pd != null || f.n_control != null ? `${f.n_pd ?? "?"}/${f.n_control ?? "?"}` : "");

  function findingHtml(f, showRegion) {
    return `<div class="finding"><div class="head"><span class="mol" data-mol="${esc(f.molecule)}">${esc(f.molecule)}</span>
      ${dirChip(f.direction)}<span class="muted">${esc(f.molecule_type || "")}</span>${showRegion ? `<span class="muted">· ${esc(f.region_text || f.region)}</span>` : ""}</div>
      <div>${esc(f.finding || "")}${f.comparison ? ` <span class="muted">(${esc(f.comparison)})</span>` : ""}</div>
      <div class="quote">“${esc(f.quote)}”</div>
      <div class="ref">${pubmed(f.pmid)} · ${esc(f.journal || "")} ${esc(f.year || "")} · ${esc(f.method || "")}${nStr(f) ? ` · n ${esc(nStr(f))}` : ""}
      ${f.approximate_region ? " · <i>region mapped approximately</i>" : ""}${f.notes ? ` · ${esc(f.notes)}` : ""}</div></div>`;
  }

  const pmLink = (q) => `https://pdmap.uni.lu/minerva/index.xhtml?id=${encodeURIComponent(PM.project || "")}&search=${encodeURIComponent(q)}`;
  function pdmapHtml(mol) {
    const m = pmGet(mol);
    if (!m) return PM.project ? `<h3 id="detail-mechanisms">PD map mechanisms</h3><p class="muted">Not curated in the PD map (${esc(PM.project)}).</p>` : "";
    return `<h3 id="detail-mechanisms">PD map mechanisms</h3>
      <div class="finding"><div><b>Pathways:</b> ${m.pw.map((p) => `<span class="chip" data-pw="${esc(p)}">${esc(p)}</span>`).join(" ")}</div>
      <div><b>Compartments:</b> ${esc(m.c.join(", ") || "—")}</div>
      <div><b>Curated reactions:</b> ${m.nr} · supported by ${m.np} papers: ${m.pm.map(pubmed).join(", ")}</div>
      <div><b>Interaction partners:</b> ${m.pa.map((x) => `<span class="mol" data-mol="${esc(x)}">${esc(x)}</span>`).join(", ") || "—"}</div>
      <div class="ref"><a href="${pmLink(mol)}" target="_blank" rel="noopener noreferrer">Open in the PD map</a> · ${esc(PM.licence || "")}</div></div>`;
  }
  function geneticsHtml(mol) {
    const g = gnGet(mol);
    if (!GN.genes) return "";
    if (!g) return `<h3 id="detail-genetics">PD genetics</h3><p class="muted">No PD GWAS locus within 1 Mb, no GWAS Catalog PD association, not on the PD gene list.</p>`;
    const tags = [g.m ? `<span class="dir mixed">known PD gene</span>` : "", g.n ? `<span class="chip">nominated by Nalls 2019</span>` : ""].join(" ");
    const loci = g.l.map(([gw, name, snp, sc]) => `<div>${esc(gw)} · locus <b>${esc(name)}</b> (${esc(snp)}) · browser evidence score ${esc(sc)}</div>`).join("");
    const cat = g.c ? `<div><b>GWAS Catalog:</b> ${g.c[0]} genome-wide significant PD associations in ${g.c[3]} studies; strongest ${esc(g.c[2] || "")} p = ${esc(Number(g.c[1]).toPrecision(2))}${g.c[5].length ? ` · ${esc(g.c[5].join(", "))}` : ""} · ${g.c[4].map(pubmed).join(", ")}</div>` : "";
    return `<h3 id="detail-genetics">PD genetics</h3><div class="finding">${tags}${loci}${cat}<div class="ref">Locus proximity is not causality. ${esc(Object.values(GN.sources || {}).join(" · "))}</div></div>`;
  }
  function parkdbHtml(mol) {
    const g = pkGet(mol);
    if (!PK.genes) return "";
    if (!g) return `<h3 id="detail-parkdb">ParkDB · prefrontal cortex RNA-seq</h3><p class="muted">Not differentially expressed in PD (padj ≥ 0.05 in both datasets) or not measured.</p>`;
    const order = ["PD", "DLB", "MSA", "PSP", "CBD", "AD"];
    const rows = order.flatMap((d) => ["D1", "D2"].map((s) => [`${s}:${d}`, g.c[`${s}:${d}`]])).filter(([, v]) => v);
    return `<h3 id="detail-parkdb">ParkDB · prefrontal cortex RNA-seq</h3>
      <p class="muted">${g.rep ? "<b>Replicated PD change</b> (padj &lt; 0.05 in both datasets, same direction)." : "PD change not replicated across the two datasets."}
      log2FC &gt; 0 = higher in disease; bulk tissue.</p><table><thead><tr><th>Contrast</th><th>n</th><th>log2FC</th><th>padj</th></tr></thead><tbody>` +
      rows.map(([k, [lfc, q]]) => `<tr><td>${esc(k.replace("D1:", "Dataset 1 · ").replace("D2:", "Dataset 2 · "))} vs control</td><td class="muted">${esc((PK.groups[k] || []).join("/"))}</td>
        <td>${q < 0.05 ? dirChip(lfc > 0 ? "increased" : "decreased") : ""} ${lfc.toFixed(2)}</td><td>${esc(q)}</td></tr>`).join("") +
      `</tbody></table><div class="ref">${esc(PK.source || "")}</div>`;
  }
  function regionPathways(lit) {
    const c = {};
    const mols = new Set(lit.map((f) => molKey(f.molecule)));
    mols.forEach((k) => { const m = PMM[k]; if (m) m.pw.forEach((p) => (c[p] = c[p] || new Set()).add(k)); });
    const top = Object.entries(c).sort((a, b) => b[1].size - a[1].size).slice(0, 10);
    if (!top.length) return "";
    return `<h3>PD map pathways of the molecules reported here</h3>` + top.map(([p, ks]) =>
      `<div class="finding"><span class="chip" data-pw="${esc(p)}">${esc(p)}</span> <span class="muted">${ks.size} molecules:</span> ${[...ks].slice(0, 12).map((k) => `<span class="mol" data-mol="${esc(k)}">${esc(k)}</span>`).join(", ")}</div>`).join("");
  }

  function showRegions(rids, title) {
    state.selection = rids;
    const prim = REG[rids[0]];
    const ev = rids.flatMap((r) => evByRegion[r] || []);
    const seen = new Set();
    const lit = rids.flatMap((r) => litByRegion[r] || []).filter((f) => !seen.has(f._i) && seen.add(f._i));
    lit.sort((a, b) => molKey(a.molecule).localeCompare(molKey(b.molecule)) || (b.year || 0) - (a.year || 0));
    const braak = rids.map((r) => REG[r].braak).find((b) => b);
    const img = ev.filter((e) => e.layer === "imaging" && /\(PDvsCN\)|free water/.test(e.feature));
    const tissue = ev.filter((e) => ["observed_pd_tissue"].includes(e.grade));
    const mech = ev.filter((e) => e.grade === "inferred_mechanism");
    const bg = ev.filter((e) => e.grade === "healthy_background" && typeof e.value === "number")
      .sort((a, b) => Math.abs(b.value) - Math.abs(a.value)).slice(0, 6);
    const geo = { surface: "cortical parcel (Desikan-Killiany)", mesh: "MASSP volumetric label", marker: "approximate marker (no volumetric label)", outside: "outside the brain" }[prim.geometry];
    let h = `<h2>${esc(title || prim.name)}</h2><dl class="kv">
      <dt>Atlas regions</dt><dd>${rids.map((r) => esc(REG[r].name)).join("; ")}</dd>
      <dt>Braak (Keo 2020)</dt><dd>${braak ? esc(braak) : "not assigned"}</dd>
      <dt>Geometry</dt><dd>${esc(geo)}</dd>
      <dt>Published findings</dt><dd>${lit.length} from ${new Set(lit.map((f) => f.pmid)).size} papers</dd></dl>`;
    if (img.length) h += `<h3>Imaging, PD vs controls</h3>` + img.map((e) =>
      `<div class="finding"><b>${esc(e.feature)}</b>: ${esc(typeof e.value === "number" ? e.value.toFixed(3) : e.value)} ${e.p != null ? `(p = ${esc(Number(e.p).toPrecision(2))})` : ""} <span class="ref">${esc(e.source)} · n ${esc(e.n)}</span></div>`).join("");
    if (lit.length) h += `<h3>Published PD tissue findings</h3>` + lit.map((f) => findingHtml(f, true)).join("");
    h += regionPathways(lit);
    const risk = [...new Set(lit.map((f) => molKey(f.molecule)))].filter((k) => isRisk(gnGet(k)));
    if (risk.length) h += `<h3>PD genetic risk genes reported here</h3><div>${risk.map((k) => `<span class="mol" data-mol="${esc(k)}">${esc(k)}</span>${gnGet(k).m ? "*" : ""}`).join(", ")}</div><p class="muted">* known PD gene; others: GWAS locus genes (GWAS Locus Browser / GWAS Catalog).</p>`;
    const cm = ev.filter((e) => e.grade === "curated_mechanism");
    if (cm.length) h += `<h3>PD map cell-type links</h3>` + cm.map((e) =>
      `<div class="finding"><b>${esc(e.feature)}</b> <span class="muted">(${esc(e.value)} molecules; ${esc(e.method)})</span><div class="ref">${esc(e.note)}</div></div>`).join("");
    if (tissue.length) h += `<h3>Biomni-PD analyses of public data</h3>` + tissue.map((e) =>
      `<div class="finding"><b>${esc(e.feature)}</b> ${e.direction ? dirChip(e.direction) : ""}<div>${esc(e.value)} ${esc(e.unit || "")}${e.q != null ? ` (q = ${esc(Number(e.q).toPrecision(2))})` : ""}</div><div class="ref">${esc(e.source)} · ${esc(e.method)}${e.note ? " · " + esc(e.note) : ""}</div></div>`).join("");
    if (mech.length) h += `<h3>Mechanism hypotheses</h3>` + mech.map((e) =>
      `<div class="finding"><b>${esc(e.feature)}</b>: ${esc(e.value)}<div class="ref">${esc(e.note)}</div></div>`).join("");
    if (bg.length) h += `<h3>Healthy-brain expression context (AHBA)</h3><p class="muted">This shows predisposition, not change in PD.</p>` + bg.map((e) =>
      `<div>${esc(e.feature)}: ${esc(e.value.toFixed(2))}</div>`).join("");
    if (!lit.length && !tissue.length && !img.length) h += `<p class="muted">No PD evidence mapped to this region yet.</p>`;
    $("panel").innerHTML = h;
    $("panel").scrollTop = 0;
  }

  function showMolecule(mol) {
    const fs = litByMol[molKey(mol)] || [];
    const byRegion = {};
    fs.forEach((f) => { const k = f.region; (byRegion[k] = byRegion[k] || []).push(f.direction); });
    let h = `<h2>${esc(mol)}</h2><p class="muted">${fs.length} published findings · ${new Set(fs.map((f) => f.pmid)).size} papers · ${Object.keys(byRegion).length} regions</p>`;
    const extra = D.evidence.filter((e) => e.grade === "observed_pd_tissue" && molKey(e.feature).includes(molKey(mol)));
    if (Object.keys(byRegion).length) h += `<h3>By region</h3><table><tbody>` + Object.entries(byRegion).sort().map(([r, ds]) =>
      `<tr><td>${esc(r)}</td><td>${[...new Set(ds)].map(dirChip).join(" ")}</td><td class="muted">${ds.length}</td></tr>`).join("") + `</tbody></table>`;
    if (extra.length) h += `<h3>Biomni-PD analyses</h3>` + extra.map((e) => `<div class="finding"><b>${esc(REG[e.region_id].name)}</b>: ${esc(e.feature)} ${e.direction ? dirChip(e.direction) : ""}<div class="ref">${esc(e.source)}</div></div>`).join("");
    h += geneticsHtml(mol);
    h += parkdbHtml(mol);
    h += pdmapHtml(mol);
    if (fs.length) h += `<h3>Published findings</h3>` + fs.map((f) => findingHtml(f, true)).join("");
    if (!fs.length && !extra.length && !pmGet(mol) && !pkGet(mol) && !gnGet(mol)) h += `<p>No region-level PD tissue findings for “${esc(mol)}” in this atlas yet.</p>`;
    $("panel").innerHTML = h;
    $("panel").scrollTop = 0;
  }

  // ---------- interactions ----------
  document.getElementById("brain").on("plotly_click", (e) => {
    const p = e.points && e.points[0];
    if (!p) return;
    const ti = p.curveNumber;
    const ci = T.cortex.indexOf(ti);
    if (ci >= 0) {
      const c = D.cortex[ci];
      let v = p.pointNumber;
      if (v == null || v < 0) { // nearest vertex fallback
        let best = Infinity;
        c.x.forEach((x, i) => { const d = (x - p.x) ** 2 + (c.y[i] - p.y) ** 2 + (c.z[i] - p.z) ** 2; if (d < best) { best = d; v = i; } });
      }
      const ix = c.vregion[v];
      if (ix >= 0) showRegions([c.regions[ix]]);
      return;
    }
    const mi = T.meshes.indexOf(ti);
    if (mi >= 0) return showRegions(D.meshes[mi].regions, D.meshes[mi].name);
    if (ti === T.markers) return showRegions([mk.rid[p.pointNumber]]);
  });
  document.body.addEventListener("click", (e) => {
    const shortcut = e.target.closest("[data-evidence-section]");
    if (shortcut) {
      e.preventDefault();
      const hash = shortcut.getAttribute("href");
      if (location.hash === hash) fromHash();
      else location.hash = hash;
      return;
    }
    const pw = e.target.closest("[data-pw]");
    if (pw) { setPathway(pw.dataset.pw); return; }
    const m = e.target.closest("[data-mol]");
    if (m) { setMolecule(m.dataset.mol); return; }
    const r = e.target.closest("[data-r]");
    if (r) showRegions([r.dataset.r]);
  });

  const mapSel = $("map");
  function fillMapSelect() {
    mapSel.innerHTML = (state.molecule ? `<option value="molecule">${esc(state.molecule)}: direction by region</option>` : "") +
      (state.pathway ? `<option value="pathway">PD map: ${esc(state.pathway)}</option>` : "") +
      (GN.genes ? `<option value="genetics">Genetics: PD risk genes with findings</option>` : "") +
      Object.entries(D.maps).map(([k, m]) => `<option value="${esc(k)}">${esc(m.label)}</option>`).join("");
    mapSel.value = state.map;
  }
  mapSel.addEventListener("change", () => {
    state.map = mapSel.value; recolor();
    if (state.map.startsWith("imaging") && $("view").value === "deep") { $("view").value = "lh_lat"; setView("lh_lat"); }
  });
  $("opacity").addEventListener("input", (e) => Plotly.restyle("brain", { opacity: e.target.value / 100 }, T.cortex));
  // views: deep = glass cortex over subcortex; cortex views = opaque surface, one hemisphere for medial views
  const VIEWS = {
    lh_lat: { hemis: ["L", "R"], eye: { x: -1.6, y: 0, z: 0.1 } }, rh_lat: { hemis: ["L", "R"], eye: { x: 1.6, y: 0, z: 0.1 } },
    lh_med: { hemis: ["L"], eye: { x: 1.6, y: 0, z: 0.1 } }, rh_med: { hemis: ["R"], eye: { x: -1.6, y: 0, z: 0.1 } },
    top: { hemis: ["L", "R"], eye: { x: 0, y: -0.01, z: 1.7 } }, bottom: { hemis: ["L", "R"], eye: { x: 0, y: 0.01, z: -1.7 } },
  };
  function setView(v) {
    const deep = v === "deep";
    const cfg = VIEWS[v] || { hemis: ["L", "R"], eye: { x: -1.15, y: 0.7, z: 0.38 } };
    const op = deep ? 0.18 : 1;
    $("opacity").value = op * 100;
    const geo = (c, k) => (deep ? c[k] : c.inflated["xyz".indexOf(k)]);
    Plotly.restyle("brain", { opacity: op, visible: D.cortex.map((c) => cfg.hemis.includes(c.hemi)),
      x: D.cortex.map((c) => geo(c, "x")), y: D.cortex.map((c) => geo(c, "y")), z: D.cortex.map((c) => geo(c, "z")) }, T.cortex);
    Plotly.restyle("brain", { visible: deep }, T.meshes.concat([T.context, T.markers]));
    $("showMarkers").checked = deep; $("showOutline").checked = deep;
    Plotly.relayout("brain", { "scene.camera": { eye: cfg.eye, up: { x: 0, y: 0, z: 1 }, center: { x: 0, y: 0, z: 0 } } });
  }
  $("view").addEventListener("change", (e) => setView(e.target.value));
  $("showMarkers").addEventListener("change", (e) => Plotly.restyle("brain", { visible: e.target.checked }, [T.markers]));
  $("showOutline").addEventListener("change", (e) => Plotly.restyle("brain", { visible: e.target.checked }, [T.context]));

  function setPathway(pw) {
    state.pathway = pw; state.map = "pathway"; $("pathway").value = pw;
    fillMapSelect(); recolor();
    const ms = Object.entries(PMM).filter(([, v]) => v.pw.includes(pw)).map(([k]) => k).sort();
    const info = (PM.pathways || {})[pw] || {};
    const reported = ms.filter((k) => litByMol[k]);
    $("panel").innerHTML = `<h2>${esc(pw)}</h2><p class="muted">PD map pathway · ${ms.length} molecules · ${info.n_pmids || 0} curated papers</p>
      <p><a href="${pmLink(pw)}" target="_blank" rel="noopener noreferrer">Open in the PD map</a></p>
      <h3>Molecules with published PD tissue findings (${reported.length})</h3><div>${reported.map((k) => `<span class="mol" data-mol="${esc(k)}">${esc(k)}</span> <span class="muted">(${litByMol[k].length})</span>`).join(", ") || "none"}</div>
      <h3>All pathway molecules</h3><div class="muted">${ms.map((k) => `<span class="mol" data-mol="${esc(k)}">${esc(k)}</span>`).join(", ")}</div>`;
    $("panel").scrollTop = 0;
  }
  function setMolecule(mol) {
    const key = molKey(mol);
    const canonical = (litByMol[key] && litByMol[key][0].molecule) || mol;
    state.molecule = canonical; state.map = "molecule";
    $("q").value = canonical; $("fText").value = canonical;
    fillMapSelect(); recolor(); showMolecule(canonical); renderTable();
  }
  const regionByName = {};
  D.regions.forEach((r) => (regionByName[r.name.toLowerCase()] = r.id));
  const meshByName = {};
  D.meshes.forEach((m) => (meshByName[m.name.toLowerCase()] = m));
  function onSearch() {
    const q = $("q").value.trim();
    if (!q) return;
    const ql = q.toLowerCase();
    if (meshByName[ql]) return showRegions(meshByName[ql].regions, meshByName[ql].name);
    if (regionByName[ql]) return showRegions([regionByName[ql]]);
    if (litByMol[molKey(q)] || pmGet(q) || pkGet(q) || gnGet(q)) return setMolecule(q);
    const pwHit = PATHWAYS.find((p) => p.toLowerCase() === ql);
    if (pwHit) return setPathway(pwHit);
    const part = D.regions.filter((r) => r.name.toLowerCase().includes(ql));
    if (part.length && part.length <= 4) return showRegions(part.map((r) => r.id));
    const molPart = Object.keys(litByMol).filter((k) => k.includes(molKey(q)));
    if (molPart.length === 1) return setMolecule(litByMol[molPart[0]][0].molecule);
    setMolecule(q);
  }
  $("q").addEventListener("change", onSearch);
  $("q").addEventListener("keydown", (e) => { if (e.key === "Enter") onSearch(); });
  $("clear").addEventListener("click", () => {
    state.molecule = null; state.map = "lit_count"; $("q").value = ""; $("fText").value = "";
    fillMapSelect(); recolor(); renderTable();
  });
  $("pathway").innerHTML = `<option value="">PD map pathway…</option>` + PATHWAYS.map((p) => `<option>${esc(p)}</option>`).join("");
  $("pathway").addEventListener("change", (e) => { if (e.target.value) setPathway(e.target.value); });
  const mols = Object.entries(litByMol).sort((a, b) => b[1].length - a[1].length).map(([, fs]) => fs[0].molecule);
  $("suggest").innerHTML = mols.map((m) => `<option value="${esc(m)}">`).join("") +
    D.meshes.map((m) => `<option value="${esc(m.name)}">`).join("") +
    Object.keys(PMM).filter((k) => !litByMol[k]).map((k) => `<option value="${esc(k)}">`).join("") +
    D.regions.filter((r) => r.geometry !== "mesh").map((r) => `<option value="${esc(r.name)}">`).join("");

  // ---------- table ----------
  const analysisRows = D.evidence.filter((e) => ["observed_pd_tissue", "inferred_mechanism"].includes(e.grade) || /free water/.test(e.feature));
  const tableRows = LIT.map((f) => ({ src: "literature", region: f.region_text || f.region, rid: f.regions[0], mol: f.molecule, type: f.molecule_type || "",
    dir: f.direction, finding: f.finding || "", method: f.method || "", n: nStr(f), source: f.pmid, year: f.year, quote: f.quote }))
    .concat(analysisRows.map((e) => ({ src: "analysis", region: REG[e.region_id].name, rid: e.region_id, mol: e.feature, type: e.layer,
      dir: e.direction || "", finding: String(e.value ?? ""), method: e.method || "", n: e.n || "", source: e.source, year: "", quote: e.note || "" })));
  const types = [...new Set(LIT.map((f) => f.molecule_type).filter(Boolean))].sort();
  $("fType").innerHTML += types.map((t) => `<option>${esc(t)}</option>`).join("");
  let filtered = tableRows;
  function renderTable() {
    const s = $("fSource").value, d = $("fDir").value, t = $("fType").value, q = $("fText").value.trim().toLowerCase();
    filtered = tableRows.filter((r) => (!s || r.src === s) && (!d || r.dir === d) && (!t || r.type === t) &&
      (!q || [r.region, r.mol, r.finding, r.method, r.quote, r.source].join(" ").toLowerCase().includes(q)));
    $("tbl").querySelector("tbody").innerHTML = filtered.slice(0, 1500).map((r) => `<tr>
      <td><span class="chip" data-r="${esc(r.rid)}">${esc(r.region)}</span></td>
      <td>${r.src === "literature" ? `<span class="mol" data-mol="${esc(r.mol)}">${esc(r.mol)}</span>` : esc(r.mol)}</td>
      <td>${esc(r.type)}</td><td>${r.dir ? dirChip(r.dir) : ""}</td><td>${esc(r.finding)}</td><td>${esc(r.method)}</td><td>${esc(r.n)}</td>
      <td>${r.src === "literature" ? pubmed(r.source) + (r.year ? ` (${esc(r.year)})` : "") : esc(r.source)}</td></tr>`).join("");
    $("tblCount").textContent = `${filtered.length} of ${tableRows.length} rows` + (filtered.length > 1500 ? " (first 1,500 shown)" : "");
  }
  ["fSource", "fDir", "fType"].forEach((id) => $(id).addEventListener("change", renderTable));
  $("fText").addEventListener("input", renderTable);
  $("dl").addEventListener("click", () => {
    const cols = ["src", "region", "mol", "type", "dir", "finding", "method", "n", "source", "year", "quote"];
    const clean = (v) => String(v ?? "").replace(/[\t\n\r]+/g, " ");
    const tsv = [cols.join("\t")].concat(filtered.map((r) => cols.map((c) => clean(r[c])).join("\t"))).join("\n");
    const a = document.createElement("a");
    a.href = URL.createObjectURL(new Blob([tsv], { type: "text/tab-separated-values" }));
    a.download = "pd_brain_atlas_findings.tsv"; a.click(); URL.revokeObjectURL(a.href);
  });

  $("stats").textContent = `${LIT.length} published findings · ${new Set(LIT.map((f) => f.pmid)).size} papers · ${Object.keys(litByMol).length} molecules · ${D.regions.length} regions`;
  $("sources").textContent = "Sources: " + Object.entries(D.sources).map(([k, v]) => `${k}: ${v}`).join(" · ");
  fillMapSelect(); recolor(); renderTable();
  // deep links: #mol=SNCA or #region=locus%20coeruleus (also updated on search)
  function fromHash() {
    const h = new URLSearchParams(location.hash.slice(1));
    if (h.get("mol")) setMolecule(h.get("mol"));
    if (h.get("pathway")) setPathway(h.get("pathway"));
    if (h.get("region")) { $("q").value = h.get("region"); onSearch(); }
    if (h.get("map") && (D.maps[h.get("map")] || ["molecule", "genetics"].includes(h.get("map")))) { state.map = h.get("map"); fillMapSelect(); recolor(); }
    if (h.get("view")) { $("view").value = h.get("view"); setView(h.get("view")); }
    if (["genetics", "parkdb", "mechanisms"].includes(h.get("section"))) {
      const target = document.getElementById(`detail-${h.get("section")}`);
      if (target) {
        const panel = $("panel");
        panel.scrollTop += target.getBoundingClientRect().top - panel.getBoundingClientRect().top - 16;
      }
    }
  }
  window.addEventListener("hashchange", fromHash);
  fromHash();
})();
