"use strict";

const $ = (s) => document.querySelector(s);
const POPS = ["AFR", "AMR", "EAS", "EUR", "SAS"];
const POP_NAMES = { AFR: "African", AMR: "Admixed American", EAS: "East Asian", EUR: "European", SAS: "South Asian" };
const PAGE = 100;

let maxBytes = 1024 ** 3;
let report = null;
let shown = PAGE;

const esc = (s) => String(s ?? "").replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
const fmtBytes = (n) => n >= 1024 ** 3 ? (n / 1024 ** 3).toFixed(2) + " GB" : n >= 1024 ** 2 ? (n / 1024 ** 2).toFixed(1) + " MB" : (n / 1024).toFixed(0) + " KB";
const pct = (x, d = 1) => (x * 100).toFixed(d) + "%";
const popVar = (p) => `var(--pop-${p})`;

function freqText(f) {
  if (f == null) return "—";
  if (f === 0) return "0";
  if (f < 0.001) return "<0.1%";
  return pct(f, f < 0.1 ? 1 : 0);
}

// ---------------------------------------------------------------- upload

async function loadStatus() {
  const res = await fetch("/api/status");
  const s = await res.json();
  maxBytes = s.max_upload_bytes;
  $("#max-size").textContent = fmtBytes(maxBytes);
  const dbs = s.databases;
  const badge = (name, meta) => meta
    ? `<span class="ok">${name} · ${Number(meta.records || 0).toLocaleString()} records · built ${esc((meta.built_at || "").slice(0, 10))}</span>`
    : `<span class="missing">${name} not installed</span>`;
  $("#db-status").innerHTML = badge("ClinVar", dbs.clinvar) + badge("1000 Genomes", dbs.kg);
}

function showError(msg) {
  const e = $("#error");
  e.textContent = msg;
  e.hidden = !msg;
}

function setProgress(fraction, text) {
  $("#progress").hidden = false;
  $("#bar-fill").style.width = Math.round(fraction * 100) + "%";
  $("#progress-text").textContent = text;
}

function upload(file) {
  showError("");
  if (file.size > maxBytes) {
    showError(`That file is ${fmtBytes(file.size)}; the limit is ${fmtBytes(maxBytes)}.`);
    return;
  }
  const params = new URLSearchParams({ filename: file.name, assembly: $("#assembly").value });
  const xhr = new XMLHttpRequest();
  xhr.open("POST", "/api/upload?" + params);
  xhr.setRequestHeader("Content-Type", "application/octet-stream");
  xhr.upload.onprogress = (e) => {
    if (e.lengthComputable) setProgress(0.4 * e.loaded / e.total, `Uploading ${fmtBytes(e.loaded)} of ${fmtBytes(e.total)}`);
  };
  xhr.onload = () => {
    let body = {};
    try { body = JSON.parse(xhr.responseText); } catch (_) { /* non-JSON error */ }
    if (xhr.status !== 202) {
      $("#progress").hidden = true;
      showError(body.detail || `Upload failed (HTTP ${xhr.status})`);
      return;
    }
    poll(body.job_id);
  };
  xhr.onerror = () => { $("#progress").hidden = true; showError("Upload failed — check your connection."); };
  setProgress(0, "Starting upload…");
  xhr.send(file);
}

async function poll(jobId) {
  for (;;) {
    const res = await fetch(`/api/jobs/${jobId}`);
    if (!res.ok) { showError("Lost track of the analysis job."); return; }
    const job = await res.json();
    if (job.state === "error") { $("#progress").hidden = true; showError(job.error); return; }
    if (job.state === "done") { await openReport(jobId); return; }
    setProgress(0.4 + 0.6 * (job.progress || 0), job.message || "Analysing…");
    await new Promise((r) => setTimeout(r, 1000));
  }
}

async function openReport(jobId) {
  const res = await fetch(`/api/jobs/${jobId}/result`);
  if (!res.ok) { history.replaceState(null, "", location.pathname); return; }
  report = await res.json();
  history.replaceState(null, "", "#report=" + jobId);
  $("#upload-view").hidden = true;
  $("#results-view").hidden = false;
  $("#progress").hidden = true;
  render();
}

// ---------------------------------------------------------------- render

function render() {
  renderFile();
  renderClinvarSummary();
  renderAncestry();
  const sel = $("#cat-filter");
  sel.innerHTML = `<option value="">All listed categories</option>` + report.categories
    .filter((c) => report.findings.some((f) => f.category === c.key))
    .map((c) => `<option value="${c.key}">${esc(c.label)}</option>`).join("");
  shown = PAGE;
  renderFindings();
  const notes = report.notes || [];
  $("#notes").hidden = notes.length === 0;
  $("#notes").innerHTML = `<h2>Notes</h2><ul>${notes.map((n) => `<li>${esc(n)}</li>`).join("")}</ul>`;
}

function renderFile() {
  const f = report.file;
  const fmtName = { vcf: "VCF", "23andme": "23andMe-style", ancestrydna: "AncestryDNA", myheritage: "MyHeritage / FTDNA" }[f.format] || f.format;
  const rows = [
    ["File", esc(f.name)],
    ["Size", fmtBytes(f.size)],
    ["Format", esc(fmtName)],
    ["Genome build", esc(f.assembly)],
    f.sample ? ["Sample", esc(f.sample)] : null,
    ["Genotype calls", f.genotypes.toLocaleString()],
    f.call_rate != null ? ["Call rate", pct(f.call_rate)] : null,
    ["No-calls", f.no_calls.toLocaleString()],
  ].filter(Boolean);
  $("#file-summary").innerHTML = `<h2>Your file</h2><dl class="kv">${rows.map(([k, v]) => `<dt>${k}</dt><dd>${v}</dd>`).join("")}</dl>`;
}

function renderClinvarSummary() {
  const cv = report.clinvar;
  const el = $("#clinvar-summary");
  if (!cv.available) {
    el.innerHTML = `<h2>ClinVar</h2><p class="muted">ClinVar database is not installed on this server.</p>`;
    return;
  }
  const checked = Object.values(cv.checked).reduce((a, b) => a + b, 0);
  const items = report.categories.map((c) => {
    const n = cv.carried[c.key] || 0;
    const listed = report.findings.some((f) => f.category === c.key);
    const label = listed ? `<button data-cat="${c.key}">${esc(c.label)}</button>` : esc(c.label);
    return `<li><span><span class="badge c-${c.key}">${n.toLocaleString()}</span> ${label}</span>
            <span class="muted num">${(cv.checked[c.key] || 0).toLocaleString()} ${cv.checked[c.key] === 1 ? "site" : "sites"} checked</span></li>`;
  }).join("");
  el.innerHTML = `<h2>ClinVar overview</h2>
    <p class="muted">${checked.toLocaleString()} of your genotypes are at ClinVar-classified sites.
      Counts show how many you carry the classified allele for.</p>
    <ul class="cat-list">${items}</ul>`;
  el.querySelectorAll("button[data-cat]").forEach((b) => b.addEventListener("click", () => {
    $("#cat-filter").value = b.dataset.cat;
    shown = PAGE;
    renderFindings();
    $("#findings").scrollIntoView({ behavior: "smooth" });
  }));
}

function renderAncestry() {
  const a = report.ancestry;
  const el = $("#ancestry");
  if (!a.available) {
    el.innerHTML = `<h2>Ancestry (1000 Genomes)</h2><p class="muted">${esc(a.reason)}</p>`;
    return;
  }
  const bar = a.populations.map((p) =>
    `<div style="width:${p.proportion * 100}%;background:${popVar(p.code)}" title="${esc(p.name)} ${pct(p.proportion)}"></div>`).join("");
  const rows = a.populations.map((p) => `
    <div class="pop-row">
      <span class="swatch" style="background:${popVar(p.code)}"></span>
      <span>${esc(p.name)} <span class="muted">(${p.code})</span></span>
      <div class="track">
        <div class="fill" style="width:${p.proportion * 100}%;background:${popVar(p.code)}"></div>
        <div class="ci" style="left:${p.ci_low * 100}%;width:${Math.max(0, p.ci_high - p.ci_low) * 100}%"></div>
      </div>
      <span class="pct">${pct(p.proportion)} <span class="muted small">±${pct((p.ci_high - p.ci_low) / 2, 1)}</span></span>
    </div>`).join("");
  el.innerHTML = `<h2>Ancestry (1000 Genomes superpopulations)</h2>
    <p class="muted">Estimated from ${a.sites_used.toLocaleString()} ancestry-informative SNPs by fitting your genotypes
      as a mixture of the five 1000 Genomes Phase 3 superpopulations. Whiskers show a 95% bootstrap interval.
      This is a continental-level estimate, not a genealogy.</p>
    <div class="ancestry-bar">${bar}</div>
    <div class="pop-rows">${rows}</div>
    ${a.note ? `<p class="muted">${esc(a.note)}</p>` : ""}`;
}

function filtered() {
  const q = $("#search").value.trim().toLowerCase();
  const cat = $("#cat-filter").value;
  const stars = Number($("#star-filter").value);
  const zyg = $("#zyg-filter").value;
  return report.findings.filter((f) =>
    (!cat || f.category === cat) &&
    f.stars >= stars &&
    (!zyg || f.zygosity === zyg) &&
    (!q || [f.gene, f.rsid, f.name, f.significance, ...f.conditions].join(" ").toLowerCase().includes(q)));
}

function freqCell(f) {
  if (!f.frequencies) return `<span class="muted">Not in 1000 Genomes</span>`;
  const fr = f.frequencies;
  const max = Math.max(...POPS.map((p) => fr[p] || 0), 0.0001);
  const bars = POPS.map((p) => {
    const h = Math.max(4, ((fr[p] || 0) / max) * 100);
    return `<div class="f" title="${POP_NAMES[p]}: ${freqText(fr[p])}"><div style="height:${h}%;background:${popVar(p)}"></div></div>`;
  }).join("");
  return `<div class="freqs">${bars}</div><div class="freq-label">Global ${freqText(fr.ALL)} · max ${POPS.reduce((a, p) => (fr[p] || 0) > (fr[a] || 0) ? p : a)} ${freqText(max)}</div>`;
}

function renderFindings() {
  const rows = filtered();
  const label = Object.fromEntries(report.categories.map((c) => [c.key, c.label]));
  $("#findings-count").textContent = `${rows.length.toLocaleString()} of ${report.findings.length.toLocaleString()} findings`;
  $("#findings tbody").innerHTML = rows.slice(0, shown).map((f) => `
    <tr>
      <td class="gene">${esc(f.gene || "—")}</td>
      <td>
        <div>${f.rsid ? `<code>${esc(f.rsid)}</code>` : ""} <span class="small">chr${esc(f.chrom)}:${f.pos} ${esc(f.ref)}→${esc(f.alt)}</span></div>
        <div class="small">${esc(f.name)}</div>
      </td>
      <td><code>${esc(f.genotype)}</code><div class="small">${f.zygosity === "homozygous" ? "2 copies" : "1 copy"} of ${esc(f.alt)}${f.opposite_strand ? " (reported on the opposite strand)" : ""}</div></td>
      <td>
        <span class="badge c-${f.category}">${esc(label[f.category])}</span>
        <div class="small">${esc(f.significance)}</div>
        <div><span class="stars" title="${esc(f.review_status)}">${"★".repeat(f.stars)}${"☆".repeat(4 - f.stars)}</span>
          <a class="small" href="${esc(f.url)}" target="_blank" rel="noopener">ClinVar ↗</a></div>
      </td>
      <td>${f.conditions.length ? esc(f.conditions.slice(0, 4).join("; ")) + (f.conditions.length > 4 ? ` <span class="small">+${f.conditions.length - 4} more</span>` : "") : `<span class="muted">—</span>`}</td>
      <td>${freqCell(f)}</td>
    </tr>`).join("") || `<tr><td colspan="6" class="muted">No findings match these filters.</td></tr>`;
  $("#more").hidden = rows.length <= shown;
}

// ---------------------------------------------------------------- export

function download(name, text, type) {
  const a = document.createElement("a");
  a.href = URL.createObjectURL(new Blob([text], { type }));
  a.download = name;
  a.click();
  setTimeout(() => URL.revokeObjectURL(a.href), 1000);
}

function exportCsv() {
  const head = ["category", "gene", "rsid", "chrom", "pos", "ref", "alt", "genotype", "zygosity", "significance",
    "review_status", "stars", "conditions", "af_all", ...POPS.map((p) => "af_" + p.toLowerCase()), "clinvar_url"];
  const cell = (v) => { const s = String(v ?? ""); return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s; };
  const lines = filtered().map((f) => [f.category, f.gene, f.rsid, f.chrom, f.pos, f.ref, f.alt, f.genotype, f.zygosity,
    f.significance, f.review_status, f.stars, f.conditions.join("; "), f.frequencies?.ALL,
    ...POPS.map((p) => f.frequencies?.[p]), f.url].map(cell).join(","));
  download("clinvar_findings.csv", [head.join(","), ...lines].join("\n"), "text/csv");
}

// ---------------------------------------------------------------- wiring

const drop = $("#drop");
$("#file").addEventListener("change", (e) => e.target.files[0] && upload(e.target.files[0]));
["dragenter", "dragover"].forEach((ev) => drop.addEventListener(ev, (e) => { e.preventDefault(); drop.classList.add("over"); }));
["dragleave", "drop"].forEach((ev) => drop.addEventListener(ev, (e) => { e.preventDefault(); drop.classList.remove("over"); }));
drop.addEventListener("drop", (e) => e.dataTransfer.files[0] && upload(e.dataTransfer.files[0]));

["#search", "#cat-filter", "#star-filter", "#zyg-filter"].forEach((s) =>
  $(s).addEventListener("input", () => { shown = PAGE; renderFindings(); }));
$("#more").addEventListener("click", () => { shown += PAGE; renderFindings(); });
$("#export-csv").addEventListener("click", exportCsv);
$("#export-json").addEventListener("click", () => download("genome_report.json", JSON.stringify(report, null, 2), "application/json"));
$("#new-upload").addEventListener("click", () => {
  history.replaceState(null, "", location.pathname);
  $("#results-view").hidden = true;
  $("#upload-view").hidden = false;
  $("#file").value = "";
});
$("#delete-report").addEventListener("click", async () => {
  if (!confirm("Permanently delete this report from the server?")) return;
  await fetch(`/api/jobs/${report.job_id}`, { method: "DELETE" });
  $("#new-upload").click();
});

loadStatus();
const m = location.hash.match(/^#report=([\w-]+)$/);
if (m) openReport(m[1]);
