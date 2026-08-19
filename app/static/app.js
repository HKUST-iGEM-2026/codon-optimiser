"use strict";
// Codon track and result rendering. No framework, no build step.
// Every number shown here comes from /api/v1; nothing is computed twice.

const EXAMPLE_CDS =
  "ATGTCAAAGCAAAAAATAAGAATAAGATTAAAAGCTTTTGATCACACAATATTAGATCAATCTGCTGAAAAAATCGTTGAAACT" +
  "GCAAAAACTTCAGGAGCTAAGGTTGTAGGTCCAGTACCATTACCAACTGAAAAAGATGTTGTTACAATATTAAGAGCGGTTCAC" +
  "AAATACAAAGATTCAAGAGAACAATTTGAAATCAGAACTCATAAGAGATTAATCGATATCGTTAATCCATCACCAAAAACTGTT" +
  "GATGCTTTAATGAGATTAAATCTTCCAGCTGGTGTTGATATAGAAATCAAACTATAA"; // rpsJ, K670_RS00880

const SCALE_LOW = [242, 239, 228];
const SCALE_HIGH = [15, 81, 50];
const BASES_PER_BLOCK = 10;
const BLOCKS_PER_LINE = 6;

const $ = (id) => document.getElementById(id);
const el = (tag, className, text) => {
  const node = document.createElement(tag);
  if (className) node.className = className;
  if (text !== undefined) node.textContent = text;
  return node;
};

let state = { mode: null, colour: "adaptiveness", result: null };

/* --- colour scale ------------------------------------------------------- */
function scaleColour(weight) {
  const t = Math.max(0, Math.min(1, weight));
  const mix = SCALE_LOW.map((low, i) => Math.round(low + (SCALE_HIGH[i] - low) * t));
  return { css: `rgb(${mix.join(",")})`, dark: t > 0.55 };
}

/* --- request ----------------------------------------------------------- */
function payload(extra) {
  const avoid = $("avoid").value.split(/[,\s]+/).map((m) => m.trim().toUpperCase()).filter(Boolean);
  return Object.assign(
    {
      sequence: $("sequence").value,
      mode: state.mode,
      seed: Number($("seed").value),
      avoid,
      upstream: $("upstream").value.trim().toUpperCase(),
      downstream: $("downstream").value.trim().toUpperCase(),
    },
    extra || {}
  );
}

// An unhandled server error answers in plain text, not JSON; parsing it blind
// turns a readable status into "Unexpected token 'I'".
async function parseBody(response) {
  const text = await response.text();
  try {
    return JSON.parse(text);
  } catch {
    return { error: "unreadable_response", detail: `HTTP ${response.status}: ${text.slice(0, 200)}` };
  }
}

async function optimise() {
  clearErrors();
  const body = payload({ permalink: $("permalink").checked });
  if (!body.sequence.trim()) return showInputError("paste a protein or coding sequence first");
  setBusy(true);
  try {
    const response = await fetch("/api/v1/optimise", {
      method: "POST",
      headers: { "content-type": "application/json" },
      body: JSON.stringify(body),
    });
    const data = await parseBody(response);
    if (response.ok) return render(data);
    showFailure(response.status, data);
  } catch (err) {
    showInputError(`request failed: ${err.message}`);
  } finally {
    setBusy(false);
  }
}

function setBusy(busy) {
  $("busy").hidden = !busy;
  $("optimise").disabled = busy;
}

function clearErrors() {
  $("input-error").textContent = "";
  $("failure").hidden = true;
}

function showInputError(message) {
  $("input-error").textContent = message;
  $("result").hidden = true;
}

function showFailure(status, data) {
  $("result").hidden = true;
  if (status === 400) return showInputError(data.detail || "invalid sequence");
  if (status === 422 && data.error === "invalid_request") {
    return showInputError(`invalid request field: ${(data.fields || []).join(", ")}`);
  }
  if (status === 429) {
    return showInputError(
      `CPU budget exhausted (${data.budget_cpu_seconds}s per window), retry in ${data.retry_after_seconds}s`
    );
  }
  $("failure").hidden = false;
  $("failure-title").textContent =
    data.error === "solver_gave_up"
      ? "No sequence returned — the search ran out of attempts"
      : "No sequence returned — hard constraints conflict";
  $("failure-detail").textContent = data.hint || data.detail || "";
  const list = $("conflict");
  list.replaceChildren();
  (data.conflict || []).forEach((c) => {
    const li = el("li");
    li.append(el("span", "coord", `${c.start}-${c.end}`), el("span", "flag", c.constraint), el("span", "", c.detail));
    list.append(li);
  });
}

/* --- rendering --------------------------------------------------------- */
function render(data) {
  state.result = data;
  clearErrors();
  $("result").hidden = false;

  const m = data.metrics;
  const cards = [
    ["CAI", m.cai.toFixed(4), m.cai_of_input === null ? "protein input" : `input CAI ${m.cai_of_input.toFixed(4)}`],
    ["GC", `${m.gc_percent}%`, `band ${data.manifest.parameters.gc_global_band_percent.join("–")}%`],
    ["GC 50nt window", `${m.gc_window_min}–${m.gc_window_max}%`, `ceiling ${data.manifest.parameters.gc_window_band_percent[1]}%`],
    ["Changed codons", `${m.changed_codons} / ${data.codons.length}`, `${m.changed_by_constraint} forced by a constraint`],
    ["Changed after codon 12", `${m.changed_codons_after_trim}`, "5′ convention of the tables"],
    ["5′ ΔG", Number(m.five_prime.dg).toFixed(2), `kcal/mol, ${m.five_prime.upstream_nt_available}nt upstream`],
  ];
  const strip = $("result-metrics");
  strip.replaceChildren();
  cards.forEach(([label, value, sub]) => {
    const card = el("div", "metric");
    card.append(el("span", "label", label), el("span", "value", value), el("span", "sub", sub));
    strip.append(card);
  });
  $("result-note").textContent =
    `mode ${data.manifest.parameters.mode} · seed ${data.manifest.parameters.seed} · ` +
    `${m.length_nt} nt · input read as ${data.input_kind} · engine ${data.manifest.engine_version} · ` +
    `tables ${data.manifest.table_version}`;

  renderTrack(data.codons);
  renderSequence(data.sequence);
  renderList($("flags"), data.flags, (f) => [`${f.start}-${f.end}`, f.flag, f.detail], "no liabilities flagged");
  renderList($("warnings"), data.warnings.map((w) => ({ text: w })), (w) => [null, null, w.text], "no warnings");
  renderManifest(data.manifest);

  const out = $("permalink-out");
  out.textContent = data.permalink ? `${location.origin}/p/${data.permalink}` : "";
}

function renderTrack(codons) {
  const track = $("track");
  track.replaceChildren();
  track.className = `track ${state.colour === "change" ? "by-change" : ""}`;
  codons.forEach((codon) => {
    const cell = el("button", `codon ${codon.status}`, codon.chosen);
    const colour = scaleColour(codon.relative_adaptiveness ?? 1);
    if (state.colour === "adaptiveness") {
      cell.style.background = colour.css;
      if (colour.dark) cell.classList.add("dark");
    }
    cell.type = "button";
    cell.dataset.index = codon.index;
    cell.setAttribute("aria-expanded", "false");
    cell.setAttribute(
      "aria-label",
      `codon ${codon.index + 1} ${codon.aa}, chosen ${codon.chosen}, original ${codon.original}, ${codon.status}`
    );
    track.append(cell);
  });
}

function renderSequence(sequence) {
  const out = $("sequence-out");
  out.replaceChildren();
  const perLine = BASES_PER_BLOCK * BLOCKS_PER_LINE;
  for (let i = 0; i < sequence.length; i += perLine) {
    const line = sequence.slice(i, i + perLine);
    const blocks = [];
    for (let j = 0; j < line.length; j += BASES_PER_BLOCK) blocks.push(line.slice(j, j + BASES_PER_BLOCK));
    out.append(
      document.createTextNode(blocks.join(" ").padEnd(perLine + BLOCKS_PER_LINE - 1, " ") + "  "),
      el("span", "pos", String(i + line.length)),
      document.createTextNode("\n")
    );
  }
}

function renderList(node, items, fields, emptyText) {
  node.replaceChildren();
  if (!items.length) return node.append(el("li", "empty", emptyText));
  items.forEach((item) => {
    const [coord, flag, detail] = fields(item);
    const li = el("li");
    if (coord) li.append(el("span", "coord", coord));
    if (flag) li.append(el("span", "flag", flag));
    li.append(el("span", "", detail));
    node.append(li);
  });
}

function renderManifest(manifest) {
  const dl = $("manifest");
  dl.replaceChildren();
  const rows = Object.entries(manifest).filter(([k]) => k !== "parameters");
  rows.push(...Object.entries(manifest.parameters));
  rows.forEach(([key, value]) => {
    dl.append(el("dt", "", key), el("dd", "", Array.isArray(value) || typeof value === "object" ? JSON.stringify(value) : String(value)));
  });
}

/* --- tooltip ----------------------------------------------------------- */
function showTooltip(cell) {
  const codon = state.result.codons[Number(cell.dataset.index)];
  const tip = $("tooltip");
  const dl = el("dl");
  const rows = [
    ["codon", `${codon.index + 1}`],
    ["amino acid", codon.aa],
    ["original", `${codon.original}  ${fmt(codon.original_frequency)}`],
    ["chosen", `${codon.chosen}  ${fmt(codon.chosen_frequency)}`],
    ["adaptiveness", fmt(codon.relative_adaptiveness)],
    ["status", codon.status],
    ["rule", codon.rule ? `${codon.rule.rule} ${codon.rule.start}-${codon.rule.end}` : "—"],
    ["local ΔG", `${codon.dg_delta > 0 ? "+" : ""}${codon.dg_delta} kcal/mol`],
  ];
  rows.forEach(([k, v]) => dl.append(el("dt", "", k), el("dd", "", v)));
  tip.replaceChildren(dl);
  if (codon.rule) tip.append(el("p", "detail", codon.rule.detail));

  tip.hidden = false;
  const box = cell.getBoundingClientRect();
  const top = box.bottom + window.scrollY + 6;
  const left = Math.min(box.left + window.scrollX, window.innerWidth - 276);
  tip.style.top = `${top}px`;
  tip.style.left = `${Math.max(8, left)}px`;
  cell.setAttribute("aria-expanded", "true");
}

function fmt(value) {
  return value === null || value === undefined ? "—" : Number(value).toFixed(3);
}

function hideTooltip() {
  $("tooltip").hidden = true;
  document.querySelectorAll('.codon[aria-expanded="true"]').forEach((c) => c.setAttribute("aria-expanded", "false"));
}

/* --- copy and export --------------------------------------------------- */
// A rendered result stays valid when Copy or Export fails, so neither may report
// through showInputError: that hides the very sequence the user is trying to take.
// The button that was pressed says what went wrong instead.
function flash(id, message) {
  const button = $(id);
  const label = button.dataset.label || (button.dataset.label = button.textContent);
  button.textContent = message;
  setTimeout(() => (button.textContent = label), 2500);
}

async function copySequence() {
  try {
    await navigator.clipboard.writeText(state.result.sequence);
    flash("copy", "Copied");
  } catch {
    // navigator.clipboard does not exist outside a secure context, which is exactly
    // what `docker run -p 8000:8000` reached over a LAN address gives you.
    flash("copy", "Blocked — select manually");
  }
}

async function exportFasta() {
  setBusy(true);  // the export re-runs the solve server-side, so it is not instant
  try {
    const response = await fetch("/api/v1/export.fasta", {
      method: "POST",
      headers: { "content-type": "application/json" },
      body: JSON.stringify(payload()),
    });
    if (!response.ok) return showFailure(response.status, await parseBody(response));
    const blob = new Blob([await response.text()], { type: "text/plain" });
    const link = el("a");
    link.href = URL.createObjectURL(blob);
    link.download = "codon_optimised.fasta";
    link.click();
    URL.revokeObjectURL(link.href);
  } catch (err) {
    flash("fasta", `Export failed: ${err.message}`);
  } finally {
    setBusy(false);
  }
}

/* --- wiring ------------------------------------------------------------ */
function setMode(mode) {
  state.mode = mode;
  document.querySelectorAll("[data-mode]").forEach((b) => b.setAttribute("aria-checked", String(b.dataset.mode === mode)));
}

function updateCounter() {
  const text = $("sequence").value.replace(/\s|^>.*$/gm, "");
  const isDna = text.length > 0 && !/[^ACGTUacgtu]/.test(text) && text.length % 3 === 0;
  const residues = isDna ? text.length / 3 : text.length;
  const nt = (residues + (isDna ? 0 : 1)) * 3;
  $("counter").textContent = `${residues} ${isDna ? "codons" : "aa"} · ${nt} / ${document.body.dataset.cap} nt output cap`;
}

document.addEventListener("DOMContentLoaded", () => {
  setMode(document.body.dataset.defaultMode);
  document.querySelectorAll("[data-mode]").forEach((b) => b.addEventListener("click", () => setMode(b.dataset.mode)));
  document.querySelectorAll("[data-colour]").forEach((b) =>
    b.addEventListener("click", () => {
      state.colour = b.dataset.colour;
      document.querySelectorAll("[data-colour]").forEach((x) => x.setAttribute("aria-checked", String(x.dataset.colour === state.colour)));
      if (state.result) renderTrack(state.result.codons);
    })
  );
  document.querySelector('[data-colour="adaptiveness"]').setAttribute("aria-checked", "true");

  $("optimise").addEventListener("click", optimise);
  $("example").addEventListener("click", () => {
    $("sequence").value = EXAMPLE_CDS;
    updateCounter();
  });
  $("sequence").addEventListener("input", updateCounter);
  $("copy").addEventListener("click", copySequence);
  $("fasta").addEventListener("click", exportFasta);

  const track = $("track");
  track.addEventListener("mouseover", (e) => {
    const cell = e.target.closest(".codon");
    if (cell) showTooltip(cell);
  });
  track.addEventListener("mouseleave", hideTooltip);
  track.addEventListener("focusin", (e) => {
    const cell = e.target.closest(".codon");
    if (cell) showTooltip(cell);
  });
  track.addEventListener("focusout", hideTooltip);
  document.addEventListener("keydown", (e) => e.key === "Escape" && hideTooltip());

  updateCounter();
  const permalink = document.body.dataset.permalink;
  if (permalink) loadPermalink(permalink);
});

async function loadPermalink(pid) {
  let response;
  try {
    response = await fetch(`/api/v1/permalink/${pid}`);
  } catch (err) {
    // Nothing is on screen yet on a /p/{id} page, so a dead fetch must not end here
    // as an unhandled rejection: the page would just stay blank.
    return showInputError(`could not load permalink: ${err.message}`);
  }
  const data = await parseBody(response);
  if (!response.ok) return showInputError(data.detail || "unknown permalink");
  $("sequence").value = data.input_kind === "cds" ? data.codons.map((c) => c.original).join("") : data.protein.replace(/\*$/, "");
  setMode(data.manifest.parameters.mode);
  $("seed").value = data.manifest.parameters.seed;
  $("avoid").value = data.manifest.parameters.avoid.join(", ");
  updateCounter();
  render(data);
  // The manifest carries context lengths, never the context itself (sequences are
  // ephemeral), so the two fields cannot be repopulated. Say so, because pressing
  // Optimise again without them silently produces a different result.
  const p = data.manifest.parameters;
  if (p.upstream_nt || p.downstream_nt) {
    $("result-note").textContent +=
      ` · assembled with ${p.upstream_nt}nt 5′ and ${p.downstream_nt}nt 3′ context, not shown: ` +
      "re-optimising from this page without re-entering it gives a different result";
  }
}
