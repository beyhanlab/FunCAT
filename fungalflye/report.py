"""
funcat.report
~~~~~~~~~~~~~~~~~
Generates a single self-contained HTML report for a completed assembly.
No internet connection required to view — everything is embedded.

Usage (from pipeline):
    from .report import generate_report
    generate_report(
        fasta        = "out/final.fasta",
        outdir       = "out",
        run_metadata = {...},
        telo_df      = df_or_None,
        confidence_tsv = "out/confidence/contig_confidence.tsv" or None,
    )

Or standalone:
    funcat report --fasta final.fasta --outdir out/
"""

import json
import math
import datetime
from pathlib import Path
from Bio import SeqIO


# ------------------------------------------------
# Data collectors
# ------------------------------------------------

def _collect_assembly_stats(fasta):
    lengths, gcs = [], []
    for r in SeqIO.parse(str(fasta), "fasta"):
        s = str(r.seq).upper()
        lengths.append(len(s))
        gc = (s.count("G") + s.count("C")) / len(s) * 100 if len(s) else 0
        gcs.append(round(gc, 2))

    if not lengths:
        return {}

    total = sum(lengths)
    sl = sorted(lengths, reverse=True)
    cumsum, n50, l50 = 0, 0, 0
    for i, L in enumerate(sl, 1):
        cumsum += L
        if cumsum >= total / 2:
            n50, l50 = L, i
            break

    return {
        "n_contigs":   len(lengths),
        "total_bp":    total,
        "largest_bp":  max(lengths),
        "n50_bp":      n50,
        "l50":         l50,
        "mean_gc":     round(sum(gcs) / len(gcs), 2),
        "lengths":     lengths,
        "gcs":         gcs,
    }


def _collect_confidence(tsv_path):
    if not tsv_path or not Path(tsv_path).exists():
        return []
    rows = []
    with open(tsv_path) as f:
        headers = None
        for line in f:
            parts = line.strip().split("\t")
            if headers is None:
                headers = parts
                continue
            rows.append(dict(zip(headers, parts)))
    return rows


def _collect_telomeres(telo_df):
    if telo_df is None:
        return []
    records = []
    for _, row in telo_df.iterrows():
        records.append({
            "contig":    row["contig"],
            "side":      row["side"],
            "telomeric": row["telomeric"],
            "hits":      int(row["hits"]) if str(row["hits"]).isdigit() else 0,
            "max_rep":   int(row["max_consecutive_repeats"]),
        })
    return records


# ------------------------------------------------
# HTML builder
# ------------------------------------------------

def _fmt_bp(bp):
    bp = int(bp)
    if bp >= 1_000_000:
        return f"{bp/1_000_000:.2f} Mb"
    if bp >= 1_000:
        return f"{bp/1_000:.1f} kb"
    return f"{bp} bp"


def _build_html(stats, confidence_rows, telo_records, metadata):
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    assembly_name = metadata.get("assembly_name", "Assembly")
    ploidy        = metadata.get("ploidy", "unknown")
    read_type     = metadata.get("read_type", "unknown")
    enhancements  = metadata.get("enhancements", {})

    # Contig length histogram buckets (log scale)
    lengths = stats.get("lengths", [])
    if lengths:
        max_len = max(lengths)
        n_bins  = 30
        bin_size = max_len / n_bins
        bins    = [0] * n_bins
        labels  = []
        for L in lengths:
            idx = min(int(L / bin_size), n_bins - 1)
            bins[idx] += 1
        for i in range(n_bins):
            lo = int(i * bin_size)
            hi = int((i + 1) * bin_size)
            labels.append(_fmt_bp((lo + hi) // 2))
        hist_data    = json.dumps(bins)
        hist_labels  = json.dumps(labels)
    else:
        hist_data = "[]"
        hist_labels = "[]"

    # GC scatter data (contig index vs GC)
    gc_data = json.dumps([
        {"x": i + 1, "y": g} for i, g in enumerate(stats.get("gcs", []))
    ])

    # Contig lengths for scatter (sorted descending)
    sorted_lens = sorted(lengths, reverse=True)
    cumsum_data = []
    running = 0
    for i, L in enumerate(sorted_lens):
        running += L
        cumsum_data.append({"x": i + 1, "y": running})
    cumsum_json = json.dumps(cumsum_data)
    total_bp = stats.get("total_bp", 1)
    n50_line = stats.get("n50_bp", 0)

    # Telomere grid
    telo_by_contig = {}
    for t in telo_records:
        key = t["contig"]
        telo_by_contig.setdefault(key, {})
        telo_by_contig[key][t["side"]] = t["telomeric"] == "YES"

    telo_html = ""
    if telo_records:
        telo_html = '<div class="telo-grid">'
        for ctg, sides in telo_by_contig.items():
            start = sides.get("start", False)
            end   = sides.get("end",   False)
            both  = start and end
            cls   = "telo-both" if both else ("telo-one" if (start or end) else "telo-none")
            label = "◀▶ complete" if both else ("◀ start" if start else ("▶ end" if end else "none"))
            short = ctg[:18] + "…" if len(ctg) > 18 else ctg
            telo_html += f'<div class="telo-card {cls}" title="{ctg}"><span class="telo-name">{short}</span><span class="telo-label">{label}</span></div>'
        telo_html += "</div>"
        telo_complete = sum(1 for s in telo_by_contig.values() if s.get("start") and s.get("end"))
        telo_summary  = f"{telo_complete} / {len(telo_by_contig)} chromosomes fully capped"
    else:
        telo_html    = '<p class="muted">Telomere analysis was not run.</p>'
        telo_summary = "—"

    # Confidence table
    if confidence_rows:
        conf_html = '<table class="conf-table"><thead><tr><th>Contig</th><th>Length</th><th>Mean cov</th><th>CV</th><th>Label</th><th>Note</th></tr></thead><tbody>'
        for r in confidence_rows:
            label = r.get("label", "")
            cls   = {"GOOD": "good", "REVIEW": "review", "FLAG": "flag"}.get(label, "")
            conf_html += (
                f'<tr class="{cls}">'
                f'<td class="mono">{r.get("contig","")}</td>'
                f'<td>{_fmt_bp(r.get("length_bp",0))}</td>'
                f'<td>{r.get("mean_coverage","")}</td>'
                f'<td>{r.get("coverage_cv","")}</td>'
                f'<td><span class="badge {cls}">{label}</span></td>'
                f'<td class="note">{r.get("reason","")}</td>'
                f'</tr>'
            )
        conf_html += "</tbody></table>"
        n_flag   = sum(1 for r in confidence_rows if r.get("label") == "FLAG")
        n_review = sum(1 for r in confidence_rows if r.get("label") == "REVIEW")
        conf_summary = f"{len(confidence_rows)} contigs scored — {n_flag} flagged, {n_review} for review"
    else:
        conf_html    = '<p class="muted">Confidence scoring was not run.</p>'
        conf_summary = "—"

    # Enhancement badges
    enh_badges = ""
    labels_map = {
        "adaptive_params":    "Adaptive parameters",
        "iterative_polish":   "Iterative polishing",
        "purge_dups":         "Purge duplicates",
        "scaffolding":        "Scaffolding",
        "confidence_scoring": "Confidence scoring",
    }
    for k, v in enhancements.items():
        name = labels_map.get(k, k)
        cls  = "enh-on" if v else "enh-off"
        enh_badges += f'<span class="enh-badge {cls}">{name}</span>'

    assembly_verdict = ""
    n_ctg = stats.get("n_contigs", 999)
    if n_ctg <= 12:
        assembly_verdict = '<div class="verdict good-verdict">✅ Assembly appears chromosome-level complete</div>'
    elif n_ctg <= 30:
        assembly_verdict = '<div class="verdict review-verdict">⚠️ Assembly is near-complete — a few extra fragments present</div>'
    else:
        assembly_verdict = '<div class="verdict flag-verdict">⚠️ Assembly is fragmented — consider tuning parameters or enabling scaffolding</div>'

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>FunCAT Report — {assembly_name}</title>
<script src="https://cdnjs.cloudflare.com/ajax/libs/Chart.js/4.4.1/chart.umd.min.js"></script>
<style>
:root {{
  --bg: #0f1117; --surface: #1a1d27; --surface2: #22263a;
  --border: #2e3348; --text: #e2e8f0; --muted: #7a83a6;
  --green: #22c55e; --amber: #f59e0b; --red: #ef4444;
  --teal: #14b8a6; --purple: #8b5cf6; --blue: #3b82f6;
  --font: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
  --mono: "SF Mono", "Fira Code", monospace;
}}
* {{ box-sizing: border-box; margin: 0; padding: 0; }}
body {{ background: var(--bg); color: var(--text); font-family: var(--font);
       font-size: 14px; line-height: 1.6; }}
.page {{ max-width: 1100px; margin: 0 auto; padding: 32px 24px; }}
.header {{ border-bottom: 1px solid var(--border); padding-bottom: 24px; margin-bottom: 32px; }}
.header h1 {{ font-size: 26px; font-weight: 600; color: var(--text); }}
.header h1 span {{ color: var(--teal); }}
.meta {{ color: var(--muted); font-size: 13px; margin-top: 6px; }}
.meta b {{ color: var(--text); }}
.grid-4 {{ display: grid; grid-template-columns: repeat(4,1fr); gap: 16px; margin-bottom: 32px; }}
.stat-card {{ background: var(--surface); border: 1px solid var(--border);
              border-radius: 12px; padding: 20px; }}
.stat-card .label {{ font-size: 11px; text-transform: uppercase; letter-spacing: .08em;
                     color: var(--muted); margin-bottom: 6px; }}
.stat-card .value {{ font-size: 24px; font-weight: 600; color: var(--teal); }}
.stat-card .sub {{ font-size: 12px; color: var(--muted); margin-top: 2px; }}
.section {{ background: var(--surface); border: 1px solid var(--border);
            border-radius: 12px; padding: 24px; margin-bottom: 24px; }}
.section h2 {{ font-size: 15px; font-weight: 600; margin-bottom: 16px;
               padding-bottom: 10px; border-bottom: 1px solid var(--border); }}
.charts-row {{ display: grid; grid-template-columns: 1fr 1fr; gap: 24px; }}
.chart-wrap {{ position: relative; height: 220px; }}
.verdict {{ padding: 12px 16px; border-radius: 8px; font-weight: 500;
            font-size: 13px; margin-bottom: 20px; }}
.good-verdict   {{ background: rgba(34,197,94,.12); border:1px solid rgba(34,197,94,.3); color: var(--green); }}
.review-verdict {{ background: rgba(245,158,11,.12); border:1px solid rgba(245,158,11,.3); color: var(--amber); }}
.flag-verdict   {{ background: rgba(239,68,68,.12);  border:1px solid rgba(239,68,68,.3);  color: var(--red);   }}
.telo-grid {{ display: flex; flex-wrap: wrap; gap: 8px; }}
.telo-card {{ display: flex; flex-direction: column; align-items: center;
              padding: 10px 14px; border-radius: 8px; min-width: 120px;
              border: 1px solid var(--border); }}
.telo-both {{ background: rgba(34,197,94,.12); border-color: rgba(34,197,94,.3); }}
.telo-one  {{ background: rgba(245,158,11,.10); border-color: rgba(245,158,11,.3); }}
.telo-none {{ background: var(--surface2); }}
.telo-name  {{ font-family: var(--mono); font-size: 11px; color: var(--muted); }}
.telo-label {{ font-size: 12px; font-weight: 500; margin-top: 4px; }}
.telo-both .telo-label {{ color: var(--green); }}
.telo-one  .telo-label {{ color: var(--amber); }}
.telo-none .telo-label {{ color: var(--muted); }}
.conf-table {{ width: 100%; border-collapse: collapse; font-size: 13px; }}
.conf-table th {{ text-align: left; padding: 8px 12px; font-weight: 500;
                  color: var(--muted); border-bottom: 1px solid var(--border);
                  font-size: 11px; text-transform: uppercase; letter-spacing: .06em; }}
.conf-table td {{ padding: 8px 12px; border-bottom: 1px solid rgba(46,51,72,.5); }}
.conf-table tr:last-child td {{ border-bottom: none; }}
.conf-table tr.flag   {{ background: rgba(239,68,68,.06); }}
.conf-table tr.review {{ background: rgba(245,158,11,.06); }}
.badge {{ display: inline-block; padding: 2px 8px; border-radius: 20px;
          font-size: 11px; font-weight: 600; letter-spacing: .04em; }}
.badge.good   {{ background: rgba(34,197,94,.15);  color: var(--green); }}
.badge.review {{ background: rgba(245,158,11,.15); color: var(--amber); }}
.badge.flag   {{ background: rgba(239,68,68,.15);  color: var(--red);   }}
.mono {{ font-family: var(--mono); font-size: 12px; }}
.note {{ color: var(--muted); font-size: 12px; max-width: 340px; }}
.enh-badge {{ display: inline-block; padding: 3px 10px; border-radius: 20px;
              font-size: 12px; margin: 3px; }}
.enh-on  {{ background: rgba(20,184,166,.15); color: var(--teal);
            border: 1px solid rgba(20,184,166,.3); }}
.enh-off {{ background: var(--surface2); color: var(--muted);
            border: 1px solid var(--border); text-decoration: line-through; }}
.muted {{ color: var(--muted); font-size: 13px; }}
.share-box {{ background: var(--surface2); border: 1px solid var(--border);
              border-radius: 8px; padding: 14px 18px; font-size: 13px;
              color: var(--muted); margin-top: 24px; }}
.share-box b {{ color: var(--text); }}
@media (max-width: 700px) {{
  .grid-4 {{ grid-template-columns: 1fr 1fr; }}
  .charts-row {{ grid-template-columns: 1fr; }}
}}
</style>
</head>
<body>
<div class="page">

  <div class="header">
    <h1>🧬 FunCAT Report — <span>{assembly_name}</span></h1>
    <div class="meta">
      Generated <b>{now}</b> &nbsp;·&nbsp;
      Ploidy <b>{ploidy}</b> &nbsp;·&nbsp;
      Reads <b>{read_type}</b> &nbsp;·&nbsp;
      FunCAT v0.3
    </div>
  </div>

  {assembly_verdict}

  <div class="grid-4">
    <div class="stat-card">
      <div class="label">Contigs</div>
      <div class="value">{stats.get("n_contigs","—")}</div>
      <div class="sub">scaffolds / chromosomes</div>
    </div>
    <div class="stat-card">
      <div class="label">Total size</div>
      <div class="value">{_fmt_bp(stats.get("total_bp",0))}</div>
      <div class="sub">{stats.get("total_bp",0):,} bp</div>
    </div>
    <div class="stat-card">
      <div class="label">N50</div>
      <div class="value">{_fmt_bp(stats.get("n50_bp",0))}</div>
      <div class="sub">L50 = {stats.get("l50","—")} contigs</div>
    </div>
    <div class="stat-card">
      <div class="label">Mean GC</div>
      <div class="value">{stats.get("mean_gc","—")}%</div>
      <div class="sub">Largest: {_fmt_bp(stats.get("largest_bp",0))}</div>
    </div>
  </div>

  <div class="section">
    <h2>Contig length distribution &amp; cumulative assembly</h2>
    <div class="charts-row">
      <div class="chart-wrap"><canvas id="histChart"></canvas></div>
      <div class="chart-wrap"><canvas id="cumulChart"></canvas></div>
    </div>
  </div>

  <div class="section">
    <h2>GC content per contig</h2>
    <div class="chart-wrap" style="height:180px"><canvas id="gcChart"></canvas></div>
  </div>

  <div class="section">
    <h2>Telomere analysis &nbsp;<span class="muted" style="font-weight:400;font-size:13px">{telo_summary}</span></h2>
    {telo_html}
  </div>

  <div class="section">
    <h2>Contig confidence &nbsp;<span class="muted" style="font-weight:400;font-size:13px">{conf_summary}</span></h2>
    {conf_html}
  </div>

  <div class="section">
    <h2>Enhancement modules</h2>
    <div>{enh_badges if enh_badges else '<span class="muted">No enhancement metadata available.</span>'}</div>
  </div>

  <div class="share-box">
    <b>This report is fully self-contained.</b> No internet connection required to view it.
    Share it by sending this single HTML file — it works offline in any modern browser.
  </div>

</div>

<script>
const CHART_DEFAULTS = {{
  color: '#7a83a6',
  borderColor: '#2e3348',
  plugins: {{ legend: {{ display: false }}, tooltip: {{ callbacks: {{}} }} }},
  scales: {{
    x: {{ ticks: {{ color: '#7a83a6', maxRotation: 0, maxTicksLimit: 8 }},
          grid:  {{ color: 'rgba(46,51,72,0.5)' }} }},
    y: {{ ticks: {{ color: '#7a83a6' }},
          grid:  {{ color: 'rgba(46,51,72,0.5)' }} }},
  }},
}};

function deepMerge(base, over) {{
  const out = Object.assign({{}}, base);
  for (const k of Object.keys(over)) {{
    if (typeof over[k] === 'object' && !Array.isArray(over[k]) && over[k] !== null)
      out[k] = deepMerge(base[k] || {{}}, over[k]);
    else out[k] = over[k];
  }}
  return out;
}}

// Contig length histogram
new Chart(document.getElementById('histChart'), {{
  type: 'bar',
  data: {{
    labels: {hist_labels},
    datasets: [{{ data: {hist_data}, backgroundColor: 'rgba(20,184,166,0.6)',
                  borderColor: '#14b8a6', borderWidth: 1 }}],
  }},
  options: deepMerge(CHART_DEFAULTS, {{
    plugins: {{ tooltip: {{ callbacks: {{ title: (i) => 'around ' + i[0].label }} }} }},
    scales: {{ y: {{ title: {{ display: true, text: 'contigs', color: '#7a83a6' }} }} }},
  }}),
}});

// Cumulative assembly (Lx plot)
const totalBp = {total_bp};
const n50Line = {n50_line};
new Chart(document.getElementById('cumulChart'), {{
  type: 'line',
  data: {{
    datasets: [{{
      data: {cumsum_json},
      borderColor: '#8b5cf6', backgroundColor: 'rgba(139,92,246,0.08)',
      borderWidth: 2, pointRadius: 0, fill: true,
      tension: 0.3,
    }}],
  }},
  options: deepMerge(CHART_DEFAULTS, {{
    parsing: {{ xAxisKey: 'x', yAxisKey: 'y' }},
    scales: {{
      x: {{ title: {{ display: true, text: 'contigs (sorted by length)', color: '#7a83a6' }} }},
      y: {{ title: {{ display: true, text: 'cumulative bp', color: '#7a83a6' }},
            ticks: {{ callback: v => v >= 1e6 ? (v/1e6).toFixed(1)+'M' : v >= 1e3 ? (v/1e3).toFixed(0)+'k' : v }} }},
    }},
    plugins: {{
      annotation: {{
        annotations: {{
          n50: {{ type: 'line', yMin: totalBp/2, yMax: totalBp/2,
                  borderColor: 'rgba(245,158,11,0.6)', borderWidth: 1,
                  borderDash: [4,4] }},
        }}
      }}
    }},
  }}),
}});

// GC scatter
new Chart(document.getElementById('gcChart'), {{
  type: 'scatter',
  data: {{
    datasets: [{{ data: {gc_data}, backgroundColor: 'rgba(59,130,246,0.5)',
                  pointRadius: 3, pointHoverRadius: 5 }}],
  }},
  options: deepMerge(CHART_DEFAULTS, {{
    parsing: {{ xAxisKey: 'x', yAxisKey: 'y' }},
    scales: {{
      x: {{ title: {{ display: true, text: 'contig index', color: '#7a83a6' }} }},
      y: {{ title: {{ display: true, text: 'GC %', color: '#7a83a6' }},
            min: 0, max: 100 }},
    }},
  }}),
}});
</script>
</body>
</html>"""

    return html


# ------------------------------------------------
# Public API
# ------------------------------------------------

def generate_report(
    fasta,
    outdir,
    run_metadata=None,
    telo_df=None,
    confidence_tsv=None,
):
    """
    Generate a self-contained HTML assembly report.

    Parameters
    ----------
    fasta           : path to final assembly FASTA
    outdir          : directory to write report into
    run_metadata    : dict with keys: assembly_name, ploidy, read_type, enhancements
    telo_df         : pandas DataFrame from qc.scan_telomeres(), or None
    confidence_tsv  : path to confidence/contig_confidence.tsv, or None

    Returns
    -------
    Path to generated HTML file.
    """

    fasta  = Path(fasta)
    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True)

    if run_metadata is None:
        run_metadata = {}

    metadata = {
        "assembly_name": run_metadata.get("assembly_name", fasta.stem),
        "ploidy":        run_metadata.get("ploidy",        "unknown"),
        "read_type":     run_metadata.get("read_type",     "unknown"),
        "enhancements":  run_metadata.get("enhancements",  {}),
    }

    print("\n[funcat] Generating HTML report...")

    stats            = _collect_assembly_stats(fasta)
    confidence_rows  = _collect_confidence(confidence_tsv)
    telo_records     = _collect_telomeres(telo_df)

    html = _build_html(stats, confidence_rows, telo_records, metadata)

    out_path = outdir / f"{metadata['assembly_name']}_funcat_report.html"
    out_path.write_text(html, encoding="utf-8")

    print(f"[funcat] ✅ Report saved: {out_path}")
    print(f"             Open in any browser — no internet required\n")

    return out_path
