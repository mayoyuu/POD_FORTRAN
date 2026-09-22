#!/usr/bin/env python3
"""Render selected ADS error-domain hours and sampled size history to one SVG.

Usage:
  python3 tools/plot_srp_ads_history.py --prefix output/orbit_srp_ads --hours 0 24 360
The input CSV ranges and principal axes are sampled estimates, not interval bounds.
"""
import argparse
import csv
import json
import math
from pathlib import Path
from xml.sax.saxutils import escape


def read_rows(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def line(x1, y1, x2, y2, color="#718096", width=1, dash=""):
    more = f' stroke-dasharray="{dash}"' if dash else ""
    return (
        f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" '
        f'y2="{y2:.2f}" stroke="{color}" stroke-width="{width}"{more}/>'
    )


def text(x, y, content, size=13, color="#263238", anchor="start"):
    return (
        f'<text x="{x:.2f}" y="{y:.2f}" fill="{color}" font-size="{size}" '
        f'text-anchor="{anchor}" font-family="Arial,sans-serif">{escape(str(content))}</text>'
    )


def scatter_panel(rows, hour, x0, y0, width, height):
    center_x = x0 + width / 2
    center_y = y0 + height / 2
    projected = []
    for row in rows:
        x = float(row["dx_km"])
        y = float(row["dy_km"])
        z = float(row["dz_km"])
        # Fixed isometric view of the three physical position errors.
        projected.append(((x - y) / math.sqrt(2), (x + y - 2 * z) / math.sqrt(6), z))
    extent = max((max(abs(px), abs(py)) for px, py, _ in projected), default=0.0)
    extent = max(extent, 1e-12)
    scale = 0.43 * min(width, height) / extent
    out = [
        f'<rect x="{x0}" y="{y0}" width="{width}" height="{height}" '
        'fill="#f8fafc" stroke="#cbd5e1"/>',
        text(x0 + 14, y0 + 22, f"hour {hour:g}: position error domain", 14),
        line(center_x - 0.43 * width, center_y, center_x + 0.43 * width, center_y, "#d3dbe3"),
        line(center_x, center_y - 0.40 * height, center_x, center_y + 0.40 * height, "#d3dbe3"),
    ]
    zmin = min((p[2] for p in projected), default=0.0)
    zmax = max((p[2] for p in projected), default=0.0)
    for px, py, pz in projected:
        fraction = (pz - zmin) / (zmax - zmin) if zmax > zmin else 0.5
        red = round(35 + 185 * fraction)
        blue = round(195 - 120 * fraction)
        color = f"#{red:02x}64{blue:02x}"
        out.append(
            f'<circle cx="{center_x + px * scale:.2f}" '
            f'cy="{center_y - py * scale:.2f}" r="2.1" fill="{color}" fill-opacity=".67"/>'
        )
    out.append(text(x0 + 12, y0 + height - 10, "isometric x/y/z projection; values in km", 11, "#64748b"))
    return "\n".join(out)


def size_panel(history, x0, y0, width, height):
    left, right = x0 + 58, x0 + width - 18
    top, bottom = y0 + 36, y0 + height - 45
    hours = [float(row["hour"]) for row in history]
    position = [float(row["pos_max_km"]) for row in history]
    axes = [[float(row[f"pos_axis{k}_km"]) for row in history] for k in (1, 2, 3)]
    hmax = max(hours) if hours else 1.0
    vmax = max(position + sum(axes, [])) if position else 1.0
    hmax = max(hmax, 1e-12)
    vmax = max(vmax, 1e-12)
    def coordinates(values):
        return " ".join(
            f"{left + (h / hmax) * (right - left):.2f},"
            f"{bottom - (v / vmax) * (bottom - top):.2f}"
            for h, v in zip(hours, values)
        )
    out = [
        f'<rect x="{x0}" y="{y0}" width="{width}" height="{height}" '
        'fill="#ffffff" stroke="#cbd5e1"/>',
        text(x0 + 14, y0 + 22, "Sampled position-domain size over time", 14),
        line(left, top, left, bottom, "#475569"),
        line(left, bottom, right, bottom, "#475569"),
        text(left, bottom + 20, "0", 11, "#64748b"),
        text(right, bottom + 20, f"{hmax:g} h", 11, "#64748b", "end"),
        text(left - 6, top + 2, f"{vmax:.3g} km", 11, "#64748b", "end"),
        text(left - 6, bottom + 2, "0", 11, "#64748b", "end"),
    ]
    series = [(position, "#be123c", "max |dr|")] + list(zip(
        axes, ("#0f766e", "#2563eb", "#8b5cf6"), ("axis 1", "axis 2", "axis 3")
    ))
    for index, (values, color, label) in enumerate(series):
        out.append(f'<polyline points="{coordinates(values)}" fill="none" '
                   f'stroke="{color}" stroke-width="2"/>')
        out.append(line(left + index * 145, y0 + height - 15,
                        left + index * 145 + 18, y0 + height - 15, color, 2))
        out.append(text(left + index * 145 + 23, y0 + height - 11, label, 11))
    return "\n".join(out)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", required=True, help="output prefix used by run_srp_ads_history")
    parser.add_argument("--hours", nargs="+", type=float, help="saved hours to show (default: start and end)")
    parser.add_argument("--output", help="SVG path (default: PREFIX_domain.svg)")
    args = parser.parse_args()
    prefix = Path(args.prefix)
    history = read_rows(str(prefix) + "_history.csv")
    if not history:
        parser.error("history CSV is empty")
    desired = args.hours if args.hours is not None else [0.0, float(history[-1]["hour"])]
    selected = []
    for hour in desired:
        if all(abs(float(row["hour"]) - hour) > 1e-8 for row in history):
            parser.error(f"hour {hour:g} is not a saved checkpoint")
        if all(abs(hour - existing) > 1e-8 for existing in selected):
            selected.append(hour)
    selected = selected[:8]
    shape = read_rows(str(prefix) + "_shape.csv")
    report_path = Path(str(prefix) + "_report.json")
    report = json.loads(report_path.read_text(encoding="utf-8")) if report_path.exists() else {}
    panel_width, panel_height = 620, 310
    columns = 2
    rows = math.ceil(len(selected) / columns)
    canvas_width = 1280
    canvas_height = 395 + rows * 330
    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{canvas_width}" '
        f'height="{canvas_height}" viewBox="0 0 {canvas_width} {canvas_height}">',
        '<rect width="100%" height="100%" fill="#f1f5f9"/>',
        text(24, 32, "SRP–ADS sampled error domain", 21),
        text(24, 54, "Ranges and axes are estimates from 512 fixed domain points.", 12, "#64748b"),
        size_panel(history, 20, 72, 1240, 290),
    ]
    if report.get("first_failed_hour") is not None:
        out.append(text(24, 381, f"First observed validation failure: hour {report['first_failed_hour']:g}", 12, "#b91c1c"))
    for index, hour in enumerate(selected):
        rows_at_hour = [row for row in shape if abs(float(row["hour"]) - hour) < 1e-8]
        out.append(scatter_panel(rows_at_hour, hour, 20 + (index % 2) * 640,
                                 395 + (index // 2) * 330, panel_width, panel_height))
    out.append("</svg>")
    destination = Path(args.output) if args.output else Path(str(prefix) + "_domain.svg")
    destination.write_text("\n".join(out), encoding="utf-8")
    print(destination)


if __name__ == "__main__":
    main()
