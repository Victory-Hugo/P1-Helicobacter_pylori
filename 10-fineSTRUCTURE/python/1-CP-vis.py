#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import ast
import os
import sys
import math
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Force non-interactive backend to avoid Qt dependency
import matplotlib.pyplot as plt


# seaborn / aquarel are optional; script still runs without them
try:
    import seaborn as sns
    _HAS_SEABORN = True
except Exception:
    _HAS_SEABORN = False

def _try_apply_aquarel_theme():
    try:
        from aquarel import load_theme
        theme = load_theme('boxy_light')
        theme.apply()
        return theme
    except Exception:
        return None

def parse_list(s: str):
    """
    Parse string representations of lists:
    1) JSON / Python literal formats: '["a","b"]' or "['a','b']"
    2) Comma-separated strings: 'a,b,c'
    """
    if s is None:
        return []
    s = s.strip()
    if not s:
        return []
    # Try literal_eval first
    if (s.startswith('[') and s.endswith(']')) or (s.startswith('(') and s.endswith(')')):
        try:
            val = ast.literal_eval(s)
            if isinstance(val, (list, tuple)):
                return [str(x).strip() for x in val]
        except Exception:
            pass
    # Fall back to comma splitting
    return [x.strip() for x in s.split(',') if x.strip()]

def read_colors_map(path_or_str: str):
    """
    Read color mapping from a two-column headerless CSV:
    column 1: donor name
    column 2: HEX color
    If the argument is not a path, try literal_eval to parse a dict.
    """
    if path_or_str and os.path.exists(path_or_str):
        dfc = pd.read_csv(path_or_str, header=None)
        if dfc.shape[1] < 2:
            raise ValueError(f"colors CSV requires at least two columns without a header; current columns={dfc.shape[1]}")
        dfc = dfc.iloc[:, :2]
        dfc.columns = ['donor', 'color']
        # Remove duplicates while keeping the first occurrence
        dfc = dfc.dropna().drop_duplicates(subset=['donor'], keep='first')
        return dict(zip(dfc['donor'].astype(str), dfc['color'].astype(str)))
    # Fallback: if not a path, try parsing as a dictionary literal
    try:
        maybe = ast.literal_eval(path_or_str)
        if isinstance(maybe, dict):
            return {str(k): str(v) for k, v in maybe.items()}
    except Exception:
        pass
    raise FileNotFoundError("Colors CSV file not found, and the argument could not be parsed as a dictionary literal.")

def main():
    parser = argparse.ArgumentParser(
        description="Make donor-wise horizontal boxplots with global normalization."
    )
    parser.add_argument("--input", required=True, help="Absolute path to the input table (e.g., 1.txt)")
    parser.add_argument("--sep", default="\\t", help="Input delimiter, default '\\t'")
    parser.add_argument("--value-vars", required=True,
                        help="List of donor column names; supports JSON/py-literal or comma-separated strings")
    parser.add_argument("--colors-csv", required=True,
                        help="Color mapping CSV without header (donor, hex)")
    parser.add_argument("--valid-recipients", required=True,
                        help="List of valid recipients; supports JSON/py-literal or comma-separated strings")
    parser.add_argument("--n-cols", type=int, default=3, help="Number of subplots per row (default 3)")
    parser.add_argument("--out", default="Chromopainter_2.pdf", help="Output PDF path")
    args = parser.parse_args()

    RECIPIENT_COLUMN = "Recipient"
    ATTRIBUTE_COLUMN = "Attribute"
    ATTRIBUTE_COLUMN_LEGACY = "\u5c5e\u6027"
    DONOR_COLUMN = "Donor"
    VALUE_COLUMN = "Contribution"

    # Parse arguments
    value_vars = parse_list(args.value_vars)
    valid_recipients = parse_list(args.valid_recipients)
    if not value_vars:
        raise SystemExit("value_vars is empty. Provide at least one donor column name.")
    if not valid_recipients:
        raise SystemExit("valid_recipients is empty. Provide at least one recipient name.")

    colors = read_colors_map(args.colors_csv)

    # Load data
    df = pd.read_csv(args.input, sep=args.sep, dtype=str)
    if ATTRIBUTE_COLUMN not in df.columns and ATTRIBUTE_COLUMN_LEGACY in df.columns:
        df = df.rename(columns={ATTRIBUTE_COLUMN_LEGACY: ATTRIBUTE_COLUMN})

    # Convert numeric columns: only convert value_vars to numeric and keep non-numeric as NaN
    for c in value_vars:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    # Validate required columns
    required_cols = {RECIPIENT_COLUMN, ATTRIBUTE_COLUMN}
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise SystemExit(f"The input table is missing required columns: {missing}. Required columns: {sorted(required_cols)}")

    # Validate presence of value_vars
    missing_vals = [v for v in value_vars if v not in df.columns]
    if missing_vals:
        print(f"Warning: the following value_vars columns are missing and will be ignored: {missing_vals}", file=sys.stderr)
        value_vars = [v for v in value_vars if v in df.columns]
        if not value_vars:
            raise SystemExit("No valid value_vars columns were found in the data.")

    # Reshape into long format
    df_long = df.melt(
        id_vars=[RECIPIENT_COLUMN, ATTRIBUTE_COLUMN],
        value_vars=value_vars,
        var_name=DONOR_COLUMN,
        value_name=VALUE_COLUMN
    )

    # Global normalization (min/max computed on valid recipient subset only)
    global_data = df_long[df_long[ATTRIBUTE_COLUMN].isin(valid_recipients)].copy()

    # Abort if subset is empty because normalization would be invalid
    if global_data.empty:
        raise SystemExit("The subset of valid recipients is empty. Ensure valid_recipients matches values in column 'Attribute'.")

    global_min = pd.to_numeric(global_data[VALUE_COLUMN], errors="coerce").min()
    global_max = pd.to_numeric(global_data[VALUE_COLUMN], errors="coerce").max()

    # Handle the degenerate case where max == min to avoid division by zero
    denom = (global_max - global_min)
    if pd.isna(global_min) or pd.isna(global_max):
        raise SystemExit("All contribution values are NaN; normalization is impossible. Check the input data.")

    if denom == 0:
        # If all values are equal, assign a constant normalized value
        df_long["global_normalized"] = 0.5
    else:
        df_long["global_normalized"] = (pd.to_numeric(df_long[VALUE_COLUMN], errors="coerce") - global_min) / denom

    # Donors to plot: must exist in both the data and the color map
    donors_in_data = set(df_long[DONOR_COLUMN].dropna().unique().tolist())
    donor_list = [d for d in colors.keys() if d in donors_in_data]
    if not donor_list:
        raise SystemExit("No donors can be plotted. Verify that donor names in the color map match the data columns.")

    # Theme and global font configuration
    theme = _try_apply_aquarel_theme()
    plt.rcParams['font.sans-serif'] = ['Arial']
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    plt.rcParams['axes.unicode_minus'] = False

    # Compute subplot layout
    n_donors = len(donor_list)
    n_cols = max(1, int(args.n_cols))
    n_rows = math.ceil(n_donors / n_cols)

    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(4 * n_cols, 2.5 * n_rows),
        squeeze=False
    )

    # Plotting loop
    for i, donor in enumerate(donor_list):
        row = i // n_cols
        col = i % n_cols
        ax = axes[row][col]

        donor_data = df_long[df_long[DONOR_COLUMN] == donor].copy()
        donor_data = donor_data[donor_data[ATTRIBUTE_COLUMN].isin(valid_recipients)]

        if donor_data.empty:
            ax.set_visible(False)
            continue

        # Fallback to matplotlib boxplot if seaborn is unavailable
        if _HAS_SEABORN:
            import seaborn as sns
            sns.boxplot(
                data=donor_data,
                x="global_normalized",
                y=ATTRIBUTE_COLUMN,
                color=colors.get(donor, None),
                width=0.5,
                showfliers=True,
                flierprops=dict(
                    marker='o',
                    markersize=3,
                    markerfacecolor=colors.get(donor, None),
                    markeredgecolor=colors.get(donor, None),
                    alpha=0.5,
                ),
                boxprops={'edgecolor': colors.get(donor, None)},
                ax=ax
            )
        else:
            # Matplotlib fallback: draw per attribute
            ycats = list(donor_data[ATTRIBUTE_COLUMN].dropna().unique())
            y_index = {v: i+1 for i, v in enumerate(ycats)}
            positions = [y_index[v] for v in donor_data[ATTRIBUTE_COLUMN]]
            grouped = donor_data.groupby(ATTRIBUTE_COLUMN)["global_normalized"].apply(list)
            bp = ax.boxplot(
                grouped.values,
                vert=False,
                patch_artist=True,
                showfliers=True
            )
            # Apply colors
            for patch in bp['boxes']:
                patch.set_facecolor(colors.get(donor, "#BBBBBB"))
                patch.set_edgecolor(colors.get(donor, "#666666"))
            for fl in bp.get('fliers', []):
                fl.set_marker('o')
                fl.set_markersize(3)
                fl.set_alpha(0.5)
                fl.set_markerfacecolor(colors.get(donor, "#888888"))
                fl.set_markeredgecolor(colors.get(donor, "#888888"))
            ax.set_yticks(range(1, len(ycats)+1))
            ax.set_yticklabels(ycats)

        ax.set_title(f"Donor: {donor}")
        ax.grid(False)
        ax.set_xlabel("Global normalized contribution")
        ax.set_ylabel("Recipient")

    # Remove unused axes
    total_subplots = n_rows * n_cols
    if n_donors < total_subplots:
        for j in range(n_donors, total_subplots):
            row = j // n_cols
            col = j % n_cols
            fig.delaxes(axes[row][col])

    if theme is not None:
        try:
            theme.apply_transforms()
        except Exception:
            pass

    plt.tight_layout()
    out_path = args.out
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    plt.savefig(out_path)
    print(f"Saved figure to: {out_path}")
    # Use plt.show() manually if desired
    # plt.show()

if __name__ == "__main__":
    main()
