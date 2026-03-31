"""
Summarize binding metrics from concatenated lifetime events CSV.

Reads a CSV produced by join_csvs.py (columns: resid_i, resid_j, start_frame,
end_frame, frames, lifetime_ns, chainID_i, chainID_j, resname_i) and produces
a summary table with one row per unique resid_i.

Output columns:
  resid_i   – The target region (chain ID or resname, e.g. "4", "CLA")
  n_events  – Number of binding events for that resid_i
  median_ms – Median lifetime across events, in milliseconds
    iqr_ms    – Interquartile range (Q3-Q1) of lifetimes, in milliseconds
    mean_ms   – Mean lifetime across events, in milliseconds
    std_ms    – Standard deviation of lifetimes, in milliseconds
  sum_ms    – Total bound time (sum of lifetime_ns), in milliseconds    
    p_value   – Optional raw Mann-Whitney U p-value versus a reference resid_i
    p_value_bonferroni – Optional Bonferroni-corrected p-value versus a reference resid_i

Usage:
  python3 summarze_metrics_binding.py -c binding_pairs.csv -o binding_pairs_summary.csv
"""

import argparse
import pandas as pd
import sys
import yaml
from scipy.stats import mannwhitneyu


def parse_args():
    parser = argparse.ArgumentParser(
        description="Summarize binding metrics per resid_i from a lifetime events CSV."
    )
    parser.add_argument(
        "-c", "--csv", required=True,
        help="Input CSV file (concatenated lifetime events)."
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output summary CSV file."
    )
    parser.add_argument(
        "-total_sim_time_microseconds", type=float, default=None,
        help="Total simulation time in microseconds. If provided, a 'probability' "
             "column is added: sum_ms / (total_sim_time_ms)."
    )
    parser.add_argument(
        "-equivalent_chains_yaml", type=str, default=None,
        help="YAML file mapping group labels to lists of equivalent chain IDs. "
             "When provided, resid_i chain IDs are replaced by their group label "
             "(e.g. chains 1,2,5,6,g,n -> LHCBM) before aggregation."
    )
    parser.add_argument(
        "-reference_resid_i", type=str, default=None,
        help="Optional resid_i value to use as a reference group for per-row "
             "Mann-Whitney U tests on lifetime_ns. When provided, the output "
             "includes raw and Bonferroni-corrected p-values."
    )
    return parser.parse_args()


def main():
    args = parse_args()
    base_columns = ["resid_i", "n_events", "median_ms", "iqr_ms", "mean_ms", "std_ms", "sum_ms"]
    if args.reference_resid_i is not None:
        base_columns.extend(["p_value", "p_value_bonferroni"])

    # Read input, treating 'NA' as a literal string (it's a resname for
    # whole-chain entries where resname_i='NA')
    df = pd.read_csv(args.csv, keep_default_na=False, na_values=[''])

    if df.empty:
        print(f"Warning: Input CSV {args.csv} has no data rows. Writing empty summary.")
        summary = pd.DataFrame(columns=base_columns)
        summary.to_csv(args.output, index=False)
        print(f"Output saved to {args.output}")
        return

    required_cols = {"resid_i", "lifetime_ns"}
    missing = required_cols - set(df.columns)
    if missing:
        print(f"Error: Input CSV missing required columns: {missing}", file=sys.stderr)
        sys.exit(1)

    # Convert lifetime_ns to float (should already be, but be safe)
    df["lifetime_ns"] = pd.to_numeric(df["lifetime_ns"], errors="coerce")

    # Drop rows where lifetime_ns could not be parsed
    df = df.dropna(subset=["lifetime_ns"])

    if df.empty:
        print(f"Warning: No valid lifetime_ns values in {args.csv}. Writing empty summary.")
        summary = pd.DataFrame(columns=base_columns)
        summary.to_csv(args.output, index=False)
        print(f"Output saved to {args.output}")
        return

    # Convert ns -> ms
    ns_to_ms = 1e-3  # 1 ns = 1e-3 microseconds

    # If equivalent chains YAML is provided, map chain IDs to group labels
    if args.equivalent_chains_yaml is not None:
        with open(args.equivalent_chains_yaml, 'r') as f:
            equiv = yaml.safe_load(f)
        # Build chain_id -> label mapping
        chain_to_label = {}
        for label, chain_list in equiv.items():
            for chain_id in chain_list:
                chain_to_label[str(chain_id)] = label
        # Replace resid_i values that match a chain ID
        df["resid_i"] = df["resid_i"].astype(str).map(
            lambda x: chain_to_label.get(x, x)
        )
        print(f"Mapped chain IDs to group labels using {args.equivalent_chains_yaml}")

    # Group by resid_i and compute summary statistics
    grouped = df.groupby("resid_i", sort=False)["lifetime_ns"]

    summary = pd.DataFrame({
        "resid_i": grouped.first().index,
        "n_events": grouped.count().values,
        "median_ms": (grouped.median() * ns_to_ms).values,
        "iqr_ms": ((grouped.quantile(0.75) - grouped.quantile(0.25)) * ns_to_ms).values,
        "mean_ms": (grouped.mean() * ns_to_ms).values,
        "std_ms": (grouped.std(ddof=0) * ns_to_ms).values,
        "sum_ms": (grouped.sum() * ns_to_ms).values,
    })

    if args.reference_resid_i is not None:
        reference_resid_i = str(args.reference_resid_i)
        grouped_lifetimes = {
            resid_i: series.to_numpy()
            for resid_i, series in grouped
        }

        if reference_resid_i not in grouped_lifetimes:
            print(
                f"Error: reference resid_i '{reference_resid_i}' not found in input CSV.",
                file=sys.stderr,
            )
            sys.exit(1)

        reference_values = grouped_lifetimes[reference_resid_i]
        n_tests = max(len(grouped_lifetimes) - 1, 1)
        raw_p_values = []
        corrected_p_values = []

        for resid_i in summary["resid_i"]:
            resid_i = str(resid_i)
            if resid_i == reference_resid_i:
                p_value = 1.0
            else:
                _, p_value = mannwhitneyu(
                    grouped_lifetimes[resid_i],
                    reference_values,
                    alternative="two-sided",
                )
            raw_p_values.append(p_value)
            corrected_p_values.append(min(p_value * n_tests, 1.0))

        summary["p_value"] = raw_p_values
        summary["p_value_bonferroni"] = corrected_p_values

    # Compute probability if total simulation time is provided
    if args.total_sim_time_microseconds is not None:
        total_sim_time_ms = args.total_sim_time_microseconds 
        summary["probability"] = (summary["sum_ms"] / total_sim_time_ms) * 100 #as percentage

    # Sort by sum_ms descending (most bound first)
    summary = summary.sort_values("sum_ms", ascending=False).reset_index(drop=True)

    # Round all two two decimal places for cleaner output
    summary["median_ms"] = summary["median_ms"].round(2)
    summary["iqr_ms"] = summary["iqr_ms"].round(2)
    summary["mean_ms"] = summary["mean_ms"].round(2)
    summary["std_ms"] = summary["std_ms"].round(2)
    summary["sum_ms"] = summary["sum_ms"].round(2)
    if "p_value" in summary.columns:
        summary["p_value"] = summary["p_value"].round(6)
        summary["p_value_bonferroni"] = summary["p_value_bonferroni"].round(6)
    if "probability" in summary.columns:
        summary["probability"] = summary["probability"].round(2)

    summary.to_csv(args.output, index=False)
    print(f"Summary ({len(summary)} groups) saved to {args.output}")


if __name__ == "__main__":
    main()
