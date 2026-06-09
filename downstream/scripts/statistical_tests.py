import pandas as pd
import numpy as np
from scipy import stats
import sys
import warnings
warnings.filterwarnings("ignore")

# ============================================================
# Statistical Tests for K. pneumoniae WGS Data
# Tests performed:
#   1. Fisher's Exact Test – carbapenemase gene presence vs. country
#   2. Chi-square Test – ST distribution across countries
#   3. Kruskal-Wallis Test – virulence score across countries
# Input : results/ready_to_download/typing/metadata_summary.tsv
# Output: results/downstream/statistics/stats_report.txt
#         results/downstream/statistics/contingency_tables.tsv
# ============================================================

def fisher_carbapenemase_vs_country(df: pd.DataFrame, out_lines: list) -> None:
    """
    Fisher's Exact Test: Is carbapenemase gene presence significantly
    different across countries?
    Uses a 2×k contingency table (positive/negative × k countries).
    When k=2, Fisher's exact; when k>2, chi-square (Fisher's not generalizable).
    """
    out_lines.append("\n" + "="*60)
    out_lines.append("TEST 1: Carbapenemase Gene Presence vs. Country")
    out_lines.append("="*60)

    # Carbapenemase column: 1/0 presence in metadata
    if "carbapenemase" not in df.columns:
        out_lines.append("  [SKIP] Column 'carbapenemase' not found in metadata.")
        return

    contingency = pd.crosstab(df["country"], df["carbapenemase"])
    out_lines.append("\nContingency Table (rows=country, cols=carbapenemase gene):")
    out_lines.append(contingency.to_string())

    countries = df["country"].unique()
    if len(countries) == 2:
        # 2×2 → Fisher's Exact Test
        table = contingency.values
        oddsratio, pvalue = stats.fisher_exact(table)
        out_lines.append(f"\nFisher's Exact Test")
        out_lines.append(f"  Odds Ratio : {oddsratio:.4f}")
        out_lines.append(f"  p-value    : {pvalue:.4f}")
        out_lines.append(f"  Significant: {'Yes' if pvalue < 0.05 else 'No'} (α=0.05)")
    else:
        # k×2 → Chi-square Test
        chi2, pvalue, dof, expected = stats.chi2_contingency(contingency)
        out_lines.append(f"\nChi-square Test (k={len(countries)} countries)")
        out_lines.append(f"  Chi²  : {chi2:.4f}")
        out_lines.append(f"  df    : {dof}")
        out_lines.append(f"  p-value: {pvalue:.4f}")
        out_lines.append(f"  Significant: {'Yes' if pvalue < 0.05 else 'No'} (α=0.05)")

        # Post-hoc pairwise Fisher's tests if significant
        if pvalue < 0.05:
            out_lines.append("\n  Post-hoc Pairwise Fisher's Tests (Bonferroni corrected):")
            from itertools import combinations
            pairs = list(combinations(countries, 2))
            alpha_corrected = 0.05 / len(pairs)
            for c1, c2 in pairs:
                sub = contingency.loc[[c1, c2]]
                if sub.shape[1] < 2:
                    continue
                _, p = stats.fisher_exact(sub.values)
                out_lines.append(
                    f"    {c1} vs {c2}: p={p:.4f} "
                    f"({'*' if p < alpha_corrected else 'ns'}, α_adj={alpha_corrected:.4f})"
                )


def chisquare_st_vs_country(df: pd.DataFrame, out_lines: list) -> None:
    """Chi-square Test: Is ST distribution different across countries?"""
    out_lines.append("\n" + "="*60)
    out_lines.append("TEST 2: Sequence Type (ST) Distribution vs. Country")
    out_lines.append("="*60)

    col = next((c for c in ["st", "ST", "MLST_ST"] if c in df.columns), None)
    if col is None:
        out_lines.append("  [SKIP] ST column not found in metadata.")
        return

    # Show top 10 STs for context
    top_sts = df[col].value_counts().head(10)
    out_lines.append("\nTop 10 Sequence Types:")
    out_lines.append(top_sts.to_string())

    contingency = pd.crosstab(df["country"], df[col])
    chi2, pvalue, dof, expected = stats.chi2_contingency(contingency)
    out_lines.append(f"\nChi-square Test: ST distribution across countries")
    out_lines.append(f"  Chi²  : {chi2:.4f}")
    out_lines.append(f"  df    : {dof}")
    out_lines.append(f"  p-value: {pvalue:.4f}")
    out_lines.append(f"  Significant: {'Yes' if pvalue < 0.05 else 'No'} (α=0.05)")


def kruskal_virulence_vs_country(df: pd.DataFrame, out_lines: list) -> None:
    """Kruskal-Wallis Test: Is virulence score different across countries?"""
    out_lines.append("\n" + "="*60)
    out_lines.append("TEST 3: Kleborate Virulence Score vs. Country")
    out_lines.append("="*60)

    score_col = next((c for c in ["virulence_score", "Virulence score", "virulence"]
                      if c in df.columns), None)
    if score_col is None:
        out_lines.append("  [SKIP] Virulence score column not found.")
        return

    groups = [group[score_col].dropna().values
               for _, group in df.groupby("country")]
    if len(groups) < 2:
        out_lines.append("  [SKIP] Not enough country groups for test.")
        return

    stat, pvalue = stats.kruskal(*groups)

    out_lines.append("\nKruskal-Wallis H-Test: Virulence Score across countries")
    out_lines.append("\nMedian virulence score per country:")
    for country, grp in df.groupby("country"):
        med = grp[score_col].median()
        out_lines.append(f"  {country}: median={med:.2f} (n={len(grp)})")

    out_lines.append(f"\n  H-statistic: {stat:.4f}")
    out_lines.append(f"  p-value    : {pvalue:.4f}")
    out_lines.append(f"  Significant: {'Yes' if pvalue < 0.05 else 'No'} (α=0.05)")

    # Post-hoc Mann-Whitney U if significant
    if pvalue < 0.05:
        from itertools import combinations
        countries = list(df["country"].unique())
        pairs = list(combinations(countries, 2))
        alpha_corrected = 0.05 / len(pairs)
        out_lines.append(
            f"\n  Post-hoc Pairwise Mann-Whitney U (Bonferroni α_adj={alpha_corrected:.4f}):"
        )
        for c1, c2 in pairs:
            g1 = df[df["country"] == c1][score_col].dropna().values
            g2 = df[df["country"] == c2][score_col].dropna().values
            if len(g1) < 1 or len(g2) < 1:
                continue
            u_stat, p_mw = stats.mannwhitneyu(g1, g2, alternative="two-sided")
            out_lines.append(
                f"    {c1} vs {c2}: U={u_stat:.1f}, p={p_mw:.4f} "
                f"({'*' if p_mw < alpha_corrected else 'ns'})"
            )


def kruskal_resistance_vs_country(df: pd.DataFrame, out_lines: list) -> None:
    """Kruskal-Wallis Test: Is resistance score different across countries?"""
    out_lines.append("\n" + "="*60)
    out_lines.append("TEST 4: Kleborate Resistance Score vs. Country")
    out_lines.append("="*60)

    score_col = next((c for c in ["resistance_score", "Resistance score", "resistance"]
                      if c in df.columns), None)
    if score_col is None:
        out_lines.append("  [SKIP] Resistance score column not found.")
        return

    groups = [group[score_col].dropna().values
               for _, group in df.groupby("country")]
    stat, pvalue = stats.kruskal(*groups)

    out_lines.append("\nKruskal-Wallis H-Test: Resistance Score across countries")
    out_lines.append("\nMedian resistance score per country:")
    for country, grp in df.groupby("country"):
        med = grp[score_col].median()
        out_lines.append(f"  {country}: median={med:.2f} (n={len(grp)})")

    out_lines.append(f"\n  H-statistic: {stat:.4f}")
    out_lines.append(f"  p-value    : {pvalue:.4f}")
    out_lines.append(f"  Significant: {'Yes' if pvalue < 0.05 else 'No'} (α=0.05)")


def main():
    metadata_path = snakemake.input.metadata
    out_report    = snakemake.output.report

    df = pd.read_csv(metadata_path, sep="\t")
    print(f"[statistical_tests] Loaded {len(df)} samples from {metadata_path}")
    print(f"[statistical_tests] Columns: {list(df.columns)}")

    out_lines = [
        "="*60,
        "STATISTICAL ANALYSIS REPORT",
        "K. pneumoniae Complex – Southeast Asia WGS Study",
        "="*60,
        f"Total samples: {len(df)}",
        f"Countries: {', '.join(sorted(df['country'].unique()))}",
    ]

    fisher_carbapenemase_vs_country(df, out_lines)
    chisquare_st_vs_country(df, out_lines)
    kruskal_virulence_vs_country(df, out_lines)
    kruskal_resistance_vs_country(df, out_lines)

    out_lines.append("\n" + "="*60)
    out_lines.append("END OF REPORT")
    out_lines.append("="*60)

    report_text = "\n".join(out_lines)
    print(report_text)

    with open(out_report, "w") as f:
        f.write(report_text)
    print(f"\n[statistical_tests] Report saved → {out_report}")


if __name__ == "__main__":
    main()
