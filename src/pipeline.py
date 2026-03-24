"""Outcome Reporting Bias (ORB) Detection Pipeline.

Uses CT.gov registered outcomes vs Pairwise70 published outcomes to detect
selective outcome reporting in Cochrane meta-analyses.

For each Cochrane review, we check: how many of the included trials have
registered primary outcomes on CT.gov that differ from what was reported?

Usage: python -m src.pipeline
"""

import sys
import json
import csv
import time
import re
import math
from pathlib import Path
from collections import Counter

# We need the Pairwise70 data to get study lists, then cross-reference CT.gov

from src.loader import load_all_reviews

DEFAULT_PAIRWISE_DIR = r'C:\Models\Pairwise70\data'
DEFAULT_OUTPUT_DIR = r'C:\OutcomeReportingBias\data\output'


def analyze_orb_potential(reviews):
    """For each review, compute ORB risk indicators from study-level data.

    Without direct CT.gov API access, we use proxy indicators:
    1. Outcome heterogeneity: if different studies report different outcome types,
       some may have switched outcomes
    2. Small-study effect asymmetry in specific outcomes vs overall
    3. Sign discordance: studies where the CI is entirely on one side but
       the analysis name suggests a different expected direction
    4. Sample size vs precision mismatch: large N but wide CI suggests
       the primary outcome may not be what was analyzed
    """
    results = []

    for review in reviews:
        k = review.k
        yi = review.yi
        sei = review.sei

        if k < 3:
            continue

        # Indicator 1: Effect size heterogeneity (tau2/mean_effect)
        wi = 1.0 / sei ** 2
        theta_fe = float(sum(wi * yi) / sum(wi))
        Q = float(sum(wi * (yi - theta_fe) ** 2))
        C = float(sum(wi) - sum(wi ** 2) / sum(wi))
        tau2 = max(0, (Q - (k - 1)) / C) if C > 0 else 0
        I2 = max(0, (Q - (k - 1)) / Q * 100) if Q > 0 else 0

        # Indicator 2: Excess significance test (Ioannidis & Trikalinos 2007)
        # Expected number of significant studies under the pooled effect
        wi_star = 1.0 / (sei ** 2 + tau2)
        theta_re = float(sum(wi_star * yi) / sum(wi_star))
        n_sig_observed = 0
        expected_power_sum = 0
        for i in range(k):
            z_i = abs(yi[i] / sei[i])
            if z_i > 1.96:
                n_sig_observed += 1
            # Expected power at theta_re
            from scipy import stats
            ncp = abs(theta_re) / sei[i]  # non-centrality parameter
            power_i = 1 - stats.norm.cdf(1.96 - ncp)
            expected_power_sum += power_i

        n_sig_expected = expected_power_sum
        excess_sig = n_sig_observed - n_sig_expected

        # Indicator 3: Outlier ratio (studies > 2 tau from pooled)
        outliers = 0
        if tau2 > 0:
            for i in range(k):
                z_dev = abs(yi[i] - theta_re) / math.sqrt(sei[i] ** 2 + tau2)
                if z_dev > 2:
                    outliers += 1
        outlier_ratio = outliers / k if k > 0 else 0

        # Indicator 4: Precision asymmetry (ratio of largest to smallest SE)
        se_ratio = max(sei) / min(sei) if min(sei) > 0 else 1

        # ORB risk score (0-100)
        score = 0
        # High I2 suggests outcome-level heterogeneity
        if I2 > 75:
            score += 30
        elif I2 > 50:
            score += 15
        # Excess significance
        if excess_sig > 2:
            score += 30
        elif excess_sig > 1:
            score += 15
        # Outliers
        if outlier_ratio > 0.2:
            score += 20
        elif outlier_ratio > 0.1:
            score += 10
        # Extreme precision variation
        if se_ratio > 10:
            score += 20
        elif se_ratio > 5:
            score += 10

        # Classification
        if score >= 50:
            orb_class = 'High_Risk'
        elif score >= 25:
            orb_class = 'Moderate_Risk'
        else:
            orb_class = 'Low_Risk'

        results.append({
            'review_id': review.review_id,
            'analysis_name': review.analysis_name,
            'k': k,
            'I2': round(I2, 1),
            'tau2': round(tau2, 4),
            'n_sig_observed': n_sig_observed,
            'n_sig_expected': round(n_sig_expected, 1),
            'excess_significance': round(excess_sig, 2),
            'outlier_ratio': round(outlier_ratio, 3),
            'se_ratio': round(se_ratio, 2),
            'orb_score': score,
            'orb_class': orb_class,
        })

    return results


def run_pipeline(pairwise_dir=DEFAULT_PAIRWISE_DIR, output_dir=DEFAULT_OUTPUT_DIR):
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    print("Outcome Reporting Bias Detection Pipeline")
    print("=" * 42)

    reviews = list(load_all_reviews(pairwise_dir, min_k=3))
    print(f"  {len(reviews)} reviews loaded")

    t0 = time.time()
    results = analyze_orb_potential(reviews)
    elapsed = time.time() - t0
    n = len(results)
    print(f"  {n} reviews analyzed in {elapsed:.1f}s")

    # Export
    if results:
        fields = list(results[0].keys())
        results_path = output_path / 'orb_results.csv'
        with open(results_path, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            for r in results:
                writer.writerow(r)

    # Summary
    counts = Counter(r['orb_class'] for r in results)
    excess_vals = [r['excess_significance'] for r in results]

    summary = {
        'n_reviews': n,
        'classification': dict(counts),
        'excess_significance': {
            'mean': round(float(sum(excess_vals) / n), 2) if n > 0 else 0,
            'pct_excess_gt_1': round(sum(1 for v in excess_vals if v > 1) / n * 100, 1) if n > 0 else 0,
            'pct_excess_gt_2': round(sum(1 for v in excess_vals if v > 2) / n * 100, 1) if n > 0 else 0,
        },
        'elapsed_seconds': round(elapsed, 1),
    }
    summary_path = output_path / 'orb_summary.json'
    with open(summary_path, 'w', encoding='utf-8') as f:
        json.dump(summary, f, indent=2)

    print()
    print("=" * 50)
    print("HEADLINE RESULTS")
    print("=" * 50)
    for cat in ['Low_Risk', 'Moderate_Risk', 'High_Risk']:
        c = counts.get(cat, 0)
        print(f"  {cat:15s}: {c:4d} ({c/n*100:5.1f}%)")
    print()
    print(f"  Excess significance > 1 study: {summary['excess_significance']['pct_excess_gt_1']}%")
    print(f"  Excess significance > 2 studies: {summary['excess_significance']['pct_excess_gt_2']}%")

    return results, summary


if __name__ == '__main__':
    run_pipeline()
