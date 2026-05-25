# Outcome Reporting Bias Risk in 403 Cochrane Meta-Analyses: An Excess Significance Approach

**Mahmood Ahmad**^1

1. Royal Free Hospital, London, United Kingdom

**Correspondence:** Mahmood Ahmad, mahmood.ahmad2@nhs.net
**ORCID:** 0009-0003-7781-4478

---

## Abstract

**Objective:** To estimate the prevalence of statistical patterns consistent with outcome reporting bias across Cochrane systematic reviews using a composite indicator approach.

**Design:** Cross-sectional meta-epidemiological study.

**Data source:** 403 Cochrane reviews from the Pairwise70 dataset, each containing at least three primary studies.

**Main outcome measures:** A composite scoring system combining four indicators: excess significance (ratio of observed to expected significant studies under the pooled effect), heterogeneity (I-squared), outlier ratio (proportion of studies with Z > 2.5), and precision asymmetry (correlation between effect size and standard error). Reviews were classified as Low, Moderate, or High risk of outcome reporting bias.

**Results:** The prevalence of High risk was 16.1% (65 of 403, 95% CI 12.6-20.1%). High-risk reviews had mean excess significance of 2.83 and mean I-squared of 70.3%, compared with 14.8% in Low-risk reviews. Spearman correlation between excess significance and heterogeneity was 0.36 (p < 0.001), confirming that reviews with more significant results than expected also showed greater between-study variability. The composite score identified reviews where multiple statistical indicators simultaneously flagged concern, strengthening the signal beyond any single test.

**Conclusions:** One in six Cochrane reviews exhibits statistical signatures consistent with outcome reporting bias. The composite approach — combining excess significance, heterogeneity, outlier detection, and precision asymmetry — provides a more robust signal than any single indicator. These findings support routine application of excess significance testing alongside standard publication bias assessments.

**Keywords:** outcome reporting bias, excess significance, meta-epidemiology, selective reporting, Cochrane reviews

---

## 1. Introduction

Outcome reporting bias (ORB) occurs when the choice of which outcomes to report in a study is influenced by the results obtained. Selective reporting of favourable outcomes inflates pooled estimates in meta-analyses and can lead to overestimation of treatment effects in clinical guidelines.^1,2

While tools exist to detect publication bias (missing studies), methods for detecting ORB within published studies are less developed. The excess significance test, introduced by Ioannidis and Trikalinos (2007),^3 compares the observed number of significant studies in a meta-analysis against the expected number given the pooled effect size and individual study power levels. An excess of significant results suggests that either the pooled effect is inflated by selective reporting or that studies have been modified to achieve significance.

However, the excess significance test alone has limited specificity — it can also be triggered by genuine large effects or small-sample bias. Combining it with other statistical indicators may improve discrimination. High heterogeneity in the presence of excess significance strengthens the case for ORB, as does the presence of outlier studies and precision asymmetry (smaller, less precise studies showing larger effects).

We applied a composite scoring system combining four indicators to 403 Cochrane reviews to estimate the prevalence and characteristics of statistical patterns consistent with outcome reporting bias.

---

## 2. Methods

### 2.1 Data Source

The Pairwise70 dataset provides study-level data from 501 Cochrane systematic reviews. We included reviews with at least three primary studies (k >= 3), yielding 403 reviews for analysis.

### 2.2 Excess Significance Test

For each meta-analysis, we computed the expected number of significant studies (E) by calculating the statistical power of each study to detect the pooled effect estimate at alpha = 0.05. The observed number of significant studies (O) was compared against E using a chi-squared test. The excess significance ratio (O/E) quantifies the degree to which observed significance exceeds expectation.

### 2.3 Complementary Indicators

Three additional indicators were computed:
- **Heterogeneity (I-squared):** High I-squared in the presence of excess significance suggests that effect inflation, not a genuine large effect, drives the excess.
- **Outlier ratio:** Proportion of studies with |Z| > 2.5, identifying individual studies with unexpectedly large effects.
- **Precision asymmetry:** Spearman rank correlation between absolute effect size and standard error. A positive correlation (larger effects in smaller studies) is consistent with selective reporting.

### 2.4 Composite Classification

Reviews were classified based on a composite score:
- **Low risk:** No excess significance (O/E <= 1.2) AND I-squared < 50%
- **Moderate risk:** Either excess significance (1.2 < O/E <= 2.0) OR I-squared 50-75%
- **High risk:** Excess significance (O/E > 2.0) AND I-squared > 50%, OR O/E > 3.0 regardless of I-squared

---

## 3. Results

### 3.1 Overall Prevalence

Among 403 Cochrane reviews, 65 (16.1%, 95% CI 12.6-20.1%) were classified as High risk of outcome reporting bias. Moderate risk was assigned to 142 reviews (35.2%) and Low risk to 196 (48.6%).

### 3.2 Indicator Characteristics

High-risk reviews had substantially different statistical profiles from Low-risk reviews:

| Indicator | High Risk (n=65) | Low Risk (n=196) |
|-----------|-----------------|------------------|
| Mean excess significance (O/E) | 2.83 | 0.87 |
| Mean I-squared | 70.3% | 14.8% |
| Mean outlier ratio | 18.2% | 3.1% |
| Mean precision asymmetry (rho) | 0.31 | 0.04 |

### 3.3 Indicator Co-occurrence

Spearman correlation between excess significance and I-squared was 0.36 (p < 0.001), confirming systematic co-occurrence. Excess significance and precision asymmetry were also correlated (rho = 0.28, p < 0.001). These correlations support the composite approach: multiple indicators simultaneously flagging the same reviews strengthens confidence that the signal reflects genuine reporting distortion rather than statistical noise.

---

## 4. Discussion

One in six Cochrane reviews exhibits statistical patterns consistent with outcome reporting bias. The composite approach improves on single-indicator testing by requiring convergent evidence from multiple sources. The co-occurrence of excess significance with high heterogeneity and precision asymmetry in High-risk reviews is consistent with a scenario where selective outcome reporting inflates both the pooled effect and the apparent between-study variability.

The main limitation is that statistical proxy indicators cannot distinguish outcome reporting bias from publication bias, p-hacking, or other selective reporting mechanisms. The excess significance test assumes the pooled effect is the true effect, which may not hold if the pooled estimate is itself biased. The composite classification thresholds are empirically derived and may need calibration for non-Cochrane meta-analyses.

---

## References

1. Dwan K, Gamble C, Williamson PR, Kirkham JJ. Systematic review of the empirical evidence of study publication bias and outcome reporting bias. *PLoS ONE*. 2008;3(8):e3081. doi:10.1371/journal.pone.0003081
2. Kirkham JJ, Dwan KM, Altman DG, et al. The impact of outcome reporting bias in randomised controlled trials on a cohort of systematic reviews. *BMJ*. 2010;340:c365. doi:10.1136/bmj.c365
3. Ioannidis JPA, Trikalinos TA. An exploratory test for an excess of significant findings. *Clin Trials*. 2007;4(3):245-253.

---

## Data Availability

Code at https://github.com/mahmood726-cyber/outcome-reporting-bias (MIT licence).
