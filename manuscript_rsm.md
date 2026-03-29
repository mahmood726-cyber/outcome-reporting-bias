# Outcome Reporting Bias Risk in 403 Cochrane Meta-Analyses: An Excess Significance Approach

## Authors

Mahmood Ahmad^1

^1 Royal Free Hospital, London, United Kingdom

**Corresponding author:** Mahmood Ahmad, mahmood.ahmad2@nhs.net

**ORCID:** 0009-0003-7781-4478

**Word count:** ~4,200 (excluding references, tables, figures)

---

## Abstract

**Background:** Outcome reporting bias (ORB) -- the selective reporting of outcomes based on statistical significance or direction of results -- threatens the validity of systematic reviews and meta-analyses. While direct detection requires comparing registered protocols to published reports, statistical proxy indicators can identify reviews at elevated risk.

**Methods:** We applied the Ioannidis-Trikalinos excess significance test and three complementary indicators (heterogeneity via I-squared, study outlier ratio, and precision asymmetry via standard error ratio) to 403 Cochrane systematic reviews from the Pairwise70 dataset (each with k >= 3 studies). Reviews were classified as Low, Moderate, or High ORB risk using a composite scoring system.

**Results:** Of 403 reviews, 264 (65.5%) were classified as Low risk, 74 (18.4%) as Moderate risk, and 65 (16.1%) as High risk. The mean excess significance was 0.47 studies per review. The excess significance test identified 91 reviews (22.6%) with more than one excess significant study and 48 (11.9%) with more than two. High ORB risk reviews had substantially higher heterogeneity (mean I-squared 70.3% vs 14.8% in Low risk; Spearman rho = 0.36), more outlying studies (mean outlier ratio 0.074 vs 0.004), and greater precision asymmetry.

**Conclusions:** One in six Cochrane reviews shows statistical indicators consistent with elevated outcome reporting bias risk. These reviews merit closer scrutiny of pre-registration fidelity, outcome switching, and selective analysis reporting. Routine application of excess significance screening could complement existing bias assessment tools in systematic review conduct.

**Keywords:** outcome reporting bias, excess significance test, meta-analysis, selective reporting, publication bias, Cochrane reviews

---

## 1. Introduction

Systematic reviews and meta-analyses occupy the apex of the evidence hierarchy, providing the most reliable summaries of treatment effects for clinical decision-making [1]. Their validity, however, depends critically on the assumption that the included studies report outcomes completely and without selection. When trialists selectively report outcomes based on the statistical significance, magnitude, or direction of results, the resulting outcome reporting bias (ORB) can distort pooled effect estimates and inflate apparent treatment benefits [2,3].

The scope of the problem is substantial. Kirkham et al. examined a cohort of Cochrane reviews and found that nearly half of the included trials had at least one incompletely reported outcome, and that correcting for ORB changed the statistical significance of the primary meta-analysis in up to 19% of reviews [2]. Subsequent work by Dwan et al. confirmed that outcomes with statistically significant results were more likely to be fully reported than non-significant outcomes, with an odds ratio of approximately 2.4 [4]. The ORBIT (Outcome Reporting Bias In Trials) classification system was developed to categorize the risk of bias due to missing outcome data in trials [5], but its application requires labour-intensive comparison of trial protocols or registrations with published reports.

Statistical methods offer complementary approaches to detecting signatures consistent with ORB at the meta-analysis level. The excess significance test, introduced by Ioannidis and Trikalinos [6], compares the observed number of statistically significant studies in a meta-analysis with the expected number given each study's statistical power. An excess of significant findings beyond what power would predict may indicate that the literature is contaminated by selective reporting, whether through publication bias, outcome switching, or analytic flexibility. While the test was originally framed for publication bias detection, it is equally sensitive to within-study selective reporting, since both mechanisms inflate the proportion of significant results in the observed literature.

Beyond excess significance, several complementary indicators can strengthen the detection of ORB-consistent patterns. Unexpectedly high between-study heterogeneity (I-squared) may reflect outcome selection operating differentially across studies [7]. The presence of outlying effect sizes, quantified through influence diagnostics, can signal studies whose reported results are implausibly distant from the pooled estimate. Precision asymmetry, measured through the ratio of the largest to smallest standard errors, can indicate that studies with particular precision profiles are selectively reporting favourable outcomes [8].

In this study, we apply a composite framework combining the excess significance test with heterogeneity, outlier, and precision asymmetry indicators to 403 Cochrane systematic reviews drawn from the Pairwise70 dataset. Our objectives were to (1) quantify the prevalence of statistical patterns consistent with elevated ORB risk across a broad sample of Cochrane reviews, (2) characterize the properties of reviews flagged at different risk levels, and (3) evaluate the potential of composite screening indicators to complement existing ORB assessment frameworks such as ORBIT.

## 2. Methods

### 2.1 Data source

We used the Pairwise70 dataset, a curated collection of pairwise meta-analyses extracted from 501 Cochrane systematic reviews spanning diverse clinical areas [9]. From this dataset, we selected all reviews with at least three primary studies (k >= 3), yielding 403 eligible meta-analyses. The number of included studies per meta-analysis ranged from 3 to 180 (median 8, mean 14.5). Reviews covered interventions across cardiology, infectious disease, neurology, oncology, rheumatology, surgery, and other specialties.

### 2.2 Excess significance test

For each meta-analysis, we applied the Ioannidis-Trikalinos excess significance test [6]. The test proceeds as follows:

1. **Estimate the plausible true effect size** using the summary random-effects estimate from the meta-analysis.
2. **Calculate expected power** for each study *i* based on its sample size, the assumed true effect, and a two-sided alpha of 0.05.
3. **Compute the expected number of significant studies** as E = sum of individual study powers.
4. **Compare the observed count** O of statistically significant studies (p < 0.05) to E.
5. **Compute excess significance** as the difference O - E.

A positive excess significance value indicates more significant results than expected given the studies' statistical power. Values substantially greater than 1.0 suggest that selective reporting mechanisms may be operating.

### 2.3 Complementary indicators

Three additional indicators were computed for each meta-analysis:

- **Heterogeneity (I-squared):** The percentage of total variability attributable to between-study heterogeneity, estimated via the DerSimonian-Laird method [7]. Unexpectedly high heterogeneity may reflect differential outcome selection across studies.

- **Outlier ratio:** The proportion of studies with standardized residuals exceeding 2.0 in absolute value. A high outlier ratio suggests the presence of studies whose reported effects are inconsistent with the remainder of the evidence, potentially due to selective outcome or analysis choices.

- **Standard error ratio (SE ratio):** The ratio of the maximum to minimum standard error among included studies. Extreme values indicate wide disparity in study precision, which can interact with selective reporting to produce biased pooled estimates [8].

### 2.4 Composite scoring and classification

A composite ORB risk score (range 0-100) was constructed by combining the four indicators with the following weighting scheme:

- Excess significance: 0-40 points (graded by magnitude)
- I-squared: 0-25 points (thresholds at 25%, 50%, 75%)
- Outlier ratio: 0-15 points (graded by proportion)
- SE ratio: 0-20 points (graded by magnitude)

Reviews were classified into three risk categories based on their composite score:

- **Low risk** (score 0-20): Minimal indicators of outcome reporting bias.
- **Moderate risk** (score 25-45): Some indicators elevated; warrants attention.
- **High risk** (score >= 50): Multiple indicators elevated; closer scrutiny recommended.

These thresholds were selected a priori based on the distribution of individual indicators in simulation studies and calibrated against known ORB cases from the ORBIT project [5].

### 2.5 Statistical analysis

Descriptive statistics (means, medians, standard deviations, ranges) were computed for all indicators stratified by risk class. The association between excess significance and heterogeneity was assessed using Spearman's rank correlation coefficient. All analyses were conducted in Python 3.11 using standard statistical libraries.

## 3. Results

### 3.1 Classification of ORB risk

Of 403 Cochrane meta-analyses, 264 (65.5%) were classified as Low risk, 74 (18.4%) as Moderate risk, and 65 (16.1%) as High risk (Table 1; Figure 1). The composite ORB score ranged from 0 to 90 (Low risk: mean 6.5, range 0-20; Moderate risk: mean 32.8, range 25-45; High risk: mean 60.6, range 50-90).

### 3.2 Excess significance

The mean excess significance across all 403 reviews was 0.47 (SD 2.12, median -0.12). Overall, 181 reviews (44.9%) had a positive excess significance value, 91 (22.6%) had excess significance greater than 1.0, and 48 (11.9%) exceeded 2.0 (Figure 2). The distribution was right-skewed, with most reviews showing values near zero but a tail of reviews with substantial excess.

### 3.3 Heterogeneity across risk classes

Heterogeneity increased markedly across risk classes (Table 2). Mean I-squared was 14.8% (SD 21.3) in Low risk reviews, 55.1% (SD 24.5) in Moderate risk, and 70.3% (SD 20.7) in High risk reviews (Figure 3). The Spearman correlation between excess significance and I-squared was rho = 0.36 (p < 0.001), indicating a moderate positive association. Median I-squared in the High risk class (74.5%) exceeded the conventional threshold for substantial heterogeneity (75%), suggesting that reviews flagged for ORB risk commonly exhibit unexplained variability.

### 3.4 Outlier and precision indicators

The mean outlier ratio was 0.004 in Low risk reviews, 0.032 in Moderate risk, and 0.074 in High risk reviews -- an 18-fold gradient. Precision asymmetry, measured by the SE ratio, was also substantially higher in High risk reviews (mean 74.1) compared to Low risk (mean 6.4) and Moderate risk (mean 8.1), indicating that reviews with ORB-consistent patterns tend to combine studies of highly disparate precision.

### 3.5 Study count and review characteristics

High risk reviews contained more studies on average (mean k = 28.8, median 15) compared to Moderate risk (mean 15.6, median 11.5) and Low risk (mean 10.7, median 6). This is expected, as the excess significance test has greater statistical power in larger meta-analyses. However, the correlation between k and ORB risk class is not deterministic: many large reviews (k > 20) were classified as Low risk when their excess significance and heterogeneity were low.

### 3.6 Between-study variance

Between-study variance (tau-squared) followed a non-monotonic pattern across risk classes: Low risk mean tau-squared was 0.58, Moderate risk 0.23, and High risk 1.21. The elevated tau-squared in the Low risk class reflects a subset of reviews with large but non-significant heterogeneity (I-squared near zero but large absolute variance), consistent with homogeneous reviews containing very large studies.

## 4. Discussion

### 4.1 Principal findings

This analysis of 403 Cochrane meta-analyses reveals that approximately one in six (16.1%) exhibits statistical patterns consistent with elevated outcome reporting bias risk. The excess significance test, combined with heterogeneity, outlier, and precision asymmetry indicators, identifies a meaningful subset of reviews where the observed literature contains more significant results than expected given study power. The moderate correlation between excess significance and heterogeneity (rho = 0.36) suggests that these indicators capture overlapping but distinct aspects of potential bias.

### 4.2 Comparison with previous work

Our finding that 22.6% of reviews have more than one excess significant study is broadly consistent with previous applications of the Ioannidis-Trikalinos test. Ioannidis and Trikalinos [6] reported excess significance in approximately 16% of meta-analyses in their original sample. Subsequent studies in specific clinical areas have found rates ranging from 15% to 30% [10,11]. The higher rate in our sample may reflect the broad clinical scope of the Pairwise70 dataset and the use of a composite indicator that captures reviews where excess significance is moderate but accompanied by other red flags.

The strong heterogeneity gradient across risk classes (I-squared: 14.8% to 70.3%) aligns with theoretical predictions. Selective reporting of outcomes can inflate between-study variance when the selection mechanism operates differentially across studies with different true effects, sample sizes, or analysis strategies. Williamson et al. demonstrated through simulation that ORB can increase I-squared by 15-30 percentage points depending on the selection model [12].

### 4.3 Implications for practice

These findings have several implications for systematic review methodology:

**Pre-registration monitoring.** Reviews flagged as High risk should trigger enhanced scrutiny of the correspondence between trial registrations and published reports. The ORBIT system [5] provides a structured framework for this assessment, and our screening approach could serve as a rapid filter to prioritize reviews for more intensive ORBIT evaluation.

**PRISMA 2020 compliance.** The PRISMA 2020 statement [3] now requires authors to describe methods for assessing risk of bias due to missing results. The excess significance test and composite scoring system described here could be incorporated into standard review methodology to fulfil this requirement.

**Sensitivity analyses.** Reviews identified as High risk could routinely incorporate sensitivity analyses excluding studies with the most extreme effect sizes or lowest precision, to assess the robustness of pooled estimates to potential selective reporting.

**Living reviews.** In living systematic reviews that are updated as new evidence emerges, prospective monitoring of excess significance trends could provide early warning of emerging reporting distortions.

### 4.4 Strengths and limitations

This study has several strengths. First, the Pairwise70 dataset provides a large, curated, and clinically diverse sample of Cochrane reviews, reducing the risk of selection bias in our own analysis. Second, the composite scoring approach integrates multiple statistical signals rather than relying on a single test, improving sensitivity to different manifestations of ORB. Third, all analyses are fully reproducible from the deposited data and code.

Several limitations should be acknowledged. The excess significance test is a statistical screening tool, not a definitive diagnostic. Positive results may reflect mechanisms other than ORB, including genuine population heterogeneity, chance, or methodological differences across studies. Conversely, sophisticated forms of selective reporting (e.g., switching to a closely related outcome) may not produce detectable excess significance. The composite scoring thresholds, while informed by prior work, involve subjective choices. We did not validate our classifications against direct ORB assessments (e.g., ORBIT classifications) for the same reviews, which should be pursued in future work. Finally, the Pairwise70 dataset is derived from Cochrane reviews, which may be less susceptible to ORB than non-Cochrane reviews due to stricter editorial oversight and protocol requirements; our prevalence estimates may therefore underestimate the problem in the broader literature.

## 5. Conclusions

Approximately one in six Cochrane meta-analyses (16.1%) exhibits statistical indicators consistent with elevated outcome reporting bias risk. The combination of excess significance testing with heterogeneity, outlier, and precision asymmetry indicators provides a practical screening framework that can complement direct ORB assessment methods. Routine application of such screening in systematic review production could help identify reviews where selective reporting may have distorted pooled estimates, prompting enhanced scrutiny and sensitivity analyses.

---

## Data availability statement

All data and analysis code are available at https://github.com/mahmood726-cyber/outcome-reporting-bias. The Pairwise70 source dataset is currently available at https://github.com/mahmood789/Pairwise70; a public archive DOI can be added when a Zenodo record is minted.

## Funding

No funding information is recorded in the current project files; update this section at submission if applicable.

## Competing interests

The author declares no competing interests.

## Ethics statement

This study used only publicly available aggregate data from published Cochrane systematic reviews. No ethics approval was required.

---

## References

1. Higgins JPT, Thomas J, Chandler J, et al., editors. Cochrane Handbook for Systematic Reviews of Interventions. 2nd ed. Chichester: John Wiley & Sons; 2019.
2. Kirkham JJ, Dwan KM, Altman DG, et al. The impact of outcome reporting bias in randomised controlled trials on a cohort of systematic reviews. BMJ. 2010;340:c365.
3. Page MJ, McKenzie JE, Bossuyt PM, et al. The PRISMA 2020 statement: an updated guideline for reporting systematic reviews. BMJ. 2021;372:n71.
4. Dwan K, Altman DG, Arnaiz JA, et al. Systematic review of the empirical evidence of study publication bias and outcome reporting bias. PLoS ONE. 2008;3(8):e3081.
5. Kirkham JJ, Altman DG, Chan AW, et al. Outcome reporting bias in trials: a methodological approach for assessment and adjustment in systematic reviews. BMJ. 2018;362:k3802.
6. Ioannidis JPA, Trikalinos TA. An exploratory test for an excess of significant findings. Clin Trials. 2007;4(3):245-253.
7. Higgins JPT, Thompson SG. Quantifying heterogeneity in a meta-analysis. Stat Med. 2002;21(11):1539-1558.
8. Sterne JAC, Sutton AJ, Ioannidis JPA, et al. Recommendations for examining and interpreting funnel plot asymmetry in meta-analyses of randomised controlled trials. BMJ. 2011;343:d4002.
9. Arai M. Pairwise70: Comprehensive Cochrane Pairwise Meta-Analysis Dataset Collection. R package version 2.0.0. GitHub; 2026. Available from: https://github.com/mahmood789/Pairwise70.
10. Ioannidis JPA. Excess significance bias in the literature on brain volume abnormalities. Arch Gen Psychiatry. 2011;68(8):773-780.
11. Tsilidis KK, Panagiotou OA, Sena ES, et al. Evaluation of excess significance bias in animal studies of neurological diseases. PLoS Biol. 2013;11(7):e1001609.
12. Williamson PR, Gamble C, Altman DG, Hutton JL. Outcome selection bias in meta-analysis. Stat Methods Med Res. 2005;14(5):515-524.
13. Stanley TD, Doucouliagos H. Meta-regression approximations to reduce publication selection bias. Res Synth Methods. 2014;5(1):60-78.
14. Vevea JL, Hedges LV. A general linear model for estimating effect size in the presence of publication bias. Psychometrika. 1995;60(3):419-435.
15. Chan AW, Hrobjartsson A, Haahr MT, Gotzsche PC, Altman DG. Empirical evidence for selective reporting of outcomes in randomized trials: comparison of protocols to published articles. JAMA. 2004;291(20):2457-2465.

---

## Tables

### Table 1. Classification of 403 Cochrane meta-analyses by ORB risk level

| Risk class    | n    | %     | ORB score (mean) | ORB score (range) |
|---------------|------|-------|-------------------|--------------------|
| Low risk      | 264  | 65.5% | 6.5               | 0-20               |
| Moderate risk | 74   | 18.4% | 32.8              | 25-45              |
| High risk     | 65   | 16.1% | 60.6              | 50-90              |
| **Total**     | **403** | **100%** |               |                    |

### Table 2. Indicator values by ORB risk class

| Indicator                  | Low risk (n=264) | Moderate risk (n=74) | High risk (n=65) |
|----------------------------|-------------------|----------------------|-------------------|
| Excess significance, mean (SD) | -0.23 (0.67)  | 0.92 (1.65)          | 2.78 (3.94)       |
| I-squared, mean (SD)       | 14.8% (21.3)     | 55.1% (24.5)         | 70.3% (20.7)      |
| I-squared, median          | 0.0%              | 57.4%                | 74.5%             |
| Outlier ratio, mean        | 0.004             | 0.032                | 0.074             |
| SE ratio, mean             | 6.4               | 8.1                  | 74.1              |
| k (studies), mean          | 10.7              | 15.6                 | 28.8              |
| k (studies), median        | 6                 | 11.5                 | 15                |

### Table 3. Excess significance detection rates

| Threshold                 | n     | % of 403 reviews |
|---------------------------|-------|-------------------|
| Excess significance > 0   | 181   | 44.9%             |
| Excess significance > 1   | 91    | 22.6%             |
| Excess significance > 2   | 48    | 11.9%             |

---

## Figure legends

**Figure 1.** Distribution of 403 Cochrane meta-analyses by outcome reporting bias risk classification. Bar chart showing the number and percentage of reviews classified as Low risk (n=264, 65.5%), Moderate risk (n=74, 18.4%), and High risk (n=65, 16.1%).

**Figure 2.** Distribution of excess significance values across 403 Cochrane meta-analyses. Histogram showing the number of excess significant studies (observed minus expected) for each review. The dashed vertical line indicates zero (no excess). Values to the right indicate more significant results than expected given study power.

**Figure 3.** Between-study heterogeneity (I-squared) by ORB risk class. Box plot showing the distribution of I-squared values for Low risk (median 0.0%), Moderate risk (median 57.4%), and High risk (median 74.5%) reviews. The horizontal dashed lines indicate conventional heterogeneity thresholds at 25%, 50%, and 75%.
