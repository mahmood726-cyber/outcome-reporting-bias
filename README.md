# Outcome Reporting Bias Risk in 403 Cochrane Meta-Analyses

What is the prevalence of statistical patterns consistent with outcome reporting bias across Cochrane systematic reviews? We applied excess significance testing and three complementary indicators to 403 Cochrane reviews from the Pairwise70 dataset, each containing at least three primary studies. A composite scoring system combining excess significance, heterogeneity, outlier ratio, and precision asymmetry classified reviews into Low, Moderate, and High risk categories. The prevalence of High risk was 16.1% (65 of 403, 95% CI 12.6-20.1%), with mean excess significance of 2.83 and mean I-squared of 70.3% compared to 14.8% in Low risk reviews. Spearman correlation between excess significance and heterogeneity was 0.36, confirming that reviews with more significant results than expected showed greater between-study variability. These findings suggest that one in six Cochrane reviews exhibits statistical signatures warranting closer scrutiny of pre-registration fidelity and outcome switching. A limitation is that statistical proxy indicators cannot distinguish outcome reporting bias from publication bias or other selective reporting mechanisms.

**Live dashboard:** <https://mahmood726-cyber.github.io/outcomereportingbias/>

## Run

Open `index.html` (or `index.html`) in any modern browser. No build step.

For local development:

```bash
python -m http.server 8000
# then open http://localhost:8000/
```

## Test

```bash
python -m pytest -q
```

The suite under `tests/` includes 1 test file(s).

## Repo layout

| Path | Purpose |
|---|---|
| `index.html` | the dashboard (main artifact) |
| `index.html` | landing page |
| `tests/` | pytest tests |
| `e156-submission/` | E156 micro-paper bundle |
| `E156-PROTOCOL.md` | project metadata (E156 entry #126) |

## License

See `LICENSE` (MIT).
