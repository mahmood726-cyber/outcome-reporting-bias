"""Tests for ORB detection pipeline."""

import math
import sys
from types import SimpleNamespace

import numpy as np
import pytest

sys.path.insert(0, str(__import__('pathlib').Path(__file__).parent.parent))

from src.pipeline import analyze_orb_potential, run_pipeline, resolve_paths


def _synthetic_review():
    return SimpleNamespace(
        review_id='CD000001',
        analysis_name='Synthetic review',
        k=5,
        yi=np.array([-0.6, -0.4, -0.3, -0.2, -0.1]),
        sei=np.array([0.1, 0.12, 0.1, 0.14, 0.11]),
    )


class TestORBAnalysis:
    def test_pipeline_runs_on_synthetic_review(self):
        results = analyze_orb_potential([_synthetic_review()])
        assert len(results) == 1
        assert results[0]['orb_class'] in ('Low_Risk', 'Moderate_Risk', 'High_Risk')
        assert 0 <= results[0]['I2'] <= 100
        assert results[0]['orb_score'] >= 0

    def test_excess_sig_calculation(self):
        """Excess significance should be finite and reasonable."""
        results = analyze_orb_potential([_synthetic_review()])
        assert math.isfinite(results[0]['excess_significance'])
        assert -20 < results[0]['excess_significance'] < 20


def test_run_pipeline_uses_repo_relative_sibling_projects(tmp_path, monkeypatch):
    projects_root = tmp_path / 'projects'
    project_root = projects_root / 'OutcomeReportingBias'
    project_root.mkdir(parents=True)

    paths = resolve_paths(project_root=project_root, projects_root=projects_root)
    paths['pairwise_dir'].mkdir(parents=True, exist_ok=True)
    monkeypatch.setattr('src.pipeline.load_all_reviews', lambda pairwise_dir, min_k=3: [_synthetic_review()])

    results, summary = run_pipeline(project_root=project_root, projects_root=projects_root)

    assert len(results) == 1
    assert summary['n_reviews'] == 1
    assert (project_root / 'data' / 'output' / 'orb_results.csv').exists()
    assert (project_root / 'data' / 'output' / 'orb_summary.json').exists()


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
