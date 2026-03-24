"""Tests for ORB detection pipeline."""
import sys, math, numpy as np, pytest
sys.path.insert(0, str(__import__('pathlib').Path(__file__).parent.parent))


class TestORBAnalysis:
    def test_pipeline_runs(self):
        from src.pipeline import analyze_orb_potential
        sys.path.insert(0, r'C:\FragilityAtlas')
        from src.loader import load_review
        reviews = []
        r = load_review(r'C:\Models\Pairwise70\data\CD000028_pub4_data.rda')
        if r:
            reviews.append(r)
        results = analyze_orb_potential(reviews)
        assert len(results) == 1
        assert results[0]['orb_class'] in ('Low_Risk', 'Moderate_Risk', 'High_Risk')
        assert 0 <= results[0]['I2'] <= 100
        assert results[0]['orb_score'] >= 0

    def test_excess_sig_calculation(self):
        """Excess significance should be finite and reasonable."""
        from src.pipeline import analyze_orb_potential
        sys.path.insert(0, r'C:\FragilityAtlas')
        from src.loader import load_review
        r = load_review(r'C:\Models\Pairwise70\data\CD000028_pub4_data.rda')
        if r:
            results = analyze_orb_potential([r])
            assert math.isfinite(results[0]['excess_significance'])
            assert -20 < results[0]['excess_significance'] < 20


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
