"""Synthetic contract tests, NOT aerodynamic or experimental observations."""
import copy
import unittest
from validation.compare_design_rankings import compare, RankingInputError


def fixture():
    def row(name, cd):
        point = dict(cl=1., cd=cd, cd_abs_bound=.0001, provenance="synthetic unit test")
        return dict(id=name, geometry_sha256="a"*64, geometry_file=name+".dat",
                    reference=point, model=dict(point, converged=True, local_supersonic=False))
    return dict(schema_version=1, reference_kind="experiment", source="SYNTHETIC",
                geometry_audited=True, conditions_audited=True, bounds_audited=True,
                geometry_audit="synthetic", conditions_audit="synthetic", bounds_audit="synthetic",
                target_cl=1., cl_tolerance=2e-5, geometries=[row("a", .02), row("b", .019)])


class DesignRankingTests(unittest.TestCase):
    def test_resolved_agreement_and_counts(self):
        pair = compare(fixture())["pairs"][0]
        self.assertEqual(pair["ranking"], "agrees")
        self.assertAlmostEqual(pair["reference"]["delta_drag_counts"], -10.)
        self.assertAlmostEqual(pair["reference"]["absolute_bound"], .0002)

    def test_disagreement(self):
        data = fixture()
        data["geometries"][1]["model"]["cd"] = .021
        self.assertEqual(compare(data)["pairs"][0]["ranking"], "disagrees")

    def test_uncertainty_prevents_false_success(self):
        data = fixture()
        data["geometries"][1]["reference"]["cd_abs_bound"] = .002
        self.assertEqual(compare(data)["pairs"][0]["ranking"], "unresolved")

    def test_audit_gates(self):
        for field in ("geometry_audited", "conditions_audited", "bounds_audited"):
            data = fixture()
            data[field] = False
            with self.assertRaises(RankingInputError):
                compare(data)

    def test_invalid_solves_and_lift_mismatch(self):
        for field, value in (("converged", False), ("local_supersonic", True),
                             ("cd", None), ("cd", float("nan")), ("cl", 1.01)):
            data = fixture()
            data["geometries"][1]["model"][field] = value
            with self.assertRaises(RankingInputError):
                compare(data)

    def test_missing_bounds_and_duplicate_geometry(self):
        data = fixture()
        data["geometries"][1]["reference"]["cd_abs_bound"] = 0
        with self.assertRaises(RankingInputError):
            compare(data)
        data = fixture()
        data["geometries"][1] = copy.deepcopy(data["geometries"][0])
        with self.assertRaises(RankingInputError):
            compare(data)


if __name__ == "__main__":
    unittest.main()
