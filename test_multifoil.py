"""Physics-based checks, run with: python -m unittest -v test_multifoil.py"""
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from multifoil import (PanelSolver, ViiflowSolver, load_case, naca4, read_dat,
                       repanel, transform, validate_geometry)


def joukowski(panels=360, epsilon=.08):
    """Circle centered at -epsilon, radius 1+epsilon, mapped by z=w+1/w."""
    radius = 1+epsilon
    theta = np.linspace(0,2*np.pi,panels+1)
    w = -epsilon+radius*np.exp(1j*theta)
    z = w+1/w
    chord = z.real.max()-z.real.min()
    xy = np.column_stack(((z.real-z.real.min())/chord,z.imag/chord))
    xy[-1] = xy[0]
    return xy,radius,chord


class PanelPhysics(unittest.TestCase):
    def test_symmetric_zero_lift_and_antisymmetric_polar(self):
        p = PanelSolver([("main",naca4("0012",201))])
        self.assertLess(abs(p.solve(0).cl),1e-10)
        self.assertAlmostEqual(p.solve(-5).cl,-p.solve(5).cl,places=10)
        r = p.solve(5)
        self.assertLess(r.residual,1e-10)
        self.assertLess(r.diagnostics["kutta_residual"],1e-10)
        self.assertIsNone(r.cd)

    def test_joukowski_exact_lift(self):
        xy,radius,chord = joukowski()
        result = PanelSolver([("exact",xy)]).solve(5)
        exact = 8*np.pi*radius/chord*np.sin(np.deg2rad(5))
        self.assertLess(abs(result.cl/exact-1),.0001)
        self.assertLess(abs(result.diagnostics["cl_pressure_total"]/exact-1),.005)
        self.assertLess(abs(result.cd_pressure_numerical),.002)

    def test_panel_refinement_and_pressure_drag(self):
        results = [PanelSolver([("main",naca4("0012",n+1))]).solve(5) for n in (80,160,320)]
        self.assertLess(abs(results[2].cl-results[1].cl),abs(results[1].cl-results[0].cl))
        self.assertLess(abs(results[2].cd_pressure_numerical),abs(results[0].cd_pressure_numerical)/3)

    def test_multi_element_order_independence_and_coupling(self):
        a = ("main",naca4("2412",201))
        b = ("flap",transform(naca4("0012",141),.3,-20,translation=(.93,-.055)))
        p = PanelSolver([a,b]).solve(4)
        q = PanelSolver([b,a]).solve(4)
        self.assertAlmostEqual(p.cl,q.cl,places=10)
        self.assertAlmostEqual(p.cm,q.cm,places=10)
        uncoupled = PanelSolver([a]).solve(4).cl+PanelSolver([b]).solve(4).cl
        self.assertGreater(abs(p.cl-uncoupled),.1)
        self.assertLess(p.residual,1e-10)

    def test_far_separated_elements(self):
        xy = naca4("0012",141)
        one = PanelSolver([("a",xy)]).solve(5).cl
        two = PanelSolver([("a",xy),("b",xy+[0,100])]).solve(5).cl
        self.assertLess(abs(two/(2*one)-1),1e-3)

    def test_mirror_symmetry_of_whole_configuration(self):
        e = [("main",naca4("2412",161)),
             ("flap",transform(naca4("0012",121),.3,-20,translation=(.93,-.055)))]
        mirrored = [(name,(xy*[1,-1])[::-1].copy()) for name,xy in e]
        a,b = PanelSolver(e).solve(4),PanelSolver(mirrored).solve(-4)
        self.assertAlmostEqual(a.cl,-b.cl,places=9)
        self.assertAlmostEqual(a.cm,-b.cm,places=9)

    def test_reference_chord_normalization(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/"case.json"
            def calc(scale):
                path.write_text(json.dumps({"reference_chord":scale,
                    "moment_reference":[scale*.25,0],
                    "elements":[{"naca":"2412","scale":scale,"panels":160}]}))
                _,e,ref = load_case(path)
                return PanelSolver(e,ref).solve(4)
            a,b = calc(1),calc(1000)
            self.assertAlmostEqual(a.cl,b.cl,places=10)
            self.assertAlmostEqual(a.cm,b.cm,places=10)

    def test_surface_and_off_body_velocity_consistency(self):
        p = PanelSolver([("main",naca4("0012",161))])
        r = p.solve(3)
        v = p.velocity(p.mid+p.n*1e-7,3)
        self.assertLess(np.max(np.abs(np.sum(v*p.n,axis=1))),1e-3)
        self.assertLess(np.max(abs(np.sum(v*p.t,axis=1)-r.surfaces[0]["ue_signed"])),1e-3)

    def test_collision_rejection(self):
        xy = naca4("0012",101)
        with self.assertRaises(ValueError):
            validate_geometry([("a",xy),("b",xy+[.1,0])])
        with self.assertRaises(ValueError):
            validate_geometry([("a",xy),("b",xy*.1+[.4,0])])

    def test_selig_orientation(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/"foil.dat"
            np.savetxt(path,naca4("2412",101)[::-1],header="NACA 2412",comments="")
            xy = read_dat(path)
            self.assertGreater(PanelSolver([("main",xy)]).solve(0).cl,0)

    def test_blunt_te_closure(self):
        xy = naca4("0012",201)
        xy[0,1] = .0005
        xy[-1,1] = -.0005
        p = PanelSolver([("main",xy)])
        r = p.solve(0)
        self.assertLess(abs(r.cl),1e-10)
        self.assertEqual(p.ns,len(xy))

    def test_examples(self):
        for case in ["two_element","three_element"]:
            _,e,ref = load_case(Path(__file__).parent/"examples"/(case+".json"))
            self.assertTrue(PanelSolver(e,ref).solve(4).converged)


if __name__ == "__main__":
    unittest.main()
