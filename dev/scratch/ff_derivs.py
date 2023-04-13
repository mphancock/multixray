import sys
from pathlib import Path
import IMP
import IMP.core
import IMP.atom
import IMP.algebra
sys.path.append(str(Path(Path.home(), "xray/src")))
from charmm import charmm_restraints


if __name__ == "__main__":
    pdb_file = Path(Path.home(), "xray/dev/scratch/one_res.pdb")
    m = IMP.Model()
    h = IMP.atom.read_pdb(str(pdb_file), m)

    pids = IMP.atom.Selection(h).get_selected_particle_indexes()
    rs = charmm_restraints(
        m=m,
        h=h,
        eps=False
    )

    xyz = IMP.core.XYZ(m, pids[0])
    rset = IMP.RestraintSet(rs, 2)
    score = rset.evaluate(calc_derivs=True)
    dff_dx = xyz.get_derivative(0)
    print(dff_dx)

    score_tot = 0
    dff_tot_dx = 0
    for r in rs:
        print(r.get_name())
        score_term = r.evaluate(calc_derivs=True)
        score_tot = score_tot + score_term

        dff_term_dx = xyz.get_derivative(0)
        dff_tot_dx = dff_tot_dx + dff_term_dx
        # print()
        print(xyz.get_derivative(0))

        break

    print(dff_tot_dx)

    # Check if you can manually reset derivatives
    print(xyz.get_derivative(0))
    da = IMP.DerivativeAccumulator(1)
    xyz.add_to_derivatives(-xyz.get_derivatives(), da)
    da = IMP.DerivativeAccumulator(2)
    xyz.add_to_derivatives(IMP.algebra.Vector3D(10,10,5), da)
    print(xyz.get_derivative(0))

