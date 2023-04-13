import IMP, IMP.atom, IMP.core
from xray_restraint import XtalRestraint
from charmm import charmm_restraints


def score(
        pdb_file,
        uc_dim,
        sg_symbol,
        cif_file,
        scale=False
):
    m = IMP.Model()
    # sel = IMP.atom.AndPDBSelector(IMP.atom.NonWaterNonHydrogenPDBSelector(), IMP.atom.ATOMPDBSelector())
    sel = IMP.atom.ATOMPDBSelector()
    h = IMP.atom.read_pdb(str(pdb_file), m, sel)

    for pid in m.get_particle_indexes():
        if IMP.atom.Atom.get_is_setup(m, pid):
            IMP.core.XYZR(m, pid).set_coordinates_are_optimized(True)

    ps = [m.get_particle(pid) for pid in m.get_particle_indexes() if
          IMP.atom.Atom.get_is_setup(m, pid)]

    rs_charmm = IMP.RestraintSet(m, 1.0)
    rs_charmm.add_restraints(charmm_restraints(m,h))
    charmm_score = rs_charmm.evaluate(False)

    rs_xtal = IMP.RestraintSet(m, 1.0)
    r_xtal = XtalRestraint(
        m=m,
        ps=ps,
        uc_dim=uc_dim,
        sg_symbol=sg_symbol,
        f_obs_file=str(cif_file),
        scale=scale
    )
    rs_xtal.add_restraint(r_xtal)
    xray_score = rs_xtal.unprotected_evaluate(None)

    return xray_score

    # scores = {"xtal": xtal_score, "mm": charmm_score}
    # print(score_input["pdb"], xtal_score, charmm_score)
    # return score_input["pdb"], scores
