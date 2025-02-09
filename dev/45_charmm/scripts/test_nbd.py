import sys
from pathlib import Path
import IMP
import IMP.atom

import IMP
import IMP.core
import IMP.algebra


# This example adds a restraint on nonbonded interactions
# after excluding a set of bonded interactions.

m = IMP.Model()
# The set of particles
ps = IMP.container.ListSingletonContainer(
    m, IMP.core.create_xyzr_particles(m, 20, 1.0))

# create a bond between two particles
bd0 = IMP.atom.Bonded.setup_particle(m, ps.get_indexes()[0])
bd1 = IMP.atom.Bonded.setup_particle(m, ps.get_indexes()[1])
IMP.atom.create_custom_bond(bd0, bd1, 2.0)

# Set up the nonbonded list for all pairs that are touching
# and let things move 3 before updating the list
nbl = IMP.container.ClosePairContainer(ps, 100.0, 3.0)

sdps = IMP.core.SoftSpherePairScore(1)
evr = IMP.container.PairsRestraint(sdps, nbl)

evr.evaluate(False)

print(nbl.get_indexes())