## \example core/simple.py
# Illustration of simple usage of the IMP library from Python.
from pathlib import Path
import IMP
import IMP.algebra
import IMP.core
import sys
import random
random.seed(0)
sys.path.append(str(Path(Path.home(), "xray/src")))
from log_statistics import LogStatistics


IMP.setup_from_argv(sys.argv, "simple example")

m = IMP.Model()

# Create two "untyped" particles
p1 = m.add_particle('p1')
p2 = m.add_particle('p2')
p3 = m.add_particle('p3')

# "Decorate" the particles with x,y,z attributes (point-like particles)
d1 = IMP.core.XYZ.setup_particle(m, p1)
d2 = IMP.core.XYZ.setup_particle(m, p2)
d3 = IMP.core.XYZ.setup_particle(m, p3)

# Use some XYZ-specific functionality (set coordinates)
for d in [d1, d2, d3]:
    d.set_coordinates(IMP.algebra.Vector3D(random.random()*10, random.random()*10, random.random()*10))
    d.set_coordinates_are_optimized(True)

print(d1, d2, d3)

rs = list()
# Harmonically restrain p1 to be zero distance from the origin
f = IMP.core.Harmonic(0.0, 1.0)
s = IMP.core.DistanceToSingletonScore(f, IMP.algebra.Vector3D(0., 0., 0.))
rs.append(IMP.core.SingletonRestraint(m, s, p1))

# Harmonically restrain p1 and p2 to be distance 5.0 apart
for pair_1, pair_2 in ((p1, p2), (p1, p3), (p2, p3)):
    f = IMP.core.Harmonic(5.0, 1.0)
    s = IMP.core.DistancePairScore(f)
    r = IMP.core.PairRestraint(m, s, (p1, p2))
    print(r.evaluate(calc_derivs=False))
    rs.append(r)

# Optimize the x,y,z coordinates of both particles with conjugate gradients
sf = IMP.core.RestraintsScoringFunction(rs, "scoring function")

o = IMP.core.ConjugateGradients(m)
o.set_scoring_function(sf)

log = LogStatistics(
    m=m,
    pids=[p1, p2],
    m_0=None,
    pids_0=None,
    restraints=rs,
    freq=1,
    silent=False
)
log.update_always()
o.add_optimizer_state(log)

# o.set_log_level(IMP.TERSE)
o.set_gradient_threshold(0)
o.optimize(max_steps=25)

# float_keys = o.get_optimized_attributes()
pid_1 = m.get_particle_indexes()[0]
float_key = d.get_xyz_keys()[0]
print(type(pid_1), type(float_key))
float_id = IMP.FloatIndex(pid_1, float_key)
o.get_scaled_derivative(float_id)

print(log.get_period())
print(log.get_number_of_updates())

# o.remove_optimizer_state(log)
print(d1, d2, d3)
del o
del m
