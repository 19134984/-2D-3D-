"""Exact D2Q9 moment audit of arXiv:2602.06686v1, Eqs. (6), (9), (30).

No numerical solver is modified by this audit. Uses rational arithmetic.
"""
from fractions import Fraction as Q

e = [(0, 0), (1, 0), (0, 1), (-1, 0), (0, -1),
     (1, 1), (-1, 1), (-1, -1), (1, -1)]
w = [Q(4, 9)] + [Q(1, 9)] * 4 + [Q(1, 36)] * 4
cs2 = Q(1, 3)
u = [Q(1, 50), Q(-3, 100)]
force = [Q(1, 10000), Q(1, 5000)]
snu, sq = Q(5, 3), Q(8, 19)  # tau_plus=0.6, magic=3/16
uf = sum(a*b for a, b in zip(u, force))

paper_source = []
for (x, y), weight in zip(e, w):
    eu = x*u[0] + y*u[1]
    ef = x*force[0] + y*force[1]
    paper_source.append(weight*((1-snu/2)*eu*ef/cs2**2
                                + (1-sq/2)*(ef-uf)/cs2))

mass_source = sum(paper_source)
momentum_source = [sum(c[k]*s for c, s in zip(e, paper_source))
                   for k in range(2)]
assert mass_source == (sq-snu)*uf/(2*cs2)
assert momentum_source == [(1-sq/2)*f for f in force]
# In Eq. (30), both regularized terms have zero first moment under Eq. (9).
# Therefore p_post = rho*u + first_moment(S), while p_pre=rho*u-F/2.
paper_increment = [f/2+s for f, s in zip(force, momentum_source)]
assert paper_increment == [(Q(3, 2)-sq/2)*f for f in force]
assert mass_source != 0
assert paper_increment != force

print('Paper Eq. (6) mass source:', mass_source, '=', float(mass_source))
print('Required mass source: 0')
print('Paper Eqs. (6)+(9)+(30) momentum gain / F:', Q(3, 2)-sq/2)
print('Required momentum gain / F: 1')

# A force-compatible, conserved-moment formulation:
# f_post = feq + P1(fneq) + (1-snu)*P2(fneq) + (1-sq)*P3(fneq)
#          + wi*(e.F)/cs2
#          + (1-snu/2)*wi*[(e.u)*(e.F)/cs2**2-u.F/cs2].
source = [weight*((x*force[0]+y*force[1])/cs2
          +(1-snu/2)*((x*u[0]+y*u[1])*(x*force[0]+y*force[1])/cs2**2-uf/cs2))
          for (x, y), weight in zip(e, w)]
assert sum(source) == 0
assert [sum(c[k]*s for c, s in zip(e, source)) for k in range(2)] == force
print('Conserved-moment source: exact zero mass and exact momentum F: PASS')
