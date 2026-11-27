import sys
sys.path.insert(0,'C:/Users/nrvig/Projects/math-commons/research/marcus-de-oliveira/zheng-n4-replay/cases-v-b')
import replay_T38_common as c
T=(0,1,5,9,14,18)
print(c.generators(T), c.orbit_rep(T))
d,rays=c.extreme_rays_of_contact_cone(T)
print(d,len(rays))
gens=c.generators(T)
print([ (j,c.f_name(j)) for j in gens])
