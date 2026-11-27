import sys, os, pickle, time
sys.path.insert(0, r"C:\Users\nrvig\Projects\math-commons\research\marcus-de-oliveira\zheng-n4-replay\cases-v-b")
sys.dont_write_bytecode = True
import replay_T38_common as M
ck = M.Checker("proto")
which = sys.argv[1]
if which == "17":
    order = [3,1,53,6,15,17,16,4]
else:
    order = [6,2,5,53,0,54,3,15]
t=time.time()
rd = M.radial_data(ck, order, "q,p,j,c,e,d,b,a,h")
print(time.time()-t)
pickle.dump(rd, open("rd%s.pkl"%which,"wb"))
print(len(rd["C0"]), len(rd["Q"]), len(rd["G"]))
