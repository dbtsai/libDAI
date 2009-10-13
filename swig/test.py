import dai

v1 = dai.Var(0,4)
v2 = dai.Var(1,3)
vs = dai.VarSet(v1,v2)
f = dai.Factor(vs)
fVec = dai.VecFactor()
fVec.append(f)
fg = dai.FactorGraph(fVec)
