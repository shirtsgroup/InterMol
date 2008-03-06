import protein
import time
import sys
N = int(sys.argv[1])
p = protein.ProteinClass("ACDEFGHIKLMNPQRSTVWY")
t1 = time.time()
for i in range(N):
 E = p.Energy()
 E = p.RotateToPhiPsi(3,-60.,-45.)
t2 = time.time()
print t2-t1
