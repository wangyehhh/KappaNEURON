from neuron.rxd.generalizedReaction import GeneralizedReaction, molecules_per_mM_um3
import sys
import os
import scipy.sparse

print(sys.argv[0])
print(os.getcwd())
n = 7
n2 = 7
m = scipy.sparse.eye(n, n2)
print(m)
class kappa(GeneralizedReaction):
    def print(self):
        print(self._indices_dict)
g = GeneralizedReaction()
k = kappa()
dic2 = k.print()