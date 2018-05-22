import numpy as np
import math
import copy
from sim import model
from matplotlib import pyplot as plt

n = 1
R = 1.4
Y = 1000

def manyObservables(n):
    results = []
    for i in range(0, n):
        o = observable()
        if o is not None:
            results.append(o)
    plt.hist(results)
    plt.show()

def observable():
    shower = makeJets(model())
    jets, beam = recombine(shower)
    energeticJet = jets[0]
    for i in range(1, len(jets)):
        if energy(jets[i].p) > energy(energeticJet.p):
            energeticJet = jets[i]
    js = energeticJet.subjets
    if(len(js) < 2):
        print ("not enough jet components")
        return None
    j1 = js[0]
    j2 = js[1]
    if transverse_momentum(j1.p, beam) < transverse_momentum(j2.p, beam):
        tmp = j1
        j1 = j2
        j2 = j1 #j1 is biggest, j2 second biggest
    for i in range(2, len(js)):
        ji = js[i]
        if transverse_momentum(ji.p, beam) > transverse_momentum(j1.p, beam):
            j2 = j1
            j1 = ji
        elif transverse_momentum(ji.p, beam) > transverse_momentum(j2.p, beam):
            j2 = ji
    #now have correct j1, j2
    sm = angular_distance(j1.p, j2.p, beam) * energy(j1.p) * energy(j2.p)
    return sm


def recombine(momentums):
    jets = []
    beam = np.zeros((3, 1)) 
    beam.shape = (3, 1)
    for a in momentums:
        beam += a.p
    while len(momentums) > 0:
        best_ij = 1000000000000
        the_i = 0
        the_j = 0
        for i in range(0, len(momentums)):
            for j in range(i + 1, len(momentums)):
                dij = eq4(momentums[i].p, momentums[j].p, beam)
                if dij < best_ij:
                    best_ij = dij
                    the_i = i
                    the_j = j
        best_iB = 1000000000000
        beam_i = 0
        for i in range(0, len(momentums)):
            diB = eq5(momentums[i].p, beam)
            if diB < best_iB:
                best_iB = diB
                beam_i = i
        #print("best interjet", best_ij, "best jet-beam", best_iB)
        if best_ij < best_iB:
            if the_i < the_j:
                t = the_i
                the_i = the_j
                the_j = t #i is now greater than j and can be popped first
            first = momentums.pop(the_i)
            second = momentums.pop(the_j)
            momentums.append(Jet.add(first, second))
            #print("W")
        else:
            #print("L")
            jets.append(momentums.pop(beam_i))
    return jets, beam


def makeJets(momentums):
    jets = []
    for m in momentums:
        jets.append(Jet.fromP(m))
    return jets

def energy(v):
    return math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)

def eq5(i, beam):
    return transverse_momentum(i, beam) ** (2 * n)

def eq4(i, j, beam):
    pi = transverse_momentum(i, beam) ** (2 * n)
    pj = transverse_momentum(j, beam) ** (2 * n)
    return min(pi, pj) * angular_distance(i, j, beam) / R


def angular_distance(i, j, beam):
    """Calculates the angular distance between i and j in relation to beam"""
    normal_b_i = np.cross(beam, i, axis=0)
    #proj_j = np.dot(normal_b_i, i) / np.linalg.norm(i) / np.linalg.norm(normal_b_i)
    d_phi = theta(normal_b_i, j)
    d_rapidity = prapidity(theta(beam, i)) - prapidity(theta(beam, j))
    angular_distance = math.sqrt(d_phi ** 2 + d_rapidity ** 2)
    return angular_distance

def transverse_momentum(v, beam):
    """Calculates the norm of the component of v that is orthogonal to beam"""
    return np.linalg.norm(v - np.matmul(v, np.transpose(beam)) * beam / np.linalg.norm(beam))


def theta(a, b):
    """Calculates the angle between two vectors"""
    a = np.dot(np.transpose(a), b) / np.linalg.norm(a) / np.linalg.norm(b)
    return math.acos(a)

def prapidity(theta):
    """Calculates pseudorapidity for a given theta"""
    return -math.log(math.tan(theta / 2))

class Jet:
    @classmethod
    def fromP(cls, momentum):
        obj = cls()
        obj.p = momentum
        obj.subjets = [obj]
        return obj
    @classmethod
    def fromTwo(cls, one, other):
        self = cls()
        self.p = one.p + other.p
        self.subjets = []
        self.subjets += one.subjets
        self.subjets += other.subjets
        return self

    def add(self, other):
        return Jet.fromTwo(self, other)

if __name__ == "__main__":
    #observable()
    manyObservables(1000)
