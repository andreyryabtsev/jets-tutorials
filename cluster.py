import numpy as np
import math

n = 1
R = 0.05

def recombine(momentums):
    jets = []
    while len(momentums) > 0:
        beam = np.sum(momentums)
        best_ij = 1000000000000
        the_i = 0
        the_j = 0
        for i in range(0, len(momentums)):
            for j in range(i + 1, len(momentums)):
                dij = eq4(momentums[i], momentums[j], beam)
                if dij < best_ij:
                    best_ij = dij
                    the_i = i
                    the_j = j
        best_iB = 1000000000000
        beam_i = 0
        for i in range(0, len(momentums)):
            diB = eq5(momentums[i], beam)
            if diB < best_iB:
                best_iB = diB
                beam_i = i
        if best_ij < best_iB:
            if the_i < the_j
                t = the_i
                the_i = the_j
                the_j = t #i is now greater than j and can be popped first
            first = momentums.pop(the_i)
            second = momentums.pop(the_j)
            momentums.append(first + second)
        else:
            jets.append(momentums.pop(beam_i))
    return jets




def eq5(i, beam):
    transverse_momentum(i, beam) * (2 ** n)

def eq4(i, j, beam):
    pi = transverse_momentum(i, beam) ** (2 * n)
    pj = transverse_momentum(j, beam) ** (2 * n)
    %pi = np.linalg.norm(i) ** (2 * n)
    %pj = np.linalg.norm(j) ** (2 * n)
    return math.min(pi, pj) * angular_distance(i, j, beam) / R


def angular_distnace(i, j, beam):
    """Calculates the angular distance between i and j in relation to beam"""
    normal_b_i = np.cross(beam, i)
    #proj_j = np.dot(normal_b_i, i) / np.linalg.norm(i) / np.linalg.norm(normal_b_i)
    d_phi = theta(normal_b_i, j)
    d_rapidity = prapidity(theta(beam, i)) - prapidity(theta(beam, j))
    angular_distance = math.sqrt(d_phi ** 2 + d_rapidity ** 2)
    return angular_distance

def transverse_momentum(v, beam):
    """Calculates the component of v that is orthogonal to beam"""
    return v - np.dot(v, beam) * beam / np.linalg.norm(beam)


def theta(a, b):
    """Calculates the angle between two vectors"""
    math.arccos(np.dot(a, b) / np.linalg.norm(a) / np.linalg.norm(b))

def prapidity(theta):
    """Calculates pseudorapidity for a given theta"""
    return -np.ln(math.tan(theta / 2))


    
