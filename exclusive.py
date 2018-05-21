import numpy as np
import math
import copy
from sim import model
from matplotlib import pyplot as plt

n = 1
d_cut = 1
R = 0.05
Y = 1000

def simulate_recombine():
    global n, R
    Ns = [-1, 0, 1]
    Rs = [0.01, 0.05, 0.1, 0.5, 1]
    showers = []
    for i in range(0, Y):
        showers.append(model())
    print("Done generating, recombining...")
    for a in range(3):
        n = Ns[a]
        for b in range(5):
            R = Rs[b]
            jets = np.zeros(Y)
            for i in range(Y):
                if i > 0 and i % 250 == 0:
                    print("Recombining, " + str(i / 10) + "% done...")
                momentums = copy.deepcopy(showers[i])
                hh = recombine(momentums)
                jets[i] = len(hh)
            print("Recombined with preset " + str(1 + 5*a + b) + "/15...")
            plt.subplot(3, 5, 1 + 5 * a + b)
            plt.hist(jets)
            plt.title('n = ' + str(n) + ', R = ' + str(R))
    plt.tight_layout()
    plt.show()


def recombine(momentums):
    while len(momentums) > 1:
        beam = np.zeros((3, 1)) 
        beam.shape = (3, 1)
        for a in momentums:
            beam += a
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
        if best_ij > d_cut and best_iB > d_cut:
            break
        #print("best interjet", best_ij, "best jet-beam", best_iB)
        if best_ij < best_iB:
            if the_i < the_j:
                t = the_i
                the_i = the_j
                the_j = t #i is now greater than j and can be popped first
            first = momentums.pop(the_i)
            second = momentums.pop(the_j)
            momentums.append(first + second)
        else:
            del momentums[beam_i]
            #momentums.append(momentums.pop(beam_i))
    return momentums




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


if __name__ == "__main__":
    simulate_recombine()
