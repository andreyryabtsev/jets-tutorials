import random, math
import numpy as np

E_CRIT = 3
T_MIN = 2
T_MAX = 6

def model(px, py, pz):
    partons = []
    P = np.array([px, py, pz])
    P.shape = (3, 1)
    original = Parton(0, 0, 0, P, random_t())
    partons.append(original)
    hadrons = []
    
    ticks = 0
    while len(partons) > 0:
        ticks += 1
        
        i = 0
        while i < len(partons):
            parton = partons[i]
            parton.t -= 1
            parton.pos += parton.P
            if parton.t <= 0:
                #splitting
                del partons[i]
                i = i - 1
                A, B = split(parton)
                A = Parton(parton.pos[0], parton.pos[1], parton.pos[2], A, random_t())
                B = Parton(parton.pos[0], parton.pos[1], parton.pos[2], B, random_t())
                if energy(A) > E_CRIT:
                    partons.insert(i, A)
                    i = i + 1
                else:
                    hadrons.append(A)
                if energy(B) > E_CRIT:
                    partons.insert(i, B)
                    i = i + 1
                else:
                    hadrons.append(B)


            i = i + 1
        if ticks % 10 == 0:
            print("t =", ticks)
            print("# of partons:", len(partons)) 
    print("# of hadrons:", len(hadrons))

def split(parton):
    """Calculate a random splitting of a parton with energy E and 3D momentum P.
    Returns two resulting partons
    """
   
    P = parton.P
    z, theta = random_pair()
    azimuthal = random.random() * math.pi

    thetaAxis = standardOrthogonal(P)
    thetaA = rotationMatrix2(thetaAxis, theta)
    thetaB = rotationMatrix2(thetaAxis, -theta)
    phiRotation = rotationMatrix2(P, azimuthal)
    matrixA = np.matmul(phiRotation, thetaA)
    matrixB = np.matmul(phiRotation, thetaB)
    
    A = z * np.matmul(matrixA, P)
    B = (1 - z) * np.matmul(matrixB, P)
    
    return A, B

def random_t():
    return T_MIN + (T_MAX - T_MIN) * random.random()

def random_pair():
    return randomInverse(1), randomInverse(math.pi / 2)


def randomInverse(endValue):
    #PDF = 1/(1+x) CDF = ln(x+1) CDF_INV = e^x - 1
    M = 0.001
    ylow = 0
    yhigh = math.log((endValue + M) / M)
    rnd = random.random() * (yhigh - ylow) + ylow
    return M * (math.exp(rnd) - 1)

def rotationMatrix(axis, angle):
    assert axis.shape == (3,1)
    axis = axis / np.linalg.norm(axis)
    return np.cos(angle) * np.identity(3) + np.sin(angle) * skewSymMatrix(axis) + (1 - np.cos(angle)) * np.tensordot(axis, axis)

def rotationMatrix2(axis, angle):
    assert axis.shape == (3, 1)
    axis = axis / np.linalg.norm(axis)
    x = axis[0]
    y = axis[1]
    z = axis[2]
    c = np.cos(angle)
    s = np.sin(angle)
    mat = np.array([[c + x ** 2 * (1 - c), x * y * (1 - c) - z * s, x * z * (1 - c) + y * s],
                     [y * x * (1 - c) + z * s, c + y ** 2 * (1 - c), y * z * (1 - c) - x * s],
                     [z * x * (1 - c) - y * s, z * y * (1 - c) + x * s, c + z ** 2 * (1 - c)]])
    mat.shape = (3, 3)
    return mat

def energy(parton):
    return np.linalg.norm(parton.P)

def skewSymMatrix(vector):
    assert vector.shape == (3,1)
    return np.array([[0, -vector[2], vector[1]], [vector[2], 0, -vector[0]], [-vector[1], vector[0], 0]])

def standardOrthogonal(vector):
    assert vector.shape == (3, 1)
    result = np.array([1, 1, -(vector[0] + vector[1]) / vector[2]])
    result.shape = (3, 1)
    return result


class Parton:
    def __init__(self, x, y, z, P, t):
        self.pos = np.array([x, y, z], dtype =float)
        self.pos.shape = (3, 1)
        self.P = P
        self.t = t



model(5, 12, 4)