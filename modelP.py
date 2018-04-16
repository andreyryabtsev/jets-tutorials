import random, math
import numpy as np
import matplotlib.pyplot as plt
import time
from mpl_toolkits.mplot3d import Axes3D

E_CRIT = 3
T_MIN = 2
T_MAX = 6
X_INI = 0
Y_INI = 0
Z_INI = 0
M = 0.09
dt = 0.01
E_INI = 12
alpha = 0.5

def draw(fig, ax, partons, deadPartons, hadrons):
    ax.clear()
    ax.scatter(0, 0, 0)
    for parton in partons:
        xs, ys, zs = getLine(parton)
        ax.plot(xs, ys, zs)

    for parton in deadPartons:
        xs, ys, zs = getLine(parton)
        ax.plot(xs, ys, zs)
    for hadron in hadrons:
        xs, ys, zs = getLine(hadron)
        ax.plot(xs, ys, zs)
        ax.scatter(hadron.pos[0][0], hadron.pos[1][0], hadron.pos[2][0])
    fig.canvas.draw()

def getLine(parton):
    xs = [parton.birth[0][0], parton.pos[0][0]]
    ys = [parton.birth[1][0], parton.pos[1][0]]
    zs = [parton.birth[2][0], parton.pos[2][0]]
    return xs, ys, zs

def model():
    #GRAPHICS
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    plt.ion()
    fig.show()
    ##


    #INITIAL PROTON COLLISION
    E1 = E_INI# * random_falling(alpha * E_INI)
    #random energy between 0 and E_INI
    #theta1 = random.random() * math.pi * 2
    #theta2 = random.random() * math.pi
    P1 = np.random.rand(3, 1)
    P1 = P1 / np.linalg.norm(P1) * E1
    P2 = -P1
    print(random_falling(alpha * E_INI))

    parton1 = Parton(X_INI, Y_INI, Z_INI, P1, random_t())
    parton2 = Parton(X_INI, Y_INI, Z_INI, P2, random_t())

    ##

    partons = []
    partons.append(parton1)
    partons.append(parton2)
    hadrons = []
    deadPartons = []
    
    timeout = 0
    ticks = 0
    while timeout > 0 or  len(partons) > 0:
        ticks += 1
        timeout -= 1
        
        for hadron in hadrons:
            hadron.pos += hadron.P

        i = 0
        while i < len(partons):
            parton = partons[i]
            parton.t -= 1
            parton.pos += parton.P
            if parton.t <= 0:
                #splitting
                del partons[i]
                deadPartons.append(parton)
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
                if len(partons) == 0:
                    timeout = 10

            i = i + 1
        time.sleep(dt)
        draw(fig, ax, partons, deadPartons, hadrons)
    plt.ioff()
    print("# of hadrons:", len(hadrons))
    draw(fig, ax, partons, deadPartons, hadrons)
    input("Done simulating, any key to close...")

def split(parton):
    """Calculate a random splitting of a parton with energy E and 3D momentum P.
    Returns two resulting partons
    """
   
    P = parton.P
    z, theta = random_pair()
    azimuthal = random.random() * math.pi

    thetaAxis = standardOrthogonal(P)
    thetaA = rotationMatrix(thetaAxis, theta)
    #thetaB = rotationMatrix(thetaAxis, -theta)
    phiRotation = rotationMatrix(P, azimuthal)
    matrixA = np.matmul(phiRotation, thetaA)
    #matrixB = np.matmul(phiRotation, thetaB)
    
    A = z * np.matmul(matrixA, P)
    B = P - A
    
    return A, B

def random_t():
    return T_MIN + (T_MAX - T_MIN) * random.random()

def random_pair():
    return randomInverse(1), randomInverse(math.pi / 2)

def random_falling(a):
    y = random.random() * (1 - math.exp(-a)) / a
    return -math.log(1 - a * y) / a
    #return -math.log(1 - a * random.random()) / a

def randomInverse(endValue):
    #PDF = 1/(1+x) CDF = ln(x+1) CDF_INV = e^x - 1
    ylow = 0
    yhigh = math.log((endValue + M) / M)
    rnd = random.random() * (yhigh - ylow) + ylow
    return M * (math.exp(rnd) - 1)

def rotationMatrix(axis, angle):
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
        self.birth = np.array([x, y, z], dtype = float)
        self.birth.shape = (3, 1)

model()
