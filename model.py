import random, math

E_CRIT = 1000

def model(E, px, py, pz):
    pass


def travel(E, px, py, pz):

    pass



def decay(E, px, py, pz):

    pass


def random_pair():
    return randomInverse(1), randomInverse(math.pi / 2)


def randomInverse(endValue):
    #PDF = 1/(1+x) CDF = ln(x+1) CDF_INV = e^x - 1
    M = 0.001
    ylow = 0
    yhigh = math.log((endValue + M) / M)
    rnd = random.random() * (yhigh - ylow) + ylow
    return M * (math.exp(rnd) - 1)


print(random_pair())
