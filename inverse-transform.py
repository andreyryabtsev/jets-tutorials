import numpy as np
import math, random
from scipy.special import erfinv

# use probability function N = exp(-(x - 5) ** 2 / 8) / math.sqrt(8 * math.pi)
# CDF is then: ( math.erf((x - 5) / 2 / math.sqrt(2) + math.erf(5 / 2 / math.sqrt(2) ) / 2
# for x = [0, 10], y = [CDF(0), CDF(10)] = [0, math.erf(5 / 2 / math.sqrt(2))]
#
# the inverse of the CDF is 2 * math.sqrt(2) * erfinv(2y - math.erf(5 / 2 / math.sqrt(2))) + 5

def CDF_inv(y):
    return 2 * math.sqrt(2) * erfinv(2 * y - math.erf(5 / 2 / math.sqrt(2))) + 5



#Find a Gaussian random number with mean = 6, sd = 2, range = [0, 10]
ylow = 0
yhigh = math.erf(5 / 2 / math.sqrt(2))
y = random.random() * (yhigh - ylow) + ylow
x = CDF_inv(y)
print(x)
