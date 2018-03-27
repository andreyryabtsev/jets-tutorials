import numpy as np
import math, random
import matplotlib.pyplot as plt

# use probability function N = math.exp(-x) 
# CDF = 1 - math.exp(-x)
# CDF_inv = x = -ln(1-y)
# range of x = [0, 5] --> range of y = 
NUM = 100

def prob(x):
    return math.exp(-x)

def CDF(x):
    return 1 - math.exp(-x)

def CDF_inv(y):
    return -math.log(1-y)

def inv_trans(low, high):
    ylow = CDF(low)
    yhigh = CDF(high)
    y = random.random() * (yhigh - ylow) + ylow
    return CDF_inv(y)

def try_ar(low, high): 
    x = random.random() * 4
    y = random.random() #since max(prob(x)) over x = [-1, 5] is 1, range of 0 to 1 works
    return x, y, prob(x)

#Find a random number from e^(-x) distribution in range 0 to 5
# Use inverse-transform
x = inv_trans(0, 5)
print("Random number with inverse-transform: " + str(x))
print("Now generating " + str(NUM) + " numbers with accept-reject...")

X = np.empty(NUM)
x_accepted = 0
x_generated = 0
for i in range(0, NUM):
    x, y, pr = try_ar(0, 5)
    while y > pr:
        x_generated += 1
        x, y, pr = try_ar(0, 5)
    x_generated += 1
    x_accepted += 1
    X[i] = x

print("Generated; acceptance rate: " + str(x_accepted / x_generated))
plt.hist(X)
plt.show()
