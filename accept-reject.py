import numpy as np
import math, random

# use probability function N = exp(-(x - 5) ** 2 / 8) / math.sqrt(8 * math.pi)

def prob(x):
    return math.exp(-(x - 5) ** 2 / 8) / math.sqrt(8 * math.pi)

def get_rnd():
    y = 100
    x = 0
    x_pr = prob(x)
    while y > x_pr:
        x = random.random() * 10
        y = random.random()
        x_pr = prob(x)
    return x


#Find a Gaussian random number with mean = 5, sd = 2, range = [0, 10]
print(get_rnd())
