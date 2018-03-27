import numpy as np
import matplotlib.pyplot as plt

def gen_gauss(mean, sd, low, high):
    n = low - 1
    while n < low or n > high:
        n = np.random.normal(mean, sd)
    return n


print(gen_gauss(5, 2, 0, 10))
