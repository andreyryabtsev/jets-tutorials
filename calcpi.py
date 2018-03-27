import numpy as np
import matplotlib.pyplot as plt
import random, sys

def find_pi(n):
    arr = np.empty(n)
    for i in range (0, n):
        arr[i] = np.sqrt((random.random() * 2 - 1) ** 2 + (random.random() * 2 - 1) ** 2) < 1
    return np.sum(arr) / n * 4


def several_pi(n, points = 1000):
    results = np.empty(n)
    for i in range (0, n):
        results[i] = find_pi(points)
    return np.std(results)

for i in range(1, 12):
    n = 2 ** i
    print(n, "calculations with 1000 points: sigma =", several_pi(n))


for i in range(2, 6):
    pts = 10 ** i
    print("10 calculations with", pts, "points: sigma =", several_pi(10, pts))
