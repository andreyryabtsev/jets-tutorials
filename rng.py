import numpy as np
import matplotlib.pyplot as plt
import random

arr = np.empty(1000)
for i in range(1, 1000):
    arr[i-1] = random.random()
plt.hist(arr)
plt.title("Random uniform numbers")
plt.xlabel("Value")
plt.ylabel("Frequency")
plt.show()
