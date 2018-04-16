from modelP import random_falling

x = 0
for i in range(1, 10000):
    r = random_falling()
    x += r

print(x/10000)
