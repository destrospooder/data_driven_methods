import numpy as np

rng = np.random.default_rng(seed = 488)

x = 2 # x0
counter = 0

while (x**2 + 3 * x - 37) > 0.01 or (x**2 + 3 * x - 37) < -0.01:
    r = (1 - (-1)) * rng.random() - 1
    while x + r < 0 or x + r > 10:
        r = (1 - (-1)) * rng.random() - 1
    x = x + r
    print(x)
    counter += 1

print('Root: ' + str(x))
print('There were ' + str(counter) + ' iterations.')