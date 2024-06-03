import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

list = []

for k in tqdm(range(10)):
    rng = np.random.default_rng(seed = k)

    x = 2 # x0
    counter = 0

    ball_size = 4

    while (x**2 + 3 * x - 37) > 0.01 or (x**2 + 3 * x - 37) < -0.01:
        r = (2 * ball_size) * rng.random() - ball_size # random number between -ball_size and ball_size, uniform distribution
        while x + r < 0 or x + r > 10:
            r = (2 * ball_size) * rng.random() - ball_size # random number between -ball_size and ball_size, uniform distribution
        x = x + r
        counter += 1

    list.append(counter)

plt.stem(list)
plt.xlabel('RNG Seed')
plt.ylabel('Number of Iterations')
plt.title('Number of Iterations vs. Seed, Mean = ' + str(np.mean(list)))
plt.show()

print('debug')
