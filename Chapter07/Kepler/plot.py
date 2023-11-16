import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

plt.style.use(['grid', 'science', 'notebook', 'mylegend'])

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

# df = pd.read_csv('ellipse_test.csv')
df = pd.read_csv('N_orbits.csv')

E0 = df['E'][0]
L0 = df['L'][0]
df['dE [%]'] = np.abs((df['E'] - E0) / E0) * 100
df['dL [%]'] = np.abs((df['L'] - L0) / L0) * 100

print(df.describe())

fig, ax = plt.subplots(1, 1)
ax.scatter(df['x'], df['y'], s=5)

ax.set_aspect('equal')
# ax.set_xscale('log')
# ax.set_yscale('log')

ax.set_title('Ellipse test')
ax.set_xlabel('x')
ax.set_ylabel('y')

# ax.legend()

fig.tight_layout()
# fig.savefig('.png', dpi=200)
plt.show()
