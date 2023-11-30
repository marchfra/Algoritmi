import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

# Constants
IMAGES_FOLDER = 'images'

# Matplotlib configuration
plt.style.use(['grid', 'science', 'notebook', 'mylegend'])
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

# Import data
df = pd.read_csv('data/data.csv')

# Manipulate data

def f(x):
	c1 = -0.5
	c2 = 0
	return 0.5 * x**2 + c1 * x + c2

# Create figure
fig, ax = plt.subplots(1, 1)
ax.scatter(df['x'], df['y'])

ax.set_title(r'')
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')

# ax.legend()

fig.tight_layout()
# fig.savefig(f'{IMAGES_FOLDER}/bvp2.png', dpi=200)
plt.show()
