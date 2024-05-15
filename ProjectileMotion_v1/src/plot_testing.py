import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import colormaps as cmaps

# Matplotlib configuration
plt.style.use(['grid', 'science', 'notebook', 'mylegend'])
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']


# Import data
df = pd.read_csv('data/testing.csv')

df = df.groupby(by=df['param'])
groups = df.groups.keys()

fig, ax = plt.subplots(1, 1)
for g in groups:
	dfg = df.get_group(g).reset_index(drop=True)[::1]
	t = dfg['t'].to_numpy()
	x = dfg['x'].to_numpy()
	ax.plot(t, x - 1 * int(g), label='param' if int(g) ==
            1 else 'no param')  # , c=dfg.index.to_numpy(), marker='.')
# plt.set_cmap(plt.get_cmap('viridis'))
ax.scatter(1.24443, 1)
ax.axvline(1.24443, ls='--', c='k', alpha=.3)

ax.legend()

fig.tight_layout()
plt.show()
