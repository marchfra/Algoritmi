import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Constants
IMAGES_FOLDER = 'images'

# Matplotlib configuration
plt.style.use(['grid', 'science', 'notebook', 'mylegend'])
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

# Import data
df = pd.read_csv('data/data.csv')

# Manipulate data
dfgroup = df.groupby(by=df['theta'])
groups = dfgroup.groups.keys()

# Create figure
fig, ax = plt.subplots(1, 1)
for i, g in enumerate(groups):
	x = dfgroup.get_group(g)['x'].to_numpy()
	y = dfgroup.get_group(g)['y'].to_numpy()
	ax.plot(x, y, color=colors[i], label=f'{g * 180 / np.pi:.0f}°')
# x = dfgroup.get_group(1)['x'].to_numpy()
# y = dfgroup.get_group(1)['y'].to_numpy()
# ax.plot(x, y)


# ax.set_xscale('log')
# ax.set_yscale('log')

ax.set_title('')
ax.set_xlabel('x [adim]')
ax.set_ylabel('y [adim]')

ax.legend(title='Initial angle')

fig.tight_layout()
# fig.savefig(f'{IMAGES_FOLDER}/figure_name.png', dpi=200)
plt.show()
