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
df = df.groupby(df['s'])
groups = df.groups.keys()

# Create figure
fig, ax = plt.subplots(1, 1)
for s in groups:
	ax.plot(df.get_group(s)['x'], df.get_group(s)['f'], label=f'{s=}')


# ax.set_xscale('log')
# ax.set_yscale('log')

ax.set_title('')
ax.set_xlabel('E')
ax.set_ylabel('')

ax.legend()

fig.tight_layout()
fig.savefig(f'{IMAGES_FOLDER}/energy_levels.png', dpi=200)
plt.show()
