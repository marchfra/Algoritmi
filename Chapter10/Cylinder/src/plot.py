import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
# from matplotlib.ticker import FormatStrFormatter

# Constants
IMAGES_FOLDER = 'images'

# Matplotlib configuration
plt.style.use(['grid', 'science', 'notebook', 'mylegend'])
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

# Import data
# df = pd.read_csv('data/data_jacobi.csv')
# df = pd.read_csv('data/data_gauss.csv')
df = pd.read_csv('data/data_sor.csv')

# Manipulate data
x = df['x'].unique()
y = df['y'].unique()
z = np.reshape(df['M'].to_numpy(), (len(x), len(y))).T

# Create figure
fig, ax = plt.subplots(1, 1)
data = ax.contourf(x, y, z, levels=500, cmap=plt.cm.RdYlBu.reversed())
contour_lines = ax.contour(x, y, z, levels=data.levels[::50], colors='k', alpha=0.5, linewidths=1)

cbar = fig.colorbar(data)
cbar.ax.set_ylabel(r'$\phi(x, y)$')
cbar.add_lines(contour_lines)

cbar.ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.2g'))
cbar.ax.invert_yaxis()

# ax.set_xscale('log')
# ax.set_yscale('log')

ax.set_title('')
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')

ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))

ax.set_aspect(1)

# ax.legend()

fig.tight_layout()
fig.savefig(f'{IMAGES_FOLDER}/cylinder_sor.png', dpi=200)
plt.show()
