import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

plt.style.use(['grid', 'science', 'notebook', 'mylegend'])

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

df = pd.read_csv('data.csv')

df = df.groupby('k')
groups = df.groups.keys()

fig, ax = plt.subplots(1, 1)
for k in groups:
	ax.plot(df.get_group(k)['x'], df.get_group(k)['phi'], label=f'{k}')

ax.set_title('Shooting')
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$\phi(x)$')

ax.legend(title=r'Value of $k$', ncols=2)

fig.tight_layout()
# fig.savefig('.png', dpi=200)


df = pd.read_csv('final.csv')
df = df.groupby('k')
groups = df.groups.keys()

fig, ax = plt.subplots(1, 1)
for k in groups:
	ax.plot(df.get_group(k)['x'], df.get_group(k)['phi'], label=f'{k/np.pi:.0f}π')

ax.set_title('Final result')
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$\phi(x)$')

ax.legend(title=r'Value of $k$', ncols=2)

fig.tight_layout()
# fig.savefig('.png', dpi=200)

plt.show()
