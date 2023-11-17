import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

plt.style.use(['grid', 'science', 'notebook', 'mylegend'])

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

dfF = pd.read_csv('forward.csv')
dfB = pd.read_csv('backward.csv')

fig, ax = plt.subplots(1, 1)
ax.plot(dfF['x'][:len(dfF)//2], dfF['psi'][:len(dfF)//2], label='Forward Integration')
ax.plot(dfB['x'][:len(dfF)//2], dfB['psi'][:len(dfF)//2], label='Backward Integration')

ax.set_yscale('log')

ax.set_title('Usual naive integration')
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$\psi(x)$')

ax.legend()

fig.tight_layout()


df = pd.read_csv('matching_shoot.csv')
df = df.groupby('E')
groups = df.groups.keys()

fig, ax = plt.subplots(1, 1)
for i, E in enumerate(groups):
	ax.plot(df.get_group(E)['x'][:800], df.get_group(E)['psi'][:800], c=colors[i], label=f'{E}')
	ax.plot(df.get_group(E)['x'][800:], df.get_group(E)['psi'][800:], c=colors[i])

ax.set_yscale('log')

ax.set_title('Matching Point Shooting')
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$\psi(x)$')

ax.legend()

fig.tight_layout()


df = pd.read_csv('res_plot.csv')

fig, ax = plt.subplots(1, 1)
ax.plot(df['E'], df['Res'])

ax.set_title('Residual')
ax.set_xlabel(r'$E$')
ax.set_ylabel(r'$Res(E)$')

fig.tight_layout()


df = pd.read_csv('final.csv')

fig, ax = plt.subplots(1, 1)
ax.plot(df['x'][:800], df['psi'][:800], c=colors[0])
ax.plot(df['x'][800:], df['psi'][800:], c=colors[0])

# ax.set_yscale('log')

ax.set_title('Final result')
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$\psi(x)$')

# ax.legend()

fig.tight_layout()
# fig.savefig('.png', dpi=200)plt.show()
plt.show()
