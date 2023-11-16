import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

plt.style.use(['grid', 'science', 'notebook', 'mylegend'])

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

df = pd.read_csv('data.csv')

dfgroup = df.groupby(by='method')
methods = dfgroup.groups.keys()

xdense = np.linspace(min(dfgroup.get_group('RK2')['t']), max(dfgroup.get_group('RK2')['t']))

fig, ax = plt.subplots(1, 1)
for m in methods:
	print(dfgroup.get_group(m).head())
	ax.plot(dfgroup.get_group(m)['t'], dfgroup.get_group(m)['E'], label=m)
# ax.plot(xdense, xdense, label='retta')

# ax.set_xscale('log')
# ax.set_yscale('log')

ax.set_title('')
ax.set_xlabel('t')
ax.set_ylabel('E')

ax.legend()

fig.tight_layout()
# fig.savefig('.png', dpi=200)
plt.show()
