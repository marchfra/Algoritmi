import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

plt.style.use(['grid', 'science', 'notebook', 'mylegend'])

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

df = pd.read_csv('data.csv')

dfgroup = df.groupby(by='h')

# print(dfgroup.get_group(0.5))

def exact(t):
	return np.exp(-0.5 * t**2)

fig, ax = plt.subplots(1, 2, figsize=(16, 6), sharey=True)
for h in list(dfgroup.groups.keys())[::-1]:
	ax[0].plot(dfgroup.get_group(h)['t'][:-1], dfgroup.get_group(h)['y'][:-1], label=f'h = {h:.3f}')
	ax[1].plot(dfgroup.get_group(h)['t'][:-1], dfgroup.get_group(h)['rel_err'][:-1], label=f'h = {h}')
ax[0].plot(np.linspace(0, 3), exact(np.linspace(0, 3)), label='exact')


# ax.set_xscale('log')
# ax.set_yscale('log')

ax[0].set_title('y(t)')
ax[1].set_title('rel_err(t)')
ax[0].set_xlabel('t')
ax[1].set_xlabel('t')
# ax.set_ylabel('')

ax[0].legend()

fig.tight_layout()
# fig.savefig('.png', dpi=200)
plt.show()
