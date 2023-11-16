import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

plt.style.use(['grid', 'science', 'notebook', 'mylegend'])

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

df = pd.read_csv('data.csv')

dfgroup = df.groupby('s')
groups = dfgroup.groups.keys()

fig, ax = plt.subplots(1, 1)
for s in groups:
	ax.plot(dfgroup.get_group(s)['r'], dfgroup.get_group(s)['phi'], label=s)
ax.scatter(20, 1, marker='x', c='tab:red')

# ax.set_xscale('log')
# ax.set_yscale('log')

ax.set_title('')
ax.set_xlabel(r'$r$')
ax.set_ylabel(r'$\phi(r)$')

ax.legend(title=r"Initial guess for $\phi'(r = 0)$", ncols=2)

fig.tight_layout()


df = pd.read_csv('residual.csv')

fig2, ax2 = plt.subplots(1, 1)
ax2.plot(df['s'], df['Res'])

ax2.set_title('')
ax2.set_xlabel(r'$s$')
ax2.set_ylabel(r'$Res(s)$')

# ax2.legend()

fig2.tight_layout()



df = pd.read_csv('final.csv')

fig3, ax3 = plt.subplots(1, 1)
# ax3.plot(df['r'], df['phi'], label=r'$\phi(r)$')
ax3.plot(df['r'], df['Phi'], label=r'$\Phi(r)$')
ax3.plot(df['r'], df['Exact'], label=r'Exact $\Phi(r)$')

ax3.set_title('')
ax3.set_xlabel(r'$r$')
ax3.set_ylabel('')

ax3.legend()

fig3.tight_layout()

# fig.savefig('.png', dpi=200)
plt.show()
