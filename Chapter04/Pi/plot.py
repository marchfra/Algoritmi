import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

plt.style.use(['grid', 'science', 'notebook', 'mylegend'])

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

df = pd.read_csv('data.csv')

xdense = np.linspace(min(df['N']), max(df['N']))
ydense = 1 / np.sqrt(xdense)

fig, ax = plt.subplots(1, 1)
ax.plot(xdense, ydense, color=colors[0], label=r'$\frac{1}{\sqrt{N}}$')
ax.scatter(df['N'], df['err'], color=colors[1], label=r'err$(N)$')
ax.scatter(df['N'], df['sigma'], color=colors[2], label=r'$\sigma(N)$')

ax.set_xscale('log')
ax.set_yscale('log')

ax.set_title(r'Error in estimating $\pi$ with Monte Carlo methods')
ax.set_xlabel('N')
ax.set_ylabel('')

ax.legend()

fig.tight_layout()
fig.savefig('error.png', dpi=200)
# plt.show()
