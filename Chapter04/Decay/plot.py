import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.style.use(['grid', 'science', 'notebook'])

df = pd.read_csv('data.csv')

df['teo'] = np.exp(-df['lambda'] * df['t'])

fig, ax = plt.subplots(1, 1)
ax.plot(df['t'], df['teo'], color='tab:green', label=r'$e^{-\lambda t}$')
ax.scatter(df['t'], df['N_norm'], marker='.', label='simulation', zorder=10)

ax.set_yscale('log')

ax.set_title('Radioactive decay')
ax.set_xlabel(r'$t$')
ax.set_ylabel(r'$N/N_0$')

ax.legend(frameon=True, fancybox=True, framealpha=0.8)

fig.tight_layout()
fig.savefig('rad_decay.png', dpi=200)
# plt.show()
