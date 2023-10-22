import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.style.use(['grid', 'science', 'notebook'])

df = pd.read_csv('random.csv')

fig, ax = plt.subplots(1, 1)
ax.scatter(np.arange(len(df['r'])), df['r'], s=1)

fig.tight_layout()
# for fmt in ['png', 'pdf']:
# 	fig.savefig(f'random.{fmt}', dpi=200)
# plt.show()

df = pd.read_csv('rand_err.csv')

fig, ax = plt.subplots(1, 1);
x = np.linspace(min(df['N']), max(df['N']))
ax.scatter(df['N'], df['k1'], label='k = 1', s=3)
ax.scatter(df['N'], df['k2'], label='k = 2', s=3)
ax.plot(x, 1 / np.sqrt(x), label='1/sqrt(N)')

ax.set_xlabel('N')

ax.legend(frameon=True, fancybox=True, framealpha=0.8)

ax.set_xscale('log')
ax.set_yscale('log')

fig.tight_layout()
fig.savefig('random.png', dpi=200)
plt.show()
