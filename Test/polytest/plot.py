import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

plt.style.use(['grid', 'science', 'notebook', 'mylegend'])

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

df = pd.read_csv('data.csv')

# df['p'] = df['x']**3 + 5 * df['x']**2 - 3 * df['x']

fig, ax = plt.subplots(1, 1)
ax.plot(df['x'], df['pol'], label='pol')
ax.plot(df['x'], df['horn'], label='horn')
# ax.plot(df['x'], df['p'], label='p')


# ax.set_xscale('log')
# ax.set_yscale('log')

ax.set_title('')
ax.set_xlabel('')
ax.set_ylabel('')

ax.legend()

fig.tight_layout()
fig.savefig('test1.png', dpi=200)
# plt.show()
