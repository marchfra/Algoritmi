import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

plt.style.use(['grid', 'science', 'notebook', 'mylegend'])

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

df = pd.read_csv('data.csv')



fig, ax = plt.subplots(1, 1)
ax.plot(df['x'], df['pol'], label='horner')
ax.plot(df['x'], df['x']**2 - 2, label='pol')


# ax.set_xscale('log')
# ax.set_yscale('log')

ax.set_title('')
ax.set_xlabel('')
ax.set_ylabel('')

ax.legend()

fig.tight_layout()
# fig.savefig('.png', dpi=200)
plt.show()
