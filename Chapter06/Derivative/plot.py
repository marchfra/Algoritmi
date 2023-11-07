import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

plt.style.use(['grid', 'science', 'notebook', 'mylegend'])

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

df = pd.read_csv('data.csv')



fig, ax = plt.subplots(1, 1)
ax.scatter(1 / df['h'], df['fd'], label='fd')
ax.scatter(1 / df['h'], df['bd'], label='bd')
ax.scatter(1 / df['h'], df['cd'], label='cd')
ax.scatter(1 / df['h'], df['hd'], label='hd')


ax.set_xscale('log')
ax.set_yscale('log')

ax.set_title('')
ax.set_xlabel('h')
ax.set_ylabel('err')

ax.legend()

fig.tight_layout()
# fig.savefig('.png', dpi=200)
plt.show()
