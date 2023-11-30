import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

plt.style.use(['grid', 'science', 'notebook', 'mylegend'])

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

df = pd.read_csv('data/data.csv')

dfgroup = df.groupby(by='method')

fig, ax = plt.subplots(1, 1)
ax.scatter(dfgroup.get_group('RK4')['x'], dfgroup.get_group('RK4')['y'], label='RK4')
ax.scatter(dfgroup.get_group('RK2')['x'], dfgroup.get_group('RK2')['y'], label='RK2')
ax.scatter(dfgroup.get_group('Euler')['x'][:20], dfgroup.get_group('Euler')['y'][:20], label='Euler')

# ax.set_xscale('log')
# ax.set_yscale('log')

ax.set_title('')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_aspect(1)

ax.legend()

fig.tight_layout()
# fig.savefig('fig1.png', dpi=200)

fig2, ax2 = plt.subplots(2, 1, sharex=True, figsize=(8,8))
ax2[0].plot(dfgroup.get_group('RK4')['t'], dfgroup.get_group('RK4')['x'], label='RK4')
# ax2[0].plot(dfgroup.get_group('RK2')['t'], dfgroup.get_group('RK2')['x'], label='RK2')
# ax2[0].plot(dfgroup.get_group('Euler')['t'][:20], dfgroup.get_group('Euler')['x'][:20], label='Euler')
ax2[1].plot(dfgroup.get_group('RK4')['t'], dfgroup.get_group('RK4')['y'], label='RK4')
# ax2[1].plot(dfgroup.get_group('RK2')['t'], dfgroup.get_group('RK2')['y'], label='RK2')
# ax2[1].plot(dfgroup.get_group('Euler')['t'][:20], dfgroup.get_group('Euler')['y'][:20], label='Euler')

# ax2.set_xscale('log')
# ax2.set_yscale('log')

# ax2.set_title('')
# ax2[0].set_xlabel('t')
ax2[1].set_xlabel('t')
ax2[0].set_ylabel('x')
ax2[1].set_ylabel('y')

ax2[0].legend()
ax2[1].legend()

fig2.tight_layout()
# fig2.savefig('fig2.png', dpi=200)

df = pd.read_csv('data/convergence.csv')

dfgroup = df.groupby(by='method')

fig3, ax3 = plt.subplots(1, 1)
rk4 = ax3.scatter(dfgroup.get_group('RK4')['dt'], dfgroup.get_group('RK4')['err'], label='RK4')
rk2 = ax3.scatter(dfgroup.get_group('RK2')['dt'], dfgroup.get_group('RK2')['err'], label='RK2')
eul = ax3.scatter(dfgroup.get_group('Euler')['dt'][:20], dfgroup.get_group('Euler')['err'][:20], label='Euler')

ax3.set_xscale('log')
ax3.set_yscale('log')

ax3.set_title('Convergence test')
ax3.set_xlabel('dt')
ax3.set_ylabel('err')

ax3.legend(handles=[eul, rk2, rk4])

fig3.tight_layout()
# fig3.savefig('convergence.png', dpi=200)
plt.show()
