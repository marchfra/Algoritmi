import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

# Constants
IMAGES_FOLDER = 'images'

# Matplotlib configuration
plt.style.use(['grid', 'science', 'notebook', 'mylegend'])
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

# Import data
df = pd.read_csv('data/data.csv')

# Manipulate data
dfRK = df[['t'] + [f'{var}RK' for var in ['th1', 'th2', 'w1', 'w2']]]
dfPV = df[['t'] + [f'{var}PV' for var in ['th1', 'th2', 'w1', 'w2']]]

dfRK.columns = [f'{var}' for var in ['t', 'th1', 'th2', 'w1', 'w2']]
dfPV.columns = [f'{var}' for var in ['t', 'th1', 'th2', 'w1', 'w2']]

# print(f'Full df:\n{df.head(10)}\n')
# print(f'RK:\n{dfRK.head(10)}\n')
# print(f'PV:\n{dfPV.head(10)}\n')

# Create figure
fig, ax = plt.subplots(2, 1, sharex=True)
mark = 'o'
ax[0].plot(dfRK['t'], dfRK['th1'], label=r'$\theta_1^{RK}$')
ax[1].plot(dfRK['t'], dfRK['th2'], label=r'$\theta_2^{RK}$')
ax[0].plot(dfPV['t'], dfPV['th1'], label=r'$\theta_1^{PV}$')
ax[1].plot(dfPV['t'], dfPV['th2'], label=r'$\theta_2^{PV}$')

# ax[0].axhline(np.pi / 4)
# ax[0].axhline(-np.pi / 4)
# ax[1].axhline(np.pi / 4)
# ax[1].axhline(-np.pi / 4)

ax[0].set_title('')
ax[1].set_xlabel(r'$t$')
ax[0].set_ylabel(r'$\theta_1$')
ax[1].set_ylabel(r'$\theta_2$')

ax[0].legend()
ax[1].legend()

fig.tight_layout()


fig, ax = plt.subplots(1, 1)
ax.plot(dfRK['t'], dfRK['th1'], label=r'$\theta_1^{RK}$', c=colors[0])
ax.plot(dfRK['t'], dfRK['th2'], label=r'$\theta_2^{RK}$', c=colors[2])
ax.plot(dfPV['t'], dfPV['th1'], label=r'$\theta_1^{PV}$', c=colors[1])
ax.plot(dfPV['t'], dfPV['th2'], label=r'$\theta_2^{PV}$', c=colors[3])

ax.set_xlim(0, 50)

ax.set_title('Coupled pendula')
ax.set_xlabel(r'$t$')
ax.set_ylabel(r'$\theta$')

ax.legend(ncols=2)

fig.tight_layout()
fig.savefig(f'{IMAGES_FOLDER}/marchisotti_coupled_pendula.png', dpi=200)
plt.show()
