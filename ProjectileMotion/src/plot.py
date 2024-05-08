import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from ticker import multiple_formatter

# Constants
IMAGES_FOLDER = 'images'
SAVE_IMAGES = False
DOWNSAMPLE = 1
DRAW_LEGENDS = False

# Matplotlib configuration
plt.style.use(['grid', 'science', 'notebook', 'mylegend'])
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']


# +---------------------------+
# | STEP 0: GETTING CONSTANTS |
# +---------------------------+

df = pd.read_csv('data/constants.csv')
L = float(df['chi'][0])
tau = float(df['tau'][0])
mu = float(df['mu'][0])
B = float(df['B'][0])
b = float(df['b'][0])

print(f'+------------------------------+')
print(f'|           CONSTANTS          |')
print(f'+------------------------------+')
print(f'|         L = {L:<7.1f} m        |')
print(f'|       tau = {tau:<7.5f} s        |')
print(f'|        mu = {mu:<7.2f} kg       |')
print(f'|         B = {B:<7.1e} kg/m     |')
print(f'+------------------------------+')


# +---------------------------------+
# | STEP 1: INTEGRATING UP TO x = L |
# +---------------------------------+

# Import data
df = pd.read_csv('data/shooting.csv')

# Manipulate data
dfgroup = df.groupby(by=df['theta'])
groups = dfgroup.groups.keys()
groups = [list(groups)[0]] + [list(groups)[-1]]

# Create figure
fig, ax = plt.subplots(1, 1)
for i, g in enumerate(list(groups)[::DOWNSAMPLE]):
	x = dfgroup.get_group(g)['x'].to_numpy()
	y = dfgroup.get_group(g)['y'].to_numpy()
	ax.plot(x, y, c=colors[i % len(colors)], label=f'{g} rad')
ax.axvline(L, c='k', ls='--', alpha=.7)

ax.set_title(rf'Shooting up to $x = ({L} + \epsilon)$ m, $B = {B}$ kg/m')
ax.set_xlabel(r'$x$ [m]')
ax.set_ylabel(r'$y$ [m]')

ax.legend(title='Initial angle', fontsize=15)

if DRAW_LEGENDS:
	ax.legend(title='Initial angle', ncol=int(
		len(groups) / DOWNSAMPLE / 4), fontsize=15)

fig.tight_layout()
if SAVE_IMAGES:
	fig.savefig(f'{IMAGES_FOLDER}/Initial shooting.png', dpi=200)


# +-----------------------------------+
# | STEP 2: FINDING ROOTS OF x(t) = L |
# +-----------------------------------+ */

# Import data
t_star = pd.read_csv('data/x_t.csv')['t*'].to_numpy()


# Create figure
fig, ax = plt.subplots(1, 1)
for i, (g, T) in enumerate(zip(list(groups)[::DOWNSAMPLE], t_star[::DOWNSAMPLE])):
	t = dfgroup.get_group(g)['t'].to_numpy()
	x = dfgroup.get_group(g)['x'].to_numpy()
	ax.plot(t, x, c=colors[i % len(colors)], label=f'{g:.2f} rad')
	ax.axvline(T, ymax=.87, c=colors[i % len(colors)], ls='--', alpha=.3)
ax.axhline(L, c='k', ls='--', alpha=.7)

ax.set_title(rf'Finding $t^*$ s.t. $x(t^*) = {L}$ m, $B = {B}$ kg/m')
ax.set_xlabel(r'$t$ [s]')
ax.set_ylabel(r'$x$ [m]')

if DRAW_LEGENDS:
	ax.legend(title='Initial angle', ncol=int(
		len(groups) / DOWNSAMPLE / 4), fontsize=15)

fig.tight_layout()
if SAVE_IMAGES:
	fig.savefig(f'{IMAGES_FOLDER}/Finding t_end.png', dpi=200)


# +--------------------------------+
# | STEP 2.5: SHOOTING UP TO x = L |
# +--------------------------------+

# Import data
df2 = pd.read_csv('data/shooting2.csv')

# Manipulate data
df2group = df2.groupby(by=df2['theta'])
groups2 = df2group.groups.keys()
groups2 = [list(groups2)[0]] + [list(groups2)[-1]]

# Create figure
fig, ax = plt.subplots(1, 1)
for i, g in enumerate(list(groups2)[::DOWNSAMPLE]):
	x = df2group.get_group(g)['x'].to_numpy()
	y = df2group.get_group(g)['y'].to_numpy()
	ax.plot(x, y, c=colors[i % len(colors)], label=f'{g:.2f} rad')
ax.axvline(L, c='k', ls='--', alpha=.7)

ax.set_title(rf'Shooting up to $x = {L}$ m, $B = {B}$ kg/m')
ax.set_xlabel(r'$x$ [m]')
ax.set_ylabel(r'$y$ [m]')

if DRAW_LEGENDS:
	ax.legend(title='Initial angle', ncol=int(
		len(groups2) / DOWNSAMPLE / 4), fontsize=15)

fig.tight_layout()
if SAVE_IMAGES:
	fig.savefig(f'{IMAGES_FOLDER}/Second shooting.png', dpi=200)


# +----------------------------+
# | STEP 3: RESIDUAL OF y(x=L) |
# +----------------------------+ */

# Import data
df_res = pd.read_csv('data/residual.csv')

# Create figure
fig, ax = plt.subplots(1, 1)
ax.plot(df_res['theta'], df_res['y*'])
ax.axhline(0, c='k', ls='--', alpha=0.7)

# def f(x):
# 	return -20 * (x - 0.787)**2 + 0.19

# ax.plot(df_res['theta'], f(df_res['theta']))

df_res = pd.read_csv('data/optimal_search.csv')
# print(df_res['theta'])
# df_res = df_res.sort_values(by=['theta'])
_theta = df_res['theta'].to_numpy()
_y = df_res['y*'].to_numpy()
ax.quiver(_theta[:-1], _y[:-1], _theta[1:] - _theta[:-1],
          _y[1:] - _y[:-1], scale_units='xy', angles='xy')
ax.plot(_theta, _y)

ax.set_title(rf'Residual plot of $y(x={L}$m$)$, $B = {B}$ kg/m')
ax.set_xlabel(r'$\theta$ [rad]')
ax.set_ylabel(r'$y(x=1)$ [m]')

# ax.xaxis.set_major_locator(plt.MultipleLocator(np.pi / 78))
# ax.xaxis.set_minor_locator(plt.MultipleLocator(np.pi / (78 * 4)))
# ax.xaxis.set_major_formatter(plt.FuncFormatter(multiple_formatter(78)))

fig.tight_layout()
if SAVE_IMAGES:
	fig.savefig(f'{IMAGES_FOLDER}/Residual.png', dpi=200)


# +---------------------------------+
# | STEP 5: INTEGRATING WITH THETA* |
# +---------------------------------+ */

# Import data
df = pd.read_csv('data/final_trajectories.csv')
df = df.groupby(by=df['theta'])
groups = df.groups.keys()

# Create figure
fig, ax = plt.subplots(1, 1)
for g in groups:
	x = df.get_group(g)['x']
	y = df.get_group(g)['y']
	ax.plot(x, y, label=f'{g:.2f} rad')
ax.axvline(L, c='k', ls='--', alpha=0.7)

ax.set_title(rf'Optimal trajectory, $B = {B}$ kg/m')
ax.set_xlabel(r'$x$ [m]')
ax.set_ylabel(r'$y$ [m]')

ax.legend(title='Launch angle')

fig.tight_layout()
if SAVE_IMAGES:
	fig.savefig(f'{IMAGES_FOLDER}/Optimal trajectories.png', dpi=200)


plt.show()
