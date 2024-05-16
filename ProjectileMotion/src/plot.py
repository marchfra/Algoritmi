import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from helper_functions import *

# Constants
IMAGES_FOLDER = 'images'
DATA_FOLDER = 'data'
SAVE_IMAGES = False
DOWNSAMPLE = 4

# Matplotlib configuration
plt.style.use(['grid', 'science', 'notebook', 'mylegend'])
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']


# +---------------------------+
# | STEP 0: GETTING CONSTANTS |
# +---------------------------+

constants = pd.read_csv('data/constants.csv')
L = float(constants['chi'][0])
tau = float(constants['tau'][0])
mu = float(constants['mu'][0])
B = float(constants['B'][0])
b = float(constants['b'][0])
V0 = float(constants['V0'][0])
h = float(constants['Y_targ'][0])


def savefig(fig: mpl.figure.Figure, figure_name: str) -> None:
	if SAVE_IMAGES:
		fig.savefig(f'{IMAGES_FOLDER}/{figure_name}.png', dpi=200)


def draw_target(ax: mpl.axes._axes.Axes, L: float, h: float) -> None:
	'''
	This function draws the target on the supplied axes.
	'''
	scale: float = 1.5
	ax.scatter(L, h, c='black', marker='o', s=scale * 350)
	ax.scatter(L, h, c='white', marker='o', s=scale * 340)
	ax.scatter(L, h, c='black', marker='o', s=scale * 200)
	ax.scatter(L, h, c='dodgerblue', marker='o', s=scale * 100)
	ax.scatter(L, h, c='red', marker='o', s=scale * 40)
	ax.scatter(L, h, c='yellow', marker='o', s=scale * 10)
	ax.scatter(L, h, c='black', marker='.', s=scale * 0.5)


def print_constants() -> None:
	print(f'+------------------------------+')
	print(f'|           CONSTANTS          |')
	print(f'+------------------------------+')
	print(f'|         L = {L:<7.1f} m        |')
	print(f'|         h = {h:<7.1f} m        |')
	print(f'|        V0 = {V0:<7.2f} m/s      |')
	print(f'|       tau = {tau:<7.5f} s        |')
	print(f'|        mu = {mu:<7.2f} kg       |')
	print(f'|         B = {B:<7.1e} kg/m     |')
	print(f'+------------------------------+')


def shooting_plot() -> None:
	# Import data
	df = pd.read_csv(f'{DATA_FOLDER}/shooting.csv')

	# Manipulate data
	df = df.groupby(by=df['theta'])
	groups = list(df.groups.keys())[::DOWNSAMPLE]

	# Create figure
	fig, ax = plt.subplots(1, 1)
	for g in groups:
		curr_df = df.get_group(g)
		x = curr_df['x'].to_numpy() * L
		y = curr_df['y'].to_numpy() * L
		ax.plot(x, y, label=f'{g:.2f} rad')
	draw_target(ax, L, h)

	ax.set_title('Shooting')
	ax.set_xlabel(r'$x$ [m]')
	ax.set_ylabel(r'$y$ [m]')

	ax.legend(title='Launch angle', ncol=2)

	fig.tight_layout()
	savefig(fig, 'shooting')


def residual_plot() -> None:
	# Import data
	df = pd.read_csv(f'{DATA_FOLDER}/residual.csv')

	# Create figure
	fig, ax = plt.subplots(1, 1)
	x = df['theta'].to_numpy()
	res = df['res'].to_numpy() * L
	ax.plot(x, res)
	ax.axhline(0, c='k', ls='--', alpha=0.7)

	ax.set_title('Residual')
	ax.set_xlabel(r'$\theta$ [rad]')
	ax.set_ylabel(r'$\Delta y$ [m]')

	fig.tight_layout()
	savefig(fig, 'shooting')


def final_plot() -> None:
	# Import data
	df = pd.read_csv(f'{DATA_FOLDER}/finalTrajectories.csv')

	# Manipulate data
	df = df.groupby(by=df['theta'])
	groups = list(df.groups.keys())

	# Create figure
	fig, ax = plt.subplots(1, 1)
	for g in groups:
		curr_df = df.get_group(g)
		x = curr_df['x'].to_numpy() * L
		y = curr_df['y'].to_numpy() * L
		ax.plot(x, y, label=f'{g:.2f} rad')
	draw_target(ax, L, h)

	ax.set_title('Final trajectories')
	ax.set_xlabel(r'$x$ [m]')
	ax.set_ylabel(r'$y$ [m]')

	ax.legend(title='Launch angle')

	fig.tight_layout()
	savefig(fig, 'final')


def analytical(x: np.ndarray, sol_number: int) -> float:
	if sol_number not in range(2):
		raise ValueError('sol_number must be either 0 or 1')

	x = x / L
	v0 = 10.0 * tau / L
	# print(f'{v0=}')
	theta = 0.5 * (sol_number * np.pi + (-1)**sol_number * np.arcsin(1 / v0**2))
	u0 = v0 * np.cos(theta)
	v0 = v0 * np.sin(theta)
	return L * (-0.5 * (x / u0)**2 + v0 / u0 * x)


def comparison_plot() -> None:
	# Import data
	df = pd.read_csv(f'{DATA_FOLDER}/noFriction.csv')

	# Manipulate data
	df = df.groupby(by=df['theta'])
	groups = list(df.groups.keys())

	# Create figure
	fig, ax = plt.subplots(1, 1)
	for i, g in enumerate(groups):
		curr_df = df.get_group(g)
		x = curr_df['x'].to_numpy() * L
		y = curr_df['y'].to_numpy() * L
		exact = analytical(x, i)
		ax.plot(x, abs(y - exact), label=f'{g:.2f} rad')
		# ax.scatter(x, exact, label=f'{g:.2f} rad')
		# ax.scatter(x, y, label=f'{g:.2f} rad')
	draw_target(ax, 10.0, 0.0)

	ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

	ax.set_title('Comparison with analytical solution')
	ax.set_xlabel(r'$x$ [m]')
	ax.set_ylabel(r'$\Delta y$ [m]')

	ax.legend(title='Launch angle')

	fig.tight_layout()
	savefig(fig, 'comparison')


def main() -> None:
	# print_constants()
	# shooting_plot()
	# residual_plot()
	# final_plot()
	comparison_plot()
	plt.show()


if __name__ == '__main__':
	main()
