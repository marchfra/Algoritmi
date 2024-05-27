import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Constants
IMAGES_FOLDER = 'images'
DATA_FOLDER = 'data'
SAVE_IMAGES = False
DOWNSAMPLE = 4

# Matplotlib configuration
plt.style.use(['grid', 'science', 'notebook', 'mylegend'])
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']


# Getting constants
constants = pd.read_csv('data/constants.csv')
L = float(constants['chi'][0])
tau = float(constants['tau'][0])
mu = float(constants['mu'][0])
B = float(constants['B'][0])
b = float(constants['b'][0])
V0 = float(constants['V0'][0])
h = float(constants['YTarg'][0])


def savefig(fig, figure_name: str) -> None:
	if SAVE_IMAGES:
		fig.savefig(f'{IMAGES_FOLDER}/{figure_name}.png', dpi=200)


def draw_target(ax, L: float, h: float) -> None:
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

	ax.legend(title='Launch angle', ncol=len(groups) / 4)

	fig.tight_layout()
	savefig(fig, 'shooting')


def residual_plot() -> None:
	# Import data
	df = pd.read_csv(f'{DATA_FOLDER}/residual.csv')

	thetas = pd.read_csv(
		f'{DATA_FOLDER}/noFriction.csv')['theta'].to_numpy()

	# Create figure
	fig, ax = plt.subplots(1, 1)
	theta = df['theta'].to_numpy()
	res = df['res'].to_numpy() * L
	ax.plot(theta, res, marker='o')
	ax.axhline(0, c='k', ls='--', alpha=0.7)

	ax.axvline(0.6877629863480418, c='tab:orange',
	           ls='--', label='analytical zero')
	ax.axvline(0.8830333404468548, c='tab:orange',
	           ls='--', label='analytical zero')
	for theta in thetas:
		ax.axvline(theta, c='tab:red', ls='--', label='numerical zero')

	ax.set_title('Residual')
	ax.set_xlabel(r'$\theta$ [rad]')
	ax.set_ylabel(r'$\Delta y$ [m]')

	fig.tight_layout()
	savefig(fig, 'residual')


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


def analytical(x: np.ndarray | float, sol_number: int) -> np.ndarray | float:
	if sol_number not in range(2):
		raise ValueError('sol_number must be either 0 or 1')

	x = x / L
	v0 = 10.0 * tau / L
	# print(f'{v0=}')
	theta = 0.5 * (sol_number * np.pi + (-1)**sol_number * np.arcsin(1 / v0**2))
	print(f'Ana theta{sol_number + 1} = {theta:.7f}')
	# theta = np.pi / 4
	u0 = v0 * np.cos(theta)
	v0 = v0 * np.sin(theta)
	return (-0.5 * (x / u0)**2 + v0 / u0 * x) * L


def comparison_plot() -> None:
	# Import data
	df = pd.read_csv(f'{DATA_FOLDER}/noFriction.csv')

	# Manipulate data
	df = df.groupby(by=df['theta'])
	groups = list(df.groups.keys())

	# Create figure
	fig, ax = plt.subplots(1, 1)
	for i, g in enumerate(groups):
		# if i == 1:
		# 	break
		print(f'\nNum theta{i + 1} = {g:.7f}')
		curr_df = df.get_group(g)
		x = curr_df['x'].to_numpy() * L
		y = curr_df['y'].to_numpy() * L
		exact = analytical(x, i)
		ax.plot(x[:-1], abs(y - exact)[:-1], label=f'{g:.2f} rad')
		# ax.plot(x, y, label='num', marker='o')
		# x = np.linspace(9.9, 10.1, 1000)
		# exact = analytical(x, i)
		# ax.plot(x, exact, label='exact')
	# ax.scatter(10.0, 0.0, c='tab:red')
	draw_target(ax, 10.0, 0.0)

	ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

	ax.set_title('Comparison with analytical solution')
	ax.set_xlabel(r'$x$ [m]')
	ax.set_ylabel(r'$\Delta y$ [m]')

	ax.legend(title='Launch angle', ncol=2)

	fig.tight_layout()
	savefig(fig, 'comparison')


def linear_plot() -> None:
	# Import data
	df = pd.read_csv(f'{DATA_FOLDER}/linear.csv')

	# Manipulate data
	x = df['x'].to_numpy()
	y = df['y'].to_numpy()

	# Create figure
	fig, ax = plt.subplots(1, 1)
	ax.plot(x, y)
	ax.scatter(-3, 4, c='k')
	ax.scatter(2, 1.5, c='k')

	ax.set_title('Linear interpolation test')
	ax.set_xlabel(r'$x$ [m]')
	ax.set_ylabel(r'$y$ [m]')

	fig.tight_layout()
	savefig(fig, 'linear')


def quadratic_plot() -> None:
	# Import data
	df = pd.read_csv(f'{DATA_FOLDER}/quadratic.csv')

	# Manipulate data
	x = df['x'].to_numpy()
	y = df['y'].to_numpy()

	# Create figure
	fig, ax = plt.subplots(1, 1)
	ax.plot(x, y)
	ax.scatter(-3, 4, c='k')
	ax.scatter(2, 1.5, c='k')
	ax.scatter(1.1, 2.3, c='k')

	ax.set_title('Quadratic interpolation test')
	ax.set_xlabel(r'$x$ [m]')
	ax.set_ylabel(r'$y$ [m]')

	fig.tight_layout()
	savefig(fig, 'linear')


def testing_plot() -> None:
	# Import data
	df = pd.read_csv(f'{DATA_FOLDER}/data.csv')

	# Create figure
	fig, ax = plt.subplots(1, 1)
	t = df['t'].to_numpy()  # * tau
	x = df['x'].to_numpy()  # * L
	y = df['y'].to_numpy()  # * L
	u = df['u'].to_numpy()  # * L / tau
	v = df['v'].to_numpy()  # * L / tau
	exact = analytical(x, 0)
	# ax.plot(x, abs(y - exact), label=r'$\pi / 4$ rad')
	ax.plot(x, y, label='numerical')
	ax.plot(x, exact, label='analytical')
	# draw_target(ax, 10.0, 0.0)

	# ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

	ax.set_title('Testing, theta = 45°')
	# ax.set_xlabel(r'$x$ [m]')
	# ax.set_ylabel(r'$\Delta y$ [m]')

	ax.legend()

	fig.tight_layout()
	savefig(fig, 'testing')


def convergence_plot() -> None:
	# Import data
	df = pd.read_csv(f'{DATA_FOLDER}/convergence.csv')

	# Manipulate data
	# df = df.groupby(by=df['theta'])
	# groups = list(df.groups.keys())

	# Create figure
	fig, ax = plt.subplots(1, 1)
	# for g in enumerate(groups):
	# dfcurr = df.get_group(g)
	ax.scatter(df['dt'], df['difference'])  # , label=f'{g:.2f} rad')

	ax.set_xscale('log')
	ax.set_yscale('log')

	# ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

	ax.set_title('Convergence')
	ax.set_xlabel(r'$dt$ [s]')
	ax.set_ylabel(r'$\Delta y$ [m]')

	# ax.legend(title='Launch angle', ncol=2)

	fig.tight_layout()
	savefig(fig, 'convergence')


def interpolation_plot() -> None:
	# Import data
	df = pd.read_csv(f'{DATA_FOLDER}/interpolationOrder.csv')

	# Manipulate data
	df = df.groupby(by=df['sol_number'])
	groups = df.groups.keys()

	# Create figure
	fig, ax = plt.subplots(1, 1)
	for g in groups:
		dfcurr = df.get_group(g)
		theta = dfcurr['theta'].to_numpy()
		theta_shift = theta - np.mean(theta)
		ax.plot(dfcurr['order'], theta_shift, marker='o', label=g)

	ax.set_xscale('log', base=2)
	# ax.set_yscale('log')

	ax.set_title('Interpolation order')
	ax.set_xlabel(r'order')
	ax.set_ylabel(r'$\theta - < \theta >$ [rad]')

	ax.legend(title='Solution number')

	fig.tight_layout()
	savefig(fig, 'interpolation')


def main() -> None:
	print_constants()
	# linear_plot()
	# quadratic_plot()
	shooting_plot()
	residual_plot()
	final_plot()
	comparison_plot()
	# interpolation_plot()
	# testing_plot()
	# convergence_plot()
	plt.show()


if __name__ == '__main__':
	main()
