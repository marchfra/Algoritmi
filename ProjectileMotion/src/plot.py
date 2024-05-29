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

	# Create figure
	fig, ax = plt.subplots(1, 1)
	theta = df['theta'].to_numpy()
	res = df['res'].to_numpy() * L
	ax.plot(theta, res, marker='o')
	ax.axhline(0, c='k', ls='--', alpha=0.7)

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


def comparison_plot() -> None:
	# Import data
	df = pd.read_csv(f'{DATA_FOLDER}/noFriction.csv')

	# Manipulate data
	df = df.groupby(by=df['theta'])
	groups = list(df.groups.keys())

	# Create figure
	fig, ax = plt.subplots(1, 1)
	for g in groups:
		curr_df = df.get_group(g)
		x = curr_df['x'].to_numpy() * L
		y = curr_df['delta_y'].to_numpy() * L
		ax.plot(x, y, label=f'{g:.2f} rad')
	draw_target(ax, 10.0, 0.0)

	ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

	ax.set_title('Comparison with analytical solution')
	ax.set_xlabel(r'$x$ [m]')
	ax.set_ylabel(r'$\Delta y$ [m]')

	ax.legend(title='Launch angle', ncol=1)

	fig.tight_layout()
	savefig(fig, 'comparison')


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


def dt_search() -> None:
	df = pd.read_csv(f'{DATA_FOLDER}/dt_search.csv')

	fig, ax = plt.subplots(1, 1)
	for i, theta in enumerate(df.columns[1:]):
		theta = df[theta].to_numpy()
		mean = np.mean(theta)
		ax.plot(df['dt'], theta - mean, marker='o', label=i)
		ax.scatter(df['dt'][4], (theta - mean)[4], marker='o', s=400)
		ax.scatter(df['dt'][4], (theta - mean)[4], marker='o', s=300, c='white')

	ax.set_xscale('log', base=2)

	ax.set_title('Time step size search')
	ax.set_xlabel(r'time step size')
	ax.set_ylabel(r'$\theta - < \theta >$ [rad]')

	ax.legend(title='Solution number')

	fig.tight_layout()
	savefig(fig, 'dt_search')


def order_search() -> None:
	df = pd.read_csv(f'{DATA_FOLDER}/order_search.csv')

	fig, ax = plt.subplots(1, 1)
	for i, theta in enumerate(df.columns[1:]):
		theta = df[theta].to_numpy()
		mean = np.mean(theta)
		ax.plot(df['order'], theta - mean, marker='o', label=i)
		# ax.scatter(df['order'][4], (theta - mean)[4], marker='o', s=400)
		# ax.scatter(df['order'][4], (theta - mean)[4], marker='o', s=300, c='white')

	# ax.set_xscale('log', base=2)

	ax.set_title('Polynomial fit order search')
	ax.set_xlabel(r'order')
	ax.set_ylabel(r'$\theta - < \theta >$ [rad]')

	ax.legend(title='Solution number')

	fig.tight_layout()
	savefig(fig, 'order_search')


def main() -> None:
	# print_constants()
	# shooting_plot()
	# dt_search()
	# order_search()
	# residual_plot()
	# final_plot()
	# comparison_plot()
	plt.show()


if __name__ == '__main__':
	main()
