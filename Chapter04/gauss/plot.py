import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def gauss(x, mu, sigma):
	norm = 1 / (sigma * np.sqrt(2*np.pi))
	expon = ((x - mu) / sigma)**2 * 0.5
	return norm * np.exp(-expon)

plt.style.use(['grid', 'science', 'notebook'])

df = pd.read_csv('data.csv')

xdense = np.linspace(-5, 5, 250)
ydense = gauss(xdense, 0, 0.5)

fig, ax = plt.subplots(1, 1)
ax.plot(xdense, ydense, color='tab:green', label='gaussian')
ax.scatter(df['x'][::16], df['y'][::16], marker='.', label='random distrib\n(downsampled by 16)', zorder=10)

ax.set_title('Acceptance/rejection method')
ax.set_xlabel('x')
ax.set_ylabel('y')

ax.legend(frameon=True, fancybox=True, framealpha=0.8)

fig.tight_layout()
fig.savefig('gauss.png', dpi=200)
# plt.show()
