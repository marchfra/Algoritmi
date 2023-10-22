import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

plt.style.use(['grid', 'science', 'notebook'])

df = pd.read_csv('data.csv')

fig, ax = plt.subplots(1, 1)
asymp = ax.axhline(np.pi / 2, linestyle='--', linewidth=1, color='red', alpha=0.5, label='asymptote')
six, = ax.plot(df['x'], df['F'], label='Si(x)')

ax.set_title('Si(x)')
ax.set_xlabel('x')
ax.set_ylabel('y')

ax.legend(handles=[six, asymp], frameon=True, fancybox=True, framealpha=0.8)

fig.tight_layout()
for fmt in ['png', 'pdf']:
	fig.savefig(f'Si(x).{fmt}', dpi=200)
# plt.show()
