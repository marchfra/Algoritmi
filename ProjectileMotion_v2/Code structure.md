# Code structure

- Constants output

- Shooting plot
	- Define $\theta$ range
	- Define initial conditions
	- Integrate while printing

- Find roots of residual
	- Define initial conditions
	- Integrate with exit condition $y(t-1) > h \wedge y(t) < h$
	- Return interpolation $- h$
