# Project #2: Realistic Projectile Motion

We consider the motion of a particle subject to drag force

$$
\begin{cases}
    \ddot{x} = F_{d, x} \\
    \ddot{y} = -g + F_{d, x}
\end{cases}
$$

where $\vec{F}_{d} = -Bv\vec{v}$ is the drag force due to air resistance (for the present calculation you can use $B = 4.0 \cdot 10^{-5} \text{ m}^{-1}$.

This force is always opposite to velocity and therefore remember to write it in vector components.

For a given inital velocity $v_0$ and distance $L$ to a target, determine the angles (if any) you must orient your cannon at in order to hit the target.

## Writing the project

The project shold be a .pdf file no more than about 10 pages long (excluding the code listed in the appendix).

The document structure should consist of:
- An abstract/introductory part where the physical system is explained and why we need to resort to numerical integration
- A test section where the numerical method is validated against known analytical/reference solution
- A model study of the problem including plots
- A final summary/discussion
- An appendix including the code used for the project, using fixed-sized fonts
