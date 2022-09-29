This program uses 1D Galerkin code to solve the heat transfer problem

$$u_t - u_{xx} = f(x,t), \quad (x,t) \in (0,1) \times (0,1)$$	

with initial and Dirichlet boundary conditions
$$
\begin{align*}
	&u(x,0) = \sin{(\pi x)},\\
	&u(0,t) = u(1,t) = 0
\end{align*}
$$
and function 

$$f(x,t) = (\pi^2 - 1)e^{-t} \sin{(\pi x)}$$.
