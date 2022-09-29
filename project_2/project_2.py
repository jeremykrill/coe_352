import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math


def forward_euler(num_t, dt, K, f, a, m):
    for j in range(num_t - 1):
        K[:,j + 1] = K[:,j] + np.linalg.solve(m, dt*(f[:,j] - np.dot(a, K[:,j])))

def backward_euler(num_t, dt, K, f, a, m):
    for j in range(num_t - 1):
        K[:,j + 1] = np.linalg.solve(m*dt + a, f[:,j+1] + np.dot(m, K[:,j])*dt)

def func(x, t):
    return ((math.pi ** 2) - 1) * math.exp(-t) * math.sin(math.pi * x)

def g_quad(x, t, h):
    f1 = func((h * (-1) * math.sqrt(1/3) + 2 * x) / 2, t)
    f2 = func((h * math.sqrt(1/3) + 2 * x) / 2, t)
    return (f1 + f2) * h / 2

if __name__ == '__main__':
    print('This program solves the heat transfer problem u_t - u_xx = (\pi^2 - 1)e^{-t}sin(\pi x).')
    print('It uses initial and Dirichlet BCs u(x,0) = sin(\pi x), u(0,t) = u(1,t) = 0.')
    discretization = input('Please input "f" for forward Euler, or "b" for backward Euler: ')
    num_t_steps = int(input('Please input the number of time steps: '))
    dt = 1 / (num_t_steps - 1)
    nodes = int(input('Please input the number of nodes: '))
    dx = 1 / (nodes - 1)

    u = np.zeros((nodes, num_t_steps))
    f = np.zeros((nodes - 1, num_t_steps))
    K = np.zeros((nodes - 2, num_t_steps))
    m = np.zeros((nodes - 2, nodes - 2))
    a = np.zeros((nodes - 2, nodes - 2))

    for i in range(nodes - 2):
      m[i][i] = 4 * dx / 6
    for i in range(nodes - 3):
      m[i][i + 1] = dx / 6
      m[i + 1][i] = dx / 6

    for i in range(nodes - 2):
      a[i][i] = 2 / dx
    for i in range(nodes - 3):
      a[i][i + 1] = -1 / dx
      a[i + 1][i] = -1 / dx

    t_mesh = np.zeros(num_t_steps)
    for i in range(1, num_t_steps):
        t_mesh[i] = t_mesh[i - 1] + 1 / (num_t_steps - 1)
    x_mesh = np.zeros(nodes)
    for i in range(1, nodes):
        x_mesh[i] = x_mesh[i - 1] + 1 / (nodes - 1)

    for j in range(num_t_steps):
      for i in range(1, nodes - 1):
        f[i,j] = g_quad(x_mesh[i], t_mesh[j], dx)
    f = np.delete(f, (0), axis=0)

    if discretization == 'f':
        forward_euler(num_t_steps, dt, K, f, a, m)
    else:
        backward_euler(num_t_steps, dt, K, f, a, m)

    u[1:-1] = K

    fig, ax = plt.subplots()
    ax.plot(x_mesh, u[:,-1], label='none')
    ttl = 'Heat transfer solution for ' + str(nodes) + ' nodes and ' + str(num_t_steps) + ' time steps at t = 1.'
    ax.set(xlabel='x', ylabel='u(x,t)', title=ttl)
    ax.grid()
    plt.show()
