import sys
import os

import numpy as np
from scipy.integrate import solve_ivp
from errno import EEXIST

def function(x,y):
    # This is the right-hand side of the first-order ordinary differential 
    # equation dx/dt = fun.
    fun = 2 * x * y + y * 2
    return fun

def EulerCromer():
    
    x = -1
    y = -1

    dx = 0.01

    ex = np.array([])
    ey = np.array([])

    while (x < 1):
        y = y + dx * function(x, y)
        x = x + dx
        ex = np.append(ex, x)
        ey = np.append(ey, y)
    answer = zip(ex,ey)
    return answer


def RungeKutta(x1, x2, h, y0):
    '''Runge-Kutta solution to differential equations'''
    
    dp = int((x2 - x1) / h)
    x = np.linspace(x1, x2, dp + 1, endpoint = True)
    y = np.ones(dp + 1,) * y0

    for i in range(0, int(dp)):
        k1 = function(x[i], y[i])
        k2 = function(x[i] + h/2, y[i] + h*k1/2)
        k3 = function(x[i] + h/2, y[i] + h*k2/2)
        k4 = function(x[i] + h, y[i] + h*k3)
              
        y[i + 1] = y[i] + h/6 * (k1 + 2*k2 + 2*k3 + k4)


    answer = zip(x,y)
    return answer

def solve(x1, x2, h):
    dp = (x2 - x1) / h
    x = np.linspace(x1, x2, int(dp + 1), endpoint = True)
    y = function(x = -1, y = -1)

    ex = np.array([])
    ey = np.array([])

    solver = solve_ivp(function, t_span = [x1, x2], y0 = [y,], method='LSODA', t_eval= x)
    answer = zip(solver.t, solver.y[0])
    return answer

if __name__ == '__main__':
        
    cwd = os.getcwd()

    try:
        os.makedirs(cwd + '/data')
    except OSError as exc:
        if exc.errno == EEXIST and os.path.isdir(cwd + '/data'):
            pass
        else: 
            raise

    filepath = cwd + '/data/' + 'test explicit' + '.txt'
    with open(filepath, 'w') as file:
        for ix, iy in RungeKutta(x1 = 0, x2 = 5, h = 0.01, y0 = 1):
            file.write(str(ix) + ',' + str(iy) + '\n')