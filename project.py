import gillespie
import numpy
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def prop(X, coeff):
    w = numpy.zeros(4)
    w[0] = coeff[0] * X[0]
    w[1] = coeff[2] * X[0] * X[1]
    w[2] = coeff[3] * X[1]
    w[3] = coeff[1] * X[0]**2
    return w

def stoch():
    return numpy.array([[1, 0], [-1, 1], [0, -1], [-1, 0]]) 

def uppg3():
    X0 = (100, 500)
    tspan = (0, 10)
    coeff = (10, 0.005, 0.01, 3)

    #With standard values
    tvec, Xarr = SSA(prop, stoch, X0, tspan, coeff)
    plt.plot(tvec,Xarr[:,0],'r-',label="Prey- OG", linewidth=5)
    plt.plot(tvec,Xarr[:,1],'r-', label="Predator - OG", linewidth=2)
    
    plt.legend()
    plt.show()


def SIR():
    