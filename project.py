from gillespie import SSA
import numpy
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def SIR_stoch():
    return numpy.array([[-1, 1, 0],
                       [0, 1, 1]])


def SIR_prop(y, coeff):
    N = y[0] + y[1] + y[2]
    w = numpy.zeros(2)

    w[0] = coeff[0] * y[1]* y[0]/N
    w[1] =  coeff[1] * y[1]

    return y

def SIR():
    x0 = (995, 5, 0)
    tspan = (0, 120)

    beta = 0.3
    gamma = 1/7
    coeff = [beta, gamma]

    tvec, xarr = SSA(SIR_prop, SIR_stoch, x0, tspan, coeff)
    plt.plot(tvec,Xarr[:,0],'b-',label="Suseptable")
    plt.plot(tvec,Xarr[:,1],'r-', label="Infected")
    plt.plot(tvec,Xarr[:,2],'g-', label="Recovered")
    
    plt.legend()
    plt.show()

SIR()

