from gillespie import SSA
import numpy
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def SIR_stoch():
    return numpy.array([[-1, 1, 0],
                       [0, -1, 1]])


def SIR_prop(y, coeff):
    N = y[0] + y[1] + y[2]
    w = numpy.zeros(2)
    w[0] = coeff[0] * y[1] / N * y[0]
    w[1] =  coeff[1] * y[1]

    return w

def SIR():
    x0 = (995, 5, 0)
    tspan = (0, 120)

    beta = 0.3
    gamma = 1/7
    coeff = [beta, gamma]

    tvec, xarr = SSA(SIR_prop, SIR_stoch, x0, tspan, coeff)
    plt.plot(tvec,xarr[:,0],'b-',label="Suseptable")
    plt.plot(tvec,xarr[:,1],'r-', label="Infected")
    plt.plot(tvec,xarr[:,2],'g-', label="Recovered")
    
    plt.legend()
    plt.show()


def SEIR_prop(y, coeff):
    w = numpy.zeros(3)
    N = y[0] + y[1] + y[2] + y[3]

    w[0] = coeff[1] * y[2] / N * y[0]
    w[1] = coeff[0] * y[1]
    w[2] = coeff[2] * y[2]

    return w


def SEIR_stoch():
    return numpy.array([[-1, 1, 0, 0],
                        [0, -1, 1, 0],
                        [0, 0, -1, 1]])

def SEIR():
    x0 = (995, 0, 5, 0)
    tspan = (0, 120)

    beta = 0.3
    gamma = 1/7
    alfa = 1/2
    coeff = [alfa, beta, gamma]

    tvec, xarr = SSA(SEIR_prop, SEIR_stoch, x0, tspan, coeff)
    plt.plot(tvec,xarr[:,0],'b-',label="Suseptable")
    plt.plot(tvec,xarr[:,1],'y-', label="Exposed")
    plt.plot(tvec,xarr[:,2],'r-', label="Infected")
    plt.plot(tvec,xarr[:,3],'g-', label="Recovered")
    
    plt.legend()
    plt.show()


def SEIRD_prop(y, coeff):
    w = numpy.zeros(4)
    N = y[0] + y[1] + y[2] + y[3] + y[4]

    w[0] = coeff[1] * y[2] / N * y[0]
    w[1] = coeff[0] * y[1]
    w[2] = coeff[2] * y[2]
    w[3] = coeff[3] * y[2]

    return w


def SEIRD_stoch():
    return numpy.array([[-1, 1, 0, 0, 0],
                        [0, -1, 1, 0, 0],
                        [0, 0, -1, 1, 0],
                        [0, 0, -1, 0, 1]])

def SEIRD():
    x0 = (995, 0, 5, 0, 0)
    tspan = (0, 120)

    beta = 0.3
    gamma = 1/7
    alfa = 1/2
    my = 0.05
    coeff = [alfa, beta, gamma, my]

    tvec, xarr = SSA(SEIRD_prop, SEIRD_stoch, x0, tspan, coeff)
    plt.plot(tvec,xarr[:,0],'b-',label="Suseptable")
    plt.plot(tvec,xarr[:,1],'y-', label="Exposed")
    plt.plot(tvec,xarr[:,2],'r-', label="Infected")
    plt.plot(tvec,xarr[:,3],'g-', label="Recovered")
    plt.plot(tvec,xarr[:,4],'k-', label="Dead")
    
    plt.legend()
    plt.show()


def SEIRDV_prop(y, coeff):
    w = numpy.zeros(5)
    N = numpy.sum(y)

    w[0] = coeff[1] * y[2] / N * y[0]
    w[1] = coeff[0] * y[1]
    w[2] = coeff[2] * y[2]
    w[3] = coeff[3] * y[2]
    if(y[0] > 0):
        w[4] = coeff[4]
    else:
        w[4] = 0

    return w


def SEIRDV_stoch():
    return numpy.array([[-1, 1, 0, 0, 0, 0],
                        [0, -1, 1, 0, 0, 0],
                        [0, 0, -1, 1, 0, 0],
                        [0, 0, -1, 0, 1, 0],
                        [-1, 0, 0, 0, 0, 1]])

def SEIRDV():
    x0 = (995, 0, 5, 0, 0, 0)
    tspan = (0, 120)

    beta = 0.3
    gamma = 1/7
    alfa = 1/2
    my = 0.05
    v = 5
    coeff = [alfa, beta, gamma, my, v]

    tvec, xarr = SSA(SEIRDV_prop, SEIRDV_stoch, x0, tspan, coeff)
    plt.plot(tvec,xarr[:,0],'b-',label="Suseptable")
    plt.plot(tvec,xarr[:,1],'y-', label="Exposed")
    plt.plot(tvec,xarr[:,2],'r-', label="Infected")
    plt.plot(tvec,xarr[:,3],'g-', label="Recovered")
    plt.plot(tvec,xarr[:,4],'k-', label="Dead")
    plt.plot(tvec,xarr[:,5],'m-', label="Vaccinated")
    
    plt.legend()
    plt.show()

SEIRDV()

