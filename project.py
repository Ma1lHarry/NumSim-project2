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

    w[0] = coeff[0] * y[1]* y[0]/N
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



def SEIR_stoch():
    return numpy.array([[-1, 1, 0, 0],
                       [0, -1, 1, 0],
                       [0, 0, -1, 1]])

def SEIR_prop(y, coeff):
    N = y[0] + y[1] + y[2] + y[3]
    w = numpy.zeros(3)

    w[0] = coeff[1] * y[2]/N * y[0]
    w[1] = coeff[0] * y[1]
    w[2] = coeff[2] * y[2]
    
    return w


def SEIR():
    x0 = (995, 0, 5, 0)
    tspan = (0, 120)

    alpha = 0.3
    beta = 0.3
    gamma = 1/7
    coeff = [alpha,beta, gamma]

    tvec, xarr = SSA(SEIR_prop, SEIR_stoch, x0, tspan, coeff)
    plt.plot(tvec,xarr[:,0],'b-',label="Suseptable")
    plt.plot(tvec,xarr[:,1],'y-', label="exposed")
    plt.plot(tvec,xarr[:,2],'r-', label="infected")
    plt.plot(tvec,xarr[:,3],'g-', label="Recovered")
    
    plt.legend()
    plt.show()

def SEIRD_stoch():
    return numpy.array([#S  E  I  R  D
                       [-1, 1, 0, 0, 0],
                       [0, -1, 1, 0, 0],
                       [0, 0, -1, 1, 0],
                       [0, 0, -1, 0, 1]
                       ])

def SEIRD_prop(y, coeff):
    N = y[0] + y[1] + y[2] + y[3] + y[4]
    w = numpy.zeros(4)

    w[0] = coeff[1] * y[2]/N * y[0]
    w[1] = coeff[0] * y[1]
    w[2] = coeff[2] * y[2]
    w[3] = coeff[3] * y[2]
    return w


def SEIRD():
    x0 = (995, 0, 5, 0, 0)
    tspan = (0, 120)

    alpha = 0.3
    beta = 0.3
    gamma = 1/7
    my = 1/100
    coeff = [alpha, beta, gamma, my]

    tvec, xarr = SSA(SEIRD_prop, SEIRD_stoch, x0, tspan, coeff)
    plt.plot(tvec,xarr[:,0],'b-',label="Suseptable")
    plt.plot(tvec,xarr[:,1],'y-', label="exposed")
    plt.plot(tvec,xarr[:,2],'r-', label="infected")
    plt.plot(tvec,xarr[:,3],'g-', label="Recovered")
    plt.plot(tvec,xarr[:,4],'k-', label="Dead")
    
    plt.legend()
    plt.show()

def SEIRDV_stoch():
    return numpy.array([#S  E  I  R  D  V  
                       [-1, 1, 0, 0, 0, 0],
                       [ 0,-1, 1, 0, 0, 0],
                       [ 0, 0,-1, 1, 0, 0],
                       [ 0, 0,-1, 0, 1, 0],
                       [-1, 0, 0, 0, 0, 1]
                       ])

def SEIRDV_prop(y, coeff):
    N = y[0] + y[1] + y[2] + y[3] + y[4] + y[5]
    w = numpy.zeros(5)

    w[0] = coeff[1] * y[2]/N * y[0]
    w[1] = coeff[0] * y[1]
    w[2] = coeff[2] * y[2]
    w[3] = coeff[3] * y[2]
    if(y[0] > 0):
        w[4] = coeff[4]
    else: 
        w[4] = 0

    return w


def SEIRDV():
    x0 = (995, 0, 5, 0, 0, 0)
    tspan = (0, 120)

    alpha = 0.3
    beta = 0.3
    gamma = 1/7
    my = 1/100
    omega = 5
    coeff = [alpha, beta, gamma, my, omega]

    tvec, xarr = SSA(SEIRDV_prop, SEIRDV_stoch, x0, tspan, coeff)
    plt.plot(tvec,xarr[:,0],'b-',label="Suseptable")
    plt.plot(tvec,xarr[:,1],'y-', label="exposed")
    plt.plot(tvec,xarr[:,2],'r-', label="infected")
    plt.plot(tvec,xarr[:,3],'g-', label="Recovered")
    plt.plot(tvec,xarr[:,4],'k-', label="Dead")
    plt.plot(tvec,xarr[:,5],'m-', label="Vacinated")
    
    plt.legend()
    plt.show()

SEIRDV()