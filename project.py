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


def SEIRV1V2M_stoch():
    return numpy.array([#S  E  I  R V1 V2  M  D  K
                       [-1, 1, 0, 0, 0, 0, 0, 0, 0],   # S  -> E   (smitta)
                       [ 0, 1, 0, 0,-1, 0, 0, 0, 0],   # V1 -> E   (genombrott dos 1)
                       [ 0, 1, 0, 0, 0,-1, 0, 0, 0],   # V2 -> E   (genombrott dos 2)
                       [ 0,-1, 1, 0, 0, 0, 0, 0, 1],   # E  -> I   (+1 kunskap)
                       [ 0, 0,-1, 1, 0, 0, 0, 0, 0],   # I  -> R
                       [-1, 0, 0, 0, 1, 0, 0, 0, 0],   # S  -> V1
                       [ 0, 0, 0,-1, 1, 0, 0, 0, 0],   # R  -> V1
                       [ 0, 0, 0, 0,-1, 1, 0, 0, 0],   # V1 -> V2
                       [ 0, 0, 0,-1, 0, 1, 0, 0, 0],   # R  -> V2
                       [ 0, 0, 0, 0,-1, 0, 1, 0, 0],   # V1 -> M
                       [ 0, 0, 0, 0, 0,-1, 1, 0, 0],   # V2 -> M
                       [ 0, 0,-1, 0, 0, 0, 0, 1, 0]    # I  -> D   (död)
                       ])

def SEIRV1V2M_prop(y, coeff):
    #     S      E      I      R     V1     V2      M      D 
    N = y[0] + y[1] + y[2] + y[3] + y[4] + y[5] + y[6] + y[7]
    w = numpy.zeros(12)
    
    w[0] = coeff[0] * y[2]/N * y[0] # S  -> E 
    w[1] = coeff[1] * y[2]/N * y[4] # V1 -> E 
    w[2] = coeff[2] * y[2]/N * y[5] # V2 -> E
    w[3] = coeff[4] * y[1]          # E  -> I
    w[4] = coeff[5] * y[2]          # I  -> R
    if(y[8] > 200):
        w[5] = coeff[6] * y[0]          # S  -> V1
        w[6] = coeff[7] * y[3]          # R  -> V1
        w[7] = coeff[8] * y[4]          # V1 -> V2
        w[8] = coeff[9] * y[3]          # R  -> V2
        
    else:
        w[5] = 0 # cant be vacinated if there are no sick
        w[6] = 0 # cant be vacinated if there are no sick
        w[7] = 0 # cant be vacinated if there are no sick
        w[8] = 0 # cant be vacinated if there are no sick
    w[9] = coeff[10] * y[4]          # V1 -> M
    w[10] = coeff[11] * y[5]        # V2 -> M
    w[11] = coeff[3] * y[2]         # I -> D

    return w


def SEIRV1V2M():
    x0 = (995, 0, 5, 0, 0, 0, 0, 0, 0)
    tspan = (0, 120)

    α  = 0.5    # infection rate S → E          (R0 = α/τ = 0.5·7 = 3.5)
    α1 = 0.2    # breakthrough V1 → E           (dose 1 gives ~60% protection: 0.4·α)
    α2 = 0.05   # breakthrough V2 → E           (dose 2 gives ~90% protection: 0.1·α)
    α3 = 0.01   # death I → D                   (IFR = α3/(τ+α3))
    β  = 1/5    # E → I, incubation ~β⁻¹ days
    τ  = 1/14    # I → R, infectious ~τ⁻¹ days
    μ1 = 0.01   # S → V1, ~1% of susceptibles vaccinated per day
    μ2 = 0.005  # R → V1, recovered get vaccinated less often
    μ3 = 1/21   # V1 → V2, dose 2 about 3 weeks after dose 1
    μ4 = 0.002  # R → V2, few recovered go straight to dose 2
    ω1 = 1/30   # V1 → M, immunity after dose 1 takes ~30 days
    ω2 = 1/14   # V2 → M, immunity after dose 2 takes ~14 days
    coeff = [α, α1, α2, α3 , β, τ, μ1, μ2, μ3, μ4, ω1, ω2]

    tvec, xarr = SSA(SEIRV1V2M_prop, SEIRV1V2M_stoch, x0, tspan, coeff)

    K_star = 200

    labels = ["Susceptible", "Exposed", "Infected", "Recovered",
              "Vaccinated dose 1", "Vaccinated dose 2", "Immune (vaccine)", "Dead"]
    colors = ['b', 'y', 'r', 'g', 'c', 'm', 'tab:orange', 'k']

    _, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(9, 8),
                                   gridspec_kw={'height_ratios': [3, 1]})

    # Population compartments (S, E, I, R, V1, V2, M, D)
    for i in range(8):
        ax1.step(tvec, xarr[:, i], color=colors[i], label=labels[i], where='post')

    # Knowledge K with the approval threshold
    ax2.step(tvec, xarr[:, 8], color='tab:purple', label="Knowledge K", where='post')
    ax2.axhline(K_star, color='gray', linestyle='--', label="K* (vaccine approved)")

    # Mark when the vaccine is approved, if it happens within tspan
    approved = numpy.where(xarr[:, 8] >= K_star)[0]
    if len(approved) > 0:
        t_vacc = tvec[approved[0]]
        for ax in (ax1, ax2):
            ax.axvline(t_vacc, color='gray', linestyle=':')
        ax1.text(t_vacc, ax1.get_ylim()[1] * 0.95, f" vaccine day {t_vacc:.0f}",
                 color='gray', va='top')

    ax1.set_ylabel("Number of individuals")
    ax1.legend(loc='center right')
    ax2.set_xlabel("Time (days)")
    ax2.set_ylabel("Cumulative cases")
    ax2.legend(loc='lower right')
    plt.tight_layout()
    plt.show()

SEIRV1V2M()
