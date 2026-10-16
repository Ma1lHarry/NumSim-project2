import math
import random
import numpy

#==================================
# Gillespies algoritm
# Numeriska metoder och Simulering
#==================================
def SSA(prop,stoch,X0,tspan,coeff):
    tvec=numpy.zeros(1)
    tvec[0]=tspan[0]
    Xarr=numpy.zeros([1,len(X0)])
    Xarr[0,:]=X0
    t=tvec[0]
    X=X0
    sMat=stoch()
    counter = 0

    while t<tspan[1]:
        counter += 1
        re=prop(X,coeff)
        a0=sum(re)
        if a0>1e-10:
            # Bestäm reaktionstid
            r1=random.random()
            tau=-numpy.log(r1)/a0
            t+=tau
            tvec=numpy.append(tvec,t)

            # Bestäm vilken reaktion
            r2=random.random()
            cre=numpy.cumsum(re)
            cre=cre/cre[-1]
            r=0

            while cre[r]<r2:
                r+=1
            X=X+sMat[r,:]
            Xarr=numpy.vstack([Xarr,X])
        else: # Alla reaktionshastigheter 0
            print('Simuleringen avslutad vid t=', t)
            return tvec, Xarr

    return tvec, Xarr