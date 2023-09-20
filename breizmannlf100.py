from scipy.integrate import quad, dblquad, nquad
import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import random

restm = 5e-4 #unit is GeV
lf = 1e+2
mu = restm * lf
s0 = mu/10

funcint = lambda p,t: (1/(2*np.pi*s0**2))*np.exp(-((p*np.cos(t)-mu)**2 + (p*np.sin(t))**2/s0**2))

funcnorm = dblquad(funcint, 0, 1, lambda x: 0, lambda x: 2*np.pi)
print(funcnorm)
#def funcinttry(p):
#    return p*(1/(2*np.pi*s0))*np.exp(-((p*np.cos(t)-mu)**2 + (p*np.sin(t))**2/s0**2))
#funcnormtry = quad(funcinttry, 0, 10)
#print(funcnormtry)

#define functions to compute the growth rate

func = lambda p, t: (1/(2*np.pi*s0**2))*np.exp(-((p*np.cos(t)-mu)**2 + (p*np.sin(t))**2/s0**2)) #define 2D distribution function
funcint = lambda p, t: p*(1/(2*np.pi*s0**2))*np.exp(-((p*np.cos(t)-mu)**2 + (p*np.sin(t))**2/s0**2)) #multiply by p before integrating over p
def g(t):
    return quad(funcint, 0, 1, args=(t,)) #integrate over p to obtain g(\theta)
#print("g at 2", g(2)[0])
#dgdtpre = lambda p, t: (numer1(p, t)*numer2(p, t))/(denom1(p, t)*denom2(p, t)) #define dg(p, \theta)/d\theta
#dgdtint = lambda p, t: p*(numer1(p, t)*numer2(p, t))/(denom1(p, t)*denom2(p, t)) #multiply by p before integrating
#define dg/d\theta p integrand
def dgdtint(p,t):
    return (6e-3)*(2*p**2*np.cos(t)*np.sin(t) - 2*p*(p*np.cos(t)-mu)*np.sin(t))*(1/(2*np.pi*s0**4))*np.exp(-((p*np.cos(t)-mu)**2 + (p*np.sin(t))**2/s0**2))
#define dg/d\theta by integrating over p
def dgdt(t):
    return quad(dgdtint, 0, 1, args=(t,)) #integrate over p to obtain dg(\theta)/d\theta

#constants
epsilon0 = 8.8e-12
e = 1.6e-19
n = 1e+3 #in m^-3, equivalent to 1e-3 per cm^3
nb = 1e-9
npl = 1e-5
m = 9.1e-31
enat = 0.30286
kgtoev = 5.62e+26
contrast = 1e-3 #nb/npl
cm3toev3 = (1.98e-14)**3
wp = np.sqrt(npl*cm3toev3*enat**2/m*kgtoev)
print(wp)
#def d2gdt2(t):
    #return quad()
def flist(t):
    return 2*np.array(g(t))*np.sin(t) + (np.cos(t) - kpar)*np.array(dgdt(t)) #Breizmann numerator defined as a function, returned as a list due to modification to solve problem with array
def f(t):
    return tuple(flist(t)) #list is now converted to tuple
def dfdtlist(t):
    return 2*np.array(g(t))*np.cos(t) + (kpar - 3*np.sin(t))*np.array(dgdt(t)) + (np.cos(t) - kpar)*np.array(d2gdt2(t))
def dfdt(t):
    return tuple(dfdtlist(t)) #list is now converted to tuple
#print(kperp, kpar, t0, g(t0)[0], dgdt(t0)[0], f(t0)[0])



#k values
kparbin = np.linspace(0.98, 1.02, 30) #range of kpar
kperbin = np.linspace(0.01, 1, 30) #range of kperp

#loop over k values
#compute the growth rate
kparvals = []
kperpvals = []
growth = []
with open('piclf100.txt', 'w') as file:
    for kpar in kparbin: #starts a loop in k1
        for kperp in kperbin: #starts a loop in k2
            k = np.sqrt(kperp**2 + kpar**2) #k in terms of its components
            muplus = (kpar/k + (kperp/k)*np.sqrt(k**2 - 1))/k #angular (Breizmann) integration upper limit as function of k_i
            muminus = (kpar/k - (kperp/k)*np.sqrt(k**2 - 1))/k #angular integration (Breizmann) lower limit as function of k_i
            t1 = np.arccos((kpar/k + (kperp/k)*np.sqrt(k**2 - 1))/k) #angular (Breizmann) integration upper limit as function of k_i
            t2 = np.arccos((kpar/k - (kperp/k)*np.sqrt(k**2 - 1))/k) #angular integration (Breizmann) lower limit as function of k_i
            t0 = (t1 + t2)/2 #mid-point between upper and lower limit
            print(kpar, kperp, muplus, muminus, t1, t2) #prints theta1 and theta2
            t1p = t1+1e-10 #shift the lower limit to force integrand to avoid zero or ultra-small values of lower limit t1 where it blows up; shift adjusted to current
            t2p = t2-1e-10
            int = lambda t: math.pi*contrast*restm*f(t)[0]/np.sqrt((np.cos(t)-np.cos(t1))*(np.cos(t2)-np.cos(t)))
            gr10 = quad(int, t1p, t0) # first half of the integral
            gr02 = quad(int, t0, t2p) #second half of the integral
            #print('kpar=', kpar, 'kperp=', kperp, 't1=', t1, 'integrand at t0=', int(t0), 'integrand at t1=', int(t1), 'first integral=', gr10[0], 'second integral=', gr02[0])
            gr = gr10 + gr02
            grnorm = gr/k**3
            data = kpar, kperp, grnorm[0]
            #kparvals.append(kpar)
            #kperpvals.append(kperp)
            #growth.append(grnorm[0])
            print(data) #print to check
            #save the growth rate
            file.write(' '.join(str(x) for x in data)+'\n')
#np.savetxt('testnew2.txt', np.column_stack[(kparvals,kperpvals,growth)])

#plt.contourf(kpar, kperp, grnorm[0])
#plt.show()
