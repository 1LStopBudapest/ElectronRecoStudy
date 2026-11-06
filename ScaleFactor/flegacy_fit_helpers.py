
import math
import mpmath
import numpy as np
from numpy.polynomial import chebyshev as Chev
from scipy.special import erf
from scipy.integrate import quad


'''
def crystalball(x, m0, sigma, alpha, n, sigma_2, tailLeft):
    
    result = []
    for m in x:
        rval=0
        t = (m-m0)/sigma
        t0 = (m-m0)/sigma_2

        absAlpha = abs(alpha)
        # if( tailLeft >= 0 ):
        if( 1 ): # force tail on right side like this...
            if (t>0):
                rval= np.exp(-0.5*t0*t0)
            elif (t > -absAlpha):
                rval= np.exp(-0.5*t*t)
            else:
                a = np.exp(-0.5*absAlpha*absAlpha)
                b = np.exp(n*(t+absAlpha))
                rval = a*b
        else:
            # rather fit high tail for n < 0
            if (t0<0):
                rval= np.exp(-0.5*t*t)
            elif (t0 < absAlpha):
                rval= np.exp(-0.5*t0*t0)
            else:
                absN = abs( n )
                a =  ((absN/absAlpha)**absN)*np.exp(-0.5*absAlpha*absAlpha)
                b= absN/absAlpha - absAlpha
                rval= a/((b + t0)**absN)

        result.append(rval)
    return result
'''

def mycrystalball(x_arr,N,alfa,n,x0,sigma):
    result = []
    N = mpmath.mpf(N)
    #alfa = mpmath.mpf(alfa)
    #n = mpmath.mpf(n)
    #x0 = mpmath.mpf(x0)
    #sigma = mpmath.mpf(sigma)

    #print N, alfa, n, x0, sigma

    alfabs = abs(alfa)
    nabs = abs(n)

    for x in x_arr:
        t0 = mpmath.mpf( (x-x0)/sigma )
        
        #D = math.sqrt(math.pi/2) * (1 + math.erf(alfabs/math.sqrt(2)))
        #C = (n/alfabs) * (1 / (n-1)) * np.exp(-( (alfabs**2) / 2))
        #N = 1 / (sigma*(C+D))
        B = mpmath.mpf(nabs/alfabs - alfabs)
        A = mpmath.power(nabs/alfabs, nabs) * mpmath.exp( -(alfabs**2) / 2 )

        try:
            if t0 > -alfa:
                inp = N*mpmath.exp(-((x-x0)**2) / (2*sigma**2) )
                result.append( float(mpmath.re(inp))  )  ##
            else:
                inp = N*A*mpmath.power(B-t0 , -n)  # was -n!!
                result.append( float(mpmath.re(inp))  )
        except TypeError as e:
            print(e)
            print(inp)

    
    return result


def bound_parameter(p0, p ):
    return p0 * (1 + 0.2*np.tanh(p-p0))

def crystalballPlusChev(x_arr, N,alfa,n,x0,sigma, *chevargs):

    #N = bound_parameter(p0[0], N)
    #alfa = bound_parameter(p0[1], alfa)
    #n = bound_parameter(p0[2], n)
    #x0 = bound_parameter(p0[3], x0)
    #sigma = bound_parameter(p0[4], sigma)


    result = []
    N = mpmath.mpf(N)
    #alfa = mpmath.mpf(alfa)
    #n = mpmath.mpf(n)
    #x0 = mpmath.mpf(x0)
    #sigma = mpmath.mpf(sigma)

    chebpars = list(chevargs)
    #for i in range(len(chebpars)):
    #    chebpars[i] = bound_parameter(p0[5+i], chebpars[i])

    alfabs = abs(alfa)
    nabs = abs(n)

    for x in x_arr:
        t0 = mpmath.mpf( (x-x0)/sigma )
        
        #D = math.sqrt(math.pi/2) * (1 + math.erf(alfabs/math.sqrt(2)))
        #C = (n/alfabs) * (1 / (n-1)) * np.exp(-( (alfabs**2) / 2))
        #N = 1 / (sigma*(C+D))
        B = mpmath.mpf(nabs/alfabs - alfabs)
        A = mpmath.power(nabs/alfabs, nabs) * mpmath.exp( -(alfabs**2) / 2 )

        try:
            if t0 > -alfa:
                inp = N*mpmath.exp(-((x-x0)**2) / (2*sigma**2) )
                result.append( float(mpmath.re(inp)) + Chev.chebval(x,chebpars ))  ##
            else:
                inp = N*A*mpmath.power(B-t0 , -n)  # was -n!!
                result.append( float(mpmath.re(inp)) + Chev.chebval(x,chebpars ))
        except TypeError as e:
            print(e)
            print(inp)

    return result

def crystalballPlusChevODR(B, x):
    return np.array(crystalballPlusChev(x, *B))

def boundedCrytalballPlusChev(p0):
    def boundedCrystalballPlusChevODR(B, x):
        return np.array(crystalballPlusChev(x, p0, *B))
    return boundedCrystalballPlusChevODR


    
# alpha > 0, n > 1 !
def crystalball_SOF(x, scale, alpha, n, xbar, sigma):
    n_over_alpha = n/abs(alpha)
    exp = np.exp(-0.5*alpha ** 2)
    A = (n_over_alpha)**n*exp
    B =  n_over_alpha - abs(alpha)
    C = n_over_alpha/(n-1)*exp
    D = np.sqrt(0.5*np.pi)*(1 + erf(abs(alpha)/np.sqrt(2)))
    N = 1/(sigma*(C + D))

    mask = (x - xbar)/sigma > -alpha

    gaussian = scale*N*np.exp(-0.5*((x[mask]-xbar)/sigma)**2)
    powerlaw = scale*N*A*(B - (x[~mask]-xbar)/sigma)**-n

    result = np.zeros_like(x)
    result[mask] = gaussian
    result[~mask] = powerlaw

    result = np.concatenate((powerlaw, gaussian))

    return result.tolist()


# this just doesn't work
def BWG(xarr, N, M, Gamma, k, sigma):
    result = []
    for E in xarr:
        def integrand(Ev):
            return k / ((Ev**2 - M**2)**2 + (M*Gamma)**2) / (sigma*np.sqrt(2*math.pi)) * np.exp( -(Ev-E)**2 / (2*(sigma**2))  )
        result.append(N * quad(integrand,-np.infty, np.infty)[0])
    return result


# this is wikipedia and makes no sense...
def BWG_voigt(xarr, N, M, Gamma, k, sigma):
    E = np.array(xarr)
    u1 = (E-M) / (math.sqrt(2)*sigma)
    u2 = (E+M) / (math.sqrt(2)*sigma)
    a = k*math.pi/(2*(sigma**2))

    def integrand(t):
        return math.exp(-(t**2)) / ( (u1-t)**2 * (u2-t)**2 + a**2 )
    H2 = a/math.pi * quad(integrand,-np.infty, np.infty)[0]

    return N*H2/(sigma**2 * 2 *math.sqrt(math.pi))

# this is correct!
def Voigt(xarr, N, M, Gamma, sigma):
    result = []
    for E in xarr:
        u1 = (E-M)/(np.sqrt(2)*sigma)
        u2 = (E+M)/(np.sqrt(2)*sigma)
        a = Gamma*M/(2*(sigma**2))
        def integrand(t):
            return np.exp(-(t**2)) / ( (u1-t)**2 * ((u2-t)**2) + a**2 )
        result.append( N * a / ( sigma**2 *2*np.sqrt(math.pi)*math.pi  ) *   quad(integrand,-np.infty, np.infty, limit=150)[0]   )
    return result

def voigtPlusChev(xarr, N, M, Gamma, sigma, *chevargs):
    chebpars = list(chevargs)
    result = []
    for E in xarr:
        u1 = (E-M)/(np.sqrt(2)*sigma)
        u2 = (E+M)/(np.sqrt(2)*sigma)
        a = Gamma*M/(2*(sigma**2))
        def integrand(t):
            return np.exp(-(t**2)) / ( (u1-t)**2 * ((u2-t)**2) + a**2 )
        result.append( N * a / ( sigma**2 *2*np.sqrt(math.pi)*math.pi  ) *   quad(integrand,-np.infty, np.infty, limit=150)[0]  + Chev.chebval(E,chebpars )  )
    return result

def voigtPlusChevODR(B, x):
    return np.array(voigtPlusChev(x, *B))


def cmsshape(xarr, alpha, beta, gamma):
    result = []
    peak = 3.1
    for x in xarr:

        erf = math.erfc((alpha - x) * beta)
        u = (x - peak)*gamma

        if(u < -70):
            u = 1e20
        elif( u>70 ):
            u = 0
        else:
            u = np.exp(-u)
        result.append(erf*u)
    return result




def mycmsshape(x, alfa,beta,A,b,c):
    result = []
    for xi in x:
        expo = float(A*mpmath.exp(b*xi+c))
        errorf = math.erfc((alfa - xi) * beta)
        result.append(expo*errorf)
    return result


def gauss(x,A,x0,sigma):
    result = []
    for xi in x:
        result.append( A*np.exp(  -(xi-x0)**2 / (2* (sigma**2)) )  )
    return result

