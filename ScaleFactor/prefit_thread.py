import os

import ROOT
from ROOT import TFile
import numpy as np
from numpy.polynomial import chebyshev as Chev
import math
import mpmath
import multiprocessing
from scipy.optimize import curve_fit
from scipy.integrate import quad
from scipy.special import erf
import warnings
from scipy.integrate import IntegrationWarning
import matplotlib.pyplot as plt
import yaml
import ctypes
import scipy.odr as odr



def get_parser():
    import argparse
    argParser = argparse.ArgumentParser(description = "Argument parser")
    argParser.add_argument('config_file',help='config_file')
    argParser.add_argument('--threads', nargs='?', const=1, type=int)
    return argParser




class Infix:
    def __init__(self, function):
        self.function = function
    def __ror__(self, other):
        return Infix(lambda x, self=self, other=other: self.function(other, x))
    def __or__(self, other):
        return self.function(other)
    def __rlshift__(self, other):
        return Infix(lambda x, self=self, other=other: self.function(other, x))
    def __rshift__(self, other):
        return self.function(other)
    def __call__(self, value1, value2):
        return self.function(value1, value2)

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
            print e
            print inp

    
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


def trapcalc(x,y):
    sum = 0
    for i in range(len(x)-1):
        if(y[i] < y[i+1]):
            sum += y[i] * (x[i+1]-x[i])  +  (y[i+1]-y[i]) * (x[i+1]-x[i])/2
        else:
            sum += y[i+1] * (x[i+1]-x[i])  +  (y[i]-y[i+1]) * (x[i+1]-x[i])/2
    return sum

def useBigger(a,b):
    return a if a >= b else b

def useSmaller(a,b):
    return a if a <= b else b

def qsum(*args):
    sm = 0
    for arg in args:
        sm += arg**2
    return np.sqrt(sm)
qs = Infix(qsum)


def getIntegralRangeSystError(integral, integrand, xrange_min, xrange_max, binwidth):
    #print "hi"
    res = useBigger( abs(integral - quad(integrand, xrange_min*0.95, xrange_max*1.05)[0]), abs(integral - quad(integrand, xrange_min*1.05, xrange_max*0.95)[0]) ) / binwidth
    #print "hi2"
    return res


def subtract_background(x, y, y_pm, chevpars, bGetError = False, chev_std = []):
    y_sub = []
    y_pm_sub = []

    if bGetError:
        parameters = []
        for i in range(100):
            rel_std = list(chev_std)
            for j in range(len(chev_std)):
                rel_std[j] = abs(chev_std[j] / chevpars[j])

            sign = np.random.choice([-1,1])
            devi = abs(np.random.normal(0, max(rel_std)))

            parameters.append(   [chevpars[k] * ( 1 + sign * devi) for k in range(len(chevpars))]   )


    for i in range(len(y)):
        bkg = Chev.chebval(x[i],chevpars)
        y_sub.append(y[i] - bkg if y[i] - bkg > 0 else 0)

        if bGetError:
            errorbar_scale = abs(y_pm[i]/y[i])  * abs(y[i] - bkg) 

            ybksyst = []
            for k in range(len(parameters)):
                ybksyst.append( Chev.chebval(x[i],parameters[k] ))
            #y_pm_sub.append( errorbar_scale |qs| np.std(ybksyst)  )    # FIXME FIXME: Make this a parameter?
            y_pm_sub.append( errorbar_scale   )

    if bGetError:
        return y_sub, y_pm_sub
    return y_sub

def get_bksub_integral(x, y, y_pm, chevpars, chev_std, int_range, repeats):

    #print chevpars
    #print chev_std
    #stri = ""
    #for i in range(len(chev_std)):
    #    stri += str(abs(chev_std[i])/chevpars[i]*100)+"% "
    #print "params: "+str(chevpars)
    #print "Cheverrors: "+stri


    x_inrange = []
    y_inrange = []
    y_pm_inrange = []


    x_inrange_systlow = []
    y_inrange_systlow = []
    y_pm_inrange_systlow = []


    x_inrange_systhigh = []
    y_inrange_systhigh = []
    y_pm_inrange_systhigh = []

    for i in range(len(x)):

        if int_range[0]*0.95 <= x[i] <= int_range[1]*1.05:
            x_inrange_systhigh.append(x[i])
            y_inrange_systhigh.append(y[i])
            y_pm_inrange_systhigh.append(y_pm[i])


        if int_range[0] <= x[i] <= int_range[1]:
            x_inrange.append(x[i])
            y_inrange.append(y[i])
            y_pm_inrange.append(y_pm[i])


        if int_range[0]*1.05 <= x[i] <= int_range[1]*0.95:
            x_inrange_systlow.append(x[i])
            y_inrange_systlow.append(y[i])
            y_pm_inrange_systlow.append(y_pm[i])
    
        if x[i] > int_range[1]*1.05:
            break


    y_sub = subtract_background(x_inrange, y_inrange, y_pm_inrange, chevpars, bGetError=False)

    y_sub_rangesyt_high =  subtract_background(x_inrange_systhigh, y_inrange_systhigh, y_pm_inrange_systhigh, chevpars, bGetError=False)
    y_sub_rangesyt_low =  subtract_background(x_inrange_systlow, y_inrange_systlow, y_pm_inrange_systlow, chevpars, bGetError=False)

    integral_w = sum(y_sub)

    integral_rangesyst_std = useBigger(abs(integral_w - sum(y_sub_rangesyt_high)), abs(integral_w - sum(y_sub_rangesyt_low)))


    integrals_for_std = []

    for i in range(repeats):

        rel_std = list(chev_std)
        for j in range(len(chev_std)):
            rel_std[j] = abs(chev_std[j] / chevpars[j])


        sign = np.random.choice([-1,1])
        devi = abs(np.random.normal(0, max(rel_std)))

        parameters = [chevpars[k] * ( 1 + sign * devi) for k in range(len(chevpars))]

        #print parameters
        ysub_i = subtract_background(x_inrange, y_inrange, y_pm_inrange, parameters, bGetError=False)
        integral_w_i = sum(ysub_i)
        #print integral_w_i
        if 0.5 * integral_w < integral_w_i < 2 * integral_w:
            integrals_for_std.append( sum(ysub_i) )

    # calculate integral with expecting a linear background
    x1 = int_range[0]
    x2 = int_range[1]
    y1 = Chev.chebval(x1,chevpars )
    y2 = Chev.chebval(x2,chevpars )
    a = (y2-y1)/(x2-x1)
    b = y1 - a*x1

    integral_linbk_w = 0
    for i in range(len(x_inrange)):
        trapstep = y_inrange[i] - (a*x_inrange[i]+b)
        integral_linbk_w += trapstep if trapstep > 0 else 0


    #print integrals_for_std
    #print integral_w, str(abs(np.sqrt(np.std(integrals_for_std)**2 + integral_w))/integral_w*100)+"%"
    return integral_w, np.sqrt(np.std(integrals_for_std)**2 + integral_w), integral_linbk_w, integral_rangesyst_std



def odr_dualfit_cb(x, y, x_pm, y_pm, chevpars, cbpars):

    p0 = []
    for cp in cbpars:
        p0.append(cp)
    for cp in chevpars:
        p0.append(cp)

    '''
    popt_sig, pcov_sig = curve_fit(crystalballPlusChev, x, y, p0, sigma=y_pm, absolute_sigma=True, maxfev=5000)
    perr_sig = np.sqrt(np.diag(pcov_sig))

    print "success"

    class myoutput:
        beta = 0
        sd_beta = 0

    ret = myoutput()
    ret.beta = popt_sig
    ret.sd_beta = perr_sig
    return ret
    '''

    linear = odr.Model(crystalballPlusChevODR)
    mydata = odr.RealData(x, y, sx=x_pm, sy=y_pm)  # 


    myodr = odr.ODR(mydata, linear, beta0=p0)
    myoutput = myodr.run()
    return myoutput

def odr_dualfit_voigt(x, y, x_pm, y_pm, chevpars, vpars):
    p0 = []
    for vp in vpars:
        p0.append(vp)
    for cp in chevpars:
        p0.append(cp)

    linear = odr.Model(voigtPlusChevODR)
    mydata = odr.RealData(x, y, sx=x_pm, sy=y_pm)

    myodr = odr.ODR(mydata, linear, beta0=p0)
    myoutput = myodr.run()
    return myoutput



def peak_fit(function, x, y, y_pm, config_file, isData, index, fix_sigma = False, fix_sigma_Delta = False):

    if isData and index in config_file['data_override_plots_xrange']:
        xrange_min = config_file['data_override_plots_xrange'][index][0]
        xrange_max = config_file['data_override_plots_xrange'][index][1]
    else:
        xrange_min = config_file['plots_xrange'][0]
        xrange_max = config_file['plots_xrange'][1]

    binwidth = abs(x[1] - x[0])
    if function == "gauss":
        bounds = ((-np.inf,-np.inf,-np.inf),(np.inf,np.inf,np.inf))  # fixme 

        popt_sig, pcov_sig = curve_fit(gauss, x, y, config_file['gauss_p0'], sigma=y_pm, absolute_sigma=True, bounds=bounds,maxfev=5000)
        perr_sig = np.sqrt(np.diag(pcov_sig))
        
        def integrand(x):
            return gauss([x], popt_sig[0],popt_sig[1],popt_sig[2])[0]
        integral = quad(integrand, xrange_min, xrange_max)

        integral_w = integral[0] / binwidth
        Delta_integral = math.sqrt(2*math.pi) * math.sqrt((popt_sig[2]**2) * (perr_sig[0]**2) + (popt_sig[0]**2) * (perr_sig[2]**2)) / binwidth

        intrange_syst_std = getIntegralRangeSystError(integral[0], integrand, xrange_min, xrange_max, binwidth)


    elif function == "crystalball":

        # force left end to come off
        if isData and index in config_file['help_cb_zeros']:
            y = list(y)
            for i in range(len(x)):
                if x[i] >= config_file[('MC_' if not isData else '') + 'help_cb_zeros'][index]:
                    break
                y[i] = 0
                

        if not isData and index in config_file['MC_override_crystalball_p0']:
            p0 = config_file['MC_override_crystalball_p0'][index]
        elif isData and index in config_file['data_override_crystalball_p0']:
            p0 = config_file['data_override_crystalball_p0'][index]
        else:
            p0 = config_file['crystalball_p0']


        if config_file["resonance"] == "Z":
            bounds = ((0,0,1,70,0),(np.inf,np.inf,np.inf,110,np.inf)) 
        else:
            bounds = ((0,0,1,2.7,0),(np.inf,np.inf,np.inf,3.3,np.inf)) 



        if not fix_sigma: # or index == 10:
            popt_sig, pcov_sig = curve_fit(mycrystalball, x, y, p0, bounds=bounds, sigma=y_pm, absolute_sigma=True,maxfev=5000)
            perr_sig = np.sqrt(np.diag(pcov_sig))
        else:
            bounds = ((0,0,1,2.7),(np.inf,np.inf,np.inf,3.3)) 
            popt_sig, pcov_sig = curve_fit(lambda x, N, alfa, n, x0: mycrystalball(x, N, alfa, n, x0, fix_sigma), x, y, p0[0:-1], bounds=bounds, sigma=y_pm, absolute_sigma=True,maxfev=5000)
            perr_sig = np.sqrt(np.diag(pcov_sig))

            popt_sig = np.append(popt_sig, [fix_sigma])
            perr_sig = np.append(perr_sig, [fix_sigma_Delta])

        if index in config_file[('MC_' if not isData else '') + 'dualfit_bounds_cb']:  
            integral_w, Delta_integral, intrange_syst_std = 1,1,1   #it's getting dual fit, skip now
        else:
            integral_w, Delta_integral, intrange_syst_std = getIntegralErrorWithMCMethod(mycrystalball, popt_sig, perr_sig, binwidth, config_file, index, isData)

    elif function == "bwg":
        bounds = ((0,85,0.01,0.00001),(np.inf,95,50,20))  

        p0 = config_file['BWG_p0']
        p0[0] = 10000000 * (max(y) / 5000)

        if isData and index in config_file['data_override_bwg_p0']:
            p0 = config_file['data_override_bwg_p0'][index]

        popt_sig, pcov_sig = curve_fit(Voigt, x, y, p0, sigma=y_pm, absolute_sigma=True, bounds=bounds,maxfev=5000)
        perr_sig = np.sqrt(np.diag(pcov_sig))
        

        if config_file['bUseBWGauss'] and index in config_file[('MC_' if not isData else '') + 'dualfit_bounds_bwg']:
            integral_w, Delta_integral, intrange_syst_std =  1,1,1  #it's getting dual fit, skip now
        else:
            integral_w, Delta_integral, intrange_syst_std = getIntegralErrorWithMCMethod(Voigt, popt_sig, perr_sig, binwidth, config_file, index, isData)

    return popt_sig, perr_sig, integral_w, Delta_integral, intrange_syst_std



def getIntegralErrorWithMCMethod(function, p0, p0err, binwidth, config_file, index, isData, bEstimateError=True):

    #xrange_min = config_file[('MC_' if not isData else '') + 'raw_integral_bounds'][index][0]
    #xrange_max = config_file[('MC_' if not isData else '') + 'raw_integral_bounds'][index][-1]

    master_integral_range = [75,105] if config_file['resonance'] == "Z" else [2.7,3.5]
    xrange_min = master_integral_range[0]
    xrange_max = master_integral_range[1]

    def integrand(x):
        return function([x], *p0)[0]
    integral = quad(integrand, xrange_min, xrange_max)
    integral_w = integral[0] / binwidth

    integral_rangesyst_std = getIntegralRangeSystError(integral[0], integrand, xrange_min, xrange_max, binwidth)

    if bEstimateError:
        # error with MC method
        parameters = []
        for i in range(config_file["integral_error_MCmethod_iterations"]): 

            parameters.append(
                [np.random.normal(p0[k], p0err[k]) for k in range(len(p0))]
            )

        integrals_for_Delta = []
        integration_warnings = 0
        for i in range(len(parameters)):
            def integrand_i(x):
                return function([x], *(parameters[i]))[0]

            integration_warning = False
            with warnings.catch_warnings(record=True) as caught_warnings:
                warnings.simplefilter("always", IntegrationWarning)
                integral_i = (quad(integrand_i, xrange_min, xrange_max))[0] / binwidth

                for warning in caught_warnings:
                    if isinstance(warning.message, IntegrationWarning):
                        integration_warning = True
                        integration_warnings += 1

                if  not np.isnan(integral_i) and not np.isinf(integral_i) and not integration_warning:
                    # ignore crazy values
                    if 1 < integral_i < 3*integral_w:
                        integrals_for_Delta.append( integral_i)

            #except IntegrationWarning as e:
            #    print ("Data " if isData else "MC ") + "Fit "+str(index)+"  IntegrationWarning "
        
        if integration_warnings > 0:
            print ("Data " if isData else "MC ") + "Fit "+str(index)+", int delta est: "+str(integration_warnings)+"x  IntegrationWarning "

        #print integral_w
        #print integrals_for_Delta, np.mean(integrals_for_Delta), np.std(integrals_for_Delta)
        Delta_integral = np.std(integrals_for_Delta)
    else:
        Delta_integral = False

    return integral_w, Delta_integral, integral_rangesyst_std






def fit_psip(index, x_in, y_in, y_pm_in, config_file, sigma, Delta_sigma):

    if config_file['ignore_mid_zeros_in_psip'][index]:
        x, y, y_pm = [], [], []
        for i in range(len(x_in)):
            if y_in[i] > 0 or i == 0 or i == len(x_in)-1:
                x.append(x_in[i])
                y.append(y_in[i])
                y_pm.append(y_pm_in[i])
    else:
        x = x_in
        y = y_in
        y_pm = y_pm_in

    binwidth = abs(x[1] - x[0])
    try:
        bounds = ((0,3.6,0.5*sigma),(np.inf,3.770,2*sigma))  
        popt_sig, pcov_sig = curve_fit(gauss, x, y, [100,3.6, sigma], sigma=y_pm, absolute_sigma=True, bounds=bounds,maxfev=5000)
        perr_sig = np.sqrt(np.diag(pcov_sig))

    except Exception as e:
        bounds = ((0,3.6),(np.inf,3.770))  
        popt_sig, pcov_sig = curve_fit(lambda x, A, mass: gauss(x, A, mass, sigma), x, y, [100,3.6], sigma=y_pm, absolute_sigma=True, bounds=bounds,maxfev=5050)
        perr_sig = np.sqrt(np.diag(pcov_sig))

        popt_sig = np.append(popt_sig, [sigma])
        perr_sig = np.append(perr_sig, [Delta_sigma])

    #print popt_sig

    def integrand(x):
        return gauss([x], popt_sig[0],popt_sig[1],popt_sig[2])[0]
    integral = quad(integrand, x[0],x[-1])

    integral_w = integral[0] / binwidth
    Delta_integral = math.sqrt(2*math.pi) * math.sqrt((popt_sig[2]**2) * (perr_sig[0]**2) + (popt_sig[0]**2) * (perr_sig[2]**2)) / binwidth

    return popt_sig, perr_sig, integral_w, Delta_integral





def makePlot(output_dir, config_file, binname, index, isData, x_axis, x_axis_pm, y_axis, y_axis_pm, x_axis_black, y_axis_black, doBkg, chevx, chevparams, chevstats, hist_integral, Delta_hist_integral, linbksub_integral, integral_range_std, gauss_settings, cb_settings, psip_settings):
    
    plt.close()

    chevorder = config_file["chevorder"]

    if isData and index in config_file['data_override_plots_xrange']:
        lower_x = config_file['data_override_plots_xrange'][index][0]
        upper_x = config_file['data_override_plots_xrange'][index][-1]
        if cb_settings['bDraw'] and config_file["bFitCrystalball"] and index in config_file['dualfit_bounds_cb']:
            if config_file['dualfit_bounds_cb'][index][0] > lower_x:
                lower_x = config_file['dualfit_bounds_cb'][index][0]
            if config_file['dualfit_bounds_cb'][index][1] < upper_x:
                upper_x = config_file['dualfit_bounds_cb'][index][1]
        x = np.linspace(lower_x,upper_x,10000)

    elif cb_settings['bDraw'] and config_file["bFitCrystalball"] and index in config_file['dualfit_bounds_cb'] and isData:
        x = np.linspace(config_file['dualfit_bounds_cb'][index][0],config_file['dualfit_bounds_cb'][index][-1],10000)
    elif cb_settings['bDraw'] and config_file["bFitCrystalball"] and index in config_file['MC_dualfit_bounds_cb'] and not isData:
        x = np.linspace(config_file['MC_dualfit_bounds_cb'][index][0],config_file['MC_dualfit_bounds_cb'][index][-1],10000)

    elif cb_settings['bDraw'] and config_file["bFitCrystalball"] and index in config_file['dualfit_bounds_bwg'] and isData:
        x = np.linspace(config_file['dualfit_bounds_bwg'][index][0],config_file['dualfit_bounds_bwg'][index][-1],10000)
    elif cb_settings['bDraw'] and config_file["bFitCrystalball"] and index in config_file['MC_dualfit_bounds_bwg'] and not isData:
        x = np.linspace(config_file['MC_dualfit_bounds_bwg'][index][0],config_file['MC_dualfit_bounds_bwg'][index][-1],10000)


    elif not isData and index in config_file['mc_override_plots_xrange']:
        x = np.linspace(config_file['mc_override_plots_xrange'][index][0],config_file['mc_override_plots_xrange'][index][-1],10000)
    elif isData and doBkg:
        x = np.linspace(config_file['background_fit_bounds'][index][0],config_file['background_fit_bounds'][index][-1],10000)
    elif not isData and doBkg:
        x = np.linspace(config_file['MC_background_fit_bounds'][index][0],config_file['MC_background_fit_bounds'][index][-1],10000)
    elif cb_settings['bDraw'] and not isData:
        x = np.linspace(config_file['MC_peak_fit_bounds_cb'][index][0],config_file['MC_peak_fit_bounds_cb'][index][-1],10000)
    else:
        x = np.linspace(config_file['plots_xrange'][0],config_file['plots_xrange'][-1],10000)

    fig = plt.figure(42,figsize=(15, 10))
    ax = fig.add_subplot(111)

    rawint_text = ("BKsub integral/bw = " if isData or doBkg else "Raw integral/bw = ")+"{:.1f}".format(hist_integral)+u"\u00B1"+("{:.1f}".format(100*Delta_hist_integral / hist_integral) if hist_integral else "?")+"%"
    if linbksub_integral:
        integral_bksyst_Delta = linbksub_integral - hist_integral
        rawint_text = rawint_text + u"\u00B1" +("{:.1f}".format(100*integral_bksyst_Delta / hist_integral) if hist_integral else "?")+"%"
        rawint_text = rawint_text + u"\u00B1" +("{:.1f}".format(100*integral_range_std / hist_integral) if hist_integral else "?")+"%"


    ax.plot(x_axis,y_axis, "sb", label = rawint_text)
    ax.errorbar(x_axis, y_axis, yerr=y_axis_pm, xerr = x_axis_pm, capsize=3, fmt=".", ecolor = "black")
    ax.plot(x_axis_black, y_axis_black, "x", color="black")
    ax.set_xlabel("m [GeV]")
    ax.set_ylabel("Events")
    ax.xaxis.set_label_coords(.9, -.1)
    ax.set_title(binname)
    goodlim_y = ax.get_ylim()

    if isData and index in config_file['data_override_plots_xrange']:
        min_y = 9999
        max_y = -1
        for i in range(len(x_axis)):
            if config_file['data_override_plots_xrange'][index][0] <= x_axis[i] <= config_file['data_override_plots_xrange'][index][-1]:
                if y_axis[i] < min_y:
                    min_y = y_axis[i]
                if y_axis[i] > max_y:
                    max_y = y_axis[i]

        goodlim_y = [min_y*0.9, max_y*1.1]



    if doBkg:
        if len (chevstats[0]) == 0:
            chi2 = "?"
        else:
            chi2 = "{:.1f}".format(chevstats[0][0] / (len(chevx) - (chevorder + 1)))
        ax.plot(x,Chev.chebval(x,chevparams),label=str(chevorder)+r"th order chev, $\chi^{2}$/dof="+chi2)


    if gauss_settings['bDraw']:
        if config_file['bUseBWGauss']:
            gauss_y = np.array(Voigt(x, gauss_settings['popt_sig'][0],gauss_settings['popt_sig'][1],gauss_settings['popt_sig'][2],gauss_settings['popt_sig'][3]) )
        else:
            gauss_y = np.array(gauss(x, gauss_settings['popt_sig'][0],gauss_settings['popt_sig'][1],gauss_settings['popt_sig'][2]) )

        if doBkg:
            gauss_y = np.array(Chev.chebval(x,gauss_settings["chevparams"])) + gauss_y

            ax.plot(x,Chev.chebval(x,gauss_settings["chevparams"]),label=str(chevorder)+r"th order chev for Gauss, $\chi^{2}$/dof="+chi2, color="magenta")

        ax.plot(x,gauss_y,color="orange",label=("Breit-Wigner-" if config_file['bUseBWGauss'] else "") + "Gauss, integral/bw = "+"{:.1f}".format(gauss_settings["integral"])+u"\u00B1"+("{:.1f}".format(100*gauss_settings["Delta_integral"] / gauss_settings["integral"]) if gauss_settings["integral"] else "?")+"%")

    if cb_settings['bDraw']:
        cb_y = np.array(mycrystalball(x, cb_settings['popt_sig'][0],cb_settings['popt_sig'][1],cb_settings['popt_sig'][2],cb_settings['popt_sig'][3],cb_settings['popt_sig'][4])  )

        if doBkg:
            cb_y =    np.array(Chev.chebval(x,chevparams)) + cb_y


        text = "Crystalball"
        if cb_settings['fix_sigma']:
            text = text + r"(fixed $\sigma=$"+"{:.0f}".format(cb_settings['fix_sigma']*1000)+" MeV)"


        if psip_settings['bDraw']:
            psip_y = np.array(gauss(x, psip_settings['popt_sig'][0],psip_settings['popt_sig'][1],psip_settings['popt_sig'][2]) )
            cb_y_plus_psip = cb_y + psip_y

            total_integral = cb_settings["integral"] + psip_settings["integral"]
            total_integral_Delta = cb_settings["Delta_integral"] |qs| psip_settings["Delta_integral"]

            fulltext=text + r"+$\psi$'(Gauss)" + ", integral/bw = "+"{:.1f}".format(total_integral)+u"\u00B1"+("{:.1f}".format(100*total_integral_Delta / total_integral) if total_integral else "?")+"%"
            ax.plot(x,cb_y_plus_psip,color="darkmagenta",label=fulltext)
        
        total_integral = cb_settings["integral"]
        total_integral_Delta = cb_settings["Delta_integral"]

        fulltext = text + ", integral/bw = "+"{:.1f}".format(total_integral)+u"\u00B1"+("{:.1f}".format(100*total_integral_Delta / total_integral) if total_integral else "?")+"%"  
        
        if doBkg:
            total_integral_bksyst_Delta = abs(cb_settings["bksyst_Delta_integral"])
            fulltext = fulltext + u"\u00B1" +("{:.1f}".format(100*total_integral_bksyst_Delta / total_integral) if total_integral else "?")+"%"

        fulltext = fulltext + u"\u00B1"+("{:.1f}".format(100* cb_settings["integral_rangesyst_std"] / total_integral) if total_integral else "?")+"%"  

        ax.plot(x,cb_y,color="crimson",label=fulltext)
        


    if isData and index in config_file['data_override_plots_xrange']:
        ax.set_xlim(config_file['data_override_plots_xrange'][index])
        ax.set_xticks(np.arange(config_file['data_override_plots_xrange'][index][0], config_file['data_override_plots_xrange'][index][1], config_file['xticks_density']), minor=True)
    else:
        ax.set_xlim(config_file['plots_xrange'])
        ax.set_xticks(np.arange(config_file['plots_xrange'][0], config_file['plots_xrange'][1], config_file['xticks_density']), minor=True)


    ax.set_ylim(goodlim_y)
    ax.legend(loc='best', title=("Data" if isData else "MC")+", "+config_file['workingpoint'])
    ax.grid(linestyle = "--")

    vlinecolor = 'gray' if index in config_file[("MC_" if not isData else "")+'dualfit_bounds_cb'] else 'k'
    master_integral_range = [75,105] if config_file['resonance'] == "Z" else [2.7,3.5]
    #ax.vlines(config_file[("MC_" if not isData else "")+"raw_integral_bounds"][index][0], goodlim_y[0],goodlim_y[1], color=vlinecolor)
    #ax.vlines(config_file[("MC_" if not isData else "")+"raw_integral_bounds"][index][1], goodlim_y[0],goodlim_y[1], color=vlinecolor)
    ax.vlines(master_integral_range[0], goodlim_y[0],goodlim_y[1], color=vlinecolor)
    ax.vlines(master_integral_range[1], goodlim_y[0],goodlim_y[1], color=vlinecolor)



    if doBkg:
        if not os.path.exists(output_dir+"/fits/chev"+str(chevorder)):
            os.makedirs(output_dir+"/fits/chev"+str(chevorder))
        plt.savefig(output_dir+"/fits/chev"+str(chevorder)+"/figure"+str(index)+".png")
    else:
        if not os.path.exists(output_dir+"/fits/"):
            os.makedirs(output_dir+"/fits/")
        plt.savefig(output_dir+"/fits/figure"+str(index)+".png")
    plt.close()



# plotdata contains: data_x_axis, data_x_axis_pm, data_y_axis, data_y_axis_pm, data_x_axis_black, data_y_axis_black, data_chevx, data_chevparams, data_chevstats, data_chevorder, data_hist_integral, data_Delta_hist_integral, data_gauss_settings, data_cb_settings, data_psip_settings
def makeCommonDataMCFitPlot(output_dir, config_file, binname, index, data_plotinfo, mc_plotinfo):
    
    plt.close()

    mc_doBkg = config_file['fit_background_on_MC']



    if index in config_file['data_override_plots_xrange']:
        data_x = np.linspace(config_file['data_override_plots_xrange'][index][0],config_file['data_override_plots_xrange'][index][-1],10000)
    elif index in config_file['dualfit_bounds_cb']:
        data_x = np.linspace(config_file['dualfit_bounds_cb'][index][0],config_file['dualfit_bounds_cb'][index][-1],10000)
    else:
        data_x = np.linspace(useBigger(config_file['background_fit_bounds'][index][0], config_file['plots_xrange'][0]), useSmaller(config_file['background_fit_bounds'][index][-1], config_file['plots_xrange'][-1]),10000)

    if index in config_file['mc_override_plots_xrange']:
        mc_x = np.linspace(config_file['mc_override_plots_xrange'][index][0],config_file['mc_override_plots_xrange'][index][-1],10000)
    elif mc_doBkg and index in config_file['MC_dualfit_bounds_cb']:
        mc_x = np.linspace(config_file['MC_dualfit_bounds_cb'][index][0],config_file['MC_dualfit_bounds_cb'][index][-1],10000)
    elif mc_doBkg:
        mc_x = np.linspace(useBigger(config_file['MC_background_fit_bounds'][index][0], config_file['plots_xrange'][0]),useSmaller(config_file['MC_background_fit_bounds'][index][-1],config_file['plots_xrange'][-1]),10000)
    else:
        mc_x = np.linspace(config_file['plots_xrange'][0],config_file['plots_xrange'][-1],10000)


    fig = plt.figure(42,figsize=(15, 10))
    ax = fig.add_subplot(111)
    ax.plot(data_plotinfo["x_axis"], data_plotinfo["y_axis"], "sb", label = "DATA, BKsub integral/bw = "+"{:.1f}".format(data_plotinfo["hist_integral"])+u"\u00B1"+("{:.1f}".format(100*data_plotinfo["Delta_hist_integral"] / data_plotinfo["hist_integral"]) if data_plotinfo["hist_integral"] else "?")+"%")
    ax.errorbar(data_plotinfo["x_axis"], data_plotinfo["y_axis"], yerr=data_plotinfo["y_axis_pm"], xerr = data_plotinfo["x_axis_pm"], capsize=3, fmt=".", ecolor = "black")

    ax.plot(mc_plotinfo["x_axis"], mc_plotinfo["y_axis"], color="darkred", marker="s", ls="", label = "MC, "+("BKsub integral/bw = " if mc_doBkg else "Raw integral/bw = ")+"{:.1f}".format(mc_plotinfo["hist_integral"])+u"\u00B1"+("{:.1f}".format(100*mc_plotinfo["Delta_hist_integral"] / mc_plotinfo["hist_integral"]) if mc_plotinfo["hist_integral"] else "?")+"%")
    ax.errorbar(mc_plotinfo["x_axis"], mc_plotinfo["y_axis"], yerr=mc_plotinfo["y_axis_pm"], xerr = mc_plotinfo["x_axis_pm"], capsize=3, fmt=".", ecolor = "black", color="red")


    ax.plot(data_plotinfo["x_axis_black"], data_plotinfo["y_axis_black"], "x", color="black")
    ax.plot(mc_plotinfo["x_axis_black"], mc_plotinfo["y_axis_black"], "x", color="black")

    ax.set_xlabel("m [GeV]")
    ax.set_ylabel("Events")
    ax.xaxis.set_label_coords(.9, -.1)
    ax.set_title(binname)
    goodlim_y = ax.get_ylim()



    if len (data_plotinfo["chevstats"][0]) == 0:
        chi2 = "?"
    else:
        chi2 = "{:.1f}".format(data_plotinfo["chevstats"][0][0] / (len(data_plotinfo["chevx"]) - (config_file["chevorder"] + 1)))
    ax.plot(data_x,Chev.chebval(data_x,data_plotinfo["chevparams"]),label=str(config_file["chevorder"])+r"th order chev, data, $\chi^{2}$/dof="+chi2)
    if mc_doBkg:
        if len (mc_plotinfo["chevstats"][0]) == 0:
            chi2 = "?"
        else:
            chi2 = "{:.1f}".format(mc_plotinfo["chevstats"][0][0] / (len(mc_plotinfo["chevx"]) - (config_file["chevorder"] + 1)))
        ax.plot(mc_x,Chev.chebval(mc_x,mc_plotinfo["chevparams"]),label=str(config_file["chevorder"])+r"th order chev, MC, $\chi^{2}$/dof="+chi2)


    if data_plotinfo["gauss_settings"]['bDraw']:
        if config_file['bUseBWGauss']:
            gauss_y = np.array(Voigt(data_x, data_plotinfo["gauss_settings"]['popt_sig'][0],data_plotinfo["gauss_settings"]['popt_sig'][1],data_plotinfo["gauss_settings"]['popt_sig'][2],data_plotinfo["gauss_settings"]['popt_sig'][3]) )
        else:
            gauss_y = np.array(gauss(data_x, data_plotinfo["gauss_settings"]['popt_sig'][0],data_plotinfo["gauss_settings"]['popt_sig'][1],data_plotinfo["gauss_settings"]['popt_sig'][2]) )

        gauss_y = np.array(Chev.chebval(data_x,data_plotinfo["chevparams"])) + gauss_y

        ax.plot(data_x,gauss_y,color="purple",label=("Breit-Wigner-" if config_file['bUseBWGauss'] else "") + "Gauss, integral/bw = "+"{:.1f}".format(data_plotinfo["gauss_settings"]["integral"])+u"\u00B1"+("{:.1f}".format(100*data_plotinfo["gauss_settings"]["Delta_integral"] / data_plotinfo["gauss_settings"]["integral"]) if data_plotinfo["gauss_settings"]["integral"] else "?")+"%")
    

    if mc_plotinfo["gauss_settings"]['bDraw']:
        if config_file['bUseBWGauss']:
            gauss_y = np.array(Voigt(mc_x, mc_plotinfo["gauss_settings"]['popt_sig'][0],mc_plotinfo["gauss_settings"]['popt_sig'][1],mc_plotinfo["gauss_settings"]['popt_sig'][2],mc_plotinfo["gauss_settings"]['popt_sig'][3]) )
        else:
            gauss_y = np.array(gauss(mc_x, mc_plotinfo["gauss_settings"]['popt_sig'][0],mc_plotinfo["gauss_settings"]['popt_sig'][1],mc_plotinfo["gauss_settings"]['popt_sig'][2]) )

        if mc_doBkg:
            gauss_y = np.array(Chev.chebval(mc_x,mc_plotinfo["chevparams"])) + gauss_y

        ax.plot(mc_x,gauss_y,color="magenta",label=("Breit-Wigner-" if config_file['bUseBWGauss'] else "") + "Gauss, integral/bw = "+"{:.1f}".format(mc_plotinfo["gauss_settings"]["integral"])+u"\u00B1"+("{:.1f}".format(100*mc_plotinfo["gauss_settings"]["Delta_integral"] / mc_plotinfo["gauss_settings"]["integral"]) if mc_plotinfo["gauss_settings"]["integral"] else "?")+"%")

    if data_plotinfo["cb_settings"]['bDraw']:
        cb_y = np.array(mycrystalball(data_x, data_plotinfo["cb_settings"]['popt_sig'][0],data_plotinfo["cb_settings"]['popt_sig'][1],data_plotinfo["cb_settings"]['popt_sig'][2],data_plotinfo["cb_settings"]['popt_sig'][3],data_plotinfo["cb_settings"]['popt_sig'][4])  )

        if data_plotinfo["psip_settings"]['bDraw']:
            psip_y = np.array(gauss(data_x, data_plotinfo["psip_settings"]['popt_sig'][0],data_plotinfo["psip_settings"]['popt_sig'][1],data_plotinfo["psip_settings"]['popt_sig'][2]) )
            cb_y = cb_y + psip_y
            total_integral = data_plotinfo["cb_settings"]["integral"] + data_plotinfo["psip_settings"]["integral"]
            total_integral_Delta = data_plotinfo["cb_settings"]["Delta_integral"] |qs| data_plotinfo["psip_settings"]["Delta_integral"]
        else:
            total_integral = data_plotinfo["cb_settings"]["integral"]
            total_integral_Delta = data_plotinfo["cb_settings"]["Delta_integral"]

        cb_y =    np.array(Chev.chebval(data_x,data_plotinfo["chevparams"])) + cb_y

        text = "Crystalball"
        if data_plotinfo["psip_settings"]['bDraw']:
            text = text + r"+$\psi$'(Gauss)"
        if data_plotinfo["cb_settings"]['fix_sigma']:
            text = text + r"(fixed $\sigma=$"+"{:.0f}".format(data_plotinfo["cb_settings"]['fix_sigma']*1000)+" MeV)"
        ax.plot(data_x,cb_y,color="orange",label=text+", integral/bw = "+"{:.1f}".format(total_integral)+u"\u00B1"+("{:.1f}".format(100*total_integral_Delta / total_integral) if total_integral else "?")+"%")


    if mc_plotinfo["cb_settings"]['bDraw']:
        cb_y = np.array(mycrystalball(mc_x, mc_plotinfo["cb_settings"]['popt_sig'][0],mc_plotinfo["cb_settings"]['popt_sig'][1],mc_plotinfo["cb_settings"]['popt_sig'][2],mc_plotinfo["cb_settings"]['popt_sig'][3],mc_plotinfo["cb_settings"]['popt_sig'][4])  )

        if mc_plotinfo["psip_settings"]['bDraw']:
            psip_y = np.array(gauss(mc_x, mc_plotinfo["psip_settings"]['popt_sig'][0],mc_plotinfo["psip_settings"]['popt_sig'][1],mc_plotinfo["psip_settings"]['popt_sig'][2]) )
            cb_y = cb_y + psip_y
            total_integral = mc_plotinfo["cb_settings"]["integral"] + mc_plotinfo["psip_settings"]["integral"]
            total_integral_Delta = mc_plotinfo["cb_settings"]["Delta_integral"] |qs| mc_plotinfo["psip_settings"]["Delta_integral"]
        else:
            total_integral = mc_plotinfo["cb_settings"]["integral"]
            total_integral_Delta = mc_plotinfo["cb_settings"]["Delta_integral"]

        if mc_doBkg:
            cb_y =    np.array(Chev.chebval(mc_x,mc_plotinfo["chevparams"])) + cb_y

        text = "Crystalball"
        if mc_plotinfo["psip_settings"]['bDraw']:
            text = text + r"+$\psi$'(Gauss)"
        if mc_plotinfo["cb_settings"]['fix_sigma']:
            text = text + r"(fixed $\sigma=$"+"{:.0f}".format(mc_plotinfo["cb_settings"]['fix_sigma']*1000)+" MeV)"
        ax.plot(mc_x,cb_y,color="crimson",label=text+", integral/bw = "+"{:.1f}".format(total_integral)+u"\u00B1"+("{:.1f}".format(100*total_integral_Delta / total_integral) if total_integral else "?")+"%")

    if config_file['resonance'] == "Z":
        xmin = useSmaller(data_x[0], mc_x[0])
        xmax = useBigger(data_x[-1], mc_x[-1])
        ax.set_xlim([60 if 60 < xmin else xmin,120 if 120 > xmax else xmax])
    ax.set_ylim(goodlim_y)
    ax.legend(loc='best', title=config_file['workingpoint'])
    ax.grid(linestyle = "--")


    if not os.path.exists(output_dir+"/fits/chev"+str(config_file["chevorder"])):
        os.makedirs(output_dir+"/fits/chev"+str(config_file["chevorder"]))
    plt.savefig(output_dir+"/fits/chev"+str(config_file["chevorder"])+"/figure"+str(index)+".png")
    plt.close()

def makeCommonSubtractedDataMCFitPlot(output_dir, config_file, binname, index, data_plotinfo, mc_plotinfo):
    
    plt.close()

    mc_doBkg = config_file['fit_background_on_MC']

    if index in config_file['data_override_plots_xrange']:
        data_x = np.linspace(config_file['data_override_plots_xrange'][index][0],config_file['data_override_plots_xrange'][index][-1],10000)
    elif index in config_file['dualfit_bounds_cb']:
        data_x = np.linspace(config_file['dualfit_bounds_cb'][index][0],config_file['dualfit_bounds_cb'][index][-1],10000)
    else:
        data_x = np.linspace(useBigger(config_file['background_fit_bounds'][index][0], config_file['plots_xrange'][0]), useSmaller(config_file['background_fit_bounds'][index][-1], config_file['plots_xrange'][-1]),10000)

    if index in config_file['mc_override_plots_xrange']:
        mc_x = np.linspace(config_file['mc_override_plots_xrange'][index][0],config_file['mc_override_plots_xrange'][index][-1],10000)
    elif mc_doBkg and index in config_file['MC_dualfit_bounds_cb']:
        mc_x = np.linspace(config_file['MC_dualfit_bounds_cb'][index][0],config_file['MC_dualfit_bounds_cb'][index][-1],10000)
    elif mc_doBkg:
        mc_x = np.linspace(useBigger(config_file['MC_background_fit_bounds'][index][0], config_file['plots_xrange'][0]),useSmaller(config_file['MC_background_fit_bounds'][index][-1],config_file['plots_xrange'][-1]),10000)
    else:
        mc_x = np.linspace(config_file['plots_xrange'][0],config_file['plots_xrange'][-1],10000)

    data_y_axis_sub, data_y_axis_sub_pm = subtract_background(data_plotinfo["x_axis"], data_plotinfo["y_axis"], data_plotinfo["y_axis_pm"], data_plotinfo["chevparams"], True, data_plotinfo["chev_std_devs"])
    if mc_doBkg:
        mc_y_axis_sub, mc_y_axis_sub_pm = subtract_background(mc_plotinfo["x_axis"], mc_plotinfo["y_axis"], mc_plotinfo["y_axis_pm"], mc_plotinfo["chevparams"], True, mc_plotinfo["chev_std_devs"])
    else:
        mc_y_axis_sub, mc_y_axis_sub_pm = mc_plotinfo["y_axis"], mc_plotinfo["y_axis_pm"]

    # TODO: Show <0 as black points?
        
    data_x_axis_sub = data_plotinfo["x_axis"]
    data_x_axis_err_sub = data_plotinfo["x_axis_pm"]
    mc_x_axis_sub = mc_plotinfo["x_axis"]
    mc_x_axis_err_sub = mc_plotinfo["x_axis_pm"]
    if config_file['background_fit_bounds'][index][0] > config_file['plots_xrange'][0]:
        # don't show points outside bk fit range
        for i in range(len(data_x_axis_sub)):
            if data_x_axis_sub[i] >= config_file['background_fit_bounds'][index][0]:
                data_x_axis_sub = data_x_axis_sub[i:]
                data_x_axis_err_sub = data_x_axis_err_sub[i:]
                data_y_axis_sub = data_y_axis_sub[i:]
                data_y_axis_sub_pm = data_y_axis_sub_pm[i:]
                break

        for i in range(len(mc_x_axis_sub)):
            if mc_x_axis_sub[i] >= config_file['background_fit_bounds'][index][0]:
                mc_x_axis_sub = mc_x_axis_sub[i:]
                mc_x_axis_err_sub = mc_x_axis_err_sub[i:]
                mc_y_axis_sub = mc_y_axis_sub[i:]
                mc_y_axis_sub_pm = mc_y_axis_sub_pm[i:]
                break

    if config_file['background_fit_bounds'][index][-1] < config_file['plots_xrange'][-1]:
        for i in range(len(data_x_axis_sub) - 1, -1, -1):
            if data_x_axis_sub[i] <= config_file['background_fit_bounds'][index][0]:
                data_x_axis_sub = data_x_axis_sub[:(i+1)]
                data_x_axis_err_sub = data_x_axis_err_sub[:(i+1)]
                data_y_axis_sub = data_y_axis_sub[:(i+1)]
                data_y_axis_sub_pm = data_y_axis_sub_pm[:(i+1)]
                break

        for i in range(len(mc_x_axis_sub) - 1, -1, -1):
            if mc_x_axis_sub[i] <= config_file['background_fit_bounds'][index][0]:
                mc_x_axis_sub = mc_x_axis_sub[:(i+1)]
                mc_x_axis_err_sub = mc_x_axis_err_sub[:(i+1)]
                mc_y_axis_sub = mc_y_axis_sub[:(i+1)]
                mc_y_axis_sub_pm = mc_y_axis_sub_pm[:(i+1)]
                break
    



    fig = plt.figure(42,figsize=(15, 10))
    ax = fig.add_subplot(111)
    ax.plot(data_x_axis_sub, data_y_axis_sub, "sb", label = "DATA, BKsub integral/bw = "+"{:.1f}".format(data_plotinfo["hist_integral"])+u"\u00B1"+("{:.1f}".format(100*data_plotinfo["Delta_hist_integral"] / data_plotinfo["hist_integral"]) if data_plotinfo["hist_integral"] else "?")+"%")
    ax.errorbar(data_x_axis_sub, data_y_axis_sub, data_y_axis_sub_pm, xerr = data_x_axis_err_sub, capsize=3, fmt=".", ecolor = "black")

    ax.plot(mc_x_axis_sub,mc_y_axis_sub, color="darkred", marker="s", ls="", label = "MC, "+("BKsub integral/bw = " if mc_doBkg else "Raw integral/bw = ")+"{:.1f}".format(mc_plotinfo["hist_integral"])+u"\u00B1"+("{:.1f}".format(100*mc_plotinfo["Delta_hist_integral"] / mc_plotinfo["hist_integral"]) if mc_plotinfo["hist_integral"] else "?")+"%")
    ax.errorbar(mc_x_axis_sub, mc_y_axis_sub, mc_y_axis_sub_pm, xerr = mc_x_axis_err_sub, capsize=3, fmt=".", ecolor = "black", color="red")


    ax.plot(data_plotinfo["x_axis_black"], data_plotinfo["y_axis_black"], "x", color="black")
    ax.plot(mc_plotinfo["x_axis_black"], mc_plotinfo["y_axis_black"], "x", color="black")

    ax.set_xlabel("m [GeV]")
    ax.set_ylabel("Events - BK")
    ax.xaxis.set_label_coords(.9, -.1)
    ax.set_title(binname)
    goodlim_y = ax.get_ylim()



    if data_plotinfo["gauss_settings"]['bDraw']:
        if config_file['bUseBWGauss']:
            gauss_y = np.array(Voigt(data_x, data_plotinfo["gauss_settings"]['popt_sig'][0],data_plotinfo["gauss_settings"]['popt_sig'][1],data_plotinfo["gauss_settings"]['popt_sig'][2],data_plotinfo["gauss_settings"]['popt_sig'][3]) )
        else:
            gauss_y = np.array(gauss(data_x, data_plotinfo["gauss_settings"]['popt_sig'][0],data_plotinfo["gauss_settings"]['popt_sig'][1],data_plotinfo["gauss_settings"]['popt_sig'][2]) )

        ax.plot(data_x,gauss_y,color="purple",label=("Breit-Wigner-" if config_file['bUseBWGauss'] else "") + "Gauss, integral/bw = "+"{:.1f}".format(data_plotinfo["gauss_settings"]["integral"])+u"\u00B1"+("{:.1f}".format(100*data_plotinfo["gauss_settings"]["Delta_integral"] / data_plotinfo["gauss_settings"]["integral"]) if data_plotinfo["gauss_settings"]["integral"] else "?")+"%")
    

    if mc_plotinfo["gauss_settings"]['bDraw']:
        if config_file['bUseBWGauss']:
            gauss_y = np.array(Voigt(mc_x, mc_plotinfo["gauss_settings"]['popt_sig'][0],mc_plotinfo["gauss_settings"]['popt_sig'][1],mc_plotinfo["gauss_settings"]['popt_sig'][2],mc_plotinfo["gauss_settings"]['popt_sig'][3]) )
        else:
            gauss_y = np.array(gauss(mc_x, mc_plotinfo["gauss_settings"]['popt_sig'][0],mc_plotinfo["gauss_settings"]['popt_sig'][1],mc_plotinfo["gauss_settings"]['popt_sig'][2]) )

        ax.plot(mc_x,gauss_y,color="magenta",label=("Breit-Wigner-" if config_file['bUseBWGauss'] else "") + "Gauss, integral/bw = "+"{:.1f}".format(mc_plotinfo["gauss_settings"]["integral"])+u"\u00B1"+("{:.1f}".format(100*mc_plotinfo["gauss_settings"]["Delta_integral"] / mc_plotinfo["gauss_settings"]["integral"]) if mc_plotinfo["gauss_settings"]["integral"] else "?")+"%")

    if data_plotinfo["cb_settings"]['bDraw']:
        cb_y = np.array(mycrystalball(data_x, data_plotinfo["cb_settings"]['popt_sig'][0],data_plotinfo["cb_settings"]['popt_sig'][1],data_plotinfo["cb_settings"]['popt_sig'][2],data_plotinfo["cb_settings"]['popt_sig'][3],data_plotinfo["cb_settings"]['popt_sig'][4])  )

        if data_plotinfo["psip_settings"]['bDraw']:
            psip_y = np.array(gauss(data_x, data_plotinfo["psip_settings"]['popt_sig'][0],data_plotinfo["psip_settings"]['popt_sig'][1],data_plotinfo["psip_settings"]['popt_sig'][2]) )
            cb_y = cb_y + psip_y
            total_integral = data_plotinfo["cb_settings"]["integral"] + data_plotinfo["psip_settings"]["integral"]
            total_integral_Delta = data_plotinfo["cb_settings"]["Delta_integral"] |qs| data_plotinfo["psip_settings"]["Delta_integral"]
        else:
            total_integral = data_plotinfo["cb_settings"]["integral"]
            total_integral_Delta = data_plotinfo["cb_settings"]["Delta_integral"]

        text = "Crystalball"
        if data_plotinfo["psip_settings"]['bDraw']:
            text = text + r"+$\psi$'(Gauss)"
        if data_plotinfo["cb_settings"]['fix_sigma']:
            text = text + r"(fixed $\sigma=$"+"{:.0f}".format(data_plotinfo["cb_settings"]['fix_sigma']*1000)+" MeV)"
        ax.plot(data_x,cb_y,color="orange",label=text+", integral/bw = "+"{:.1f}".format(total_integral)+u"\u00B1"+("{:.1f}".format(100*total_integral_Delta / total_integral) if total_integral else "?")+"%")


    if mc_plotinfo["cb_settings"]['bDraw']:
        cb_y = np.array(mycrystalball(mc_x, mc_plotinfo["cb_settings"]['popt_sig'][0],mc_plotinfo["cb_settings"]['popt_sig'][1],mc_plotinfo["cb_settings"]['popt_sig'][2],mc_plotinfo["cb_settings"]['popt_sig'][3],mc_plotinfo["cb_settings"]['popt_sig'][4])  )

        if mc_plotinfo["psip_settings"]['bDraw']:
            psip_y = np.array(gauss(mc_x, mc_plotinfo["psip_settings"]['popt_sig'][0],mc_plotinfo["psip_settings"]['popt_sig'][1],mc_plotinfo["psip_settings"]['popt_sig'][2]) )
            cb_y = cb_y + psip_y
            total_integral = mc_plotinfo["cb_settings"]["integral"] + mc_plotinfo["psip_settings"]["integral"]
            total_integral_Delta = mc_plotinfo["cb_settings"]["Delta_integral"] |qs| mc_plotinfo["psip_settings"]["Delta_integral"]
        else:
            total_integral = mc_plotinfo["cb_settings"]["integral"]
            total_integral_Delta = mc_plotinfo["cb_settings"]["Delta_integral"]

        text = "Crystalball"
        if mc_plotinfo["psip_settings"]['bDraw']:
            text = text + r"+$\psi$'(Gauss)"
        if mc_plotinfo["cb_settings"]['fix_sigma']:
            text = text + r"(fixed $\sigma=$"+"{:.0f}".format(mc_plotinfo["cb_settings"]['fix_sigma']*1000)+" MeV)"
        ax.plot(mc_x,cb_y,color="crimson",label=text+", integral/bw = "+"{:.1f}".format(total_integral)+u"\u00B1"+("{:.1f}".format(100*total_integral_Delta / total_integral) if total_integral else "?")+"%")



    if config_file['resonance'] == "Z":
        xmin = useSmaller(data_x[0], mc_x[0])
        xmax = useBigger(data_x[-1], mc_x[-1])
        ax.set_xlim([60 if 60 < xmin else xmin,120 if 120 > xmax else xmax])
    ax.set_ylim(goodlim_y)
    ax.legend(loc='best', title=config_file['workingpoint'])
    ax.grid(linestyle = "--")


    if not os.path.exists(output_dir+"/fits/chev"+str(config_file["chevorder"])):
        os.makedirs(output_dir+"/fits/chev"+str(config_file["chevorder"]))
    plt.savefig(output_dir+"/fits/chev"+str(config_file["chevorder"])+"/figure"+str(index)+"_bksub.png")
    plt.close()



def TnPEfficiency(mutex, shared_data_raw_integrals, shared_data_raw_integral_sigmas, shared_data_raw_integral_sigmas_rangesyst, shared_data_raw_integral_sigmas_bksyst, shared_data_gauss_integrals, shared_data_gauss_integral_sigmas, shared_data_gauss_integral_sigmas_rangesyst, shared_data_gauss_integral_sigmas_bksyst, shared_data_integrals_cb, shared_data_integral_cb_sigmas, shared_data_integral_cb_sigmas_rangesyst, shared_data_integral_cb_sigmas_bksyst, shared_mc_raw_integrals, shared_mc_raw_integral_sigmas, shared_mc_raw_integral_sigmas_rangesyst, shared_mc_raw_integral_sigmas_bksyst, shared_mc_gauss_integrals, shared_mc_gauss_integral_sigmas, shared_mc_gauss_integral_sigmas_rangesyst, shared_mc_gauss_integral_sigmas_bksyst, shared_mc_integrals_cb, shared_mc_integral_cb_sigmas, shared_mc_integral_cb_sigmas_rangesyst, shared_mc_integral_cb_sigmas_bksyst, shared_mc_alt_raw_integrals, shared_mc_alt_raw_integral_sigmas, shared_mc_alt_raw_integral_sigmas_rangesyst, shared_mc_alt_raw_integral_sigmas_bksyst, shared_mc_alt_gauss_integrals, shared_mc_alt_gauss_integral_sigmas, shared_mc_alt_gauss_integral_sigmas_rangesyst, shared_mc_alt_gauss_integral_sigmas_bksyst, shared_mc_alt_integrals_cb, shared_mc_alt_integral_cb_sigmas, shared_mc_alt_integral_cb_sigmas_rangesyst, shared_mc_alt_integral_cb_sigmas_bksyst,
                       datafile, mcfile, mcfile_alt, bins, binnames, indexlist, binindex_range_to_fit, workdir, config_file, skipped_lowpt_bins):


    first_one = True
    for index in indexlist:
        fit_in_this_bin = binindex_range_to_fit[0] <= index < binindex_range_to_fit[1]

        use_fixed_sigma_data = False
        use_fixed_sigma_data_Delta = False
        if index in config_file['fix_sigma_from_previous_bin_data']:
            if first_one:
                print "Error: Use different number of threads with fixed sigma from previous bin in bin "+str(index)
            else:
                print "Fixing sigma for bin "+str(index)+" to "+str(sigma_cb_fit)
                use_fixed_sigma_data = sigma_cb_fit
                use_fixed_sigma_data_Delta = use_fixed_sigma_data_Delta

        first_one = False

        # for JPSi:  raw, CB
        # for Z:  raw, BWG, CB
        data_plotinfo, data_raw_integral, data_raw_integral_sigma, data_gauss_integral, data_gauss_integral_sigma, data_integral_cb, data_integral_cb_sigma, sigma_cb_fit, Delta_sigma_cb_fit = TnPEfficiencyInBin(True, datafile, workdir, "dataplots/", bins[index], binnames[index], index, fit_in_this_bin, config_file, use_fixed_sigma_data, use_fixed_sigma_data_Delta)

        


        if index in config_file['fix_sigma_from_previous_bin_data']:
            print "Fitted sigma after fixing it down: "+str(sigma_cb_fit)

        # for JPsi: raw, fixed width fit CB
        # for Z: raw, BWG, CB
        mc_plotinfo, mc_raw_integral, mc_raw_integral_sigma, mc_gauss_integral, mc_gauss_integral_sigma, mc_integral_cb, mc_integral_cb_sigma, _, _ = TnPEfficiencyInBin(False, mcfile, workdir, "mcplots/", bins[index], binnames[index], index, fit_in_this_bin, config_file, sigma_cb_fit if config_file['resonance'] == "JPsi" else False, Delta_sigma_cb_fit if config_file['resonance'] == "JPsi" else False)

        if fit_in_this_bin:
            makeCommonDataMCFitPlot(workdir+"/commonplots/", config_file, binnames[index], index, data_plotinfo, mc_plotinfo)
            makeCommonSubtractedDataMCFitPlot(workdir+"/commonplots/", config_file, binnames[index], index, data_plotinfo, mc_plotinfo)

            with mutex:
                shared_data_raw_integrals[index - 4*skipped_lowpt_bins] = data_raw_integral
                shared_data_raw_integral_sigmas[index - 4*skipped_lowpt_bins] = data_raw_integral_sigma
                shared_data_raw_integral_sigmas_bksyst[index - 4*skipped_lowpt_bins] = data_plotinfo["linkb_hist_integral_Delta"]
                shared_data_raw_integral_sigmas_rangesyst[index - 4*skipped_lowpt_bins] = data_plotinfo["hist_integral_rangesyst_std"]

                shared_data_gauss_integrals[index - 4*skipped_lowpt_bins] = data_gauss_integral
                shared_data_gauss_integral_sigmas[index - 4*skipped_lowpt_bins] = data_gauss_integral_sigma
                shared_data_gauss_integral_sigmas_bksyst[index - 4*skipped_lowpt_bins] = data_plotinfo["gauss_settings"]["bksyst_Delta_integral"] if "bksyst_Delta_integral" in data_plotinfo["gauss_settings"] else 0
                shared_data_gauss_integral_sigmas_rangesyst[index - 4*skipped_lowpt_bins] = data_plotinfo["gauss_settings"]["integral_rangesyst_std"] if "integral_rangesyst_std" in data_plotinfo["gauss_settings"] else 0

                if data_plotinfo["psip_settings"]["bDraw"] and config_file["include_psip_result"]: 
                    data_integral_cb += data_plotinfo["psip_settings"]["integral"]
                    data_integral_cb_sigma = data_integral_cb_sigma |qs| data_plotinfo["psip_settings"]["Delta_integral"]

                shared_data_integrals_cb[index - 4*skipped_lowpt_bins] = data_integral_cb 
                shared_data_integral_cb_sigmas[index - 4*skipped_lowpt_bins] = data_integral_cb_sigma
                shared_data_integral_cb_sigmas_bksyst[index - 4*skipped_lowpt_bins] = data_plotinfo["cb_settings"]["bksyst_Delta_integral"] if "bksyst_Delta_integral" in data_plotinfo["cb_settings"] else 0
                shared_data_integral_cb_sigmas_rangesyst[index - 4*skipped_lowpt_bins] = data_plotinfo["cb_settings"]["integral_rangesyst_std"] if "integral_rangesyst_std" in data_plotinfo["cb_settings"] else 0

                shared_mc_raw_integrals[index - 4*skipped_lowpt_bins] = mc_raw_integral
                shared_mc_raw_integral_sigmas[index - 4*skipped_lowpt_bins] = mc_raw_integral_sigma
                shared_mc_raw_integral_sigmas_bksyst[index - 4*skipped_lowpt_bins] = mc_plotinfo["linkb_hist_integral_Delta"]
                shared_mc_raw_integral_sigmas_rangesyst[index - 4*skipped_lowpt_bins] = mc_plotinfo["hist_integral_rangesyst_std"]

                shared_mc_gauss_integrals[index - 4*skipped_lowpt_bins] = mc_gauss_integral
                shared_mc_gauss_integral_sigmas[index - 4*skipped_lowpt_bins] = mc_gauss_integral_sigma
                shared_mc_gauss_integral_sigmas_bksyst[index - 4*skipped_lowpt_bins] = mc_plotinfo["gauss_settings"]["bksyst_Delta_integral"] if "bksyst_Delta_integral" in mc_plotinfo["gauss_settings"] else 0
                shared_mc_gauss_integral_sigmas_rangesyst[index - 4*skipped_lowpt_bins] = mc_plotinfo["gauss_settings"]["integral_rangesyst_std"] if "integral_rangesyst_std" in mc_plotinfo["gauss_settings"] else 0

                if mc_plotinfo["psip_settings"]["bDraw"] and config_file["include_psip_result"]: 
                    mc_integral_cb += mc_plotinfo["psip_settings"]["integral"]
                    mc_integral_cb_sigma = mc_integral_cb_sigma |qs| mc_plotinfo["psip_settings"]["Delta_integral"]

                shared_mc_integrals_cb[index - 4*skipped_lowpt_bins] = mc_integral_cb
                shared_mc_integral_cb_sigmas[index - 4*skipped_lowpt_bins] = mc_integral_cb_sigma
                shared_mc_integral_cb_sigmas_bksyst[index - 4*skipped_lowpt_bins] = mc_plotinfo["cb_settings"]["bksyst_Delta_integral"] if "bksyst_Delta_integral" in mc_plotinfo["cb_settings"] else 0
                shared_mc_integral_cb_sigmas_rangesyst[index - 4*skipped_lowpt_bins] = mc_plotinfo["cb_settings"]["integral_rangesyst_std"] if "integral_rangesyst_std" in mc_plotinfo["cb_settings"] else 0
 
        if mcfile_alt:
            mc_alt_plotinfo, mc_alt_raw_integral, mc_alt_raw_integral_sigma, mc_alt_gauss_integral, mc_alt_gauss_integral_sigma, mc_alt_integral_cb, mc_alt_integral_cb_sigma,  _, _ = TnPEfficiencyInBin(False, mcfile_alt, workdir, "mc_alt_plots/", bins[index], binnames[index], index, fit_in_this_bin, config_file, sigma_cb_fit if config_file['resonance'] == "JPsi" else False, Delta_sigma_cb_fit if config_file['resonance'] == "JPsi" else False)



            if fit_in_this_bin:
                makeCommonDataMCFitPlot(workdir+"/commonplots_alt/", config_file, binnames[index], index, data_plotinfo, mc_alt_plotinfo)
                makeCommonSubtractedDataMCFitPlot(workdir+"/commonplots_alt/", config_file, binnames[index], index, data_plotinfo, mc_alt_plotinfo)

                with mutex:
                    shared_mc_alt_raw_integrals[index - 4*skipped_lowpt_bins] = mc_alt_raw_integral
                    shared_mc_alt_raw_integral_sigmas[index - 4*skipped_lowpt_bins] = mc_alt_raw_integral_sigma
                    shared_mc_alt_raw_integral_sigmas_bksyst[index - 4*skipped_lowpt_bins] = mc_alt_plotinfo["linkb_hist_integral_Delta"]
                    shared_mc_alt_raw_integral_sigmas_rangesyst[index - 4*skipped_lowpt_bins] = mc_alt_plotinfo["hist_integral_rangesyst_std"]

                    shared_mc_alt_gauss_integrals[index - 4*skipped_lowpt_bins] = mc_alt_gauss_integral
                    shared_mc_alt_gauss_integral_sigmas[index - 4*skipped_lowpt_bins] = mc_alt_gauss_integral_sigma
                    shared_mc_alt_gauss_integral_sigmas_bksyst[index - 4*skipped_lowpt_bins] = mc_alt_plotinfo["gauss_settings"]["bksyst_Delta_integral"] if "bksyst_Delta_integral" in mc_alt_plotinfo["gauss_settings"] else 0
                    shared_mc_alt_gauss_integral_sigmas_rangesyst[index - 4*skipped_lowpt_bins] = mc_alt_plotinfo["gauss_settings"]["integral_rangesyst_std"] if "integral_rangesyst_std" in mc_alt_plotinfo["gauss_settings"] else 0

                    if mc_alt_plotinfo["psip_settings"]["bDraw"] and config_file["include_psip_result"]:   
                        mc_alt_integral_cb += mc_alt_plotinfo["psip_settings"]["integral"]
                        mc_alt_integral_cb_sigma = mc_alt_integral_cb_sigma |qs| mc_alt_plotinfo["psip_settings"]["Delta_integral"]

                    shared_mc_alt_integrals_cb[index - 4*skipped_lowpt_bins] = mc_alt_integral_cb
                    shared_mc_alt_integral_cb_sigmas[index - 4*skipped_lowpt_bins] = mc_alt_integral_cb_sigma
                    shared_mc_alt_integral_cb_sigmas_bksyst[index - 4*skipped_lowpt_bins] = mc_alt_plotinfo["cb_settings"]["bksyst_Delta_integral"] if "bksyst_Delta_integral" in mc_alt_plotinfo["cb_settings"] else 0
                    shared_mc_alt_integral_cb_sigmas_rangesyst[index - 4*skipped_lowpt_bins] = mc_alt_plotinfo["cb_settings"]["integral_rangesyst_std"] if "integral_rangesyst_std" in mc_alt_plotinfo["cb_settings"] else 0



def TnPEfficiencyInBin(isData, inputfile_name, workdir, plotdir, bin, binname, index, bDoFit, config_file, fix_sigma, fix_sigma_Delta):

    doBkg = isData or config_file['fit_background_on_MC']

    master_integral_range = [75,105] if config_file['resonance'] == "Z" else [2.7,3.5]

    inputfile = ROOT.TFile.Open(inputfile_name, "read")
    background_fit_bounds = config_file[('MC_' if not isData else '') + 'background_fit_bounds']
    peak_fit_bounds = config_file[('MC_' if not isData else '') + 'peak_fit_bounds']
    peak_fit_bounds_cb = config_file[('MC_' if not isData else '') + 'peak_fit_bounds_cb']
    background_fit_bounds_cmsshape = config_file[('MC_' if not isData else '') + 'background_fit_bounds_cmsshape']

    hist = inputfile.Get(bin)  

    if not ((index in config_file['skip_rebin_data'] and isData) or (index in config_file['skip_rebin_mc'] and not isData)):
        hist.Rebin()


    x_axis = []
    x_axis_pm = []
    x_axis_black = []
    y_axis_black = []
    y_axis = []
    y_axis_pm = []

    x_axis_unprocessed = []
    x_axis_unprocessed_pm = []
    y_axis_unprocessed = []
    y_axis_unprocessed_pm = []
    for i in range(hist.GetNbinsX()):
        binval = hist.GetBinCenter(i)

        if binval < 1.0 or binval > config_file['overflow']: # overflow bins
            continue

        bincontent = hist.GetBinContent(i)
        x_axis_unprocessed.append(binval)
        x_axis_unprocessed_pm.append(abs(binval - hist.GetBinLowEdge(i)))
        y_axis_unprocessed.append(bincontent)
        y_axis_unprocessed_pm.append(hist.GetBinError(i) if hist.GetBinError(i) > 0 else 0.00001)
        if bincontent > 0:
            x_axis.append(binval)
            x_axis_pm.append(abs(binval - hist.GetBinLowEdge(i)))
            y_axis.append(bincontent)
            y_axis_pm.append(hist.GetBinError(i))
        else:
            x_axis_black.append(binval)
            y_axis_black.append(bincontent)
        

    dualfit_x = []
    dualfit_x_pm = []
    dualfit_y = []
    dualfit_y_pm = []

    fitx = []
    fity = []
    chevx = []
    chevy = []
    chevy_pm = []
    peakx = []
    peaky = []
    peaky_pm = []
    peak_cb_x = []
    peak_cb_y = []
    peak_cb_y_pm = []

    psip_x = []
    psip_y = []
    psip_y_pm = []

    if isData or config_file['FilterZerosFromMCFitInput']:
        for i in range(len(x_axis)):
            x_above_psip =              (background_fit_bounds[index][4] < x_axis[i] < background_fit_bounds[index][5]) if len(background_fit_bounds[index]) > 5 else False
            x_above_unknown_resonance = (background_fit_bounds[index][6] < x_axis[i] < background_fit_bounds[index][7]) if len(background_fit_bounds[index]) > 7 else False
            if ((background_fit_bounds_cmsshape[index][0] < x_axis[i] < background_fit_bounds_cmsshape[index][1]) or (background_fit_bounds_cmsshape[index][2] < x_axis[i] < background_fit_bounds_cmsshape[index][3])):
                fitx.append(x_axis[i])
                fity.append(y_axis[i])

            if peak_fit_bounds[index][0] < x_axis[i] < peak_fit_bounds[index][1]:
                peakx.append(x_axis[i])
                peaky.append(y_axis[i])
                peaky_pm.append(y_axis_pm[i])

            if peak_fit_bounds_cb[index][0] < x_axis[i] < peak_fit_bounds_cb[index][1]:
                peak_cb_x.append(x_axis[i])
                peak_cb_y.append(y_axis[i])
                peak_cb_y_pm.append(y_axis_pm[i])

            if ((background_fit_bounds[index][0] < x_axis[i] < background_fit_bounds[index][1]) or (background_fit_bounds[index][2] < x_axis[i] < background_fit_bounds[index][3])) or x_above_psip or x_above_unknown_resonance:
                chevx.append(x_axis[i])
                chevy.append(y_axis[i])
                chevy_pm.append(y_axis_pm[i])

            if index in config_file['psip_fit_bounds']:
                if config_file['psip_fit_bounds'][index][0] < x_axis[i] < config_file['psip_fit_bounds'][index][1]:
                    psip_x.append(x_axis[i])
                    psip_y.append(y_axis[i])
                    psip_y_pm.append(y_axis_pm[i])

            if index in config_file[('MC_' if not isData else '') + 'dualfit_bounds_cb']:
                if config_file[('MC_' if not isData else '') + 'dualfit_bounds_cb'][index][0] < x_axis[i] < config_file[('MC_' if not isData else '') + 'dualfit_bounds_cb'][index][1]:
                    dualfit_x.append(x_axis[i])
                    dualfit_x_pm.append(x_axis_pm[i])
                    dualfit_y.append(y_axis[i])
                    dualfit_y_pm.append(y_axis_pm[i])

    else:
        for i in range(len(x_axis_unprocessed)):
            x_above_psip =              (background_fit_bounds[index][4] < x_axis_unprocessed[i] < background_fit_bounds[index][5]) if len(background_fit_bounds[index]) > 5 else False
            x_above_unknown_resonance = (background_fit_bounds[index][6] < x_axis_unprocessed[i] < background_fit_bounds[index][7]) if len(background_fit_bounds[index]) > 7 else False
            if ((background_fit_bounds_cmsshape[index][0] < x_axis_unprocessed[i] < background_fit_bounds_cmsshape[index][1]) or (background_fit_bounds_cmsshape[index][2] < x_axis_unprocessed[i] < background_fit_bounds_cmsshape[index][3])):
                fitx.append(x_axis_unprocessed[i])
                fity.append(y_axis_unprocessed[i])

            if peak_fit_bounds[index][0] < x_axis_unprocessed[i] < peak_fit_bounds[index][1]:
                peakx.append(x_axis_unprocessed[i])
                peaky.append(y_axis_unprocessed[i])
                peaky_pm.append(y_axis_unprocessed_pm[i])

            if peak_fit_bounds_cb[index][0] < x_axis_unprocessed[i] < peak_fit_bounds_cb[index][1]:
                peak_cb_x.append(x_axis_unprocessed[i])
                peak_cb_y.append(y_axis_unprocessed[i])
                peak_cb_y_pm.append(y_axis_unprocessed_pm[i])

            if ((background_fit_bounds[index][0] < x_axis_unprocessed[i] < background_fit_bounds[index][1]) or (background_fit_bounds[index][2] < x_axis_unprocessed[i] < background_fit_bounds[index][3])) or x_above_psip or x_above_unknown_resonance:
                chevx.append(x_axis_unprocessed[i])
                chevy.append(y_axis_unprocessed[i])
                chevy_pm.append(y_axis_unprocessed_pm[i])

            if index in config_file['psip_fit_bounds']:
                if config_file['psip_fit_bounds'][index][0] < x_axis_unprocessed[i] < config_file['psip_fit_bounds'][index][1]:
                    psip_x.append(x_axis_unprocessed[i])
                    psip_y.append(y_axis_unprocessed[i])
                    psip_y_pm.append(y_axis_unprocessed_pm[i])

            if index in config_file[('MC_' if not isData else '') + 'dualfit_bounds_cb']:
                if config_file[('MC_' if not isData else '') + 'dualfit_bounds_cb'][index][0] < x_axis[i] < config_file[('MC_' if not isData else '') + 'dualfit_bounds_cb'][index][1]:
                    dualfit_x.append(x_axis[i])
                    dualfit_x_pm.append(x_axis_pm[i])
                    dualfit_y.append(y_axis[i])
                    dualfit_y_pm.append(y_axis_pm[i])


    fig = plt.figure(99,figsize=(15, 10))
    ax = fig.add_subplot(111)
    ax.plot(x_axis,y_axis, "sb")
    ax.errorbar(x_axis, y_axis, yerr=y_axis_pm, xerr=x_axis_pm , capsize=3, fmt=".", ecolor = "black")
    ax.plot(x_axis_black, y_axis_black, "x", color="black")
    ax.set_xlabel("m [GeV]")
    ax.set_ylabel("Events")
    ax.xaxis.set_label_coords(.9, -.1)




    if isData and index in config_file['data_override_plots_xrange']:
        ax.set_xlim(config_file['data_override_plots_xrange'][index])
    else:
        ax.set_xlim(config_file['plots_xrange'])

    ax.set_title(binname)
    if not os.path.exists(workdir+plotdir):
        os.makedirs(workdir+plotdir)
    plt.savefig(workdir+plotdir+"figure"+str(index)+".png")
    plt.close()

    if not bDoFit:
        return False, False, False, False, False, False, False, False, False

    # another sanity check
    if len(x_axis) < 6:
        print ("Warning: Data " if isData else "Warning: MC ") + "Figure "+str(index)+" is likely empty, skip from fits"
        index += 1
        return False, False, False, False, False, False, False, False, False

    print ("Data " if isData else "MC ") + "Fit "+str(index)+", "+binname
    if doBkg:
        weights = []
        for sigma in chevy_pm:
            weights.append(1/sigma)

        chevparams, chevstats = Chev.chebfit(chevx, chevy, config_file['chevorder'], full=True, w=weights)

        A = Chev.chebvander(chevx, config_file['chevorder'])
        dof = len(chevx) - (config_file['chevorder'] + 1)  # degrees of freedom
        resid_var = chevstats[0][0] / dof
        cov = resid_var * np.linalg.inv(np.dot(A.T, A))
        chev_std_devs = np.sqrt(np.diag(cov))

        chevparams_gauss = list(chevparams)


        if config_file["bFitGauss"]:
            gauss_y_sub, gauss_y_pm_sub = subtract_background(peakx, peaky, peaky_pm, chevparams_gauss, True, chev_std_devs)

            try:
                gauss_popt_sig, gauss_perr_sig, gauss_integral, gauss_Delta_integral, gauss_integral_rangesyst_std = peak_fit("gauss" if not config_file['bUseBWGauss'] else "bwg", peakx, gauss_y_sub, gauss_y_pm_sub, config_file, isData, index)
            except Exception as er:
                print ("Data " if isData else "MC ") + "Gauss Fit "+str(index)+": "+str(er)
                gauss_popt_sig, gauss_perr_sig, gauss_integral, gauss_Delta_integral, gauss_integral_rangesyst_std = False, False, False, False, False

            if config_file['bUseBWGauss'] and index in config_file[('MC_' if not isData else '') + 'dualfit_bounds_bwg'] and gauss_integral:
                odrout = odr_dualfit_voigt(dualfit_x, dualfit_y, dualfit_x_pm, dualfit_y_pm, chevparams_gauss, gauss_popt_sig)
                gauss_popt_sig = odrout.beta[:4]
                gauss_perr_sig = odrout.sd_beta[:4]
                chevparams_gauss = odrout.beta[4:]  

                gauss_integral, gauss_Delta_integral, gauss_integral_rangesyst_std = getIntegralErrorWithMCMethod(Voigt, gauss_popt_sig, gauss_perr_sig, abs(peakx[1] - gauss_y_sub[0]), config_file, index, isData)


            gauss_settings = {
                "bDraw": bool(gauss_integral),
                "popt_sig": gauss_popt_sig,
                "perr_sig": gauss_perr_sig,
                "integral": gauss_integral,
                "Delta_integral": gauss_Delta_integral,
                #"bksyst_Delta_integral": integral_w_linbk - cb_integral,
                "integral_rangesyst_std": gauss_integral_rangesyst_std,
                "chevparams": chevparams_gauss,
            }
        else:
            gauss_popt_sig, gauss_perr_sig, gauss_integral, gauss_Delta_integral = False, False, False, False
            gauss_settings = {
                "bDraw": False,
            }


        cb_y_sub, cb_y_pm_sub = subtract_background(peak_cb_x, peak_cb_y, peak_cb_y_pm, chevparams, True, chev_std_devs)
        psip_y_sub, psip_y_pm_sub = subtract_background(psip_x, psip_y, psip_y_pm, chevparams, True, chev_std_devs)
        peak_cb_x_forsub = list(peak_cb_x)

        if index in config_file['slice_out_cb_zeros']:
            while(True):
                for i in range(len(peak_cb_x)):
                    if config_file['slice_out_cb_zeros'][index][0] < peak_cb_x[i] < config_file['slice_out_cb_zeros'][index][1] and cb_y_sub[i] <= 0:
                        cb_y_sub[i] = (cb_y_sub[i-1] + cb_y_sub[i+1])/2
                        print "averaged to "+str(cb_y_sub[i])+"at "+str(peak_cb_x[i])
                        break
                else: #wasnobreak
                    break

        if index in config_file['slice_out_points']:
            peak_slice = []
            peak_pm_slice = []
            peak_x_slice = []
            for i in range(len(peak_cb_x)):
                if not (config_file['slice_out_points'][index][0] <= peak_cb_x[i] <= config_file['slice_out_points'][index][1]):
                    peak_slice.append(cb_y_sub[i])
                    peak_pm_slice.append(cb_y_pm_sub[i])
                    peak_x_slice.append(peak_cb_x[i])
            
            peak_cb_x_forsub = peak_x_slice
            cb_y_sub = peak_slice
            cb_y_pm_sub = peak_pm_slice
        


        if config_file["bFitCrystalball"]:
            try:
                cb_popt_sig, cb_perr_sig, cb_integral, cb_Delta_integral, cb_integral_rangesyst_std = peak_fit("crystalball", peak_cb_x_forsub, cb_y_sub, cb_y_pm_sub, config_file, isData, index, fix_sigma, fix_sigma_Delta)    ### change: added fix_sigma for bad bin
            except Exception as er:
                print ("Data " if isData else "MC ") + "CB Fit "+str(index)+": "+str(er)
                cb_popt_sig, cb_perr_sig, cb_integral, cb_Delta_integral, cb_integral_rangesyst_std = False, False, False, False, False

            if index in config_file[('MC_' if not isData else '') + 'dualfit_bounds_cb'] and cb_integral:
                odrout = odr_dualfit_cb(dualfit_x, dualfit_y, dualfit_x_pm, dualfit_y_pm, chevparams, cb_popt_sig)
                cb_popt_sig = odrout.beta[:5]
                cb_perr_sig = odrout.sd_beta[:5]
                chevparams = odrout.beta[5:]
                #chev_std_devs = odrout.sd_beta[5:]     # FIXME: not currect, but...

                #test = []
                #for i in range(len(chevparams)):
                #    test.append(chev_std_devs[i]/abs(chevparams[i])*100)
                #print test

                cb_integral, cb_Delta_integral, cb_integral_rangesyst_std = getIntegralErrorWithMCMethod(mycrystalball, cb_popt_sig, cb_perr_sig, abs(peak_cb_x[1] - peak_cb_x[0]), config_file, index, isData)


            if not index in config_file["data_skip_bksyst_estimation"] and cb_integral:
                # estimate error from just linear background...
                # area = integral(sig+bk) - integral(linear bk)
                # first part:
                fullparams = list(cb_popt_sig)
                for i in range(len(chevparams)):
                    fullparams.append(chevparams[i])
                integral_w_sigbk, _, _ = getIntegralErrorWithMCMethod(crystalballPlusChev, fullparams, False, abs(peak_cb_x[1] - peak_cb_x[0]), config_file, index, isData, bEstimateError=False)

                # second part:
                y1 = Chev.chebval(master_integral_range[0],chevparams )
                y2 = Chev.chebval(master_integral_range[1],chevparams )
                rangelen = master_integral_range[1] - master_integral_range[0]
                # assume y1 < y2
                if y1 > y2:
                    temp = y1
                    y1 = y2
                    y2 = temp
                area_linbk = rangelen * y1 + (rangelen * (y2 - y1))/2

                integral_w_linbk = integral_w_sigbk - area_linbk / abs(peak_cb_x[1] - peak_cb_x[0])
            else:
                integral_w_linbk = cb_integral

            cb_settings = {
                "bDraw": bool(cb_integral),
                "popt_sig": cb_popt_sig,
                "perr_sig": cb_perr_sig,
                "integral": cb_integral,
                "Delta_integral": cb_Delta_integral,
                "bksyst_Delta_integral": integral_w_linbk - cb_integral,
                "integral_rangesyst_std": cb_integral_rangesyst_std,
                "fix_sigma": fix_sigma  
            }

            psip_popt_sig, psip_perr_sig, psip_integral, psip_Delta_integral = False, False, False, False
            if index in config_file['psip_fit_bounds'] and cb_integral:
                try:
                    psip_popt_sig, psip_perr_sig, psip_integral, psip_Delta_integral = fit_psip(index, psip_x, psip_y_sub, psip_y_pm_sub, config_file, cb_popt_sig[4], cb_perr_sig[4]) 
                except Exception as er:
                    print ("Data " if isData else "MC ") + "Psi' Fit "+str(index)+": "+str(er)
        
                psip_settings = {
                    "bDraw": bool(psip_integral),
                    "popt_sig": psip_popt_sig,
                    "perr_sig": psip_perr_sig,
                    "integral": psip_integral,
                    "Delta_integral": psip_Delta_integral,
                }
            else:
                psip_settings = {
                    "bDraw": False,
                }
        else:
            cb_popt_sig, cb_perr_sig, cb_integral, cb_Delta_integral = False, False, False, False
            cb_settings = {
                "bDraw": False,
            }
            psip_settings = {
                "bDraw": False,
            }


        bksub_integral, bksub_integral_pm, linbksub_integral, rangesyst_std = get_bksub_integral(x_axis, y_axis, y_axis_pm, chevparams, chev_std_devs, master_integral_range, config_file["integral_error_MCmethod_iterations"])

        if index in config_file["data_skip_bksyst_estimation"]:
            linbksub_integral = bksub_integral

        makePlot(workdir+plotdir, config_file, binname, index, isData, x_axis, x_axis_pm, y_axis, y_axis_pm, x_axis_black, y_axis_black, doBkg, chevx, chevparams, chevstats, bksub_integral, bksub_integral_pm, linbksub_integral, rangesyst_std, gauss_settings, cb_settings, psip_settings)
    else:
        if config_file["bFitGauss"]:
            try:
                gauss_popt_sig, gauss_perr_sig, gauss_integral, gauss_Delta_integral, gauss_integral_rangesyst_std = peak_fit("gauss" if not config_file['bUseBWGauss'] else "bwg", peakx, peaky, peaky_pm, config_file, isData, index)
            except Exception as er:
                print ("Data " if isData else "MC ") + "Fit "+str(index)+": "+str(er)
                gauss_popt_sig, gauss_perr_sig, gauss_integral, gauss_Delta_integral, gauss_integral_rangesyst_std = False, False, False, False

            gauss_settings = {
                "bDraw": bool(gauss_integral),
                "popt_sig": gauss_popt_sig,
                "perr_sig": gauss_perr_sig,
                "integral": gauss_integral,
                "Delta_integral": gauss_Delta_integral,
                "integral_rangesyst_std": gauss_integral_rangesyst_std,
            }
        else:
            gauss_popt_sig, gauss_perr_sig, gauss_integral, gauss_Delta_integral = False, False, False, False
            gauss_settings = {
                "bDraw": False,
            }


        if config_file["bFitCrystalball"]:

            try:
                cb_popt_sig, cb_perr_sig, cb_integral, cb_Delta_integral, cb_integral_rangesyst_std = peak_fit("crystalball", peak_cb_x, peak_cb_y, peak_cb_y_pm, config_file, isData, index, fix_sigma, fix_sigma_Delta)
            except Exception as er:
                print ("Data " if isData else "MC ") + "Fit "+str(index)+": "+str(er)
                cb_popt_sig, cb_perr_sig, cb_integral, cb_Delta_integral, cb_integral_rangesyst_std = False, False, False, False, False

            cb_settings = {
                "bDraw": bool(cb_integral),
                "popt_sig": cb_popt_sig,
                "perr_sig": cb_perr_sig,
                "integral": cb_integral,
                "Delta_integral": cb_Delta_integral,
                "integral_rangesyst_std": cb_integral_rangesyst_std,
                "fix_sigma": fix_sigma
            }

            psip_popt_sig, psip_perr_sig, psip_integral, psip_Delta_integral = False, False, False, False
            if index in config_file['psip_fit_bounds'] and cb_integral and not fix_sigma:
                try:
                    psip_popt_sig, psip_perr_sig, psip_integral, psip_Delta_integral = fit_psip(index, psip_x, psip_y, psip_y_pm, config_file, cb_popt_sig[4], cb_perr_sig[4])

                    cb_integral += psip_integral
                    cb_Delta_integral = cb_Delta_integral + psip_Delta_integral
                except Exception as er:
                    print ("Data " if isData else "MC ") + "Psi' Fit "+str(index)+": "+str(er)

                psip_settings = {
                    "bDraw": bool(psip_integral),
                    "popt_sig": psip_popt_sig,
                    "perr_sig": psip_perr_sig,
                    "integral": psip_integral,
                    "Delta_integral": psip_Delta_integral,
                }
            else:
                psip_settings = {
                    "bDraw": False,
                }
        else:
            cb_popt_sig, cb_perr_sig, cb_integral, cb_Delta_integral = False, False, False, False
            cb_settings = {
                "bDraw": False,
            }
            psip_settings = {
                "bDraw": False,
            }


        # find out the range
        histxaxis = hist.GetXaxis()
        lowbinSystDown = -1
        lowbin = -1
        lowbinSystUp = -1
        highbin = -1
        highbinSystDown = -1
        highbinSystUp = -1
        for k in range(1,histxaxis.GetNbins() + 1):
            if lowbinSystDown == -1 and histxaxis.GetBinCenter(k) >= master_integral_range[0]*0.95:
                lowbinSystDown = k
            if lowbin == -1 and histxaxis.GetBinCenter(k) >= master_integral_range[0]:
                lowbin = k
            if lowbinSystUp == -1 and histxaxis.GetBinCenter(k) >= master_integral_range[0]*1.05:
                lowbinSystUp = k

            if highbinSystDown == -1 and histxaxis.GetBinCenter(k) > master_integral_range[1]*0.95:
                highbinSystDown = k - 1
            if highbin == -1 and histxaxis.GetBinCenter(k) > master_integral_range[1]:
                highbin = k - 1
            if highbinSystUp == -1 and histxaxis.GetBinCenter(k) > master_integral_range[1]*1.05:
                highbinSystUp = k - 1
                break

        hist_integral_error = ctypes.c_double(0)
        hist_integral = hist.IntegralAndError(lowbin, highbin, hist_integral_error)
        hist_integral_error = hist_integral_error.value


        hist_integral_low_error = ctypes.c_double(0)
        hist_integral_low  = hist.IntegralAndError(lowbinSystUp,   highbinSystDown, hist_integral_low_error)
        hist_integral_high = hist.IntegralAndError(lowbinSystDown, highbinSystUp,   hist_integral_low_error)

        hist_integral_rangesyst_std = useBigger(abs(hist_integral - hist_integral_high), abs(hist_integral - hist_integral_low))

        makePlot(workdir+plotdir, config_file, binname, index, isData, x_axis, x_axis_pm, y_axis, y_axis_pm, x_axis_black, y_axis_black, doBkg, False, False, False, hist_integral, hist_integral_error, 0, hist_integral_rangesyst_std, gauss_settings, cb_settings, psip_settings)

    if doBkg:
        raw_integral = bksub_integral
        raw_integral_sigma = bksub_integral_pm
        linkb_hist_integral_Delta = linbksub_integral - raw_integral
        raw_integral_rangesyst_std = rangesyst_std
    else: 
        raw_integral = hist_integral
        raw_integral_sigma = hist_integral_error
        linkb_hist_integral_Delta = False
        raw_integral_rangesyst_std = hist_integral_rangesyst_std


    if cb_integral:
        return_sigma_for_fixed_width = cb_popt_sig[4]
        return_Delta_sigma_for_fixed_width = cb_perr_sig[4]
    else:
        return_sigma_for_fixed_width = False
        return_Delta_sigma_for_fixed_width = False


    plotinfo = {
        "x_axis": x_axis,
        "x_axis_pm": x_axis_pm,
        "y_axis": y_axis,
        "y_axis_pm": y_axis_pm,
        "x_axis_black": x_axis_black,
        "y_axis_black": y_axis_black,
        "chevx": chevx if doBkg else False,
        "chevparams": chevparams if doBkg else False,
        "chevstats": chevstats if doBkg else False,
        "chev_std_devs": chev_std_devs if doBkg else False,
        "hist_integral": raw_integral,
        "Delta_hist_integral": raw_integral_sigma,
        "linkb_hist_integral_Delta":  linkb_hist_integral_Delta,
        "hist_integral_rangesyst_std": raw_integral_rangesyst_std,
        "gauss_settings": gauss_settings,
        "cb_settings": cb_settings,
        "psip_settings": psip_settings,
    }

    inputfile.Close()
    return plotinfo, raw_integral, raw_integral_sigma, gauss_integral, gauss_Delta_integral, cb_integral, cb_Delta_integral, return_sigma_for_fixed_width, return_Delta_sigma_for_fixed_width



def get_efficiency_error(P, F, DP, DF):
    return math.sqrt(  (  (  F/((P+F)**2)   )*DP   )**2    +    (  (  -P/((P+F)**2)     )*DF   )**2   )


def write_mparray(name, mparray, outfile):
    outfile.write(name+": [")
    for i in range(len(mparray) - 1):
        outfile.write(str(mparray[i])+", ")
    outfile.write(str(mparray[-1]) + "]\n")


def main():

    options = get_parser().parse_args()

    with open(options.config_file, 'r') as file:
        config_file = yaml.safe_load(file)

    for i in range(len(options.config_file)-1, 0, -1):
        if options.config_file[i] == '/':
            workdir = options.config_file[:i]+"/"
            break


    #for i in range(len(config_file['raw_integral_bounds'])):
    #    if i in config_file['MC_dualfit_bounds_cb']:
    #        config_file['MC_raw_integral_bounds'][i] = config_file['MC_dualfit_bounds_cb'][i]
    #    if i in config_file['dualfit_bounds_cb']:
    #        config_file['raw_integral_bounds'][i] = config_file['dualfit_bounds_cb'][i]


    datafile = workdir + config_file['datafile_name']
    mcfile = workdir + config_file['mcfile_name']
    if config_file['mcalt_name']:
        mcfile_alt = workdir + config_file['mcalt_name']
    else:
        mcfile_alt = False

    ptbins = config_file['ptbins']
    ptbins_num = config_file['ptbins_num']
    ptbins_center = []
    for i in range(1,len(ptbins_num),1):
        ptbins_center.append((ptbins_num[i] + ptbins_num[i-1])/2.0 )

    binw = []
    for i in range(len(ptbins_center)):
        binw.append(ptbins_center[i] - ptbins_num[i])

    bins = []
    binnames = []

    binindex_range_to_fit = []

    for i in range(len(ptbins) - 1):
        if len(binindex_range_to_fit) == 0 and ptbins_num[i] >= config_file['fit_in_pt_range'][0]:
            binindex_range_to_fit.append(len(bins))

        if len(binindex_range_to_fit) == 1 and ptbins_num[i] >= config_file['fit_in_pt_range'][1]:
            binindex_range_to_fit.append(len(bins))

        exec('bins.append("bin%s_abs(el_sc_eta)_0p00To1p44_el_pt_%sTo%s_Pass") ' % (str(3*i).zfill(2), ptbins[i], ptbins[i+1]))
        exec('bins.append("bin%s_abs(el_sc_eta)_0p00To1p44_el_pt_%sTo%s_Fail") ' % (str(3*i).zfill(2), ptbins[i], ptbins[i+1]))
        exec('bins.append("bin%s_abs(el_sc_eta)_1p57To2p50_el_pt_%sTo%s_Pass") ' % (str(3*i + 2).zfill(2), ptbins[i], ptbins[i+1]))
        exec('bins.append("bin%s_abs(el_sc_eta)_1p57To2p50_el_pt_%sTo%s_Fail") ' % (str(3*i + 2).zfill(2), ptbins[i], ptbins[i+1]))

        exec('binnames.append(r"Bin%s, $0 < \| \eta \| < 1.4442$,  $%s < p_{T} < %s$, PASS") ' % (str(2*i), ptbins[i][:-3], ptbins[i+1][:-3]))
        exec('binnames.append(r"Bin%s, $0 < \| \eta \| < 1.4442$,  $%s < p_{T} < %s$, FAIL") ' % (str(2*i), ptbins[i][:-3], ptbins[i+1][:-3]))
        exec('binnames.append(r"Bin%s, $1.566 < \| \eta \| < 2.5$,  $%s < p_{T} < %s$, PASS") ' % (str(2*i + 1), ptbins[i][:-3], ptbins[i+1][:-3]))
        exec('binnames.append(r"Bin%s, $1.566 < \| \eta \| < 2.5$,  $%s < p_{T} < %s$, FAIL") ' % (str(2*i + 1), ptbins[i][:-3], ptbins[i+1][:-3]))

    if len(binindex_range_to_fit) == 1 :
        binindex_range_to_fit.append(len(bins))


    ptbins_center_fit = ptbins_center[(binindex_range_to_fit[0]/4):(binindex_range_to_fit[1]/4)]
    binw_fit = binw[(binindex_range_to_fit[0]/4):(binindex_range_to_fit[1]/4)]

    if options.threads is None:
        threads = multiprocessing.cpu_count()

    elif options.threads == -1:
        threads = multiprocessing.cpu_count()
    elif options.threads < 1:
        print "Ain't the sharpest tool, are we?"
        quit()

    else:
        threads = options.threads


    # get number of fitted bins
    skipped_lowpt_bins = 0
    num_fitted_pt_bins = 0
    for pt in config_file['ptbins_num']:
        if pt < config_file['fit_in_pt_range'][0]:
            skipped_lowpt_bins += 1
        elif pt < config_file['fit_in_pt_range'][1]:
            num_fitted_pt_bins += 1
        else:
            break


    # decide how many bins per thread
    if num_fitted_pt_bins*4 <= threads:
        threads = num_fitted_pt_bins*4
        files_per_job = 1
        extrafiles = 0
    else:
        files_per_job = int(num_fitted_pt_bins*4 / threads)
        extrafiles = num_fitted_pt_bins*4 % threads


    indexlists = []
    binindex = skipped_lowpt_bins * 4
    for i in range(threads):
        indexlists.append([])
        for j in range(files_per_job):
            indexlists[-1].append(binindex)
            binindex += 1

        if extrafiles > 0:
            indexlists[-1].append(binindex)
            binindex += 1
            extrafiles -= 1
    
    def cyclic_index(index, threads):
        while(index >= threads):
            index -= threads
        return index
        
    for i in range(skipped_lowpt_bins):
        for j in range(4):
            indexlists[cyclic_index(j, threads)].append(4*i + j)

    free_worker = 0
    for i in range(binindex, len(bins), 1):
        indexlists[cyclic_index(free_worker, threads)].append(i)
        free_worker += 1


    mutex = multiprocessing.Lock()
    shared_data_raw_integrals = multiprocessing.Array('d',range(num_fitted_pt_bins * 4))
    shared_data_raw_integral_sigmas = multiprocessing.Array('d',range(num_fitted_pt_bins * 4))
    shared_data_raw_integral_sigmas_rangesyst = multiprocessing.Array('d',range(num_fitted_pt_bins * 4))
    shared_data_raw_integral_sigmas_bksyst = multiprocessing.Array('d',range(num_fitted_pt_bins * 4))
    shared_data_gauss_integrals = multiprocessing.Array('d',range(num_fitted_pt_bins * 4))
    shared_data_gauss_integral_sigmas = multiprocessing.Array('d',range(num_fitted_pt_bins * 4))
    shared_data_gauss_integral_sigmas_rangesyst = multiprocessing.Array('d',range(num_fitted_pt_bins * 4))
    shared_data_gauss_integral_sigmas_bksyst = multiprocessing.Array('d',range(num_fitted_pt_bins * 4))
    shared_data_integrals_cb = multiprocessing.Array('d',range(num_fitted_pt_bins * 4))
    shared_data_integral_cb_sigmas = multiprocessing.Array('d',range(num_fitted_pt_bins * 4))
    shared_data_integral_cb_sigmas_rangesyst = multiprocessing.Array('d',range(num_fitted_pt_bins * 4))
    shared_data_integral_cb_sigmas_bksyst = multiprocessing.Array('d',range(num_fitted_pt_bins * 4))
    shared_mc_raw_integrals = multiprocessing.Array('d',range(num_fitted_pt_bins * 4))
    shared_mc_raw_integral_sigmas = multiprocessing.Array('d',range(num_fitted_pt_bins * 4))
    shared_mc_raw_integral_sigmas_rangesyst = multiprocessing.Array('d',range(num_fitted_pt_bins * 4))
    shared_mc_raw_integral_sigmas_bksyst = multiprocessing.Array('d',range(num_fitted_pt_bins * 4))
    shared_mc_gauss_integrals = multiprocessing.Array('d',range(num_fitted_pt_bins * 4))
    shared_mc_gauss_integral_sigmas = multiprocessing.Array('d',range(num_fitted_pt_bins * 4))
    shared_mc_gauss_integral_sigmas_rangesyst = multiprocessing.Array('d',range(num_fitted_pt_bins * 4))
    shared_mc_gauss_integral_sigmas_bksyst = multiprocessing.Array('d',range(num_fitted_pt_bins * 4))
    shared_mc_integrals_cb = multiprocessing.Array('d',range(num_fitted_pt_bins * 4))
    shared_mc_integral_cb_sigmas = multiprocessing.Array('d',range(num_fitted_pt_bins * 4))
    shared_mc_integral_cb_sigmas_rangesyst = multiprocessing.Array('d',range(num_fitted_pt_bins * 4))
    shared_mc_integral_cb_sigmas_bksyst = multiprocessing.Array('d',range(num_fitted_pt_bins * 4))
    shared_mc_alt_raw_integrals = multiprocessing.Array('d',range(num_fitted_pt_bins * 4))
    shared_mc_alt_raw_integral_sigmas = multiprocessing.Array('d',range(num_fitted_pt_bins * 4))
    shared_mc_alt_raw_integral_sigmas_rangesyst = multiprocessing.Array('d',range(num_fitted_pt_bins * 4))
    shared_mc_alt_raw_integral_sigmas_bksyst = multiprocessing.Array('d',range(num_fitted_pt_bins * 4))
    shared_mc_alt_gauss_integrals = multiprocessing.Array('d',range(num_fitted_pt_bins * 4))
    shared_mc_alt_gauss_integral_sigmas = multiprocessing.Array('d',range(num_fitted_pt_bins * 4))
    shared_mc_alt_gauss_integral_sigmas_rangesyst = multiprocessing.Array('d',range(num_fitted_pt_bins * 4))
    shared_mc_alt_gauss_integral_sigmas_bksyst = multiprocessing.Array('d',range(num_fitted_pt_bins * 4))
    shared_mc_alt_integrals_cb = multiprocessing.Array('d',range(num_fitted_pt_bins * 4))
    shared_mc_alt_integral_cb_sigmas = multiprocessing.Array('d',range(num_fitted_pt_bins * 4))
    shared_mc_alt_integral_cb_sigmas_rangesyst = multiprocessing.Array('d',range(num_fitted_pt_bins * 4))
    shared_mc_alt_integral_cb_sigmas_bksyst = multiprocessing.Array('d',range(num_fitted_pt_bins * 4))



    print "Opened all input resources"

    proc = {}
    for i in range(threads):
        proc[i] = multiprocessing.Process(target=TnPEfficiency, args=(mutex, shared_data_raw_integrals, shared_data_raw_integral_sigmas, shared_data_raw_integral_sigmas_rangesyst, shared_data_raw_integral_sigmas_bksyst, shared_data_gauss_integrals, shared_data_gauss_integral_sigmas, shared_data_gauss_integral_sigmas_rangesyst, shared_data_gauss_integral_sigmas_bksyst, shared_data_integrals_cb, shared_data_integral_cb_sigmas, shared_data_integral_cb_sigmas_rangesyst, shared_data_integral_cb_sigmas_bksyst, shared_mc_raw_integrals, shared_mc_raw_integral_sigmas, shared_mc_raw_integral_sigmas_rangesyst, shared_mc_raw_integral_sigmas_bksyst, shared_mc_gauss_integrals, shared_mc_gauss_integral_sigmas, shared_mc_gauss_integral_sigmas_rangesyst, shared_mc_gauss_integral_sigmas_bksyst, shared_mc_integrals_cb, shared_mc_integral_cb_sigmas, shared_mc_integral_cb_sigmas_rangesyst, shared_mc_integral_cb_sigmas_bksyst, shared_mc_alt_raw_integrals, shared_mc_alt_raw_integral_sigmas, shared_mc_alt_raw_integral_sigmas_rangesyst, shared_mc_alt_raw_integral_sigmas_bksyst, shared_mc_alt_gauss_integrals, shared_mc_alt_gauss_integral_sigmas, shared_mc_alt_gauss_integral_sigmas_rangesyst, shared_mc_alt_gauss_integral_sigmas_bksyst, shared_mc_alt_integrals_cb, shared_mc_alt_integral_cb_sigmas, shared_mc_alt_integral_cb_sigmas_rangesyst, shared_mc_alt_integral_cb_sigmas_bksyst, 
                       datafile, mcfile, mcfile_alt, bins, binnames, indexlists[i], binindex_range_to_fit, workdir, config_file, skipped_lowpt_bins))
        proc[i].start()

    for i in range(threads):
        proc[i].join()
    

    print "Making efficiency plots"


    ## FIXME FIXME FIXME
    # A BWG-re kulon bksyst kell ha valaha nezzuk meg!! Most csak syntax miatt bent van a CB syst

    if config_file['resonance'] == "Z":
        makeEfficiencyPlot(workdir+"efficiency-raw", shared_data_raw_integrals, shared_data_raw_integral_sigmas, shared_data_raw_integral_sigmas_bksyst, shared_data_raw_integral_sigmas_rangesyst, shared_mc_raw_integrals, shared_mc_raw_integral_sigmas, shared_mc_raw_integral_sigmas_bksyst, shared_mc_raw_integral_sigmas_rangesyst, ptbins_center_fit, binw_fit)
        makeEfficiencyPlot(workdir+"efficiency-bwg", shared_data_gauss_integrals, shared_data_gauss_integral_sigmas, shared_data_gauss_integral_sigmas_bksyst, shared_data_gauss_integral_sigmas_rangesyst, shared_mc_gauss_integrals, shared_mc_gauss_integral_sigmas, shared_mc_gauss_integral_sigmas_bksyst, shared_mc_gauss_integral_sigmas_rangesyst, ptbins_center_fit, binw_fit)
        makeEfficiencyPlot(workdir+"efficiency-crystalball", shared_data_integrals_cb, shared_data_integral_cb_sigmas, shared_data_integral_cb_sigmas_bksyst, shared_data_integral_cb_sigmas_rangesyst, shared_mc_integrals_cb, shared_mc_integral_cb_sigmas, shared_mc_integral_cb_sigmas_bksyst, shared_mc_integral_cb_sigmas_rangesyst, ptbins_center_fit, binw_fit)

        if config_file['mcalt_name']:
            makeEfficiencyPlot(workdir+"efficiency-raw-alt", shared_data_raw_integrals, shared_data_raw_integral_sigmas, shared_data_raw_integral_sigmas_bksyst, shared_data_raw_integral_sigmas_rangesyst, shared_mc_alt_raw_integrals, shared_mc_alt_raw_integral_sigmas, shared_mc_alt_raw_integral_sigmas_bksyst, shared_mc_alt_raw_integral_sigmas_rangesyst, ptbins_center_fit, binw_fit)
            makeEfficiencyPlot(workdir+"efficiency-bwg-alt", shared_data_gauss_integrals, shared_data_gauss_integral_sigmas, shared_data_gauss_integral_sigmas_bksyst, shared_data_gauss_integral_sigmas_rangesyst, shared_mc_alt_gauss_integrals, shared_mc_alt_gauss_integral_sigmas, shared_mc_alt_gauss_integral_sigmas_bksyst, shared_mc_alt_gauss_integral_sigmas_rangesyst, ptbins_center_fit, binw_fit)
            makeEfficiencyPlot(workdir+"efficiency-crystalball-alt", shared_data_integrals_cb, shared_data_integral_cb_sigmas, shared_data_integral_cb_sigmas_bksyst, shared_data_integral_cb_sigmas_rangesyst, shared_mc_alt_integrals_cb, shared_mc_alt_integral_cb_sigmas, shared_mc_alt_integral_cb_sigmas_bksyst, shared_mc_alt_integral_cb_sigmas_rangesyst, ptbins_center_fit, binw_fit)

    # For JPsi
    else:
        makeEfficiencyPlot(workdir+"efficiency-raw", shared_data_raw_integrals, shared_data_raw_integral_sigmas, shared_data_raw_integral_sigmas_bksyst, shared_data_raw_integral_sigmas_rangesyst, shared_mc_raw_integrals, shared_mc_raw_integral_sigmas, shared_mc_raw_integral_sigmas_bksyst, shared_mc_raw_integral_sigmas_rangesyst, ptbins_center_fit, binw_fit)
        makeEfficiencyPlot(workdir+"efficiency-crystalball-raw", shared_data_integrals_cb, shared_data_integral_cb_sigmas, shared_data_integral_cb_sigmas_bksyst, shared_data_integral_cb_sigmas_rangesyst, shared_mc_raw_integrals, shared_mc_raw_integral_sigmas, shared_mc_raw_integral_sigmas_bksyst, shared_mc_raw_integral_sigmas_rangesyst, ptbins_center_fit, binw_fit)
        makeEfficiencyPlot(workdir+"efficiency-crystalball-fixedwidth", shared_data_integrals_cb, shared_data_integral_cb_sigmas, shared_data_integral_cb_sigmas_bksyst, shared_data_integral_cb_sigmas_rangesyst, shared_mc_integrals_cb, shared_mc_integral_cb_sigmas, shared_mc_integral_cb_sigmas_bksyst, shared_mc_integral_cb_sigmas_rangesyst, ptbins_center_fit, binw_fit)


    with open(workdir+"/results_"+config_file["resonance"]+"_chev"+str(config_file['chevorder'])+".yaml", "w") as out:
        write_mparray("shared_data_raw_integrals", shared_data_raw_integrals, out)
        write_mparray("shared_data_raw_integral_sigmas", shared_data_raw_integral_sigmas, out)
        write_mparray("shared_data_raw_integral_sigmas_bksyst", shared_data_raw_integral_sigmas_bksyst, out)
        write_mparray("shared_data_raw_integral_sigmas_rangesyst", shared_data_raw_integral_sigmas_rangesyst, out)
        write_mparray("shared_data_gauss_integrals", shared_data_gauss_integrals, out)
        write_mparray("shared_data_gauss_integral_sigmas", shared_data_gauss_integral_sigmas, out)
        write_mparray("shared_data_gauss_integral_sigmas_bksyst", shared_data_gauss_integral_sigmas_bksyst, out)
        write_mparray("shared_data_gauss_integral_sigmas_rangesyst", shared_data_gauss_integral_sigmas_rangesyst, out)
        write_mparray("shared_data_integrals_cb", shared_data_integrals_cb, out)
        write_mparray("shared_data_integral_cb_sigmas", shared_data_integral_cb_sigmas, out)
        write_mparray("shared_data_integral_cb_sigmas_bksyst", shared_data_integral_cb_sigmas_bksyst, out)
        write_mparray("shared_data_integral_cb_sigmas_rangesyst", shared_data_integral_cb_sigmas_rangesyst, out)
        write_mparray("shared_mc_raw_integrals", shared_mc_raw_integrals, out)
        write_mparray("shared_mc_raw_integral_sigmas",shared_mc_raw_integral_sigmas , out)
        write_mparray("shared_mc_raw_integral_sigmas_bksyst",shared_mc_raw_integral_sigmas_bksyst , out)
        write_mparray("shared_mc_raw_integral_sigmas_rangesyst",shared_mc_raw_integral_sigmas_rangesyst , out)
        write_mparray("shared_mc_gauss_integrals", shared_mc_gauss_integrals, out)
        write_mparray("shared_mc_gauss_integral_sigmas",shared_mc_gauss_integral_sigmas , out)
        write_mparray("shared_mc_gauss_integral_sigmas_bksyst",shared_mc_gauss_integral_sigmas_bksyst , out)
        write_mparray("shared_mc_gauss_integral_sigmas_rangesyst",shared_mc_gauss_integral_sigmas_rangesyst , out)
        write_mparray("shared_mc_integrals_cb", shared_mc_integrals_cb, out)
        write_mparray("shared_mc_integral_cb_sigmas", shared_mc_integral_cb_sigmas, out)
        write_mparray("shared_mc_integral_cb_sigmas_bksyst", shared_mc_integral_cb_sigmas_bksyst, out)
        write_mparray("shared_mc_integral_cb_sigmas_rangesyst", shared_mc_integral_cb_sigmas_rangesyst, out)
        write_mparray("shared_mc_alt_raw_integrals", shared_mc_alt_raw_integrals, out)
        write_mparray("shared_mc_alt_raw_integral_sigmas", shared_mc_alt_raw_integral_sigmas, out)
        write_mparray("shared_mc_alt_raw_integral_sigmas_bksyst", shared_mc_alt_raw_integral_sigmas_bksyst, out)
        write_mparray("shared_mc_alt_raw_integral_sigmas_rangesyst", shared_mc_alt_raw_integral_sigmas_rangesyst, out)
        write_mparray("shared_mc_alt_gauss_integrals", shared_mc_alt_gauss_integrals, out)
        write_mparray("shared_mc_alt_gauss_integral_sigmas", shared_mc_alt_gauss_integral_sigmas, out)
        write_mparray("shared_mc_alt_gauss_integral_sigmas_bksyst", shared_mc_alt_gauss_integral_sigmas_bksyst, out)
        write_mparray("shared_mc_alt_gauss_integral_sigmas_rangesyst", shared_mc_alt_gauss_integral_sigmas_rangesyst, out)
        write_mparray("shared_mc_alt_integrals_cb", shared_mc_alt_integrals_cb, out)
        write_mparray("shared_mc_alt_integral_cb_sigmas", shared_mc_alt_integral_cb_sigmas, out)
        write_mparray("shared_mc_alt_integral_cb_sigmas_bksyst", shared_mc_alt_integral_cb_sigmas_bksyst, out)
        write_mparray("shared_mc_alt_integral_cb_sigmas_rangesyst", shared_mc_alt_integral_cb_sigmas_rangesyst, out)

        write_mparray("ptbins_center_fit", ptbins_center_fit, out)
        write_mparray("binw_fit", binw_fit, out)



def makeEfficiencyPlot(filename, data_integrals, data_integral_sigmas, data_integral_sigmas_bksyst, data_integral_sigmas_rangesyst, mc_integrals, mc_integral_sigmas, mc_integral_sigmas_bksyst, mc_integral_sigmas_rangesyst, ptbins_center_fit, binw_fit):
    data_eff_barrel = []
    data_eff_barrel_error = []
    data_eff_barrel_error_total = []
    mc_eff_barrel = []
    mc_eff_barrel_error = []
    mc_eff_barrel_error_total = []
    data_eff_endcap = []
    data_eff_endcap_error = []
    data_eff_endcap_error_total = []
    mc_eff_endcap = []
    mc_eff_endcap_error = []
    mc_eff_endcap_error_total = []

    scale_factor_barrel = []
    scale_factor_barrel_error = []
    scale_factor_x_barrel = []
    scale_factor_x_barrel_binw = []
    scale_factor_barrel_black = []
    scale_factor_barrel_black_x = []
    scale_factor_barrel_black_x_binw = []

    scale_factor_endcap = []
    scale_factor_endcap_error = []
    scale_factor_x_endcap = []
    scale_factor_x_endcap_binw = []
    scale_factor_endcap_black = []
    scale_factor_endcap_black_x = []
    scale_factor_endcap_black_x_binw = []

    pt_bin_number = 0
    for i in range(0,len(data_integrals),4):
        
        if data_integrals[i] and data_integrals[i+1]:
            data_barrel_pass = data_integrals[i]
            data_barrel_pass_error = data_integral_sigmas[i] if data_integral_sigmas[i] and data_integral_sigmas[i]/data_integrals[i] < 1 else False
            data_barrel_pass_error_bksyst = abs(data_integral_sigmas_bksyst[i])
            data_barrel_pass_error_rangesyst = abs(data_integral_sigmas_rangesyst[i])

            data_barrel_fail = data_integrals[i+1] 
            data_barrel_fail_error = data_integral_sigmas[i+1] if data_integral_sigmas[i+1] and data_integral_sigmas[i+1]/data_integrals[i+1] < 1 else False
            data_barrel_fail_error_bksyst = abs(data_integral_sigmas_bksyst[i+1])
            data_barrel_fail_error_rangesyst = abs(data_integral_sigmas_rangesyst[i+1])

            data_eff_barrel.append(data_barrel_pass / (data_barrel_pass + data_barrel_fail))
            data_eff_barrel_error.append( get_efficiency_error(data_barrel_pass, data_barrel_fail, data_barrel_pass_error, data_barrel_fail_error) if data_barrel_pass_error and data_barrel_fail_error else 0  )
            data_eff_barrel_error_total.append( get_efficiency_error(data_barrel_pass, data_barrel_fail, qsum(data_barrel_pass_error, data_barrel_pass_error_bksyst, data_barrel_pass_error_rangesyst), qsum(data_barrel_fail_error, data_barrel_fail_error_bksyst,data_barrel_fail_error_rangesyst)) if data_barrel_pass_error and data_barrel_fail_error else 0  )

        else:
            data_eff_barrel.append(-1)
            data_eff_barrel_error.append( 0.001 )
            data_eff_barrel_error_total.append( 0.001 )


        if data_integrals[i+2] and data_integrals[i+3]:
            data_endcap_pass = data_integrals[i+2]
            data_endcap_pass_error = data_integral_sigmas[i+2] if data_integral_sigmas[i+2]  and data_integral_sigmas[i+2]/data_integrals[i+2] < 1 else False
            data_endcap_pass_error_bksyst = abs(data_integral_sigmas_bksyst[i+2])
            data_endcap_pass_error_rangesyst = abs(data_integral_sigmas_rangesyst[i+2])

            data_endcap_fail = data_integrals[i+3]
            data_endcap_fail_error = data_integral_sigmas[i+3] if data_integral_sigmas[i+3]  and data_integral_sigmas[i+3]/data_integrals[i+3] < 1 else False
            data_endcap_fail_error_bksyst = abs(data_integral_sigmas_bksyst[i+3])
            data_endcap_fail_error_rangesyst = abs(data_integral_sigmas_rangesyst[i+3])

            data_eff_endcap.append(data_endcap_pass / (data_endcap_pass + data_endcap_fail))
            data_eff_endcap_error.append( get_efficiency_error(data_endcap_pass, data_endcap_fail, data_endcap_pass_error, data_endcap_fail_error) if data_endcap_pass_error and data_endcap_fail_error else 0  )
            data_eff_endcap_error_total.append( get_efficiency_error(data_endcap_pass, data_endcap_fail, qsum(data_endcap_pass_error, data_endcap_pass_error_bksyst, data_endcap_pass_error_rangesyst), qsum(data_endcap_fail_error, data_endcap_fail_error_bksyst,data_endcap_fail_error_rangesyst)) if data_endcap_pass_error and data_endcap_fail_error else 0  )
        else:
            data_eff_endcap.append(-1)
            data_eff_endcap_error.append( 0.001 )
            data_eff_endcap_error_total.append( 0.001 )



        if mc_integrals[i] and mc_integrals[i+1]:
            mc_barrel_pass = mc_integrals[i]
            mc_barrel_pass_error = mc_integral_sigmas[i] if mc_integral_sigmas[i] and mc_integral_sigmas[i]/mc_integrals[i] < 1 else False
            mc_barrel_pass_error_bksyst = abs(mc_integral_sigmas_bksyst[i])
            mc_barrel_pass_error_rangesyst = abs(mc_integral_sigmas_rangesyst[i])

            mc_barrel_fail = mc_integrals[i+1]
            mc_barrel_fail_error = mc_integral_sigmas[i+1] if mc_integral_sigmas[i+1]  and mc_integral_sigmas[i+1]/mc_integrals[i+1] < 1 else False
            mc_barrel_fail_error_bksyst = abs(mc_integral_sigmas_bksyst[i+1])
            mc_barrel_fail_error_rangesyst = abs(mc_integral_sigmas_rangesyst[i+1])

            mc_eff_barrel.append(mc_barrel_pass / (mc_barrel_pass + mc_barrel_fail))
            mc_eff_barrel_error.append( get_efficiency_error(mc_barrel_pass, mc_barrel_fail, mc_barrel_pass_error, mc_barrel_fail_error) if mc_barrel_pass_error and mc_barrel_fail_error else 0  )
            mc_eff_barrel_error_total.append( get_efficiency_error(mc_barrel_pass, mc_barrel_fail, qsum(mc_barrel_pass_error, mc_barrel_pass_error_bksyst, mc_barrel_pass_error_rangesyst), qsum(mc_barrel_fail_error, mc_barrel_fail_error_bksyst,mc_barrel_fail_error_rangesyst)) if mc_barrel_pass_error and mc_barrel_fail_error else 0  )
        
        else:
            mc_eff_barrel.append(-1)
            mc_eff_barrel_error.append( 0.001 )
            mc_eff_barrel_error_total.append( 0.001 )



        if mc_integrals[i+2] and mc_integrals[i+3]:
            mc_endcap_pass = mc_integrals[i+2]
            mc_endcap_pass_error = mc_integral_sigmas[i+2] if mc_integral_sigmas[i+2] and mc_integral_sigmas[i+2]/mc_integrals[i+2] < 1 else False
            mc_endcap_pass_error_bksyst = abs(mc_integral_sigmas_bksyst[i+2])
            mc_endcap_pass_error_rangesyst = abs(mc_integral_sigmas_rangesyst[i+2])

            mc_endcap_fail = mc_integrals[i+3]
            mc_endcap_fail_error = mc_integral_sigmas[i+3] if mc_integral_sigmas[i+3] and mc_integral_sigmas[i+3]/mc_integrals[i+3] < 1 else False
            mc_endcap_fail_error_bksyst = abs(mc_integral_sigmas_bksyst[i+3])
            mc_endcap_fail_error_rangesyst = abs(mc_integral_sigmas_rangesyst[i+3])

            mc_eff_endcap.append(mc_endcap_pass / (mc_endcap_pass + mc_endcap_fail))
            mc_eff_endcap_error.append( get_efficiency_error(mc_endcap_pass, mc_endcap_fail, mc_endcap_pass_error, mc_endcap_fail_error) if mc_endcap_pass_error and mc_endcap_fail_error else 0  )
            mc_eff_endcap_error_total.append( get_efficiency_error(mc_endcap_pass, mc_endcap_fail, qsum(mc_endcap_pass_error, mc_endcap_pass_error_bksyst,mc_endcap_pass_error_rangesyst), qsum(mc_endcap_fail_error, mc_endcap_fail_error_bksyst,mc_endcap_fail_error_rangesyst)) if mc_endcap_pass_error and mc_endcap_fail_error else 0  )
        else:
            mc_eff_endcap.append(-1)
            mc_eff_endcap_error.append( 0.001 )
            #mc_eff_endcap_error_up.append( 0.001 )
            mc_eff_endcap_error_total.append( 0.001 )





        if data_eff_barrel[-1] > 0 and mc_eff_barrel[-1] > 0:

            if data_eff_barrel_error[-1] > 0 and mc_eff_barrel_error[-1]:
                scale_factor_barrel.append(data_eff_barrel[-1]/mc_eff_barrel[-1])
                scale_factor_x_barrel.append(ptbins_center_fit[pt_bin_number])
                scale_factor_x_barrel_binw.append(binw_fit[pt_bin_number])
                scale_factor_barrel_error.append( ((data_eff_barrel_error_total[-1]/data_eff_barrel[-1]) |qs| (mc_eff_barrel_error_total[-1]/mc_eff_barrel[-1])) * scale_factor_barrel[-1]   )
            else:
                scale_factor_barrel_black.append(data_eff_barrel[-1]/mc_eff_barrel[-1])
                scale_factor_barrel_black_x.append(ptbins_center_fit[pt_bin_number])
                scale_factor_barrel_black_x_binw.append(binw_fit[pt_bin_number])

        if data_eff_endcap[-1] > 0 and mc_eff_endcap[-1] > 0:

            if data_eff_endcap_error[-1] > 0 and mc_eff_endcap_error[-1] > 0:
                scale_factor_endcap.append(data_eff_endcap[-1]/mc_eff_endcap[-1])
                scale_factor_x_endcap.append(ptbins_center_fit[pt_bin_number])
                scale_factor_x_endcap_binw.append(binw_fit[pt_bin_number])
                scale_factor_endcap_error.append( ((data_eff_endcap_error_total[-1]/data_eff_endcap[-1]) |qs| (mc_eff_endcap_error_total[-1]/mc_eff_endcap[-1])  ) * scale_factor_endcap[-1]   )
            else:
                scale_factor_endcap_black.append(data_eff_endcap[-1]/mc_eff_endcap[-1])
                scale_factor_endcap_black_x.append(ptbins_center_fit[pt_bin_number])
                scale_factor_endcap_black_x_binw.append(binw_fit[pt_bin_number])

        pt_bin_number += 1



    axs = []
    bottom_axs = []
    fig = plt.figure(tight_layout=True, figsize=(15, 10))
    gs = plt.GridSpec(2, 2, height_ratios=[2, 1],  figure=fig)# hspace=0,
    for i in range(2):
        axs.append(fig.add_subplot(gs[i]))
        bottom_axs.append(fig.add_subplot(gs[i+2], sharex=axs[-1]))

    ax_barrel = axs[0]
    ax_endcap = axs[1]
    bottom_ax_barrel = bottom_axs[0]
    bottom_ax_endcap = bottom_axs[1]

    ax_barrel.set_xlabel(r"$p_{T}$ [GeV]")
    ax_barrel.set_ylabel("Efficiency")
    ax_barrel.set_title(r"$|\eta | < 1.4442$")
    ax_endcap.set_title(r"$1.566 < |\eta | < 2.50$")
    bottom_ax_barrel.set_xlabel(r"$p_{T}$ [GeV]")
    bottom_ax_barrel.set_ylabel("Data/MC")

    #print len(ptbins_center_fit), len(data_eff_barrel), len(data_eff_barrel_error), len(binw_fit)
    #print ptbins_center_fit, binw_fit

    #ax_barrel.plot(ptbins_center_fit, data_eff_barrel, "bs")
    ax_barrel.errorbar(ptbins_center_fit, data_eff_barrel, yerr=data_eff_barrel_error_total, xerr=binw_fit, capsize=3, fmt=".", color="blue", ecolor = "blue", label="Data")
    #ax_barrel.plot(ptbins_center_fit, mc_eff_barrel, "rs")
    ax_barrel.errorbar(ptbins_center_fit, mc_eff_barrel, yerr=mc_eff_barrel_error_total, xerr=binw_fit, capsize=3, fmt=".", color="crimson", ecolor = "crimson", label="MC")

    bottom_ax_barrel.errorbar(scale_factor_x_barrel, scale_factor_barrel, yerr=scale_factor_barrel_error , xerr=scale_factor_x_barrel_binw, capsize=3, fmt=".", color="black", ecolor = "black")
    bottom_ax_barrel.errorbar(scale_factor_barrel_black_x, scale_factor_barrel_black, xerr=scale_factor_barrel_black_x_binw, fmt="x", color="black", ecolor="black")


    #ax_endcap.plot(ptbins_center_fit, data_eff_endcap, "bs")
    ax_endcap.errorbar(ptbins_center_fit, data_eff_endcap, yerr=data_eff_endcap_error_total, xerr=binw_fit, capsize=3, fmt=".", color="blue", ecolor = "blue", label="Data")
    #ax_endcap.plot(ptbins_center_fit, mc_eff_endcap, "rs")
    ax_endcap.errorbar(ptbins_center_fit, mc_eff_endcap, yerr=mc_eff_endcap_error_total, xerr=binw_fit, capsize=3, fmt=".", color="crimson", ecolor = "crimson", label="MC")

    bottom_ax_endcap.errorbar(scale_factor_x_endcap, scale_factor_endcap, yerr=scale_factor_endcap_error, xerr=scale_factor_x_endcap_binw, capsize=3, fmt=".", color="black", ecolor = "black")
    bottom_ax_endcap.errorbar(scale_factor_endcap_black_x, scale_factor_endcap_black, xerr=scale_factor_endcap_black_x_binw, fmt="x", color="black", ecolor="black")

    ax_barrel.legend(loc="best")

    ax_barrel.set_ylim([0,1.02])
    ax_endcap.set_ylim([0,1.02])
    ax_barrel.set_yticks([0,0.2,0.4,0.6,0.8,1.0])
    ax_barrel.set_yticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0], minor=True)
    ax_barrel.grid(ls="--", which="major")
    ax_endcap.set_yticks([0,0.2,0.4,0.6,0.8,1.0])
    ax_endcap.set_yticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0], minor=True)
    ax_endcap.grid(ls="--", which="major")

    current_xlim = ax_barrel.get_xlim()
    ax_barrel.set_xticks([3,5,12,20,50,100])
    ax_endcap.set_xticks([3,5,12,20,50,100])
    ax_barrel.set_xlim(current_xlim)
    ax_endcap.set_xlim(current_xlim)

    ec_l = bottom_ax_barrel.get_ylim()
    if ec_l[0] > 1:
        bottom_ax_barrel.set_ylim([0.98,ec_l[1]])
    if ec_l[1] < 1:
        bottom_ax_barrel.set_ylim([ec_l[0],1.02])

    ec_l = bottom_ax_endcap.get_ylim()
    if ec_l[0] > 1:
        bottom_ax_endcap.set_ylim([0.98,ec_l[1]])
    if ec_l[1] < 1:
        bottom_ax_endcap.set_ylim([ec_l[0],1.02])

    #bottom_ax_barrel.hlines(1.0, bottom_ax_barrel.get_xlim()[0], bottom_ax_barrel.get_xlim()[1], linestyles="dashed", color="crimson", alpha=0.7)
    #bottom_ax_endcap.hlines(1.0, bottom_ax_endcap.get_xlim()[0], bottom_ax_endcap.get_xlim()[1], linestyles="dashed", color="crimson")

    bottom_ax_barrel.axhline(1.0, ls="--", color="crimson", alpha=0.7)
    bottom_ax_endcap.axhline(1.0, ls="--", color="crimson", alpha=0.7)
    plt.savefig(filename)













if __name__ == "__main__":
    main()




