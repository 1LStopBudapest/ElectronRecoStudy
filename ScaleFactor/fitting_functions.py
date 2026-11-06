import numpy as np
from numpy.polynomial import Chebyshev
from scipy.stats import crystalball
from scipy.stats import chi2 as get_p_value
import math
from iminuit import Minuit
from iminuit.cost import ExtendedBinnedNLL


from numpy.typing import ArrayLike


from futils import qs

############################
### --- Model Functions
############################

# actual one used by template fit
def cms_shape(xarr, alpha, beta, gamma, peak):
    result = []
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
    return np.array(result)

# Actual one used by template fit
def normalized_background(x, alpha, beta, gamma):

    # the term with peak can be rewritten as
    # exp(peak*gamma) * exp(-x*gamma)
    # so it is just a constant multiplicative factor independent of x - drops out after normalization
    # so for normalized cmsshape, it does not matter - do not use as fit param, fix at 90.

    shape = cms_shape(x, alpha, beta, gamma, 90)
    # Prevent numerical problems
    shape = np.clip(shape, 0, None)
    integral = np.sum(shape)
    if integral <= 0 or not np.isfinite(integral):
        return np.ones_like(x) / len(x)
    return shape / integral

# for Template fit
def background_cdf(x, alpha, beta, gamma):
    pdf = normalized_background(x, alpha, beta, gamma)
    return np.cumsum(pdf)


def baker_cousins_residual(observed, expected):
    observed = np.asarray(observed, dtype=float)
    expected = np.asarray(expected, dtype=float)

    residual = np.zeros_like(observed)

    mask = observed > 0

    residual[mask] = np.sign(observed[mask] - expected[mask]) * np.sqrt(
        2 * (
            observed[mask] * np.log(observed[mask] / expected[mask])
            - (observed[mask] - expected[mask])
        )
    )

    residual[~mask] = -np.sqrt(2 * expected[~mask])

    return residual





############################
### --- Extended Binned NLL for analytic function
############################



def do_minuit_fit(input_x, input_y, p0 = {}, limits = {}):

    iminuit_x = np.array(input_x)
    iminuit_y = np.array(input_y)
    dx = iminuit_x[1]-iminuit_x[0]
    edges = np.concatenate(([iminuit_x[0]-dx/2], iminuit_x+dx/2))

    bin_width = iminuit_x[1] - iminuit_x[0]

    xmin = edges[0]
    xmax = edges[-1]
    def background_cdf(x, c1, c2, c3):
        poly = Chebyshev([1, c1, c2, c3])
        pint = poly.integ()
        t = 2*(x-xmin)/(xmax-xmin) - 1
        norm = pint(1) - pint(-1)
        return (pint(t) - pint(-1)) / norm


    def model_cdf(edges,
                nsig,
                nbkg,
                mu,
                sigma,
                beta,
                m,
                c1,
                c2,
                c3):
        sig = crystalball.cdf(
            edges,
            beta,
            m,
            loc=mu,
            scale=sigma,
        )
        bkg = background_cdf(
            edges,
            c1,
            c2,
            c3,
        )
        return nsig*sig + nbkg*bkg

    cost = ExtendedBinnedNLL(iminuit_y, edges, model_cdf)


    m = Minuit(
        cost,
        nsig=2200 if not "nsig" in p0 else p0["nsig"],
        nbkg=2800 if not 'nbkg' in p0 else p0['nbkg'],
        mu=3.10 if not 'mu' in p0 else p0['mu'],
        sigma=0.07 if not 'sigma' in p0 else p0['sigma'],
        beta=2.0 if not 'beta' in p0 else p0['beta'],
        m=3.0 if not 'm' in p0 else p0['m'],
        c1=0.0 if not 'c1' in p0 else p0['c1'],
        c2=0.0 if not 'c2' in p0 else p0['c2'],
        c3=0.0 if not 'c3' in p0 else p0['c3'],
    )

    #print((0, None) if not "nsig" in limits else (limits["nsig"][0], limits["nsig"][1]))


    m.limits["nsig"] = (0, None) if not "nsig" in limits else (limits["nsig"][0], limits["nsig"][1])
    m.limits["nbkg"] = (0, None) if not 'nbkg' in limits else (limits['nbkg'][0], limits['nbkg'][1])
    m.limits["sigma"] = (0.001, 0.20) if not 'sigma' in limits else (limits['sigma'][0], limits['sigma'][1])
    m.limits["beta"] = (0.5, 10) if not 'beta' in limits else (limits['beta'][0], limits['beta'][1])
    m.limits["m"] = (1.1, 20) if not 'm' in limits else (limits['m'][0], limits['m'][1])

    m.migrad()
    m.hesse()


    # objects for later plotting
    x_sample = np.linspace(edges[0],edges[-1]+bin_width/2,10000)

    def signal_pdf(x, mu, sigma, beta, m):
        return crystalball.pdf(x, beta, m, loc=mu, scale=sigma)

    def background_pdf(x, c1, c2, c3):
        poly = Chebyshev([1, c1, c2, c3])
        t = 2*(x-xmin)/(xmax-xmin) - 1
        norm = poly.integ()(1) - poly.integ()(-1)
        return (2/(xmax-xmin)) * poly(t) / norm


    pars = m.values

    signal = pars["nsig"] * signal_pdf(
        x_sample,
        pars["mu"],
        pars["sigma"],
        pars["beta"],
        pars["m"],
    )

    background = pars["nbkg"] * background_pdf(
        x_sample,
        pars["c1"],
        pars["c2"],
        pars["c3"],
    )

    total = signal + background



    # 4 Get Binned Signal

    def get_signal_binned(edges, nsig, mu, sigma, beta, m):
        cdf = crystalball.cdf(
            edges,
            beta,
            m,
            loc=mu,
            scale=sigma,
        )

        return nsig * np.diff(cdf)

    def signal_binned_errors(edges, minuit):
        names = ["nsig", "mu", "sigma", "beta", "m"]

        values = np.array([minuit.values[name] for name in names], dtype=float)
        errors = np.array([minuit.errors[name] for name in names], dtype=float)

        covariance = np.array(
            [[minuit.covariance[a, b] for b in names] for a in names],
            dtype=float,
        )

        # Numerical Jacobian:
        # J[i, j] = d(signal_bin_i) / d(parameter_j)

        n_bins = len(edges) - 1
        n_pars = len(names)

        J = np.zeros((n_bins, n_pars))

        for j in range(n_pars):
            # Use a small step relative to the fitted parameter error
            h = errors[j] * 1e-3

            # Protect against zero parameter errors
            if h == 0:
                h = max(abs(values[j]) * 1e-6, 1e-8)

            p_plus = values.copy()
            p_minus = values.copy()

            p_plus[j] += h
            p_minus[j] -= h

            f_plus = get_signal_binned(edges, *p_plus)
            f_minus = get_signal_binned(edges, *p_minus)

            J[:, j] = (f_plus - f_minus) / (2 * h)

        # Covariance matrix of the binned signal
        signal_covariance = J @ covariance @ J.T

        # Variance of each bin
        signal_variance = np.diag(signal_covariance)

        # Numerical roundoff can occasionally give tiny negative values
        signal_errors = np.sqrt(np.maximum(signal_variance, 0))

        return signal_errors, signal_covariance


    signal_binned = get_signal_binned(
        edges,
        pars["nsig"],
        pars["mu"],
        pars["sigma"],
        pars["beta"],
        pars["m"],
    )

    uncert_signal_binned, signal_binned_covariance = signal_binned_errors(
        edges,
        m,
    )




    # expected counts/bin from the fitted CDF
    expected = np.diff(model_cdf(edges, *pars))


    model_deviance = abs(iminuit_y - expected)

    # pearson pull would be
    #sigma = np.sqrt(np.maximum(expected, 1e-12))
    #pull = (y - expected) / sigma

    pull = baker_cousins_residual(iminuit_y, expected)



    # 5 Uncertainty Band
    pcov = m.covariance

    y_fit_curve = total * bin_width

    uncert_band = None
    uncert_expected = None
    if pcov is not None:

        
        # 5.1 Uncertainty Band (Smooth Sampling)
        n_samples = 100
        vopts = np.random.multivariate_normal(m.values, pcov, n_samples)
        sampled_ydata = []
        for v in vopts:
            if v[0] <= 0: continue
            if v[1] <= 0: continue
            if v[3] <= 0: continue
            signal_i = v[0] * signal_pdf(
                x_sample,
                v[2],
                v[3],
                v[4],
                v[5],
            )
            background_i = v[1] * background_pdf(
                x_sample,
                v[6],
                v[7],
                v[8],
            )
            total_i = signal_i + background_i
            y_sample = total_i * bin_width
            sampled_ydata.append(y_sample)
        sampled_ydata = np.vstack(sampled_ydata)
        model_uncert = np.nanstd(sampled_ydata, axis=0)

        uncert_band = (y_fit_curve - model_uncert, y_fit_curve + model_uncert)


        # 5.2 Uncertainty Band (Sampling for Binned Expected)
        n_samples = 5000

        samples = np.random.multivariate_normal(
            np.array(m.values),
            np.array(pcov),
            n_samples
        )

        sampled_expected = []

        for v in samples:

            nsig, nbkg, mu, sigma, beta, shape_m, c1, c2, c3 = v

            if nsig <= 0 or nbkg <= 0 or sigma <= 0:
                continue

            expected_i = np.diff(
                model_cdf(
                    edges,
                    nsig,
                    nbkg,
                    mu,
                    sigma,
                    beta,
                    shape_m,
                    c1,
                    c2,
                    c3,
                )
            )

            sampled_expected.append(expected_i)

        sampled_expected = np.asarray(sampled_expected)

        uncert_expected = np.std(sampled_expected, axis=0)


    # distinction: signal, background, and total are smooth curves, evaluated at any x. They are therefore derived from the PDF
    # the expected is binned, it is the actual expected counts in each data bin from the fit. Derived from the CDF

    minuit_info = {
        "m": m,
        "signal": signal,
        "background": background,
        "expected": expected,
        "signal_binned": signal_binned,
        "bin_edges": edges,
        "uncert_signal_binned": uncert_signal_binned,
        "cov_signal_binned": signal_binned_covariance,
        "pull": pull,
        "uncert_band": uncert_band,
        "uncert_expected": uncert_expected,
        "uncert_model_deviance": model_deviance,
        "x_sample": x_sample
    }

    return minuit_info












############################
### --- Template Fitting
############################








def do_template_fit(x_data, y_data, y_mc, config_file, index):
    from iminuit.cost import Template


    x = np.asarray(x_data)
    data = np.asarray(y_data, dtype=float)
    mc = np.asarray(y_mc, dtype=float)

    dx = x[1] - x[0]
    edges = np.r_[
        x - dx/2,
        x[-1] + dx/2
    ]

    # ------------------------------------------------------------
    # Signal template
    # ------------------------------------------------------------
    #
    # Template() expects the raw MC histogram, as a np.histogram:
    # tuple of histogram values, and bin edges
    signal_template_nph = (np.asarray(mc), np.asarray(edges))
    signal_template_norm = mc / np.sum(mc)


    # ------------------------------------------------------------
    # Background model
    # ------------------------------------------------------------
    #
    # normalized_background(x, alpha, beta, gamma) returns
    # sum(background) = 1.
    #
    def background_model(x, nbkg, alpha, beta, gamma):

        return nbkg * background_cdf(    #normalized_background(
            x,
            alpha,
            beta,
            gamma
        )



    # ------------------------------------------------------------
    # Expected counts function
    # ------------------------------------------------------------
    #
    def expected_counts(
        x: ArrayLike, 
        nsig: float,
        nbkg: float,
        alpha: float,
        beta: float,
        gamma: float
    ):
        background_template_norm = normalized_background(
            x, alpha, beta, gamma
        )
        return nsig * signal_template_norm + nbkg * background_template_norm
        


    # ------------------------------------------------------------
    # Limits from config
    # ------------------------------------------------------------
    n_total = np.sum(data)

    if index in config_file["template_fit_p0"].keys():


        if isinstance(config_file["template_fit_p0"][index][0], list):
            nsig_over_mc = config_file["template_fit_p0"][index][0][1]


            nsig_p0 = nsig_over_mc * sum(y_mc)
            if nsig_p0 > n_total: nsig_p0 = 0.95*n_total
            nbkg_p0 = n_total - nsig_p0
        else:
            nsig_p0  = config_file["template_fit_p0"][index][0] if config_file["template_fit_p0"][index][0] else 0.5*n_total
            nbkg_p0  = config_file["template_fit_p0"][index][1] if config_file["template_fit_p0"][index][1] else 0.5*n_total



        alpha_p0 = config_file["template_fit_p0"][index][2]
        beta_p0  = config_file["template_fit_p0"][index][3]
        gamma_p0 = config_file["template_fit_p0"][index][4]


    else:
        nsig_p0=0.5 * n_total
        nbkg_p0=0.5 * n_total
        alpha_p0=60
        beta_p0=0.05
        gamma_p0=0.1


    # ------------------------------------------------------------
    # Build iminuit Template cost
    # ------------------------------------------------------------
    #
    # Component 0: signal MC template
    # Component 1: parametric background
    #
    # Template() expects the raw MC histogram, as a np.histogram,
    # which is a tuple of edges and values.
    # Template will handle normalization internally
    # And also incorporates MC stats into likelihood
    c = Template(
        data, edges,
        (
            signal_template_nph[0],
            background_model,
        ),
    )


    # ------------------------------------------------------------
    # Minuit
    # ------------------------------------------------------------
    m = Minuit(
        c,
        x0=nsig_p0,
        x1_nbkg=nbkg_p0,
        x1_alpha=alpha_p0,
        x1_beta=beta_p0,
        x1_gamma=gamma_p0,
    )

    # Limits
    if index in config_file["template_fit_limits"].keys():

        if config_file["template_fit_limits"][index][0] and isinstance( config_file["template_fit_limits"][index][0][0], str ):  # nsig/mc split, devi by nsig value [1]

            nsig_up = (nsig_over_mc + config_file["template_fit_limits"][index][0][1]) * sum(y_mc)
            if nsig_up > n_total: nsig_up = n_total

            nsig_down = (nsig_over_mc - config_file["template_fit_limits"][index][0][1]) * sum(y_mc)
            if nsig_up > n_total: nsig_up = n_total

            m.limits["x0"] = (nsig_down, nsig_up)
            m.limits["x1_nbkg"]  = (0, None)
        else:
            m.limits["x0"]  = config_file["template_fit_limits"][index][0] if config_file["template_fit_limits"][index][0] else (0, None)
            m.limits["x1_nbkg"]  = config_file["template_fit_limits"][index][1] if config_file["template_fit_limits"][index][1] else (0, None)

        m.limits["x1_alpha"] = config_file["template_fit_limits"][index][2] if config_file["template_fit_limits"][index][2] else (25, 60)
        m.limits["x1_beta"]  = config_file["template_fit_limits"][index][3] if config_file["template_fit_limits"][index][3] else (0.001, 0.15)
        m.limits["x1_gamma"] = config_file["template_fit_limits"][index][4] if config_file["template_fit_limits"][index][4] else (-2,2)
    else:
        m.limits["x0"] = (0, None)
        m.limits["x1_nbkg"] = (0, None)
        m.limits["x1_alpha"] = (25, 60)
        m.limits["x1_beta"] = (0.001, 0.15)
        m.limits["x1_gamma"] = (-2,2)


    #m.errordef = Minuit.LIKELIHOOD
    m.migrad()
    m.hesse()

    pars = np.array([
        m.values["x0"],
        m.values["x1_nbkg"],
        m.values["x1_alpha"],
        m.values["x1_beta"],
        m.values["x1_gamma"],
    ], dtype=float)

    # Plot it

    signal = (
        m.values["x0"] * signal_template_norm
    )

    background = (
        m.values["x1_nbkg"]
        * normalized_background(
            x,
            m.values["x1_alpha"],
            m.values["x1_beta"],
            m.values["x1_gamma"]
        )
    )

    expected = signal + background

    
    # FIXME use this as the global def for baker_cousisns?
    def baker_cousins_pull(data, expected):
        """
        Baker-Cousins signed likelihood residual.

        data     : observed bin counts
        expected : fitted expected bin counts
        """
        data = np.asarray(data, dtype=float)
        expected = np.asarray(expected, dtype=float)
        expected = np.maximum(expected, 1e-12)
        pull = np.zeros_like(data)
        # n > 0
        mask = data > 0
        pull[mask] = (
            np.sign(data[mask] - expected[mask])
            * np.sqrt(
                2.0 * (
                    expected[mask]
                    - data[mask]
                    + data[mask]
                    * np.log(data[mask] / expected[mask])
                )
            )
        )
        # n == 0
        mask_zero = data == 0
        pull[mask_zero] = -np.sqrt(
            2.0 * expected[mask_zero]
        )
        return pull

    pull = baker_cousins_pull(
        data,
        expected,
    )


    # chi2
    pearson_chi2 = np.sum((data - expected)**2 / expected)
    ndof = len(data) - m.nfit
    chi2_over_ndof = pearson_chi2 / ndof


    # Poisson deviance
    mask = data > 0
    deviance = 2 * (
        np.sum(expected - data)
        + np.sum(
            data[mask] * np.log(data[mask] / expected[mask])
        )
    )
    deviance_over_ndof = deviance / ndof



    p_value = get_p_value.sf(deviance, ndof)




    uncert_expected = None
    '''
    pcov = np.asarray(m.covariance)
    if pcov is not None:

        n_samples = 5000

        samples = np.random.multivariate_normal(
            pars,
            pcov,
            size=n_samples
        )

        sampled_expected = []

        for v in samples:

            nsig, nbkg, alpha, beta, gamma = v

            # Respect physical constraints
            if nsig < 0 or nbkg < 0:
                continue

            expected_i = expected_counts(
                x,
                nsig,
                nbkg,
                alpha,
                beta,
                gamma
            )

            sampled_expected.append(expected_i)

        sampled_expected = np.asarray(sampled_expected)

        uncert_expected = np.std(
            sampled_expected,
            axis=0
        )
    '''

    # FIXME: Maybe not? Check what this gives. Consider something else.
    uncert_signal = (
        m.errors["x0"]
        * signal_template_norm
    )



    plotinfo = {
        "m": m,
        "pull": pull,
        "x": x,
        "data": data,
        "mc": mc,
        "uncert_mc": np.sqrt(mc),
        "expected": expected,
        "signal": signal,
        "background": background,
        "cc": (sum(y_data) - m.values["x1_nbkg"], np.sqrt(sum(y_data)) |qs| m.errors["x1_nbkg"]),
        "chi2": pearson_chi2,
        "deviance": deviance,
        "ndof": ndof,
        "p_value": p_value,
        "uncert_expected": uncert_expected,
        "uncert_signal": uncert_signal
    }

    return plotinfo











def do_template_fitBNLL(x_data, y_data, y_mc, config_file, index):
    x = x_data
    data = np.asarray(y_data, dtype=float)
    mc = np.asarray(y_mc, dtype=float)


    # Normalize MC template to unit area
    signal_template = mc / np.sum(mc)


    def expected_counts(
        x,
        x0,
        x1_nbkg,
        x1_alpha,
        x1_beta,
        x1_gamma
    ):
        background_template = normalized_background(
            x, x1_alpha, x1_beta, x1_gamma
        )
        return (
            x0 * signal_template
            + x1_nbkg * background_template
        )


    def poisson_nll(
        x0,
        x1_nbkg,
        x1_alpha,
        x1_beta,
        x1_gamma
    ):
        expected = expected_counts(
            x,
            x0,
            x1_nbkg,
            x1_alpha,
            x1_beta,
            x1_gamma
        )
        # Avoid log(0)
        expected = np.clip(expected, 1e-12, None)
        nll = np.sum(
            expected - data * np.log(expected)
        )
        return nll


    n_total = np.sum(data)

    if index in config_file["template_fit_p0"].keys():


        if isinstance(config_file["template_fit_p0"][index][0], list):
            nsig_over_mc = config_file["template_fit_p0"][index][0][1]


            nsig_p0 = nsig_over_mc * sum(y_mc)
            if nsig_p0 > n_total: nsig_p0 = 0.95*n_total
            nbkg_p0 = n_total - nsig_p0
        else:
            nsig_p0  = config_file["template_fit_p0"][index][0] if config_file["template_fit_p0"][index][0] else 0.5*n_total
            nbkg_p0  = config_file["template_fit_p0"][index][1] if config_file["template_fit_p0"][index][1] else 0.5*n_total



        alpha_p0 = config_file["template_fit_p0"][index][2]
        beta_p0  = config_file["template_fit_p0"][index][3]
        gamma_p0 = config_file["template_fit_p0"][index][4]


    else:
        nsig_p0=0.5 * n_total
        nbkg_p0=0.5 * n_total
        alpha_p0=60
        beta_p0=0.05
        gamma_p0=0.1






    m = Minuit(
        poisson_nll,
        x0=nsig_p0,
        x1_nbkg=nbkg_p0,
        x1_alpha=alpha_p0,
        x1_beta=beta_p0,
        x1_gamma=gamma_p0
    )

    # Limits


    if index in config_file["template_fit_limits"].keys():

        if config_file["template_fit_limits"][index][0] and isinstance( config_file["template_fit_limits"][index][0][0], str ):  # nsig/mc split, devi by nsig value [1]

            nsig_up = (nsig_over_mc + config_file["template_fit_limits"][index][0][1]) * sum(y_mc)
            if nsig_up > n_total: nsig_up = n_total

            nsig_down = (nsig_over_mc - config_file["template_fit_limits"][index][0][1]) * sum(y_mc)
            if nsig_up > n_total: nsig_up = n_total

            m.limits["x0"] = (nsig_down, nsig_up)
            m.limits["x1_nbkg"]  = (0, None)
        else:
            m.limits["x0"]  = config_file["template_fit_limits"][index][0] if config_file["template_fit_limits"][index][0] else (0, None)
            m.limits["x1_nbkg"]  = config_file["template_fit_limits"][index][1] if config_file["template_fit_limits"][index][1] else (0, None)

        m.limits["x1_alpha"] = config_file["template_fit_limits"][index][2] if config_file["template_fit_limits"][index][2] else (25, 60)
        m.limits["x1_beta"]  = config_file["template_fit_limits"][index][3] if config_file["template_fit_limits"][index][3] else (0.001, 0.15)
        m.limits["x1_gamma"] = config_file["template_fit_limits"][index][4] if config_file["template_fit_limits"][index][4] else (-2,2)
    else:
        m.limits["x0"] = (0, None)
        m.limits["x1_nbkg"] = (0, None)
        m.limits["x1_alpha"] = (25, 60)
        m.limits["x1_beta"] = (0.001, 0.15)
        m.limits["x1_gamma"] = (-2,2)


    m.errordef = Minuit.LIKELIHOOD
    m.migrad()
    m.hesse()

    pars = np.array([
        m.values["x0"],
        m.values["x1_nbkg"],
        m.values["x1_alpha"],
        m.values["x1_beta"],
        m.values["x1_gamma"],
    ], dtype=float)

    # Plot it

    expected = expected_counts(
        x,
        m.values["x0"],
        m.values["x1_nbkg"],
        m.values["x1_alpha"],
        m.values["x1_beta"],
        m.values["x1_gamma"]
    )

    signal = (
        m.values["x0"] * signal_template
    )

    background = (
        m.values["x1_nbkg"]
        * normalized_background(
            x,
            m.values["x1_alpha"],
            m.values["x1_beta"],
            m.values["x1_gamma"]
        )
    )


    
    # FIXME use this as the global def for baker_cousisns?
    def baker_cousins_pull(data, expected):
        """
        Baker-Cousins signed likelihood residual.

        data     : observed bin counts
        expected : fitted expected bin counts
        """
        data = np.asarray(data, dtype=float)
        expected = np.asarray(expected, dtype=float)
        expected = np.maximum(expected, 1e-12)
        pull = np.zeros_like(data)
        # n > 0
        mask = data > 0
        pull[mask] = (
            np.sign(data[mask] - expected[mask])
            * np.sqrt(
                2.0 * (
                    expected[mask]
                    - data[mask]
                    + data[mask]
                    * np.log(data[mask] / expected[mask])
                )
            )
        )
        # n == 0
        mask_zero = data == 0
        pull[mask_zero] = -np.sqrt(
            2.0 * expected[mask_zero]
        )
        return pull

    pull = baker_cousins_pull(
        data,
        expected,
    )


    # chi2
    pearson_chi2 = np.sum((data - expected)**2 / expected)
    ndof = len(data) - m.nfit
    chi2_over_ndof = pearson_chi2 / ndof


    # Poisson deviance
    mask = data > 0
    deviance = 2 * (
        np.sum(expected - data)
        + np.sum(
            data[mask] * np.log(data[mask] / expected[mask])
        )
    )
    deviance_over_ndof = deviance / ndof



    p_value = get_p_value.sf(deviance, ndof)




    uncert_expected = None
    '''
    pcov = np.asarray(m.covariance)
    if pcov is not None:

        n_samples = 5000

        samples = np.random.multivariate_normal(
            pars,
            pcov,
            size=n_samples
        )

        sampled_expected = []

        for v in samples:

            x0, x1_nbkg, x1_alpha, x1_beta, x1_gamma = v

            # Respect physical constraints
            if x0 < 0 or x1_nbkg < 0:
                continue

            expected_i = expected_counts(
                x,
                x0,
                x1_nbkg,
                x1_alpha,
                x1_beta,
                x1_gamma
            )

            sampled_expected.append(expected_i)

        sampled_expected = np.asarray(sampled_expected)

        uncert_expected = np.std(
            sampled_expected,
            axis=0
        )
    '''

    uncert_signal = (
        m.errors["x0"]
        * signal_template
    )


    plotinfo = {
        "m": m,
        "pull": pull,
        "x": x,
        "data": data,
        "mc": mc,
        "uncert_mc": np.sqrt(mc),
        "expected": expected,
        "signal": signal,
        "background": background,
        "cc": (sum(y_data) - m.values["x1_nbkg"], np.sqrt(sum(y_data)) |qs| m.errors["x1_nbkg"]),
        "chi2": pearson_chi2,
        "deviance": deviance,
        "ndof": ndof,
        "p_value": p_value,
        "uncert_expected": uncert_expected,
        "uncert_signal": uncert_signal
    }


    return plotinfo