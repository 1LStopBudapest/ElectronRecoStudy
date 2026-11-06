import numpy as np
import os
import numpy as np
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import yaml



def get_parser():
    import argparse
    argParser = argparse.ArgumentParser(description = "Argument parser")
    argParser.add_argument('zdata16pre',help='config_file')
    argParser.add_argument('zdata16post',help='config_file')
    argParser.add_argument('zdata17',help='config_file')
    argParser.add_argument('zdata18',help='config_file')
    argParser.add_argument('jpsi18',help='config_file')
    argParser.add_argument('--fname', nargs='?', const="", default="", type=str)
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

def qsum(*args):
    sm = 0
    for arg in args:
        sm += arg**2
    return np.sqrt(sm)
qs = Infix(qsum)




def main():



    options = get_parser().parse_args()


    with open(options.zdata16pre, 'r') as file:
        z16pre = yaml.safe_load(file)

    with open(options.zdata16post, 'r') as file:
        z16post = yaml.safe_load(file)

    with open(options.zdata17, 'r') as file:
        z17 = yaml.safe_load(file)

    with open(options.zdata18, 'r') as file:
        z18 = yaml.safe_load(file)

    with open(options.jpsi18, 'r') as file:
        j18 = yaml.safe_load(file)



    fname = "" if not options.fname else options.fname + "-"


    makeEfficiencyPlot(f"SF/all_years-{fname}definitive.png", 

                        z16pre["shared_data_minuit_fr_eff"],
                        z16pre["shared_data_minuit_fr_eff_stat_err"],
                        z16pre["shared_mc_cutcount_fr_eff"],
                        z16pre["shared_mc_cutcount_fr_eff_total_err"],
                        z16pre["ptbins_center_fit"],
                        z16pre["binw_fit"],

                        z16post["shared_data_minuit_fr_eff"],
                        z16post["shared_data_minuit_fr_eff_stat_err"],
                        z16post["shared_mc_cutcount_fr_eff"],
                        z16post["shared_mc_cutcount_fr_eff_total_err"],
                        z16post["ptbins_center_fit"],
                        z16post["binw_fit"],

                        z17["shared_data_minuit_fr_eff"],
                        z17["shared_data_minuit_fr_eff_stat_err"],
                        z17["shared_mc_cutcount_fr_eff"],
                        z17["shared_mc_cutcount_fr_eff_total_err"],
                        z17["ptbins_center_fit"],
                        z17["binw_fit"],

                        z18["shared_data_minuit_fr_eff"],
                        z18["shared_data_minuit_fr_eff_stat_err"],
                        z18["shared_mc_cutcount_fr_eff"],
                        z18["shared_mc_cutcount_fr_eff_total_err"],
                        z18["ptbins_center_fit"],
                        z18["binw_fit"],

                        j18["shared_data_minuit_eff"][(j18["skipped_lowpt_bins"])*2:],
                        j18["shared_data_minuit_eff_total_err"][(j18["skipped_lowpt_bins"])*2:],
                        j18["shared_mc_cutcount_fr_eff"][(j18["skipped_lowpt_bins"])*2:],
                        j18["shared_mc_cutcount_fr_eff_total_err"][(j18["skipped_lowpt_bins"])*2:],
                        j18["ptbins_center_fit"],
                        j18["binw_fit"],

                        )




def getSFerrorbar(data, mc, data_err, mc_err, rho=0.5):
    sf = data / mc
    data_relerr = data_err/data
    mc_relerr = mc_err/mc
    rel_err = np.sqrt(  data_relerr**2  + mc_relerr**2   - 2*rho*data_relerr*mc_relerr  )
    return rel_err * sf


def makeEfficiencyPlot(filename, 

    z16pre_data_eff, z16pre_data_eff_error, z16pre_mc_eff, z16pre_mc_eff_error,
    z16pre_ptbins_center_fit, z16pre_binw_fit, 

    z16post_data_eff, z16post_data_eff_error, z16post_mc_eff, z16post_mc_eff_error,
    z16post_ptbins_center_fit, z16post_binw_fit, 

    z17_data_eff, z17_data_eff_error, z17_mc_eff, z17_mc_eff_error,
    z17_ptbins_center_fit, z17_binw_fit, 


    z18_data_eff, z18_data_eff_error, z18_mc_eff, z18_mc_eff_error,
    z18_ptbins_center_fit, z18_binw_fit, 


    jpsi_data_eff, jpsi_data_eff_error, jpsi_mc_eff, jpsi_mc_eff_error,
    jpsi_ptbins_center_fit, jpsi_binw_fit,


    ):




    def getSF(
        z_data_eff, z_data_eff_error, z_mc_eff, z_mc_eff_error,
        z_ptbins_center_fit, z_binw_fit, 
        ):
        
        z_data_eff_barrel = z_data_eff[0::2]
        z_data_eff_barrel_error_total = z_data_eff_error[0::2]
        z_mc_eff_barrel = z_mc_eff[0::2]
        z_mc_eff_barrel_error_total = z_mc_eff_error[0::2]
        z_data_eff_endcap = z_data_eff[1::2]
        z_data_eff_endcap_error_total = z_data_eff_error[1::2]
        z_mc_eff_endcap = z_mc_eff[1::2]
        z_mc_eff_endcap_error_total = z_mc_eff_error[1::2]

        z_scale_factor_barrel = []
        z_scale_factor_barrel_error = []
        z_scale_factor_x_barrel = []
        z_scale_factor_x_barrel_binw = []
        z_scale_factor_barrel_black = []
        z_scale_factor_barrel_black_x = []
        z_scale_factor_barrel_black_x_binw = []

        z_scale_factor_endcap = []
        z_scale_factor_endcap_error = []
        z_scale_factor_x_endcap = []
        z_scale_factor_x_endcap_binw = []
        z_scale_factor_endcap_black = []
        z_scale_factor_endcap_black_x = []
        z_scale_factor_endcap_black_x_binw = []


        pt_bin_number = 0
        for i in range(len(z_data_eff_barrel)):
            if z_data_eff_barrel[i] > 0 and z_mc_eff_barrel[i] > 0:

                if z_data_eff_barrel_error_total[i] > 0 and z_mc_eff_barrel_error_total[i]:
                    z_scale_factor_barrel.append(z_data_eff_barrel[i]/z_mc_eff_barrel[i])
                    z_scale_factor_x_barrel.append(z_ptbins_center_fit[pt_bin_number])
                    z_scale_factor_x_barrel_binw.append(z_binw_fit[pt_bin_number])
                    z_scale_factor_barrel_error.append( getSFerrorbar(z_data_eff_barrel[i], z_mc_eff_barrel[i], z_data_eff_barrel_error_total[i],  z_mc_eff_barrel_error_total[i] )   )  
                else:
                    z_scale_factor_barrel_black.append(z_data_eff_barrel[i]/z_mc_eff_barrel[i])
                    z_scale_factor_barrel_black_x.append(z_ptbins_center_fit[pt_bin_number])
                    z_scale_factor_barrel_black_x_binw.append(z_binw_fit[pt_bin_number])

            if z_data_eff_endcap[i] > 0 and z_mc_eff_endcap[i] > 0:

                if z_data_eff_endcap_error_total[i] > 0 and z_mc_eff_endcap_error_total[i] > 0:
                    z_scale_factor_endcap.append(z_data_eff_endcap[i]/z_mc_eff_endcap[i])
                    z_scale_factor_x_endcap.append(z_ptbins_center_fit[pt_bin_number])
                    z_scale_factor_x_endcap_binw.append(z_binw_fit[pt_bin_number])
                    z_scale_factor_endcap_error.append(  getSFerrorbar(z_data_eff_endcap[i], z_mc_eff_endcap[i], z_data_eff_endcap_error_total[i],  z_mc_eff_endcap_error_total[i])    )
                else:
                    z_scale_factor_endcap_black.append(z_data_eff_endcap[i]/z_mc_eff_endcap[i])
                    z_scale_factor_endcap_black_x.append(z_ptbins_center_fit[pt_bin_number])
                    z_scale_factor_endcap_black_x_binw.append(z_binw_fit[pt_bin_number])

            pt_bin_number += 1

        return z_scale_factor_barrel, z_scale_factor_barrel_error, z_scale_factor_x_barrel, z_scale_factor_x_barrel_binw, z_scale_factor_barrel_black, z_scale_factor_barrel_black_x, z_scale_factor_barrel_black_x_binw, z_scale_factor_endcap, z_scale_factor_endcap_error, z_scale_factor_x_endcap, z_scale_factor_x_endcap_binw, z_scale_factor_endcap_black, z_scale_factor_endcap_black_x, z_scale_factor_endcap_black_x_binw




    z_16pre_scale_factor_barrel, z_16pre_scale_factor_barrel_error, z_16pre_scale_factor_x_barrel, z_16pre_scale_factor_x_barrel_binw, z_16pre_scale_factor_barrel_black, z_16pre_scale_factor_barrel_black_x, z_16pre_scale_factor_barrel_black_x_binw, z_16pre_scale_factor_endcap, z_16pre_scale_factor_endcap_error, z_16pre_scale_factor_x_endcap, z_16pre_scale_factor_x_endcap_binw, z_16pre_scale_factor_endcap_black, z_16pre_scale_factor_endcap_black_x, z_16pre_scale_factor_endcap_black_x_binw = getSF(z16pre_data_eff, z16pre_data_eff_error, z16pre_mc_eff, z16pre_mc_eff_error, z16pre_ptbins_center_fit, z16pre_binw_fit)

    z_16post_scale_factor_barrel, z_16post_scale_factor_barrel_error, z_16post_scale_factor_x_barrel, z_16post_scale_factor_x_barrel_binw, z_16post_scale_factor_barrel_black, z_16post_scale_factor_barrel_black_x, z_16post_scale_factor_barrel_black_x_binw, z_16post_scale_factor_endcap, z_16post_scale_factor_endcap_error, z_16post_scale_factor_x_endcap, z_16post_scale_factor_x_endcap_binw, z_16post_scale_factor_endcap_black, z_16post_scale_factor_endcap_black_x, z_16post_scale_factor_endcap_black_x_binw = getSF(z16post_data_eff, z16post_data_eff_error, z16post_mc_eff, z16post_mc_eff_error, z16post_ptbins_center_fit, z16post_binw_fit)

    z_17_scale_factor_barrel, z_17_scale_factor_barrel_error, z_17_scale_factor_x_barrel, z_17_scale_factor_x_barrel_binw, z_17_scale_factor_barrel_black, z_17_scale_factor_barrel_black_x, z_17_scale_factor_barrel_black_x_binw, z_17_scale_factor_endcap, z_17_scale_factor_endcap_error, z_17_scale_factor_x_endcap, z_17_scale_factor_x_endcap_binw, z_17_scale_factor_endcap_black, z_17_scale_factor_endcap_black_x, z_17_scale_factor_endcap_black_x_binw = getSF(z17_data_eff, z17_data_eff_error, z17_mc_eff, z17_mc_eff_error, z17_ptbins_center_fit, z17_binw_fit)

    z_18_scale_factor_barrel, z_18_scale_factor_barrel_error, z_18_scale_factor_x_barrel, z_18_scale_factor_x_barrel_binw, z_18_scale_factor_barrel_black, z_18_scale_factor_barrel_black_x, z_18_scale_factor_barrel_black_x_binw, z_18_scale_factor_endcap, z_18_scale_factor_endcap_error, z_18_scale_factor_x_endcap, z_18_scale_factor_x_endcap_binw, z_18_scale_factor_endcap_black, z_18_scale_factor_endcap_black_x, z_18_scale_factor_endcap_black_x_binw = getSF(z18_data_eff, z18_data_eff_error, z18_mc_eff, z18_mc_eff_error, z18_ptbins_center_fit, z18_binw_fit)


    j_18_scale_factor_barrel, j_18_scale_factor_barrel_error, j_18_scale_factor_x_barrel, j_18_scale_factor_x_barrel_binw, j_18_scale_factor_barrel_black, j_18_scale_factor_barrel_black_x, j_18_scale_factor_barrel_black_x_binw, j_18_scale_factor_endcap, j_18_scale_factor_endcap_error, j_18_scale_factor_x_endcap, j_18_scale_factor_x_endcap_binw, j_18_scale_factor_endcap_black, j_18_scale_factor_endcap_black_x, j_18_scale_factor_endcap_black_x_binw = getSF(jpsi_data_eff, jpsi_data_eff_error, jpsi_mc_eff, jpsi_mc_eff_error, jpsi_ptbins_center_fit, jpsi_binw_fit)

    

    axs = []
    bottom_axs = []
    fig = plt.figure(tight_layout=True, figsize=(15, 10))
    gs = plt.GridSpec(2, 2, height_ratios=[3, 1],  figure=fig)# hspace=0,
    for i in range(2):
        axs.append(fig.add_subplot(gs[i]))
        bottom_axs.append(fig.add_subplot(gs[i+2], sharex=axs[-1]))

    ax_barrel = axs[0]
    ax_endcap = axs[1]
    bottom_ax_barrel = bottom_axs[0]
    bottom_ax_endcap = bottom_axs[1]


    ax_barrel.set_xlabel(r"$p_{T}$ [GeV]", fontsize=14)
    ax_barrel.set_ylabel("Data/MC", fontsize=15)
    ax_endcap.set_xlabel(r"$p_{T}$ [GeV]", fontsize=14)
    ax_endcap.set_ylabel("Data/MC", fontsize=14)
    ax_barrel.set_title(r"$|\eta | < 1.4442$", fontsize=15)
    ax_endcap.set_title(r"$1.566 < |\eta | < 2.50$", fontsize=15)
    bottom_ax_barrel.set_xlabel(r"$p_{T}$ [GeV]", fontsize=14)
    bottom_ax_barrel.set_ylabel("Year/UL2018", fontsize=15)
    bottom_ax_endcap.set_xlabel(r"$p_{T}$ [GeV]", fontsize=14)
    bottom_ax_endcap.set_ylabel("Year/UL2018", fontsize=14)



    ratio17_barrel = []
    ratio17_endcap = []
    ratio16post_barrel = []
    ratio16post_endcap = []
    ratio16pre_barrel = []
    ratio16pre_endcap = []

    for i in range(len(z_18_scale_factor_barrel)):
        ratio17_barrel.append(   z_17_scale_factor_barrel[i] / z_18_scale_factor_barrel[i]     )
        ratio17_endcap.append(   z_17_scale_factor_endcap[i] / z_18_scale_factor_endcap[i]     )
        ratio16post_barrel.append(   z_16post_scale_factor_barrel[i] / z_18_scale_factor_barrel[i]     )
        ratio16post_endcap.append(   z_16post_scale_factor_endcap[i] / z_18_scale_factor_endcap[i]     )
        ratio16pre_barrel.append(   z_16pre_scale_factor_barrel[i] / z_18_scale_factor_barrel[i]     )
        ratio16pre_endcap.append(   z_16pre_scale_factor_endcap[i] / z_18_scale_factor_endcap[i]     )





    # stagger the points in the bins on the x axis  binc - 25%, binc -12.5%, binc, binc+12.5%, binc+25%
    for i in range(len(z_16pre_scale_factor_x_barrel)):
        z_16pre_scale_factor_x_barrel[i] -= z_16pre_scale_factor_x_barrel_binw[i] / 2
        z_16pre_scale_factor_x_endcap[i] -= z_16pre_scale_factor_x_barrel_binw[i] / 2
        z_16post_scale_factor_x_barrel[i] -= z_16post_scale_factor_x_barrel_binw[i] / 4
        z_16post_scale_factor_x_endcap[i] -= z_16post_scale_factor_x_barrel_binw[i] / 4
        z_18_scale_factor_x_barrel[i] += z_18_scale_factor_x_barrel_binw[i] / 4
        z_18_scale_factor_x_endcap[i] += z_18_scale_factor_x_barrel_binw[i] / 4

    for i in range(len(j_18_scale_factor_x_barrel)):
        j_18_scale_factor_x_barrel[i] += j_18_scale_factor_x_barrel_binw[i] / 2
        j_18_scale_factor_x_endcap[i] += j_18_scale_factor_x_barrel_binw[i] / 2

    z_16pre_scale_factor_x_barrel_binw = [[binw*0.5 for binw in z_16pre_scale_factor_x_barrel_binw], [binw*1.5 for binw in z_16pre_scale_factor_x_barrel_binw]]
    z_16pre_scale_factor_x_endcap_binw = [[binw*0.5 for binw in z_16pre_scale_factor_x_endcap_binw], [binw*1.5 for binw in z_16pre_scale_factor_x_endcap_binw]]

    z_16post_scale_factor_x_barrel_binw = [[binw*0.75 for binw in z_16post_scale_factor_x_barrel_binw], [binw*1.25 for binw in z_16post_scale_factor_x_barrel_binw]]
    z_16post_scale_factor_x_endcap_binw = [[binw*0.75 for binw in z_16post_scale_factor_x_endcap_binw], [binw*1.25 for binw in z_16post_scale_factor_x_endcap_binw]]


    z_18_scale_factor_x_barrel_binw = [[binw*1.25 for binw in z_18_scale_factor_x_barrel_binw], [binw*0.75 for binw in z_18_scale_factor_x_barrel_binw]]
    z_18_scale_factor_x_endcap_binw = [[binw*1.25 for binw in z_18_scale_factor_x_endcap_binw], [binw*0.75 for binw in z_18_scale_factor_x_endcap_binw]]

    j_18_scale_factor_x_barrel_binw = [[binw*1.5 for binw in j_18_scale_factor_x_barrel_binw], [binw*0.5 for binw in j_18_scale_factor_x_barrel_binw]]
    j_18_scale_factor_x_endcap_binw = [[binw*1.5 for binw in j_18_scale_factor_x_endcap_binw], [binw*0.5 for binw in j_18_scale_factor_x_endcap_binw]]


    # 16pre
    ax_barrel.errorbar(z_16pre_scale_factor_x_barrel, z_16pre_scale_factor_barrel, yerr=z_16pre_scale_factor_barrel_error, xerr=z_16pre_scale_factor_x_barrel_binw, capsize=3, fmt=".", color="black", ecolor = "black", label="UL2016preVFP")
    ax_barrel.errorbar(z_16pre_scale_factor_barrel_black_x, z_16pre_scale_factor_barrel_black, xerr=z_16pre_scale_factor_barrel_black_x_binw, fmt="x", color="chocolate", ecolor="black")

    ax_endcap.errorbar(z_16pre_scale_factor_x_endcap, z_16pre_scale_factor_endcap, yerr=z_16pre_scale_factor_endcap_error , xerr=z_16pre_scale_factor_x_endcap_binw, capsize=3, fmt=".", color="black", ecolor = "black")
    ax_endcap.errorbar(z_16pre_scale_factor_endcap_black_x, z_16pre_scale_factor_endcap_black, xerr=z_16pre_scale_factor_endcap_black_x_binw, fmt="x", color="chocolate", ecolor="black")


    ax_barrel.errorbar(z_16post_scale_factor_x_barrel, z_16post_scale_factor_barrel, yerr=z_16post_scale_factor_barrel_error, xerr=z_16post_scale_factor_x_barrel_binw, capsize=3, fmt=".", color="blue", ecolor = "blue", label="UL2016postVFP")
    ax_barrel.errorbar(z_16post_scale_factor_barrel_black_x, z_16post_scale_factor_barrel_black, xerr=z_16post_scale_factor_barrel_black_x_binw, fmt="x", color="chocolate", ecolor="blue")

    ax_endcap.errorbar(z_16post_scale_factor_x_endcap, z_16post_scale_factor_endcap, yerr=z_16post_scale_factor_endcap_error, xerr=z_16post_scale_factor_x_endcap_binw, capsize=3, fmt=".", color="blue", ecolor = "blue")
    ax_endcap.errorbar(z_16post_scale_factor_endcap_black_x, z_16post_scale_factor_endcap_black, xerr=z_16post_scale_factor_endcap_black_x_binw, fmt="x", color="chocolate", ecolor="blue")


    ax_barrel.errorbar(z_17_scale_factor_x_barrel, z_17_scale_factor_barrel, yerr=z_17_scale_factor_barrel_error, xerr=z_17_scale_factor_x_barrel_binw, capsize=3, fmt=".", color="red", ecolor = "red", label="UL2017")
    ax_barrel.errorbar(z_17_scale_factor_barrel_black_x, z_17_scale_factor_barrel_black, xerr=z_17_scale_factor_barrel_black_x_binw, fmt="x", color="chocolate", ecolor="red")

    ax_endcap.errorbar(z_17_scale_factor_x_endcap, z_17_scale_factor_endcap, yerr=z_17_scale_factor_endcap_error, xerr=z_17_scale_factor_x_endcap_binw, capsize=3, fmt=".", color="red", ecolor = "red")
    ax_endcap.errorbar(z_17_scale_factor_endcap_black_x, z_17_scale_factor_endcap_black, xerr=z_17_scale_factor_endcap_black_x_binw, fmt="x", color="chocolate", ecolor="red")


    ax_barrel.errorbar(z_18_scale_factor_x_barrel, z_18_scale_factor_barrel, yerr=z_18_scale_factor_barrel_error, xerr=z_18_scale_factor_x_barrel_binw, capsize=3, fmt=".", color="magenta", ecolor = "magenta", label="UL2018")
    ax_barrel.errorbar(z_18_scale_factor_barrel_black_x, z_18_scale_factor_barrel_black, xerr=z_18_scale_factor_barrel_black_x_binw, fmt="x", color="chocolate", ecolor="magenta")

    ax_endcap.errorbar(z_18_scale_factor_x_endcap, z_18_scale_factor_endcap, yerr=z_18_scale_factor_endcap_error, xerr=z_18_scale_factor_x_endcap_binw, capsize=3, fmt=".", color="magenta", ecolor = "magenta")
    ax_endcap.errorbar(z_18_scale_factor_endcap_black_x, z_18_scale_factor_endcap_black, xerr=z_18_scale_factor_endcap_black_x_binw, fmt="x", color="chocolate", ecolor="magenta")



    ax_barrel.errorbar(j_18_scale_factor_x_barrel, j_18_scale_factor_barrel, yerr=j_18_scale_factor_barrel_error, xerr=j_18_scale_factor_x_barrel_binw, capsize=3, fmt=".", color="orange", ecolor = "orange", label="JPsi UL2018")
    ax_barrel.errorbar(j_18_scale_factor_barrel_black_x, j_18_scale_factor_barrel_black, xerr=j_18_scale_factor_barrel_black_x_binw, fmt="x", color="chocolate", ecolor="orange")

    ax_endcap.errorbar(j_18_scale_factor_x_endcap, j_18_scale_factor_endcap, yerr=j_18_scale_factor_endcap_error, xerr=j_18_scale_factor_x_endcap_binw, capsize=3, fmt=".", color="orange", ecolor = "orange")
    ax_endcap.errorbar(j_18_scale_factor_endcap_black_x, j_18_scale_factor_endcap_black, xerr=j_18_scale_factor_endcap_black_x_binw, fmt="x", color="chocolate", ecolor="orange")


    bottom_ax_barrel.errorbar(z_16pre_scale_factor_x_barrel, ratio16pre_barrel, xerr=z_16pre_scale_factor_x_barrel_binw, capsize=3, fmt=".", color="black", ecolor = "black")
    bottom_ax_endcap.errorbar(z_16pre_scale_factor_x_endcap, ratio16pre_endcap, xerr=z_16pre_scale_factor_x_endcap_binw, capsize=3, fmt=".", color="black", ecolor = "black")

    bottom_ax_barrel.errorbar(z_16post_scale_factor_x_barrel, ratio16post_barrel, xerr=z_16post_scale_factor_x_barrel_binw, capsize=3, fmt=".", color="blue", ecolor = "blue")
    bottom_ax_endcap.errorbar(z_16post_scale_factor_x_endcap, ratio16post_endcap, xerr=z_16post_scale_factor_x_endcap_binw, capsize=3, fmt=".", color="blue", ecolor = "blue")

    bottom_ax_barrel.errorbar(z_17_scale_factor_x_barrel, ratio17_barrel, xerr=z_17_scale_factor_x_barrel_binw, capsize=3, fmt=".", color="red", ecolor = "red")
    bottom_ax_endcap.errorbar(z_17_scale_factor_x_endcap, ratio17_endcap, xerr=z_17_scale_factor_x_endcap_binw, capsize=3, fmt=".", color="red", ecolor = "red")


    ax_barrel.set_xlim([1,102])
    ax_endcap.set_xlim([1,102])

    ax_barrel.set_xticks([3,5,12,20,50,100])
    ax_endcap.set_xticks([3,5,12,20,50,100])

    ax_barrel.grid(ls="--", which="major")
    ax_endcap.grid(ls="--", which="major")


    bottom_ax_barrel.grid(ls="--", which="major")
    bottom_ax_endcap.grid(ls="--", which="major")

    ax_barrel.legend(loc="best", fontsize=17)
    
    ec_l = ax_barrel.get_ylim()
    if ec_l[0] > 1:
        ax_barrel.set_ylim([0.98,ec_l[1]])
    if ec_l[1] < 1:
        ax_barrel.set_ylim([ec_l[0],1.02])

    ec_l = ax_endcap.get_ylim()
    if ec_l[0] > 1:
        ax_endcap.set_ylim([0.98,ec_l[1]])
    if ec_l[1] < 1:
        ax_endcap.set_ylim([ec_l[0],1.02])
    
    #bottom_ax_barrel.hlines(1.0, bottom_ax_barrel.get_xlim()[0], bottom_ax_barrel.get_xlim()[1], linestyles="dashed", color="crimson", alpha=0.7)
    #bottom_ax_endcap.hlines(1.0, bottom_ax_endcap.get_xlim()[0], bottom_ax_endcap.get_xlim()[1], linestyles="dashed", color="crimson")


    ax_barrel.axhline(1.0, ls="--", color="crimson", alpha=0.7)
    ax_endcap.axhline(1.0, ls="--", color="crimson", alpha=0.7)
    plt.savefig(filename)




if __name__ == "__main__":
    main()

