import numpy as np
import os
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import yaml



def get_parser():
    import argparse
    argParser = argparse.ArgumentParser(description = "Argument parser")
    argParser.add_argument('jpsidata',help='result_file')
    argParser.add_argument('zdata',help='result_file')
    argParser.add_argument('--fname1', nargs='?', const="JPsi", default="JPsi", type=str)
    argParser.add_argument('--fname2', nargs='?', const="Z18", default="Z18", type=str)
    argParser.add_argument('--import1', nargs='?', const="False", default="False", type=str)
    argParser.add_argument('--import2', nargs='?', const="False", default="False", type=str)
    argParser.add_argument('--fr1', nargs='?', const="False", default="False", type=str)
    argParser.add_argument('--fr2', nargs='?', const="True", default="True", type=str)
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

    with open(options.jpsidata, 'r') as file:
        j = yaml.safe_load(file)

    with open(options.zdata, 'r') as file:
        z = yaml.safe_load(file)


    use_import_for_first =  True if options.import1 == "True" else False
    use_import_for_second = True if options.import2 == "True" else False

    use_fr_for_first =  True if options.fr1 == "True" else False
    use_fr_for_second = True if options.fr2 == "True" else False



    first_data_eff = "shared_data_minuit_eff"
    first_data_eff_error = "shared_data_minuit_eff_total_err"
    first_mc_eff = "shared_mc_cutcount_eff"
    first_mc_eff_error = "shared_mc_cutcount_eff_total_err"
    if use_import_for_first:
        first_data_eff = "data_imported_eff"
        first_data_eff_error = "data_imported_eff_total_err"
        first_mc_eff = "shared_mc_cutcount_fr_eff"
        first_mc_eff_error = "shared_mc_cutcount_fr_eff_total_err"
    if use_fr_for_first:
        first_data_eff = "shared_data_minuit_fr_eff"
        first_data_eff_error = "shared_data_minuit_fr_eff_stat_err"
        first_mc_eff = "shared_mc_cutcount_fr_eff"
        first_mc_eff_error = "shared_mc_cutcount_fr_eff_total_err"



    second_data_eff = "shared_data_minuit_eff"
    second_data_eff_error = "shared_data_minuit_eff_total_err"
    second_mc_eff = "shared_mc_cutcount_eff"
    second_mc_eff_error = "shared_mc_cutcount_eff_total_err"
    if use_import_for_second:
        second_data_eff = "data_imported_eff"
        second_data_eff_error = "data_imported_eff_total_err"
        second_mc_eff = "shared_mc_cutcount_fr_eff"
        second_mc_eff_error = "shared_mc_cutcount_fr_eff_total_err"
    if use_fr_for_second:
        second_data_eff = "shared_data_minuit_fr_eff"
        second_data_eff_error = "shared_data_minuit_fr_eff_stat_err"
        second_mc_eff = "shared_mc_cutcount_fr_eff"
        second_mc_eff_error = "shared_mc_cutcount_fr_eff_total_err"



    makeEfficiencyPlot("SF/%svs%s_definitive.png" % (options.fname1, options.fname2), 
                        j[first_data_eff][(j["skipped_lowpt_bins"])*2:],
                        j[first_data_eff_error][(j["skipped_lowpt_bins"])*2:],
                        j[first_mc_eff][(j["skipped_lowpt_bins"])*2:],
                        j[first_mc_eff_error][(j["skipped_lowpt_bins"])*2:],
                        j["ptbins_center_fit"],
                        j["binw_fit"],
                        z[second_data_eff][(z["skipped_lowpt_bins"])*2:],
                        z[second_data_eff_error][(z["skipped_lowpt_bins"])*2:],
                        z[second_mc_eff][(z["skipped_lowpt_bins"])*2:],
                        z[second_mc_eff_error][(z["skipped_lowpt_bins"])*2:],
                        z["ptbins_center_fit"],
                        z["binw_fit"],
                        options.fname1,
                        options.fname2
                        )



def cutoff_errorbar_at_one(errors, efficiencies):
    newerrors = [[],[]]
    use = False
    for error, eff in zip(errors, efficiencies):
        newerrors[0].append(error)
        if error + eff > 1:
            use = True
            print("Truncated errorbar at 1! FIXME!")
            newerrors[1].append(1 - eff)
        else:
            newerrors[1].append(error)

    if use:
        return newerrors
    return errors






def getSFerrorbar(data, mc, data_err, mc_err, rho=0.5):
    sf = data / mc
    data_relerr = data_err/data
    mc_relerr = mc_err/mc
    rel_err = np.sqrt(  data_relerr**2  + mc_relerr**2   - 2*rho*data_relerr*mc_relerr  )
    return rel_err * sf



def makeEfficiencyPlot(filename, 
    jpsi_data_eff, jpsi_data_eff_error, jpsi_mc_eff, jpsi_mc_eff_error,
    jpsi_ptbins_center_fit, jpsi_binw_fit,

    z_data_eff, z_data_eff_error, z_mc_eff, z_mc_eff_error,
    z_ptbins_center_fit, z_binw_fit, 
    
    fname1, fname2):


    print (fname1, fname2)


    jpsi_data_eff_barrel = jpsi_data_eff[0::2]
    jpsi_data_eff_barrel_error_total = jpsi_data_eff_error[0::2]
    jpsi_mc_eff_barrel = jpsi_mc_eff[0::2]
    jpsi_mc_eff_barrel_error_total = jpsi_mc_eff_error[0::2]
    jpsi_data_eff_endcap = jpsi_data_eff[1::2]
    jpsi_data_eff_endcap_error_total = jpsi_data_eff_error[1::2]
    jpsi_mc_eff_endcap = jpsi_mc_eff[1::2]
    jpsi_mc_eff_endcap_error_total = jpsi_mc_eff_error[1::2]

    jpsi_scale_factor_barrel = []
    jpsi_scale_factor_barrel_error = []
    jpsi_scale_factor_x_barrel = []
    jpsi_scale_factor_x_barrel_binw = []
    jpsi_scale_factor_barrel_black = []
    jpsi_scale_factor_barrel_black_x = []
    jpsi_scale_factor_barrel_black_x_binw = []

    jpsi_scale_factor_endcap = []
    jpsi_scale_factor_endcap_error = []
    jpsi_scale_factor_x_endcap = []
    jpsi_scale_factor_x_endcap_binw = []
    jpsi_scale_factor_endcap_black = []
    jpsi_scale_factor_endcap_black_x = []
    jpsi_scale_factor_endcap_black_x_binw = []

    pt_bin_number = 0
    for i in range(len(jpsi_data_eff_barrel)):
        if jpsi_data_eff_barrel[i] > 0 and jpsi_mc_eff_barrel[i] > 0:

            if jpsi_data_eff_barrel_error_total[i] > 0 and jpsi_mc_eff_barrel_error_total[i]:
                jpsi_scale_factor_barrel.append(jpsi_data_eff_barrel[i]/jpsi_mc_eff_barrel[i])
                jpsi_scale_factor_x_barrel.append(jpsi_ptbins_center_fit[pt_bin_number])
                jpsi_scale_factor_x_barrel_binw.append(jpsi_binw_fit[pt_bin_number])
                jpsi_scale_factor_barrel_error.append( getSFerrorbar(jpsi_data_eff_barrel[i], jpsi_mc_eff_barrel[i], jpsi_data_eff_barrel_error_total[i],  jpsi_mc_eff_barrel_error_total[i]   )  )
            else:
                jpsi_scale_factor_barrel_black.append(jpsi_data_eff_barrel[i]/jpsi_mc_eff_barrel[i])
                jpsi_scale_factor_barrel_black_x.append(jpsi_ptbins_center_fit[pt_bin_number])
                jpsi_scale_factor_barrel_black_x_binw.append(jpsi_binw_fit[pt_bin_number])

        if jpsi_data_eff_endcap[i] > 0 and jpsi_mc_eff_endcap[i] > 0:

            if jpsi_data_eff_endcap_error_total[i] > 0 and jpsi_mc_eff_endcap_error_total[i] > 0:
                jpsi_scale_factor_endcap.append(jpsi_data_eff_endcap[i]/jpsi_mc_eff_endcap[i])
                jpsi_scale_factor_x_endcap.append(jpsi_ptbins_center_fit[pt_bin_number])
                jpsi_scale_factor_x_endcap_binw.append(jpsi_binw_fit[pt_bin_number])
                jpsi_scale_factor_endcap_error.append( getSFerrorbar(jpsi_data_eff_endcap[i], jpsi_mc_eff_endcap[i], jpsi_data_eff_endcap_error_total[i],  jpsi_mc_eff_endcap_error_total[i])     )
            else:
                jpsi_scale_factor_endcap_black.append(jpsi_data_eff_endcap[i]/jpsi_mc_eff_endcap[i])
                jpsi_scale_factor_endcap_black_x.append(jpsi_ptbins_center_fit[pt_bin_number])
                jpsi_scale_factor_endcap_black_x_binw.append(jpsi_binw_fit[pt_bin_number])

        pt_bin_number += 1




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

    ax_barrel.set_xlabel(r"$p_{T}$ [GeV]", fontsize=14)
    ax_endcap.set_xlabel(r"$p_{T}$ [GeV]", fontsize=14)
    ax_barrel.set_ylabel("Efficiency", fontsize=15)
    ax_barrel.set_title(r"$|\eta | < 1.4442$", fontsize=15)
    ax_endcap.set_title(r"$1.566 < |\eta | < 2.50$", fontsize=15)
    bottom_ax_barrel.set_xlabel(r"$p_{T}$ [GeV]", fontsize=14)
    bottom_ax_endcap.set_xlabel(r"$p_{T}$ [GeV]", fontsize=14)
    bottom_ax_barrel.set_ylabel("Data/MC", fontsize=15)


    # FIXME: EZT RENDESEN, HOGY A SF ERRORBARBA IS SZAMITSON!!!!?

    jpsi_data_eff_barrel_error_total = cutoff_errorbar_at_one(jpsi_data_eff_barrel_error_total, jpsi_data_eff_barrel)
    jpsi_mc_eff_barrel_error_total = cutoff_errorbar_at_one(jpsi_mc_eff_barrel_error_total, jpsi_mc_eff_barrel)
    z_data_eff_barrel_error_total = cutoff_errorbar_at_one(z_data_eff_barrel_error_total, z_data_eff_barrel)
    z_mc_eff_barrel_error_total = cutoff_errorbar_at_one(z_mc_eff_barrel_error_total, z_mc_eff_barrel)
    jpsi_data_eff_endcap_error_total = cutoff_errorbar_at_one(jpsi_data_eff_endcap_error_total, jpsi_data_eff_endcap)
    jpsi_mc_eff_endcap_error_total = cutoff_errorbar_at_one(jpsi_mc_eff_endcap_error_total, jpsi_mc_eff_endcap)
    z_data_eff_endcap_error_total = cutoff_errorbar_at_one(z_data_eff_endcap_error_total, z_data_eff_endcap)
    z_mc_eff_endcap_error_total = cutoff_errorbar_at_one(z_mc_eff_endcap_error_total, z_mc_eff_endcap)


    


    #ax_barrel.plot(ptbins_center_fit, data_eff_barrel, "bs")
    ax_barrel.errorbar(jpsi_ptbins_center_fit, jpsi_data_eff_barrel, yerr=jpsi_data_eff_barrel_error_total, xerr=jpsi_binw_fit, capsize=3, fmt=".", color="blue", ecolor = "blue", label=fname1+" Data")
    #ax_barrel.plot(ptbins_center_fit, mc_eff_barrel, "rs")
    ax_barrel.errorbar(jpsi_ptbins_center_fit, jpsi_mc_eff_barrel, yerr=jpsi_mc_eff_barrel_error_total, xerr=jpsi_binw_fit, capsize=3, fmt=".", color="crimson", ecolor = "crimson", label=fname1+" MC")


    ax_barrel.errorbar(z_ptbins_center_fit, z_data_eff_barrel, yerr=z_data_eff_barrel_error_total, xerr=z_binw_fit, capsize=3, fmt=".", color="cyan", ecolor = "cyan", label=fname2+" Data")
    ax_barrel.errorbar(z_ptbins_center_fit, z_mc_eff_barrel, yerr=z_mc_eff_barrel_error_total, xerr=z_binw_fit, capsize=3, fmt=".", color="magenta", ecolor = "magenta", label=fname2+" MC")


    bottom_ax_barrel.errorbar(jpsi_scale_factor_x_barrel, jpsi_scale_factor_barrel, yerr=jpsi_scale_factor_barrel_error, xerr=jpsi_scale_factor_x_barrel_binw, capsize=3, fmt=".", color="black", ecolor = "black", label = fname1 + " SF")
    bottom_ax_barrel.errorbar(jpsi_scale_factor_barrel_black_x, jpsi_scale_factor_barrel_black, xerr=jpsi_scale_factor_barrel_black_x_binw, fmt="x", color="black", ecolor="black")


    bottom_ax_barrel.errorbar(z_scale_factor_x_barrel, z_scale_factor_barrel, yerr=z_scale_factor_barrel_error, xerr=z_scale_factor_x_barrel_binw, capsize=3, fmt=".", color="chocolate", ecolor = "chocolate", label = fname2 + " SF")
    bottom_ax_barrel.errorbar(z_scale_factor_barrel_black_x, z_scale_factor_barrel_black, xerr=z_scale_factor_barrel_black_x_binw, fmt="x", color="chocolate", ecolor="chocolate")

    #ax_endcap.plot(ptbins_center_fit, data_eff_endcap, "bs")
    ax_endcap.errorbar(jpsi_ptbins_center_fit, jpsi_data_eff_endcap, yerr=jpsi_data_eff_endcap_error_total, xerr=jpsi_binw_fit, capsize=3, fmt=".", color="blue", ecolor = "blue", label=fname1+" Data")
    #ax_endcap.plot(ptbins_center_fit, mc_eff_endcap, "rs")
    ax_endcap.errorbar(jpsi_ptbins_center_fit, jpsi_mc_eff_endcap, yerr=jpsi_mc_eff_endcap_error_total, xerr=jpsi_binw_fit, capsize=3, fmt=".", color="crimson", ecolor = "crimson", label=fname1+" MC")

    ax_endcap.errorbar(z_ptbins_center_fit, z_data_eff_endcap, yerr=z_data_eff_endcap_error_total, xerr=z_binw_fit, capsize=3, fmt=".", color="cyan", ecolor = "cyan", label=fname2+" Data")
    ax_endcap.errorbar(z_ptbins_center_fit, z_mc_eff_endcap, yerr=z_mc_eff_endcap_error_total, xerr=z_binw_fit, capsize=3, fmt=".", color="magenta", ecolor = "magenta", label=fname2+" MC")


    bottom_ax_endcap.errorbar(jpsi_scale_factor_x_endcap, jpsi_scale_factor_endcap, yerr=jpsi_scale_factor_endcap_error, xerr=jpsi_scale_factor_x_endcap_binw, capsize=3, fmt=".", color="black", ecolor = "black")
    bottom_ax_endcap.errorbar(jpsi_scale_factor_endcap_black_x, jpsi_scale_factor_endcap_black, xerr=jpsi_scale_factor_endcap_black_x_binw, fmt="x", color="black", ecolor="black")


    bottom_ax_endcap.errorbar(z_scale_factor_x_endcap, z_scale_factor_endcap, yerr=z_scale_factor_endcap_error, xerr=z_scale_factor_x_endcap_binw, capsize=3, fmt=".", color="chocolate", ecolor = "chocolate")
    bottom_ax_endcap.errorbar(z_scale_factor_endcap_black_x, z_scale_factor_endcap_black, xerr=z_scale_factor_endcap_black_x_binw, fmt="x", color="chocolate", ecolor="chocolate")

    ax_barrel.legend(loc="best", fontsize=17)
    bottom_ax_barrel.legend(loc="best", fontsize=17)

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

