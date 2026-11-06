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
import matplotlib
matplotlib.use('Agg') # non-interactive backend
import matplotlib.pyplot as plt
#plt.rcParams["text.usetex"] = True
from matplotlib.gridspec import GridSpec
from matplotlib.offsetbox import TextArea, VPacker, HPacker, AnchoredOffsetbox
from matplotlib.colors import LinearSegmentedColormap
import yaml
import ctypes
import scipy.odr as odr

from iminuit import Minuit


from fitting_functions import do_minuit_fit

from fitting_functions import do_template_fitBNLL as do_template_fit

from futils import qs, overlapIndicesOfLists, trapcalc, useBigger, useSmaller



def get_parser():
    import argparse
    argParser = argparse.ArgumentParser(description = "Argument parser")
    argParser.add_argument('config_file',help='config_file')
    argParser.add_argument('--threads', nargs='?', const=1, type=int)
    argParser.add_argument('--bins', nargs='?', const=1, type=str)
    argParser.add_argument('--effi', nargs='?', const="", type=str)
    return argParser






def get_minuit_data_cc(minuit_y, m):
    totsum = sum(minuit_y)
    return totsum - m.values['nbkg'], np.sqrt(totsum) |qs| m.errors['nbkg']




def plot_two_with_text(left_plot_func, right_plot_func, text_info, doBkg, output_dir, bin_number, config_file, minuit_infotext = ""):
    """
    Create a canvas with:
    - Left plot
    - Right plot
    - Small text info area on the right

    Args:
        left_plot_func: function that takes an Axes object and plots on it (left plot)
        right_plot_func: function that takes an Axes object and plots on it (right plot)
        text_info: string to display in the text area

    """
    
    if minuit_infotext:
        fig, axes = plt.subplots(1, 4, figsize=(15, 5), gridspec_kw={'width_ratios': [1, 8, 8, 2]}) #, constrained_layout=True)
        ax_minuit, ax_left, ax_right, ax_text = axes

        ax_minuit.axis('off')
        ax_minuit.text(0.5, 0.5, minuit_infotext, ha='left', va='center', wrap=True, fontsize=8)
    else:
        # Create figure with custom width ratios: left 4, right 4, text 1
        fig, axes = plt.subplots(1, 3, figsize=(15, 5), gridspec_kw={'width_ratios': [4, 4, 1]}) #, constrained_layout=True)
    
        ax_left, ax_right, ax_text = axes
    
    # Hide axes for text area
    ax_text.axis('off')
    
    # Plot on left if function provided
    if left_plot_func is not None:
        left_plot_func(ax_left)
    
    # Plot on right if function provided
    if right_plot_func is not None:
        right_plot_func(ax_right)
    
    # Add text to text area
    ax_text.text(0.5, 0.5, text_info, ha='center', va='center', wrap=True, fontsize=16)
    
    for ax in [ax_left, ax_right]:
        ax.set_xlabel("mass [GeV]") #, fontsize=18, labelpad=3)
        ax.set_ylabel("Events") #, fontsize=18)

        ax.xaxis.set_label_coords(1.0, -0.08)  
        ax.xaxis.label.set_horizontalalignment('right')

    plt.tight_layout()

    if doBkg:
        if not os.path.exists(output_dir+"/fits/chev"+str(config_file["chevorder"])):
            os.makedirs(output_dir+"/fits/chev"+str(config_file["chevorder"]))
        plt.savefig(output_dir+"/fits/chev"+str(config_file["chevorder"])+"/bin"+str(bin_number)+".png")
    else:
        if not os.path.exists(output_dir+"/fits/"):
            os.makedirs(output_dir+"/fits/")
        plt.savefig(output_dir+"/fits/bin"+str(bin_number)+".png")
    plt.close()




def plot_two_with_pull_and_text(left_plot_func, right_plot_func, text_info, doBkg, output_dir, bin_number, config_file, rangeplot="", bTemplate=False, minuit_infotext = ""):
    """
    Create a canvas with:
    - Left plot
    - Right plot
    - Small text info area on the right

    Args:
        left_plot_func: function that takes an Axes object and plots on it (left plot)
        right_plot_func: function that takes an Axes object and plots on it (right plot)
        text_info: string to display in the text area

    """
    


    if minuit_infotext:

        fig = plt.figure(figsize=(21, 8)) #, constrained_layout=True)
        gs = GridSpec(
            2, 4,
            width_ratios=[1, 4, 4, 1],
            height_ratios=[4, 1],
            #hspace=0.05,
            #wspace=0.25,
        )

        ax_minuit = fig.add_subplot(gs[:, 0])  # spans both rows

        ax_left = fig.add_subplot(gs[0, 1])
        ax_left_pull = fig.add_subplot(gs[1, 1], sharex=ax_left)

        ax_right = fig.add_subplot(gs[0, 2])
        ax_right_pull = fig.add_subplot(gs[1, 2], sharex=ax_right)

        ax_text = fig.add_subplot(gs[:, 3])   # spans both rows

        ax_minuit.axis('off')
        ax_minuit.text(0.001, 0.999, minuit_infotext, ha='left', va='top', wrap=True, fontsize=11)





        if isinstance(text_info, str):
            ax_text.text(0.5, 0.5, text_info, ha='center', va='center', wrap=True, fontsize=16)
        else:


            lines = []
            for line_tup in text_info:
                if isinstance(line_tup, list):
                    children = []
                    for hline_tup in line_tup:
                        children.append( TextArea(hline_tup[0], textprops=hline_tup[1])   )

                elif isinstance(line_tup, tuple):
                    children = [TextArea(line_tup[0], textprops=line_tup[1])]
                elif isinstance(line_tup, str):
                    children = [TextArea(line_tup, textprops={"color": "black", "fontsize": 16})]
                    

                line = HPacker(
                    children=children,
                    align="baseline",
                    sep=0,
                )
                lines.append(line)

            text_box = VPacker(
                children=lines,
                align="center",
                sep=0,
            )

            anchored_text = AnchoredOffsetbox(
                loc="upper left",
                child=text_box,
                frameon=False,
                bbox_to_anchor=(0.05, 0.95),
                bbox_transform=ax_text.transAxes,
            )

            ax_text.add_artist(anchored_text)


    else:
        fig = plt.figure(figsize=(15, 7))
        gs = GridSpec(
            2, 3,
            width_ratios=[4, 4, 1],
            height_ratios=[4, 1],
            #hspace=0.05,
            #wspace=0.25,
        )

        ax_left = fig.add_subplot(gs[0, 0])
        ax_left_pull = fig.add_subplot(gs[1, 0], sharex=ax_left)

        ax_right = fig.add_subplot(gs[0, 1])
        ax_right_pull = fig.add_subplot(gs[1, 1], sharex=ax_right)

        ax_text = fig.add_subplot(gs[:, 2])   # spans both rows
        # Add text to text area
        ax_text.text(0.5, 0.5, text_info, ha='center', va='center', wrap=True, fontsize=16)

    







    # Hide axes for text area
    ax_text.axis("off")

    
    # Plot on left if function provided
    if left_plot_func is not None:
        left_plot_func(ax_left, ax_left_pull)
    
    # Plot on right if function provided
    if right_plot_func is not None:
        right_plot_func(ax_right, ax_right_pull)
    
    
    for ax in [ax_left, ax_right]:
        if not ax.get_xlabel():
            ax.set_xlabel("mass [GeV]") #, fontsize=18, labelpad=3)
        if not ax.get_ylabel():
            ax.set_ylabel("Events") #, fontsize=18)

        ax.xaxis.set_label_coords(1.0, -0.08)  
        ax.xaxis.label.set_horizontalalignment('right')

        #plt.setp(ax.get_xticklabels(), visible=False)


    for ax_pull in [ax_right_pull, ax_left_pull]:
        ax_pull.set_ylim(-3, 3)
        ax_pull.set_ylabel("Pull (Baker-Cousins)")
        ax_pull.grid(ls="--", alpha=0.4)

    #fig.suptitle("")

    plt.tight_layout()



    if rangeplot:
        rangeplot = "/"+rangeplot



    fitdir = "/fits/chev"+str(config_file["chevorder"])+rangeplot if not bTemplate else "/templatefit"+rangeplot


    if doBkg:
        if not os.path.exists(output_dir+fitdir):
            os.makedirs(output_dir+fitdir)
        plt.savefig(output_dir+fitdir+"/bin"+str(bin_number)+".png")
    else:
        if not os.path.exists(output_dir+"/fits/"):
            os.makedirs(output_dir+"/fits/")
        plt.savefig(output_dir+"/fits/bin"+str(bin_number)+".png")
    plt.close()








def minuit_data_plot(config_file, workdir, binnames, data_pass_plotinfo, data_fail_plotinfo, pass_index, fail_index, minuit_efficiency, minuit_efficiency_error, minuit_efficiency_range_error, rangeplot=""):
    # Minuit plot


    which_range = "_"+rangeplot if rangeplot else ""    # "", "_down", or "_up"


    def base_plot(ax, ax_pull, plotinfo, index, sProbeType):

        if config_file['bUseMinuit']:
            hist_integral = plotinfo["minuit_histintegral"][0]
            hist_integral_error = plotinfo["minuit_histintegral"][1]
            hist_integral_rangesyst_error = np.std([hist_integral, plotinfo["minuit_down_histintegral"][0], plotinfo["minuit_up_histintegral"][0] ])

            rawint_text = "Data, integral - bkg = \n{:.1f}".format(hist_integral)+u"\u00B1"+("{:.1f}".format(100*hist_integral_error / hist_integral) if hist_integral else "?")+"%"+r"$^{\mathrm{stat}}$"

            if not which_range:
                rawint_text = rawint_text + u"\u00B1" +("{:.1f}".format(100*hist_integral_rangesyst_error / hist_integral) if hist_integral else "?")+"%"+r"$^{\mathrm{range}}$"

        else:
            hist_integral = plotinfo["hist_integral"][0]

            rawint_text = "BKsub integral/bw = \n{:.1f}".format(hist_integral)+u"\u00B1"+("{:.1f}".format(100*plotinfo["Delta_hist_integral"] / hist_integral) if hist_integral else "?")+"%"+r"$^{\mathrm{stat}}$"

            if plotinfo["linkb_hist_integral_Delta"]:
                rawint_text = rawint_text + u"\u00B1" +("{:.1f}".format(100*plotinfo["linkb_hist_integral_Delta"] / hist_integral) if hist_integral else "?")+"%"+r"$^{\mathrm{bk}}$"

        m = plotinfo['minuit_info'+which_range]["m"]
        xcent = np.array(plotinfo["minuit_x"+which_range])
        y = np.array(plotinfo["minuit_y"+which_range])
        bin_width = xcent[1]-xcent[0]
        edges = np.concatenate(([xcent[0]-bin_width/2], xcent+bin_width/2))

        x = plotinfo['minuit_info'+which_range]["x_sample"]

        x_axis = plotinfo["x_axis"]
        y_axis = plotinfo["y_axis"]
        x_axis_black = plotinfo["x_axis_black"]
        y_axis_black = plotinfo["y_axis_black"]
        x_axis_pm = plotinfo["x_axis_pm"]
        y_axis_pm = plotinfo["y_axis_pm"]
        ax.plot(x_axis,y_axis, "sb", label = rawint_text)
        ax.errorbar(x_axis, y_axis, yerr=y_axis_pm, xerr = x_axis_pm, capsize=3, fmt=".", ecolor = "black")
        ax.plot(x_axis_black, y_axis_black, "x", color="black")


        goodlim_y = ax.get_ylim()

        if index in config_file['data_override_plots_xrange']:
            min_y = 9999
            max_y = -1
            for i in range(len(x_axis)):
                if config_file['data_override_plots_xrange'][index][0] <= x_axis[i] <= config_file['data_override_plots_xrange'][index][-1]:
                    if y_axis[i] < min_y:
                        min_y = y_axis[i]
                    if y_axis[i] > max_y:
                        max_y = y_axis[i]

            goodlim_y = [min_y*0.9, max_y*1.1]


        signal = plotinfo['minuit_info'+which_range]["signal"]
        background = plotinfo['minuit_info'+which_range]["background"]
        total = signal + background

        pull = plotinfo['minuit_info'+which_range]["pull"]

        ax_pull.axhline(0, color='black', lw=1)
        ax_pull.axhline(1, color='gray', ls='--', lw=1)
        ax_pull.axhline(-1, color='gray', ls='--', lw=1)

        ax_pull.errorbar(
            xcent,
            pull,
            fmt='o',
            color='black',
            markersize=4,
        )





        #ax.plot(
        #    x,
        #    signal * bin_width,
        #    color="green",
        #    ls=":",
        #    lw=2,
        #    label="Signal",
        #)

        ax.plot(
            x,
            background * bin_width,
            color="blue",
            ls="--",
            lw=2,
            label=str(config_file["chevorder"])+"th order chev",
        )

        ax.plot(
            x,
            total * bin_width,
            color="red",
            lw=2,
            label="Total fit",
        )

        # 5. Uncertainty Band (from Sampling)
        if plotinfo['minuit_info'+which_range]["uncert_band"] is not None:
            ax.fill_between(x, plotinfo['minuit_info'+which_range]["uncert_band"][0], plotinfo['minuit_info'+which_range]["uncert_band"][1],
                                color="hotpink", alpha=0.3) #, label='Fit Uncertainty')


        if index in config_file['data_override_plots_xrange']:
            ax.set_xlim(config_file['data_override_plots_xrange'][index])
            ax.set_xticks(np.arange(config_file['data_override_plots_xrange'][index][0], config_file['data_override_plots_xrange'][index][1], config_file['xticks_density']), minor=True)
        else:
            ax.set_xlim(config_file['plots_xrange'])
            ax.set_xticks(np.arange(config_file['plots_xrange'][0], config_file['plots_xrange'][1], config_file['xticks_density']), minor=True)


        #ax.set_ylim(goodlim_y)


        fulltext = r"Signal (CB)" + ", integral/bw = "+"\n{:.1f}".format(m.values["nsig"])+u"\u00B1"+("{:.1f}".format(100*m.errors["nsig"] / m.values["nsig"]) if m.values["nsig"] else "?")+"%"+r"$^{\mathrm{st}}$"


        if not rangeplot:
            # nominal plot gets range syst
            integral_down = plotinfo["minuit_info_down"]["m"].values["nsig"]
            integral_up = plotinfo["minuit_info_up"]["m"].values["nsig"]

            #greater_diff = m.values["nsig"] - integral_down if abs(m.values["nsig"] - integral_down) > abs(m.values["nsig"] - integral_up) else integral_up
            #integral_rangesyst_std = abs(m.values["nsig"] - greater_diff)

            integral_rangesyst_std = np.std([m.values["nsig"],integral_down,integral_up])

            fulltext += u"\u00B1"+("{:.1f}".format(100*integral_rangesyst_std / m.values["nsig"]) if m.values["nsig"] else "?")+"%"+r"$^{\mathrm{range}}$"


        ax.stairs(
            plotinfo["minuit_info"+which_range]["signal_binned"],
            plotinfo["minuit_info"+which_range]["bin_edges"],
            #where="mid",
            color="teal",
            linestyle="--",
            label=fulltext,
            zorder=1,
            linewidth=2,
        )


        display_binname = binnames[index][:4] + ", DATA " + sProbeType + " PROBES\n"
        display_binname += binnames[index][6:-4] + " WP: " + config_file['workingpoint']
        ax.set_title(display_binname, fontsize=11)

        leg = ax.legend(loc='best', fontsize=10)
        leg.get_title().set_fontsize(17)
        leg.get_frame().set_alpha(0.42)

        ax.grid(linestyle = "--")

        #vlinecolor = 'gray'
        #ax.vlines(master_integral_range[0], goodlim_y[0],goodlim_y[1], color=vlinecolor)
        #ax.vlines(master_integral_range[1], goodlim_y[0],goodlim_y[1], color=vlinecolor)


    def pass_plot(ax, ax_pull):
        return base_plot(ax, ax_pull, data_pass_plotinfo, pass_index, "PASSING")

    def fail_plot(ax, ax_pull):
        return base_plot(ax, ax_pull, data_fail_plotinfo, fail_index, "FAILING")

    # build text info:
    text_info = config_file["resonance"] + " DATA\n\n"
    #text_info += str(binnames[pass_index][:-4])+"\n"

    #text_info += "C&C -BK in\n[%g,%g] GeV:\n%s=%.3f\n" % (master_integral_range[0], master_integral_range[1], r'$\eta$', data_efficiency)
    #text_info += u"\u00B1"+("{:.1f}".format(100*data_efficiency_stat_error / data_efficiency) if data_efficiency else "?")+"%"+r"$^{\mathrm{stat}}$"+"\n"
    #text_info += u"\u00B1"+("{:.1f}".format(100*data_efficiency_rangesyst_error / data_efficiency) if data_efficiency else "?")+"%"+r"$^{\mathrm{range}}$"+"\n\n"



    text_info += "CB+chev3 fit:\n%s=%.1f%%\n" % (r'$\eta$', minuit_efficiency*100)
    text_info += u"\u00B1"+("{:.1f}".format(100*minuit_efficiency_error) if minuit_efficiency else "?")+"%"+r"$^{\mathrm{stat,(abs)}}$"+"\n"

    if not rangeplot:
        text_info += u"\u00B1"+("{:.1f}".format(100*minuit_efficiency_range_error) if minuit_efficiency else "?")+"%"+r"$^{\mathrm{range,(abs)}}$"+"\n"


    plot_two_with_pull_and_text(pass_plot, fail_plot, text_info, True, workdir+"dataplots/", str(binnames[pass_index][3]), config_file, rangeplot=rangeplot)




def minuit_template_plot(config_file, workdir, binnames, data_pass_plotinfo, data_fail_plotinfo, pass_index, fail_index, minuit_efficiency_tuple, mc_efficiency_tuple, rangeplot="", minuit_range_efficiencies = (), mc_range_efficiencies = ()):
    # Minuit plot


    which_range = "_"+rangeplot if rangeplot else ""


    def base_plot(ax, ax_pull, plotinfo, index, sProbeType):

        hist_integral = plotinfo["cc"][0]
        hist_integral_error = plotinfo["cc"][1]

        rawint_text = "Data, integral - bkg = \n{:.1f}".format(hist_integral)+u"\u00B1"+("{:.1f}".format(100*hist_integral_error / hist_integral) if hist_integral else "?")+"%"+r"$^{\mathrm{stat}}$"

        if not which_range:
            hist_integral_rangesyst_error = np.std([hist_integral, plotinfo["cc_down"][0], plotinfo["cc_up"][0] ])

            rawint_text = rawint_text + u"\u00B1" +("{:.1f}".format(100*hist_integral_rangesyst_error / hist_integral) if hist_integral else "?")+"%"+r"$^{\mathrm{range}}$"





        m = plotinfo["m"]

        x = plotinfo["x"]
        data = plotinfo["data"]
        mc = plotinfo["mc"]
        expected = plotinfo["expected"]
        signal = plotinfo["signal"]
        background = plotinfo["background"]
        pull = plotinfo["pull"]




        mctext = r"MC" + ", integral/bw = "+"\n{:.1f}".format(plotinfo["mc_hist_integral"][0])+u"\u00B1"+("{:.1f}".format(100*plotinfo["mc_hist_integral"][1] / plotinfo["mc_hist_integral"][0]) if plotinfo["mc_hist_integral"][0] else "?")+"%"+r"$^{\mathrm{st}}$"


        ax_mc = ax.twinx()
        ax_mc.set_ylabel("Events (MC)", color="darkorange")
        ax_mc.tick_params(axis="y", labelcolor="darkorange")
        mc_handle, = ax_mc.step(
            x,
            mc,
            where="mid",
            color="darkorange",
            label=mctext,
            zorder=1
        )


        data_handle = ax.errorbar(
            x,
            data,
            yerr=np.sqrt(data),
            fmt="o",
            color="black",
            label=rawint_text,
            zorder=2
        )


        ax_mc.set_zorder(0)
        ax.set_zorder(1)
        ax.patch.set_visible(False)

        fulltext = r"Fit signal" + ", integral/bw = "+"\n{:.1f}".format(m.values["x0"])+u"\u00B1"+("{:.1f}".format(100*m.errors["x0"] / m.values["x0"]) if m.values["x0"] else "?")+"%"+r"$^{\mathrm{st}}$"

        if not rangeplot:
            # nominal plot gets range syst
            integral_down = plotinfo["m_down"].values["x0"]
            integral_up = plotinfo["m_up"].values["x0"]

            #greater_diff = m.values["x0"] - integral_down if abs(m.values["x0"] - integral_down) > abs(m.values["x0"] - integral_up) else integral_up
            #integral_rangesyst_std = abs(m.values["x0"] - greater_diff)

            integral_rangesyst_std = np.std([m.values["x0"],integral_down,integral_up])

            fulltext += u"\u00B1"+("{:.1f}".format(100*integral_rangesyst_std / m.values["x0"]) if m.values["x0"] else "?")+"%"+r"$^{\mathrm{range}}$"

        fit_handle, = ax.step(
            x,
            expected,
            where="mid",
            color="red",
            label="Fit total",
            zorder=3
        )


        #if len(plotinfo["mc_hist_integral"]) > 2 and not rangeplot:
        #    mctext += u"\u00B1"+("{:.1f}".format(100*plotinfo["mc_hist_integral"][2] / plotinfo["mc_hist_integral"][0]) if plotinfo["mc_hist_integral"][0] else "?")+"%"+r"$^{\mathrm{range}}$"

        signal_handle, = ax.step(
            x,
            signal,
            where="mid",
            color="teal",
            linestyle="--",
            label=fulltext,
            zorder=1
        )

        bk_handle, = ax.step(
            x,
            background,
            where="mid",
            color="blue",
            linestyle="--",
            label="Fit bkg (CMSshape)",
            zorder=1
        )

        ax_pull.errorbar(
            x,
            pull,
            fmt="o",
            color="black",
            markersize=4,
        )

        ax_pull.axhline(
            0,
            color="black",
            linewidth=1,
        )

        ax_pull.axhline(
            1,
            color="gray",
            linestyle="--",
        )

        ax_pull.axhline(
            -1,
            color="gray",
            linestyle="--",
        )

        ax_pull.axhline(
            2,
            color="gray",
            linestyle=":")

        ax_pull.axhline(
            -2,
            color="gray",
            linestyle=":"
        )

        ax_pull.set_ylabel("BC pull")
        ax_pull.set_xlabel("pair mass [GeV]")
        ax_pull.set_ylim(-5, 5)



        ax.set_xlabel("pair mass [GeV]")
        ax.set_ylabel("Events (Data)")


        '''
        goodlim_y = ax.get_ylim()

        if index in config_file['data_override_plots_xrange']:
            min_y = 9999
            max_y = -1
            for i in range(len(x)):
                if config_file['data_override_plots_xrange'][index][0] <= x[i] <= config_file['data_override_plots_xrange'][index][-1]:
                    if y_axis[i] < min_y:
                        min_y = y_axis[i]
                    if y_axis[i] > max_y:
                        max_y = y_axis[i]

            goodlim_y = [min_y*0.9, max_y*1.1]
        '''




        # FIXME FIXME FIXME: How to draw uncertainty band?
        '''
        # 5. Uncertainty Band (Sampling)
        pcov = m.covariance
        if pcov is not None:
            n_samples = 100
            vopts = np.random.multivariate_normal(m.values, pcov, n_samples)
            sampled_ydata = []
            for v in vopts:
                if v[0] <= 0: continue
                if v[1] <= 0: continue
                if v[3] <= 0: continue
                signal_i = v[0] * signal_pdf(
                    x,
                    v[2],
                    v[3],
                    v[4],
                    v[5],
                )
                background_i = v[1] * background_pdf(
                    x,
                    v[6],
                    v[7],
                    v[8],
                )
                total_i = signal_i + background_i
                y_sample = total_i * bin_width
                sampled_ydata.append(y_sample)
            sampled_ydata = np.vstack(sampled_ydata)
            model_uncert = np.nanstd(sampled_ydata, axis=0)
            ax.fill_between(x, y_fit_curve - model_uncert, y_fit_curve + model_uncert,
                                color="hotpink", alpha=0.3) #, label='Fit Uncertainty')
        '''

        '''
        if index in config_file['data_override_plots_xrange']:
            ax.set_xlim(config_file['data_override_plots_xrange'][index])
            ax.set_xticks(np.arange(config_file['data_override_plots_xrange'][index][0], config_file['data_override_plots_xrange'][index][1], config_file['xticks_density']), minor=True)
        else:
            ax.set_xlim(config_file['plots_xrange'])
            ax.set_xticks(np.arange(config_file['plots_xrange'][0], config_file['plots_xrange'][1], config_file['xticks_density']), minor=True)


        ax.set_ylim(goodlim_y)
        '''

        display_binname = binnames[index][:4] + ", DATA " + sProbeType + " PROBES\n"
        display_binname += binnames[index][6:-4] + " WP: " + config_file['workingpoint']
        ax.set_title(display_binname, fontsize=11)

        leg = ax.legend(loc='best', fontsize=10, handles = [
            data_handle,
            mc_handle,
            fit_handle,
            signal_handle,
            bk_handle
        ])
        leg.get_title().set_fontsize(17)
        leg.get_frame().set_alpha(0.42)

        '''
        ax.grid(linestyle = "--")
        # Get the left-axis limits
        ymin, ymax = ax.get_ylim()
        # Get the right-axis limits
        mc_ymin, mc_ymax = ax_mc.get_ylim()
        # Choose the left-axis tick positions
        left_ticks = ax.get_yticks()
        # Convert their positions to the corresponding right-axis values
        right_ticks = (
            mc_ymin
            + (left_ticks - ymin) / (ymax - ymin) * (mc_ymax - mc_ymin)
        )
        ax_mc.set_yticks(right_ticks)
        '''

    def pass_plot(ax, ax_pull):
        return base_plot(ax, ax_pull, data_pass_plotinfo, pass_index, "PASSING")

    def fail_plot(ax, ax_pull):
        return base_plot(ax, ax_pull, data_fail_plotinfo, fail_index, "FAILING")


    # build text info:
    text_info = []
    text_info.append( config_file["collection"] + "\n" + config_file["workingpoint"] + "\n\n"  )
    text_info.append( config_file["resonance"] + " DATA")
    #text_info += str(binnames[pass_index][:-4])+"\n"

    #text_info += "C&C -BK in\n[%g,%g] GeV:\n%s=%.3f\n" % (master_integral_range[0], master_integral_range[1], r'$\eta$', data_efficiency)
    #text_info += u"\u00B1"+("{:.1f}".format(100*data_efficiency_stat_error / data_efficiency) if data_efficiency else "?")+"%"+r"$^{\mathrm{stat}}$"+"\n"
    #text_info += u"\u00B1"+("{:.1f}".format(100*data_efficiency_rangesyst_error / data_efficiency) if data_efficiency else "?")+"%"+r"$^{\mathrm{range}}$"+"\n\n"



    minuit_efficiency = minuit_efficiency_tuple[0]
    minuit_efficiency_error  = minuit_efficiency_tuple[1]


    line = "Template fit:\n%s=%.1f%%\n" % (r'$\eta$', minuit_efficiency*100)
    line += u"\u00B1"+("{:.1f}".format(100*minuit_efficiency_error) if minuit_efficiency else "?")+"%"+r"$^{\mathrm{stat,(abs)}}$"+"\n" 


    if not rangeplot:
        minuit_efficiency_range_error = minuit_efficiency_tuple[2]
        line +=   u"\u00B1"+("{:.1f}".format(100*minuit_efficiency_range_error) if minuit_efficiency else "?")+"%"+r"$^{\mathrm{range,(abs)}}$"+"\n" 

    text_info.append( ( line , {"color": "teal", "fontsize": 16}   ) )


    text_info.append( "\n\n")

    text_info.append( ( "MC CC" , {"color": "darkorange", "fontsize": 16}   ) )



    mc_efficiency = mc_efficiency_tuple[0]
    mc_efficiency_stat_error = mc_efficiency_tuple[1]


    text_info.append( ( "%s=%.1f%%" % (r'$\eta$', mc_efficiency*100) , {"color": "darkorange", "fontsize": 16}   ) )
    text_info.append( ( u"\u00B1"+("{:.1f}".format(100*mc_efficiency_stat_error) if mc_efficiency else "?")+"%"+r"$^{\mathrm{stat,(abs)}}$" , {"color": "darkorange", "fontsize": 16}   ) )

    if not rangeplot:
        mc_efficiency_rangesyst_error = mc_efficiency_tuple[2]
        text_info.append( ( u"\u00B1"+("{:.1f}".format(100*mc_efficiency_rangesyst_error) if mc_efficiency else "?")+"%"+r"$^{\mathrm{range,(abs)}}$" , {"color": "darkorange", "fontsize": 16}   ) )


    
    sf = minuit_efficiency / mc_efficiency
    sf_error_stat = getSFerrorbar(minuit_efficiency, mc_efficiency, minuit_efficiency_error, mc_efficiency_stat_error)



    text_info.append("\n\n")
    text_info.append("SF=%.2f" % sf)
    text_info.append(u"\u00B1"+("{:.1f}".format(100*sf_error_stat/sf) if sf else "?")+"%"+r"$^{\mathrm{stat}}$")

    if not rangeplot:
        sf_error_range = np.std( [sf,  minuit_range_efficiencies[0]/mc_range_efficiencies[0], minuit_range_efficiencies[1]/mc_range_efficiencies[1] ])
        text_info.append(u"\u00B1"+("{:.1f}".format(100*sf_error_range/sf) if sf else "?")+"%"+r"$^{\mathrm{range}}$")


    minuit_infotext = ""
    if config_file["bDebugReportParams"]:
        minuit_infotext = str(data_pass_plotinfo["m"]) + "\n\n" + str(data_fail_plotinfo["m"])


        def get_minuit_infotext(plotinfo):

            m = plotinfo["m"]
            lines = []
            lines.append("Migrad (Poisson NLL):")
            lines.append(f"FCN = {m.fval:.4g}")
            lines.append(f"EDM = {m.fmin.edm:.3g}")
            lines.append(f"Nfcn = {m.nfcn}")
            #lines.append("──────────")
            lines.append(f"ndof = {plotinfo["ndof"]}")
            lines.append(f"Pearson chi2/ndof = {plotinfo["chi2"]/plotinfo["ndof"]:.3g}")
            lines.append(f"Poisson devi/ndof = {plotinfo["deviance"]/plotinfo["ndof"]:.3g}")
            #lines.append(f"p value = {plotinfo["p_value"]:.3g}")
            #lines.append("──────────")
            lines.append("✓ Valid Minimum" if m.fmin.is_valid else "✗ INVALID minimum")

            lines.append("⚠ ABOVE EDM threshold (goal x 10)" if m.fmin.is_above_max_edm else "✓ Below EDM threshold (goal x 10)")

            lines.append("⚠ SOME parameters at limit" if m.fmin.has_parameters_at_limit else "✓ No parameters at limit" )

            lines.append("⚠ Reached call limit" if m.fmin.has_reached_call_limit else "✓ Below call limit")

            lines.append("✗ Hesse FAILED" if m.fmin.hesse_failed else "✓ Hesse ok")

            if m.fmin.has_accurate_covar:
                cov_status = "✓ Covariance accurate"
            elif m.fmin.has_made_posdef_covar:
                cov_status = "⚠ Covariance FORCED pos. def."
            elif m.fmin.has_posdef_covar:
                cov_status = "✓ Covariance pos. def."
            else:
                cov_status = "✗ Covariance NOT pos. def."

            lines.append(cov_status)

            #lines.append("──────────")
            for name in m.parameters:
                value = m.values[name]
                error = m.errors[name]
                lower, upper = m.limits[name]

                at_limit = ""
                if lower is not None:
                    if abs(value - lower) < 1e-5*value:
                        at_limit = "⚠ "
                if upper is not None:
                    if abs(value - upper) < 1e-5*value:
                        at_limit = "⚠ "


                dims = {
                    "alpha": " GeV",
                    "beta": " 1/GeV",
                    "gamma": " 1/GeV"
                }

                dimension = ""
                if name in dims.keys():
                    dimension = dims[name]


                lines.append(f"{at_limit}{name:6s} = ({value:.3g} ± {error:.2g}){dimension}")


            lines.append(f"nsig/MCint = {m.values["x0"]/plotinfo["mc_hist_integral"][0]:.3g}")


            #cov_array = np.asarray(m.covariance)
            #corr_array = np.asarray(m.covariance.correlation())


            minuit_infotext = "\n".join(lines)
            return minuit_infotext


        pass_minuit_infotext = get_minuit_infotext(data_pass_plotinfo)
        fail_minuit_infotext = get_minuit_infotext(data_fail_plotinfo)

        minuit_infotext = "PASSING - "+pass_minuit_infotext+"\n──────────\nFAILING - "+fail_minuit_infotext


    plot_two_with_pull_and_text(pass_plot, fail_plot, text_info, True, workdir+"dataplots/", str(binnames[pass_index][3]), config_file, rangeplot=rangeplot, bTemplate=True, minuit_infotext = minuit_infotext)








'''
def getIntegralRangeSystError(integral, integrand, xrange_min, xrange_max, binwidth):
    res = useBigger( abs(integral - quad(integrand, xrange_min*0.95, xrange_max*1.05)[0]), abs(integral - quad(integrand, xrange_min*1.05, xrange_max*0.95)[0]) ) / binwidth
    return res
'''

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
            if y[i] > 0:
                errorbar_scale = abs(y_pm[i]/y[i])  * abs(y[i] - bkg) 

                ybksyst = []
                for k in range(len(parameters)):
                    ybksyst.append( Chev.chebval(x[i],parameters[k] ))
                #y_pm_sub.append( errorbar_scale |qs| np.std(ybksyst)  )    # FIXME FIXME: Make this a parameter?
                y_pm_sub.append( errorbar_scale   )
            else:
                y_pm_sub.append(1)

    if bGetError:
        return y_sub, y_pm_sub
    return y_sub

def get_bksub_integral(x, y, y_pm, chevpars, chev_std, int_range, fit_range, repeats, resonance):

    #print chevpars
    #print chev_std
    #stri = ""
    #for i in range(len(chev_std)):
    #    stri += str(abs(chev_std[i])/chevpars[i]*100)+"% "
    #print ("params: "+str(chevpars))
    #print ("Cheverrors: "+stri)


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

    integral_w_low = sum(y_sub_rangesyt_low)
    integral_w_high = sum(y_sub_rangesyt_high)


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


    ## BKSYST ##
    # calculate integral with expecting a linear background
    # a) try to find local minimum to the left of the first integral range
    locmin_i = 0
    iIsFirstBin = 0
    for i in range(len(x)):
        if x[i] < fit_range[0]: continue
        if x[i] >= int_range[0]: break

        if iIsFirstBin == 0:
            locmin_i = i
            iIsFirstBin += 1
            continue

        if y[i] < y[locmin_i]:
            locmin_i = i
            iIsFirstBin += 1

    if iIsFirstBin > 1 and resonance == "Z":    # FIXME also test for higher end?
        x1 = x[locmin_i]
        x2 = int_range[1]
        y1 = y[locmin_i]
        
        for i in range(len(x)):
            if x[i] >= int_range[1]:
                y2 = y[i]
                break

    else:
        # b) no local minimum, use previous version
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

    linbk_params = [[x1, x2], [y1, y2]]

    #print integrals_for_std
    #print integral_w, str(abs(np.sqrt(np.std(integrals_for_std)**2 + integral_w))/integral_w*100)+"%"
    return [integral_w, integral_w_low, integral_w_high], np.sqrt(np.std(integrals_for_std)**2 + integral_w), integral_linbk_w, linbk_params



def odr_dualfit_cb(x, y, x_pm, y_pm, chevpars, cbpars):

    p0 = []
    for cp in cbpars:
        p0.append(cp)
    for cp in chevpars:
        p0.append(cp)

    '''
    popt_sig, pcov_sig = curve_fit(crystalballPlusChev, x, y, p0, sigma=y_pm, absolute_sigma=True, maxfev=5000)
    perr_sig = np.sqrt(np.diag(pcov_sig))

    print( "success")

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


def curve_dualfit_cb(x, y, y_pm, p0 ):
    lower = []
    upper = []
    for par in p0:
        lower.append(par*0.95 if par > 0 else par*1.05)
        upper.append(par*1.05 if par > 1 else par*0.95)

    popt, pcov = curve_fit(crystalballPlusChev, x, y, p0=p0, sigma=y_pm, absolute_sigma=True, bounds=(lower, upper) )

    cov_signal = pcov[:5,:5]
    perr = np.sqrt(np.diag(pcov))


    cb_popt_sig = popt[:5]
    cb_perr_sig = perr[:5]
    chevparams = popt[5:]
    chev_std_devs = perr[5:] 
    return cb_popt_sig, cb_perr_sig, chevparams, chev_std_devs, cov_signal







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
            integral_w, Delta_integral = 1,1   #it's getting dual fit, skip now
        else:
            integral_w, Delta_integral = getIntegralErrorWithMCMethod(mycrystalball, popt_sig, perr_sig, binwidth, config_file, index, isData)

    elif function == "bwg":
        bounds = ((0,85,0.01,0.00001),(np.inf,95,50,20))  

        p0 = config_file['BWG_p0']
        p0[0] = 10000000 * (max(y) / 5000)

        if isData and index in config_file['data_override_bwg_p0']:
            p0 = config_file['data_override_bwg_p0'][index]

        popt_sig, pcov_sig = curve_fit(Voigt, x, y, p0, sigma=y_pm, absolute_sigma=True, bounds=bounds,maxfev=5000)
        perr_sig = np.sqrt(np.diag(pcov_sig))
        

        if config_file['bUseBWGauss'] and index in config_file[('MC_' if not isData else '') + 'dualfit_bounds_bwg']:
            integral_w, Delta_integral =  1,1  #it's getting dual fit, skip now
        else:
            integral_w, Delta_integral = getIntegralErrorWithMCMethod(Voigt, popt_sig, perr_sig, binwidth, config_file, index, isData)

    return popt_sig, perr_sig, integral_w, Delta_integral



def getIntegralErrorWithMCMethod(function, p0, p0err, binwidth, config_file, index, isData, bEstimateError=True):

    #xrange_min = config_file[('MC_' if not isData else '') + 'raw_integral_bounds'][index][0]
    #xrange_max = config_file[('MC_' if not isData else '') + 'raw_integral_bounds'][index][-1]

    master_integral_range = config_file["master_integral_range"]
    xrange_min = master_integral_range[0]
    xrange_max = master_integral_range[1]

    def integrand(x):
        return function([x], *p0)[0]
    integral = quad(integrand, xrange_min, xrange_max)
    integral_w = integral[0] / binwidth

    integral_w_low = (quad(integrand, xrange_min*0.95, xrange_max*0.95))[0] / binwidth
    integral_w_high = (quad(integrand, xrange_min*1.05, xrange_max*1.05))[0] / binwidth
    #integral_rangesyst_std = getIntegralRangeSystError(integral[0], integrand, xrange_min, xrange_max, binwidth)

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



    return [integral_w,integral_w_low, integral_w_high], Delta_integral




def is_valid_cb(p, p_nom):
    N, alpha, n, x0, sigma = p
    if N <= 0 or n <= 0 or sigma < 0.001: return False
    # caps to avoid overflows
    if N > 1e3 * p_nom[0]:   # e.g., 1000x nominal
        return False
    if sigma > 1.0:        
        return False
    if n > 1000 * p_nom[2]: # CB tail exponent too large
        return False
    return True


def getIntegralErrorWithMCMethodCov(function, p_nom, cov, binwidth, config_file, index):

    master_integral_range = config_file["master_integral_range"]
    xrange_min = master_integral_range[0]
    xrange_max = master_integral_range[1]

    def get_integral(p, xmin, xmax):
        def integrand(x):
            return function([x], *p)[0]
        integral = quad(integrand, xmin, xmax)
        integral_w = integral[0] / binwidth
        return integral_w


    eigvals, eigvecs = np.linalg.eigh(cov)
    eps = 1e-12  # tiny floor
    eigvals[eigvals < eps] = eps

    cov_fixed = np.dot(eigvecs, np.dot(np.diag(eigvals), eigvecs.T))
    cov_fixed = (cov_fixed + cov_fixed.T) / 2.0

    try:
        np.linalg.cholesky(cov_fixed)
    except np.linalg.LinAlgError:
        print("Data fit "+str(index)+": Covariance is NOT PSD: Cholesky failed, fallback")
        cov_fixed = np.diag(np.diag(cov_fixed))


    samples = np.random.multivariate_normal(p_nom, cov_fixed, size=config_file["integral_error_MCmethod_iterations"])
    mask = np.array([is_valid_cb(p, p_nom) for p in samples])
    samples = samples[mask]


    


    I_nominal = get_integral(p_nom, xrange_min, xrange_max)

    integrals_for_Delta = []
    integration_warnings = 0
    for p in samples:
        integration_warning = False
        with warnings.catch_warnings(record=True) as caught_warnings:
            warnings.simplefilter("always", IntegrationWarning)
            integral_i = get_integral(p, xrange_min, xrange_max)

            for warning in caught_warnings:
                if isinstance(warning.message, IntegrationWarning):
                    integration_warning = True
                    integration_warnings += 1

            if  not np.isnan(integral_i) and not np.isinf(integral_i) and not integration_warning:
                # ignore crazy values
                if 1 < integral_i < 3*I_nominal:
                    integrals_for_Delta.append( integral_i)

        
    if integration_warnings > 0:
        print ("Data Fit "+str(index)+", int delta est: "+str(integration_warnings)+"x  IntegrationWarning ")


    I_err = np.std(integrals_for_Delta, ddof=1)

    I_low = get_integral(p_nom, xrange_min*0.95, xrange_max*0.95)
    I_high = get_integral(p_nom, xrange_min*1.05, xrange_max*1.05)

    return [I_nominal, I_low, I_high], I_err






def getIntegralFuncErrorFromLinBK(chevparams, binwidth, config_file, linbk_params):

    # integral of underlying chev
    master_integral_range = config_file["master_integral_range"]
    xrange_min = master_integral_range[0]
    xrange_max = master_integral_range[1]

    def integrand(x):
        return Chev.chebval(x,chevparams)
    integral = quad(integrand, xrange_min, xrange_max)
    integral_w = integral[0] / binwidth

    #print("Integral under chev: "+str(integral_w))

    # integral with lin
    x1 = linbk_params[0][0]
    x2 = linbk_params[0][1]
    y1 = linbk_params[1][0]
    y2 = linbk_params[1][1]

    integral_linbk = 0.5 * (y1 + y2) * abs(x2 - x1)  / binwidth
    #print("Integral under lin: "+str(integral_linbk))


    return abs(integral_w - integral_linbk)



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





def getSampledXForPlotting(config_file, index, isData, doBkg, cb_settings):



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

    return x


# UNUSED!
def makePlot(output_dir, config_file, binname, index, isData, x_axis, x_axis_pm, y_axis, y_axis_pm, x_axis_black, y_axis_black, doBkg, chevx, chevparams, chevstats, hist_integral, Delta_hist_integral, hist_fullrangeintegral, Delta_hist_fullrangeintegral, linbksub_integral, linbk_params, integral_range_std, gauss_settings, cb_settings, psip_settings):
    
    plt.close()

    master_integral_range = config_file["master_integral_range"]

    chevorder = config_file["chevorder"]
    if index in config_file['chevorder_in_bin']:
        chevorder = config_file['chevorder_in_bin'][index]

    x = getSampledXForPlotting(config_file, index, isData, doBkg, cb_settings)

    fig = plt.figure(42,figsize=(15, 10))
    ax = fig.add_subplot(111)

    rawint_text = ("BKsub integral/bw = " if isData or doBkg else "Raw integral/bw = ")+"\n{:.1f}".format(hist_integral)+u"\u00B1"+("{:.1f}".format(100*Delta_hist_integral / hist_integral) if hist_integral else "?")+"%"+r"$^{\mathrm{st}}$"


    if linbksub_integral:
        integral_bksyst_Delta = abs(linbksub_integral - hist_integral)
        rawint_text = rawint_text + u"\u00B1" +("{:.1f}".format(100*integral_bksyst_Delta / hist_integral) if hist_integral else "?")+"%"+r"$^{\mathrm{bk}}$"
    rawint_text = rawint_text + u"\u00B1" +("{:.1f}".format(100*integral_range_std / hist_integral) if hist_integral else "?")+"%"+r"$^{\mathrm{r}}$"

    if not (isData or doBkg) and hist_fullrangeintegral:
        rawint_text += "\nFull range raw integral/bw = \n"+"{:.1f}".format(hist_fullrangeintegral)+u"\u00B1"+("{:.1f}".format(100*Delta_hist_fullrangeintegral / hist_fullrangeintegral) if hist_fullrangeintegral else "?")+"%"+r"$^{\mathrm{st}}$"

    ax.plot(x_axis,y_axis, "sb", label = rawint_text)
    ax.errorbar(x_axis, y_axis, yerr=y_axis_pm, xerr = x_axis_pm, capsize=3, fmt=".", ecolor = "black")
    ax.plot(x_axis_black, y_axis_black, "x", color="black")
    ax.set_xlabel("m [GeV]", fontsize=18, labelpad=3)
    ax.set_ylabel("Events", fontsize=18)
    ax.xaxis.set_label_coords(0.9, -0.05)


    passing_or_failing = binname[-4:]
    binname = binname[:-4]
    binname += "WP: "+config_file['workingpoint']

    ax.set_title(binname, fontsize = 20)
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
        ax.plot(x,Chev.chebval(x,chevparams),label=str(chevorder)+r"th order chev") #, $\chi^{2}$/dof="+chi2)


        # Draw the line BK estimate
        if (not index in config_file["data_skip_bksyst_estimation"] and isData) or (not index in config_file["mc_skip_bksyst_estimation"] and not isData):
            
            ax.plot(linbk_params[0], linbk_params[1], color="gray", ls ="--")


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

        fulltext = text + ", integral/bw = "+"\n{:.1f}".format(total_integral)+u"\u00B1"+("{:.1f}".format(100*total_integral_Delta / total_integral) if total_integral else "?")+"%"+r"$^{\mathrm{st}}$"
        
        if doBkg and cb_settings["bksyst_Delta_integral"]:
            total_integral_bksyst_Delta = abs(cb_settings["bksyst_Delta_integral"])
            fulltext = fulltext + u"\u00B1" +("{:.1f}".format(100*total_integral_bksyst_Delta / total_integral) if total_integral else "?")+"%"+r"$^{\mathrm{bk}}$"

        fulltext = fulltext + u"\u00B1"+("{:.1f}".format(100* cb_settings["integral_rangesyst_std"] / total_integral) if total_integral else "?")+"%"+r"$^{\mathrm{r}}$"

        ax.plot(x,cb_y,color="crimson",label=fulltext)

        


    if isData and index in config_file['data_override_plots_xrange']:
        ax.set_xlim(config_file['data_override_plots_xrange'][index])
        ax.set_xticks(np.arange(config_file['data_override_plots_xrange'][index][0], config_file['data_override_plots_xrange'][index][1], config_file['xticks_density']), minor=True)
    else:
        ax.set_xlim(config_file['plots_xrange'])
        ax.set_xticks(np.arange(config_file['plots_xrange'][0], config_file['plots_xrange'][1], config_file['xticks_density']), minor=True)


    ax.set_ylim(goodlim_y)
    leg = ax.legend(loc='best', title=("Data" if isData else "MC")+", "+passing_or_failing, fontsize=17)
    leg.get_title().set_fontsize(17)
    leg.get_frame().set_alpha(0.42)

    ax.grid(linestyle = "--")



    vlinecolor = 'gray' if index in config_file[("MC_" if not isData else "")+'dualfit_bounds_cb'] else 'k'
    
    #ax.vlines(config_file[("MC_" if not isData else "")+"raw_integral_bounds"][index][0], goodlim_y[0],goodlim_y[1], color=vlinecolor)
    #ax.vlines(config_file[("MC_" if not isData else "")+"raw_integral_bounds"][index][1], goodlim_y[0],goodlim_y[1], color=vlinecolor)
    ax.vlines(master_integral_range[0], goodlim_y[0],goodlim_y[1], color=vlinecolor)
    ax.vlines(master_integral_range[1], goodlim_y[0],goodlim_y[1], color=vlinecolor)





    if doBkg:
        if not os.path.exists(output_dir+"/fits/chev"+str(config_file["chevorder"])):
            os.makedirs(output_dir+"/fits/chev"+str(config_file["chevorder"]))
        plt.savefig(output_dir+"/fits/chev"+str(config_file["chevorder"])+"/figure"+str(index)+".png")
    else:
        if not os.path.exists(output_dir+"/fits/"):
            os.makedirs(output_dir+"/fits/")
        plt.savefig(output_dir+"/fits/figure"+str(index)+".png")
    plt.close()


# UNUSED
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


    chevorder = config_file["chevorder"]
    if index in config_file['chevorder_in_bin']:
        chevorder = config_file['chevorder_in_bin'][index]


    if len (data_plotinfo["chevstats"][0]) == 0:
        chi2 = "?"
    else:
        chi2 = "{:.1f}".format(data_plotinfo["chevstats"][0][0] / (len(data_plotinfo["chevx"]) - (chevorder + 1)))
    ax.plot(data_x,Chev.chebval(data_x,data_plotinfo["chevparams"]),label=str(chevorder)+r"th order chev, data, $\chi^{2}$/dof="+chi2)
    if mc_doBkg:
        if len (mc_plotinfo["chevstats"][0]) == 0:
            chi2 = "?"
        else:
            chi2 = "{:.1f}".format(mc_plotinfo["chevstats"][0][0] / (len(mc_plotinfo["chevx"]) - (chevorder + 1)))
        ax.plot(mc_x,Chev.chebval(mc_x,mc_plotinfo["chevparams"]),label=str(chevorder)+r"th order chev, MC, $\chi^{2}$/dof="+chi2)


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




def TnPEfficiency(mutex, effi, effierr, sea,
                       datafile, mcfile, mcfile_alt, bins, binnames, indexlist, binindex_range_to_fit, fitted_figures, workdir, config_file, skipped_lowpt_bins):



    bCB = False


    master_integral_range = config_file["master_integral_range"]
    first_one = True

    for i in range(0, len(indexlist), 2):
        pass_index = indexlist[i]
        fail_index = indexlist[i+1]

        fit_in_this_bin = binindex_range_to_fit[0] <= pass_index < binindex_range_to_fit[1]
        # olyan nincs hogy egyiket fit masikat nem, eleg pass-t nezni

        if fitted_figures != "all":
            fit_in_this_bin = pass_index in fitted_figures



        if "bUseTemplateFit" in config_file.keys() and config_file["bUseTemplateFit"]:


            plotinfo_window_pass, plotinfo_window_down_pass, plotinfo_window_up_pass, plotinfo_fr_pass = TnPTemplateEfficiencyInBin(datafile, mcfile, workdir, "dataplots/", bins[pass_index], binnames[pass_index], pass_index, config_file)
            plotinfo_window_fail, plotinfo_window_down_fail, plotinfo_window_up_fail, plotinfo_fr_fail = TnPTemplateEfficiencyInBin(datafile, mcfile, workdir, "dataplots/", bins[fail_index], binnames[fail_index], fail_index, config_file)

            # in mass window
            # histintegral - bk


            data_pass = plotinfo_window_pass["cc"][0]
            data_fail = plotinfo_window_fail["cc"][0]
            data_efficiency = data_pass / (data_pass + data_fail)

            data_efficiency_rangeup = plotinfo_window_up_pass["cc"][0] / (plotinfo_window_up_pass["cc"][0] + plotinfo_window_up_fail["cc"][0])
            data_efficiency_rangedown = plotinfo_window_down_pass["cc"][0] / (plotinfo_window_down_pass["cc"][0] + plotinfo_window_down_fail["cc"][0])


            data_efficiency_stat_error = get_efficiency_error(data_pass, data_fail, plotinfo_window_pass["cc"][1], plotinfo_window_fail["cc"][1])
            data_efficiency_rangesyst_error = np.std([data_efficiency, data_efficiency_rangeup, data_efficiency_rangedown])
            data_efficiency_total_error = data_efficiency_stat_error |qs| data_efficiency_rangesyst_error


            # from template fit (window)

            minuit_efficiency = plotinfo_window_pass["m"].values["x0"] / (plotinfo_window_pass["m"].values["x0"] + plotinfo_window_fail["m"].values["x0"])
            minuit_efficiency_error = get_efficiency_error(plotinfo_window_pass["m"].values["x0"], plotinfo_window_fail["m"].values["x0"], plotinfo_window_pass["m"].errors["x0"], plotinfo_window_fail["m"].errors["x0"])

            minuit_efficiency_down = plotinfo_window_down_pass["m"].values["x0"] / (plotinfo_window_down_pass["m"].values["x0"] + plotinfo_window_down_fail["m"].values["x0"])
            minuit_efficiency_up = plotinfo_window_up_pass["m"].values["x0"] / (plotinfo_window_up_pass["m"].values["x0"] + plotinfo_window_up_fail["m"].values["x0"])


            minuit_efficiency_down_error = get_efficiency_error(plotinfo_window_down_pass["m"].values["x0"], plotinfo_window_down_fail["m"].values["x0"], plotinfo_window_down_pass["m"].errors["x0"], plotinfo_window_down_fail["m"].errors["x0"])
            minuit_efficiency_up_error = get_efficiency_error(plotinfo_window_up_pass["m"].values["x0"], plotinfo_window_up_fail["m"].values["x0"], plotinfo_window_up_pass["m"].errors["x0"], plotinfo_window_up_fail["m"].errors["x0"])

            minuit_efficiency_range_error = np.std([minuit_efficiency, minuit_efficiency_down, minuit_efficiency_up])

            minuit_efficiency_total_error = minuit_efficiency_error |qs| minuit_efficiency_range_error


            # full mass range
            # histintegral - bk


            data_pass_fr = plotinfo_fr_pass["cc"][0]
            data_fail_fr = plotinfo_fr_fail["cc"][0]
            data_efficiency_fr = data_pass_fr / (data_pass_fr + data_fail_fr)

            data_efficiency_fr_error = get_efficiency_error(data_pass_fr, data_fail_fr, plotinfo_fr_pass["cc"][1], plotinfo_fr_fail["cc"][1])



            # from template fit (full range)

            minuit_fr_efficiency = plotinfo_fr_pass["m"].values["x0"] / (plotinfo_fr_pass["m"].values["x0"] + plotinfo_fr_fail["m"].values["x0"])
            minuit_fr_efficiency_error = get_efficiency_error(plotinfo_fr_pass["m"].values["x0"], plotinfo_fr_fail["m"].values["x0"], plotinfo_fr_pass["m"].errors["x0"], plotinfo_fr_fail["m"].errors["x0"])


            # monte carlo

            mc_pass = plotinfo_window_pass["mc_hist_integral"]
            mc_fail = plotinfo_window_fail["mc_hist_integral"]
            if not isinstance(mc_pass, tuple) or not isinstance(mc_fail, tuple):
                raise NotImplementedError


            mc_pass_down = plotinfo_window_down_pass["mc_hist_integral"]
            mc_fail_down = plotinfo_window_down_fail["mc_hist_integral"]
            if not isinstance(mc_pass_down, tuple) or not isinstance(mc_fail_down, tuple):
                raise NotImplementedError


            mc_pass_up = plotinfo_window_up_pass["mc_hist_integral"]
            mc_fail_up = plotinfo_window_up_fail["mc_hist_integral"]
            if not isinstance(mc_pass_up, tuple) or not isinstance(mc_fail_up, tuple):
                raise NotImplementedError
            
            mc_efficiency = get_efficiency(mc_pass[0],  mc_fail[0])
            mc_efficiency_stat_error = get_efficiency_error(mc_pass[0], mc_fail[0], mc_pass[1], mc_fail[1])

            mc_eff_down = get_efficiency(mc_pass_down[0],  mc_fail_down[0])
            mc_eff_up = get_efficiency(mc_pass_up[0],  mc_fail_up[0]) 

            mc_eff_down_stat_error = get_efficiency_error(mc_pass_down[0], mc_fail_down[0], mc_pass_down[1], mc_fail_down[1])
            mc_eff_up_stat_error = get_efficiency_error(mc_pass_up[0], mc_fail_up[0], mc_pass_up[1], mc_fail_up[1])


            mc_efficiency_rangesyst_error = np.std([mc_efficiency , mc_eff_down, mc_eff_up  ])
            mc_efficiency_total_error = mc_efficiency_stat_error |qs| mc_efficiency_rangesyst_error


            # fullrange mc
            mc_pass_fr = plotinfo_fr_pass["mc_hist_integral"]
            mc_fail_fr = plotinfo_fr_fail["mc_hist_integral"]
            if not isinstance(mc_pass_fr, tuple) or not isinstance(mc_fail_fr, tuple):
                raise NotImplementedError

            mc_efficiency_fr = get_efficiency(mc_pass_fr[0], mc_fail_fr[0])
            mc_efficiency_fr_total_error = get_efficiency_error(mc_pass_fr[0], mc_fail_fr[0], mc_pass_fr[1], mc_fail_fr[1])


            # template plotting

            minuit_template_plot(config_file, workdir, binnames, plotinfo_window_pass, plotinfo_window_fail, pass_index, fail_index, (minuit_efficiency, minuit_efficiency_error, minuit_efficiency_range_error), (mc_efficiency, mc_efficiency_stat_error, mc_efficiency_rangesyst_error), rangeplot="", minuit_range_efficiencies=(minuit_efficiency_down, minuit_efficiency_up), mc_range_efficiencies = (mc_eff_down, mc_eff_up))
            minuit_template_plot(config_file, workdir, binnames, plotinfo_window_down_pass, plotinfo_window_down_fail, pass_index, fail_index, (minuit_efficiency_down, minuit_efficiency_down_error), (mc_eff_down, mc_eff_down_stat_error), rangeplot="down")
            minuit_template_plot(config_file, workdir, binnames, plotinfo_window_up_pass, plotinfo_window_up_fail, pass_index, fail_index, (minuit_efficiency_up, minuit_efficiency_up_error), (mc_eff_up, mc_eff_up_stat_error), rangeplot="up")
            minuit_template_plot(config_file, workdir, binnames, plotinfo_fr_pass, plotinfo_fr_fail, pass_index, fail_index, (minuit_fr_efficiency, minuit_fr_efficiency_error), (mc_efficiency_fr, mc_efficiency_fr_total_error), rangeplot="fr")

            efficiency_loss_histogram_from_templatefit(config_file, workdir, pass_index//2, binnames[pass_index], plotinfo_fr_pass["x"], plotinfo_fr_pass["signal"], plotinfo_fr_pass["uncert_signal"], plotinfo_fr_fail["signal"], plotinfo_fr_fail["uncert_signal"], plotinfo_fr_pass["mc"], plotinfo_fr_pass["uncert_mc"], plotinfo_fr_fail["mc"], plotinfo_fr_fail["uncert_mc"], plotinfo_window_pass["x"], plotinfo_window_pass["signal"], plotinfo_window_pass["uncert_signal"], plotinfo_window_fail["signal"], plotinfo_window_fail["uncert_signal"], plotinfo_window_pass["mc"], plotinfo_window_pass["uncert_mc"], plotinfo_window_fail["mc"], plotinfo_window_fail["uncert_mc"])




        else:  # Regular Function Fitting

        
            pass_use_fixed_sigma_data = False
            pass_use_fixed_sigma_data_Delta = False
            fail_use_fixed_sigma_data = False
            fail_use_fixed_sigma_data_Delta = False
            if pass_index in config_file['fix_sigma_from_previous_bin_data']:
                if first_one:
                    print ("Error: Use different number of threads with fixed sigma from previous bin in bin "+str(pass_index))
                else:
                    print ("Fixing sigma for bin "+str(pass_index)+" to "+str(data_pass_plotinfo["sigma_cb_fit"]))
                    use_fixed_sigma_data = data_pass_plotinfo["sigma_cb_fit"]
                    use_fixed_sigma_data_Delta = data_pass_plotinfo["Delta_sigma_cb_fit"]
            if fail_index in config_file['fix_sigma_from_previous_bin_data']:
                if first_one:
                    print ("Error: Use different number of threads with fixed sigma from previous bin in bin "+str(fail_index))
                else:
                    print ("Fixing sigma for bin "+str(fail_index)+" to "+str(data_fail_plotinfo["sigma_cb_fit"]))
                    use_fixed_sigma_data = data_fail_plotinfo["sigma_cb_fit"]
                    use_fixed_sigma_data_Delta = data_fail_plotinfo["Delta_sigma_cb_fit"]

            first_one = False

            
            if effi is None:
                data_pass_plotinfo = TnPEfficiencyInBin(True, datafile, workdir, "dataplots/", bins[pass_index], binnames[pass_index], pass_index, fit_in_this_bin, config_file, pass_use_fixed_sigma_data, pass_use_fixed_sigma_data_Delta)
                data_fail_plotinfo = TnPEfficiencyInBin(True, datafile, workdir, "dataplots/", bins[fail_index], binnames[fail_index], fail_index, fit_in_this_bin, config_file, fail_use_fixed_sigma_data, fail_use_fixed_sigma_data_Delta)


                # Handle if plotinfo is False, i.e. error was found in sanity checks
                # This is needed to let the code run for std properly
                if data_pass_plotinfo is False:
                    print(f"Data Fit {pass_index} (Pass) completely failed, skipping bin with pass-fail indexes {pass_index}, {fail_index}")
                    continue
                if data_fail_plotinfo is False:
                    print(f"Data Fit {fail_index} (Fail) completely failed, skipping bin with pass-fail indexes {pass_index}, {fail_index}")
                    continue




                # get Data efficiency
                # bksub c&c


                if config_file['bUseMinuit']:
                    data_pass = data_pass_plotinfo['minuit_histintegral'][0]
                    data_fail = data_fail_plotinfo['minuit_histintegral'][0]
                    data_efficiency = data_pass / (data_pass + data_fail)

                    data_efficiency_rangeup = data_pass_plotinfo['minuit_up_histintegral'][0] / (data_pass_plotinfo['minuit_up_histintegral'][0] + data_fail_plotinfo['minuit_up_histintegral'][0])
                    data_efficiency_rangedown = data_pass_plotinfo['minuit_down_histintegral'][0] / (data_pass_plotinfo['minuit_down_histintegral'][0] + data_fail_plotinfo['minuit_down_histintegral'][0])


                    data_efficiency_stat_error = get_efficiency_error(data_pass, data_fail, data_pass_plotinfo['minuit_histintegral'][1], data_fail_plotinfo['minuit_histintegral'][1])
                    data_efficiency_rangesyst_error = np.std([data_efficiency, data_efficiency_rangeup, data_efficiency_rangedown])
                    data_efficiency_total_error = data_efficiency_stat_error |qs| data_efficiency_rangesyst_error

                else:
                    data_pass = data_pass_plotinfo["hist_integral"]
                    data_fail = data_fail_plotinfo["hist_integral"]
                    if not isinstance(data_pass, list) or not isinstance(data_fail, list):
                        raise NotImplementedError
                    
                    data_efficiency = data_pass[0] / (data_pass[0] + data_fail[0])
                    pass_integral_error_stat_plus_bksyst = data_pass_plotinfo["Delta_hist_integral"] |qs| data_pass_plotinfo["linkb_hist_integral_Delta"]
                    fail_integral_error_stat_plus_bksyst = data_fail_plotinfo["Delta_hist_integral"] |qs| data_fail_plotinfo["linkb_hist_integral_Delta"]
                    
                    data_efficiency_stat_error = get_efficiency_error(data_pass[0], data_fail[0], pass_integral_error_stat_plus_bksyst, fail_integral_error_stat_plus_bksyst)
                    data_efficiency_rangesyst_error = np.std([data_efficiency , data_pass[1] / (data_pass[1] + data_fail[1]) , data_pass[2] / (data_pass[2] + data_fail[2])])
                    data_efficiency_total_error = data_efficiency_stat_error |qs| data_efficiency_rangesyst_error

                # crystalball
                bCB = data_pass_plotinfo["cb_settings"]["bDraw"] and data_fail_plotinfo["cb_settings"]["bDraw"]
                if bCB:
                    cb_pass = data_pass_plotinfo["cb_settings"]["integral"]
                    cb_fail = data_fail_plotinfo["cb_settings"]["integral"]

                    if not isinstance(cb_pass, list) or not isinstance(cb_fail, list):
                        raise NotImplementedError

                    cb_efficiency = cb_pass[0]  / (cb_pass[0] + cb_fail[0])
                    pass_integral_error_stat_plus_bksyst = data_pass_plotinfo["cb_settings"]["Delta_integral"] |qs| data_pass_plotinfo["cb_settings"]["bksyst_Delta_integral"]
                    fail_integral_error_stat_plus_bksyst = data_fail_plotinfo["cb_settings"]["Delta_integral"] |qs| data_fail_plotinfo["cb_settings"]["bksyst_Delta_integral"]

                    cb_efficiency_stat_error = get_efficiency_error(cb_pass[0], cb_fail[0], pass_integral_error_stat_plus_bksyst, fail_integral_error_stat_plus_bksyst)
                    cb_efficiency_rangesyst_error = np.std([cb_efficiency , cb_pass[1] / (cb_pass[1] + cb_fail[1]) , cb_pass[2] / (cb_pass[2] + cb_fail[2])])
                    cb_efficiency_total_error = cb_efficiency_stat_error |qs| cb_efficiency_rangesyst_error








                ## Make Data plot
                if data_pass_plotinfo["minuit_info"] is not None:
                    # minuit out plot




                    minuit_efficiency = data_pass_plotinfo["minuit_info"]["m"].values["nsig"] / (data_pass_plotinfo["minuit_info"]["m"].values["nsig"] + data_fail_plotinfo["minuit_info"]["m"].values["nsig"])
                    minuit_efficiency_error = get_efficiency_error(data_pass_plotinfo["minuit_info"]["m"].values["nsig"], data_fail_plotinfo["minuit_info"]["m"].values["nsig"], data_pass_plotinfo["minuit_info"]["m"].errors["nsig"], data_fail_plotinfo["minuit_info"]["m"].errors["nsig"])

                    minuit_efficiency_down = data_pass_plotinfo["minuit_info_down"]["m"].values["nsig"] / (data_pass_plotinfo["minuit_info_down"]["m"].values["nsig"] + data_fail_plotinfo["minuit_info_down"]["m"].values["nsig"])
                    minuit_efficiency_up = data_pass_plotinfo["minuit_info_up"]["m"].values["nsig"] / (data_pass_plotinfo["minuit_info_up"]["m"].values["nsig"] + data_fail_plotinfo["minuit_info_up"]["m"].values["nsig"])


                    minuit_efficiency_down_error = get_efficiency_error(data_pass_plotinfo["minuit_info_down"]["m"].values["nsig"], data_fail_plotinfo["minuit_info_down"]["m"].values["nsig"], data_pass_plotinfo["minuit_info_down"]["m"].errors["nsig"], data_fail_plotinfo["minuit_info_down"]["m"].errors["nsig"])
                    minuit_efficiency_up_error = get_efficiency_error(data_pass_plotinfo["minuit_info_up"]["m"].values["nsig"], data_fail_plotinfo["minuit_info_up"]["m"].values["nsig"], data_pass_plotinfo["minuit_info_up"]["m"].errors["nsig"], data_fail_plotinfo["minuit_info_up"]["m"].errors["nsig"])

                    minuit_efficiency_range_error = np.std([minuit_efficiency, minuit_efficiency_down, minuit_efficiency_up])

                    minuit_efficiency_total_error = minuit_efficiency_error |qs| minuit_efficiency_range_error


                    # Minuit plot

                    minuit_data_plot(config_file, workdir, binnames, data_pass_plotinfo, data_fail_plotinfo, pass_index, fail_index, minuit_efficiency, minuit_efficiency_error, minuit_efficiency_range_error, rangeplot="")
                    minuit_data_plot(config_file, workdir, binnames, data_pass_plotinfo, data_fail_plotinfo, pass_index, fail_index, minuit_efficiency_down, minuit_efficiency_down_error, 0, rangeplot="down")
                    minuit_data_plot(config_file, workdir, binnames, data_pass_plotinfo, data_fail_plotinfo, pass_index, fail_index, minuit_efficiency_up, minuit_efficiency_up_error, 0, rangeplot="up")




                else:
                    # older odr out plot
                    def base_plot(ax, plotinfo, index, sProbeType):
                        hist_integral = plotinfo["hist_integral"][0]

                        rawint_text = "BKsub integral/bw = \n{:.1f}".format(hist_integral)+u"\u00B1"+("{:.1f}".format(100*plotinfo["Delta_hist_integral"] / hist_integral) if hist_integral else "?")+"%"+r"$^{\mathrm{stat}}$"

                        if plotinfo["linkb_hist_integral_Delta"]:
                            rawint_text = rawint_text + u"\u00B1" +("{:.1f}".format(100*plotinfo["linkb_hist_integral_Delta"] / hist_integral) if hist_integral else "?")+"%"+r"$^{\mathrm{bk}}$"

                        x = getSampledXForPlotting(config_file, index, isData=True, doBkg=True, cb_settings=plotinfo["cb_settings"])

                        x_axis = plotinfo["x_axis"]
                        y_axis = plotinfo["y_axis"]
                        x_axis_black = plotinfo["x_axis_black"]
                        y_axis_black = plotinfo["y_axis_black"]
                        x_axis_pm = plotinfo["x_axis_pm"]
                        y_axis_pm = plotinfo["y_axis_pm"]
                        ax.plot(x_axis,y_axis, "sb", label = rawint_text)
                        ax.errorbar(x_axis, y_axis, yerr=y_axis_pm, xerr = x_axis_pm, capsize=3, fmt=".", ecolor = "black")
                        ax.plot(x_axis_black, y_axis_black, "x", color="black")


                        goodlim_y = ax.get_ylim()

                        if index in config_file['data_override_plots_xrange']:
                            min_y = 9999
                            max_y = -1
                            for i in range(len(x_axis)):
                                if config_file['data_override_plots_xrange'][index][0] <= x_axis[i] <= config_file['data_override_plots_xrange'][index][-1]:
                                    if y_axis[i] < min_y:
                                        min_y = y_axis[i]
                                    if y_axis[i] > max_y:
                                        max_y = y_axis[i]

                            goodlim_y = [min_y*0.9, max_y*1.1]


                        chevstats = plotinfo["chevstats"]
                        chevx = plotinfo["chevx"]
                        chevparams = plotinfo["chevparams"]

                        chevorder = config_file["chevorder"]
                        if index in config_file['chevorder_in_bin']:
                            chevorder = config_file['chevorder_in_bin'][index]

                        if len (chevstats[0]) == 0:
                            chi2 = "?"
                        else:
                            chi2 = "{:.1f}".format(chevstats[0][0] / (len(chevx) - (chevorder + 1)))
                        ax.plot(x,Chev.chebval(x,chevparams),label=str(chevorder)+r"th order chev")


                        # Draw the line BK estimate
                        if (not index in config_file["data_skip_bksyst_estimation"]):
                            
                            ax.plot(plotinfo["linbk_params"][0], plotinfo["linbk_params"][1], color="gray", ls ="--")


                        if plotinfo["gauss_settings"]['bDraw']:
                            if config_file['bUseBWGauss']:
                                gauss_y = np.array(Voigt(x, plotinfo["gauss_settings"]['popt_sig'][0],plotinfo["gauss_settings"]['popt_sig'][1],plotinfo["gauss_settings"]['popt_sig'][2],plotinfo["gauss_settings"]['popt_sig'][3]) )
                            else:
                                gauss_y = np.array(gauss(x, plotinfo["gauss_settings"]['popt_sig'][0],plotinfo["gauss_settings"]['popt_sig'][1],plotinfo["gauss_settings"]['popt_sig'][2]) )

                            gauss_y = np.array(Chev.chebval(x,plotinfo["gauss_settings"]["chevparams"])) + gauss_y
                            ax.plot(x,Chev.chebval(x,plotinfo["gauss_settings"]["chevparams"]),label=str(chevorder)+r"th order chev for Gauss, $\chi^{2}$/dof="+chi2, color="magenta")

                            ax.plot(x,gauss_y,color="orange",label=("Breit-Wigner-" if config_file['bUseBWGauss'] else "") + "Gauss, integral/bw = "+"{:.1f}".format(plotinfo["gauss_settings"]["integral"])+u"\u00B1"+("{:.1f}".format(100*plotinfo["gauss_settings"]["Delta_integral"] / plotinfo["gauss_settings"]["integral"]) if plotinfo["gauss_settings"]["integral"] else "?")+"%")


                        if plotinfo["cb_settings"]['bDraw']:
                            cb_y = np.array(mycrystalball(x, plotinfo["cb_settings"]['popt_sig'][0],plotinfo["cb_settings"]['popt_sig'][1],plotinfo["cb_settings"]['popt_sig'][2],plotinfo["cb_settings"]['popt_sig'][3],plotinfo["cb_settings"]['popt_sig'][4])  )

                            # bk lift:
                            cb_y =  np.array(Chev.chebval(x,chevparams)) + cb_y


                            text = "Crystalball"
                            if plotinfo["cb_settings"]['fix_sigma']:
                                text = text + r"(fixed $\sigma=$"+"{:.0f}".format(plotinfo["cb_settings"]['fix_sigma']*1000)+" MeV)"


                            if plotinfo["psip_settings"]['bDraw']:
                                psip_y = np.array(gauss(x, plotinfo["psip_settings"]['popt_sig'][0],plotinfo["psip_settings"]['popt_sig'][1],plotinfo["psip_settings"]['popt_sig'][2]) )
                                cb_y_plus_psip = cb_y + psip_y

                                total_integral = plotinfo["cb_settings"]["integral"][0] + plotinfo["psip_settings"]["integral"][0]
                                total_integral_Delta = plotinfo["cb_settings"]["Delta_integral"] |qs| plotinfo["psip_settings"]["Delta_integral"]

                                fulltext=text + r"+$\psi$'(Gauss)" + ", integral/bw = "+"{:.1f}".format(total_integral)+u"\u00B1"+("{:.1f}".format(100*total_integral_Delta / total_integral) if total_integral else "?")+"%"
                                ax.plot(x,cb_y_plus_psip,color="darkmagenta",label=fulltext)
                            
                            total_integral = plotinfo["cb_settings"]["integral"][0]
                            total_integral_Delta = plotinfo["cb_settings"]["Delta_integral"]

                            fulltext = text + ", integral/bw = "+"\n{:.1f}".format(total_integral)+u"\u00B1"+("{:.1f}".format(100*total_integral_Delta / total_integral) if total_integral else "?")+"%"+r"$^{\mathrm{st}}$"
                            
                            if plotinfo["cb_settings"]["bksyst_Delta_integral"]:
                                total_integral_bksyst_Delta = abs(plotinfo["cb_settings"]["bksyst_Delta_integral"])
                                fulltext = fulltext + u"\u00B1" +("{:.1f}".format(100*total_integral_bksyst_Delta / total_integral) if total_integral else "?")+"%"+r"$^{\mathrm{bk}}$"

                            #fulltext = fulltext + u"\u00B1"+("{:.1f}".format(100* plotinfo["cb_settings"]["integral_rangesyst_std"] / total_integral) if total_integral else "?")+"%"+r"$^{\mathrm{r}}$"

                            ax.plot(x,cb_y,color="crimson",label=fulltext)

                            


                        if index in config_file['data_override_plots_xrange']:
                            ax.set_xlim(config_file['data_override_plots_xrange'][index])
                            ax.set_xticks(np.arange(config_file['data_override_plots_xrange'][index][0], config_file['data_override_plots_xrange'][index][1], config_file['xticks_density']), minor=True)
                        else:
                            ax.set_xlim(config_file['plots_xrange'])
                            ax.set_xticks(np.arange(config_file['plots_xrange'][0], config_file['plots_xrange'][1], config_file['xticks_density']), minor=True)


                        ax.set_ylim(goodlim_y)


                        display_binname = binnames[index][:4] + ", DATA " + sProbeType + " PROBES\n"
                        display_binname += binnames[index][6:-4] + " WP: " + config_file['workingpoint']
                        ax.set_title(display_binname, fontsize=11)

                        leg = ax.legend(loc='best', fontsize=10)
                        leg.get_title().set_fontsize(17)
                        leg.get_frame().set_alpha(0.42)

                        ax.grid(linestyle = "--")

                        vlinecolor = 'gray'
                        ax.vlines(master_integral_range[0], goodlim_y[0],goodlim_y[1], color=vlinecolor)
                        ax.vlines(master_integral_range[1], goodlim_y[0],goodlim_y[1], color=vlinecolor)


                    def pass_plot(ax):
                        return base_plot(ax, data_pass_plotinfo, pass_index, "PASSING")

                    def fail_plot(ax):
                        return base_plot(ax, data_fail_plotinfo, fail_index, "FAILING")

                    # build text info:
                    text_info = config_file["resonance"] + " DATA\n\n"
                    #text_info += str(binnames[pass_index][:-4])+"\n"
                    text_info += "C&C -BK in\n[%g,%g] GeV:\n%s=%.3f\n" % (master_integral_range[0], master_integral_range[1], r'$\eta$', data_efficiency)
                    text_info += u"\u00B1"+("{:.1f}".format(100*data_efficiency_stat_error / data_efficiency) if data_efficiency else "?")+"%"+r"$^{\mathrm{stat}}$"+"\n"
                    text_info += u"\u00B1"+("{:.1f}".format(100*data_efficiency_rangesyst_error / data_efficiency) if data_efficiency else "?")+"%"+r"$^{\mathrm{range}}$"+"\n\n"


                    text_info += "CB fit in\n[%g,%g] GeV:\n%s=%.3f\n" % (master_integral_range[0], master_integral_range[1], r'$\eta$', cb_efficiency)
                    text_info += u"\u00B1"+("{:.1f}".format(100*cb_efficiency_stat_error / cb_efficiency) if cb_efficiency else "?")+"%"+r"$^{\mathrm{stat}}$"+"\n"
                    text_info += u"\u00B1"+("{:.1f}".format(100*cb_efficiency_rangesyst_error / cb_efficiency) if cb_efficiency else "?")+"%"+r"$^{\mathrm{range}}$"+"\n"


                    plot_two_with_text(pass_plot, fail_plot, text_info, True, workdir+"dataplots/", str(binnames[pass_index][3]), config_file)


            ## MC

            mc_pass_plotinfo = TnPEfficiencyInBin(False, mcfile, workdir, "mcplots/", bins[pass_index], binnames[pass_index], pass_index, fit_in_this_bin, config_file, data_pass_plotinfo["sigma_cb_fit"] if config_file['resonance'] == "JPsi" else False, data_pass_plotinfo["Delta_sigma_cb_fit"] if config_file['resonance'] == "JPsi" else False)
            mc_fail_plotinfo = TnPEfficiencyInBin(False, mcfile, workdir, "mcplots/", bins[fail_index], binnames[fail_index], fail_index, fit_in_this_bin, config_file, data_fail_plotinfo["sigma_cb_fit"] if config_file['resonance'] == "JPsi" else False, data_fail_plotinfo["Delta_sigma_cb_fit"] if config_file['resonance'] == "JPsi" else False)


            # get MC efficiency
            # always cut and count
            # range
            mc_pass = mc_pass_plotinfo["hist_integral"]
            mc_fail = mc_fail_plotinfo["hist_integral"]
            if not isinstance(mc_pass, list) or not isinstance(mc_fail, list):
                raise NotImplementedError
            
            mc_efficiency = mc_pass[0] / (mc_pass[0] + mc_fail[0])
            mc_efficiency_stat_error = get_efficiency_error(mc_pass[0], mc_fail[0], mc_pass_plotinfo["Delta_hist_integral"], mc_fail_plotinfo["Delta_hist_integral"])
            mc_efficiency_rangesyst_error = np.std([mc_efficiency , mc_pass[1] / (mc_pass[1] + mc_fail[1]) , mc_pass[2] / (mc_pass[2] + mc_fail[2])])
            mc_efficiency_total_error = mc_efficiency_stat_error |qs| mc_efficiency_rangesyst_error

            # fullrange
            mc_pass_fr = mc_pass_plotinfo["hist_fullrangeintegral"]
            mc_fail_fr = mc_fail_plotinfo["hist_fullrangeintegral"]

            mc_efficiency_fr = mc_pass_fr / (mc_pass_fr + mc_fail_fr)
            mc_efficiency_fr_total_error = get_efficiency_error(mc_pass_fr, mc_fail_fr, mc_pass_plotinfo["Delta_hist_fullrangeintegral"], mc_fail_plotinfo["Delta_hist_fullrangeintegral"])


            bRangeSetToFullRangeOnly = False
            if (config_file["bUseFullRangeCutAndCountMC"] and master_integral_range[0] == mc_pass_plotinfo["fullrange_meaning"][0] and master_integral_range[1] == mc_pass_plotinfo["fullrange_meaning"][1]):
                bRangeSetToFullRangeOnly = True



            ## Make MC plot
            def base_plot_mc(ax, plotinfo, index, sProbeType):
                hist_integral = plotinfo["hist_integral"][0]
                rawint_text = ("Integral/bw = ")+"\n{:.1f}".format(hist_integral)+u"\u00B1"+("{:.1f}".format(100*plotinfo["Delta_hist_integral"] / hist_integral) if hist_integral else "?")+"%"+r"$^{\mathrm{stat}}$"

                if plotinfo["hist_fullrangeintegral"] and not bRangeSetToFullRangeOnly:
                    rawint_text += "\nFR Integral/bw = \n"+"{:.1f}".format(plotinfo["hist_fullrangeintegral"])+u"\u00B1"+("{:.1f}".format(100*plotinfo["Delta_hist_fullrangeintegral"] / plotinfo["hist_fullrangeintegral"]) if plotinfo["hist_fullrangeintegral"] else "?")+"%"+r"$^{\mathrm{stat}}$"

                x_axis = plotinfo["x_axis"]
                y_axis = plotinfo["y_axis"]
                x_axis_black = plotinfo["x_axis_black"]
                y_axis_black = plotinfo["y_axis_black"]
                x_axis_pm = plotinfo["x_axis_pm"]
                y_axis_pm = plotinfo["y_axis_pm"]
                ax.plot(x_axis,y_axis, "sb", label = rawint_text)
                ax.errorbar(x_axis, y_axis, yerr=y_axis_pm, xerr = x_axis_pm, capsize=3, fmt=".", ecolor = "black")
                ax.plot(x_axis_black, y_axis_black, "x", color="black")

                display_binname = binnames[index][:4] + ", MC " + sProbeType + " PROBES\n"
                display_binname += binnames[index][6:-4] + " WP: " + config_file['workingpoint']
                ax.set_title(display_binname, fontsize=11)

                leg = ax.legend(loc='best', fontsize=10)
                leg.get_frame().set_alpha(0.42)

                ax.grid(linestyle = "--")

                vlinecolor = 'gray'
                goodlim_y = ax.get_ylim()
                ax.vlines(master_integral_range[0], goodlim_y[0],goodlim_y[1], color=vlinecolor)
                ax.vlines(master_integral_range[1], goodlim_y[0],goodlim_y[1], color=vlinecolor)

            def pass_plot_mc(ax):
                return base_plot_mc(ax, mc_pass_plotinfo, pass_index, "PASSING")

            def fail_plot_mc(ax):
                return base_plot_mc(ax, mc_fail_plotinfo, fail_index, "FAILING")

            # build text info:
            text_info = config_file["collection"] +"\n" + config_file["resonance"] + " MC SIMULATION\n\n"
            #text_info += str(binnames[pass_index][:-4])+"\n"

            if not bRangeSetToFullRangeOnly:
                text_info += "C&C in\n[%g, %g] GeV:\n%s=%.1f" % (master_integral_range[0], master_integral_range[1], r'$\eta$', mc_efficiency*100) + "%\n"
                text_info += u"\u00B1"+("{:.1f}".format(100*mc_efficiency_stat_error) if mc_efficiency else "?")+"%"+r"$^{\mathrm{stat},(abs)}$"+"\n"
                text_info += u"\u00B1"+("{:.1f}".format(100*mc_efficiency_rangesyst_error) if mc_efficiency else "?")+"%"+r"$^{\mathrm{range},(abs)}$"+"\n\n"


            if not config_file["bUseFullRangeCutAndCountMC"]:
                text_info += "(Info)\n"
            text_info += "C&C in\n[%g, %g] GeV:\n%s=%.1f" % (mc_pass_plotinfo["fullrange_meaning"][0], mc_pass_plotinfo["fullrange_meaning"][1], r'$\eta$', mc_efficiency_fr*100) + "%\n"
            text_info += u"\u00B1"+("{:.1f}".format(100*mc_efficiency_fr_total_error) if mc_efficiency_fr else "?")+"%"+r"$^{\mathrm{stat},(abs)}$"


            plot_two_with_text(pass_plot_mc, fail_plot_mc, text_info, False, workdir+"mcplots/", str(binnames[pass_index][3]), config_file)


            efficiency_loss_histogram_from_minuitfit(config_file, workdir, pass_index//2, binnames[pass_index], 
                                                        data_pass_plotinfo["minuit_x"], 
                                                        data_pass_plotinfo["minuit_info"]["signal_binned"], 
                                                        data_pass_plotinfo["minuit_info"]["uncert_signal_binned"], 
                                                        data_fail_plotinfo["minuit_x"], 
                                                        data_fail_plotinfo["minuit_info"]["signal_binned"], 
                                                        data_fail_plotinfo["minuit_info"]["uncert_signal_binned"], 

                                                        mc_pass_plotinfo["x_axis_unprocessed"],
                                                        mc_pass_plotinfo["y_axis_unprocessed"],
                                                        mc_pass_plotinfo["y_axis_unprocessed_pm"],
                                                        mc_fail_plotinfo["x_axis_unprocessed"],
                                                        mc_fail_plotinfo["y_axis_unprocessed"],
                                                        mc_fail_plotinfo["y_axis_unprocessed_pm"]
                                                    )





        # save results for summarizer scripts
        if fit_in_this_bin:
            #makeCommonDataMCFitPlot(workdir+"/commonplots/", config_file, binnames[index], index, data_plotinfo, mc_plotinfo)
            #makeCommonSubtractedDataMCFitPlot(workdir+"/commonplots/", config_file, binnames[index], index, data_plotinfo, mc_plotinfo)

            with mutex:
                if effi is None:
                    sea.shared_data_cutcount_eff[pass_index // 2] = data_efficiency
                    sea.shared_data_cutcount_eff_total_err[pass_index // 2] = data_efficiency_total_error
                    sea.shared_data_cutcount_eff_stat_err[pass_index // 2] = data_efficiency_stat_error
                    sea.shared_data_cutcount_eff_syst_err[pass_index // 2] = data_efficiency_rangesyst_error

                    if bCB:
                        sea.shared_data_cbfit_eff[pass_index // 2] = cb_efficiency
                        sea.shared_data_cbfit_eff_total_err[pass_index // 2] = cb_efficiency_total_error
                        sea.shared_data_cbfit_eff_stat_err[pass_index // 2] = cb_efficiency_stat_error
                        sea.shared_data_cbfit_eff_syst_err[pass_index // 2] = cb_efficiency_rangesyst_error

                    sea.shared_data_minuit_eff[pass_index // 2] = minuit_efficiency 
                    sea.shared_data_minuit_eff_stat_err[pass_index // 2] = minuit_efficiency_error
                    sea.shared_data_minuit_eff_syst_err[pass_index // 2] = minuit_efficiency_range_error
                    sea.shared_data_minuit_eff_total_err[pass_index // 2] = minuit_efficiency_total_error

                    if "bUseTemplateFit" in config_file.keys() and config_file["bUseTemplateFit"]:
                        sea.shared_data_minuit_fr_eff[pass_index // 2] = minuit_fr_efficiency 
                        sea.shared_data_minuit_fr_eff_stat_err[pass_index // 2] = minuit_fr_efficiency_error

                sea.shared_mc_cutcount_eff[pass_index // 2] = mc_efficiency
                sea.shared_mc_cutcount_eff_total_err[pass_index // 2] = mc_efficiency_total_error
                sea.shared_mc_cutcount_eff_stat_err[pass_index // 2] = mc_efficiency_stat_error
                sea.shared_mc_cutcount_eff_syst_err[pass_index // 2] = mc_efficiency_rangesyst_error

                sea.shared_mc_cutcount_fr_eff[pass_index // 2] = mc_efficiency_fr
                sea.shared_mc_cutcount_fr_eff_total_err[pass_index // 2] = mc_efficiency_fr_total_error


                


        '''   FIXME: RE-IMPLEMENT FOR NEW GROUPING
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

                    shared_mc_alt_raw_fullrangeintegrals[index - 4*skipped_lowpt_bins] = mc_plotinfo["hist_fullrangeintegral"]
                    shared_mc_alt_raw_fullrangeintegral_sigmas[index - 4*skipped_lowpt_bins] = mc_plotinfo["Delta_hist_fullrangeintegral"]

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
        '''



def TnPTemplateEfficiencyInBin(datafile, mcfile, workdir, plotdir, bin, binname, index, config_file):

    datf = ROOT.TFile.Open(datafile, "read")
    mcf = ROOT.TFile.Open(mcfile, "read")

    hist_data = datf.Get(bin)
    hist_mc = mcf.Get(bin)

    
    x_data = []
    x_mc = []
    y_data = []
    y_mc = []


    x_data_up = []
    x_mc_up = []
    y_data_up = []
    y_mc_up = []


    x_data_down = []
    x_mc_down = []
    y_data_down = []
    y_mc_down = []

    x_data_fr = []
    x_mc_fr = []
    y_data_fr = []
    y_mc_fr = []




    fitbound_low = config_file["template_fit_bounds"][index][0] if index in config_file["template_fit_bounds"] else config_file["template_fit_master_bounds"][0]
    fitbound_high = config_file["template_fit_bounds"][index][1] if index in config_file["template_fit_bounds"] else config_file["template_fit_master_bounds"][1]


    down_lowshift = 2 if not index in config_file['minuit_rangedown_shift'] else config_file['minuit_rangedown_shift'][index][0]
    down_highshift = 2 if not index in config_file['minuit_rangedown_shift'] else config_file['minuit_rangedown_shift'][index][1]


    up_lowshift = 2 if not index in config_file['minuit_rangeup_shift'] else config_file['minuit_rangeup_shift'][index][0]
    up_highshift = 2 if not index in config_file['minuit_rangeup_shift'] else config_file['minuit_rangeup_shift'][index][1]


    # Get MC cc
    # find out the range
    histxaxis = hist_mc.GetXaxis()
    lowbinSystDown = -1
    lowbin = -1
    lowbinSystUp = -1
    highbin = -1
    highbinSystDown = -1
    highbinSystUp = -1
    for k in range(1,histxaxis.GetNbins() + 1):
        if lowbinSystDown == -1 and histxaxis.GetBinCenter(k) >= fitbound_low - down_lowshift:
            lowbinSystDown = k
        if lowbin == -1 and histxaxis.GetBinCenter(k) >= fitbound_low:
            lowbin = k
        if lowbinSystUp == -1 and histxaxis.GetBinCenter(k) >= fitbound_low + down_highshift:
            lowbinSystUp = k

        if highbinSystDown == -1 and histxaxis.GetBinCenter(k) > fitbound_high - up_lowshift:
            highbinSystDown = k - 1
        if highbin == -1 and histxaxis.GetBinCenter(k) > fitbound_high:
            highbin = k - 1
        if highbinSystUp == -1 and histxaxis.GetBinCenter(k) > fitbound_high - up_highshift:
            highbinSystUp = k - 1
            break

    hist_mc_integral_error = ctypes.c_double(0)
    hist_mc_integral = hist_mc.IntegralAndError(lowbin, highbin, hist_mc_integral_error)
    hist_mc_integral_error = hist_mc_integral_error.value

    hist_mc_fullrangeintegral_error = ctypes.c_double(0)
    hist_mc_fullrangeintegral = hist_mc.IntegralAndError(1, histxaxis.GetNbins(), hist_mc_fullrangeintegral_error)
    hist_mc_fullrangeintegral_error = hist_mc_fullrangeintegral_error.value


    hist_mc_integral_low_error = ctypes.c_double(0)
    hist_mc_integral_high_error = ctypes.c_double(0)
    hist_mc_integral_low  = hist_mc.IntegralAndError(lowbinSystDown,   highbinSystDown, hist_mc_integral_low_error)
    hist_mc_integral_high = hist_mc.IntegralAndError(lowbinSystUp, highbinSystUp,   hist_mc_integral_high_error)

    hist_mc_integral_low_error = hist_mc_integral_low_error.value
    hist_mc_integral_high_error = hist_mc_integral_high_error.value





    skipped_first_bins = 0
    for i in range(hist_data.GetNbinsX()):
        binval = hist_data.GetBinCenter(i)


        if binval < 1.0 or binval > config_file['overflow']: # overflow bins
            continue

        bincontent = hist_data.GetBinContent(i)


        # bincontent at 30 GeV might be there as 0 - don't throw off the fit
        if i < 3 and bincontent < 1:
            skipped_first_bins += 1
            continue

        x_data_fr.append(binval)
        y_data_fr.append(bincontent)




        if fitbound_low - down_lowshift < binval < fitbound_high - down_highshift:
            x_data_down.append(binval)
            y_data_down.append(bincontent)


        if fitbound_low + up_lowshift < binval < fitbound_high + up_highshift:
            x_data_up.append(binval)
            y_data_up.append(bincontent)


        if fitbound_low < binval < fitbound_high:
            x_data.append(binval)
            y_data.append(bincontent)


    for i in range(hist_mc.GetNbinsX()):
        binval = hist_mc.GetBinCenter(i)

        if binval < 1.0 or binval > config_file['overflow']: # overflow bins
            continue

        bincontent = hist_mc.GetBinContent(i)


        if skipped_first_bins > 0:
            skipped_first_bins -= 1
            continue

        x_mc_fr.append(binval)
        y_mc_fr.append(bincontent)



        down_lowshift = 2 if not index in config_file['minuit_rangedown_shift'] else config_file['minuit_rangedown_shift'][index][0]
        down_highshift = 2 if not index in config_file['minuit_rangedown_shift'] else config_file['minuit_rangedown_shift'][index][1]

        if fitbound_low - down_lowshift < binval < fitbound_high - down_highshift:
            x_mc_down.append(binval)
            y_mc_down.append(bincontent)

        up_lowshift = 2 if not index in config_file['minuit_rangeup_shift'] else config_file['minuit_rangeup_shift'][index][0]
        up_highshift = 2 if not index in config_file['minuit_rangeup_shift'] else config_file['minuit_rangeup_shift'][index][1]

        if fitbound_low + up_lowshift < binval < fitbound_high + up_highshift:
            x_mc_up.append(binval)
            y_mc_up.append(bincontent)


        if fitbound_low < binval < fitbound_high:
            x_mc.append(binval)
            y_mc.append(bincontent)



    assert np.allclose(x_data, x_mc)
    assert np.allclose(x_data_fr, x_mc_fr)


    # set range fit

    plotinfo_window = do_template_fit(np.asarray(x_data), np.asarray(y_data), np.asarray(y_mc), config_file, index)
    plotinfo_window_down = do_template_fit(np.asarray(x_data_down), np.asarray(y_data_down), np.asarray(y_mc_down), config_file, index)
    plotinfo_window_up = do_template_fit(np.asarray(x_data_up), np.asarray(y_data_up), np.asarray(y_mc_up), config_file, index)

    plotinfo_window["mc_hist_integral"] = (hist_mc_integral, hist_mc_integral_error, np.std([hist_mc_integral, hist_mc_integral_low, hist_mc_integral_high]))
    plotinfo_window["mc_hist_integral_range"] = (hist_mc_integral_low, hist_mc_integral_high)
    plotinfo_window_down["mc_hist_integral"] = (hist_mc_integral_low, hist_mc_integral_low_error)
    plotinfo_window_up["mc_hist_integral"] = (hist_mc_integral_high, hist_mc_integral_high_error)




    plotinfo_window["m_up"] = plotinfo_window_up["m"]
    plotinfo_window["cc_up"] = plotinfo_window_up["cc"]
    plotinfo_window["m_down"] = plotinfo_window_down["m"]
    plotinfo_window["cc_down"] = plotinfo_window_down["cc"]


    plotinfo_fr     = do_template_fit(np.asarray(x_data_fr), np.asarray(y_data_fr), np.asarray(y_mc_fr), config_file, index)
    plotinfo_fr["mc_hist_integral"] = (hist_mc_fullrangeintegral, hist_mc_fullrangeintegral_error)



    if config_file['bDebugReportParams']:
        print(f"Fit {index} mass window fit m: ")
        test_s = str(plotinfo_window["m"])
        print(test_s)
        print("\n\n")

        print(f"Fit {index} FULL RANGE fit m: ")
        print(plotinfo_fr["m"])


    # So the plan is to have a function, call for both shortrange and fullrange, it does the fits and the plot variables, returns a plotinfo, we pass those out and outside do the side-by-side pass and fail plot

    return plotinfo_window, plotinfo_window_down, plotinfo_window_up, plotinfo_fr



def TnPEfficiencyInBin(isData, inputfile_name, workdir, plotdir, bin, binname, index, bDoFit, config_file, fix_sigma, fix_sigma_Delta):

    doBkg = isData or config_file['fit_background_on_MC']

    master_integral_range = config_file["master_integral_range"]    # else [2.7,3.5]
    master_integral_systlow_range = [71, 100]  if config_file['resonance'] == "Z" else [2.55 ,3.33]
    master_integral_systhigh_range = [78, 110]  if config_file['resonance'] == "Z" else [2.83 ,3.67]

    inputfile = ROOT.TFile.Open(inputfile_name, "read")
    background_fit_bounds = config_file[('MC_' if not isData else '') + 'background_fit_bounds']
    peak_fit_bounds = config_file[('MC_' if not isData else '') + 'peak_fit_bounds']
    peak_fit_bounds_cb = config_file[('MC_' if not isData else '') + 'peak_fit_bounds_cb']
    background_fit_bounds_cmsshape = config_file[('MC_' if not isData else '') + 'background_fit_bounds_cmsshape']

    hist = inputfile.Get(bin)  

    if not ((index in config_file['skip_rebin_data'] and isData) or (index in config_file['skip_rebin_mc'] and not isData)):
        hist.Rebin()


    minuit_info = None
    minuit_info_down = None
    minuit_info_up = None

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
        y_axis_unprocessed_pm.append(hist.GetBinError(i) if hist.GetBinError(i) > 0 else 1)
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


    minuit_x = []
    minuit_x_down = []
    minuit_x_up = []
    minuit_y = []
    minuit_y_up = []
    minuit_y_down = []

    minuit_histintegral = None
    minuit_histintegral_error = None
    
    minuit_up_histintegral = None
    minuit_up_histintegral_error = None
    minuit_down_histintegral = None
    minuit_down_histintegral_error = None



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

    if (config_file['FilterZerosFromMCFitInput'] and not isData) or (config_file['FilterZerosFromDataFitInput'] and isData):
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

            if config_file['bUseMinuit'] and index in config_file['minuit_bounds']:
                if config_file['minuit_bounds'][index][0] < x_axis[i] < config_file['minuit_bounds'][index][1]:
                    minuit_x.append(x_axis[i])
                    minuit_y.append(y_axis[i])


                down_lowshift = 0.2 if not index in config_file['minuit_rangedown_shift'] else config_file['minuit_rangedown_shift'][index][0]
                down_highshift = 0.2 if not index in config_file['minuit_rangedown_shift'] else config_file['minuit_rangedown_shift'][index][1]


                if config_file['minuit_bounds'][index][0] - down_lowshift < x_axis[i] < config_file['minuit_bounds'][index][1] - down_highshift:
                    minuit_x_down.append(x_axis[i])
                    minuit_y_down.append(y_axis[i])

                up_lowshift = 0.2 if not index in config_file['minuit_rangeup_shift'] else config_file['minuit_rangeup_shift'][index][0]
                up_highshift = 0.2 if not index in config_file['minuit_rangeup_shift'] else config_file['minuit_rangeup_shift'][index][1]

                if config_file['minuit_bounds'][index][0] + up_lowshift < x_axis[i] < config_file['minuit_bounds'][index][1] + up_highshift:
                    minuit_x_up.append(x_axis[i])
                    minuit_y_up.append(y_axis[i])

    else:
        y_axis_unprocessed_pm_zerofix = []   # set zeros to 1 so that divisions work
        for y, y_pm in zip(y_axis_unprocessed, y_axis_unprocessed_pm):
            assert y >= 0
            y_axis_unprocessed_pm_zerofix.append(y_pm if y > 0 else 1)


        for i in range(len(x_axis_unprocessed)):
            x_above_psip =              (background_fit_bounds[index][4] < x_axis_unprocessed[i] < background_fit_bounds[index][5]) if len(background_fit_bounds[index]) > 5 else False
            x_above_unknown_resonance = (background_fit_bounds[index][6] < x_axis_unprocessed[i] < background_fit_bounds[index][7]) if len(background_fit_bounds[index]) > 7 else False
            if ((background_fit_bounds_cmsshape[index][0] < x_axis_unprocessed[i] < background_fit_bounds_cmsshape[index][1]) or (background_fit_bounds_cmsshape[index][2] < x_axis_unprocessed[i] < background_fit_bounds_cmsshape[index][3])):
                fitx.append(x_axis_unprocessed[i])
                fity.append(y_axis_unprocessed[i])

            if peak_fit_bounds[index][0] < x_axis_unprocessed[i] < peak_fit_bounds[index][1]:
                peakx.append(x_axis_unprocessed[i])
                peaky.append(y_axis_unprocessed[i])
                peaky_pm.append(y_axis_unprocessed_pm_zerofix[i])

            if peak_fit_bounds_cb[index][0] < x_axis_unprocessed[i] < peak_fit_bounds_cb[index][1]:
                peak_cb_x.append(x_axis_unprocessed[i])
                peak_cb_y.append(y_axis_unprocessed[i])
                peak_cb_y_pm.append(y_axis_unprocessed_pm_zerofix[i])

            if ((background_fit_bounds[index][0] < x_axis_unprocessed[i] < background_fit_bounds[index][1]) or (background_fit_bounds[index][2] < x_axis_unprocessed[i] < background_fit_bounds[index][3])) or x_above_psip or x_above_unknown_resonance:
                chevx.append(x_axis_unprocessed[i])
                chevy.append(y_axis_unprocessed[i])
                chevy_pm.append(y_axis_unprocessed_pm_zerofix[i])

            if index in config_file['psip_fit_bounds']:
                if config_file['psip_fit_bounds'][index][0] < x_axis_unprocessed[i] < config_file['psip_fit_bounds'][index][1]:
                    psip_x.append(x_axis_unprocessed[i])
                    psip_y.append(y_axis_unprocessed[i])
                    psip_y_pm.append(y_axis_unprocessed_pm_zerofix[i])

            if index in config_file[('MC_' if not isData else '') + 'dualfit_bounds_cb']:
                if config_file[('MC_' if not isData else '') + 'dualfit_bounds_cb'][index][0] < x_axis_unprocessed[i] < config_file[('MC_' if not isData else '') + 'dualfit_bounds_cb'][index][1]:
                    dualfit_x.append(x_axis_unprocessed[i])
                    dualfit_x_pm.append(x_axis_unprocessed_pm[i])
                    dualfit_y.append(y_axis_unprocessed[i])
                    dualfit_y_pm.append(y_axis_unprocessed_pm_zerofix[i])

            if config_file['bUseMinuit'] and index in config_file['minuit_bounds']:
                if config_file['minuit_bounds'][index][0] < x_axis_unprocessed[i] < config_file['minuit_bounds'][index][1]:
                    minuit_x.append(x_axis_unprocessed[i])
                    minuit_y.append(y_axis_unprocessed[i])


                down_lowshift = 0.2 if not index in config_file['minuit_rangedown_shift'] else config_file['minuit_rangedown_shift'][index][0]
                down_highshift = 0.2 if not index in config_file['minuit_rangedown_shift'] else config_file['minuit_rangedown_shift'][index][1]

                if config_file['minuit_bounds'][index][0] - down_lowshift < x_axis_unprocessed[i] < config_file['minuit_bounds'][index][1] - down_highshift:
                    minuit_x_down.append(x_axis_unprocessed[i])
                    minuit_y_down.append(y_axis_unprocessed[i])


                up_lowshift = 0.2 if not index in config_file['minuit_rangeup_shift'] else config_file['minuit_rangeup_shift'][index][0]
                up_highshift = 0.2 if not index in config_file['minuit_rangeup_shift'] else config_file['minuit_rangeup_shift'][index][1]

                if config_file['minuit_bounds'][index][0] + up_lowshift < x_axis_unprocessed[i] < config_file['minuit_bounds'][index][1] + up_highshift:
                    minuit_x_up.append(x_axis_unprocessed[i])
                    minuit_y_up.append(y_axis_unprocessed[i])

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

        curr_ylim = ax.get_ylim()
        max_y = 0
        min_y = 99999
        for i in range(len(x_axis)):
            if config_file['data_override_plots_xrange'][index][0] <= x_axis[i]:
                if y_axis[i] > max_y: 
                    max_y = y_axis[i]
                if y_axis[i] < min_y: 
                    min_y = y_axis[i]
            if config_file['data_override_plots_xrange'][index][1] < x_axis[i]:
                break

        if max_y < curr_ylim[1]:
            ax.set_ylim(curr_ylim[0], max_y)
            if min_y > curr_ylim[0]:
                ax.set_ylim(0.95*min_y, max_y*1.05)

    else:
        ax.set_xlim(config_file['plots_xrange'])

    ax.set_title(binname)
    if not os.path.exists(workdir+plotdir):
        os.makedirs(workdir+plotdir)
    plt.savefig(workdir+plotdir+"figure"+str(index)+".png")
    plt.close()

    if not bDoFit:
        return False

    # another sanity check
    if len(x_axis) < 3:
        print (("Warning: Data " if isData else "Warning: MC ") + "Figure "+str(index)+" is likely empty, skip from fits")
        index += 1
        return False


    # determine do fit
    bDoFit = isData or config_file["fit_MC"]


    print (("Data Fit " if isData else "MC c&c") +str(index)+", "+binname)
    if doBkg and bDoFit:
        weights = []
        for sigma in chevy_pm:
            weights.append(1/sigma)

        chevorder = config_file["chevorder"]
        if index in config_file['chevorder_in_bin']:
            chevorder = config_file['chevorder_in_bin'][index]

        chevparams, chevstats = Chev.chebfit(chevx, chevy, chevorder, full=True, w=weights)

        A = Chev.chebvander(chevx, chevorder)
        dof = len(chevx) - (chevorder + 1)  # degrees of freedom
        resid_var = chevstats[0][0] / dof
        cov = resid_var * np.linalg.inv(np.dot(A.T, A))
        chev_std_devs = np.sqrt(np.diag(cov))

        chevparams_gauss = list(chevparams)



        if config_file["bFitGauss"]:
            gauss_y_sub, gauss_y_pm_sub = subtract_background(peakx, peaky, peaky_pm, chevparams_gauss, True, chev_std_devs)

            try:
                gauss_popt_sig, gauss_perr_sig, gauss_integral, gauss_Delta_integral = peak_fit("gauss" if not config_file['bUseBWGauss'] else "bwg", peakx, gauss_y_sub, gauss_y_pm_sub, config_file, isData, index)
            except Exception as er:
                print ("Data " if isData else "MC ") + "Gauss Fit "+str(index)+": "+str(er)
                gauss_popt_sig, gauss_perr_sig, gauss_integral, gauss_Delta_integral = False, False, False, False

            if config_file['bUseBWGauss'] and index in config_file[('MC_' if not isData else '') + 'dualfit_bounds_bwg'] and gauss_integral:
                odrout = odr_dualfit_voigt(dualfit_x, dualfit_y, dualfit_x_pm, dualfit_y_pm, chevparams_gauss, gauss_popt_sig)
                gauss_popt_sig = odrout.beta[:4]
                gauss_perr_sig = odrout.sd_beta[:4]
                chevparams_gauss = odrout.beta[4:]  


                cov_full = odrout.cov_beta * odrout.res_var
                cov_signal = cov_full[:4, :4]

                gauss_integral, gauss_Delta_integral = getIntegralErrorWithMCMethodCov(Voigt, gauss_popt_sig, cov_signal, abs(peakx[1] - gauss_y_sub[0]), config_file, index)



            gauss_settings = {
                "bDraw": bool(gauss_integral),
                "popt_sig": gauss_popt_sig,
                "perr_sig": gauss_perr_sig,
                "integral": gauss_integral,
                "Delta_integral": gauss_Delta_integral,
                #"bksyst_Delta_integral": integral_w_linbk - cb_integral,
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
                        print ("averaged to "+str(cb_y_sub[i])+"at "+str(peak_cb_x[i]))
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
                #else:
                #    print("Slicing point at "+str(peak_cb_x[i]))
            
            peak_cb_x_forsub = peak_x_slice
            cb_y_sub = peak_slice
            cb_y_pm_sub = peak_pm_slice
        


        if config_file["bFitCrystalball"]:
            try:
                cb_popt_sig, cb_perr_sig, cb_integral, cb_Delta_integral = peak_fit("crystalball", peak_cb_x_forsub, cb_y_sub, cb_y_pm_sub, config_file, isData, index, fix_sigma, fix_sigma_Delta)    ### change: added fix_sigma for bad bin
            except Exception as er:
                print ("Data " if isData else "MC ") + "CB Fit "+str(index)+": "+str(er)
                cb_popt_sig, cb_perr_sig, cb_integral, cb_Delta_integral = False, False, False, False

            if index in config_file[('MC_' if not isData else '') + 'dualfit_bounds_cb'] and cb_integral:
                odrout = odr_dualfit_cb(dualfit_x, dualfit_y, dualfit_x_pm, dualfit_y_pm, chevparams, cb_popt_sig)
                cb_popt_sig = odrout.beta[:5]
                cb_perr_sig = odrout.sd_beta[:5]
                chevparams = odrout.beta[5:]
                #chev_std_devs = odrout.sd_beta[5:]     # FIXME: not currect, but... PUT BACK NOW?

                #test = []
                #for i in range(len(chevparams)):
                #    test.append(chev_std_devs[i]/abs(chevparams[i])*100)
                #print test

                cov_full = odrout.cov_beta * odrout.res_var
                cov_signal = cov_full[:5, :5]

                #cb_integral, cb_Delta_integral = getIntegralErrorWithMCMethodCov(mycrystalball, cb_popt_sig, cov_signal, abs(peak_cb_x[1] - peak_cb_x[0]), config_file, index)
                #print("ODR: ",cb_perr_sig)
                #print("ODR: ",cb_Delta_integral/cb_integral)

                # triple fit...
                #cb_popt_sig, cb_perr_sig, chevparams, chev_std_devs, cov_signal = curve_dualfit_cb(dualfit_x, dualfit_y, dualfit_y_pm, odrout.beta)
                #print("C: ",cb_perr_sig)
                #print("C: ",cb_Delta_integral/cb_integral)



                cb_integral, cb_Delta_integral = getIntegralErrorWithMCMethodCov(mycrystalball, cb_popt_sig, cov_signal, abs(peak_cb_x[1] - peak_cb_x[0]), config_file, index)


                #print("OG Signal yield = %.1f ± %.1f"
                #    % (cb_integral[0], cb_Delta_integral))
                #print("Signal yield = %.1f ± %.1f"
                #    % (m.values["nsig"], m.errors["nsig"]))


            if not index in config_file["data_skip_bksyst_estimation"] and cb_integral:
                # estimate error from just linear background...
                # area = integral(sig+bk) - integral(linear bk)
                # first part:
                fullparams = list(cb_popt_sig)
                for i in range(len(chevparams)):
                    fullparams.append(chevparams[i])
                integral_w_sigbk, _ = getIntegralErrorWithMCMethod(crystalballPlusChev, fullparams, False, abs(peak_cb_x[1] - peak_cb_x[0]), config_file, index, isData, bEstimateError=False)

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
                "bksyst_Delta_integral": 0,       # integral_w_linbk - cb_integral, # placeholder!
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


        # Note for future me: The raw hist integral and stuff MUST be here after fitting the peak, because it relies on the
        # chev background, which must come after the dual fitting
        # from minuit bk is further down

        
        if index in config_file[('MC_' if not isData else '') + 'dualfit_bounds_cb']:
            fit_range = config_file[('MC_' if not isData else '') + 'dualfit_bounds_cb'][index]
        else:
            fit_range = [config_file[('MC_' if not isData else '') + 'background_fit_bounds'][index][0],  config_file[('MC_' if not isData else '') + 'background_fit_bounds'][index][-1]]


        bksub_integral, bksub_integral_pm, linbksub_integral, linbk_params = get_bksub_integral(x_axis, y_axis, y_axis_pm, chevparams, chev_std_devs, master_integral_range, fit_range, config_file["integral_error_MCmethod_iterations"], config_file["resonance"])

        Delta_funcInt_from_linbk = getIntegralFuncErrorFromLinBK(chevparams, abs(x_axis[1] - x_axis[0]), config_file, linbk_params)


        hist_fullrangeintegral, hist_fullrangeintegral_error = False, False


        cb_settings["bksyst_Delta_integral"] = Delta_funcInt_from_linbk

        if index in config_file["data_skip_bksyst_estimation"]:
            linbksub_integral = False
            cb_settings["bksyst_Delta_integral"] = False



        #makePlot(workdir+plotdir, config_file, binname, index, isData, x_axis, x_axis_pm, y_axis, y_axis_pm, x_axis_black, y_axis_black, doBkg, chevx, chevparams, chevstats, bksub_integral, bksub_integral_pm, 0, 0, linbksub_integral, linbk_params, rangesyst_std, gauss_settings, cb_settings, psip_settings)

        if config_file["bUseMinuit"]:
            # iminuit test for data, minuit fit
            p0 = {}
            if index in config_file['minuit_p0']:
                p0 = config_file['minuit_p0'][index]

            limits = {}
            if index in config_file['minuit_limits']:
                limits = config_file['minuit_limits'][index]



            minuit_info = do_minuit_fit(minuit_x, minuit_y, p0, limits)

            if config_file['bDebugReportParams']:
                print(f"Fit {index} values: ",minuit_info["m"].values)


            minuit_info_down = do_minuit_fit(minuit_x_down, minuit_y_down, p0, limits)
            minuit_info_up = do_minuit_fit(minuit_x_up, minuit_y_up, p0, limits)


            minuit_histintegral, minuit_histintegral_error = get_minuit_data_cc(minuit_y, minuit_info["m"])
            minuit_up_histintegral, minuit_up_histintegral_error = get_minuit_data_cc(minuit_y_up, minuit_info_up["m"])
            minuit_down_histintegral, minuit_down_histintegral_error = get_minuit_data_cc(minuit_y_down, minuit_info_down["m"])

            #print(m.values["nsig"], minuit_histintegral)


    else:

        if config_file["bFitGauss"] and bDoFit:
            try:
                gauss_popt_sig, gauss_perr_sig, gauss_integral, gauss_Delta_integral, gauss_integral_rangesyst_std = peak_fit("gauss" if not config_file['bUseBWGauss'] else "bwg", peakx, peaky, peaky_pm, config_file, isData, index)
            except Exception as er:
                print ("Data " if isData else "MC ") + "Fit "+str(index)+": "+str(er)
                gauss_popt_sig, gauss_perr_sig, gauss_integral, gauss_Delta_integral, gauss_integral_rangesyst_std = False, False, False, False

            gauss_settings = {
                "bDraw": bool(gauss_integral[0]),
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


        if config_file["bFitCrystalball"] and bDoFit:

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
            if lowbinSystDown == -1 and histxaxis.GetBinCenter(k) >= master_integral_systlow_range[0]:
                lowbinSystDown = k
            if lowbin == -1 and histxaxis.GetBinCenter(k) >= master_integral_range[0]:
                lowbin = k
            if lowbinSystUp == -1 and histxaxis.GetBinCenter(k) >= master_integral_systhigh_range[0]:
                lowbinSystUp = k

            if highbinSystDown == -1 and histxaxis.GetBinCenter(k) > master_integral_systlow_range[1]:
                highbinSystDown = k - 1
            if highbin == -1 and histxaxis.GetBinCenter(k) > master_integral_range[1]:
                highbin = k - 1
            if highbinSystUp == -1 and histxaxis.GetBinCenter(k) > master_integral_systhigh_range[1]:
                highbinSystUp = k - 1
                break

        hist_integral_error = ctypes.c_double(0)
        hist_integral = hist.IntegralAndError(lowbin, highbin, hist_integral_error)
        hist_integral_error = hist_integral_error.value

        hist_fullrangeintegral_error = ctypes.c_double(0)
        hist_fullrangeintegral = hist.IntegralAndError(1, histxaxis.GetNbins(), hist_fullrangeintegral_error)
        hist_fullrangeintegral_error = hist_fullrangeintegral_error.value


        hist_integral_low_error = ctypes.c_double(0)
        hist_integral_low  = hist.IntegralAndError(lowbinSystDown,   highbinSystDown, hist_integral_low_error)
        hist_integral_high = hist.IntegralAndError(lowbinSystUp, highbinSystUp,   hist_integral_low_error)


        #hist_integral_rangesyst_std = useBigger(abs(hist_integral - hist_integral_high), abs(hist_integral - hist_integral_low))

        #makePlot(workdir+plotdir, config_file, binname, index, isData, x_axis, x_axis_pm, y_axis, y_axis_pm, x_axis_black, y_axis_black, doBkg, False, False, False, hist_integral, hist_integral_error, hist_fullrangeintegral, hist_fullrangeintegral_error, 0, [], hist_integral_rangesyst_std, gauss_settings, cb_settings, psip_settings)

    if doBkg:
        raw_integral = bksub_integral
        raw_integral_sigma = bksub_integral_pm
        linkb_hist_integral_Delta = abs(linbksub_integral - raw_integral[0]) if linbksub_integral else False
    else: 
        raw_integral = [hist_integral, hist_integral_low, hist_integral_high]
        raw_integral_sigma = hist_integral_error
        linkb_hist_integral_Delta = False


    if cb_integral:
        return_sigma_for_fixed_width = cb_popt_sig[4]
        return_Delta_sigma_for_fixed_width = cb_perr_sig[4]
    else:
        return_sigma_for_fixed_width = False
        return_Delta_sigma_for_fixed_width = False


    fullrange_meaning = (hist.GetXaxis().GetXmin(), hist.GetXaxis().GetXmax())



    plotinfo = {
        "x_axis": x_axis,
        "x_axis_pm": x_axis_pm,
        "y_axis": y_axis,
        "y_axis_pm": y_axis_pm,
        "x_axis_black": x_axis_black,
        "x_axis_unprocessed": x_axis_unprocessed,
        "x_axis_unprocessed_pm": x_axis_unprocessed_pm,
        "y_axis_unprocessed": y_axis_unprocessed,
        "y_axis_unprocessed_pm": y_axis_unprocessed_pm,
        "dualfit_x": dualfit_x,
        "dualfit_y": dualfit_y,
        "y_axis_black": y_axis_black,
        "chevx": chevx if doBkg else False,
        "chevparams": chevparams if doBkg else False,
        "chevstats": chevstats if doBkg else False,
        "chev_std_devs": chev_std_devs if doBkg else False,
        "hist_integral": raw_integral,
        "Delta_hist_integral": raw_integral_sigma,
        "hist_fullrangeintegral": hist_fullrangeintegral,
        "Delta_hist_fullrangeintegral": hist_fullrangeintegral_error,
        "fullrange_meaning": fullrange_meaning,
        "linkb_hist_integral_Delta":  linkb_hist_integral_Delta,
        "linbk_params": linbk_params if doBkg else None,
        "gauss_settings": gauss_settings,
        "cb_settings": cb_settings,
        "psip_settings": psip_settings,
        "sigma_cb_fit": return_sigma_for_fixed_width,
        "Delta_sigma_cb_fit": return_Delta_sigma_for_fixed_width,
        "minuit_info": minuit_info,
        "minuit_info_down": minuit_info_down,
        "minuit_info_up": minuit_info_up,
        "minuit_x": minuit_x,
        "minuit_y": minuit_y,
        "minuit_x_down": minuit_x_down,
        "minuit_y_down": minuit_y_down,
        "minuit_x_up": minuit_x_up,
        "minuit_y_up": minuit_y_up,
        "minuit_histintegral": [minuit_histintegral, minuit_histintegral_error],
        "minuit_down_histintegral": [minuit_down_histintegral, minuit_down_histintegral_error],
        "minuit_up_histintegral": [minuit_up_histintegral, minuit_up_histintegral_error],

    }

    inputfile.Close()
    #return plotinfo, raw_integral, raw_integral_sigma, gauss_integral, gauss_Delta_integral, cb_integral, cb_Delta_integral, return_sigma_for_fixed_width, return_Delta_sigma_for_fixed_width
    return plotinfo



def efficiency_loss_histogram_from_templatefit(config_file, workdir, index, binname, x, signal_pass, uncert_signal_pass, signal_fail, uncert_signal_fail, mc_pass, uncert_mc_pass, mc_fail, uncert_mc_fail, w_x, w_signal_pass, w_uncert_signal_pass, w_signal_fail, w_uncert_signal_fail, w_mc_pass, w_uncert_mc_pass, w_mc_fail, w_uncert_mc_fail):



    x_data = []
    y_data = []
    y_data_pm = []
    y_mc = []
    y_mc_pm = []
    x_mc = []


    
    x_sf = []
    y_sf = []
    y_sf_pm = []


    for p, p_pm, f, f_pm, p_mc, p_mc_pm, f_mc, f_mc_pm, xi in zip(signal_pass, uncert_signal_pass, signal_fail, uncert_signal_fail, mc_pass, uncert_mc_pass, mc_fail, uncert_mc_fail, x):
        if p+f > 0:
            data_eff = p / (p + f)
            data_eff_pm = get_efficiency_error(p, f, p_pm, f_pm)
            y_data.append(data_eff)
            y_data_pm.append( data_eff_pm  )
            x_data.append(xi)
        if p_mc+f_mc > 0:
            mc_eff = p_mc / (p_mc + f_mc)
            mc_eff_pm = get_efficiency_error(p_mc, f_mc, p_mc_pm, f_mc_pm)
            y_mc.append(mc_eff)
            y_mc_pm.append( mc_eff_pm )
            x_mc.append(xi)

        if p+f > 0 and p_mc+f_mc > 0:
            if mc_eff != 0:
                sf_pm = getSFerrorbar(data_eff, mc_eff, data_eff_pm, mc_eff_pm   ) 
                if sf_pm > 2:  # leave off huge points that distort scale?
                    continue
                x_sf.append(xi)
                y_sf.append(data_eff / mc_eff)
                y_sf_pm.append( sf_pm )


    w_y_data = []
    w_y_data_pm = []
    w_x_data = []


    w_x_sf = []
    w_y_sf = []
    w_y_sf_pm = []

    for p, p_pm, f, f_pm, p_mc, p_mc_pm, f_mc, f_mc_pm, xi in zip(w_signal_pass, w_uncert_signal_pass, w_signal_fail, w_uncert_signal_fail, w_mc_pass, w_uncert_mc_pass, w_mc_fail, w_uncert_mc_fail, w_x):
        if p+f > 0:
            data_eff = p / (p + f)
            data_eff_pm = get_efficiency_error(p, f, p_pm, f_pm)

            w_y_data.append(data_eff)
            w_y_data_pm.append( data_eff_pm)
            w_x_data.append(xi)

            if p_mc + f_mc > 0:
                mc_eff = p_mc / (p_mc + f_mc)
                mc_eff_pm = get_efficiency_error(p_mc, f_mc, p_mc_pm, f_mc_pm)

                if mc_eff != 0:

                    sf_pm = getSFerrorbar(data_eff, mc_eff, data_eff_pm, mc_eff_pm) 
                    if sf_pm > 2:  # leave off huge points that distort scale?
                        continue
                    w_x_sf.append(xi)
                    w_y_sf.append(data_eff / mc_eff)
                    w_y_sf_pm.append(sf_pm )

    efficiency_loss_histogram_plotting(config_file, workdir, index, binname, x_data, y_data, y_data_pm, x_mc, y_mc, y_mc_pm, x_sf, y_sf, y_sf_pm, w_x_data, w_y_data, w_y_data_pm, w_x_sf, w_y_sf, w_y_sf_pm)



def efficiency_loss_histogram_from_minuitfit(config_file, workdir, index, binname, x_data_pass, signal_pass, uncert_signal_pass, x_data_fail, signal_fail, uncert_signal_fail, x_mc_pass, mc_pass, uncert_mc_pass, x_mc_fail, mc_fail, uncert_mc_fail):


    # The fit boundaries are NOT necessarily the same 
    # Find the overlapping region where both fits exist

    '''
    print(f"INFO For Index {index}")
    print(f"signal_pass:{signal_pass}")
    print(f"uncert_signal_pass:{uncert_signal_pass}")
    print(f"signal_fail:{signal_fail}")
    print(f"ucnert_signal_fail:{uncert_signal_fail}")
    print(f"mc_pass:{mc_pass}")
    print(f"uncert_mc_pass:{uncert_mc_pass}")
    print(f"mc_fail:{mc_fail}")
    print(f"uncert_mc_fail:{uncert_mc_fail}")
    '''






    x_data_pass_start, x_data_pass_end, x_data_fail_start, x_data_fail_end = overlapIndicesOfLists(x_data_pass, x_data_fail)


    x_data = []
    y_data = []
    y_data_pm = []
    y_mc = []
    y_mc_pm = []
    x_mc = []


    for i_pass, i_fail in zip(range(x_data_pass_start, x_data_pass_end, 1), range(x_data_fail_start, x_data_fail_end, 1) ):
        p = signal_pass[i_pass]
        f = signal_fail[i_fail]

        p_pm = uncert_signal_pass[i_pass]
        f_pm = uncert_signal_fail[i_fail]

        if p+f > 0:
            y_data.append(p/(p+f))
            x_data.append( x_data_pass[i_pass]  )
            y_data_pm.append(   get_efficiency_error(p, f, p_pm, f_pm))



    x_mc_pass_start, x_mc_pass_end, x_mc_fail_start, x_mc_fail_end = overlapIndicesOfLists(x_mc_pass, x_mc_fail)


    for i_pass, i_fail in zip(range(x_mc_pass_start, x_mc_pass_end, 1), range(x_mc_fail_start, x_mc_fail_end, 1) ):
        p = mc_pass[i_pass]
        f = mc_fail[i_fail]

        p_pm = uncert_mc_pass[i_pass]
        f_pm = uncert_mc_fail[i_fail]

        if p+f > 0:
            y_mc.append(p/(p+f))
            x_mc.append(x_mc_pass[i_pass])
            y_mc_pm.append(   get_efficiency_error(p, f, p_pm, f_pm))



    x_data_common = []
    y_data_common = []
    y_data_pm_common = []
    y_mc_common = []
    y_mc_pm_common = []
    x_mc_common = []


    for i in range(len(x_data)):
        if x_data[i] not in x_mc:
            continue

        mc_i = x_mc.index(x_data[i])

        x_data_common.append(x_data[i])
        y_data_common.append(y_data[i])
        y_data_pm_common.append(y_data_pm[i])
        x_mc_common.append(x_mc[mc_i])
        y_mc_common.append(y_mc[mc_i])
        y_mc_pm_common.append(y_mc_pm[mc_i])



    x_sf = []
    y_sf = []
    y_sf_pm = []

    for i in range(len(x_data_common)):

        if y_mc_common[i] > 0:
            x_sf.append(x_data_common[i])
            y_sf.append(y_data_common[i] / y_mc_common[i])
            y_sf_pm.append( getSFerrorbar(y_data_common[i], y_mc_common[i], y_data_pm_common[i], y_mc_pm_common[i])    )



    efficiency_loss_histogram_plotting(config_file, workdir, index, binname, x_data, y_data, y_data_pm, x_mc, y_mc, y_mc_pm, x_sf, y_sf, y_sf_pm, None, None, None, None, None, None)


def efficiency_loss_histogram_plotting(config_file, workdir, index, binname, x_data, y_data, y_data_pm, x_mc, y_mc, y_mc_pm, x_sf, y_sf, y_sf_pm, w_x_data, w_y_data, w_y_data_pm, w_x_sf, w_y_sf, w_y_sf_pm):

    # Create figure with two subplots, height ratio 2:1
    fig, (ax, ax_sf) = plt.subplots(
        2, 1,
        figsize=(8, 7),
        sharex=True,
        gridspec_kw={"height_ratios": [2, 1], "hspace": 0.05}
    )

    ax.set_title(binname)

    # Top plot
    ax.errorbar(x_mc, y_mc, yerr=cropErrorbars(y_mc_pm, y_mc), fmt='o', label="MC [FR]", color="orange")
    ax.errorbar(x_data, y_data, yerr=cropErrorbars(y_data_pm, y_data), fmt='o', label="Data [FR]", color = "green")

    if w_y_data is not None:
        ax.errorbar(w_x_data, w_y_data, yerr=w_y_data_pm, fmt='o', label="Data [from window]", color="navy")

    ax.set_ylabel("Efficiency")
    ax.legend()

    # FIXME FIXME FIXME ÍGY NEM JÓ, A SF ERRORBAR IS ASYMM HA ILYEN ASYMM IZÉ KELL.......

    # Bottom plot: SF = Data / MC
    ax_sf.errorbar(x_sf, y_sf, yerr=y_sf_pm, fmt='o', color="green")

    if w_y_sf is not None:
        ax_sf.errorbar(w_x_sf, w_y_sf, yerr=w_y_sf_pm, fmt='o', color="navy")

    ax_sf.set_xlabel("mass [GeV]")
    ax_sf.set_ylabel("SF")

    # reference line at SF = 1
    ax_sf.axhline(1.0, color="gray", linestyle="--", linewidth=1)

    fig.tight_layout()

    if not os.path.exists(workdir + "/commonplots/binned_effi/"):
        os.makedirs(workdir + "/commonplots/binned_effi/")

    fig.savefig(
        workdir + f"/commonplots/binned_effi/bin{index}.png",
        bbox_inches="tight",
        pad_inches=0.02,
    )

    plt.close()

    ############################################################
    # Rebinned scatter plot


    if not os.path.exists(workdir+f"/energy_loss_profiles/bin{index}"):
        # No energy loss TProfiles exist, just skip
        return 

    
    proot_file = ROOT.TFile.Open(workdir+f"/energy_loss_profiles/bin{index}/profile_{config_file["collection"].split(" ", 1)[0]}_bin{index}.root")

    rkeys = proot_file.GetListOfKeys()
    if rkeys.GetEntries() != 1:
        raise RuntimeError(f"Expected exactly 1 object in TProfile root file, found {rkeys.GetEntries()}")

    rkey = rkeys.At(0)
    profile = rkey.ReadObj()

    # Check that it was found
    if not profile:
        raise RuntimeError("Could not find TProfile in ROOT file")
    if not profile.InheritsFrom("TProfile"):
        raise TypeError(f"Expected TProfile, got {profile.ClassName()}")


    rbsf_x = []
    rbsf_x_pm = []
    rbsf_y = []
    rbsf_y_pm = []

    for x, y, y_pm in zip(x_sf, y_sf, y_sf_pm):
        bin_number = profile.FindBin(x)
        energyloss = profile.GetBinContent(bin_number)
        energyloss_err = profile.GetBinError(bin_number)

        # FIXME band aid?
        if not energyloss_err > 0:
            continue 


        rbsf_x.append(energyloss)
        rbsf_x_pm.append(energyloss_err)
        rbsf_y.append(y)
        rbsf_y_pm.append(y_pm)

        #print(bin_number, x, y, y_pm, energyloss, energyloss_err)


    fig, ax = plt.subplots(figsize=(10, 6))

    ax.errorbar(
        rbsf_x, rbsf_y,
        xerr=rbsf_x_pm,
        yerr=rbsf_y_pm,
        fmt='o',
        color="black" if w_x_data is None else "green",
        ecolor='black' if w_x_data is None else "green",
        capsize=3,
        markersize=5,
        linestyle='none'
    )

    ax.set_xlabel(r'$\left\langle \frac{(\text{Reco probe eT}) - (\text{Gen probe eT})}{(\text{Reco probe eT})}     \right\rangle$')
    ax.set_ylabel(r'Scale factor (Data eff/MC eff)')
    

    if w_x_data is not None:


        rbsf_x = []
        rbsf_x_pm = []
        rbsf_y = []
        rbsf_y_pm = []

        for x, y, y_pm in zip(w_x_sf, w_y_sf, w_y_sf_pm):
            bin_number = profile.FindBin(x)
            energyloss = profile.GetBinContent(bin_number)
            energyloss_err = profile.GetBinError(bin_number)

            # FIXME band aid?
            if not energyloss_err > 0:
                continue 

            rbsf_x.append(energyloss)
            rbsf_x_pm.append(energyloss_err)
            rbsf_y.append(y)
            rbsf_y_pm.append(y_pm)

        ax.errorbar(
            rbsf_x, rbsf_y,
            xerr=rbsf_x_pm,
            yerr=rbsf_y_pm,
            fmt='o',
            color="navy",
            ecolor="navy",
            capsize=3,
            markersize=5,
            linestyle='none'
        )


    # reference line at SF = 1
    ax.axhline(1.0, color="gray", linestyle="--", linewidth=1)
    ax.axvline(0.0, color="gray", linestyle="--", linewidth=1)

    fig.savefig(
        workdir + f"/commonplots/binned_effi/rb_bin{index}.png",
        bbox_inches="tight",
        pad_inches=0.02,
    )
    plt.close()


    ############################################################
    # Combined plot: efficiency + SF on top, rebinned SF below
    fig = plt.figure(figsize=(10, 15))

    # Outer grid: top group + bottom plot
    gs_outer = fig.add_gridspec(
        2, 1,
        height_ratios=[3, 2],
        hspace=0.1
    )

    # Inner grid: efficiency + SF with almost no spacing
    gs_top = gs_outer[0].subgridspec(
        2, 1,
        height_ratios=[2, 1],
        hspace=0.05
    )

    ax = fig.add_subplot(gs_top[0])
    ax_sf = fig.add_subplot(gs_top[1], sharex=ax)

    # Bottom: rebinned SF
    ax_rb = fig.add_subplot(gs_outer[1])

    # Top: efficiency
    ax.set_title(binname)

    ax.errorbar(
        x_mc, y_mc,
        yerr=cropErrorbars(y_mc_pm, y_mc),
        fmt='o', label="MC [FR]", color="orange"
    )
    ax.errorbar(
        x_data, y_data,
        yerr=cropErrorbars(y_data_pm, y_data),
        fmt='o', label="Data [FR]", color="green"
    )

    if w_y_data is not None:
        ax.errorbar(
            w_x_data, w_y_data,
            yerr=w_y_data_pm,
            fmt='o', label="Data [from window]", color="navy"
        )

    ax.set_ylabel("Efficiency")
    ax.legend()
    ax.tick_params(labelbottom=False)

    # Middle: SF


    '''
    ax_sf.errorbar(
        x_sf, y_sf, yerr=y_sf_pm,
        fmt='o', color="green"
    )
    '''
    cmap = LinearSegmentedColormap.from_list(
        "green_gradient", ["lightgreen", "darkgreen"]
    )
    for i in range(len(x_sf)):
        color = cmap(i / (len(x_sf) - 1))

        ax_sf.errorbar(
            x_sf[i], y_sf[i],
            yerr=y_sf_pm[i],
            fmt='o',          # Points only
            color=color,      # Marker and error bar color
            capsize=3,
            markersize=5,
            linestyle='none'  # No connecting lines
        )




    if w_y_sf is not None:
        cmap2 = LinearSegmentedColormap.from_list(
            "blue_gradient", ["lightblue", "darkblue"]
        )

        for i in range(len(w_x_sf)):
            color = cmap2(i / (len(w_x_sf) - 1))

            ax_sf.errorbar(
                w_x_sf[i], w_y_sf[i],
                yerr=w_y_sf_pm[i],
                fmt='o',          # Points only
                color=color,      # Marker and error bar color
                capsize=3,
                markersize=5,
                linestyle='none'  # No connecting lines
            )


    ax_sf.set_xlabel("mass [GeV]")
    ax_sf.set_ylabel("SF")
    ax_sf.axhline(1.0, color="gray", linestyle="--", linewidth=1)

    # Bottom: rebinned SF vs energy loss
    # Recalculate arrays because the existing ones may have been
    # overwritten by the window-based values.

    def get_rebinned_points(xs, ys, yerrs):
        rb_x, rb_xerr, rb_y, rb_yerr = [], [], [], []

        for x, y, y_pm in zip(xs, ys, yerrs):
            bin_number = profile.FindBin(x)
            energyloss = profile.GetBinContent(bin_number)
            energyloss_err = profile.GetBinError(bin_number)

            if not energyloss_err > 0:
                continue

            rb_x.append(energyloss)
            rb_xerr.append(energyloss_err)
            rb_y.append(y)
            rb_yerr.append(y_pm)

        return rb_x, rb_xerr, rb_y, rb_yerr

    rb_x, rb_xerr, rb_y, rb_yerr = get_rebinned_points(
        x_sf, y_sf, y_sf_pm
    )




    for i in range(len(rb_x)):
        color = cmap(i / (len(rb_x) - 1))

        ax_rb.errorbar(
            rb_x[i], rb_y[i],
            xerr=rb_xerr[i],
            yerr=rb_yerr[i],
            fmt='o',          # Points only
            color=color,      # Marker and error bar color
            capsize=3,
            markersize=5,
            linestyle='none'  # No connecting lines
        )


    if w_x_data is not None:
        rb_x, rb_xerr, rb_y, rb_yerr = get_rebinned_points(
            w_x_sf, w_y_sf, w_y_sf_pm
        )

        for i in range(len(rb_x)):
            color = cmap2(i / (len(rb_x) - 1))

            ax_rb.errorbar(
                rb_x[i], rb_y[i],
                xerr=rb_xerr[i],
                yerr=rb_yerr[i],
                fmt='o',          # Points only
                color=color,      # Marker and error bar color
                capsize=3,
                markersize=5,
                linestyle='none'  # No connecting lines
            )

    ax_rb.set_xlabel(
        r'$\left\langle \frac{(\mathrm{Reco\ probe\ eT})'
        r' - (\mathrm{Gen\ probe\ eT})}'
        r'{(\mathrm{Reco\ probe\ eT})} \right\rangle$'
    )
    ax_rb.set_ylabel("Scale factor (Data eff/MC eff)")


    # reference line at SF = 1
    ax_rb.axhline(1.0, color="gray", linestyle="--", linewidth=1)
    ax_rb.axvline(0.0, color="gray", linestyle="--", linewidth=1)

    fig.savefig(
        workdir + f"/commonplots/binned_effi/combined_bin{index}.png",
        bbox_inches="tight",
        pad_inches=0.02,
    )

    plt.close(fig)











def get_efficiency(P, F):
    return P/(P+F)


def get_efficiency_error(P, F, DP, DF):
    if P + F == 0:
        print("[WARN]: P+F Zero division in get_efficiency_error, returning 0.001")
        return 0.001

    return math.sqrt(  (  (  F/((P+F)**2)   )*DP   )**2    +    (  (  -P/((P+F)**2)     )*DF   )**2   )


def write_mparray(name, mparray, outfile):
    outfile.write(name+": [")
    for i in range(len(mparray) - 1):
        outfile.write(str(mparray[i])+", ")
    outfile.write(str(mparray[-1]) + "]\n")




def read_effi(file):
    data = []
    with open(file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue

            row = list(map(float, line.split()))

            # Skip rows with first two columns = 1.44400 and 1.56600
            if row[0] == 1.44400 and row[1] == 1.56600:
                continue

            data.append(row)

    data = np.array(data)

    return data[:, 4], data[:, 5]


def main():

    options = get_parser().parse_args()



    with open(options.config_file, 'r') as file:
        config_file = yaml.safe_load(file)


    for i in range(len(options.config_file)-1, 0, -1):
        if options.config_file[i] == '/':
            workdir = options.config_file[:i]+"/"
            break



    effi = None
    effierr = None
    if options.effi:
        effi, effierr = read_effi(options.effi)
        print("Info: egammaEffi file supplied for data template fit results - data will NOT be fit or plotted")
    else:
        if "effi" in config_file.keys():
            effi, effierr = read_effi(workdir+"/"+config_file["effi"])
            print("Info: egammaEffi file supplied in config file for data template fit results - data will NOT be fit or plotted")






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

        exec(r'binnames.append(r"Bin%s, $0 < \| \eta \| < 1.4442$,  $%s < p_{T} < %s$, PASS") ' % (str(2*i), ptbins[i][:-3], ptbins[i+1][:-3]))
        exec(r'binnames.append(r"Bin%s, $0 < \| \eta \| < 1.4442$,  $%s < p_{T} < %s$, FAIL") ' % (str(2*i), ptbins[i][:-3], ptbins[i+1][:-3]))
        exec(r'binnames.append(r"Bin%s, $1.566 < \| \eta \| < 2.5$,  $%s < p_{T} < %s$, PASS") ' % (str(2*i + 1), ptbins[i][:-3], ptbins[i+1][:-3]))
        exec(r'binnames.append(r"Bin%s, $1.566 < \| \eta \| < 2.5$,  $%s < p_{T} < %s$, FAIL") ' % (str(2*i + 1), ptbins[i][:-3], ptbins[i+1][:-3]))

    if len(binindex_range_to_fit) == 1 :
        binindex_range_to_fit.append(len(bins))


    ptbins_center_fit = ptbins_center[(binindex_range_to_fit[0]//4):(binindex_range_to_fit[1]//4)]
    binw_fit = binw[(binindex_range_to_fit[0]//4):(binindex_range_to_fit[1]//4)]

    if options.threads is None:
        threads = multiprocessing.cpu_count()

    elif options.threads == -1:
        threads = multiprocessing.cpu_count()
    elif options.threads < 1:
        print ("Ain't the sharpest tool, are we?")
        quit()

    else:
        threads = options.threads


    if options.bins is not None:
        try:
            fitted_figures = [int(x) for x in options.bins.split(",")]
        except Exception:
            print("Insensible bins")
            quit()
    else:
        fitted_figures = "all"

    '''  old threading
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
    '''



    # get number of processed bins
    skipped_lowpt_bins = 0
    num_fitted_pt_bins = 0
    for pt in config_file['ptbins_num']:
        if pt < config_file['fit_in_pt_range'][0]:
            skipped_lowpt_bins += 1
        elif pt < config_file['fit_in_pt_range'][1]:
            num_fitted_pt_bins += 1
        else:
            break


    bins_to_process = range(skipped_lowpt_bins*4, (skipped_lowpt_bins+num_fitted_pt_bins)*4 )


    indexlists = [[] for _ in range(threads)]



    # This used to go until len(bins_to_process). 
    # The problem is if the have skipped lowpt bins, then this gets shorter than needed
    # Should work with bins_to_process[-1]
    # After all, we go until the last bin, and here we include the skipped lowpt ones as well to just get their distribution plot. Only the fit is skipped.
    _worker = 0
    for i in range(0, bins_to_process[-1], 2):    
        indexlists[_worker].append(i)
        indexlists[_worker].append(i+1)
        _worker += 1
        if _worker >= len(indexlists): _worker = 0



    def allocate_multiprocessing_array():
        return multiprocessing.Array('d',range((bins_to_process[-1] + 1) // 2))


    mutex = multiprocessing.Lock()

    class SharedEffArraysHolderClass:
        def __init__(self):
            self.shared_data_cutcount_eff = allocate_multiprocessing_array() 
            self.shared_data_cutcount_eff_total_err = allocate_multiprocessing_array() 
            self.shared_data_cutcount_eff_stat_err = allocate_multiprocessing_array() 
            self.shared_data_cutcount_eff_syst_err = allocate_multiprocessing_array() 
            self.shared_data_cbfit_eff = allocate_multiprocessing_array() 
            self.shared_data_cbfit_eff_total_err = allocate_multiprocessing_array() 
            self.shared_data_cbfit_eff_stat_err = allocate_multiprocessing_array() 
            self.shared_data_cbfit_eff_syst_err = allocate_multiprocessing_array() 
            self.shared_data_minuit_eff = allocate_multiprocessing_array() 
            self.shared_data_minuit_eff_total_err = allocate_multiprocessing_array() 
            self.shared_data_minuit_eff_stat_err = allocate_multiprocessing_array() 
            self.shared_data_minuit_eff_syst_err = allocate_multiprocessing_array() 
            self.shared_data_minuit_fr_eff = allocate_multiprocessing_array() 
            self.shared_data_minuit_fr_eff_stat_err = allocate_multiprocessing_array() 
            self.shared_mc_cutcount_eff = allocate_multiprocessing_array() 
            self.shared_mc_cutcount_eff_total_err = allocate_multiprocessing_array() 
            self.shared_mc_cutcount_eff_stat_err = allocate_multiprocessing_array() 
            self.shared_mc_cutcount_eff_syst_err = allocate_multiprocessing_array() 
            self.shared_mc_cutcount_fr_eff = allocate_multiprocessing_array() 
            self.shared_mc_cutcount_fr_eff_total_err = allocate_multiprocessing_array()
            self.shared_zeros = multiprocessing.Array('d',range(num_fitted_pt_bins * 4))


    sea = SharedEffArraysHolderClass()

    print( "Opened all input resources")

    proc = {}
    for i in range(threads):
        proc[i] = multiprocessing.Process(target=TnPEfficiency, args=(mutex, effi, effierr, sea, 
                       datafile, mcfile, mcfile_alt, bins, binnames, indexlists[i], binindex_range_to_fit, fitted_figures, workdir, config_file, skipped_lowpt_bins))
        proc[i].start()

    for i in range(threads):
        proc[i].join()
    

    print( "Making efficiency plots")


    if config_file['resonance'] == "Z":
        if effi is not None:
            # effi might have more bins (up to 500 GeV)
            data_imported_eff = effi[:len(sea.shared_mc_cutcount_eff)]
            data_imported_eff_error = effierr[:len(sea.shared_mc_cutcount_eff)]


            makeEfficiencyPlot(workdir+"efficiency-template-cutcount", data_imported_eff, data_imported_eff_error, sea.shared_mc_cutcount_eff, sea.shared_mc_cutcount_eff_total_err, ptbins_center_fit, binw_fit, skipped_lowpt_bins, config_file, mclabel = "(grey range)", datalabel="(templ)")

    # For JPsi
    else:
        if config_file['bUseMinuit']:
            makeEfficiencyPlot(workdir+"efficiency-minuit-mcfr", sea.shared_data_minuit_eff, sea.shared_data_minuit_eff_total_err, sea.shared_mc_cutcount_fr_eff, sea.shared_mc_cutcount_fr_eff_total_err, ptbins_center_fit, binw_fit, skipped_lowpt_bins, config_file, mclabel = "(cc fr)", datalabel="(minuit)")
        else:
            makeEfficiencyPlot(workdir+"efficiency-cutcount", sea.shared_data_cutcount_eff, sea.shared_data_cutcount_eff_total_err, sea.shared_mc_cutcount_eff, sea.shared_mc_cutcount_eff_total_err, ptbins_center_fit, binw_fit, skipped_lowpt_bins, config_file)
            makeEfficiencyPlot(workdir+"efficiency-cutcountfr", sea.shared_data_cutcount_eff, sea.shared_data_cutcount_eff_total_err, sea.shared_mc_cutcount_fr_eff, sea.shared_mc_cutcount_fr_eff_total_err, ptbins_center_fit, binw_fit, skipped_lowpt_bins, config_file)
            makeEfficiencyPlot(workdir+"efficiency-crystalball-cc", sea.shared_data_cbfit_eff, sea.shared_data_cbfit_eff_total_err, sea.shared_mc_cutcount_eff, sea.shared_mc_cutcount_eff_total_err, ptbins_center_fit, binw_fit, skipped_lowpt_bins, config_file)
            makeEfficiencyPlot(workdir+"efficiency-crystalball-ccfr", sea.shared_data_cbfit_eff, sea.shared_data_cbfit_eff_total_err, sea.shared_mc_cutcount_fr_eff, sea.shared_mc_cutcount_fr_eff_total_err, ptbins_center_fit, binw_fit, skipped_lowpt_bins, config_file)


    with open(workdir+"/results_"+config_file["resonance"]+"_chev"+str(config_file["chevorder"])+".yaml", "w") as out:
        if effi is not None:
            write_mparray("data_imported_eff", data_imported_eff, out)
            write_mparray("data_imported_eff_total_err", data_imported_eff_error, out)


        for name, value in vars(sea).items():
            write_mparray(name, value, out)

        write_mparray("ptbins_center_fit", ptbins_center_fit, out)
        write_mparray("binw_fit", binw_fit, out)
        out.write("skipped_lowpt_bins: "+str(skipped_lowpt_bins))



    # Gridding
    print("Gridding...")
    if config_file["resonance"] == "Z":
        from gridderWindowFr import create_combined_images
        create_combined_images(workdir+"/dataplots/templatefit/fr/", workdir+"/dataplots/templatefit/", workdir+"/dataplots/templatefit/gridded/")
    else:
        from gridderJpsi import create_combined_images
        create_combined_images(workdir+f"/dataplots/fits/chev{config_file["chevorder"]}", workdir+"/mcplots/fits/", workdir+"commonplots/")



def getSFerrorbar(data, mc, data_err, mc_err, rho=0.5):
    if data == 0:
        print("[WARN] Zero data division in getSFerrorbar, returning 0.001")
        return 0.001
    if mc == 0:
        print("[WARN] Zero mc division in getSFerrorbar, returning 0.001")
        return 0.001

    sf = data / mc
    data_relerr = data_err/data
    mc_relerr = mc_err/mc
    rel_err = np.sqrt(  data_relerr**2  + mc_relerr**2   - 2*rho*data_relerr*mc_relerr  )
    return rel_err * sf




def makeEfficiencyPlot(filename, data_efficiencies, data_eff_errorbars, mc_efficiencies, mc_eff_errorbars, ptbins_center_fit, binw_fit, skipped_lowpt_bins, config_file, mclabel = "", datalabel = ""):

    rho = 0.5

    binstart_barrel = 0 + skipped_lowpt_bins*2
    binstart_endcap = 1 + skipped_lowpt_bins*2

    data_eff_barrel = data_efficiencies[binstart_barrel::2]
    data_eff_barrel_error_total = data_eff_errorbars[binstart_barrel::2]
    mc_eff_barrel = mc_efficiencies[binstart_barrel::2]
    mc_eff_barrel_error_total = mc_eff_errorbars[binstart_barrel::2]
    data_eff_endcap = data_efficiencies[binstart_endcap::2]
    data_eff_endcap_error_total = data_eff_errorbars[binstart_endcap::2]
    mc_eff_endcap = mc_efficiencies[binstart_endcap::2]
    mc_eff_endcap_error_total = mc_eff_errorbars[binstart_endcap::2]

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
    for i in range(len(data_eff_barrel)):


        if data_eff_barrel[i] > 0 and mc_eff_barrel[i] > 0:

            if data_eff_barrel_error_total[i] > 0 and mc_eff_barrel_error_total[i]:
                scale_factor_barrel.append(data_eff_barrel[i]/mc_eff_barrel[i])
                scale_factor_x_barrel.append(ptbins_center_fit[pt_bin_number])
                scale_factor_x_barrel_binw.append(binw_fit[pt_bin_number])
                scale_factor_barrel_error.append( getSFerrorbar(data_eff_barrel[i], mc_eff_barrel[i], data_eff_barrel_error_total[i],  mc_eff_barrel_error_total[i], rho)   )
            else:
                scale_factor_barrel_black.append(data_eff_barrel[i]/mc_eff_barrel[i])
                scale_factor_barrel_black_x.append(ptbins_center_fit[pt_bin_number])
                scale_factor_barrel_black_x_binw.append(binw_fit[pt_bin_number])

        if data_eff_endcap[i] > 0 and mc_eff_endcap[i] > 0:

            if data_eff_endcap_error_total[i] > 0 and mc_eff_endcap_error_total[i] > 0:
                scale_factor_endcap.append(data_eff_endcap[i]/mc_eff_endcap[i])
                scale_factor_x_endcap.append(ptbins_center_fit[pt_bin_number])
                scale_factor_x_endcap_binw.append(binw_fit[pt_bin_number])
                scale_factor_endcap_error.append( getSFerrorbar(data_eff_endcap[i], mc_eff_endcap[i], data_eff_endcap_error_total[i],  mc_eff_endcap_error_total[i], rho)   )
            else:
                scale_factor_endcap_black.append(data_eff_endcap[i]/mc_eff_endcap[i])
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

    ax_barrel.set_xlabel(r"$p_{T}$ [GeV]", fontsize=14)
    ax_endcap.set_xlabel(r"$p_{T}$ [GeV]", fontsize=14)
    ax_barrel.set_ylabel("Efficiency", fontsize=15)
    ax_barrel.set_title(r"$|\eta | < 1.4442$", fontsize=15)
    ax_endcap.set_title(r"$1.566 < |\eta | < 2.50$", fontsize=15)
    bottom_ax_barrel.set_xlabel(r"$p_{T}$ [GeV]", fontsize=14)
    bottom_ax_endcap.set_xlabel(r"$p_{T}$ [GeV]", fontsize=14)
    bottom_ax_barrel.set_ylabel("Data/MC", fontsize=15)

    #print len(ptbins_center_fit), len(data_eff_barrel), len(data_eff_barrel_error), len(binw_fit)
    #print ptbins_center_fit, binw_fit

    #ax_barrel.plot(ptbins_center_fit, data_eff_barrel, "bs")
    ax_barrel.errorbar(ptbins_center_fit, data_eff_barrel, yerr=data_eff_barrel_error_total, xerr=binw_fit, capsize=3, fmt=".", color="blue", ecolor = "blue", label="Data"+datalabel)
    #ax_barrel.plot(ptbins_center_fit, mc_eff_barrel, "rs")
    ax_barrel.errorbar(ptbins_center_fit, mc_eff_barrel, yerr=mc_eff_barrel_error_total, xerr=binw_fit, capsize=3, fmt=".", color="crimson", ecolor = "crimson", label="MC"+mclabel)

    bottom_ax_barrel.errorbar(scale_factor_x_barrel, scale_factor_barrel, yerr=scale_factor_barrel_error , xerr=scale_factor_x_barrel_binw, capsize=3, fmt=".", color="black", ecolor = "black")
    bottom_ax_barrel.errorbar(scale_factor_barrel_black_x, scale_factor_barrel_black, xerr=scale_factor_barrel_black_x_binw, fmt="x", color="black", ecolor="black")


    #ax_endcap.plot(ptbins_center_fit, data_eff_endcap, "bs")
    ax_endcap.errorbar(ptbins_center_fit, data_eff_endcap, yerr=data_eff_endcap_error_total, xerr=binw_fit, capsize=3, fmt=".", color="blue", ecolor = "blue", label="Data"+datalabel)
    #ax_endcap.plot(ptbins_center_fit, mc_eff_endcap, "rs")
    ax_endcap.errorbar(ptbins_center_fit, mc_eff_endcap, yerr=mc_eff_endcap_error_total, xerr=binw_fit, capsize=3, fmt=".", color="crimson", ecolor = "crimson", label="MC"+mclabel)

    bottom_ax_endcap.errorbar(scale_factor_x_endcap, scale_factor_endcap, yerr=scale_factor_endcap_error, xerr=scale_factor_x_endcap_binw, capsize=3, fmt=".", color="black", ecolor = "black")
    bottom_ax_endcap.errorbar(scale_factor_endcap_black_x, scale_factor_endcap_black, xerr=scale_factor_endcap_black_x_binw, fmt="x", color="black", ecolor="black")

    legend = ax_barrel.legend(loc="best", fontsize=17, title=config_file["resonance"]+" "+str(config_file["year"]))
    legend.get_title().set_fontsize(17)

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

    fig.suptitle(config_file["collection"]+", WP: "+config_file["workingpoint"], fontsize=16)

    plt.savefig(filename)










def cropErrorbars(errorbars, values, lim_lower=0, lim_upper = 1):
    lower = []
    upper = []

    for eb, val in zip(errorbars, values):
        if not (lim_lower <= val <= lim_upper):
            raise Exception("cropErrorbars: value is not in limits")
        
        lower.append(  eb if val - eb > lim_lower else val - lim_lower )
        upper.append(  eb if val + eb < lim_upper else lim_upper - val )

    return (lower, upper)





if __name__ == "__main__":
    main()




