import numpy as np
import os
import numpy as np
import math
import matplotlib.pyplot as plt
import yaml



def get_parser():
    import argparse
    argParser = argparse.ArgumentParser(description = "Argument parser")
    argParser.add_argument('jpsidata',help='config_file')
    argParser.add_argument('zdata',help='config_file')
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

    makeEfficiencyPlot("SF/rawint_effs.png", 
                        j["shared_data_raw_integrals"],
                        j["shared_data_raw_integral_sigmas"],
                        j["shared_data_raw_integral_sigmas_bksyst"],
                        j["shared_data_raw_integral_sigmas_rangesyst"],
                        j["shared_mc_raw_integrals"],
                        j["shared_mc_raw_integral_sigmas"],
                        j["shared_mc_raw_integral_sigmas_bksyst"],
                        j["shared_mc_raw_integral_sigmas_rangesyst"],
                        j["ptbins_center_fit"],
                        j["binw_fit"],
                        z["shared_data_raw_integrals"],
                        z["shared_data_raw_integral_sigmas"],
                        z["shared_data_raw_integral_sigmas_bksyst"],
                        z["shared_data_raw_integral_sigmas_rangesyst"],
                        z["shared_mc_raw_integrals"],
                        z["shared_mc_raw_integral_sigmas"],
                        z["shared_mc_raw_integral_sigmas_bksyst"],
                        z["shared_mc_raw_integral_sigmas_rangesyst"],
                        z["ptbins_center_fit"],
                        z["binw_fit"]
                        )


    makeEfficiencyPlot("SF/cbint_jpsim_effs.png", 
                        j["shared_data_integrals_cb"],
                        j["shared_data_integral_cb_sigmas"],
                        j["shared_data_integral_cb_sigmas_bksyst"],
                        j["shared_data_integral_cb_sigmas_rangesyst"],
                        j["shared_mc_raw_integrals"],
                        j["shared_mc_raw_integral_sigmas"],
                        j["shared_mc_raw_integral_sigmas_bksyst"],
                        j["shared_mc_raw_integral_sigmas_rangesyst"],
                        j["ptbins_center_fit"],
                        j["binw_fit"],
                        z["shared_data_integrals_cb"],
                        z["shared_data_integral_cb_sigmas"],
                        z["shared_data_integral_cb_sigmas_bksyst"],
                        z["shared_data_integral_cb_sigmas_rangesyst"],
                        z["shared_mc_integrals_cb"],
                        z["shared_mc_integral_cb_sigmas"],
                        z["shared_mc_integral_cb_sigmas_bksyst"],
                        z["shared_mc_integral_cb_sigmas_rangesyst"],
                        z["ptbins_center_fit"],
                        z["binw_fit"]
                        )





def get_efficiency_error(P, F, DP, DF):
    return math.sqrt(  (  (  F/((P+F)**2)   )*DP   )**2    +    (  (  -P/((P+F)**2)     )*DF   )**2   )


def makeEfficiencyPlot(filename, jpsi_data_integrals, jpsi_data_integral_sigmas, jpsi_data_integral_sigmas_bksyst, jpsi_data_integral_sigmas_rangesyst, jpsi_mc_integrals, jpsi_mc_integral_sigmas, jpsi_mc_integral_sigmas_bksyst, jpsi_mc_integral_sigmas_rangesyst, jpsi_ptbins_center_fit, jpsi_binw_fit, z_data_integrals, z_data_integral_sigmas, z_data_integral_sigmas_bksyst, z_data_integral_sigmas_rangesyst, z_mc_integrals, z_mc_integral_sigmas, z_mc_integral_sigmas_bksyst, z_mc_integral_sigmas_rangesyst, z_ptbins_center_fit, z_binw_fit):

    jpsi_data_eff_barrel = []
    jpsi_data_eff_barrel_error_total = []
    jpsi_mc_eff_barrel = []
    jpsi_mc_eff_barrel_error_total = []
    jpsi_data_eff_endcap = []
    jpsi_data_eff_endcap_error_total = []
    jpsi_mc_eff_endcap = []
    jpsi_mc_eff_endcap_error_total = []

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
    for i in range(0,len(jpsi_data_integrals),4):
        
        if jpsi_data_integrals[i] and jpsi_data_integrals[i+1]:
            data_barrel_pass = jpsi_data_integrals[i]
            data_barrel_pass_error = jpsi_data_integral_sigmas[i] if jpsi_data_integral_sigmas[i] and jpsi_data_integral_sigmas[i]/jpsi_data_integrals[i] < 1 else False
            data_barrel_pass_error_bksyst = abs(jpsi_data_integral_sigmas_bksyst[i]) 
            data_barrel_pass_error_rangesyst = abs(jpsi_data_integral_sigmas_rangesyst[i])

            data_barrel_fail = jpsi_data_integrals[i+1] 
            data_barrel_fail_error = jpsi_data_integral_sigmas[i+1] if jpsi_data_integral_sigmas[i+1] and jpsi_data_integral_sigmas[i+1]/jpsi_data_integrals[i+1] < 1 else False
            data_barrel_fail_error_bksyst = abs(jpsi_data_integral_sigmas_bksyst[i+1]) 
            data_barrel_fail_error_rangesyst = abs(jpsi_data_integral_sigmas_rangesyst[i+1])

            jpsi_data_eff_barrel.append(data_barrel_pass / (data_barrel_pass + data_barrel_fail))
            jpsi_data_eff_barrel_error_total.append( get_efficiency_error(data_barrel_pass, data_barrel_fail, qsum(data_barrel_pass_error, data_barrel_pass_error_bksyst, data_barrel_pass_error_rangesyst), qsum(data_barrel_fail_error, data_barrel_fail_error_bksyst,data_barrel_fail_error_rangesyst)) if data_barrel_pass_error and data_barrel_fail_error else 0    )

        else:
            jpsi_data_eff_barrel.append(-1)
            jpsi_data_eff_barrel_error_total.append( 0.001 )

        if jpsi_data_integrals[i+2] and jpsi_data_integrals[i+3]:
            data_endcap_pass = jpsi_data_integrals[i+2]
            data_endcap_pass_error = jpsi_data_integral_sigmas[i+2] if jpsi_data_integral_sigmas[i+2]  and jpsi_data_integral_sigmas[i+2]/jpsi_data_integrals[i+2] < 1 else False
            data_endcap_pass_error_bksyst = abs(jpsi_data_integral_sigmas_bksyst[i+2]) 
            data_endcap_pass_error_rangesyst = abs(jpsi_data_integral_sigmas_rangesyst[i+2])
            data_endcap_fail = jpsi_data_integrals[i+3]
            data_endcap_fail_error = jpsi_data_integral_sigmas[i+3] if jpsi_data_integral_sigmas[i+3]  and jpsi_data_integral_sigmas[i+3]/jpsi_data_integrals[i+3] < 1 else False
            data_endcap_fail_error_bksyst = abs(jpsi_data_integral_sigmas_bksyst[i+3]) 
            data_endcap_fail_error_rangesyst = abs(jpsi_data_integral_sigmas_rangesyst[i+3])

            jpsi_data_eff_endcap.append(data_endcap_pass / (data_endcap_pass + data_endcap_fail))
            jpsi_data_eff_endcap_error_total.append( get_efficiency_error(data_endcap_pass, data_endcap_fail, qsum(data_endcap_pass_error, data_endcap_pass_error_bksyst, data_endcap_pass_error_rangesyst), qsum(data_endcap_fail_error, data_endcap_fail_error_bksyst,data_endcap_fail_error_rangesyst)) if data_endcap_pass_error and data_endcap_fail_error else 0  )
        
        else:
            jpsi_data_eff_endcap.append(-1)
            jpsi_data_eff_endcap_error_total.append( 0.001 )




        if jpsi_mc_integrals[i] and jpsi_mc_integrals[i+1]:
            mc_barrel_pass = jpsi_mc_integrals[i]
            mc_barrel_pass_error = jpsi_mc_integral_sigmas[i] if jpsi_mc_integral_sigmas[i] and jpsi_mc_integral_sigmas[i]/jpsi_mc_integrals[i] < 1 else False
            mc_barrel_pass_error_bksyst = abs(jpsi_mc_integral_sigmas_bksyst[i]) 
            mc_barrel_pass_error_rangesyst = abs(jpsi_mc_integral_sigmas_rangesyst[i])
            mc_barrel_fail = jpsi_mc_integrals[i+1]
            mc_barrel_fail_error = jpsi_mc_integral_sigmas[i+1] if jpsi_mc_integral_sigmas[i+1]  and jpsi_mc_integral_sigmas[i+1]/jpsi_mc_integrals[i+1] < 1 else False
            mc_barrel_fail_error_bksyst = abs(jpsi_mc_integral_sigmas_bksyst[i+1]) 
            mc_barrel_fail_error_rangesyst = abs(jpsi_mc_integral_sigmas_rangesyst[i+1])

            jpsi_mc_eff_barrel.append(mc_barrel_pass / (mc_barrel_pass + mc_barrel_fail))
            jpsi_mc_eff_barrel_error_total.append( get_efficiency_error(mc_barrel_pass, mc_barrel_fail, qsum(mc_barrel_pass_error, mc_barrel_pass_error_bksyst, mc_barrel_pass_error_rangesyst), qsum(mc_barrel_fail_error, mc_barrel_fail_error_bksyst,mc_barrel_fail_error_rangesyst)) if mc_barrel_pass_error and mc_barrel_fail_error else 0 )
        else:
            jpsi_mc_eff_barrel.append(-1)
            jpsi_mc_eff_barrel_error_total.append( 0.001 )


        if jpsi_mc_integrals[i+2] and jpsi_mc_integrals[i+3]:
            mc_endcap_pass = jpsi_mc_integrals[i+2]
            mc_endcap_pass_error = jpsi_mc_integral_sigmas[i+2] if jpsi_mc_integral_sigmas[i+2] and jpsi_mc_integral_sigmas[i+2]/jpsi_mc_integrals[i+2] < 1 else False
            mc_endcap_pass_error_bksyst = abs(jpsi_mc_integral_sigmas_bksyst[i+2]) 
            mc_endcap_pass_error_rangesyst = abs(jpsi_mc_integral_sigmas_rangesyst[i+2])
            mc_endcap_fail = jpsi_mc_integrals[i+3]
            mc_endcap_fail_error = jpsi_mc_integral_sigmas[i+3] if jpsi_mc_integral_sigmas[i+3] and jpsi_mc_integral_sigmas[i+3]/jpsi_mc_integrals[i+3] < 1 else False
            mc_endcap_fail_error_bksyst = abs(jpsi_mc_integral_sigmas_bksyst[i+3])
            mc_endcap_fail_error_rangesyst = abs(jpsi_mc_integral_sigmas_rangesyst[i+3])

            jpsi_mc_eff_endcap.append(mc_endcap_pass / (mc_endcap_pass + mc_endcap_fail))
            jpsi_mc_eff_endcap_error_total.append( get_efficiency_error(mc_endcap_pass, mc_endcap_fail, qsum(mc_endcap_pass_error, mc_endcap_pass_error_bksyst,mc_endcap_pass_error_rangesyst), qsum(mc_endcap_fail_error, mc_endcap_fail_error_bksyst,mc_endcap_fail_error_rangesyst)) if mc_endcap_pass_error and mc_endcap_fail_error else 0   )
        else:
            jpsi_mc_eff_endcap.append(-1)
            jpsi_mc_eff_endcap_error_total.append( 0.001 )




        if jpsi_data_eff_barrel[-1] > 0 and jpsi_mc_eff_barrel[-1] > 0:

            if jpsi_data_eff_barrel_error_total[-1] > 0 and jpsi_mc_eff_barrel_error_total[-1]:
                jpsi_scale_factor_barrel.append(jpsi_data_eff_barrel[-1]/jpsi_mc_eff_barrel[-1])
                jpsi_scale_factor_x_barrel.append(jpsi_ptbins_center_fit[pt_bin_number])
                jpsi_scale_factor_x_barrel_binw.append(jpsi_binw_fit[pt_bin_number])
                jpsi_scale_factor_barrel_error.append( ((jpsi_data_eff_barrel_error_total[-1]/jpsi_data_eff_barrel[-1]) |qs| (jpsi_mc_eff_barrel_error_total[-1]/jpsi_mc_eff_barrel[-1]) ) * jpsi_scale_factor_barrel[-1]    )   # HALOOOO
            else:
                jpsi_scale_factor_barrel_black.append(jpsi_data_eff_barrel[-1]/jpsi_mc_eff_barrel[-1])
                jpsi_scale_factor_barrel_black_x.append(jpsi_ptbins_center_fit[pt_bin_number])
                jpsi_scale_factor_barrel_black_x_binw.append(jpsi_binw_fit[pt_bin_number])

        if jpsi_data_eff_endcap[-1] > 0 and jpsi_mc_eff_endcap[-1] > 0:

            if jpsi_data_eff_endcap_error_total[-1] > 0 and jpsi_mc_eff_endcap_error_total[-1] > 0:
                jpsi_scale_factor_endcap.append(jpsi_data_eff_endcap[-1]/jpsi_mc_eff_endcap[-1])
                jpsi_scale_factor_x_endcap.append(jpsi_ptbins_center_fit[pt_bin_number])
                jpsi_scale_factor_x_endcap_binw.append(jpsi_binw_fit[pt_bin_number])
                jpsi_scale_factor_endcap_error.append( ((jpsi_data_eff_endcap_error_total[-1]/jpsi_data_eff_endcap[-1]) |qs| (jpsi_mc_eff_endcap_error_total[-1]/jpsi_mc_eff_endcap[-1]) ) * jpsi_scale_factor_endcap[-1]    )
            else:
                jpsi_scale_factor_endcap_black.append(jpsi_data_eff_endcap[-1]/jpsi_mc_eff_endcap[-1])
                jpsi_scale_factor_endcap_black_x.append(jpsi_ptbins_center_fit[pt_bin_number])
                jpsi_scale_factor_endcap_black_x_binw.append(jpsi_binw_fit[pt_bin_number])

        pt_bin_number += 1


    z_data_eff_barrel = []
    z_data_eff_barrel_error_total = []
    z_mc_eff_barrel = []
    z_mc_eff_barrel_error_total = []
    z_data_eff_endcap = []
    z_data_eff_endcap_error_total = []
    z_mc_eff_endcap = []
    z_mc_eff_endcap_error_total = []

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
    for i in range(0,len(z_data_integrals),4):
        
        if z_data_integrals[i] and z_data_integrals[i+1]:
            data_barrel_pass = z_data_integrals[i]
            data_barrel_pass_error = z_data_integral_sigmas[i] if z_data_integral_sigmas[i] and z_data_integral_sigmas[i]/z_data_integrals[i] < 1 else False
            data_barrel_pass_error_bksyst = abs(z_data_integral_sigmas_bksyst[i]) 
            data_barrel_pass_error_rangesyst = abs(z_data_integral_sigmas_rangesyst[i])
            data_barrel_fail = z_data_integrals[i+1] 
            data_barrel_fail_error = z_data_integral_sigmas[i+1] if z_data_integral_sigmas[i+1] and z_data_integral_sigmas[i+1]/z_data_integrals[i+1] < 1 else False
            data_barrel_fail_error_bksyst = abs(z_data_integral_sigmas_bksyst[i+1]) 
            data_barrel_fail_error_rangesyst = abs(z_data_integral_sigmas_rangesyst[i+1])


            z_data_eff_barrel.append(data_barrel_pass / (data_barrel_pass + data_barrel_fail))
            z_data_eff_barrel_error_total.append(get_efficiency_error(data_barrel_pass, data_barrel_fail, qsum(data_barrel_pass_error, data_barrel_pass_error_bksyst, data_barrel_pass_error_rangesyst), qsum(data_barrel_fail_error, data_barrel_fail_error_bksyst,data_barrel_fail_error_rangesyst)) if data_barrel_pass_error and data_barrel_fail_error else 0  )
        else:
            z_data_eff_barrel.append(-1)
            z_data_eff_barrel_error_total.append( 0.001 )

        if z_data_integrals[i+2] and z_data_integrals[i+3]:
            data_endcap_pass = z_data_integrals[i+2]
            data_endcap_pass_error = z_data_integral_sigmas[i+2] if z_data_integral_sigmas[i+2]  and z_data_integral_sigmas[i+2]/z_data_integrals[i+2] < 1 else False
            data_endcap_pass_error_bksyst = abs(z_data_integral_sigmas_bksyst[i+2]) 
            data_endcap_pass_error_rangesyst = abs(z_data_integral_sigmas_rangesyst[i+2])
            data_endcap_fail = z_data_integrals[i+3]
            data_endcap_fail_error = z_data_integral_sigmas[i+3] if z_data_integral_sigmas[i+3]  and z_data_integral_sigmas[i+3]/z_data_integrals[i+3] < 1 else False
            data_endcap_fail_error_bksyst = abs(z_data_integral_sigmas_bksyst[i+3])
            data_endcap_fail_error_rangesyst = abs(z_data_integral_sigmas_rangesyst[i+3]) 

            z_data_eff_endcap.append(data_endcap_pass / (data_endcap_pass + data_endcap_fail))
            z_data_eff_endcap_error_total.append(  get_efficiency_error(data_endcap_pass, data_endcap_fail, qsum(data_endcap_pass_error, data_endcap_pass_error_bksyst, data_endcap_pass_error_rangesyst), qsum(data_endcap_fail_error, data_endcap_fail_error_bksyst,data_endcap_fail_error_rangesyst)) if data_endcap_pass_error and data_endcap_fail_error else 0 )
        else:
            z_data_eff_endcap.append(-1)
            z_data_eff_endcap_error_total.append( 0.001 )


        if z_mc_integrals[i] and z_mc_integrals[i+1]:
            mc_barrel_pass = z_mc_integrals[i]
            mc_barrel_pass_error = z_mc_integral_sigmas[i] if z_mc_integral_sigmas[i] and z_mc_integral_sigmas[i]/z_mc_integrals[i] < 1 else False
            mc_barrel_pass_error_bksyst = abs(z_mc_integral_sigmas_bksyst[i]) 
            mc_barrel_pass_error_rangesyst = abs(z_mc_integral_sigmas_rangesyst[i])
            mc_barrel_fail = z_mc_integrals[i+1]
            mc_barrel_fail_error = z_mc_integral_sigmas[i+1] if z_mc_integral_sigmas[i+1]  and z_mc_integral_sigmas[i+1]/z_mc_integrals[i+1] < 1 else False
            mc_barrel_fail_error_bksyst = abs(z_mc_integral_sigmas_bksyst[i+1]) 
            mc_barrel_fail_error_rangesyst = abs(z_mc_integral_sigmas_rangesyst[i+1])

            z_mc_eff_barrel.append(mc_barrel_pass / (mc_barrel_pass + mc_barrel_fail))
            z_mc_eff_barrel_error_total.append(get_efficiency_error(mc_barrel_pass, mc_barrel_fail, qsum(mc_barrel_pass_error, mc_barrel_pass_error_bksyst, mc_barrel_pass_error_rangesyst), qsum(mc_barrel_fail_error, mc_barrel_fail_error_bksyst,mc_barrel_fail_error_rangesyst)) if mc_barrel_pass_error and mc_barrel_fail_error else 0   )
        else:
            z_mc_eff_barrel.append(-1)
            z_mc_eff_barrel_error_total.append( 0.001 )



        if z_mc_integrals[i+2] and z_mc_integrals[i+3]:
            mc_endcap_pass = z_mc_integrals[i+2]
            mc_endcap_pass_error = z_mc_integral_sigmas[i+2] if z_mc_integral_sigmas[i+2] and z_mc_integral_sigmas[i+2]/z_mc_integrals[i+2] < 1 else False
            mc_endcap_pass_error_bksyst = abs(z_mc_integral_sigmas_bksyst[i+2]) 
            mc_endcap_pass_error_rangesyst = abs(z_mc_integral_sigmas_rangesyst[i+2])
            mc_endcap_fail = z_mc_integrals[i+3]
            mc_endcap_fail_error = z_mc_integral_sigmas[i+3] if z_mc_integral_sigmas[i+3] and z_mc_integral_sigmas[i+3]/z_mc_integrals[i+3] < 1 else False
            mc_endcap_fail_error_bksyst = abs(z_mc_integral_sigmas_bksyst[i+3]) 
            mc_endcap_fail_error_rangesyst = abs(z_mc_integral_sigmas_rangesyst[i+3])

            z_mc_eff_endcap.append(mc_endcap_pass / (mc_endcap_pass + mc_endcap_fail))
            z_mc_eff_endcap_error_total.append( get_efficiency_error(mc_endcap_pass, mc_endcap_fail, qsum(mc_endcap_pass_error, mc_endcap_pass_error_bksyst,mc_endcap_pass_error_rangesyst), qsum(mc_endcap_fail_error, mc_endcap_fail_error_bksyst,mc_endcap_fail_error_rangesyst)) if mc_endcap_pass_error and mc_endcap_fail_error else 0  )
        else:
            z_mc_eff_endcap.append(-1)
            z_mc_eff_endcap_error_total.append( 0.001 )

        if z_data_eff_barrel[-1] > 0 and z_mc_eff_barrel[-1] > 0:

            if z_data_eff_barrel_error_total[-1] > 0 and z_mc_eff_barrel_error_total[-1]:
                z_scale_factor_barrel.append(z_data_eff_barrel[-1]/z_mc_eff_barrel[-1])
                z_scale_factor_x_barrel.append(z_ptbins_center_fit[pt_bin_number])
                z_scale_factor_x_barrel_binw.append(z_binw_fit[pt_bin_number])
                z_scale_factor_barrel_error.append( ((z_data_eff_barrel_error_total[-1]/z_data_eff_barrel[-1]) |qs| (z_mc_eff_barrel_error_total[-1]/z_mc_eff_barrel[-1]) ) * z_scale_factor_barrel[-1]    )  # HALOOOO!!!
            else:
                z_scale_factor_barrel_black.append(z_data_eff_barrel[-1]/z_mc_eff_barrel[-1])
                z_scale_factor_barrel_black_x.append(z_ptbins_center_fit[pt_bin_number])
                z_scale_factor_barrel_black_x_binw.append(z_binw_fit[pt_bin_number])

        if z_data_eff_endcap[-1] > 0 and z_mc_eff_endcap[-1] > 0:

            if z_data_eff_endcap_error_total[-1] > 0 and z_mc_eff_endcap_error_total[-1] > 0:
                z_scale_factor_endcap.append(z_data_eff_endcap[-1]/z_mc_eff_endcap[-1])
                z_scale_factor_x_endcap.append(z_ptbins_center_fit[pt_bin_number])
                z_scale_factor_x_endcap_binw.append(z_binw_fit[pt_bin_number])
                z_scale_factor_endcap_error.append( ((z_data_eff_endcap_error_total[-1]/z_data_eff_endcap[-1]) |qs| (z_mc_eff_endcap_error_total[-1]/z_mc_eff_endcap[-1]) ) * z_scale_factor_endcap[-1]    )
            else:
                z_scale_factor_endcap_black.append(z_data_eff_endcap[-1]/z_mc_eff_endcap[-1])
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

    ax_barrel.set_xlabel(r"$p_{T}$ [GeV]")
    ax_barrel.set_ylabel("Efficiency")
    ax_barrel.set_title(r"$|\eta | < 1.4442$")
    ax_endcap.set_title(r"$1.566 < |\eta | < 2.50$")
    bottom_ax_barrel.set_xlabel(r"$p_{T}$ [GeV]")
    bottom_ax_barrel.set_ylabel("Data/MC")


    #ax_barrel.plot(ptbins_center_fit, data_eff_barrel, "bs")
    ax_barrel.errorbar(jpsi_ptbins_center_fit, jpsi_data_eff_barrel, yerr=jpsi_data_eff_barrel_error_total, xerr=jpsi_binw_fit, capsize=3, fmt=".", color="blue", ecolor = "blue", label="JPsi Data")
    #ax_barrel.plot(ptbins_center_fit, mc_eff_barrel, "rs")
    ax_barrel.errorbar(jpsi_ptbins_center_fit, jpsi_mc_eff_barrel, yerr=jpsi_mc_eff_barrel_error_total, xerr=jpsi_binw_fit, capsize=3, fmt=".", color="crimson", ecolor = "crimson", label="JPsi MC")


    ax_barrel.errorbar(z_ptbins_center_fit, z_data_eff_barrel, yerr=z_data_eff_barrel_error_total, xerr=z_binw_fit, capsize=3, fmt=".", color="cyan", ecolor = "cyan", label="Z Data")
    ax_barrel.errorbar(z_ptbins_center_fit, z_mc_eff_barrel, yerr=z_mc_eff_barrel_error_total, xerr=z_binw_fit, capsize=3, fmt=".", color="magenta", ecolor = "magenta", label="Z MC")


    bottom_ax_barrel.errorbar(jpsi_scale_factor_x_barrel, jpsi_scale_factor_barrel, yerr=jpsi_scale_factor_barrel_error, xerr=jpsi_scale_factor_x_barrel_binw, capsize=3, fmt=".", color="black", ecolor = "black")
    bottom_ax_barrel.errorbar(jpsi_scale_factor_barrel_black_x, jpsi_scale_factor_barrel_black, xerr=jpsi_scale_factor_barrel_black_x_binw, fmt="x", color="black", ecolor="black")


    bottom_ax_barrel.errorbar(z_scale_factor_x_barrel, z_scale_factor_barrel, yerr=z_scale_factor_barrel_error, xerr=z_scale_factor_x_barrel_binw, capsize=3, fmt=".", color="chocolate", ecolor = "chocolate")
    bottom_ax_barrel.errorbar(z_scale_factor_barrel_black_x, z_scale_factor_barrel_black, xerr=z_scale_factor_barrel_black_x_binw, fmt="x", color="chocolate", ecolor="chocolate")

    #ax_endcap.plot(ptbins_center_fit, data_eff_endcap, "bs")
    ax_endcap.errorbar(jpsi_ptbins_center_fit, jpsi_data_eff_endcap, yerr=jpsi_data_eff_endcap_error_total, xerr=jpsi_binw_fit, capsize=3, fmt=".", color="blue", ecolor = "blue", label="JPsi Data")
    #ax_endcap.plot(ptbins_center_fit, mc_eff_endcap, "rs")
    ax_endcap.errorbar(jpsi_ptbins_center_fit, jpsi_mc_eff_endcap, yerr=jpsi_mc_eff_endcap_error_total, xerr=jpsi_binw_fit, capsize=3, fmt=".", color="crimson", ecolor = "crimson", label="JPsi MC")

    ax_endcap.errorbar(z_ptbins_center_fit, z_data_eff_endcap, yerr=z_data_eff_endcap_error_total, xerr=z_binw_fit, capsize=3, fmt=".", color="cyan", ecolor = "cyan", label="Z Data")
    ax_endcap.errorbar(z_ptbins_center_fit, z_mc_eff_endcap, yerr=z_mc_eff_endcap_error_total, xerr=z_binw_fit, capsize=3, fmt=".", color="magenta", ecolor = "magenta", label="Z MC")


    bottom_ax_endcap.errorbar(jpsi_scale_factor_x_endcap, jpsi_scale_factor_endcap, yerr=jpsi_scale_factor_endcap_error, xerr=jpsi_scale_factor_x_endcap_binw, capsize=3, fmt=".", color="black", ecolor = "black")
    bottom_ax_endcap.errorbar(jpsi_scale_factor_endcap_black_x, jpsi_scale_factor_endcap_black, xerr=jpsi_scale_factor_endcap_black_x_binw, fmt="x", color="black", ecolor="black")


    bottom_ax_endcap.errorbar(z_scale_factor_x_endcap, z_scale_factor_endcap, yerr=z_scale_factor_endcap_error, xerr=z_scale_factor_x_endcap_binw, capsize=3, fmt=".", color="chocolate", ecolor = "chocolate")
    bottom_ax_endcap.errorbar(z_scale_factor_endcap_black_x, z_scale_factor_endcap_black, xerr=z_scale_factor_endcap_black_x_binw, fmt="x", color="chocolate", ecolor="chocolate")

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

