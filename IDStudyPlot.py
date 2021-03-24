import os, sys
import ROOT

sys.path.append('../')
from Helper.PlotHelper import *

def get_parser():
    ''' Argument parser.
    '''
    import argparse
    argParser = argparse.ArgumentParser(description = "Argument parser")
    argParser.add_argument('--sample',           action='store',                     type=str,            default='FullSim100',                                help="Which sample?" )
    argParser.add_argument('--channel',           action='store',                     type=str,            default='All',                                help="Which particle" )
    return argParser

options = get_parser().parse_args()



sample  = options.sample
channel = options.channel

def plotEff(hnum, hden, xtitle, name, legTitle, sample, xmin = 0, ymin =0, xmax=50, ymax=1.2):
    heff = ROOT.TGraphAsymmErrors()
    heff.BayesDivide(hnum,hden)

    heff.SetLineColor(ROOT.kBlack)
    heff.SetLineWidth(2)
    heff.SetMarkerSize(0.8)
    heff.SetMarkerStyle(20)
    heff.SetMarkerColor(ROOT.kBlack)

    leg = ROOT.TLegend(0.5, 0.8, 0.9, 0.9)
    leg.AddEntry(heff, legTitle ,"p")
    c = ROOT.TCanvas('c', '', 600, 800)
    fr = c.DrawFrame(xmin, ymin, xmax, ymax)
    fr.GetYaxis().SetTitle('Eff')
    fr.GetYaxis().SetTitleSize(0.05)
    fr.GetYaxis().SetTitleOffset(0.7)
    fr.GetYaxis().SetLabelSize(0.03)
    fr.GetXaxis().SetTitle(xtitle)
    fr.GetXaxis().SetTitleSize(0.05)
    fr.GetXaxis().SetTitleOffset(0.9)
    fr.GetXaxis().SetLabelSize(0.03)
    heff.Draw("P")
    leg.Draw("SAME")
    c.SaveAs("EffPlots/"+sample+"/"+name+".png")
    c.Close()
    

def plotStackedEff(hnum, hden, colors, xtitle, name, legTitle, sample, xmin = 0, ymin =0, xmax=50, ymax=1.1, legPos = "tr"):
    if( len(hnum) > len(colors)):
        print("Error in plotStackedEff - arrays are incompatible: ",len(hnum), len(colors) )
        return
    
    heff = {}
    for i in range(len(hnum)):
        heff[i] = ROOT.TGraphAsymmErrors()
        heff[i].BayesDivide(hnum[i],hden)

        heff[i].SetLineColor(colors[i])
        heff[i].SetLineWidth(2)
        heff[i].SetMarkerSize(0.8)
        heff[i].SetMarkerStyle(20)
        heff[i].SetMarkerColor(colors[i])
        heff[i].SetName("Efficiency in sample "+sample)

    if legPos == "tr":
        leg = ROOT.TLegend(0.5, 0.8, 0.9, 0.9)
    elif legPos == "br":
        leg = ROOT.TLegend(0.5, 0.1, 0.9, 0.2)
    elif legPos == "tl":
        leg = ROOT.TLegend(0.1, 0.8, 0.5, 0.9)
    elif legPos == "bl":
        leg = ROOT.TLegend(0.1, 0.1, 0.5, 0.2)
    elif legPos == "tc":
        leg = ROOT.TLegend(0.3, 0.8, 0.7, 0.9)
    elif legPos == "bc":
        leg = ROOT.TLegend(0.3, 0.1, 0.7, 0.2)
    else:
        print("Error: Invalid TLegend position.")
        return
    
    for i in range(len(heff)):    
        leg.AddEntry(heff[i], legTitle[i] ,"p")

    c = ROOT.TCanvas('c', '', 600, 800)
    fr = c.DrawFrame(xmin, ymin, xmax, ymax)
    fr.GetYaxis().SetTitle('Eff')
    fr.GetYaxis().SetTitleSize(0.05)
    fr.GetYaxis().SetTitleOffset(0.7)
    fr.GetYaxis().SetLabelSize(0.03)
    fr.GetXaxis().SetTitle(xtitle)
    fr.GetXaxis().SetTitleSize(0.05)
    fr.GetXaxis().SetTitleOffset(0.9)
    fr.GetXaxis().SetLabelSize(0.03)
    for i in range(len(heff)):
        heff[i].Draw("P,SAME")
    leg.Draw("SAME")
    c.SaveAs("IDStudyPlots/"+sample+"/"+name+".png")
    c.Close()




f = ROOT.TFile.Open('root_files/IDStudy_Sample'+sample+'.root')

h_pass = f.Get("RecoPhotonEta_")
h_total = f.Get("TruePhotonEta_")
plotEff(h_pass, h_total, 'photon eta [GeV]', 'photon_eta-eff', "Reco photon", sample,-3,0,3)
h_pass = f.Get("IsoTrackEta_")
h_total = f.Get("TrueEleEta_")
plotEff(h_pass, h_total, 'isotrack matched ele eta [GeV]', 'isoTrack_eta-eff', "matched isoTrack", sample,-3,0,3)


h_pass = {}
colors = [ROOT.kBlack,ROOT.kRed,ROOT.kBlue,ROOT.kGreen,ROOT.kViolet,ROOT.kYellow]
legTitleEle = ["Nofilter ele","Veto Ele","Loose ele","Med ele","Tight ele"]
legTitleMuon = ["Nofilter muon","Loose muon","Med mon","Tight muon"]

h_pass[0] = f.Get("NofilterElePt_")
h_pass[1] = f.Get("VetoElePt_")
h_pass[2] = f.Get("LooseElePt_")
h_pass[3] = f.Get("MediumElePt_")
h_pass[4] = f.Get("TightElePt_")

h_total = f.Get("TrueElePt_")


plotStackedEff(h_pass, h_total, colors, 'p_{T} [GeV]', 'electron_common_pT-eff', legTitleEle, sample,0,0,50)

h_pass[0] = f.Get("NofilterEleEta_")
h_pass[1] = f.Get("VetoEleEta_")
h_pass[2] = f.Get("LooseEleEta_")
h_pass[3] = f.Get("MediumEleEta_")
h_pass[4] = f.Get("TightEleEta_")

h_total = f.Get("TrueEleEta_")

plotStackedEff(h_pass, h_total, colors, '\eta', 'electron_common_eta-eff', legTitleEle, sample,-3,0,3)


'''
h_pass = {}
h_pass[0] = f.Get("NofilterMuonPt_")
h_pass[1] = f.Get("LooseMuonPt_")
h_pass[2] = f.Get("MediumMuonPt_")
h_pass[3] = f.Get("TightMuonPt_")
h_total = f.Get("TrueMuonPt_")

plotStackedEff(h_pass, h_total, colors, 'p_{T} [GeV]', 'muon_common_pT-eff', legTitleMuon, sample,0,0,50,1.1,"tr")

h_pass[0] = f.Get("NofilterMuonEta_")
h_pass[1] = f.Get("LooseMuonEta_")
h_pass[2] = f.Get("MediumMuonEta_")
h_pass[3] = f.Get("TightMuonEta_")
h_total = f.Get("TrueMuonEta_")

plotStackedEff(h_pass, h_total, colors, '\eta', 'muon_common_eta-eff', legTitleMuon, sample,-3,0,3,1.1,"tr")
'''
