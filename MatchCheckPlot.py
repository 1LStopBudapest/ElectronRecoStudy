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
    c.SaveAs("RecoFixPlots/"+sample+"/"+name+".png")
    c.Close()
    

def plotStackedEff(hnum, hden, colors, xtitle, name, legTitle, sample, xmin = 0, ymin =0, xmax=50, ymax=1.2):
    if( not (len(hnum) == len(hden)) or (not (len(hnum) == len(colors)))):
        print("Error in plotStackedEff - arrays are incompatible.")
        return
    
    heff = {}
    for i in range(len(hnum)):
        heff[i] = ROOT.TGraphAsymmErrors()
        heff[i].BayesDivide(hnum[i],hden[i])

        heff[i].SetLineColor(colors[i])
        heff[i].SetLineWidth(2)
        heff[i].SetMarkerSize(0.8)
        heff[i].SetMarkerStyle(20)
        heff[i].SetMarkerColor(colors[i])

    leg = ROOT.TLegend(0.5, 0.8, 0.9, 0.9)
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
    c.SaveAs("RecoFixPlots/"+sample+"/"+name+".png")
    c.Close()

def PlotStacked1D(h, title, colors,drawOption="hist", islogy=False, canvasX=600, canvasY=800):
    hname = {}
    htitle = {}
    sname = {}
    
    for i in range(len(h)):
        h[i].SetLineColor(colors[i])
        h[i].SetLineWidth(2)
        h[i].SetMarkerSize(0.8)
        h[i].SetMarkerStyle(20)
        h[i].SetMarkerColor(colors[i])
        hname[i] = h[i].GetName()
        htitle[i] = h[i].GetTitle()
        sname[i] = hname[i].replace(htitle[i]+"_", "")

    leg = ROOT.TLegend(0.5, 0.85, 0.9, 0.9)
    for i in range(len(h)):
        leg.AddEntry(h[i], sname[i] ,"l")

    style1D(h[0], islogy)
    
    c = ROOT.TCanvas('c', '', canvasX, canvasY)
    c.cd()
    h[0].Draw(drawOption)
    for i in range(1,len(h),1):
        h[i].Draw(drawOption+"SAME")
    leg.Draw("SAME")
    if islogy:ROOT.gPad.SetLogy()
    c.SaveAs("./MatchCheckPlots/"+title+".png")
    c.Close()


f = ROOT.TFile.Open('root_files/MatchCheck_Sample'+sample+'.root')

Plot1D( f.Get("AllPhotonEta_") ,'./MatchCheckPlots')
Plot1D( f.Get("MatchedPhotonEta_") ,'./MatchCheckPlots')
Plot1D( f.Get("AllIsoTrackEta_") ,'./MatchCheckPlots')
Plot1D( f.Get("MatchedIsoTrackEta_") ,'./MatchCheckPlots')
Plot1D( f.Get("DoubleMatchedIsoTrackEta_") ,'./MatchCheckPlots')
Plot1D( f.Get("AllElectronEta_") ,'./MatchCheckPlots')
Plot1D( f.Get("MatchedElectronEta_") ,'./MatchCheckPlots')



colors = [ROOT.kBlack,ROOT.kRed,ROOT.kBlue,ROOT.kPink,ROOT.kOrange]
h = {}
h[0] = f.Get("AllPhotonEta_")
h[1] = f.Get("MatchedPhotonEta_")
h[2] = f.Get("MatchedPhotonGenEta_")
PlotStacked1D(h, 'common_photon_distr_eta', colors)
h = {}
h[0] = f.Get("AllIsoTrackEta_")
h[1] = f.Get("MatchedIsoTrackEta_")
h[2] = f.Get("MatchedIsoTrackGenEta_")
h[3] = f.Get("DoubleMatchedIsoTrackEta_")
h[4] = f.Get("DoubleMatchedIsoTrackGenEta_")
PlotStacked1D(h, 'common_iso_distr_eta', colors)

h = {}
h[0] = f.Get("AllElectronEta_")
h[1] = f.Get("MatchedElectronEta_")
h[2] = f.Get("MatchedElectronGenEta_")
PlotStacked1D(h, 'common_ele_distr_eta',colors)



'''
if ("Electron" == channel or "All" == channel):
    h_pass = f.Get("RecoElePt_")
    h_total = f.Get("TrueElePt_")
    plotEff(h_pass, h_total, 'ele p_{T} [GeV]', 'electron_pT-eff', "Reco ele", sample)

    h_pass = f.Get("RecoEleEta_")
    h_total = f.Get("TrueEleEta_")
    plotEff(h_pass, h_total, 'ele eta [GeV]', 'electron_eta-eff', "Reco ele", sample,-3,0,3)
if("Photon" == channel or "All" == channel):
    h_pass = f.Get("RecoPhotonPt_")
    h_total = f.Get("TruePhotonPt_")
    plotEff(h_pass, h_total, 'photon p_{T} [GeV]', 'photon_pT-eff', "Reco photon", sample)

    h_pass = f.Get("RecoPhotonEta_")
    h_total = f.Get("TruePhotonEta_")
    plotEff(h_pass, h_total, 'photon eta [GeV]', 'photon_eta-eff', "Reco photon", sample,-3,0,3)
if("IsoTrack" == channel or "All" == channel):
    h_pass = f.Get("IsoTrackPt_")
    h_total = f.Get("TrueElePt_")
    plotEff(h_pass, h_total, 'isotrack matched ele p_{T} [GeV]', 'isoTrack_pT-eff', "matched isoTrack", sample)

    h_pass = f.Get("IsoTrackEta_")
    h_total = f.Get("TrueEleEta_")
    plotEff(h_pass, h_total, 'isotrack matched ele eta [GeV]', 'isoTrack_eta-eff', "matched isoTrack", sample,-3,0,3)
if("RecoFix" == channel or "All" == channel):
    h_pass = f.Get("RecoFixedElePt_")
    h_total = f.Get("TrueElePt_")
    plotEff(h_pass, h_total, 'ele p_{T} with reco fix [GeV]', 'electron_recofix_pT-eff', "Reco fixed ele", sample)

    h_pass = f.Get("RecoFixedEleEta_")
    h_total = f.Get("TrueEleEta_")
    plotEff(h_pass, h_total, 'ele eta with reco fix[GeV]', 'electron_recofix_eta-eff', "Reco fixed ele", sample,-3,0,3)
if("CommonReco" == channel or "All" == channel):
    h_pass = {}
    h_total = {}
    colors = [ROOT.kBlack,ROOT.kRed,ROOT.kBlue]
    legTitle = ["Reco fixed ele","matched isoTrack","Reco ele"]
    h_pass[0] = f.Get("RecoFixedElePt_")
    h_pass[1] = f.Get("IsoTrackPt_")
    h_pass[2] = f.Get("RecoElePt_")

    h_total[0] = f.Get("TrueElePt_")
    h_total[1] = f.Get("TrueElePt_")
    h_total[2] = f.Get("TrueElePt_")

    plotStackedEff(h_pass, h_total, colors, 'p_{T} [GeV]', 'electron_common_pT-eff', legTitle, sample,0,0,50)

    h_pass = {}
    h_total = {}
    colors = [ROOT.kBlack,ROOT.kRed,ROOT.kBlue]
    legTitle = ["Reco fixed ele","matched isoTrack","Reco ele"]
    h_pass[0] = f.Get("RecoFixedEleEta_")
    h_pass[1] = f.Get("IsoTrackEta_")
    h_pass[2] = f.Get("RecoEleEta_")

    h_total[0] = f.Get("TrueEleEta_")
    h_total[1] = f.Get("TrueEleEta_")
    h_total[2] = f.Get("TrueEleEta_")

    plotStackedEff(h_pass, h_total, colors, 'p_{T} [GeV]', 'electron_common_eta-eff', legTitle, sample,-3,0,3)
if("testing" == channel):
    h = f.Get("momID_")
    hname = h.GetName()
    htitle = h.GetTitle()
    sname = hname.replace(htitle+"_", "")
    outputdirpath = os.path.join("./","1DPlots",sname)
    if not os.path.exists(outputdirpath):
        os.makedirs(outputdirpath)

    leg = ROOT.TLegend(0.5, 0.85, 0.9, 0.9)
    leg.AddEntry(h, sname ,"l")
    
    c = ROOT.TCanvas('c', '', 600, 800)
    fr = c.DrawFrame(0, 0, 50, 16000)
    fr.GetYaxis().SetTitle('Eff')
    fr.GetYaxis().SetTitleSize(0.05)
    fr.GetYaxis().SetTitleOffset(0.7)
    fr.GetYaxis().SetLabelSize(0.03)
    fr.GetXaxis().SetTitle("momID")
    fr.GetXaxis().SetTitleSize(0.05)
    fr.GetXaxis().SetTitleOffset(0.9)
    fr.GetXaxis().SetLabelSize(0.03)
    h.Draw("hist")
    leg.Draw("SAME")
    c.SaveAs(outputdirpath+"/"+htitle+".png")
    c.Close()
if("testing2" == channel):
    g = ROOT.TFile.Open('root_files/RecoFix_Sample'+sample+'m.root')
    h_pass = {}
    h_total = {}
    colors = [ROOT.kBlack,ROOT.kRed]
    legTitle = ["Reco ele FullSim100","Reco ele FullSim100m"]
    h_pass[0] = f.Get("RecoElePt_")
    h_pass[1] = g.Get("RecoElePt_")

    h_total[0] = f.Get("TrueElePt_")
    h_total[1] = g.Get("TrueElePt_")

    plotStackedEff(h_pass, h_total, colors, 'p_{T} [GeV]', 'electron_beforeafter_pT-eff', legTitle, sample,0,0,50)
'''