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
    argParser.add_argument('--threaded',           action='store',                     type=str,            default='1',                                help="Was it the threaded code?" )
    return argParser

options = get_parser().parse_args()



sample  = options.sample
channel = options.channel
bThreaded = options.threaded

def Plot2D(h, dir, xmin = 0, ymin =0, xmax=1, ymax=1, drawOption="hist", islogz=False, canvasX=600, canvasY=800):
    hname = h.GetName()
    htitle = h.GetTitle()
    sname = hname.replace(htitle+"_", "")
    outputdirpath = os.path.join(dir,"2DPlots/final",sname)
    if not os.path.exists(outputdirpath):
        os.makedirs(outputdirpath)

    leg = ROOT.TLegend(0.5, 0.85, 0.9, 0.9)
    leg.AddEntry(h, sname ,"l")

    #style2D(h, islogy)
    
    c = ROOT.TCanvas('c', '', canvasX, canvasY)
    if(islogz):
        c.SetLogz()
    fr = c.DrawFrame(xmin, ymin, xmax, ymax)
    c.cd()
    h.Draw(drawOption)
    leg.Draw("SAME")
    c.SaveAs(outputdirpath+"/"+htitle+".png")
    c.Close()



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
    if( not (len(hnum) == len(hden)) or ( len(hnum) > len(colors))):
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

def plotStackedEffAndRatio(hnum, hden, colors, xtitle, name, legTitle, sample, xmin = 0, ymin =-0.2, xmax=50, ymax=1.2, yratiomin = 0.5, yratiomax = 13):
    if( not (len(hnum) == len(hden)) or ( len(hnum) > len(colors))):
        print("Error in plotStackedEffAndRatio - arrays are incompatible.")
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

    # Use: 0 is reco fixed, 1 is simple reco
    #xmax = max(ROOT.TMath.MaxElement(hnum.GetN(),hnum.GetX()), ROOT.TMath.MaxElement(hden.GetN(),hden.GetX()))
    lxnum = [x for x in heff[0].GetX()]
    lynum = [x for x in heff[0].GetY()]
    lxden = [x for x in heff[1].GetX()]
    lyden = [x for x in heff[1].GetY()]


    Nbins = hnum[0].GetNbinsX()
    binwidth = hnum[0].GetBinWidth(1)
    #hRatio = ROOT.TH1F("hRatio", "hRatio", Nbins, 0, Nbins*binwidth)
    hRatio = ROOT.TH1F("hRatio", "hRatio", Nbins, hnum[0].GetBin(1), hnum[0].GetBin(Nbins))

    #print("New plot")
    #print("binwidth: ",binwidth)
    #print(lxnum)
    #print(lxden)

    for b in range(0, Nbins):
        x = binwidth/2 + b*binwidth
        #print("x is",x)
        if x in lxnum:
            id1 = lxnum.index(x)
            tnum = lynum[id1]
        else: tnum = -1.0
        if x in lxden:
            id1 = lxden.index(x)
            tden = lyden[id1]
        else: tden = -1.0
        #ratio = -1.0 if (tnum==-1 or hden==-1 or tden==0) else tnum/tden 
        ratio = -1.0 if (tnum==-1 or tden==0) else tnum/tden 
        #print(tnum, tden, ratio)
        hRatio.SetBinContent(b+1, ratio)


    hRatio.GetYaxis().SetTitle("fixed / simple")
    hRatio.GetYaxis().SetRangeUser(yratiomin,yratiomax)
    hRatio.SetTitle("")
    hRatio.SetStats(0)
    hRatio.SetLineColor(ROOT.kRed)
    hRatio.SetMarkerColor(ROOT.kRed)
    hRatio.SetMarkerSize(0.8)
    hRatio.SetMarkerStyle(20)
    hRatio.GetYaxis().SetTitleSize(0.08)
    hRatio.GetYaxis().SetTitleOffset(0.5)
    hRatio.GetYaxis().SetLabelSize(0.07)
    hRatio.GetXaxis().SetTitle(xtitle)
    hRatio.GetXaxis().SetTitleSize(0.1)
    hRatio.GetXaxis().SetTitleOffset(0.9)
    hRatio.GetXaxis().SetLabelSize(0.07)

    hRatioFrame = hRatio.Clone("RatioFrame")
    for b in range(1, hRatioFrame.GetNbinsX() + 1):
        hRatioFrame.SetBinContent(b, 1.0)
    hRatioFrame.SetLineStyle(2)
    hRatioFrame.SetLineWidth(2)
    hRatioFrame.SetLineColor(ROOT.kGreen)


    leg = ROOT.TLegend(0.5, 0.8, 0.9, 0.9)
    for i in range(len(heff)):    
        leg.AddEntry(heff[i], legTitle[i] ,"p")

    c = ROOT.TCanvas('c', '', 600, 800)
    p1 = ROOT.TPad("p1", "p1", 0, 0.3, 1, 1.0)
    p1.SetBottomMargin(0) # Upper and lower plot are joined
    p1.Draw()             # Draw the upper pad: p1
    p1.cd()
    fr = p1.DrawFrame(xmin, ymin, xmax, ymax)
    #fr = p1.DrawFrame(0., 0.0, Nbins*binwidth, 1.2)
    #fr = c.DrawFrame(xmin, ymin, xmax, ymax)
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
    #hRatio.Draw("SAME")
    #for i in range(1, hRatioFrame.GetNbinsX() + 1):
    #    print(hRatio.GetBinContent(i))

    c.cd()
    p2 = ROOT.TPad("p2", "p2", 0, 0.01, 1, 0.3)
    p2.SetTopMargin(0)
    p2.SetBottomMargin(0.2)
    p2.Draw()
    p2.cd()
    hRatio.Draw("P")
    hRatioFrame.Draw("HISTsame")
    c.SaveAs("RecoFixPlots/"+sample+"/"+name+".png")
    c.Close()

if(bThreaded):
    f = ROOT.TFile.Open('root_files/RecoFixThreadedSTACK_Sample'+sample+'.root')
else:   
    f = ROOT.TFile.Open('root_files/RecoFix_Sample'+sample+'.root')

if ("Electron" == channel or "All" == channel):
    h_pass = f.Get("RecoElePt_")
    h_total = f.Get("TrueElePt_")
    plotEff(h_pass, h_total, 'ele p_{T} [GeV]', 'electron_pT-eff', "Reco ele", sample)

    h_pass = f.Get("RecoEleEta_")
    h_total = f.Get("TrueEleEta_")
    plotEff(h_pass, h_total, 'ele eta [GeV]', 'electron_eta-eff', "Reco ele", sample,-3,0,3)

    Plot1D(h_pass ,'./RecoFixPlots')
    Plot1D(h_total ,'./RecoFixPlots')

    h_pass = f.Get("RecoEleVtxx_")
    h_total = f.Get("TrueEleVtxx_")
    plotEff(h_pass, h_total, 'ele vtx x [mm or cm have to check]', 'electron_vtxx-eff', "Reco ele", sample,-3,0,3)

    h_pass = f.Get("RecoEleVtxy_")
    h_total = f.Get("TrueEleVtxy_")
    plotEff(h_pass, h_total, 'ele vtx y [mm or cm have to check]', 'electron_vtxy-eff', "Reco ele", sample,-3,0,3)

    h_pass = f.Get("RecoEleVtxz_")
    h_total = f.Get("TrueEleVtxz_")
    plotEff(h_pass, h_total, 'ele vtx z [mm or cm have to check]', 'electron_vtxz-eff', "Reco ele", sample,-3,0,3)


if("Photon" == channel or "All" == channel):
    h_pass = f.Get("RecoPhotonPt_")
    h_total = f.Get("TrueElePt_")
    Plot1D(h_pass ,'./RecoFixPlots')
    plotEff(h_pass, h_total, 'photon p_{T} [GeV]', 'photon_pT-eff', "Reco photon", sample)

    h_pass = f.Get("RecoPhotonEta_")
    h_total = f.Get("TrueEleEta_")
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
    colors = [ROOT.kBlue,ROOT.kBlack,ROOT.kRed,ROOT.kOrange,ROOT.kGreen]
    legTitle = ["Reco fixed ele","Reco ele","matched isoTrack","Reco photon","Reco pho + isoTrack"]
    h_pass[0] = f.Get("RecoFixedElePt_")
    h_pass[1] = f.Get("RecoElePt_")
    h_pass[2] = f.Get("IsoTrackPt_")
    h_pass[3] = f.Get("RecoPhotonPt_")
    #h_pass[4] = f.Get("RecoPhotonPlusIsoTrackPt_")

    h_total[0] = f.Get("TrueElePt_")
    h_total[1] = f.Get("TrueElePt_")
    h_total[2] = f.Get("TrueElePt_")
    h_total[3] = f.Get("TrueElePt_")
    #h_total[4] = f.Get("TrueElePt_")

    plotStackedEffAndRatio(h_pass, h_total, colors, 'p_{T} [GeV]', 'electron_common_pT-eff', legTitle, sample,0,-0.05,50,1.2,0.5,2)

    h_pass = {}
    h_total = {}
    legTitle = ["Reco fixed ele","Reco ele","matched isoTrack","Reco photon","Reco pho + isoTrack"]
    h_pass[0] = f.Get("RecoFixedEleEta_")
    h_pass[1] = f.Get("RecoEleEta_")
    h_pass[2] = f.Get("IsoTrackEta_")
    h_pass[3] = f.Get("RecoPhotonEta_")
    #h_pass[4] = f.Get("RecoPhotonPlusIsoTrackEta_")

    h_total[0] = f.Get("TrueEleEta_")
    h_total[1] = f.Get("TrueEleEta_")
    h_total[2] = f.Get("TrueEleEta_")
    h_total[3] = f.Get("TrueEleEta_")
    #h_total[4] = f.Get("TrueEleEta_")

    plotStackedEffAndRatio(h_pass, h_total, colors, 'p_{T} [GeV]', 'electron_common_eta-eff', legTitle, sample,-3,-0.05,3,1.2,0.5,2)

    h_pass = {}
    h_total = {}
    colors = [ROOT.kBlue,ROOT.kBlack,ROOT.kRed]
    legTitle = ["Reco fixed ele","Reco ele","matched isoTrack"]
    h_pass[0] = f.Get("RecoFixedEleVtxx_")
    h_pass[1] = f.Get("RecoEleVtxx_")
    h_pass[2] = f.Get("IsoTrackVtxx_")

    h_total[0] = f.Get("TrueEleVtxx_")
    h_total[1] = f.Get("TrueEleVtxx_")
    h_total[2] = f.Get("TrueEleVtxx_")

    plotStackedEff(h_pass, h_total, colors, 'vtx x [cm or mm I dont remember now]', 'electron_common_vtxx-eff', legTitle, sample,-3,0,3)


    h_pass = {}
    h_total = {}
    colors = [ROOT.kBlue,ROOT.kBlack,ROOT.kRed]
    legTitle = ["Reco fixed ele","Reco ele","matched isoTrack"]
    h_pass[0] = f.Get("RecoFixedEleVtxy_")
    h_pass[1] = f.Get("RecoEleVtxy_")
    h_pass[2] = f.Get("IsoTrackVtxy_")

    h_total[0] = f.Get("TrueEleVtxy_")
    h_total[1] = f.Get("TrueEleVtxy_")
    h_total[2] = f.Get("TrueEleVtxy_")

    plotStackedEff(h_pass, h_total, colors, 'vtx y [cm or mm I dont remember now]', 'electron_common_vtxy-eff', legTitle, sample,-3,0,3)

    h_pass = {}
    h_total = {}
    colors = [ROOT.kBlue,ROOT.kBlack,ROOT.kRed]
    legTitle = ["Reco fixed ele","Reco ele","matched isoTrack"]
    h_pass[0] = f.Get("RecoFixedEleVtxz_")
    h_pass[1] = f.Get("RecoEleVtxz_")
    h_pass[2] = f.Get("IsoTrackVtxz_")

    h_total[0] = f.Get("TrueEleVtxz_")
    h_total[1] = f.Get("TrueEleVtxz_")
    h_total[2] = f.Get("TrueEleVtxz_")

    plotStackedEff(h_pass, h_total, colors, 'vtx z [cm or mm I dont remember now]', 'electron_common_vtxz-eff', legTitle, sample,-3,0,3)

    h_pass = {}
    h_total = {}
    colors = [ROOT.kBlue,ROOT.kBlack,ROOT.kRed]
    legTitle = ["Reco fixed ele","Reco ele","matched isoTrack"]
    h_pass[0] = f.Get("RecoFixedEleDxy_")
    h_pass[1] = f.Get("RecoEleDxy_")
    h_pass[2] = f.Get("IsoTrackDxy_")

    h_total[0] = f.Get("TrueEleDxy_")
    h_total[1] = f.Get("TrueEleDxy_")
    h_total[2] = f.Get("TrueEleDxy_")

    plotStackedEff(h_pass, h_total, colors, 'dxy [cm]', 'electron_common_dxy-eff', legTitle, sample,-3,0,3)

    h_pass = {}
    h_total = {}
    colors = [ROOT.kBlue,ROOT.kBlack,ROOT.kRed]
    legTitle = ["Reco fixed ele","Reco ele","matched isoTrack"]
    h_pass[0] = f.Get("RecoFixedEleDz_")
    h_pass[1] = f.Get("RecoEleDz_")
    h_pass[2] = f.Get("IsoTrackDz_")

    h_total[0] = f.Get("TrueEleDz_")
    h_total[1] = f.Get("TrueEleDz_")
    h_total[2] = f.Get("TrueEleDz_")

    plotStackedEff(h_pass, h_total, colors, 'dz [cm]', 'electron_common_dz-eff', legTitle, sample,-3,0,3)


    h_pass = {}
    h_total = {}
    colors = [ROOT.kBlue,ROOT.kBlack,ROOT.kRed,ROOT.kOrange,ROOT.kGreen]
    legTitle = ["Reco fixed ele","Reco ele","matched isoTrack","Reco photon","Reco pho + isoTrack"]
    h_pass[0] = f.Get("RecoFixedEleDzAbs_")
    h_pass[1] = f.Get("RecoEleDzAbs_")
    h_pass[2] = f.Get("IsoTrackDzAbs_")
    h_pass[3] = f.Get("RecoPhotonDzAbs_")
    h_pass[4] = f.Get("RecoPhotonPlusIsoTrackDzAbs_")

    h_total[0] = f.Get("TrueEleDzAbs_")
    h_total[1] = f.Get("TrueEleDzAbs_")
    h_total[2] = f.Get("TrueEleDzAbs_")
    h_total[3] = f.Get("TrueEleDzAbs_")
    h_total[4] = f.Get("TrueEleDzAbs_")

    plotStackedEffAndRatio(h_pass, h_total, colors, 'abs(dz) [cm]', 'electron_common_dzAbs-eff', legTitle, sample,0,-0.05,20,1.2,0.5,2)

    h_pass = {}
    h_total = {}
    colors = [ROOT.kBlue,ROOT.kBlack,ROOT.kRed,ROOT.kOrange,ROOT.kGreen]
    legTitle = ["Reco fixed ele","Reco ele","matched isoTrack","Reco photon","Reco pho + isoTrack"]
    h_pass[0] = f.Get("RecoFixedEleDxyAbs_")
    h_pass[1] = f.Get("RecoEleDxyAbs_")
    h_pass[2] = f.Get("IsoTrackDxyAbs_")
    h_pass[3] = f.Get("RecoPhotonDxyAbs_")
    h_pass[4] = f.Get("RecoPhotonPlusIsoTrackDxyAbs_")

    h_total[0] = f.Get("TrueEleDxyAbs_")
    h_total[1] = f.Get("TrueEleDxyAbs_")
    h_total[2] = f.Get("TrueEleDxyAbs_")
    h_total[3] = f.Get("TrueEleDxyAbs_")
    h_total[4] = f.Get("TrueEleDzAbs_")

    plotStackedEffAndRatio(h_pass, h_total, colors, 'abs(dxy) [cm]', 'electron_common_dxyAbs-eff', legTitle, sample,0,-0.05,20)


    Plot2D( f.Get("dRControlPlot2D_") ,'./RecoFixPlots',xmin = 0, ymin =0, xmax=1, ymax=1, drawOption="COLZ",islogz=True)


    h_pass = {}
    h_total = {}
    legTitle = ["Reco fixed ele","Reco ele"]
    h_pass[0] = f.Get("RejectBackgroundRecoFixPt_")
    h_pass[1] = f.Get("RejectBackgroundPt_")

    h_total[0] = f.Get("TrueBackgroundPt_")
    h_total[1] = f.Get("TrueBackgroundPt_")

    plotStackedEffAndRatio(h_pass, h_total, colors, 'p_{T} [GeV]', 'background_common_pT', legTitle, sample,0,-0.05,50,1.2,0,1.2)

    h_pass = {}
    h_total = {}
    legTitle = ["Reco fixed ele","Reco ele"]
    h_pass[0] = f.Get("RejectBackgroundRecoFixEta_")
    h_pass[1] = f.Get("RejectBackgroundEta_")

    h_total[0] = f.Get("TrueBackgroundEta_")
    h_total[1] = f.Get("TrueBackgroundEta_")


    plotStackedEffAndRatio(h_pass, h_total, colors, '\eta ', 'background_common_eta', legTitle, sample,-3,-0.05,3,1.2,0,1.2)
    #plotStackedEff(h_pass, h_total, colors, '\eta ', 'background_common_eta', legTitle, sample,-3,0,3)

    h_pass = {}
    h_total = {}
    legTitle = ["Reco fixed ele","Reco ele"]
    h_pass[0] = f.Get("RejectBackgroundRecoFixDxyAbs_")
    h_pass[1] = f.Get("RejectBackgroundDxyAbs_")

    h_total[0] = f.Get("TrueBackgroundDxyAbs_")
    h_total[1] = f.Get("TrueBackgroundDxyAbs_")

    plotStackedEffAndRatio(h_pass, h_total, colors, 'abs(dxy) [cm] ', 'background_common_dxyAbs', legTitle, sample, 0,-0.05,20,1.2,0,1.2)

    h_pass = {}
    h_total = {}
    legTitle = ["Reco fixed ele","Reco ele"]
    h_pass[0] = f.Get("RejectBackgroundRecoFixDzAbs_")
    h_pass[1] = f.Get("RejectBackgroundDzAbs_")

    h_total[0] = f.Get("TrueBackgroundDzAbs_")
    h_total[1] = f.Get("TrueBackgroundDzAbs_")

    plotStackedEffAndRatio(h_pass, h_total, colors, 'abs(dz) [cm]', 'background_common_dzAbs', legTitle, sample, 0,-0.05,20,1.2,0,1.2)

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
