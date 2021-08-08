import os, sys
import ROOT
import numpy as np

sys.path.append('../')
from Helper.PlotHelper import *

def get_parser():
    ''' Argument parser.
    '''
    import argparse
    argParser = argparse.ArgumentParser(description = "Argument parser")
    argParser.add_argument('--sample',           action='store',                     type=str,            default='UL17_Full99mm',                                help="Which sample?" )
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

def Plot2DEff(hnum, hden, dir, xmin = 0, ymin =0, xmax=1, ymax=1, drawOption="hist", islogz=False, canvasX=1400, canvasY=800):
    for i in range(1,hnum.GetBin(hnum.GetNbinsX(), hnum.GetNbinsY()) +1, 1):
        denom = hden.GetBinContent(i)
        if(denom == 0):
            hnum.SetBinContent(i, 0  )
        else:
            hnum.SetBinContent(i, hnum.GetBinContent(i) / denom )


    hname = hnum.GetName()
    htitle = hnum.GetTitle()
    sname = hname.replace(htitle+"_", "")
    outputdirpath = os.path.join(dir,"2DPlots/final",sname)
    if not os.path.exists(outputdirpath):
        os.makedirs(outputdirpath)

    leg = ROOT.TLegend(0.5, 0.85, 0.9, 0.9)
    leg.AddEntry(hnum, sname ,"l")

    ROOT.gStyle.SetPaintTextFormat('4.1f')
    ROOT.gStyle.SetOptStat(0)
    c = ROOT.TCanvas('c', '', canvasX, canvasY)
    if(islogz):
        c.SetLogz()
    fr = c.DrawFrame(xmin, ymin, xmax, ymax)
    c.cd()
    hnum.Draw(drawOption)
    leg.Draw("SAME")
    c.SaveAs(outputdirpath+"/"+htitle+".png")
    c.Close()



def plotEff(hnum, hden, xtitle, name, legTitle, sample, xmin = 0, ymin =0, xmax=50, ymax=1.2, ytitle='Eff'):
    heff = ROOT.TGraphAsymmErrors()
    heff.BayesDivide(hnum,hden)

    heff.SetLineColor(ROOT.kBlack)
    heff.SetLineWidth(2)
    heff.SetMarkerSize(0.8)
    heff.SetMarkerStyle(20)
    heff.SetMarkerColor(ROOT.kBlack)

    leg = ROOT.TLegend(0.5, 0.7, 0.9, 0.9)
    leg.AddEntry(heff, legTitle ,"p")
    c = ROOT.TCanvas('c', '', 600, 800)
    fr = c.DrawFrame(xmin, ymin, xmax, ymax)
    fr.GetYaxis().SetTitle(ytitle)
    fr.GetYaxis().SetTitleSize(0.05)
    fr.GetYaxis().SetTitleOffset(0.7)
    fr.GetYaxis().SetLabelSize(0.03)
    fr.GetXaxis().SetTitle(xtitle)
    fr.GetXaxis().SetTitleSize(0.05)
    fr.GetXaxis().SetTitleOffset(0.9)
    fr.GetXaxis().SetLabelSize(0.03)
    heff.Draw("P")
    leg.Draw("SAME")
    c.SaveAs("RecoFixPlots/"+sample+"/"+name+"_"+sample+".png")
    c.Close()
    

def plotStackedEff(hnum, hden, colors, xtitle, name, legTitle, sample, xmin = 0, ymin =0, xmax=50, ymax=1.2, ytitle='Efficiency'):
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
    fr.GetYaxis().SetTitle(ytitle)
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
    c.SaveAs("RecoFixPlots/"+sample+"/"+name+"_"+sample+".png")
    c.Close()

def plotStackedDistr(hnum, hden, colors, xtitle, name, legTitle, sample, xmin = 0, xmax=50, ytitle='Events'):


    colors = [ROOT.kBlue,ROOT.kCyan,ROOT.kRed,ROOT.kOrange,ROOT.kGreen]
    if( not (len(hnum) == len(hden)) or ( len(hnum) > len(colors))):
        print("Error in plotStackedDistr - arrays are incompatible.")
        return
    
    for i in range(len(hnum)):
        hnum[i].SetLineColor(colors[i])
        hnum[i].SetLineWidth(2)
        hnum[i].SetMarkerSize(0.8)
        hnum[i].SetMarkerStyle(20)
        hnum[i].SetMarkerColor(colors[i])

        hden[i].SetLineColor(ROOT.kBlack)
        hden[i].SetLineWidth(2)
        hden[i].SetMarkerSize(0.8)
        hden[i].SetMarkerStyle(20)
        hden[i].SetMarkerColor(ROOT.kBlack)

    histmax = hden[0].GetMaximum()
    ymin = 0
    ymax = histmax / 0.95


    if(len(legTitle) == len(hnum)):
        leg = ROOT.TLegend(0.5, 0.7, 0.9, 0.9)
        leg.AddEntry(hden[0], "Gen electron" ,"p")
        for i in range(len(hnum)):    
            leg.AddEntry(hnum[i], legTitle[i] ,"p")

    c = ROOT.TCanvas('c', '', 600, 800)
    fr = c.DrawFrame(xmin, ymin, xmax, ymax)
    fr.GetYaxis().SetTitle(ytitle)
    fr.GetYaxis().SetTitleSize(0.05)
    fr.GetYaxis().SetTitleOffset(1.0)
    fr.GetYaxis().SetLabelSize(0.03)
    fr.GetXaxis().SetTitle(xtitle)
    fr.GetXaxis().SetTitleSize(0.05)
    fr.GetXaxis().SetTitleOffset(0.9)
    fr.GetXaxis().SetLabelSize(0.03)

    hden[0].Draw("P,SAME")
    for i in range(len(hnum)):
        hnum[i].Draw("P,SAME")

    
    if(len(legTitle) == len(hnum)):    
        leg.Draw("SAME")



    tex = ROOT.TLatex(0.45,0.91,"CMS Simulation work in progress")
    tex.SetNDC()
    tex.SetTextFont(42)
    tex.SetTextSize(0.035)
    tex.SetLineWidth(2)
    tex.Draw()

    c.SaveAs("RecoFixPlots/"+sample+"/"+name+"_"+sample+".png")
    c.Close()









def plotStackedEffAndRatio(hnum, hden, colors, xtitle, name, legTitle, sample, xmin = 0, ymin =-0.2, xmax=50, ymax=1.2, yratiomin = 0.5, yratiomax = 13, ytitle='Efficiency', legendpos='topright'):
    if( len(hnum) > len(hden) or ( len(hnum) > len(colors))):
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
        heff[i].SetTitle("Sample: "+str(sample))

    # Use: 0 is reco fixed, 1 is simple reco
    #xmax = max(ROOT.TMath.MaxElement(hnum.GetN(),hnum.GetX()), ROOT.TMath.MaxElement(hden.GetN(),hden.GetX()))
    lxnum = [x for x in heff[0].GetX()]
    lynum = [x for x in heff[0].GetY()]
    lxden = [x for x in heff[1].GetX()]
    lyden = [x for x in heff[1].GetY()]


    ratioBinning = []
    Nbins = hnum[0].GetNbinsX()
    for i in range (Nbins + 1):
        ratioBinning.append(hnum[0].GetBinLowEdge(i+1))
    #binwidth = hnum[0].GetBinLowEdge(i+1)
    hRatio = ROOT.TH1F("hRatio", "hRatio", Nbins, np.array(ratioBinning))
    #hRatio = ROOT.TH1F("hRatio", "hRatio", Nbins, 0, Nbins*binwidth)
    #hRatio = ROOT.TH1F("hRatio", "hRatio", Nbins, hnum[0].GetBin(1), hnum[0].GetBin(Nbins))

    for i in range (Nbins):
        numerator = hnum[0].GetBinContent(i+1)
        denominator = hnum[1].GetBinContent(i+1)

        if(denominator == 0):
            ratio = -1
        else:
            ratio = numerator / denominator
        hRatio.SetBinContent(i+1, ratio)


    hRatio.GetYaxis().SetTitle("extended / standard")

    #hRatio.GetYaxis().SetRangeUser(yratiomin,yratiomax)
    #if(yratiomax == 0):
    maxim = hRatio.GetMaximum()
    if(maxim < 1):
        maxim = 1.05
    else:
        maxim /= 0.95
    hRatio.GetYaxis().SetRangeUser(yratiomin, maxim)

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

    if(len(legTitle) == len(heff)):
        if(legendpos=='bottomright'):
            leg = ROOT.TLegend(0.5, 0.02, 0.9, 0.22)
        else:
            leg = ROOT.TLegend(0.5, 0.7, 0.9, 0.9)
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
    fr.GetYaxis().SetTitle(ytitle)
    fr.GetYaxis().SetTitleSize(0.05)
    fr.GetYaxis().SetTitleOffset(0.7)
    fr.GetYaxis().SetLabelSize(0.03)
    fr.GetXaxis().SetTitle(xtitle)
    fr.GetXaxis().SetTitleSize(0.05)
    fr.GetXaxis().SetTitleOffset(0.9)
    fr.GetXaxis().SetLabelSize(0.03)
    for i in range(len(heff)):
        heff[i].Draw("P,SAME")

    if(len(legTitle) == len(heff)):
        leg.Draw("SAME")
    #hRatio.Draw("SAME")
    #for i in range(1, hRatioFrame.GetNbinsX() + 1):
    #    print(hRatio.GetBinContent(i))

    tex = ROOT.TLatex(0.45,0.91,"CMS Simulation work in progress")
    tex.SetNDC()
    tex.SetTextFont(42)
    tex.SetTextSize(0.035)
    tex.SetLineWidth(2)
    tex.Draw()

    c.cd()
    p2 = ROOT.TPad("p2", "p2", 0, 0.01, 1, 0.3)
    p2.SetTopMargin(0)
    p2.SetBottomMargin(0.2)
    p2.Draw()
    p2.cd()
    hRatio.Draw("P")
    hRatioFrame.Draw("HISTsame")
    c.SaveAs("RecoFixPlots/"+sample+"/"+name+"_"+sample+".png")
    c.Close()




def plotStackedEffAndRatio2(hnum, hden, colors, xtitle, name, legTitle, sample, xmin = 0, ymin =-0.2, xmax=50, ymax=1.2, yratiomin = 0.5, yratiomax = 13, ytitle='Efficiency', legendpos='topright'):
    if( len(hnum) > len(hden) or ( len(hnum) > len(colors))):
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
        heff[i].SetTitle("Sample: "+str(sample))

    # Use: 0 is reco fixed, 1 is simple reco
    #xmax = max(ROOT.TMath.MaxElement(hnum.GetN(),hnum.GetX()), ROOT.TMath.MaxElement(hden.GetN(),hden.GetX()))
    lxnum = [x for x in heff[0].GetX()]
    lynum = [x for x in heff[0].GetY()]
    lxden = [x for x in heff[1].GetX()]
    lyden = [x for x in heff[1].GetY()]


    ratioBinning = []
    Nbins = hnum[0].GetNbinsX()
    for i in range (Nbins + 1):
        ratioBinning.append(hnum[0].GetBinLowEdge(i+1))
    #binwidth = hnum[0].GetBinLowEdge(i+1)
    hRatio = ROOT.TH1F("hRatio", "hRatio", Nbins, np.array(ratioBinning))
    hRatio2 = ROOT.TH1F("hRatio2", "hRatio2", Nbins, np.array(ratioBinning))
    #hRatio = ROOT.TH1F("hRatio", "hRatio", Nbins, 0, Nbins*binwidth)
    #hRatio = ROOT.TH1F("hRatio", "hRatio", Nbins, hnum[0].GetBin(1), hnum[0].GetBin(Nbins))

    for i in range (Nbins):
        numerator = hnum[0].GetBinContent(i+1)
        denominator = hnum[1].GetBinContent(i+1)

        numerator2 = hnum[2].GetBinContent(i+1)

        if(denominator == 0):
            ratio = -1
            ratio2 = -1
        else:
            ratio = numerator / denominator
            ratio2 = numerator2 / denominator
        hRatio.SetBinContent(i+1, ratio)
        hRatio2.SetBinContent(i+1, ratio2)


    hRatio.GetYaxis().SetTitle("extended / standard")
    hRatio2.GetYaxis().SetTitle("extended / standard")

    #hRatio.GetYaxis().SetRangeUser(yratiomin,yratiomax)
    #if(yratiomax == 0):
    maxim = hRatio.GetMaximum()
    maxim2= hRatio2.GetMaximum()
    if(maxim < maxim2):
        maxim = maxim2 
    if(maxim < 1):
        maxim = 1.05
    else:
        maxim /= 0.95
    hRatio.GetYaxis().SetRangeUser(yratiomin, maxim)

    hRatio.SetTitle("")
    hRatio.SetStats(0)
    hRatio.SetLineColor(colors[0])
    hRatio.SetMarkerColor(colors[0])
    hRatio.SetMarkerSize(0.8)
    hRatio.SetMarkerStyle(20)
    hRatio.GetYaxis().SetTitleSize(0.08)
    hRatio.GetYaxis().SetTitleOffset(0.5)
    hRatio.GetYaxis().SetLabelSize(0.07)
    hRatio.GetXaxis().SetTitle(xtitle)
    hRatio.GetXaxis().SetTitleSize(0.1)
    hRatio.GetXaxis().SetTitleOffset(0.9)
    hRatio.GetXaxis().SetLabelSize(0.07)

    hRatio2.SetTitle("")
    hRatio2.SetStats(0)
    hRatio2.SetLineColor(colors[2])
    hRatio2.SetMarkerColor(colors[2])
    hRatio2.SetMarkerSize(0.8)
    hRatio2.SetMarkerStyle(20)
    hRatio2.GetYaxis().SetTitleSize(0.08)
    hRatio2.GetYaxis().SetTitleOffset(0.5)
    hRatio2.GetYaxis().SetLabelSize(0.07)
    hRatio2.GetXaxis().SetTitle(xtitle)
    hRatio2.GetXaxis().SetTitleSize(0.1)
    hRatio2.GetXaxis().SetTitleOffset(0.9)
    hRatio2.GetXaxis().SetLabelSize(0.07)

    hRatioFrame = hRatio.Clone("RatioFrame")
    for b in range(1, hRatioFrame.GetNbinsX() + 1):
        hRatioFrame.SetBinContent(b, 1.0)
    hRatioFrame.SetLineStyle(2)
    hRatioFrame.SetLineWidth(2)
    hRatioFrame.SetLineColor(ROOT.kGreen)

    if(len(legTitle) == len(heff)):
        if(legendpos=='bottomright'):
            leg = ROOT.TLegend(0.5, 0.02, 0.9, 0.22)
        else:
            leg = ROOT.TLegend(0.5, 0.7, 0.9, 0.9)
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
    fr.GetYaxis().SetTitle(ytitle)
    fr.GetYaxis().SetTitleSize(0.05)
    fr.GetYaxis().SetTitleOffset(0.7)
    fr.GetYaxis().SetLabelSize(0.03)
    fr.GetXaxis().SetTitle(xtitle)
    fr.GetXaxis().SetTitleSize(0.05)
    fr.GetXaxis().SetTitleOffset(0.9)
    fr.GetXaxis().SetLabelSize(0.03)
    for i in range(len(heff)):
        heff[i].Draw("P,SAME")

    if(len(legTitle) == len(heff)):
        leg.Draw("SAME")
    #hRatio.Draw("SAME")
    #for i in range(1, hRatioFrame.GetNbinsX() + 1):
    #    print(hRatio.GetBinContent(i))

    tex = ROOT.TLatex(0.45,0.91,"CMS Simulation work in progress")
    tex.SetNDC()
    tex.SetTextFont(42)
    tex.SetTextSize(0.035)
    tex.SetLineWidth(2)
    tex.Draw()

    c.cd()
    p2 = ROOT.TPad("p2", "p2", 0, 0.01, 1, 0.3)
    p2.SetTopMargin(0)
    p2.SetBottomMargin(0.2)
    p2.Draw()
    p2.cd()
    hRatio.Draw("P")
    hRatio2.Draw("P,SAME")
    hRatioFrame.Draw("HISTsame")
    c.SaveAs("RecoFixPlots/"+sample+"/"+name+"_"+sample+".png")
    c.Close()





if(bThreaded):
    f = ROOT.TFile.Open('root_files/RecoFixThreadedSTACK_Sample'+sample+'.root')
else:   
    f = ROOT.TFile.Open('root_files/RecoFix_Sample'+sample+'.root')

if ("Electron" == channel or "All" == channel):
    '''
    h_pass = f.Get("RecoElePt_")
    h_total = f.Get("TrueElePt_")
    plotEff(h_pass, h_total, 'ele p_{T} [GeV]', 'electron_pT-eff', "Reco ele", sample)
    '''
    hpassv = {}
    hpassv[0] = f.Get("RecoFixedElePt_")
    hpassv[1] = f.Get("RecoElePt_")
    hpassv[2] = f.Get("RecoPhotonPlusIsoTrackPt_")
    htotv = {}
    htotv[0] = f.Get("TrueElePt_")
    htotv[1] = f.Get("TrueElePt_")
    htotv[2] = f.Get("TrueElePt_")

    colors = [ROOT.kBlue,ROOT.kBlack,ROOT.kRed,ROOT.kOrange,ROOT.kGreen]
    legTitle = ["Extended reco ele","Standard reco ele","matched isoTrack","Reco photon","Reco pho + isoTrack"]

    plotStackedEffAndRatio(hpassv, htotv, colors, 'p_{T} [GeV]', 'TESTING_pT-eff', legTitle, sample,0,-0.05,100,1.2,-1.1,2)

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
    #h_pass = f.Get("RecoPhotonPt_")
    #h_total = f.Get("TrueElePt_")
    #Plot1D(h_pass ,'./RecoFixPlots')
    #Plot1D(h_total ,'./RecoFixPlots')
    #plotEff(h_pass, h_total, 'photon p_{T} [GeV]', 'photon_pT-eff', "Reco photon", sample)

    #h_pass = f.Get("RecoPhotonEta_")
    #h_total = f.Get("TrueEleEta_")
    #plotEff(h_pass, h_total, 'photon eta [GeV]', 'photon_eta-eff', "Reco photon", sample,-3,0,3)

    Plot1D(f.Get("AllTruePhotonPt_") ,'./RecoFixPlots/'+sample)
    Plot1D(f.Get("RecoPhotonPt_") ,'./RecoFixPlots/'+sample)
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
    #colors = [ROOT.kBlue,ROOT.kBlack,ROOT.kRed,ROOT.kOrange,ROOT.kGreen,ROOT.kAzure + 10]
    colors = [ROOT.kBlue,ROOT.kBlack,ROOT.kAzure + 10]
    #legTitle = ["Extended reco ele","Standard reco ele","matched isoTrack","Reco photon","Reco pho + isoTrack","Extended+Exrapol"]
    legTitle = ["Extended reco ele","Standard reco ele","Extended+Exrapol"]
    h_pass[0] = f.Get("RecoFixedElePt_")
    h_pass[1] = f.Get("RecoElePt_")
    #h_pass[2] = f.Get("IsoTrackPt_")
    #h_pass[3] = f.Get("RecoPhotonPt_")
    #h_pass[4] = f.Get("RecoPhotonPlusIsoTrackPt_")
    h_pass[2] = f.Get("RecoFixedExtrapolElePt_")

    h_total[0] = f.Get("TrueElePt_")
    h_total[1] = f.Get("TrueElePt_")
    h_total[2] = f.Get("TrueElePt_")
    #h_total[3] = f.Get("TrueElePt_")
    #h_total[4] = f.Get("TrueElePt_")
    #h_total[5] = f.Get("TrueElePt_")


    plotStackedEffAndRatio2(h_pass, h_total, colors, 'p_{T} [GeV]', 'electron_common_pT-eff', legTitle, sample,0,-0.05,100,1.2,0.9,2)
    #plotStackedEff(h_pass, h_total, colors, 'p_{T} [GeV]', 'electron_common_pT-eff', legTitle, sample,0,-0.05,100,1.2)
    plotStackedDistr(h_pass, h_total, colors, 'p_{T} [GeV]', 'electron_common_pT-distr', legTitle, sample,0,40)


    Plot2DEff( f.Get("RecoFixedExtrapolPtDxy2D_"), f.Get("TrueElePtDxy2D_") ,'./RecoFixPlots/'+sample,xmin = 0, ymin =0, xmax=100, ymax=15, drawOption="COLZ,text",islogz=True)
    Plot2DEff( f.Get("RecoFixedPtDxy2D_"), f.Get("TrueElePtDxy2D_") ,'./RecoFixPlots/'+sample,xmin = 0, ymin =0, xmax=100, ymax=15, drawOption="COLZ,text",islogz=True)
    
    h_pass = {}
    h_total = {}
    h_pass[0] = f.Get("RecoFixedEleEta_")
    h_pass[1] = f.Get("RecoEleEta_")
    #h_pass[2] = f.Get("IsoTrackEta_")
    #h_pass[3] = f.Get("RecoPhotonEta_")
    #h_pass[4] = f.Get("RecoPhotonPlusIsoTrackEta_")
    h_pass[2] = f.Get("RecoFixedExtrapolEleEta_")

    h_total[0] = f.Get("TrueEleEta_")
    h_total[1] = f.Get("TrueEleEta_")
    h_total[2] = f.Get("TrueEleEta_")
    #h_total[3] = f.Get("TrueEleEta_")
    #h_total[4] = f.Get("TrueEleEta_")

    plotStackedEffAndRatio2(h_pass, h_total, colors, '\eta', 'electron_common_eta-eff', legTitle, sample,-3,-0.05,3,1.2,0.9,2)
    #plotStackedEff(h_pass, h_total, colors, '\eta', 'electron_common_eta-eff', legTitle, sample,-3,-0.05,3,1.2)
    plotStackedDistr(h_pass, h_total, colors, '\eta', 'electron_common_eta-distr', legTitle, sample,-3,3)

    h_pass = {}
    h_total = {}
    h_pass[0] = f.Get("RecoFixedElePhi_")
    h_pass[1] = f.Get("RecoElePhi_")
    #h_pass[2] = f.Get("IsoTrackPhi_")
    #h_pass[3] = f.Get("RecoPhotonPhi_")
    #h_pass[4] = f.Get("RecoPhotonPlusIsoTrackPhi_")
    h_pass[2] = f.Get("RecoFixedExtrapolElePhi_")

    h_total[0] = f.Get("TrueElePhi_")
    h_total[1] = f.Get("TrueElePhi_")
    h_total[2] = f.Get("TrueElePhi_")
    #h_total[3] = f.Get("TrueElePhi_")
    #h_total[4] = f.Get("TrueElePhi_")

    plotStackedEffAndRatio2(h_pass, h_total, colors, '\phi', 'electron_common_phi-eff', legTitle, sample,-4,-0.05,4,1.2,0.9,2)


    legTitle = []
    plotStackedDistr(h_pass, h_total, colors, '\phi', 'electron_common_phi-distr', legTitle, sample,-4,4)


    h_pass = {}
    h_total = {}
    colors = [ROOT.kBlue,ROOT.kBlack,ROOT.kRed]
    legTitle = ["Extended reco ele","Standard reco ele","matched isoTrack"]
    h_pass[0] = f.Get("RecoFixedEleVtxx_")
    h_pass[1] = f.Get("RecoEleVtxx_")
    h_pass[2] = f.Get("IsoTrackVtxx_")

    h_total[0] = f.Get("TrueEleVtxx_")
    h_total[1] = f.Get("TrueEleVtxx_")
    h_total[2] = f.Get("TrueEleVtxx_")

    plotStackedEff(h_pass, h_total, colors, 'vtx x [cm]', 'electron_common_vtxx-eff', legTitle, sample,-3,0,3)


    h_pass = {}
    h_total = {}
    colors = [ROOT.kBlue,ROOT.kBlack,ROOT.kRed]
    legTitle = ["Extended reco ele","Standard reco ele","matched isoTrack"]
    h_pass[0] = f.Get("RecoFixedEleVtxy_")
    h_pass[1] = f.Get("RecoEleVtxy_")
    h_pass[2] = f.Get("IsoTrackVtxy_")

    h_total[0] = f.Get("TrueEleVtxy_")
    h_total[1] = f.Get("TrueEleVtxy_")
    h_total[2] = f.Get("TrueEleVtxy_")

    plotStackedEff(h_pass, h_total, colors, 'vtx y [cm]', 'electron_common_vtxy-eff', legTitle, sample,-3,0,3)

    h_pass = {}
    h_total = {}
    colors = [ROOT.kBlue,ROOT.kBlack,ROOT.kRed]
    legTitle = ["Extended reco ele","Standard reco ele","matched isoTrack"]
    h_pass[0] = f.Get("RecoFixedEleVtxz_")
    h_pass[1] = f.Get("RecoEleVtxz_")
    h_pass[2] = f.Get("IsoTrackVtxz_")

    h_total[0] = f.Get("TrueEleVtxz_")
    h_total[1] = f.Get("TrueEleVtxz_")
    h_total[2] = f.Get("TrueEleVtxz_")

    plotStackedEff(h_pass, h_total, colors, 'vtx z [cm]', 'electron_common_vtxz-eff', legTitle, sample,-3,0,3)


    h_pass = {}
    h_total = {}
    colors = [ROOT.kBlue,ROOT.kBlack,ROOT.kRed]
    legTitle = ["Extended reco ele","Standard reco ele","matched isoTrack"]
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
    legTitle = ["Extended reco ele","Standard reco ele","matched isoTrack"]
    h_pass[0] = f.Get("RecoFixedEleDz_")
    h_pass[1] = f.Get("RecoEleDz_")
    h_pass[2] = f.Get("IsoTrackDz_")

    h_total[0] = f.Get("TrueEleDz_")
    h_total[1] = f.Get("TrueEleDz_")
    h_total[2] = f.Get("TrueEleDz_")

    plotStackedEff(h_pass, h_total, colors, 'dz [cm]', 'electron_common_dz-eff', legTitle, sample,-3,0,3)


    #colors = [ROOT.kBlue,ROOT.kBlack,ROOT.kRed,ROOT.kOrange,ROOT.kGreen,ROOT.kAzure + 10]
    colors = [ROOT.kBlue,ROOT.kBlack,ROOT.kAzure + 10]
    #legTitle = ["Extended reco ele","Standard reco ele","matched isoTrack","Reco photon","Reco pho + isoTrack","Extended+Exrapol"]
    legTitle = ["Extended reco ele","Standard reco ele","Extended+Exrapol"]



    h_pass = {}
    h_total = {}
    #colors = [ROOT.kBlue,ROOT.kBlack,ROOT.kRed,ROOT.kOrange,ROOT.kGreen]
    #legTitle = ["Extended reco ele","Standard reco ele","matched isoTrack","Reco photon","Reco pho + isoTrack"]
    h_pass[0] = f.Get("RecoFixedEleDzAbs_")
    h_pass[1] = f.Get("RecoEleDzAbs_")
    #h_pass[2] = f.Get("IsoTrackDzAbs_")
    #h_pass[3] = f.Get("RecoPhotonDzAbs_")
    #h_pass[4] = f.Get("RecoPhotonPlusIsoTrackDzAbs_")
    h_pass[2] = f.Get("RecoFixedExtrapolEleDzAbs_")

    h_total[0] = f.Get("TrueEleDzAbs_")
    h_total[1] = f.Get("TrueEleDzAbs_")
    h_total[2] = f.Get("TrueEleDzAbs_")
    #h_total[3] = f.Get("TrueEleDzAbs_")
    #h_total[4] = f.Get("TrueEleDzAbs_")

    plotStackedEffAndRatio2(h_pass, h_total, colors, 'abs(dz) [cm]', 'electron_common_dzAbs-eff', legTitle, sample,0,-0.05,15,1.2,0.9,2)

    h_pass = {}
    h_total = {}
    h_pass[0] = f.Get("RecoFixedEleDxyAbs_")
    h_pass[1] = f.Get("RecoEleDxyAbs_")
    #h_pass[2] = f.Get("IsoTrackDxyAbs_")
    #h_pass[3] = f.Get("RecoPhotonDxyAbs_")
    #h_pass[4] = f.Get("RecoPhotonPlusIsoTrackDxyAbs_")
    h_pass[2] = f.Get("RecoFixedExtrapolEleDzAbs_")

    h_total[0] = f.Get("TrueEleDxyAbs_")
    h_total[1] = f.Get("TrueEleDxyAbs_")
    h_total[2] = f.Get("TrueEleDxyAbs_")
    #h_total[3] = f.Get("TrueEleDxyAbs_")
    #h_total[4] = f.Get("TrueEleDxyAbs_")

    plotStackedEffAndRatio2(h_pass, h_total, colors, 'abs(dxy) [cm]', 'electron_common_dxyAbs-eff', legTitle, sample,0,-0.05,15)


    Plot2D( f.Get("dRControlPlot2D_") ,'./RecoFixPlots',xmin = 0, ymin =0, xmax=1, ymax=1, drawOption="COLZ",islogz=True)

    # background
    h_pass = {}
    h_total = {}
    legTitle = ["Extended reco ele","Standard reco ele"]
    h_pass[0] = f.Get("RejectBackgroundRecoFixPt_")
    h_pass[1] = f.Get("RejectBackgroundPt_")

    h_total[0] = f.Get("TrueBackgroundPt_")
    h_total[1] = f.Get("TrueBackgroundPt_")

    plotStackedEffAndRatio(h_pass, h_total, colors, 'p_{T} [GeV]', 'background_common_pT', legTitle, sample,0,-0.05,100,1.2,0,1.2,ytitle='Rejection (1-Eff)',legendpos='bottomright')

    h_pass = {}
    h_total = {}
    legTitle = ["Extended reco ele","Standard reco ele"]
    h_pass[0] = f.Get("RejectBackgroundRecoFixEta_")
    h_pass[1] = f.Get("RejectBackgroundEta_")

    h_total[0] = f.Get("TrueBackgroundEta_")
    h_total[1] = f.Get("TrueBackgroundEta_")


    plotStackedEffAndRatio(h_pass, h_total, colors, '\eta ', 'background_common_eta', legTitle, sample,-3,-0.05,3,1.2,0,1.2,ytitle='Rejection (1-Eff)',legendpos='bottomright')
    #plotStackedEff(h_pass, h_total, colors, '\eta ', 'background_common_eta', legTitle, sample,-3,0,3)

    h_pass = {}
    h_total = {}
    legTitle = ["Extended reco ele","Standard reco ele"]
    h_pass[0] = f.Get("RejectBackgroundRecoFixPhi_")
    h_pass[1] = f.Get("RejectBackgroundPhi_")

    h_total[0] = f.Get("TrueBackgroundPhi_")
    h_total[1] = f.Get("TrueBackgroundPhi_")


    plotStackedEffAndRatio(h_pass, h_total, colors, '\phi ', 'background_common_phi', legTitle, sample,-4,-0.05,4,1.2,0,1.2,ytitle='Rejection (1-Eff)',legendpos='bottomright')
    #plotStackedEff(h_pass, h_total, colors, '\eta ', 'background_common_eta', legTitle, sample,-3,0,3)

    h_pass = {}
    h_total = {}
    legTitle = ["Extended reco ele","Standard reco ele"]
    h_pass[0] = f.Get("RejectBackgroundRecoFixDxyAbs_")
    h_pass[1] = f.Get("RejectBackgroundDxyAbs_")

    h_total[0] = f.Get("TrueBackgroundDxyAbs_")
    h_total[1] = f.Get("TrueBackgroundDxyAbs_")

    plotStackedEffAndRatio(h_pass, h_total, colors, 'abs(dxy) [cm] ', 'background_common_dxyAbs', legTitle, sample, 0,-0.05,15,1.2,0,1.2,ytitle='Rejection (1-Eff)',legendpos='bottomright')

    h_pass = {}
    h_total = {}
    legTitle = ["Extended reco ele","Standard reco ele"]
    h_pass[0] = f.Get("RejectBackgroundRecoFixDzAbs_")
    h_pass[1] = f.Get("RejectBackgroundDzAbs_")

    h_total[0] = f.Get("TrueBackgroundDzAbs_")
    h_total[1] = f.Get("TrueBackgroundDzAbs_")

    plotStackedEffAndRatio(h_pass, h_total, colors, 'abs(dz) [cm]', 'background_common_dzAbs', legTitle, sample, 0,-0.05,15,1.2,0,1.2,ytitle='Rejection (1-Eff)',legendpos='bottomright')


    # ROC


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
