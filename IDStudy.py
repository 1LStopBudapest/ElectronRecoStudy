import os, sys
import ROOT
import types
import math


#from RegSel import RegSel

sys.path.append('../')
from Helper.HistInfo import HistInfo
from Helper.MCWeight import MCWeight
from Helper.TreeVarSel import TreeVarSel
from TriggerStudy.TrigVarSel import TrigVarSel
from Sample.SampleChain import SampleChain
from Helper.VarCalc import *
from Helper.PlotHelper import *
from Sample.FileList_2016 import samples as samples_2016

def get_parser():
    ''' Argument parser.
    '''
    import argparse
    argParser = argparse.ArgumentParser(description = "Argument parser")
    argParser.add_argument('--sample',           action='store',                     type=str,            default='FullSim100',                                help="Which sample?" )
    argParser.add_argument('--year',             action='store',                     type=int,            default=2016,                                             help="Which year?" )
    argParser.add_argument('--startfile',        action='store',                     type=int,            default=0,                                                help="start from which root file like 0th or 10th etc?" )
    argParser.add_argument('--nfiles',           action='store',                     type=int,            default=-1,                                               help="No of files to run. -1 means all files" )
    argParser.add_argument('--nevents',           action='store',                    type=int,            default=-1,                                               help="No of events to run. -1 means all events" )
    argParser.add_argument('--channel',           action='store',                    type=str,            default='Electron',                                       help="Which lepton?" )
    argParser.add_argument('--region',            action='store',                    type=str,            default='mesurement',                                     help="Which lepton?" )    
    argParser.add_argument('--pJobs',             action='store',                    type=bool,            default=False,                                           help="using GPU parallel program or not" )

    return argParser

options = get_parser().parse_args()



samples  = options.sample
channel = options.channel
nEvents = options.nevents
year = options.year

isData = True if ('Run' in samples or 'Data' in samples) else False

lepOpt = 'Ele' if 'Electron' in channel else 'Mu'

DataLumi = 1.0

if year==2016:
    samplelist = samples_2016
    DataLumi = SampleChain.luminosity_2016
elif year==2017:
    samplelist = samples_2017
    DataLumi = SampleChain.luminosity_2017
else:
    samplelist = samples_2018
    DataLumi = SampleChain.luminosity_2018


# csak fullsimre
if 'FullSim0' != samples and 'FullSim5' != samples and 'FullSim10' != samples and 'FullSim100' != samples and 'FullSim100m' != samples:
    print("Sorry for now, may be expanded later.")
    quit()

sample = samples

# Make samplechain
ch = SampleChain(sample, options.startfile, options.nfiles, year).getchain()

def dR(eta_gen, eta_reco, phi_gen, phi_reco):
    return math.sqrt((eta_gen-eta_reco)**2 + (phi_gen-phi_reco)**2 )

def dRcut(dist):
    return dist < 0.1

def hasMomRecursive(i, pdgid, ch):
    if( ch.GenPart_genPartIdxMother[i] == -1):
        return False
    elif(abs(ch.GenPart_pdgId[ch.GenPart_genPartIdxMother[i]]) == pdgid):
        return True
    else:
        return hasMomRecursive( ch.GenPart_genPartIdxMother[i], pdgid, ch )



# HISTOS #
# Save histos:
hfile = ROOT.TFile( 'root_files/IDStudy_Sample'+sample+'.root', 'RECREATE')
ptBinning = [0, 3.5, 5, 12, 20, 30, 50, 80, 200]
histext = ''
histos = {}
histos['TrueElePt'] = HistInfo(hname = 'TrueElePt', sample = histext, binning = [50, 0, 100], histclass = ROOT.TH1F).make_hist()
histos['NofilterElePt'] = HistInfo(hname = 'NofilterElePt', sample = histext, binning = [50, 0, 100], histclass = ROOT.TH1F).make_hist()
histos['VetoElePt'] = HistInfo(hname = 'VetoElePt', sample = histext, binning = [50, 0, 100], histclass = ROOT.TH1F).make_hist()
histos['LooseElePt'] = HistInfo(hname = 'LooseElePt', sample = histext, binning = [50, 0, 100], histclass = ROOT.TH1F).make_hist()
histos['MediumElePt'] = HistInfo(hname = 'MediumElePt', sample = histext, binning = [50, 0, 100], histclass = ROOT.TH1F).make_hist()
histos['TightElePt'] = HistInfo(hname = 'TightElePt', sample = histext, binning = [50, 0, 100], histclass = ROOT.TH1F).make_hist()

histos['TrueMuonPt'] = HistInfo(hname = 'TrueMuonPt', sample = histext, binning = [50, 0, 100], histclass = ROOT.TH1F).make_hist()
histos['NofilterMuonPt'] = HistInfo(hname = 'NofilterMuonPt', sample = histext, binning = [50, 0, 100], histclass = ROOT.TH1F).make_hist()
#histos['VetoMuonPt'] = HistInfo(hname = 'VetoMuonPt', sample = histext, binning = [50, 0, 100], histclass = ROOT.TH1F).make_hist()
histos['LooseMuonPt'] = HistInfo(hname = 'LooseMuonPt', sample = histext, binning = [50, 0, 100], histclass = ROOT.TH1F).make_hist()
histos['MediumMuonPt'] = HistInfo(hname = 'MediumMuonPt', sample = histext, binning = [50, 0, 100], histclass = ROOT.TH1F).make_hist()
histos['TightMuonPt'] = HistInfo(hname = 'TightMuonPt', sample = histext, binning = [50, 0, 100], histclass = ROOT.TH1F).make_hist()


histos['TrueEleEta'] = HistInfo(hname = 'TrueEleEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['NofilterEleEta'] = HistInfo(hname = 'NofilterEleEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['VetoEleEta'] = HistInfo(hname = 'VetoEleEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['LooseEleEta'] = HistInfo(hname = 'LooseEleEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['MediumEleEta'] = HistInfo(hname = 'MediumEleEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['TightEleEta'] = HistInfo(hname = 'TightEleEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()

histos['TrueMuonEta'] = HistInfo(hname = 'TrueMuonEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['NofilterMuonEta'] = HistInfo(hname = 'NofilterMuonEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
#histos['VetoMuonEta'] = HistInfo(hname = 'VetoMuonEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['LooseMuonEta'] = HistInfo(hname = 'LooseMuonEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['MediumMuonEta'] = HistInfo(hname = 'MediumMuonEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['TightMuonEta'] = HistInfo(hname = 'TightMuonEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()

histos['IsoTrackEta'] = HistInfo(hname = 'IsoTrackEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['IsoTrackPt'] = HistInfo(hname = 'IsoTrackPt', sample = histext, binning = [50, 0, 100], histclass = ROOT.TH1F).make_hist()

histos['TruePhotonPt'] = HistInfo(hname = 'TruePhotonPt', sample = histext, binning = [50, 0, 100], histclass = ROOT.TH1F).make_hist()
histos['RecoPhotonPt'] = HistInfo(hname = 'RecoPhotonPt', sample = histext, binning = [50, 0, 100], histclass = ROOT.TH1F).make_hist()
histos['TruePhotonEta'] = HistInfo(hname = 'TruePhotonEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['RecoPhotonEta'] = HistInfo(hname = 'RecoPhotonEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
effs = {}


n_entries = ch.GetEntries() #num entries in tree
nevtcut = n_entries -1 if nEvents == - 1 else nEvents - 1

# RECO #
print("Calculating efficiencies with different IDs.")
for ientry in range(n_entries): #loop over events
    if ientry > nevtcut: break
    if ientry % (nevtcut/10)==0 : print 'processing ', ientry,'th event'
    ch.GetEntry(ientry)

    for i in range(ch.nGenPart):
    
        # electron, pt eta cuts, status 1, with a stop mother somewhere in mother history
        if( abs(ch.GenPart_pdgId[i]) == 11 and ch.GenPart_pt[i] >= 5 and abs(ch.GenPart_eta[i]) <= 2.5 and ch.GenPart_status[i] == 1 and hasMomRecursive(i, 1000006, ch)):
            histos['TrueElePt'].Fill(ch.GenPart_pt[i])
            histos['TrueEleEta'].Fill(ch.GenPart_eta[i])

            # check if there is a reco ele
            eledist = [100000,100000,100000,100000,100000]
            eleIdx = [-1,-1,-1,-1,-1]
            idx0 = [-1,-1,-1,-1,-1]
            for j in range(ch.nElectron):
                #eleVID(ch.Electron_vidNestedWPBitmap[j], 0, removedCuts=['pfRelIso03_all'])
                eleVID(ch.Electron_vidNestedWPBitmap[j], 1, removedCuts=['pfRelIso03_all'])
                eleVID(ch.Electron_vidNestedWPBitmap[j], 2, removedCuts=['pfRelIso03_all'])
                eleVID(ch.Electron_vidNestedWPBitmap[j], 3, removedCuts=['pfRelIso03_all'])
                eleVID(ch.Electron_vidNestedWPBitmap[j], 4, removedCuts=['pfRelIso03_all'])

                dist0 = dR(ch.GenPart_eta[i], ch.Electron_eta[j], ch.GenPart_phi[i], ch.Electron_phi[j] ) 
                for iID in range(len(eledist)):
                    #if( dist0 < eledist[iID] and ch.Electron_cutBased_Fall17_V1[j] >= iID):
                    if dist0 < eledist[iID] and eleVID(ch.Electron_vidNestedWPBitmap[j], iID, removedCuts=['pfRelIso03_all']):
                        eledist[iID] = dist0
                        idx0[iID] = j
            for iID in range(len(eledist)):
                if( dRcut(eledist[iID]) ):
                    if   0 == iID:
                        histos['NofilterElePt'].Fill(ch.GenPart_pt[i])
                        histos['NofilterEleEta'].Fill(ch.GenPart_eta[i])
                    elif 1 == iID:
                        histos['VetoElePt'].Fill(ch.GenPart_pt[i])
                        histos['VetoEleEta'].Fill(ch.GenPart_eta[i])
                    elif 2 == iID:
                        histos['LooseElePt'].Fill(ch.GenPart_pt[i])
                        histos['LooseEleEta'].Fill(ch.GenPart_eta[i])
                    elif 3 == iID:
                        histos['MediumElePt'].Fill(ch.GenPart_pt[i])
                        histos['MediumEleEta'].Fill(ch.GenPart_eta[i])
                    elif 4 == iID:
                        histos['TightElePt'].Fill(ch.GenPart_pt[i])
                        histos['TightEleEta'].Fill(ch.GenPart_eta[i])
                    else:
                        print("Itt valami nincs rendben. (Ele)")
                    eleIdx[iID] = idx0[iID]

            # match isotracks
            isodist = 100000
            isoIdx = -1
            idx0 = -1
            for j in range(ch.nIsoTrack):
                dist0 = dR(ch.GenPart_eta[i], ch.IsoTrack_eta[j], ch.GenPart_phi[i], ch.IsoTrack_phi[j] ) 
                if( dist0 < isodist):
                    isodist = dist0
                    idx0 = j
            if( dRcut(isodist) ):
                histos['IsoTrackPt'].Fill(ch.GenPart_pt[i])
                histos['IsoTrackEta'].Fill(ch.GenPart_eta[i])
                isoIdx = idx0

        # muon, pt eta cuts, status 1, with a stop mother somewhere in mother history
        if( abs(ch.GenPart_pdgId[i]) == 13 and ch.GenPart_pt[i] >= 3.5 and abs(ch.GenPart_eta[i]) <= 2.4 and ch.GenPart_status[i] == 1 and hasMomRecursive(i, 1000006, ch)):
            histos['TrueMuonPt'].Fill(ch.GenPart_pt[i])
            histos['TrueMuonEta'].Fill(ch.GenPart_eta[i])

            # check if there is a reco muon
            muondist = [100000,100000,100000,100000]
            muonIdx = [-1,-1,-1,-1]
            idx0 = [-1,-1,-1,-1]
            for j in range(ch.nMuon):
                dist0 = dR(ch.GenPart_eta[i], ch.Muon_eta[j], ch.GenPart_phi[i], ch.Muon_phi[j] ) 

                if( dist0 < muondist[0]):
                    muondist[0] = dist0
                    idx0[0] = j
                if( dist0 < muondist[1] and ch.Muon_looseId[j]):
                    muondist[1] = dist0
                    idx0[1] = j
                if( dist0 < muondist[2] and ch.Muon_mediumId[j]):
                    muondist[2] = dist0
                    idx0[2] = j
                if( dist0 < muondist[3] and ch.Muon_tightId[j]):
                    muondist[3] = dist0
                    idx0[3] = j

            for iID in range(len(muondist)):
                if( dRcut(muondist[iID]) ):
                    if   0 == iID:
                        histos['NofilterMuonPt'].Fill(ch.GenPart_pt[i])
                        histos['NofilterMuonEta'].Fill(ch.GenPart_eta[i])
                    elif 1 == iID:
                        histos['LooseMuonPt'].Fill(ch.GenPart_pt[i])
                        histos['LooseMuonEta'].Fill(ch.GenPart_eta[i])
                    elif 2 == iID:
                        histos['MediumMuonPt'].Fill(ch.GenPart_pt[i])
                        histos['MediumMuonEta'].Fill(ch.GenPart_eta[i])
                    elif 3 == iID:
                        histos['TightMuonPt'].Fill(ch.GenPart_pt[i])
                        histos['TightMuonEta'].Fill(ch.GenPart_eta[i])
                    else:
                        print("Itt valami nincs rendben (Muon).")
                    muonIdx[iID] = idx0[iID]

        # photon    
        if( abs(ch.GenPart_pdgId[i]) == 22 ): # true photon
            histos['TruePhotonPt'].Fill(ch.GenPart_pt[i])
            histos['TruePhotonEta'].Fill(ch.GenPart_eta[i])
            
            dist = 100000
            phoIdx = -1
            for j in range(ch.nPhoton):
                dist0 = dR(ch.GenPart_eta[i], ch.Photon_eta[j], ch.GenPart_phi[i], ch.Photon_phi[j] ) 
                if( dist0 < dist):
                    dist = dist0
                    phoIdx = j
            if( dist < 0.1):
                histos['RecoPhotonPt'].Fill(ch.GenPart_pt[i])
                histos['RecoPhotonEta'].Fill(ch.GenPart_eta[i])


# Save histos:
hfile.Write()
