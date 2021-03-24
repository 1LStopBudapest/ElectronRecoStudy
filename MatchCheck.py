import os, sys
import ROOT
import types
import math


#from RegSel import RegSel

sys.path.append('../')
from Helper.HistInfo import HistInfo
from Helper.MCWeight import MCWeight
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


def hasMom(i, pdgid, ch):
    if( ch.GenPart_genPartIdxMother[i] == -1):
        return False
    elif(abs(ch.GenPart_pdgId[ch.GenPart_genPartIdxMother[i]]) == pdgid):
        return True
    return False




# HISTOS #
# Save histos:
hfile = ROOT.TFile( 'root_files/MatchCheck_Sample'+sample+'.root', 'RECREATE')
ptBinning = [0, 3.5, 5, 12, 20, 30, 50, 80, 200]
histext = ''
histos = {}
histos['TrueElePt'] = HistInfo(hname = 'TrueElePt', sample = histext, binning = [50, 0, 100], histclass = ROOT.TH1F).make_hist()
histos['RecoElePt'] = HistInfo(hname = 'RecoElePt', sample = histext, binning = [50, 0, 100], histclass = ROOT.TH1F).make_hist()
histos['IsoTrackPt'] = HistInfo(hname = 'IsoTrackPt', sample = histext, binning = [50, 0, 100], histclass = ROOT.TH1F).make_hist()
histos['RecoFixedElePt'] = HistInfo(hname = 'RecoFixedElePt', sample = histext, binning = [50, 0, 100], histclass = ROOT.TH1F).make_hist()
histos['TrueEleEta'] = HistInfo(hname = 'TrueEleEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['RecoEleEta'] = HistInfo(hname = 'RecoEleEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['IsoTrackEta'] = HistInfo(hname = 'IsoTrackEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['RecoFixedEleEta'] = HistInfo(hname = 'RecoFixedEleEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()

histos['TruePhotonPt'] = HistInfo(hname = 'TruePhotonPt', sample = histext, binning = [50, 0, 100], histclass = ROOT.TH1F).make_hist()
histos['RecoPhotonPt'] = HistInfo(hname = 'RecoPhotonPt', sample = histext, binning = [50, 0, 100], histclass = ROOT.TH1F).make_hist()
histos['TruePhotonEta'] = HistInfo(hname = 'TruePhotonEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['RecoPhotonEta'] = HistInfo(hname = 'RecoPhotonEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()


histos['AllPhotonEta'] = HistInfo(hname = 'AllPhotonEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['MatchedPhotonEta'] = HistInfo(hname = 'MatchedPhotonEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['MatchedPhotonGenEta'] = HistInfo(hname = 'MatchedPhotonGenEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['AllIsoTrackEta'] = HistInfo(hname = 'AllIsoTrackEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['MatchedIsoTrackEta'] = HistInfo(hname = 'MatchedIsoTrackEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['MatchedIsoTrackGenEta'] = HistInfo(hname = 'MatchedIsoTrackGenEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['AllElectronEta'] = HistInfo(hname = 'AllElectronEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['MatchedElectronEta'] = HistInfo(hname = 'MatchedElectronEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['MatchedElectronGenEta'] = HistInfo(hname = 'MatchedElectronGenEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['DoubleMatchedIsoTrackEta'] = HistInfo(hname = 'DoubleMatchedIsoTrackEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['DoubleMatchedIsoTrackGenEta'] = HistInfo(hname = 'DoubleMatchedIsoTrackGenEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
effs = {}


n_entries = ch.GetEntries() #num entries in tree
nevtcut = n_entries -1 if nEvents == - 1 else nEvents - 1

# ele has larger impact parameter
# reco eff is low
# try improving

'''
AllCutEle = HistInfo(hname = 'All Ele Pt', sample = histext, binning = [50, 0, 100], histclass = ROOT.TH1F).make_hist()
AllCutEleGoodMom = HistInfo(hname = 'All Ele with >= 0 momID Pt', sample = histext, binning = [50, 0, 100], histclass = ROOT.TH1F).make_hist()
EleMomDist = HistInfo(hname = 'momID', sample = histext, binning = [30, 0, 30], histclass = ROOT.TH1F).make_hist()
finds = {}
ifi = 0
# MOMPDGID TESTING
print("Calculating efficiencies with and without recofix.")
for ientry in range(n_entries): #loop over events
    if ientry > nevtcut: break
    if ientry % (nevtcut/10)==0 : print 'processing ', ientry,'th event'
    ch.GetEntry(ientry)

    for i in range(ch.nGenPart):
    
        # electron, pt eta cuts, status 1
        if( abs(ch.GenPart_pdgId[i]) == 11 and ch.GenPart_pt[i] >= 5 and abs(ch.GenPart_eta[i]) <= 2.5 and ch.GenPart_status[i] == 1):
            AllCutEle.Fill(ch.GenPart_pt[i])

            # with stored mom
            if( not ch.GenPart_genPartIdxMother[i] == -1):
                found = False
                for j in range(len(finds)):
                    if(finds[j] == abs(ch.GenPart_pdgId[ch.GenPart_genPartIdxMother[i]])):
                        found = True
                if(not found):
                    finds[ifi] = abs(ch.GenPart_pdgId[ch.GenPart_genPartIdxMother[i]])
                    ifi += 1

                AllCutEleGoodMom.Fill(ch.GenPart_pt[i])
                if(ch.GenPart_pdgId[ch.GenPart_genPartIdxMother[i]] < 100):
                    EleMomDist.Fill(ch.GenPart_pdgId[ch.GenPart_genPartIdxMother[i]] )


print(finds)
hfile.Write()
quit()
'''

# RECO #
print("Calculating efficiencies with and without recofix.")
for ientry in range(n_entries): #loop over events
    if ientry > nevtcut: break
    if ientry % (nevtcut/10)==0 : print 'processing ', ientry,'th event'
    ch.GetEntry(ientry)

    there_is_good_ele = False
    # Good: SUSY jel esemenyekben, ahol van egy true electron 5 GeV pT folott es 2.5 |eta|-n belul
    for i in range(ch.nGenPart):
        # electron, pt eta cuts, status 1, with a stop mother somewhere in mother history
        if( abs(ch.GenPart_pdgId[i]) == 11 and ch.GenPart_pt[i] >= 5 and abs(ch.GenPart_eta[i]) <= 2.5 and ch.GenPart_status[i] == 1 and hasMomRecursive(i, 1000006, ch)):
            there_is_good_ele = True
            break
    if( there_is_good_ele):
        for i in range(ch.nElectron):
            histos['AllElectronEta'].Fill(ch.Electron_eta[i])

            eledist = 100000
            eleIdx = -1
            for j in range(ch.nGenPart):
            # electron, pt eta cuts, status 1, with a stop mother somewhere in mother history
                if abs(ch.GenPart_pdgId[j]) == 11 and ch.GenPart_pt[j] >= 5 and abs(ch.GenPart_eta[j]) <= 2.5 and ch.GenPart_status[j] == 1 and hasMomRecursive(j, 1000006, ch):
                    dist0 = dR(ch.GenPart_eta[j], ch.Electron_eta[i], ch.GenPart_phi[j], ch.Electron_phi[i])
                    if dRcut(dist0) and dist0 < eledist:
                        eledist = dist0
                        eleIdx = j
            
            if(eleIdx != -1):
                histos['MatchedElectronEta'].Fill(ch.Electron_eta[i])
                histos['MatchedElectronGenEta'].Fill(ch.GenPart_eta[eleIdx])


        for i in range(ch.nPhoton):
            histos['AllPhotonEta'].Fill(ch.Photon_eta[i])

            phodist = 100000
            phoIdx = -1
            for j in range(ch.nGenPart):
            # electron, pt eta cuts, status 1, with a stop mother somewhere in mother history
                if abs(ch.GenPart_pdgId[j]) == 11 and ch.GenPart_pt[j] >= 5 and abs(ch.GenPart_eta[j]) <= 2.5 and ch.GenPart_status[j] == 1 and hasMomRecursive(j, 1000006, ch):
                    dist0 = dR(ch.GenPart_eta[j], ch.Photon_eta[i], ch.GenPart_phi[j], ch.Photon_phi[i])
                    if dRcut(dist0) and dist0 < phodist:
                        phodist = dist0
                        phoIdx = j
            
            if(phoIdx != -1):
                histos['MatchedPhotonEta'].Fill(ch.Photon_eta[i])
                histos['MatchedPhotonGenEta'].Fill(ch.GenPart_eta[phoIdx])

        for i in range(ch.nIsoTrack):
            histos['AllIsoTrackEta'].Fill(ch.IsoTrack_eta[i])

            isodist = 100000
            isoIdx = -1
            for j in range(ch.nGenPart):
            # electron, pt eta cuts, status 1, with a stop mother somewhere in mother history
                if abs(ch.GenPart_pdgId[j]) == 11 and ch.GenPart_pt[j] >= 5 and abs(ch.GenPart_eta[j]) <= 2.5 and ch.GenPart_status[j] == 1 and hasMomRecursive(j, 1000006, ch):
                    dist0 = dR(ch.GenPart_eta[j], ch.IsoTrack_eta[i], ch.GenPart_phi[j], ch.IsoTrack_phi[i])
                    if( dRcut(dist0) and dist0 < isodist):
                        isodist = dist0
                        isoIdx = j
            
            if(isoIdx != -1):
                histos['MatchedIsoTrackEta'].Fill(ch.IsoTrack_eta[i])
                histos['MatchedIsoTrackGenEta'].Fill(ch.GenPart_eta[isoIdx])

            # ugyanezen esemenyekben minden isoTrack eloszlasa, amelyhez kozel van egy true elektron es egy reco photon is (1db abra)
            iso2dist = 100000
            iso2Idx = -1

            for j in range(ch.nGenPart):
                for k in range(ch.nPhoton):
                    if abs(ch.GenPart_pdgId[j]) == 11 and ch.GenPart_pt[j] >= 5 and abs(ch.GenPart_eta[j]) <= 2.5 and ch.GenPart_status[j] == 1 and hasMomRecursive(j, 1000006, ch):
                        dist0e = dR(ch.GenPart_eta[j], ch.IsoTrack_eta[i], ch.GenPart_phi[j], ch.IsoTrack_phi[i])
                        dist0p = dR(ch.Photon_eta[k], ch.IsoTrack_eta[i], ch.Photon_phi[k], ch.IsoTrack_phi[i])
                        if( dRcut(dist0e) and dRcut(dist0p) and dist0e < iso2dist):
                            iso2dist = dist0e
                            iso2Idx = j
            
            if(iso2Idx != -1):
                histos['DoubleMatchedIsoTrackEta'].Fill(ch.IsoTrack_eta[i])
                histos['DoubleMatchedIsoTrackGenEta'].Fill(ch.GenPart_eta[iso2Idx])


# Save histos:
hfile.Write()

'''

    for i in range(ch.nGenPart):
        if( abs(ch.GenPart_pdgId[i]) == 11):
            histos['AllElectronEta'].Fill(ch.GenPart_eta[i])

        # electron, pt eta cuts, status 1, with a stop mother somewhere in mother history
        if( abs(ch.GenPart_pdgId[i]) == 11 and ch.GenPart_pt[i] >= 5 and abs(ch.GenPart_eta[i]) <= 2.5 and ch.GenPart_status[i] == 1 and hasMomRecursive(i, 1000006, ch)):
            histos['TrueElePt'].Fill(ch.GenPart_pt[i])
            histos['TrueEleEta'].Fill(ch.GenPart_eta[i])

            # check if there is a reco ele / isotrack:
            eledist = 100000
            eleIdx = -1
            idx0 = -1
            for j in range(ch.nElectron):
                dist0 = dR(ch.GenPart_eta[i], ch.Electron_eta[j], ch.GenPart_phi[i], ch.Electron_phi[j] ) 
                if( dist0 < eledist):
                    eledist = dist0
                    idx0 = j
            if( dRcut(eledist) ):
                histos['RecoElePt'].Fill(ch.GenPart_pt[i])
                histos['RecoEleEta'].Fill(ch.GenPart_eta[i])

                histos['RecoFixedElePt'].Fill(ch.GenPart_pt[i])
                histos['RecoFixedEleEta'].Fill(ch.GenPart_eta[i])
                eleIdx = idx0

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

            # attempt reco fix
            if( eleIdx == -1 ):
                # reco ele not found - check if there is a reco-ed photon
                there_is_reco_photon = False
                for j in range(ch.nPhoton):
                    if( dR(ch.GenPart_eta[i], ch.Photon_eta[j], ch.GenPart_phi[i], ch.Photon_phi[j] ) < 0.1):
                        there_is_reco_photon = True
                        break
                
                if( there_is_reco_photon):
                    # Recover efficiency by looking for photon, but this has big background
                    # Check if there is an isotrack for true ele -> then the isotrack is the detected ele
                    if( not isoIdx == -1):
                        histos['RecoFixedElePt'].Fill(ch.GenPart_pt[i])
                        histos['RecoFixedEleEta'].Fill(ch.GenPart_eta[i])

                        histos['MatchedIsoTrackEta'].Fill(ch.IsoTrack_eta[isoIdx])
        



        # photon
        if( abs(ch.GenPart_pdgId[i]) == 22 ): # true photon
            histos['TruePhotonPt'].Fill(ch.GenPart_pt[i])
            histos['TruePhotonEta'].Fill(ch.GenPart_eta[i])
            histos['AllPhotonEta'].Fill(ch.GenPart_eta[i])

            if(hasMom(i, 11, ch)):
                histos['MatchedPhotonEta'].Fill(ch.GenPart_eta[i])
            
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


'''
