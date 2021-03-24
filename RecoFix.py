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
if 'FullSim0' != samples and 'FullSim5' != samples and 'FullSim10' != samples and 'FullSim100' != samples and 'FullSim100m' != samples and 'FullSim100m2' != samples:
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



def getdxy(i, ch):
    #genLV = ROOT.TLorentzVector()
    #genLV.SetPtEtaPhiM(ch.GenPart_pt[i],ch.GenPart_eta[i], ch.GenPart_phi[i], ch.GenPart_mass[i])
    px = ch.GenPart_pt[i] * cos(ch.GenPart_phi[i])
    py = ch.GenPart_pt[i] * sin(ch.GenPart_phi[i])
    return ((ch.GenPart_vy[i] - ch.GenVtx_y)  * px - (ch.GenPart_vx[i] - ch.GenVtx_x) * py) / ch.GenPart_pt[i];

def getdz(i, ch):
    px = ch.GenPart_pt[i] * cos(ch.GenPart_phi[i])
    py = ch.GenPart_pt[i] * sin(ch.GenPart_phi[i])
    pz = ch.GenPart_pt[i] / math.tan(ch.GenPart_phi[i])
    return (ch.GenPart_vz[i] - ch.GenVtx_z) - ((ch.GenPart_vx[i] - ch.GenVtx_x) * px + (ch.GenPart_vy[i] - ch.GenVtx_y)*  py) / ch.GenPart_pt[i] * (pz / ch.GenPart_pt[i]);


# HISTOS #
# Save histos:
hfile = ROOT.TFile( 'root_files/RecoFix_Sample'+sample+'.root', 'RECREATE')
ptBinning = [0, 3.5, 5, 12, 20, 30, 50, 80, 200]
histext = ''
histos = {}
histos['TrueElePt'] = HistInfo(hname = 'TrueElePt', sample = histext, binning = [50, 0, 100], histclass = ROOT.TH1F).make_hist()
histos['TrueEleEta'] = HistInfo(hname = 'TrueEleEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['TrueEleVtxx'] = HistInfo(hname = 'TrueEleVtxx', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['TrueEleVtxy'] = HistInfo(hname = 'TrueEleVtxy', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['TrueEleVtxz'] = HistInfo(hname = 'TrueEleVtxz', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['TrueEleDxy'] = HistInfo(hname = 'TrueEleDxy', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
histos['TrueEleDz'] = HistInfo(hname = 'TrueEleDz', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
histos['TrueEleDxyAbs'] = HistInfo(hname = 'TrueEleDxyAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()
histos['TrueEleDzAbs'] = HistInfo(hname = 'TrueEleDzAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()


histos['RecoElePt'] = HistInfo(hname = 'RecoElePt', sample = histext, binning = [50, 0, 100], histclass = ROOT.TH1F).make_hist()
histos['RecoEleEta'] = HistInfo(hname = 'RecoEleEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['RecoEleVtxx'] = HistInfo(hname = 'RecoEleVtxx', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['RecoEleVtxy'] = HistInfo(hname = 'RecoEleVtxy', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['RecoEleVtxz'] = HistInfo(hname = 'RecoEleVtxz', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['RecoEleDxy'] = HistInfo(hname = 'RecoEleDxy', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
histos['RecoEleDz'] = HistInfo(hname = 'RecoEleDz', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
histos['RecoEleDxyAbs'] = HistInfo(hname = 'RecoEleDxyAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()
histos['RecoEleDzAbs'] = HistInfo(hname = 'RecoEleDzAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()


histos['IsoTrackPt'] = HistInfo(hname = 'IsoTrackPt', sample = histext, binning = [50, 0, 100], histclass = ROOT.TH1F).make_hist()
histos['IsoTrackEta'] = HistInfo(hname = 'IsoTrackEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['IsoTrackVtxx'] = HistInfo(hname = 'IsoTrackVtxx', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['IsoTrackVtxy'] = HistInfo(hname = 'IsoTrackVtxy', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['IsoTrackVtxz'] = HistInfo(hname = 'IsoTrackVtxz', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['IsoTrackDxy'] = HistInfo(hname = 'IsoTrackDxy', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
histos['IsoTrackDz'] = HistInfo(hname = 'IsoTrackDz', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
histos['IsoTrackDxyAbs'] = HistInfo(hname = 'IsoTrackDxyAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()
histos['IsoTrackDzAbs'] = HistInfo(hname = 'IsoTrackDzAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()


histos['RecoFixedElePt'] = HistInfo(hname = 'RecoFixedElePt', sample = histext, binning = [50, 0, 100], histclass = ROOT.TH1F).make_hist()
histos['RecoFixedEleEta'] = HistInfo(hname = 'RecoFixedEleEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['RecoFixedEleVtxx'] = HistInfo(hname = 'RecoFixedEleVtxx', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['RecoFixedEleVtxy'] = HistInfo(hname = 'RecoFixedEleVtxy', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['RecoFixedEleVtxz'] = HistInfo(hname = 'RecoFixedEleVtxz', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['RecoFixedEleDxy'] = HistInfo(hname = 'RecoFixedEleDxy', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
histos['RecoFixedEleDz'] = HistInfo(hname = 'RecoFixedEleDz', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
histos['RecoFixedEleDxyAbs'] = HistInfo(hname = 'RecoFixedEleDxyAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()
histos['RecoFixedEleDzAbs'] = HistInfo(hname = 'RecoFixedEleDzAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()

histos['TruePhotonPt'] = HistInfo(hname = 'TruePhotonPt', sample = histext, binning = [50, 0, 100], histclass = ROOT.TH1F).make_hist()
histos['TruePhotonEta'] = HistInfo(hname = 'TruePhotonEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['TruePhotonVtxx'] = HistInfo(hname = 'TruePhotonVtxx', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['TruePhotonVtxy'] = HistInfo(hname = 'TruePhotonVtxy', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['TruePhotonVtxz'] = HistInfo(hname = 'TruePhotonVtxz', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['TruePhotonDxy'] = HistInfo(hname = 'TruePhotonDxy', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
histos['TruePhotonDz'] = HistInfo(hname = 'TruePhotonDz', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
histos['TruePhotonDxyAbs'] = HistInfo(hname = 'TruePhotonDxyAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()
histos['TruePhotonDzAbs'] = HistInfo(hname = 'TruePhotonDzAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()

histos['RecoPhotonPt'] = HistInfo(hname = 'RecoPhotonPt', sample = histext, binning = [50, 0, 100], histclass = ROOT.TH1F).make_hist()
histos['RecoPhotonEta'] = HistInfo(hname = 'RecoPhotonEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['RecoPhotonVtxx'] = HistInfo(hname = 'RecoPhotonVtxx', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['RecoPhotonVtxy'] = HistInfo(hname = 'RecoPhotonVtxy', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['RecoPhotonVtxz'] = HistInfo(hname = 'RecoPhotonVtxz', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
histos['RecoPhotonDxy'] = HistInfo(hname = 'RecoPhotonDxy', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
histos['RecoPhotonDz'] = HistInfo(hname = 'RecoPhotonDz', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
histos['RecoPhotonDxyAbs'] = HistInfo(hname = 'RecoPhotonDxyAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()
histos['RecoPhotonDzAbs'] = HistInfo(hname = 'RecoPhotonDzAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()

effs = {}


n_entries = ch.GetEntries() #num entries in tree
nevtcut = n_entries -1 if nEvents == - 1 else nEvents - 1

# ele has larger impact parameter
# reco eff is low
# try improving

# RECO #
print("Calculating efficiencies with and without recofix.")
for ientry in range(n_entries): #loop over events
    if ientry > nevtcut: break
    if ientry % (nevtcut/10)==0 : print 'processing ', ientry,'th event'
    ch.GetEntry(ientry)


    for i in range(ch.nGenPart):

        # electron, pt eta cuts, status 1, with a stop mother somewhere in mother history
        if( abs(ch.GenPart_pdgId[i]) == 11 and ch.GenPart_pt[i] >= 5 and abs(ch.GenPart_eta[i]) <= 2.5 and ch.GenPart_status[i] == 1 and hasMomRecursive(i, 1000006, ch)):

            GenPart_dxy = getdxy(i, ch)
            GenPart_dz = getdz(i, ch)

            histos['TrueElePt'].Fill(ch.GenPart_pt[i])
            histos['TrueEleEta'].Fill(ch.GenPart_eta[i])
            histos['TrueEleVtxx'].Fill(ch.GenPart_vx[i])
            histos['TrueEleVtxy'].Fill(ch.GenPart_vy[i])
            histos['TrueEleVtxz'].Fill(ch.GenPart_vz[i])
            histos['TrueEleDxy'].Fill( GenPart_dxy)
            histos['TrueEleDz'].Fill( GenPart_dz)
            histos['TrueEleDxyAbs'].Fill( abs(GenPart_dxy))
            histos['TrueEleDzAbs'].Fill( abs(GenPart_dz))

            # check if there is a reco ele / isotrack:
            eledist = 100000
            eleIdx = -1
            idx0 = -1
            for j in range(ch.nElectron):

                # remove isolation
                #eleVID(ch.Electron_vidNestedWPBitmap[j], 0, removedCuts=['pfRelIso03_all'])

                dist0 = dR(ch.GenPart_eta[i], ch.Electron_eta[j], ch.GenPart_phi[i], ch.Electron_phi[j] ) 
                if( dist0 < eledist):
                    eledist = dist0
                    idx0 = j
            if( dRcut(eledist) ):
                histos['RecoElePt'].Fill(ch.GenPart_pt[i])
                histos['RecoEleEta'].Fill(ch.GenPart_eta[i])
                histos['RecoEleVtxx'].Fill(ch.GenPart_vx[i])
                histos['RecoEleVtxy'].Fill(ch.GenPart_vy[i])
                histos['RecoEleVtxz'].Fill(ch.GenPart_vz[i])
                histos['RecoEleDxy'].Fill( GenPart_dxy)
                histos['RecoEleDz'].Fill( GenPart_dz)
                histos['RecoEleDxyAbs'].Fill( abs(GenPart_dxy))
                histos['RecoEleDzAbs'].Fill( abs(GenPart_dz))

                histos['RecoFixedElePt'].Fill(ch.GenPart_pt[i])
                histos['RecoFixedEleEta'].Fill(ch.GenPart_eta[i])
                histos['RecoFixedEleVtxx'].Fill(ch.GenPart_vx[i])
                histos['RecoFixedEleVtxy'].Fill(ch.GenPart_vy[i])
                histos['RecoFixedEleVtxz'].Fill(ch.GenPart_vz[i])
                histos['RecoFixedEleDxy'].Fill(GenPart_dxy)
                histos['RecoFixedEleDz'].Fill( GenPart_dz)
                histos['RecoFixedEleDxyAbs'].Fill(abs(GenPart_dxy))
                histos['RecoFixedEleDzAbs'].Fill( abs(GenPart_dz))

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
                histos['IsoTrackVtxx'].Fill(ch.GenPart_vx[i])
                histos['IsoTrackVtxy'].Fill(ch.GenPart_vy[i])
                histos['IsoTrackVtxz'].Fill(ch.GenPart_vz[i])
                histos['IsoTrackDxy'].Fill( GenPart_dxy)
                histos['IsoTrackDz'].Fill( GenPart_dz)
                histos['IsoTrackDxyAbs'].Fill( abs(GenPart_dxy))
                histos['IsoTrackDzAbs'].Fill( abs(GenPart_dz))
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
                        histos['RecoFixedEleVtxx'].Fill(ch.GenPart_vx[i])
                        histos['RecoFixedEleVtxy'].Fill(ch.GenPart_vy[i])
                        histos['RecoFixedEleVtxz'].Fill(ch.GenPart_vz[i])
                        histos['RecoFixedEleDxy'].Fill( GenPart_dxy)
                        histos['RecoFixedEleDz'].Fill( GenPart_dz)
                        histos['RecoFixedEleDxyAbs'].Fill(abs(GenPart_dxy))
                        histos['RecoFixedEleDzAbs'].Fill( abs(GenPart_dz))
        



        # photon
        if( abs(ch.GenPart_pdgId[i]) == 22 ): # true photon
            GenPart_dxy = getdxy(i, ch)
            GenPart_dz = getdz(i, ch)

            histos['TruePhotonPt'].Fill(ch.GenPart_pt[i])
            histos['TruePhotonEta'].Fill(ch.GenPart_eta[i])
            histos['TruePhotonVtxx'].Fill(ch.GenPart_vx[i])
            histos['TruePhotonVtxy'].Fill(ch.GenPart_vy[i])
            histos['TruePhotonVtxz'].Fill(ch.GenPart_vz[i])
            histos['TruePhotonDxy'].Fill( GenPart_dxy)
            histos['TruePhotonDz'].Fill( GenPart_dz)
            histos['TruePhotonDxyAbs'].Fill( abs(GenPart_dxy))
            histos['TruePhotonDzAbs'].Fill( abs(GenPart_dz))
            
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
                histos['RecoPhotonVtxx'].Fill(ch.GenPart_vx[i])
                histos['RecoPhotonVtxy'].Fill(ch.GenPart_vy[i])
                histos['RecoPhotonVtxz'].Fill(ch.GenPart_vz[i])
                histos['RecoPhotonDxy'].Fill( GenPart_dxy)
                histos['RecoPhotonDz'].Fill( GenPart_dz)
                histos['RecoPhotonDxyAbs'].Fill( abs(GenPart_dxy))
                histos['RecoPhotonDzAbs'].Fill( abs(GenPart_dz))



# Save histos:
hfile.Write()

'''
Legacy code
# RECO FIXING #

# for true ele, find ele -> efficiency for ele reco:
# TODO: 100 or 10 mm instead
print("Attempting reco fix.")
print 'Running over total events: ', nevtcut+1
for ientry in range(n_entries): #loop over events
    if ientry > nevtcut: break
    if ientry % (nevtcut/10)==0 : print 'processing ', ientry,'th event'
    ch.GetEntry(ientry)

    for i in range(ch.nGenPart):
        if( abs(ch.GenPart_pdgId[i]) == 11 ): # true ele
            # check if there is a reco ele:
            there_is_reco_ele = False
            for j in range(ch.nElectron):
                if( dR(ch.GenPart_eta[i], ch.Electron_eta[j], ch.GenPart_phi[i], ch.Electron_phi[j] ) < 0.1):
                    there_is_reco_ele = True
                    break
        
            if( there_is_reco_ele):
                continue
            
            # not found - check if there is a reco-ed photon

            there_is_reco_photon = False
            for j in range(ch.nPhoton):
                if( dR(ch.GenPart_eta[i], ch.Photon_eta[j], ch.GenPart_phi[i], ch.Photon_phi[j] ) < 0.1):
                    there_is_reco_photon = True
                    break
            
            if( not there_is_reco_photon):
                continue

            # Recover efficiency by looking for photon, but this has big background
            # Check if there is an isotrack for true ele -> then the isotrack is the detected ele
            dist = 100000
            IsoTrackIdx = -1
            for j in range(ch.nIsoTrack):
                dist0 = dR(ch.GenPart_eta[i], ch.IsoTrack_eta[j], ch.GenPart_phi[i], ch.IsoTrack_phi[j])
                if( dist0 < dist):
                    dist = dist0
                    IsoTrackIdx = j
            
            if(dist < 0.1):
                print("Found one!")
'''

