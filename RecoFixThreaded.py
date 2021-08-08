import os, sys
import ROOT
import types
import math
import multiprocessing
import threading

#from RegSel import RegSel
import numpy as np

sys.path.append('../')
from Helper.HistInfo import HistInfo
from Helper.MCWeight import MCWeight
from TriggerStudy.TrigVarSel import TrigVarSel
from Sample.SampleChain import SampleChain
from Helper.VarCalc import *
from Helper.PlotHelper import *
from Sample.FileList_2016 import samples as samples_2016

from isoTrackPhotonPair import *

def get_parser():
    ''' Argument parser.
    '''
    import argparse
    argParser = argparse.ArgumentParser(description = "Argument parser")
    argParser.add_argument('--sample',           action='store',                     type=str,            default='UL17_Full99mm',                                help="Which sample?" )
    argParser.add_argument('--year',             action='store',                     type=int,            default=2016,                                             help="Which year?" )
    argParser.add_argument('--startfile',        action='store',                     type=int,            default=0,                                                help="start from which root file like 0th or 10th etc?" )
    argParser.add_argument('--nfiles',           action='store',                     type=int,            default=-1,                                               help="No of files to run. -1 means all files" )
    argParser.add_argument('--nevents',           action='store',                    type=int,            default=-1,                                               help="No of events to run. -1 means all events" )
    argParser.add_argument('--Nthreads',           action='store',                    type=int,            default=-1,                                               help="How many CPU?" )
    argParser.add_argument('--channel',           action='store',                    type=str,            default='Electron',                                       help="Which lepton?" )
    argParser.add_argument('--region',            action='store',                    type=str,            default='mesurement',                                     help="Which lepton?" )    
    argParser.add_argument('--pJobs',             action='store',                    type=bool,            default=False,                                           help="using GPU parallel program or not" )
    argParser.add_argument('--vertex',             action='store',                    type=str,            default="True",                                           help="is there vertex info" )
    argParser.add_argument('--extrapolate',             action='store',                    type=str,            default="True",                                           help="extrapolate isotrack to ECAL surface" )
    return argParser

options = get_parser().parse_args()

Nthreads = options.Nthreads
cpu_available = 8
#os.cpu_count()

if(Nthreads == -1):
    Nthreads = cpu_available
elif( Nthreads < 1):
    print("Error: invalid thread number")
    exit()



samples  = options.sample
channel = options.channel
nEvents = options.nevents
year = options.year
if( options.vertex == "True"):
    bUseVertex = True
else:
    bUseVertex = False

if( options.extrapolate == "True"):
    bExtrapolate = True
else:
    bExtrapolate = False



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


sample = samples


def dR(eta_gen, eta_reco, phi_gen, phi_reco):
    phipart = phi_gen-phi_reco
    if( phipart < -math.pi):
        phipart += 2*math.pi
    elif( phipart > math.pi):
        phipart -= 2*math.pi
    return math.sqrt((eta_gen-eta_reco)**2 + (phipart)**2 )

def dRcut(dist, cut = 0.1):
    return dist < cut

def hasMomRecursive(i, pdgid, ch):
    if( ch.GenPart_genPartIdxMother[i] == -1):
        return False
    elif(abs(ch.GenPart_pdgId[ch.GenPart_genPartIdxMother[i]]) == pdgid):
        return True
    else:
        return hasMomRecursive( ch.GenPart_genPartIdxMother[i], pdgid, ch )


def extrapolateTrack(pt, eta, phi, charge, x0, y0, z0):
    MAX_R=1.29
    MAX_Z=2.935
            
    R = pt / (0.3 * 3.8 )
    
    xC = x0/100.0 - R* math.cos(phi - charge * math.pi/2.0)
    yC = y0/100.0 - R* math.sin(phi - charge * math.pi/2.0)

    # calculate x,y intersection of track
    a = (-1* yC) / xC
    rb = MAX_R
    RC2 = xC**2 + yC**2
    b = (RC2 - R**2 + rb**2) / (2*xC)

    qa = a**2 + 1
    qb = 2*a*b
    qc = b**2 - rb**2
    disc = qb**2 - 4*qa*qc
    #print "disc:", disc
    
    y,x,y_other,x_other = 0,0,0,0
    if( disc > 0):
        #print "+"
        # barrel can be hit, solution exists
        y1 = (-qb + math.sqrt(disc)) / (2*qa)
        y2 = (-qb - math.sqrt(disc)) / (2*qa)
        x1 = b + y1*a
        x2 = b + y2*a
        
        if(phi > 0):
            y = y1
            x = x1
            y_other = y2
            x_other = x2
        else:
            y = y2
            x = x2
            y_other = y1
            x_other = x1

        #Z_ECAL = z0/100.0 + R * math.sinh(eta) *(math.acos(1 - 1.29*1.29 / (2*R*R)));
        yangle=math.asin((y-yC)/R)
        xangle=math.acos((x-xC)/R)
        angle=math.atan2(y-yC,x-xC)
        phi0=phi - charge*math.pi/2.0
        AngleY=yangle - phi0
        AngleX=xangle - phi0
        Angle=angle - phi0
        Z_ECAL_X = z0/100.0 + AngleX * (R*math.sinh(eta) / charge)
        Z_ECAL_Y = z0/100.0 + AngleY * (R*math.sinh(eta) / charge)
        Z_ECAL = z0/100.0 + Angle * (R*math.sinh(eta) / charge)
                
        if(Z_ECAL > MAX_Z):
            Z_ECAL = MAX_Z
        if(Z_ECAL < -MAX_Z):
            Z_ECAL = -MAX_Z
    else:
        # Barrel cannot be hit, endcap is hit (electron spirals out)
        if(eta > 0):
            Z_ECAL = MAX_Z
        else:
            Z_ECAL = -MAX_Z
# -charge valtoztatva +ra
    X_ECAL = xC + R * math.cos(charge * (Z_ECAL-z0/100.0)/(R * math.sinh(eta)) + phi - charge * math.pi/2.0)
    Y_ECAL = yC + R * math.sin(charge * (Z_ECAL-z0/100.0)/(R * math.sinh(eta)) + phi - charge * math.pi/2.0)
  
    D_ECAL = math.sqrt(X_ECAL*X_ECAL+Y_ECAL*Y_ECAL)

    etaSC = math.asinh(Z_ECAL/D_ECAL)
    if(Y_ECAL > 0):
        phiSC = math.acos(X_ECAL/D_ECAL)
    else:
        phiSC = -1.0 * math.acos(X_ECAL/D_ECAL)
        

    return etaSC, phiSC, x, y, x_other, y_other, xC, yC






def phiToXY(phi):
    rb = 129.0 / 100.0
    x = rb * math.cos(phi)
    y = rb * math.sin(phi)
    return x,y


def GetGenPartSC(i, ch):
    return extrapolateTrack(ch.GenPart_pt[i], ch.GenPart_eta[i], ch.GenPart_phi[i], -1*math.copysign(1,ch.GenPart_pdgId[i]), ch.GenPart_vx[i], ch.GenPart_vy[i], ch.GenPart_vz[i])

# here, but not used for efficiency, no need to always recalculate the SC
# maybe use decorator in the future?
def dRFromSC(eta_reco, phi_reco, i, ch):
    etaSC, phiSC = GetGenPartSC(i, ch)
    return dR(etaSC, eta_reco, phiSC, phi_reco)




def getdxy(i, ch):
    #genLV = ROOT.TLorentzVector()
    #genLV.SetPtEtaPhiM(ch.GenPart_pt[i],ch.GenPart_eta[i], ch.GenPart_phi[i], ch.GenPart_mass[i])
    px = ch.GenPart_pt[i] * cos(ch.GenPart_phi[i])
    py = ch.GenPart_pt[i] * sin(ch.GenPart_phi[i])
    return ((ch.GenPart_vy[i] - ch.GenVtx_y)  * px - (ch.GenPart_vx[i] - ch.GenVtx_x) * py) / ch.GenPart_pt[i]

def getdz(i, ch):
    px = ch.GenPart_pt[i] * cos(ch.GenPart_phi[i])
    py = ch.GenPart_pt[i] * sin(ch.GenPart_phi[i])
    pz = ch.GenPart_pt[i] * math.sinh(ch.GenPart_eta[i])
    return (ch.GenPart_vz[i] - ch.GenVtx_z) - ((ch.GenPart_vx[i] - ch.GenVtx_x) * px + (ch.GenPart_vy[i] - ch.GenVtx_y)*  py) / ch.GenPart_pt[i] * (pz / ch.GenPart_pt[i])



# RECO #
def RecoAndRecoFix(threadID, it0, it1, nevtcut, ch_common):

    # Output histos - standalone for each Process
    ptBinning = [0, 3.5, 5, 7.5, 10, 12, 14, 16, 18, 20, 30, 40, 50, 100]
    histext = ''
    hfile = ROOT.TFile( 'root_files/RecoFixThreaded_Sample'+sample+'_Thread-'+str(threadID)+'.root', 'RECREATE')
    histos = {}
    histos['TrueElePt'] = HistInfo(hname = 'TrueElePt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()
    histos['TrueEleEta'] = HistInfo(hname = 'TrueEleEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['TrueElePhi'] = HistInfo(hname = 'TrueElePhi', sample = histext, binning = [30, -4, 4], histclass = ROOT.TH1F).make_hist()
    histos['TrueEleVtxx'] = HistInfo(hname = 'TrueEleVtxx', sample = histext, binning = [40, -1, 1], histclass = ROOT.TH1F).make_hist()
    histos['TrueEleVtxy'] = HistInfo(hname = 'TrueEleVtxy', sample = histext, binning = [40, -1, 1], histclass = ROOT.TH1F).make_hist()
    histos['TrueEleVtxz'] = HistInfo(hname = 'TrueEleVtxz', sample = histext, binning = [50, -50, 50], histclass = ROOT.TH1F).make_hist()
    histos['TrueEleDxy'] = HistInfo(hname = 'TrueEleDxy', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
    histos['TrueEleDz'] = HistInfo(hname = 'TrueEleDz', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
    histos['TrueEleDxyAbs'] = HistInfo(hname = 'TrueEleDxyAbs', sample = histext, binning = [40, 0, 15], histclass = ROOT.TH1F).make_hist()
    histos['TrueEleDzAbs'] = HistInfo(hname = 'TrueEleDzAbs', sample = histext, binning = [40, 0, 15], histclass = ROOT.TH1F).make_hist()

    histos['RecoElePt'] = HistInfo(hname = 'RecoElePt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()
    histos['RecoEleEta'] = HistInfo(hname = 'RecoEleEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['RecoElePhi'] = HistInfo(hname = 'RecoElePhi', sample = histext, binning = [30, -4, 4], histclass = ROOT.TH1F).make_hist()
    histos['RecoEleVtxx'] = HistInfo(hname = 'RecoEleVtxx', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['RecoEleVtxy'] = HistInfo(hname = 'RecoEleVtxy', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['RecoEleVtxz'] = HistInfo(hname = 'RecoEleVtxz', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['RecoEleDxy'] = HistInfo(hname = 'RecoEleDxy', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
    histos['RecoEleDz'] = HistInfo(hname = 'RecoEleDz', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
    histos['RecoEleDxyAbs'] = HistInfo(hname = 'RecoEleDxyAbs', sample = histext, binning = [40, 0, 15], histclass = ROOT.TH1F).make_hist()
    histos['RecoEleDzAbs'] = HistInfo(hname = 'RecoEleDzAbs', sample = histext, binning = [40, 0, 15], histclass = ROOT.TH1F).make_hist()

    histos['IsoTrackPt'] = HistInfo(hname = 'IsoTrackPt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()
    histos['IsoTrackEta'] = HistInfo(hname = 'IsoTrackEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['IsoTrackPhi'] = HistInfo(hname = 'IsoTrackPhi', sample = histext, binning = [30, -4, 4], histclass = ROOT.TH1F).make_hist()
    histos['IsoTrackVtxx'] = HistInfo(hname = 'IsoTrackVtxx', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['IsoTrackVtxy'] = HistInfo(hname = 'IsoTrackVtxy', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['IsoTrackVtxz'] = HistInfo(hname = 'IsoTrackVtxz', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['IsoTrackDxy'] = HistInfo(hname = 'IsoTrackDxy', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
    histos['IsoTrackDz'] = HistInfo(hname = 'IsoTrackDz', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
    histos['IsoTrackDxyAbs'] = HistInfo(hname = 'IsoTrackDxyAbs', sample = histext, binning = [40, 0, 15], histclass = ROOT.TH1F).make_hist()
    histos['IsoTrackDzAbs'] = HistInfo(hname = 'IsoTrackDzAbs', sample = histext, binning = [40, 0, 15], histclass = ROOT.TH1F).make_hist()

    histos['RecoFixedElePt'] = HistInfo(hname = 'RecoFixedElePt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()
    histos['RecoFixedEleEta'] = HistInfo(hname = 'RecoFixedEleEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['RecoFixedElePhi'] = HistInfo(hname = 'RecoFixedElePhi', sample = histext, binning = [30, -4, 4], histclass = ROOT.TH1F).make_hist()
    histos['RecoFixedEleVtxx'] = HistInfo(hname = 'RecoFixedEleVtxx', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['RecoFixedEleVtxy'] = HistInfo(hname = 'RecoFixedEleVtxy', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['RecoFixedEleVtxz'] = HistInfo(hname = 'RecoFixedEleVtxz', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['RecoFixedEleDxy'] = HistInfo(hname = 'RecoFixedEleDxy', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
    histos['RecoFixedEleDz'] = HistInfo(hname = 'RecoFixedEleDz', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
    histos['RecoFixedEleDxyAbs'] = HistInfo(hname = 'RecoFixedEleDxyAbs', sample = histext, binning = [40, 0, 15], histclass = ROOT.TH1F).make_hist()
    histos['RecoFixedEleDzAbs'] = HistInfo(hname = 'RecoFixedEleDzAbs', sample = histext, binning = [40, 0, 15], histclass = ROOT.TH1F).make_hist()

    histos['RecoFixedExtrapolElePt'] = HistInfo(hname = 'RecoFixedExtrapolElePt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()
    histos['RecoFixedExtrapolEleEta'] = HistInfo(hname = 'RecoFixedExtrapolEleEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['RecoFixedExtrapolElePhi'] = HistInfo(hname = 'RecoFixedExtrapolElePhi', sample = histext, binning = [30, -4, 4], histclass = ROOT.TH1F).make_hist()
    histos['RecoFixedExtrapolEleVtxx'] = HistInfo(hname = 'RecoFixedExtrapolEleVtxx', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['RecoFixedExtrapolEleVtxy'] = HistInfo(hname = 'RecoFixedExtrapolEleVtxy', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['RecoFixedExtrapolEleVtxz'] = HistInfo(hname = 'RecoFixedExtrapolEleVtxz', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['RecoFixedExtrapolEleDxy'] = HistInfo(hname = 'RecoFixedExtrapolEleDxy', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
    histos['RecoFixedExtrapolEleDz'] = HistInfo(hname = 'RecoFixedExtrapolEleDz', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
    histos['RecoFixedExtrapolEleDxyAbs'] = HistInfo(hname = 'RecoFixedExtrapolEleDxyAbs', sample = histext, binning = [40, 0, 15], histclass = ROOT.TH1F).make_hist()
    histos['RecoFixedExtrapolEleDzAbs'] = HistInfo(hname = 'RecoFixedExtrapolEleDzAbs', sample = histext, binning = [40, 0, 15], histclass = ROOT.TH1F).make_hist()


    histos['RecoFixedExtrapolPtDxy2D'] = ROOT.TH2F('RecoFixedExtrapolPtDxy2D_','RecoFixedExtrapolPtDxy2D', len(ptBinning) -1, np.array(ptBinning), 40,0,15)
    histos['RecoFixedPtDxy2D'] = ROOT.TH2F('RecoFixedPtDxy2D_','RecoFixedPtDxy2D', len(ptBinning) -1, np.array(ptBinning), 40,0,15)


    histos['TruePhotonPt'] = HistInfo(hname = 'TruePhotonPt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()
    histos['TruePhotonEta'] = HistInfo(hname = 'TruePhotonEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['TruePhotonPhi'] = HistInfo(hname = 'TruePhotonPhi', sample = histext, binning = [30, -4, 4], histclass = ROOT.TH1F).make_hist()
    histos['TruePhotonVtxx'] = HistInfo(hname = 'TruePhotonVtxx', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['TruePhotonVtxy'] = HistInfo(hname = 'TruePhotonVtxy', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['TruePhotonVtxz'] = HistInfo(hname = 'TruePhotonVtxz', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['TruePhotonDxy'] = HistInfo(hname = 'TruePhotonDxy', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
    histos['TruePhotonDz'] = HistInfo(hname = 'TruePhotonDz', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
    histos['TruePhotonDxyAbs'] = HistInfo(hname = 'TruePhotonDxyAbs', sample = histext, binning = [40, 0, 15], histclass = ROOT.TH1F).make_hist()
    histos['TruePhotonDzAbs'] = HistInfo(hname = 'TruePhotonDzAbs', sample = histext, binning = [40, 0, 15], histclass = ROOT.TH1F).make_hist()

    histos['RecoPhotonPt'] = HistInfo(hname = 'RecoPhotonPt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()
    histos['RecoPhotonEta'] = HistInfo(hname = 'RecoPhotonEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['RecoPhotonPhi'] = HistInfo(hname = 'RecoPhotonPhi', sample = histext, binning = [30, -4, 4], histclass = ROOT.TH1F).make_hist()
    histos['RecoPhotonVtxx'] = HistInfo(hname = 'RecoPhotonVtxx', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['RecoPhotonVtxy'] = HistInfo(hname = 'RecoPhotonVtxy', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['RecoPhotonVtxz'] = HistInfo(hname = 'RecoPhotonVtxz', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['RecoPhotonDxy'] = HistInfo(hname = 'RecoPhotonDxy', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
    histos['RecoPhotonDz'] = HistInfo(hname = 'RecoPhotonDz', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
    histos['RecoPhotonDxyAbs'] = HistInfo(hname = 'RecoPhotonDxyAbs', sample = histext, binning = [40, 0, 15], histclass = ROOT.TH1F).make_hist()
    histos['RecoPhotonDzAbs'] = HistInfo(hname = 'RecoPhotonDzAbs', sample = histext, binning = [40, 0, 15], histclass = ROOT.TH1F).make_hist()


    histos['RecoPhotonPlusIsoTrackPt'] = HistInfo(hname = 'RecoPhotonPlusIsoTrackPt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()
    histos['RecoPhotonPlusIsoTrackEta'] = HistInfo(hname = 'RecoPhotonPlusIsoTrackEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['RecoPhotonPlusIsoTrackPhi'] = HistInfo(hname = 'RecoPhotonPlusIsoTrackPhi', sample = histext, binning = [30, -4, 4], histclass = ROOT.TH1F).make_hist()
    histos['RecoPhotonPlusIsoTrackDxyAbs'] = HistInfo(hname = 'RecoPhotonPlusIsoTrackDxyAbs', sample = histext, binning = [40, 0, 15], histclass = ROOT.TH1F).make_hist()
    histos['RecoPhotonPlusIsoTrackDzAbs'] = HistInfo(hname = 'RecoPhotonPlusIsoTrackDzAbs', sample = histext, binning = [40, 0, 15], histclass = ROOT.TH1F).make_hist()

    histos['dRControlPlot2D'] = HistInfo(hname = 'dRControlPlot2D', sample = histext, binning = [[40, 0, 2],[40, 0, 2]], histclass = ROOT.TH2F).make_hist2D()

    histos['TrueBackgroundPt'] = HistInfo(hname = 'TrueBackgroundPt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()
    histos['TrueBackgroundEta'] = HistInfo(hname = 'TrueBackgroundEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['TrueBackgroundPhi'] = HistInfo(hname = 'TrueBackgroundPhi', sample = histext, binning = [30, -4, 4], histclass = ROOT.TH1F).make_hist()
    histos['TrueBackgroundDxyAbs'] = HistInfo(hname = 'TrueBackgroundDxyAbs', sample = histext, binning = [40, 0, 15], histclass = ROOT.TH1F).make_hist()
    histos['TrueBackgroundDzAbs'] = HistInfo(hname = 'TrueBackgroundDzAbs', sample = histext, binning = [40, 0, 15], histclass = ROOT.TH1F).make_hist()

    histos['RejectBackgroundPt'] = HistInfo(hname = 'RejectBackgroundPt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()
    histos['RejectBackgroundEta'] = HistInfo(hname = 'RejectBackgroundEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['RejectBackgroundPhi'] = HistInfo(hname = 'RejectBackgroundPhi', sample = histext, binning = [30, -4, 4], histclass = ROOT.TH1F).make_hist()
    histos['RejectBackgroundDxyAbs'] = HistInfo(hname = 'RejectBackgroundDxyAbs', sample = histext, binning = [40, 0, 15], histclass = ROOT.TH1F).make_hist()
    histos['RejectBackgroundDzAbs'] = HistInfo(hname = 'RejectBackgroundDzAbs', sample = histext, binning = [40, 0, 15], histclass = ROOT.TH1F).make_hist()

    histos['RejectBackgroundRecoFixPt'] = HistInfo(hname = 'RejectBackgroundRecoFixPt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()
    histos['RejectBackgroundRecoFixEta'] = HistInfo(hname = 'RejectBackgroundRecoFixEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['RejectBackgroundRecoFixPhi'] = HistInfo(hname = 'RejectBackgroundRecoFixPhi', sample = histext, binning = [30, -4, 4], histclass = ROOT.TH1F).make_hist()
    histos['RejectBackgroundRecoFixDxyAbs'] = HistInfo(hname = 'RejectBackgroundRecoFixDxyAbs', sample = histext, binning = [40, 0, 15], histclass = ROOT.TH1F).make_hist()
    histos['RejectBackgroundRecoFixDzAbs'] = HistInfo(hname = 'RejectBackgroundRecoFixDzAbs', sample = histext, binning = [40, 0, 15], histclass = ROOT.TH1F).make_hist()

    histos['AllTruePhotonPt'] = HistInfo(hname = 'AllTruePhotonPt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()
    histos['AllRecoPhotonPt'] = HistInfo(hname = 'AllRecoPhotonPt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()

    extrapolate_log = open("logs/extrapolate-log-t"+str(threadID),"w")
    # RECO CALCULATION #
    for ientry in range(it0,it1,1): 
        if ientry > nevtcut: break
        if ientry == it0 : print 'Process-'+str(threadID)+' starting at ', ientry,'th event'
        if (ientry % (nevtcut/10)==0 and ientry != it0 and ientry != it1 - 1) : print 'Process-'+str(threadID)+' processing ', ientry,'th event'
        if ientry == it1 - 1: print 'Process-'+str(threadID)+' finishing at ', ientry,'th event'

        ch = ch_common
        ch.GetEntry(ientry)


        for i in range(ch.nPhoton):
            Fill1D(histos['AllRecoPhotonPt'], ch.Photon_pt[i])


        # create isotrack - Photon pairs for extended reco
        aPhotonIsoTrackPairs = []
        aPhotonIsoTrackExtrapolPairs = []
        for i in range(ch.nIsoTrack):
            fTrackPhoDist = 100000
            fTrackPhoDistExtrapol = 100000
            dist0 = 100000
            for j in range(ch.nPhoton):
                
                gencomp_eta = ch.IsoTrack_eta[i]
                gencomp_phi = ch.IsoTrack_phi[i]
                dist0 = dR(gencomp_eta, ch.Photon_eta[j], gencomp_phi, ch.Photon_phi[j] )

                if( dist0 < fTrackPhoDist):
                    fTrackPhoDist = dist0

                if( bExtrapolate):
                    gencomp_eta = ch.IsoTrack_eta[i] + ch.IsoTrack_deltaEta[i]
                    gencomp_phi = ch.IsoTrack_phi[i] + ch.IsoTrack_deltaPhi[i]
                    dist0 = dR(gencomp_eta, ch.Photon_eta[j], gencomp_phi, ch.Photon_phi[j] )

                    if( dist0 < fTrackPhoDistExtrapol):
                        fTrackPhoDistExtrapol = dist0

            
            if( dRcut(fTrackPhoDist, 0.2) ):
                aPhotonIsoTrackPairs.append( isoTrackPhotonPair(ch, i, bExtrapolate) )

            if( dRcut(fTrackPhoDistExtrapol, 0.2) ):
                aPhotonIsoTrackExtrapolPairs.append( isoTrackPhotonPair(ch, i, bExtrapolate) )
        '''
        if(len(aPhotonIsoTrackPairs) < len(aPhotonIsoTrackExtrapolPairs)):
            print 'SUCCESS '+str( len(aPhotonIsoTrackPairs))+ ' '+str( len(aPhotonIsoTrackExtrapolPairs))
        if(len(aPhotonIsoTrackPairs) > len(aPhotonIsoTrackExtrapolPairs)):
            print 'FAIL '+str( len(aPhotonIsoTrackPairs))+ ' '+str( len(aPhotonIsoTrackExtrapolPairs))
        '''


        ## RECO STUDY ##
        for i in range(ch.nGenPart):

            ##  photon pt test
            if(abs(ch.GenPart_pdgId[i]) == 22):
                Fill1D(histos['AllTruePhotonPt'], ch.GenPart_pt[i])


            # electron, pt eta cuts, status 1, with a stop mother somewhere in mother history
            # pt >= 15, |eta| < 1
            # and abs(ch.GenPart_vx[i]) < 0.2 and abs(ch.GenPart_vy[i]) < 0.2
            if( abs(ch.GenPart_pdgId[i]) == 11 and ch.GenPart_pt[i] >= 5  and abs(ch.GenPart_eta[i]) <= 2.5 and ch.GenPart_status[i] == 1 and hasMomRecursive(i, 1000006, ch) ):

                if(bUseVertex):
                    GenPart_dxy = getdxy(i, ch)
                    GenPart_dz = getdz(i, ch)

                Fill1D(histos['TrueElePt'],ch.GenPart_pt[i])
                Fill1D(histos['TrueEleEta'],ch.GenPart_eta[i])
                Fill1D(histos['TrueElePhi'],ch.GenPart_phi[i])

                if(bUseVertex):
                    Fill1D(histos['TrueEleVtxx'],ch.GenPart_vx[i])
                    Fill1D(histos['TrueEleVtxy'],ch.GenPart_vy[i])
                    Fill1D(histos['TrueEleVtxz'],ch.GenPart_vz[i])
                    Fill1D(histos['TrueEleDxy'], GenPart_dxy)
                    Fill1D(histos['TrueEleDz'], GenPart_dz)
                    Fill1D(histos['TrueEleDxyAbs'], abs(GenPart_dxy))
                    Fill1D(histos['TrueEleDzAbs'], abs(GenPart_dz))

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
                    eleIdx = idx0
                    Fill1D(histos['RecoElePt'], ch.GenPart_pt[i])
                    Fill1D(histos['RecoEleEta'], ch.GenPart_eta[i])
                    Fill1D(histos['RecoElePhi'], ch.GenPart_phi[i])

                    if(bUseVertex):
                        Fill1D(histos['RecoEleVtxx'], ch.GenPart_vx[i])
                        Fill1D(histos['RecoEleVtxy'], ch.GenPart_vy[i])
                        Fill1D(histos['RecoEleVtxz'], ch.GenPart_vz[i])
                        Fill1D(histos['RecoEleDxy'], GenPart_dxy)
                        Fill1D(histos['RecoEleDz'], GenPart_dz)
                        Fill1D(histos['RecoEleDxyAbs'], abs(GenPart_dxy))
                        Fill1D(histos['RecoEleDzAbs'], abs(GenPart_dz))

                    Fill1D(histos['RecoFixedElePt'], ch.GenPart_pt[i])
                    Fill1D(histos['RecoFixedEleEta'], ch.GenPart_eta[i])
                    Fill1D(histos['RecoFixedElePhi'], ch.GenPart_phi[i])

                    Fill1D(histos['RecoFixedExtrapolElePt'], ch.GenPart_pt[i])
                    Fill1D(histos['RecoFixedExtrapolEleEta'], ch.GenPart_eta[i])
                    Fill1D(histos['RecoFixedExtrapolElePhi'], ch.GenPart_phi[i])

                    

                    if(bUseVertex):
                        Fill1D(histos['RecoFixedEleVtxx'], ch.GenPart_vx[i])
                        Fill1D(histos['RecoFixedEleVtxy'], ch.GenPart_vy[i])
                        Fill1D(histos['RecoFixedEleVtxz'], ch.GenPart_vz[i])
                        Fill1D(histos['RecoFixedEleDxy'], GenPart_dxy)
                        Fill1D(histos['RecoFixedEleDz'], GenPart_dz)
                        Fill1D(histos['RecoFixedEleDxyAbs'], abs(GenPart_dxy))
                        Fill1D(histos['RecoFixedEleDzAbs'], abs(GenPart_dz))

                        Fill1D(histos['RecoFixedExtrapolEleVtxx'], ch.GenPart_vx[i])
                        Fill1D(histos['RecoFixedExtrapolEleVtxy'], ch.GenPart_vy[i])
                        Fill1D(histos['RecoFixedExtrapolEleVtxz'], ch.GenPart_vz[i])
                        Fill1D(histos['RecoFixedExtrapolEleDxy'], GenPart_dxy)
                        Fill1D(histos['RecoFixedExtrapolEleDz'], GenPart_dz)
                        Fill1D(histos['RecoFixedExtrapolEleDxyAbs'], abs(GenPart_dxy))
                        Fill1D(histos['RecoFixedExtrapolEleDzAbs'], abs(GenPart_dz))

                        histos['RecoFixedExtrapolPtDxy2D'].Fill( ch.GenPart_pt[i], abs(GenPart_dxy) )
                        histos['RecoFixedPtDxy2D'].Fill( ch.GenPart_pt[i], abs(GenPart_dxy) )
                        



                # Extended reco 2.0
                # match not found electrons to isotrack-photon pairs


                if( eleIdx == -1):
                    # reco ele not found - check if there is a reco-ed photon
                    # Recover efficiency by looking for photon, but this has big background
                    # Check if there is an isotrack for true ele -> then the isotrack is the detected ele


                    fPairDist = 100000
                    iPairIdx = -1
                    for j in range(len(aPhotonIsoTrackPairs)):
                        
                        dist0 = dR(ch.GenPart_eta[i], aPhotonIsoTrackPairs[j].eta, ch.GenPart_phi[i], aPhotonIsoTrackPairs[j].phi ) 

                        if( dist0 < fPairDist):
                            fPairDist = dist0
                            iPairIdx = j

                    
                    if( dRcut(fPairDist, 0.2) ):
                        Fill1D(histos['RecoFixedElePt'], ch.GenPart_pt[i])
                        Fill1D(histos['RecoFixedEleEta'], ch.GenPart_eta[i])
                        Fill1D(histos['RecoFixedElePhi'], ch.GenPart_phi[i])

                        if(bUseVertex):
                            Fill1D(histos['RecoFixedEleVtxx'], ch.GenPart_vx[i])
                            Fill1D(histos['RecoFixedEleVtxy'], ch.GenPart_vy[i])
                            Fill1D(histos['RecoFixedEleVtxz'], ch.GenPart_vz[i])
                            Fill1D(histos['RecoFixedEleDxy'], GenPart_dxy)
                            Fill1D(histos['RecoFixedEleDz'], GenPart_dz)
                            Fill1D(histos['RecoFixedEleDxyAbs'], abs(GenPart_dxy))
                            Fill1D(histos['RecoFixedEleDzAbs'], abs(GenPart_dz))

                            histos['RecoFixedPtDxy2D'].Fill( ch.GenPart_pt[i], abs(GenPart_dxy) )
        
                        # these are not filled when ele is found, just here:
                        Fill1D(histos['RecoPhotonPlusIsoTrackPt'], ch.GenPart_pt[i])
                        Fill1D(histos['RecoPhotonPlusIsoTrackEta'], ch.GenPart_eta[i])
                        Fill1D(histos['RecoPhotonPlusIsoTrackPhi'], ch.GenPart_phi[i])
                        if(bUseVertex):
                            Fill1D(histos['RecoPhotonPlusIsoTrackDxyAbs'], abs(GenPart_dxy))
                            Fill1D(histos['RecoPhotonPlusIsoTrackDzAbs'], abs(GenPart_dz))


                    if(bExtrapolate):
                        fPairDist = 100000
                        iPairIdx = -1
                        for j in range(len(aPhotonIsoTrackExtrapolPairs)):
                            
                            dist0 = dR(ch.GenPart_eta[i], aPhotonIsoTrackExtrapolPairs[j].eta, ch.GenPart_phi[i], aPhotonIsoTrackExtrapolPairs[j].phi ) 

                            if( dist0 < fPairDist):
                                fPairDist = dist0
                                iPairIdx = j

                        
                        if( dRcut(fPairDist, 0.2) ):
                            Fill1D(histos['RecoFixedExtrapolElePt'], ch.GenPart_pt[i])
                            Fill1D(histos['RecoFixedExtrapolEleEta'], ch.GenPart_eta[i])
                            Fill1D(histos['RecoFixedExtrapolElePhi'], ch.GenPart_phi[i])

                            if(bUseVertex):
                                Fill1D(histos['RecoFixedExtrapolEleVtxx'], ch.GenPart_vx[i])
                                Fill1D(histos['RecoFixedExtrapolEleVtxy'], ch.GenPart_vy[i])
                                Fill1D(histos['RecoFixedExtrapolEleVtxz'], ch.GenPart_vz[i])
                                Fill1D(histos['RecoFixedExtrapolEleDxy'], GenPart_dxy)
                                Fill1D(histos['RecoFixedExtrapolEleDz'], GenPart_dz)
                                Fill1D(histos['RecoFixedExtrapolEleDxyAbs'], abs(GenPart_dxy))
                                Fill1D(histos['RecoFixedExtrapolEleDzAbs'], abs(GenPart_dz))

                                histos['RecoFixedExtrapolPtDxy2D'].Fill( ch.GenPart_pt[i], abs(GenPart_dxy) )
        


                # just isotracks
                isodist = 100000
                isoIdx = -1
                idx0 = -1
                
                for j in range(ch.nIsoTrack):
                    dist0 = dR(ch.GenPart_eta[i], ch.IsoTrack_eta[j], ch.GenPart_phi[i], ch.IsoTrack_phi[j] ) 
                    if( dist0 < isodist):
                        isodist = dist0
                        idx0 = j
                if( dRcut(isodist,0.2) ):
                    Fill1D(histos['IsoTrackPt'], ch.GenPart_pt[i])
                    Fill1D(histos['IsoTrackEta'], ch.GenPart_eta[i])
                    Fill1D(histos['IsoTrackPhi'], ch.GenPart_phi[i])

                    if(bUseVertex):
                        Fill1D(histos['IsoTrackVtxx'], ch.GenPart_vx[i])
                        Fill1D(histos['IsoTrackVtxy'], ch.GenPart_vy[i])
                        Fill1D(histos['IsoTrackVtxz'], ch.GenPart_vz[i])
                        Fill1D(histos['IsoTrackDxy'], GenPart_dxy)
                        Fill1D(histos['IsoTrackDz'], GenPart_dz)
                        Fill1D(histos['IsoTrackDxyAbs'], abs(GenPart_dxy))
                        Fill1D(histos['IsoTrackDzAbs'], abs(GenPart_dz))
                    isoIdx = idx0


                # just photon
                phodist = 100000
                iPhoIdx = -1
                for j in range(ch.nPhoton):
                    dist0 = dR(ch.GenPart_eta[i], ch.Photon_eta[j], ch.GenPart_phi[i], ch.Photon_phi[j] ) 
                    if( dist0 < phodist):
                        phodist = dist0
                        iPhoIdx = j
                if( dRcut(phodist,0.2) ):
                    Fill1D(histos['RecoPhotonPt'], ch.GenPart_pt[i])
                    Fill1D(histos['RecoPhotonEta'], ch.GenPart_eta[i])
                    Fill1D(histos['RecoPhotonPhi'], ch.GenPart_phi[i])
                    if(bUseVertex):
                        Fill1D(histos['RecoPhotonVtxx'], ch.GenPart_vx[i])
                        Fill1D(histos['RecoPhotonVtxy'], ch.GenPart_vy[i])
                        Fill1D(histos['RecoPhotonVtxz'], ch.GenPart_vz[i])
                        Fill1D(histos['RecoPhotonDxy'], GenPart_dxy)
                        Fill1D(histos['RecoPhotonDz'], GenPart_dz)
                        Fill1D(histos['RecoPhotonDxyAbs'], abs(GenPart_dxy))
                        Fill1D(histos['RecoPhotonDzAbs'], abs(GenPart_dz))

                histos['dRControlPlot2D'].Fill(isodist,phodist)


        # background calculation
        # create background
        background = []
        for iTR in range(ch.nIsoTrack):
            min_dRele = 1000
            minIEle = -1
            min_dRbackgnd = 1000
            minIBackgnd = -1

            for iMC in range(ch.nGenPart):
                temp_dR = dR(ch.GenPart_eta[iMC], ch.IsoTrack_eta[iTR], ch.GenPart_phi[iMC], ch.IsoTrack_phi[iTR])
                if abs(ch.GenPart_pdgId[iMC]) == 11:
                    if( temp_dR < min_dRele):
                        min_dRele = temp_dR
                        minIEle = iMC
                else:
                    if(temp_dR < min_dRbackgnd):
                        min_dRbackgnd = temp_dR
                        minIBackgnd = iMC

            # FIXME: Ide nem kellene dR cut?????
            #if( minIBackgnd != -1):
            if( dRcut(min_dRbackgnd,0.2)):
                if(not dRcut(min_dRele,0.2) ):
                    background.append(minIBackgnd)

                    # this is the denom:
                    histos['TrueBackgroundPt'].Fill( ch.GenPart_pt[minIBackgnd])
                    histos['TrueBackgroundEta'].Fill( ch.GenPart_eta[minIBackgnd])
                    histos['TrueBackgroundPhi'].Fill( ch.GenPart_phi[minIBackgnd])
                    if(bUseVertex):
                        histos['TrueBackgroundDxyAbs'].Fill( abs(getdxy(minIBackgnd, ch)))
                        histos['TrueBackgroundDzAbs'].Fill( abs(getdz(minIBackgnd, ch)))

        # To get background rejection: redo the matching, but work on the true background as input sample:
        for iBK in range(len(background)):

            i = background[iBK]
            
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
            # not dRcut: background REJECTION efficiency
            if( not dRcut(eledist) ):
                histos['RejectBackgroundPt'].Fill( ch.GenPart_pt[i])
                histos['RejectBackgroundEta'].Fill( ch.GenPart_eta[i])
                histos['RejectBackgroundPhi'].Fill( ch.GenPart_phi[i])
                if(bUseVertex):
                    histos['RejectBackgroundDxyAbs'].Fill( abs(getdxy(i, ch)))
                    histos['RejectBackgroundDzAbs'].Fill( abs(getdz(i, ch)))
            else:
                eleIdx = idx0 


            # attempt reco fix
            if( eleIdx == -1 ):
                # match isotracks
                isodist = 100000
                isoIdx = -1
                idx0 = -1
                for j in range(ch.nIsoTrack):
                    dist0 = dR(ch.GenPart_eta[i], ch.IsoTrack_eta[j], ch.GenPart_phi[i], ch.IsoTrack_phi[j] ) 
                    if( dist0 < isodist):
                        isodist = dist0
                        idx0 = j
                if( dRcut(isodist,0.2) ):
                    isoIdx = idx0

                # reco ele not found - check if there is a reco-ed photon
                there_is_reco_photon = False
                for j in range(ch.nPhoton):
                    dist0 = dR(ch.GenPart_eta[i], ch.Photon_eta[j], ch.GenPart_phi[i], ch.Photon_phi[j] )
                    if( dRcut(dist0, 0.2) ):
                        there_is_reco_photon = True
                        break
                
                if( not there_is_reco_photon or isoIdx == -1):
                    histos['RejectBackgroundRecoFixPt'].Fill( ch.GenPart_pt[i])
                    histos['RejectBackgroundRecoFixEta'].Fill( ch.GenPart_eta[i])
                    histos['RejectBackgroundRecoFixPhi'].Fill( ch.GenPart_phi[i])
                    if(bUseVertex):
                        histos['RejectBackgroundRecoFixDxyAbs'].Fill( abs(getdxy(i, ch)))
                        histos['RejectBackgroundRecoFixDzAbs'].Fill( abs(getdz(i, ch)))
    # Save histos:
    extrapolate_log.close()
    hfile.Write()


# Make samplechain
ch = SampleChain(sample, options.startfile, options.nfiles, year).getchain()

n_entries = ch.GetEntries() #num entries in tree
nevtcut = n_entries -1 if nEvents == - 1 else nEvents - 1


if __name__ == '__main__':
    print "RecoFixThreaded calculating efficiencies with and without recofix on",Nthreads,"threads."
    proc = {}
    for i in range(Nthreads):
        it0 = 0 + i * n_entries / Nthreads ;
        it1 = 0 + (i+1) * n_entries / Nthreads ;
        proc[i] = multiprocessing.Process(target=RecoAndRecoFix, args=(i, it0, it1, nevtcut, ch))
        proc[i].start()
        #print 'Started thread',i
else:
    print("Something is serously wrong.")
    quit()

for i in range(Nthreads):
        proc[i].join()

print("All threads finished. Stacking results...")
# Stack partial results


ptBinning = [0, 3.5, 5, 7.5, 10, 12, 14, 16, 18, 20, 30, 40, 50, 100]
histext = ''
hfile = ROOT.TFile( 'root_files/RecoFixThreadedSTACK_Sample'+sample+'.root', 'RECREATE')
savehistos = {}
savehistos['TrueElePt'] = HistInfo(hname = 'TrueElePt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()
savehistos['TrueEleEta'] = HistInfo(hname = 'TrueEleEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['TrueElePhi'] = HistInfo(hname = 'TrueElePhi', sample = histext, binning = [30, -4, 4], histclass = ROOT.TH1F).make_hist()
savehistos['TrueEleVtxx'] = HistInfo(hname = 'TrueEleVtxx', sample = histext, binning = [40, -1, 1], histclass = ROOT.TH1F).make_hist()
savehistos['TrueEleVtxy'] = HistInfo(hname = 'TrueEleVtxy', sample = histext, binning = [40, -1, 1], histclass = ROOT.TH1F).make_hist()
savehistos['TrueEleVtxz'] = HistInfo(hname = 'TrueEleVtxz', sample = histext, binning = [50, -50, 50], histclass = ROOT.TH1F).make_hist()
savehistos['TrueEleDxy'] = HistInfo(hname = 'TrueEleDxy', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
savehistos['TrueEleDz'] = HistInfo(hname = 'TrueEleDz', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
savehistos['TrueEleDxyAbs'] = HistInfo(hname = 'TrueEleDxyAbs', sample = histext, binning = [40, 0, 15], histclass = ROOT.TH1F).make_hist()
savehistos['TrueEleDzAbs'] = HistInfo(hname = 'TrueEleDzAbs', sample = histext, binning = [40, 0, 15], histclass = ROOT.TH1F).make_hist()

savehistos['RecoElePt'] = HistInfo(hname = 'RecoElePt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()
savehistos['RecoEleEta'] = HistInfo(hname = 'RecoEleEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['RecoElePhi'] = HistInfo(hname = 'RecoElePhi', sample = histext, binning = [30, -4, 4], histclass = ROOT.TH1F).make_hist()
savehistos['RecoEleVtxx'] = HistInfo(hname = 'RecoEleVtxx', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['RecoEleVtxy'] = HistInfo(hname = 'RecoEleVtxy', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['RecoEleVtxz'] = HistInfo(hname = 'RecoEleVtxz', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['RecoEleDxy'] = HistInfo(hname = 'RecoEleDxy', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
savehistos['RecoEleDz'] = HistInfo(hname = 'RecoEleDz', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
savehistos['RecoEleDxyAbs'] = HistInfo(hname = 'RecoEleDxyAbs', sample = histext, binning = [40, 0, 15], histclass = ROOT.TH1F).make_hist()
savehistos['RecoEleDzAbs'] = HistInfo(hname = 'RecoEleDzAbs', sample = histext, binning = [40, 0, 15], histclass = ROOT.TH1F).make_hist()

savehistos['IsoTrackPt'] = HistInfo(hname = 'IsoTrackPt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()
savehistos['IsoTrackEta'] = HistInfo(hname = 'IsoTrackEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['IsoTrackPhi'] = HistInfo(hname = 'IsoTrackPhi', sample = histext, binning = [30, -4, 4], histclass = ROOT.TH1F).make_hist()
savehistos['IsoTrackVtxx'] = HistInfo(hname = 'IsoTrackVtxx', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['IsoTrackVtxy'] = HistInfo(hname = 'IsoTrackVtxy', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['IsoTrackVtxz'] = HistInfo(hname = 'IsoTrackVtxz', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['IsoTrackDxy'] = HistInfo(hname = 'IsoTrackDxy', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
savehistos['IsoTrackDz'] = HistInfo(hname = 'IsoTrackDz', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
savehistos['IsoTrackDxyAbs'] = HistInfo(hname = 'IsoTrackDxyAbs', sample = histext, binning = [40, 0, 15], histclass = ROOT.TH1F).make_hist()
savehistos['IsoTrackDzAbs'] = HistInfo(hname = 'IsoTrackDzAbs', sample = histext, binning = [40, 0, 15], histclass = ROOT.TH1F).make_hist()

savehistos['RecoFixedElePt'] = HistInfo(hname = 'RecoFixedElePt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()
savehistos['RecoFixedEleEta'] = HistInfo(hname = 'RecoFixedEleEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['RecoFixedElePhi'] = HistInfo(hname = 'RecoFixedElePhi', sample = histext, binning = [30, -4, 4], histclass = ROOT.TH1F).make_hist()
savehistos['RecoFixedEleVtxx'] = HistInfo(hname = 'RecoFixedEleVtxx', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['RecoFixedEleVtxy'] = HistInfo(hname = 'RecoFixedEleVtxy', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['RecoFixedEleVtxz'] = HistInfo(hname = 'RecoFixedEleVtxz', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['RecoFixedEleDxy'] = HistInfo(hname = 'RecoFixedEleDxy', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
savehistos['RecoFixedEleDz'] = HistInfo(hname = 'RecoFixedEleDz', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
savehistos['RecoFixedEleDxyAbs'] = HistInfo(hname = 'RecoFixedEleDxyAbs', sample = histext, binning = [40, 0, 15], histclass = ROOT.TH1F).make_hist()
savehistos['RecoFixedEleDzAbs'] = HistInfo(hname = 'RecoFixedEleDzAbs', sample = histext, binning = [40, 0, 15], histclass = ROOT.TH1F).make_hist()


savehistos['RecoFixedExtrapolElePt'] = HistInfo(hname = 'RecoFixedExtrapolElePt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()
savehistos['RecoFixedExtrapolEleEta'] = HistInfo(hname = 'RecoFixedExtrapolEleEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['RecoFixedExtrapolElePhi'] = HistInfo(hname = 'RecoFixedExtrapolElePhi', sample = histext, binning = [30, -4, 4], histclass = ROOT.TH1F).make_hist()
savehistos['RecoFixedExtrapolEleVtxx'] = HistInfo(hname = 'RecoFixedExtrapolEleVtxx', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['RecoFixedExtrapolEleVtxy'] = HistInfo(hname = 'RecoFixedExtrapolEleVtxy', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['RecoFixedExtrapolEleVtxz'] = HistInfo(hname = 'RecoFixedExtrapolEleVtxz', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['RecoFixedExtrapolEleDxy'] = HistInfo(hname = 'RecoFixedExtrapolEleDxy', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
savehistos['RecoFixedExtrapolEleDz'] = HistInfo(hname = 'RecoFixedExtrapolEleDz', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
savehistos['RecoFixedExtrapolEleDxyAbs'] = HistInfo(hname = 'RecoFixedExtrapolEleDxyAbs', sample = histext, binning = [40, 0, 15], histclass = ROOT.TH1F).make_hist()
savehistos['RecoFixedExtrapolEleDzAbs'] = HistInfo(hname = 'RecoFixedExtrapolEleDzAbs', sample = histext, binning = [40, 0, 15], histclass = ROOT.TH1F).make_hist()


savehistos['RecoFixedExtrapolPtDxy2D'] = ROOT.TH2F('RecoFixedExtrapolPtDxy2D_','RecoFixedExtrapolPtDxy2D', len(ptBinning) -1, np.array(ptBinning), 40,0,15)
savehistos['RecoFixedPtDxy2D'] = ROOT.TH2F('RecoFixedPtDxy2D_','RecoFixedPtDxy2D', len(ptBinning) -1, np.array(ptBinning), 40,0,15)

savehistos['TruePhotonPt'] = HistInfo(hname = 'TruePhotonPt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()
savehistos['TruePhotonEta'] = HistInfo(hname = 'TruePhotonEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['TruePhotonPhi'] = HistInfo(hname = 'TruePhotonPhi', sample = histext, binning = [30, -4, 4], histclass = ROOT.TH1F).make_hist()
savehistos['TruePhotonVtxx'] = HistInfo(hname = 'TruePhotonVtxx', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['TruePhotonVtxy'] = HistInfo(hname = 'TruePhotonVtxy', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['TruePhotonVtxz'] = HistInfo(hname = 'TruePhotonVtxz', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['TruePhotonDxy'] = HistInfo(hname = 'TruePhotonDxy', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
savehistos['TruePhotonDz'] = HistInfo(hname = 'TruePhotonDz', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
savehistos['TruePhotonDxyAbs'] = HistInfo(hname = 'TruePhotonDxyAbs', sample = histext, binning = [40, 0, 15], histclass = ROOT.TH1F).make_hist()
savehistos['TruePhotonDzAbs'] = HistInfo(hname = 'TruePhotonDzAbs', sample = histext, binning = [40, 0, 15], histclass = ROOT.TH1F).make_hist()

savehistos['RecoPhotonPt'] = HistInfo(hname = 'RecoPhotonPt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()
savehistos['RecoPhotonEta'] = HistInfo(hname = 'RecoPhotonEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['RecoPhotonPhi'] = HistInfo(hname = 'RecoPhotonPhi', sample = histext, binning = [30, -4, 4], histclass = ROOT.TH1F).make_hist()
savehistos['RecoPhotonVtxx'] = HistInfo(hname = 'RecoPhotonVtxx', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['RecoPhotonVtxy'] = HistInfo(hname = 'RecoPhotonVtxy', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['RecoPhotonVtxz'] = HistInfo(hname = 'RecoPhotonVtxz', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['RecoPhotonDxy'] = HistInfo(hname = 'RecoPhotonDxy', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
savehistos['RecoPhotonDz'] = HistInfo(hname = 'RecoPhotonDz', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
savehistos['RecoPhotonDxyAbs'] = HistInfo(hname = 'RecoPhotonDxyAbs', sample = histext, binning = [40, 0, 15], histclass = ROOT.TH1F).make_hist()
savehistos['RecoPhotonDzAbs'] = HistInfo(hname = 'RecoPhotonDzAbs', sample = histext, binning = [40, 0, 15], histclass = ROOT.TH1F).make_hist()


savehistos['RecoPhotonPlusIsoTrackPt'] = HistInfo(hname = 'RecoPhotonPlusIsoTrackPt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()
savehistos['RecoPhotonPlusIsoTrackEta'] = HistInfo(hname = 'RecoPhotonPlusIsoTrackEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['RecoPhotonPlusIsoTrackPhi'] = HistInfo(hname = 'RecoPhotonPlusIsoTrackPhi', sample = histext, binning = [30, -4, 4], histclass = ROOT.TH1F).make_hist()
savehistos['RecoPhotonPlusIsoTrackDxyAbs'] = HistInfo(hname = 'RecoPhotonPlusIsoTrackDxyAbs', sample = histext, binning = [40, 0, 15], histclass = ROOT.TH1F).make_hist()
savehistos['RecoPhotonPlusIsoTrackDzAbs'] = HistInfo(hname = 'RecoPhotonPlusIsoTrackDzAbs', sample = histext, binning = [40, 0, 15], histclass = ROOT.TH1F).make_hist()

savehistos['dRControlPlot2D'] = HistInfo(hname = 'dRControlPlot2D', sample = histext, binning = [[40, 0, 2],[40, 0, 2]], histclass = ROOT.TH2F).make_hist2D()

savehistos['TrueBackgroundPt'] = HistInfo(hname = 'TrueBackgroundPt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()
savehistos['TrueBackgroundEta'] = HistInfo(hname = 'TrueBackgroundEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['TrueBackgroundPhi'] = HistInfo(hname = 'TrueBackgroundPhi', sample = histext, binning = [30, -4, 4], histclass = ROOT.TH1F).make_hist()
savehistos['TrueBackgroundDxyAbs'] = HistInfo(hname = 'TrueBackgroundDxyAbs', sample = histext, binning = [40, 0, 15], histclass = ROOT.TH1F).make_hist()
savehistos['TrueBackgroundDzAbs'] = HistInfo(hname = 'TrueBackgroundDzAbs', sample = histext, binning = [40, 0, 15], histclass = ROOT.TH1F).make_hist()

savehistos['RejectBackgroundPt'] = HistInfo(hname = 'RejectBackgroundPt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()
savehistos['RejectBackgroundEta'] = HistInfo(hname = 'RejectBackgroundEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['RejectBackgroundPhi'] = HistInfo(hname = 'RejectBackgroundPhi', sample = histext, binning = [30, -4, 4], histclass = ROOT.TH1F).make_hist()
savehistos['RejectBackgroundDxyAbs'] = HistInfo(hname = 'RejectBackgroundDxyAbs', sample = histext, binning = [40, 0, 15], histclass = ROOT.TH1F).make_hist()
savehistos['RejectBackgroundDzAbs'] = HistInfo(hname = 'RejectBackgroundDzAbs', sample = histext, binning = [40, 0, 15], histclass = ROOT.TH1F).make_hist()

savehistos['RejectBackgroundRecoFixPt'] = HistInfo(hname = 'RejectBackgroundRecoFixPt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()
savehistos['RejectBackgroundRecoFixEta'] = HistInfo(hname = 'RejectBackgroundRecoFixEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['RejectBackgroundRecoFixPhi'] = HistInfo(hname = 'RejectBackgroundRecoFixPhi', sample = histext, binning = [30, -4, 4], histclass = ROOT.TH1F).make_hist()
savehistos['RejectBackgroundRecoFixDxyAbs'] = HistInfo(hname = 'RejectBackgroundRecoFixDxyAbs', sample = histext, binning = [40, 0, 15], histclass = ROOT.TH1F).make_hist()
savehistos['RejectBackgroundRecoFixDzAbs'] = HistInfo(hname = 'RejectBackgroundRecoFixDzAbs', sample = histext, binning = [40, 0, 15], histclass = ROOT.TH1F).make_hist()


savehistos['AllTruePhotonPt'] = HistInfo(hname = 'AllTruePhotonPt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()
savehistos['AllRecoPhotonPt'] = HistInfo(hname = 'AllRecoPhotonPt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()
for threadID in range(Nthreads):
    savedfiles = ROOT.TFile.Open('root_files/RecoFixThreaded_Sample'+sample+'_Thread-'+str(threadID)+'.root')

    for key in savehistos.keys():
        savehistos[key].Add( savedfiles.Get( key+'_' ))


print("Saving stacked root file...")
# Save endresult:
# FIXME: DONT SAVE
hfile.Write()
