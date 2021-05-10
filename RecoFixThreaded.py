import os, sys
import ROOT
import types
import math
import multiprocessing
import threading

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
    argParser.add_argument('--sample',           action='store',                     type=str,            default='UL17_Full99mm',                                help="Which sample?" )
    argParser.add_argument('--year',             action='store',                     type=int,            default=2016,                                             help="Which year?" )
    argParser.add_argument('--startfile',        action='store',                     type=int,            default=0,                                                help="start from which root file like 0th or 10th etc?" )
    argParser.add_argument('--nfiles',           action='store',                     type=int,            default=-1,                                               help="No of files to run. -1 means all files" )
    argParser.add_argument('--nevents',           action='store',                    type=int,            default=-1,                                               help="No of events to run. -1 means all events" )
    argParser.add_argument('--Nthreads',           action='store',                    type=int,            default=-1,                                               help="How many CPU?" )
    argParser.add_argument('--channel',           action='store',                    type=str,            default='Electron',                                       help="Which lepton?" )
    argParser.add_argument('--region',            action='store',                    type=str,            default='mesurement',                                     help="Which lepton?" )    
    argParser.add_argument('--pJobs',             action='store',                    type=bool,            default=False,                                           help="using GPU parallel program or not" )

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
    R = pt / (0.3 * 3.8)
    xC = x0/100.0 + R* math.cos(phi - charge * 3.14159265359/2.0)
    yC = y0/100.0 + R* math.sin(phi - charge * 3.14159265359/2.0)


    # calculate x,y intersection of track
    a = (-1* yC) / xC
    rb = 129.0 / 100.0
    RC2 = xC**2 + yC**2
    b = (RC2 - R**2 + rb**2) / (2*xC)

    qa = a**2 + 1
    qb = 2*a*b
    qc = b**2 - rb
    disc = b**2 - 4*qa*qc

    y,x,y_other,x_other = 0,0,0,0
    if( disc > 0):
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
        Z_ECAL = z0/100.0 - (math.asin((y-yC)/R) - phi - charge*3.14159265359/2.0) * (R*math.sinh(eta) / charge)

        if(Z_ECAL > 2.1):
            Z_ECAL = 2.1
        if(Z_ECAL < -2.1):
            Z_ECAL = -2.1
    else:
        # Barrel cannot be hit, endcap is hit (electron spirals out)
        if(eta > 0):
            Z_ECAL = 2.1
        else:
            Z_ECAL = -2.1

    X_ECAL = xC + R * math.cos(-charge * (Z_ECAL-z0/100.0)/(R * math.sinh(eta)) + charge * 3.14159265359/2.0)
    Y_ECAL = yC + R * math.sin(-charge * (Z_ECAL-z0/100.0)/(R * math.sinh(eta)) + charge * 3.14159265359/2.0)
  
    D_ECAL = math.sqrt(X_ECAL*X_ECAL+Y_ECAL*Y_ECAL)

    etaSC = math.asinh(Z_ECAL/D_ECAL)
    if(Y_ECAL > 0):
        phiSC = math.cos(X_ECAL/D_ECAL)
    else:
        phiSC = -1.0 * math.cos(X_ECAL/D_ECAL)


    return etaSC, phiSC, x, y, x_other, y_other

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
    return ((ch.GenPart_vy[i] - ch.GenVtx_y)  * px - (ch.GenPart_vx[i] - ch.GenVtx_x) * py) / ch.GenPart_pt[i];

def getdz(i, ch):
    px = ch.GenPart_pt[i] * cos(ch.GenPart_phi[i])
    py = ch.GenPart_pt[i] * sin(ch.GenPart_phi[i])
    pz = ch.GenPart_pt[i] / math.tan(ch.GenPart_phi[i])
    return (ch.GenPart_vz[i] - ch.GenVtx_z) - ((ch.GenPart_vx[i] - ch.GenVtx_x) * px + (ch.GenPart_vy[i] - ch.GenVtx_y)*  py) / ch.GenPart_pt[i] * (pz / ch.GenPart_pt[i]);


# ele has larger impact parameter
# reco eff is low
# try improving

# RECO #
def RecoAndRecoFix(threadID, it0, it1, nevtcut, ch_common):

# Output histo - own for each Process
    ptBinning = [0, 3.5, 5, 7.5, 10, 12, 14, 16, 18, 20, 30, 40, 50, 100]
    histext = ''
    hfile = ROOT.TFile( 'root_files/RecoFixThreaded_Sample'+sample+'_Thread-'+str(threadID)+'.root', 'RECREATE')
    histos = {}
    histos['TrueElePt'] = HistInfo(hname = 'TrueElePt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()
    histos['TrueEleEta'] = HistInfo(hname = 'TrueEleEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['TrueEleVtxx'] = HistInfo(hname = 'TrueEleVtxx', sample = histext, binning = [40, -1, 1], histclass = ROOT.TH1F).make_hist()
    histos['TrueEleVtxy'] = HistInfo(hname = 'TrueEleVtxy', sample = histext, binning = [40, -1, 1], histclass = ROOT.TH1F).make_hist()
    histos['TrueEleVtxz'] = HistInfo(hname = 'TrueEleVtxz', sample = histext, binning = [50, -50, 50], histclass = ROOT.TH1F).make_hist()
    histos['TrueEleDxy'] = HistInfo(hname = 'TrueEleDxy', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
    histos['TrueEleDz'] = HistInfo(hname = 'TrueEleDz', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
    histos['TrueEleDxyAbs'] = HistInfo(hname = 'TrueEleDxyAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()
    histos['TrueEleDzAbs'] = HistInfo(hname = 'TrueEleDzAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()

    histos['RecoElePt'] = HistInfo(hname = 'RecoElePt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()
    histos['RecoEleEta'] = HistInfo(hname = 'RecoEleEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['RecoEleVtxx'] = HistInfo(hname = 'RecoEleVtxx', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['RecoEleVtxy'] = HistInfo(hname = 'RecoEleVtxy', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['RecoEleVtxz'] = HistInfo(hname = 'RecoEleVtxz', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['RecoEleDxy'] = HistInfo(hname = 'RecoEleDxy', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
    histos['RecoEleDz'] = HistInfo(hname = 'RecoEleDz', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
    histos['RecoEleDxyAbs'] = HistInfo(hname = 'RecoEleDxyAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()
    histos['RecoEleDzAbs'] = HistInfo(hname = 'RecoEleDzAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()

    histos['IsoTrackPt'] = HistInfo(hname = 'IsoTrackPt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()
    histos['IsoTrackEta'] = HistInfo(hname = 'IsoTrackEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['IsoTrackVtxx'] = HistInfo(hname = 'IsoTrackVtxx', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['IsoTrackVtxy'] = HistInfo(hname = 'IsoTrackVtxy', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['IsoTrackVtxz'] = HistInfo(hname = 'IsoTrackVtxz', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['IsoTrackDxy'] = HistInfo(hname = 'IsoTrackDxy', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
    histos['IsoTrackDz'] = HistInfo(hname = 'IsoTrackDz', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
    histos['IsoTrackDxyAbs'] = HistInfo(hname = 'IsoTrackDxyAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()
    histos['IsoTrackDzAbs'] = HistInfo(hname = 'IsoTrackDzAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()

    histos['RecoFixedElePt'] = HistInfo(hname = 'RecoFixedElePt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()
    histos['RecoFixedEleEta'] = HistInfo(hname = 'RecoFixedEleEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['RecoFixedEleVtxx'] = HistInfo(hname = 'RecoFixedEleVtxx', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['RecoFixedEleVtxy'] = HistInfo(hname = 'RecoFixedEleVtxy', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['RecoFixedEleVtxz'] = HistInfo(hname = 'RecoFixedEleVtxz', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['RecoFixedEleDxy'] = HistInfo(hname = 'RecoFixedEleDxy', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
    histos['RecoFixedEleDz'] = HistInfo(hname = 'RecoFixedEleDz', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
    histos['RecoFixedEleDxyAbs'] = HistInfo(hname = 'RecoFixedEleDxyAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()
    histos['RecoFixedEleDzAbs'] = HistInfo(hname = 'RecoFixedEleDzAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()

    histos['TruePhotonPt'] = HistInfo(hname = 'TruePhotonPt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()
    histos['TruePhotonEta'] = HistInfo(hname = 'TruePhotonEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['TruePhotonVtxx'] = HistInfo(hname = 'TruePhotonVtxx', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['TruePhotonVtxy'] = HistInfo(hname = 'TruePhotonVtxy', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['TruePhotonVtxz'] = HistInfo(hname = 'TruePhotonVtxz', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['TruePhotonDxy'] = HistInfo(hname = 'TruePhotonDxy', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
    histos['TruePhotonDz'] = HistInfo(hname = 'TruePhotonDz', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
    histos['TruePhotonDxyAbs'] = HistInfo(hname = 'TruePhotonDxyAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()
    histos['TruePhotonDzAbs'] = HistInfo(hname = 'TruePhotonDzAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()

    histos['RecoPhotonPt'] = HistInfo(hname = 'RecoPhotonPt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()
    histos['RecoPhotonEta'] = HistInfo(hname = 'RecoPhotonEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['RecoPhotonVtxx'] = HistInfo(hname = 'RecoPhotonVtxx', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['RecoPhotonVtxy'] = HistInfo(hname = 'RecoPhotonVtxy', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['RecoPhotonVtxz'] = HistInfo(hname = 'RecoPhotonVtxz', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['RecoPhotonDxy'] = HistInfo(hname = 'RecoPhotonDxy', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
    histos['RecoPhotonDz'] = HistInfo(hname = 'RecoPhotonDz', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
    histos['RecoPhotonDxyAbs'] = HistInfo(hname = 'RecoPhotonDxyAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()
    histos['RecoPhotonDzAbs'] = HistInfo(hname = 'RecoPhotonDzAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()

    histos['RecoPhotonPlusIsoTrackDxyAbs'] = HistInfo(hname = 'RecoPhotonPlusIsoTrackDxyAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()
    histos['RecoPhotonPlusIsoTrackDzAbs'] = HistInfo(hname = 'RecoPhotonPlusIsoTrackDzAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()

    histos['dRControlPlot2D'] = HistInfo(hname = 'dRControlPlot2D', sample = histext, binning = [[40, 0, 2],[40, 0, 2]], histclass = ROOT.TH2F).make_hist2D()

    histos['TrueBackgroundPt'] = HistInfo(hname = 'TrueBackgroundPt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()
    histos['TrueBackgroundEta'] = HistInfo(hname = 'TrueBackgroundEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['TrueBackgroundDxyAbs'] = HistInfo(hname = 'TrueBackgroundDxyAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()
    histos['TrueBackgroundDzAbs'] = HistInfo(hname = 'TrueBackgroundDzAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()

    histos['RejectBackgroundPt'] = HistInfo(hname = 'RejectBackgroundPt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()
    histos['RejectBackgroundEta'] = HistInfo(hname = 'RejectBackgroundEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['RejectBackgroundDxyAbs'] = HistInfo(hname = 'RejectBackgroundDxyAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()
    histos['RejectBackgroundDzAbs'] = HistInfo(hname = 'RejectBackgroundDzAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()

    histos['RejectBackgroundRecoFixPt'] = HistInfo(hname = 'RejectBackgroundRecoFixPt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()
    histos['RejectBackgroundRecoFixEta'] = HistInfo(hname = 'RejectBackgroundRecoFixEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
    histos['RejectBackgroundRecoFixDxyAbs'] = HistInfo(hname = 'RejectBackgroundRecoFixDxyAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()
    histos['RejectBackgroundRecoFixDzAbs'] = HistInfo(hname = 'RejectBackgroundRecoFixDzAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()

    histos['AllTruePhotonPt'] = HistInfo(hname = 'AllTruePhotonPt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()
    histos['AllRecoPhotonPt'] = HistInfo(hname = 'AllRecoPhotonPt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()

    interpolate_log = open("logs/interpolate-log-t"+str(threadID),"w")
    # RECO CALCULATION #
    for ientry in range(it0,it1,1): #loop over events
        if ientry > nevtcut: break
        if ientry == it0 : print 'Process-'+str(threadID)+' starting at ', ientry,'th event'
        if (ientry % (nevtcut/10)==0 and ientry != it0 and ientry != it1 - 1) : print 'Process-'+str(threadID)+' processing ', ientry,'th event'
        if ientry == it1 - 1: print 'Process-'+str(threadID)+' finishing at ', ientry,'th event'

        ch = ch_common
        ch.GetEntry(ientry)


        for i in range(ch.nPhoton):
            Fill1D(histos['AllRecoPhotonPt'], ch.Photon_pt[i])


        for i in range(ch.nGenPart):

            ## REMOVEME photon pt test
            if(abs(ch.GenPart_pdgId[i]) == 22):
                Fill1D(histos['AllTruePhotonPt'], ch.GenPart_pt[i])




            # electron, pt eta cuts, status 1, with a stop mother somewhere in mother history
            # FIXME: TESTING NOW
            if( abs(ch.GenPart_pdgId[i]) == 11 and ch.GenPart_pt[i] >= 15  and abs(ch.GenPart_eta[i]) <= 1 and ch.GenPart_status[i] == 1 and hasMomRecursive(i, 1000006, ch) and abs(ch.GenPart_vx[i]) < 0.2 and abs(ch.GenPart_vy[i]) < 0.2):

                GenPart_dxy = getdxy(i, ch)
                GenPart_dz = getdz(i, ch)

                Fill1D(histos['TrueElePt'],ch.GenPart_pt[i])
                Fill1D(histos['TrueEleEta'],ch.GenPart_eta[i])
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
                    Fill1D(histos['RecoEleVtxx'], ch.GenPart_vx[i])
                    Fill1D(histos['RecoEleVtxy'], ch.GenPart_vy[i])
                    Fill1D(histos['RecoEleVtxz'], ch.GenPart_vz[i])
                    Fill1D(histos['RecoEleDxy'], GenPart_dxy)
                    Fill1D(histos['RecoEleDz'], GenPart_dz)
                    Fill1D(histos['RecoEleDxyAbs'], abs(GenPart_dxy))
                    Fill1D(histos['RecoEleDzAbs'], abs(GenPart_dz))

                    Fill1D(histos['RecoFixedElePt'], ch.GenPart_pt[i])
                    Fill1D(histos['RecoFixedEleEta'], ch.GenPart_eta[i])
                    Fill1D(histos['RecoFixedEleVtxx'], ch.GenPart_vx[i])
                    Fill1D(histos['RecoFixedEleVtxy'], ch.GenPart_vy[i])
                    Fill1D(histos['RecoFixedEleVtxz'], ch.GenPart_vz[i])
                    Fill1D(histos['RecoFixedEleDxy'], GenPart_dxy)
                    Fill1D(histos['RecoFixedEleDz'], GenPart_dz)
                    Fill1D(histos['RecoFixedEleDxyAbs'], abs(GenPart_dxy))
                    Fill1D(histos['RecoFixedEleDzAbs'], abs(GenPart_dz))


                    etaSC_i, phiSC_i, calc_x, calc_y, calc2_x, calc2_y = GetGenPartSC(i, ch)
                    tangent_x, tangent_y = phiToXY(ch.GenPart_phi[i])
                    sc_x, sc_y = phiToXY(phiSC_i)
                    eSCx, eSCy = phiToXY(ch.Electron_phi[eleIdx])
                    print >> interpolate_log, "RECO ELE FOUND"
                    print >> interpolate_log, "Gen Ele P_T, x, y, z, charge:     ", ch.GenPart_pt[i],  ch.GenPart_vx[i], ch.GenPart_vy[i], ch.GenPart_vz[i], -1*math.copysign(1,ch.GenPart_pdgId[i])
                    print >> interpolate_log, "Gen Ele  eta, phi:                ", ch.GenPart_eta[i], ch.GenPart_phi[i]
                    print >> interpolate_log, "SC eta, phi:                      ", etaSC_i, phiSC_i
                    print >> interpolate_log, "Delta eta, Delta phi :            ", ch.GenPart_eta[i] - etaSC_i, ch.GenPart_phi[i] - phiSC_i
                    print >> interpolate_log, "Tangent crossing   (x,y), phi:    ", tangent_x, tangent_y, ch.GenPart_phi[i]
                    print >> interpolate_log, "Interpl cross SEL (x1,y1), phi:   ", calc_x, calc_y, math.atan2(calc_y, calc_x) 
                    print >> interpolate_log, "Interpl cross other (x2,y2), phi: ", calc2_x, calc2_y, math.atan2(calc2_y, calc2_x) 
                    print >> interpolate_log, "  RecoEle x y     (xSC,ySC),phi:  ", eSCx, eSCy, ch.Electron_phi[eleIdx]
                    print >> interpolate_log, " "


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
                    Fill1D(histos['IsoTrackPt'], ch.GenPart_pt[i])
                    Fill1D(histos['IsoTrackEta'], ch.GenPart_eta[i])
                    Fill1D(histos['IsoTrackVtxx'], ch.GenPart_vx[i])
                    Fill1D(histos['IsoTrackVtxy'], ch.GenPart_vy[i])
                    Fill1D(histos['IsoTrackVtxz'], ch.GenPart_vz[i])
                    Fill1D(histos['IsoTrackDxy'], GenPart_dxy)
                    Fill1D(histos['IsoTrackDz'], GenPart_dz)
                    Fill1D(histos['IsoTrackDxyAbs'], abs(GenPart_dxy))
                    Fill1D(histos['IsoTrackDzAbs'], abs(GenPart_dz))
                    isoIdx = idx0

                # attempt reco fix
                if( eleIdx == -1 ):
                    # reco ele not found - check if there is a reco-ed photon

                    # when matching ele with photon, use extrapolated track:
                    etaSC_i, phiSC_i, calc_x, calc_y, calc2_x, calc2_y = GetGenPartSC(i, ch)
                    tangent_x, tangent_y = phiToXY(ch.GenPart_phi[i])
                    sc_x, sc_y = phiToXY(phiSC_i)
                    if( ch.nPhoton > 0):
                        print >> interpolate_log, "Gen Ele P_T, x, y, z, charge:     ", ch.GenPart_pt[i],  ch.GenPart_vx[i], ch.GenPart_vy[i], ch.GenPart_vz[i], -1*math.copysign(1,ch.GenPart_pdgId[i])
                        print >> interpolate_log, "Gen Ele  eta, phi:                ", ch.GenPart_eta[i], ch.GenPart_phi[i]
                        print >> interpolate_log, "SC eta, phi:                      ", etaSC_i, phiSC_i
                        print >> interpolate_log, "Delta eta, Delta phi :            ", ch.GenPart_eta[i] - etaSC_i, ch.GenPart_phi[i] - phiSC_i
                        print >> interpolate_log, "Tangent crossing   (x,y), phi:    ", tangent_x, tangent_y, ch.GenPart_phi[i]
                        print >> interpolate_log, "Interpl cross SEL (x1,y1), phi:   ", calc_x, calc_y, math.atan2(calc_y, calc_x) 
                        print >> interpolate_log, "Interpl cross other (x2,y2), phi: ", calc2_x, calc2_y, math.atan2(calc2_y, calc2_x) 

                    #print >> interpolate_log, "Supercluster   (yECAL,xECAL):", sc_y, sc_x
                    there_is_reco_photon = False
                    for j in range(ch.nPhoton):
                        dist0 = dR(etaSC_i, ch.Photon_eta[j], phiSC_i, ch.Photon_phi[j] )
                        pSCx, pSCy = phiToXY(ch.Photon_phi[j])
                        #print >> interpolate_log, "  Photoncand-"+str(j)+" eta, phi, pt: ", ch.Photon_eta[j], ch.Photon_phi[j], ch.Photon_pt[j]
                        print >> interpolate_log, "  Photocand-"+str(j)+" (xSC,ySC),phi:      ", pSCx, pSCy, ch.Photon_phi[j]
                        if( dRcut(dist0, 0.2) ):
                            there_is_reco_photon = True
                            break
                    
                    if( ch.nPhoton > 0):
                        print >> interpolate_log, "Photon found: ", there_is_reco_photon
                        print >> interpolate_log, " "

                    if( there_is_reco_photon):
                        # Recover efficiency by looking for photon, but this has big background
                        # Check if there is an isotrack for true ele -> then the isotrack is the detected ele
                        if( not isoIdx == -1):
                            Fill1D(histos['RecoFixedElePt'], ch.GenPart_pt[i])
                            Fill1D(histos['RecoFixedEleEta'], ch.GenPart_eta[i])
                            Fill1D(histos['RecoFixedEleVtxx'], ch.GenPart_vx[i])
                            Fill1D(histos['RecoFixedEleVtxy'], ch.GenPart_vy[i])
                            Fill1D(histos['RecoFixedEleVtxz'], ch.GenPart_vz[i])
                            Fill1D(histos['RecoFixedEleDxy'], GenPart_dxy)
                            Fill1D(histos['RecoFixedEleDz'], GenPart_dz)
                            Fill1D(histos['RecoFixedEleDxyAbs'], abs(GenPart_dxy))
                            Fill1D(histos['RecoFixedEleDzAbs'], abs(GenPart_dz))
            
                            # these are not filled when ele is found, just here:
                            Fill1D(histos['RecoPhotonPlusIsoTrackDxyAbs'], abs(GenPart_dxy))
                            Fill1D(histos['RecoPhotonPlusIsoTrackDzAbs'], abs(GenPart_dz))

                # photon reco eff
                phodist = 100000
                phoIdx = -1
                for j in range(ch.nPhoton):
                    dist0 = dR(ch.GenPart_eta[i], ch.Photon_eta[j], ch.GenPart_phi[i], ch.Photon_phi[j] ) 
                    if( dist0 < phodist):
                        phodist = dist0
                        phoIdx = j
                if( dRcut(phodist,0.2) ):
                    Fill1D(histos['RecoPhotonPt'], ch.GenPart_pt[i])
                    Fill1D(histos['RecoPhotonEta'], ch.GenPart_eta[i])
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
                    histos['RejectBackgroundRecoFixDxyAbs'].Fill( abs(getdxy(i, ch)))
                    histos['RejectBackgroundRecoFixDzAbs'].Fill( abs(getdz(i, ch)))
    # Save histos:
    interpolate_log.close()
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
savehistos['TrueEleVtxx'] = HistInfo(hname = 'TrueEleVtxx', sample = histext, binning = [40, -1, 1], histclass = ROOT.TH1F).make_hist()
savehistos['TrueEleVtxy'] = HistInfo(hname = 'TrueEleVtxy', sample = histext, binning = [40, -1, 1], histclass = ROOT.TH1F).make_hist()
savehistos['TrueEleVtxz'] = HistInfo(hname = 'TrueEleVtxz', sample = histext, binning = [50, -50, 50], histclass = ROOT.TH1F).make_hist()
savehistos['TrueEleDxy'] = HistInfo(hname = 'TrueEleDxy', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
savehistos['TrueEleDz'] = HistInfo(hname = 'TrueEleDz', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
savehistos['TrueEleDxyAbs'] = HistInfo(hname = 'TrueEleDxyAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()
savehistos['TrueEleDzAbs'] = HistInfo(hname = 'TrueEleDzAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()

savehistos['RecoElePt'] = HistInfo(hname = 'RecoElePt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()
savehistos['RecoEleEta'] = HistInfo(hname = 'RecoEleEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['RecoEleVtxx'] = HistInfo(hname = 'RecoEleVtxx', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['RecoEleVtxy'] = HistInfo(hname = 'RecoEleVtxy', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['RecoEleVtxz'] = HistInfo(hname = 'RecoEleVtxz', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['RecoEleDxy'] = HistInfo(hname = 'RecoEleDxy', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
savehistos['RecoEleDz'] = HistInfo(hname = 'RecoEleDz', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
savehistos['RecoEleDxyAbs'] = HistInfo(hname = 'RecoEleDxyAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()
savehistos['RecoEleDzAbs'] = HistInfo(hname = 'RecoEleDzAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()

savehistos['IsoTrackPt'] = HistInfo(hname = 'IsoTrackPt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()
savehistos['IsoTrackEta'] = HistInfo(hname = 'IsoTrackEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['IsoTrackVtxx'] = HistInfo(hname = 'IsoTrackVtxx', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['IsoTrackVtxy'] = HistInfo(hname = 'IsoTrackVtxy', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['IsoTrackVtxz'] = HistInfo(hname = 'IsoTrackVtxz', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['IsoTrackDxy'] = HistInfo(hname = 'IsoTrackDxy', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
savehistos['IsoTrackDz'] = HistInfo(hname = 'IsoTrackDz', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
savehistos['IsoTrackDxyAbs'] = HistInfo(hname = 'IsoTrackDxyAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()
savehistos['IsoTrackDzAbs'] = HistInfo(hname = 'IsoTrackDzAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()

savehistos['RecoFixedElePt'] = HistInfo(hname = 'RecoFixedElePt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()
savehistos['RecoFixedEleEta'] = HistInfo(hname = 'RecoFixedEleEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['RecoFixedEleVtxx'] = HistInfo(hname = 'RecoFixedEleVtxx', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['RecoFixedEleVtxy'] = HistInfo(hname = 'RecoFixedEleVtxy', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['RecoFixedEleVtxz'] = HistInfo(hname = 'RecoFixedEleVtxz', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['RecoFixedEleDxy'] = HistInfo(hname = 'RecoFixedEleDxy', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
savehistos['RecoFixedEleDz'] = HistInfo(hname = 'RecoFixedEleDz', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
savehistos['RecoFixedEleDxyAbs'] = HistInfo(hname = 'RecoFixedEleDxyAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()
savehistos['RecoFixedEleDzAbs'] = HistInfo(hname = 'RecoFixedEleDzAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()

savehistos['TruePhotonPt'] = HistInfo(hname = 'TruePhotonPt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()
savehistos['TruePhotonEta'] = HistInfo(hname = 'TruePhotonEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['TruePhotonVtxx'] = HistInfo(hname = 'TruePhotonVtxx', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['TruePhotonVtxy'] = HistInfo(hname = 'TruePhotonVtxy', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['TruePhotonVtxz'] = HistInfo(hname = 'TruePhotonVtxz', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['TruePhotonDxy'] = HistInfo(hname = 'TruePhotonDxy', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
savehistos['TruePhotonDz'] = HistInfo(hname = 'TruePhotonDz', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
savehistos['TruePhotonDxyAbs'] = HistInfo(hname = 'TruePhotonDxyAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()
savehistos['TruePhotonDzAbs'] = HistInfo(hname = 'TruePhotonDzAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()

savehistos['RecoPhotonPt'] = HistInfo(hname = 'RecoPhotonPt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()
savehistos['RecoPhotonEta'] = HistInfo(hname = 'RecoPhotonEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['RecoPhotonVtxx'] = HistInfo(hname = 'RecoPhotonVtxx', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['RecoPhotonVtxy'] = HistInfo(hname = 'RecoPhotonVtxy', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['RecoPhotonVtxz'] = HistInfo(hname = 'RecoPhotonVtxz', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['RecoPhotonDxy'] = HistInfo(hname = 'RecoPhotonDxy', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
savehistos['RecoPhotonDz'] = HistInfo(hname = 'RecoPhotonDz', sample = histext, binning = [40, -6, 6], histclass = ROOT.TH1F).make_hist()
savehistos['RecoPhotonDxyAbs'] = HistInfo(hname = 'RecoPhotonDxyAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()
savehistos['RecoPhotonDzAbs'] = HistInfo(hname = 'RecoPhotonDzAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()

savehistos['RecoPhotonPlusIsoTrackDxyAbs'] = HistInfo(hname = 'RecoPhotonPlusIsoTrackDxyAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()
savehistos['RecoPhotonPlusIsoTrackDzAbs'] = HistInfo(hname = 'RecoPhotonPlusIsoTrackDzAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()

savehistos['dRControlPlot2D'] = HistInfo(hname = 'dRControlPlot2D', sample = histext, binning = [[40, 0, 2],[40, 0, 2]], histclass = ROOT.TH2F).make_hist2D()

savehistos['TrueBackgroundPt'] = HistInfo(hname = 'TrueBackgroundPt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()
savehistos['TrueBackgroundEta'] = HistInfo(hname = 'TrueBackgroundEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['TrueBackgroundDxyAbs'] = HistInfo(hname = 'TrueBackgroundDxyAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()
savehistos['TrueBackgroundDzAbs'] = HistInfo(hname = 'TrueBackgroundDzAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()

savehistos['RejectBackgroundPt'] = HistInfo(hname = 'RejectBackgroundPt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()
savehistos['RejectBackgroundEta'] = HistInfo(hname = 'RejectBackgroundEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['RejectBackgroundDxyAbs'] = HistInfo(hname = 'RejectBackgroundDxyAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()
savehistos['RejectBackgroundDzAbs'] = HistInfo(hname = 'RejectBackgroundDzAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()

savehistos['RejectBackgroundRecoFixPt'] = HistInfo(hname = 'RejectBackgroundRecoFixPt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()
savehistos['RejectBackgroundRecoFixEta'] = HistInfo(hname = 'RejectBackgroundRecoFixEta', sample = histext, binning = [30, -3, 3], histclass = ROOT.TH1F).make_hist()
savehistos['RejectBackgroundRecoFixDxyAbs'] = HistInfo(hname = 'RejectBackgroundRecoFixDxyAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()
savehistos['RejectBackgroundRecoFixDzAbs'] = HistInfo(hname = 'RejectBackgroundRecoFixDzAbs', sample = histext, binning = [40, 0, 20], histclass = ROOT.TH1F).make_hist()


savehistos['AllTruePhotonPt'] = HistInfo(hname = 'AllTruePhotonPt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()
savehistos['AllRecoPhotonPt'] = HistInfo(hname = 'AllRecoPhotonPt', sample = histext, binning = ptBinning, binopt = 'var', histclass = ROOT.TH1F).make_hist()
for threadID in range(Nthreads):
    savedfiles = ROOT.TFile.Open('root_files/RecoFixThreaded_Sample'+sample+'_Thread-'+str(threadID)+'.root')

    for key in savehistos.keys():
        savehistos[key].Add( savedfiles.Get( key+'_' ))


print("Saving stacked root file...")
# Save endresult:
hfile.Write()