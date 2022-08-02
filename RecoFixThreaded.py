import os, sys
import ROOT
import types
import math
import multiprocessing
import threading

#from RegSel import RegSel
import numpy as np
import Histograms

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
    argParser.add_argument('--sample',           action='store',                     type=str,            default='Sig_Displaced_350_335',                                help="Which sample?" )
    argParser.add_argument('--signal',           action='store',                    type=str,            default="True",                                       help="calculate signal things" )
    argParser.add_argument('--background',             action='store',                    type=str,            default="True",                        help="calculate BK things")
    argParser.add_argument('--year',             action='store',                     type=str,            default="2016PostVFP",                                             help="Which year?" )
    argParser.add_argument('--startfile',        action='store',                     type=int,            default=0,                                                help="start from which root file like 0th or 10th etc?" )
    argParser.add_argument('--nfiles',           action='store',                     type=int,            default=-1,                                               help="No of files to run. -1 means all files" )
    argParser.add_argument('--nevents',           action='store',                    type=int,            default=-1,                                               help="No of events to run. -1 means all events" )
    argParser.add_argument('--Nthreads',           action='store',                    type=int,            default=-1,                                               help="How many CPU?" )
    argParser.add_argument('--channel',           action='store',                    type=str,            default='Electron',                                       help="Which lepton?" )
    argParser.add_argument('--ptcut',             action='store',                    type=float,            default=3.5,                                           help="pt lower cut")
    argParser.add_argument('--ptcutHigh',             action='store',                    type=float,            default=-1,                                           help="pt upper cut")
    argParser.add_argument('--vertex',             action='store',                    type=str,            default="True",                                           help="is there vertex info" )
    argParser.add_argument('--extrapolate',             action='store',                    type=str,            default="False",                                           help="extrapolate isotrack to ECAL surface" )
    argParser.add_argument('--v9',             action='store',                    type=str,            default="True",                                           help="use v9 nanoAOD")
    argParser.add_argument('--binning',             action='store',                    type=int,            default=0,                                        help="which set of binning to use")
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
binningSet = options.binning
if( options.signal == "True"):
    bSignal = True
else:
    bSignal = False

if( options.background == "True"):
    bBackground = True
else:
    bBackground = False


if( options.vertex == "True"):
    bUseVertex = True
else:
    bUseVertex = False

if( options.extrapolate == "True"):
    bExtrapolate = True
else:
    bExtrapolate = False

if( options.v9 == "True"):
    bUseV9 = True
else:
    bUseV9 = False




isData = True if ('Run' in samples or 'Data' in samples) else False

lepOpt = 'Ele' if 'Electron' in channel else 'Mu'
ptcut = int(options.ptcut) if options.ptcut.is_integer() else options.ptcut
ptcutHigh = 999999 if -1 == options.ptcutHigh else (   int(options.ptcutHigh) if options.ptcutHigh.is_integer() else options.ptcutHigh  )  


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

'''
def hasMomRecursive(i, pdgid, ch):
    if( ch.GenPart_genPartIdxMother[i] == -1):
        return False
    elif(abs(ch.GenPart_pdgId[ch.GenPart_genPartIdxMother[i]]) == pdgid):
        return True
    else:
        return hasMomRecursive( ch.GenPart_genPartIdxMother[i], pdgid, ch )
'''


def hasMomRecursive(i, pdgid, ch):
    particle = i
    particle_j = particle
    depth = 0
    while(depth <= ch.nGenPart):
        momidx = ch.GenPart_genPartIdxMother[particle_j]
        if not (0 <= momidx < ch.nGenPart):
            return False
        elif momidx == particle or momidx == particle_j:
            return False
        elif(abs(ch.GenPart_pdgId[momidx]) == pdgid):
            return True
        particle_j = momidx
        depth += 1
    return False


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



def getMiniIsoType(miniRelIso, pt):
    # checks what cut the provided miniRelIso isolation passes
    # -1: fail, 0: passes loose, 1: passes tight
    # cuts calculated based on cuts in https://github.com/1LStopBudapest/Helper/blob/master/TreeVarSel.py#L282 for normal iso
    # by weighting with cone area, calculated from different cone sizes (0.3 -> 0.2), weight: 0.42942652
    absIso = miniRelIso*pt
    if pt > 5 and pt < 25:
        if absIso < 2.1471326:
            return 1
        elif absIso < 8.5885305:
            return 0
    elif pt >= 25:
        if miniRelIso < 0.0858853:
            return 1
        elif miniRelIso < 0.3435412:
            return 0
    return -1



# RECO #
def RecoAndRecoFix(threadID, it0, it1, nevtcut, ch_common, ptcut, ptcutHigh, binningSet):

    # Output histos - standalone for each Process
    hfile = ROOT.TFile( 'root_files/RecoFixThreaded_Sample'+sample+'_Thread-'+str(threadID)+'.root', 'RECREATE')
    histos = Histograms.SpawnHistograms(binningSet)

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
            Fill1D(histos['AllRecoPhoton_Pt'], ch.Photon_pt[i])

        # LowPtEle statistics
        if(bUseV9):
            for i in range(ch.nLowPtElectron):
                if(ch.LowPtElectron_pt[i] >= ptcut):
                    Fill1D(histos['LowPtEle_Pt'], ch.LowPtElectron_pt[i])
                    Fill1D(histos['LowPtEle_Eta'], ch.LowPtElectron_eta[i])
                    Fill1D(histos['LowPtEle_Phi'], ch.LowPtElectron_phi[i])
                    Fill1D(histos['LowPtEle_DxyAbs'], abs(ch.LowPtElectron_dxy[i]))
                    Fill1D(histos['LowPtEle_DzAbs'], abs(ch.LowPtElectron_dz[i]))
                    
                    Fill1D(histos['LowPtEleBDTID'], ch.LowPtElectron_ID[i])
                    Fill1D(histos['LowPtEleBDTembeddedID'], ch.LowPtElectron_embeddedID[i])
                    histos['LowPtEleBDT2D'].Fill(ch.LowPtElectron_ID[i], ch.LowPtElectron_embeddedID[i])

                    histos['LowPtEleBDT2D'].Fill(ch.LowPtElectron_ID[i], ch.LowPtElectron_embeddedID[i])

                    histos['LowPtEle_BDT_RelIso_2D'].Fill(ch.LowPtElectron_ID[i], ch.LowPtElectron_miniPFRelIso_all[i])
                    histos['LowPtEle_BDT_AbsIso_2D'].Fill(ch.LowPtElectron_ID[i], ch.LowPtElectron_miniPFRelIso_all[i] * ch.LowPtElectron_pt[i])
                    if(ch.LowPtElectron_pt[i] < 3.5):
                        histos['LowPtEle_BDT_AbsIso_2D_ptmin_3.5'].Fill(ch.LowPtElectron_ID[i], ch.LowPtElectron_miniPFRelIso_all[i] * ch.LowPtElectron_pt[i])
                    elif(ch.LowPtElectron_pt[i] < 10):
                        histos['LowPtEle_BDT_AbsIso_2D_pt3.5_10'].Fill(ch.LowPtElectron_ID[i], ch.LowPtElectron_miniPFRelIso_all[i] * ch.LowPtElectron_pt[i])
                    elif(ch.LowPtElectron_pt[i] < 25):
                        histos['LowPtEle_BDT_AbsIso_2D_pt10_25'].Fill(ch.LowPtElectron_ID[i], ch.LowPtElectron_miniPFRelIso_all[i] * ch.LowPtElectron_pt[i])

                    
                    miniIsoID = getMiniIsoType(ch.LowPtElectron_miniPFRelIso_all[i],ch.LowPtElectron_pt[i] )
                    if( miniIsoID == 0 or miniIsoID == 1):
                        histos['LowPtEle_BDT_RelIso_2D_Loose'].Fill(ch.LowPtElectron_ID[i], ch.LowPtElectron_miniPFRelIso_all[i])
                    if(miniIsoID == 1):
                        histos['LowPtEle_BDT_RelIso_2D_Tight'].Fill(ch.LowPtElectron_ID[i], ch.LowPtElectron_miniPFRelIso_all[i])



        # create isotrack - Photon pairs for extended reco
        aPhotonIsoTrackPairs = []
        aPhotonIsoTrackPairs_CBLoose = []
        aPhotonIsoTrackPairs_CBMedium = []
        aPhotonIsoTrackPairs_CBTight = []
        aPhotonIsoTrackPairs_MVA90 = []
        aPhotonIsoTrackPairs_MVA80 = []
        aPhotonIsoTrackExtrapolPairs = []
        for i in range(ch.nIsoTrack):
            # [Fail,CBLoose,CBMed,CBTight,MVA90,MVA80]
            fTrackPhoDist = [100000,100000,100000,100000,100000,100000,100000]
            fTrackPhoDistExtrapol = 100000
            dist0 = 100000
            for j in range(ch.nPhoton):
                
                gencomp_eta = ch.IsoTrack_eta[i]
                gencomp_phi = ch.IsoTrack_phi[i]
                dist0 = dR(gencomp_eta, ch.Photon_eta[j], gencomp_phi, ch.Photon_phi[j] )

                if( bExtrapolate):
                    gencomp_eta = ch.IsoTrack_eta[i] + ch.IsoTrack_deltaEta[i]
                    gencomp_phi = ch.IsoTrack_phi[i] + ch.IsoTrack_deltaPhi[i]
                    dist0 = dR(gencomp_eta, ch.Photon_eta[j], gencomp_phi, ch.Photon_phi[j] )


                if( dist0 < fTrackPhoDist[0]):
                    fTrackPhoDist[0] = dist0

                for iCut in range(1,4,1):
                    if(ch.Photon_cutBased[j] >= iCut and dist0 < fTrackPhoDist[iCut]):
                        fTrackPhoDist[iCut] = dist0

                if( ch.Photon_mvaID_WP90[j] and dist0 < fTrackPhoDist[4]):
                    fTrackPhoDist[4] = dist0
                if( ch.Photon_mvaID_WP80[j] and dist0 < fTrackPhoDist[5]):
                    fTrackPhoDist[5] = dist0

            # FIXME!!!!!! LIKE THIS, EXTRAPOLATING WILL NOT BE RECOGNIZED IN THE ID, BUT ALSO WON'T GIVE AN ERRORMESSAGE ABOUT IT!!!
            if( dRcut(fTrackPhoDist[0], 0.2) ):
                aPhotonIsoTrackPairs.append( isoTrackPhotonPair(ch, i, j, bExtrapolate) )
            if( dRcut(fTrackPhoDist[1], 0.2) ):
                aPhotonIsoTrackPairs_CBLoose.append( isoTrackPhotonPair(ch, i, j, bExtrapolate) )
            if( dRcut(fTrackPhoDist[2], 0.2) ):
                aPhotonIsoTrackPairs_CBMedium.append( isoTrackPhotonPair(ch, i, j, bExtrapolate) )
            if( dRcut(fTrackPhoDist[3], 0.2) ):
                aPhotonIsoTrackPairs_CBTight.append( isoTrackPhotonPair(ch, i, j, bExtrapolate) )
            if( dRcut(fTrackPhoDist[4], 0.2) ):
                aPhotonIsoTrackPairs_MVA90.append( isoTrackPhotonPair(ch, i, j, bExtrapolate) )
            if( dRcut(fTrackPhoDist[5], 0.2) ):
                aPhotonIsoTrackPairs_MVA80.append( isoTrackPhotonPair(ch, i, j, bExtrapolate) )

            #if( dRcut(fTrackPhoDistExtrapol, 0.2) ):
            #    aPhotonIsoTrackExtrapolPairs.append( isoTrackPhotonPair(ch, i, j, bExtrapolate) )




        ## RECO STUDY ##

        bdtcuts = [0,1,2,2.5,3,4,5]
        if(bSignal):
            for i in range(ch.nGenPart):

                ##  photon pt test
                if(abs(ch.GenPart_pdgId[i]) == 22):
                    Fill1D(histos['AllTruePhoton_Pt'], ch.GenPart_pt[i])


                # electron, pt eta cuts, status 1, with a stop mother somewhere in mother history
                # pt >= 15, |eta| < 1
                # and abs(ch.GenPart_vx[i]) < 0.2 and abs(ch.GenPart_vy[i]) < 0.2
                # hasMomRecursive(i, 1000006, ch) stop, changed to W (24)
                if(  abs(ch.GenPart_pdgId[i]) == 11 and ptcut <= ch.GenPart_pt[i] <= ptcutHigh  and abs(ch.GenPart_eta[i]) <= 2.5 and ch.GenPart_status[i] == 1 and hasMomRecursive(i, 24, ch)):

                    if(bUseVertex):
                        GenPart_dxy = getdxy(i, ch)
                        GenPart_dz = getdz(i, ch)

                    Fill1D(histos['TrueEle_Pt'],ch.GenPart_pt[i])
                    Fill1D(histos['TrueEle_Eta'],ch.GenPart_eta[i])
                    Fill1D(histos['TrueEle_Phi'],ch.GenPart_phi[i])

                    if(bUseVertex):
                        Fill1D(histos['TrueEle_Dxy'], GenPart_dxy)
                        Fill1D(histos['TrueEle_Dz'], GenPart_dz)
                        Fill1D(histos['TrueEle_DxyAbs'], abs(GenPart_dxy))
                        Fill1D(histos['TrueEle_DzAbs'], abs(GenPart_dz))

                        # WATCHME: overflow added to last pt bin here
                        ovpt = ch.GenPart_pt[i]
                        ovdxy = abs(GenPart_dxy)
                        if( ovpt > 100):
                            ovpt = 99
                        if( ovdxy > 15):
                            ovdxy = 14.8
                        histos['TrueElePtDxy2D'].Fill( ovpt, ovdxy )




                    #########################
                    ## SIGNAL RECO METHODS ##
                    #########################



                    # Standard Reco:
                    # check if there is a reco ele / isotrack:
                    bRecoAsStandard = False
                    abRecoAsStandardCutBased = [False, False, False, False]
                    abRecoAsStandardMVA = [False, False, False]
                    eledist = 100000
                    eleIdx = -1
                    idx0 = -1
                    for j in range(ch.nElectron):

                        dist0 = dR(ch.GenPart_eta[i], ch.Electron_eta[j], ch.GenPart_phi[i], ch.Electron_phi[j] ) 
                        if( dist0 < eledist):
                            eledist = dist0
                            idx0 = j
                    if( dRcut(eledist) ):
                        bRecoAsStandard = True
                        eleIdx = idx0
                    # save info:
                    foundStandardElectron = eleIdx # save this for overlap study: if -1, no overlap, if not -1, there was an overlap


                    # Standard Reco with cutbased
                    for iCut in range(4):
                        eledist = 100000
                        eleIdx = -1
                        idx0 = -1
                        for j in range(ch.nElectron):
                            #if( ch.Electron_cutBased[j] >> iCut & 1):
                            #if(ch.Electron_cutBased[j] > iCut):   <-- this one contains isolation

                            if(eleVID( ch.Electron_vidNestedWPBitmap[j], iCut + 1, removedCuts=['pfRelIso03_all'])): #cutbased id: 0:fail, 1:veto, 2:loose, 3:medium, 4:tight

                                dist0 = dR(ch.GenPart_eta[i], ch.Electron_eta[j], ch.GenPart_phi[i], ch.Electron_phi[j] ) 
                                if( dist0 < eledist):
                                    eledist = dist0
                                    idx0 = j
                        if( dRcut(eledist) ):
                            abRecoAsStandardCutBased[iCut] = True



                    # Standard Reco with MVAs
                    # check if there is a reco ele / isotrack:
                    eledist = 100000
                    eleIdx = -1
                    idx0 = -1
                    for j in range(ch.nElectron):
                        if( ch.Electron_mvaFall17V2noIso_WPL[j]):
                            dist0 = dR(ch.GenPart_eta[i], ch.Electron_eta[j], ch.GenPart_phi[i], ch.Electron_phi[j] ) 
                            if( dist0 < eledist):
                                eledist = dist0
                                idx0 = j
                    if( dRcut(eledist) ):
                        abRecoAsStandardMVA[0] = True

                    eledist = 100000
                    eleIdx = -1
                    idx0 = -1
                    for j in range(ch.nElectron):
                        if( ch.Electron_mvaFall17V2noIso_WP90[j]):
                            dist0 = dR(ch.GenPart_eta[i], ch.Electron_eta[j], ch.GenPart_phi[i], ch.Electron_phi[j] ) 
                            if( dist0 < eledist):
                                eledist = dist0
                                idx0 = j
                    if( dRcut(eledist) ):
                        abRecoAsStandardMVA[1] = True

                    eledist = 100000
                    eleIdx = -1
                    idx0 = -1
                    for j in range(ch.nElectron):
                        if( ch.Electron_mvaFall17V2noIso_WP80[j]):
                            dist0 = dR(ch.GenPart_eta[i], ch.Electron_eta[j], ch.GenPart_phi[i], ch.Electron_phi[j] ) 
                            if( dist0 < eledist):
                                eledist = dist0
                                idx0 = j
                    if( dRcut(eledist) ):
                        abRecoAsStandardMVA[2] = True




                    # V9: LowPt Electron search included
                    bRecoAsLowPt = False
                    abRecoAsLowPtBDT = [False] * len(bdtcuts)
                    bdtidx = [-1] * len(bdtcuts)
                    bdtdist = [100000] * len(bdtcuts)
                    if(bUseV9):
                        eledist = 100000
                        eleIdx = -1
                        idx0 = -1
                        for j in range(ch.nLowPtElectron):

                            dist0 = dR(ch.GenPart_eta[i], ch.LowPtElectron_eta[j], ch.GenPart_phi[i], ch.LowPtElectron_phi[j] ) 
                            if( dist0 < eledist):
                                eledist = dist0
                                idx0 = j
                            for k in range(len(bdtcuts)):
                                if(bdtdist[k] > dist0 and ch.LowPtElectron_ID[j] > bdtcuts[k]):
                                    bdtdist[k] = dist0
                                    bdtidx[k] = j

                        if( dRcut(eledist) ):
                            bRecoAsLowPt = True
                            eleIdx = idx0 
                        # save info:
                        foundLowPtElectron = eleIdx  # save this for later overlap studies with extended

                        for k in range(len(bdtcuts)):
                            if(dRcut(bdtdist[k])):
                                abRecoAsLowPtBDT[k] = True



                    # IPP reco 2.0
                    # match to isotrack-photon pairs -> then the isotrack is the detected ele


                    bRecoAsIPP = False
                    abRecoAsIPPCutBased = [False,False,False]
                    abRecoAsIPPMVA = [False,False]
                    fPairDist = 100000
                    iPairIdx = -1
                    for j in range(len(aPhotonIsoTrackPairs)):
                        dist0 = dR(ch.GenPart_eta[i], aPhotonIsoTrackPairs[j].GetEta(), ch.GenPart_phi[i], aPhotonIsoTrackPairs[j].GetPhi() ) 
                        if( dist0 < fPairDist):
                            fPairDist = dist0
                            iPairIdx = j
                    if( dRcut(fPairDist, 0.2) ):
                        bRecoAsIPP = True

                    fPairDist = 100000
                    iPairIdx = -1
                    for j in range(len(aPhotonIsoTrackPairs_CBLoose)):
                        dist0 = dR(ch.GenPart_eta[i], aPhotonIsoTrackPairs_CBLoose[j].GetEta(), ch.GenPart_phi[i], aPhotonIsoTrackPairs_CBLoose[j].GetPhi() ) 
                        if( dist0 < fPairDist):
                            fPairDist = dist0
                            iPairIdx = j
                    if( dRcut(fPairDist, 0.2) ):
                        abRecoAsIPPCutBased[0] = True

                    fPairDist = 100000
                    iPairIdx = -1
                    for j in range(len(aPhotonIsoTrackPairs_CBMedium)):
                        dist0 = dR(ch.GenPart_eta[i], aPhotonIsoTrackPairs_CBMedium[j].GetEta(), ch.GenPart_phi[i], aPhotonIsoTrackPairs_CBMedium[j].GetPhi() ) 
                        if( dist0 < fPairDist):
                            fPairDist = dist0
                            iPairIdx = j
                    if( dRcut(fPairDist, 0.2) ):
                        abRecoAsIPPCutBased[1] = True


                    fPairDist = 100000
                    iPairIdx = -1
                    for j in range(len(aPhotonIsoTrackPairs_CBTight)):
                        dist0 = dR(ch.GenPart_eta[i], aPhotonIsoTrackPairs_CBTight[j].GetEta(), ch.GenPart_phi[i], aPhotonIsoTrackPairs_CBTight[j].GetPhi() ) 
                        if( dist0 < fPairDist):
                            fPairDist = dist0
                            iPairIdx = j
                    if( dRcut(fPairDist, 0.2) ):
                        abRecoAsIPPCutBased[2] = True


                    fPairDist = 100000
                    iPairIdx = -1
                    for j in range(len(aPhotonIsoTrackPairs_MVA90)):
                        dist0 = dR(ch.GenPart_eta[i], aPhotonIsoTrackPairs_MVA90[j].GetEta(), ch.GenPart_phi[i], aPhotonIsoTrackPairs_MVA90[j].GetPhi() ) 
                        if( dist0 < fPairDist):
                            fPairDist = dist0
                            iPairIdx = j
                    if( dRcut(fPairDist, 0.2) ):
                        abRecoAsIPPMVA[0] = True

                    fPairDist = 100000
                    iPairIdx = -1
                    for j in range(len(aPhotonIsoTrackPairs_MVA80)):
                        dist0 = dR(ch.GenPart_eta[i], aPhotonIsoTrackPairs_MVA80[j].GetEta(), ch.GenPart_phi[i], aPhotonIsoTrackPairs_MVA80[j].GetPhi() ) 
                        if( dist0 < fPairDist):
                            fPairDist = dist0
                            iPairIdx = j
                    if( dRcut(fPairDist, 0.2) ):
                        abRecoAsIPPMVA[1] = True

                    



                    ############################
                    ## FILL SIGNAL HISTOGRAMS ##
                    ############################


                    if( bRecoAsStandard):
                        # Fill just standard
                        Fill1D(histos['RecoStandardEle_Pt'], ch.GenPart_pt[i])
                        Fill1D(histos['RecoStandardEle_Eta'], ch.GenPart_eta[i])
                        Fill1D(histos['RecoStandardEle_Phi'], ch.GenPart_phi[i])

                        if(bUseVertex):
                            Fill1D(histos['RecoStandardEle_Dxy'], GenPart_dxy)
                            Fill1D(histos['RecoStandardEle_Dz'], GenPart_dz)
                            Fill1D(histos['RecoStandardEle_DxyAbs'], abs(GenPart_dxy))
                            Fill1D(histos['RecoStandardEle_DzAbs'], abs(GenPart_dz))

                        # Fill standard || LowPt -> moved to end


                        # Fill extended (standard || track/pho pair) (with/without extrapol)
                        Fill1D(histos['RecoExtendedEle_Pt'], ch.GenPart_pt[i])
                        Fill1D(histos['RecoExtendedEle_Eta'], ch.GenPart_eta[i])
                        Fill1D(histos['RecoExtendedEle_Phi'], ch.GenPart_phi[i])

                        Fill1D(histos['RecoExtended+ExtrapolEle_Pt'], ch.GenPart_pt[i])
                        Fill1D(histos['RecoExtended+ExtrapolEle_Eta'], ch.GenPart_eta[i])
                        Fill1D(histos['RecoExtended+ExtrapolEle_Phi'], ch.GenPart_phi[i])

                        if(bUseVertex):
                            Fill1D(histos['RecoExtendedEle_Dxy'], GenPart_dxy)
                            Fill1D(histos['RecoExtendedEle_Dz'], GenPart_dz)
                            Fill1D(histos['RecoExtendedEle_DxyAbs'], abs(GenPart_dxy))
                            Fill1D(histos['RecoExtendedEle_DzAbs'], abs(GenPart_dz))

                            Fill1D(histos['RecoExtended+ExtrapolEle_Dxy'], GenPart_dxy)
                            Fill1D(histos['RecoExtended+ExtrapolEle_Dz'], GenPart_dz)
                            Fill1D(histos['RecoExtended+ExtrapolEle_DxyAbs'], abs(GenPart_dxy))
                            Fill1D(histos['RecoExtended+ExtrapolEle_DzAbs'], abs(GenPart_dz))


                            if( ovpt > 100):
                                ovpt = 99
                            if( ovdxy > 15):
                                ovdxy = 14.8
                            histos['ExtendedRecoExtrapolPtDxy2D'].Fill( ovpt, ovdxy )
                            histos['ExtendedRecoPtDxy2D'].Fill( ovpt, ovdxy )




                        # Fill Extended (standard || track/pho pair) || lowpt
                        Fill1D(histos['RecoExtended||LowPtEle_Pt'], ch.GenPart_pt[i])
                        Fill1D(histos['RecoExtended||LowPtEle_Eta'], ch.GenPart_eta[i])
                        Fill1D(histos['RecoExtended||LowPtEle_Phi'], ch.GenPart_phi[i])

                        # extended+lowpt+extrapol doesnt exist yet
                        #Fill1D(histos['ExtendedRecoWithLowPtExtrapolEle_Pt'], ch.GenPart_pt[i])
                        #Fill1D(histos['ExtendedRecoWithLowPtExtrapolEle_Eta'], ch.GenPart_eta[i])
                        #Fill1D(histos['ExtendedRecoWithLowPtExtrapolEle_Phi'], ch.GenPart_phi[i])
                        if(bUseVertex):
                            Fill1D(histos['RecoExtended||LowPtEle_Dxy'], GenPart_dxy)
                            Fill1D(histos['RecoExtended||LowPtEle_Dz'], GenPart_dz)
                            Fill1D(histos['RecoExtended||LowPtEle_DxyAbs'], abs(GenPart_dxy))
                            Fill1D(histos['RecoExtended||LowPtEle_DzAbs'], abs(GenPart_dz))


                    if(abRecoAsStandardCutBased[0]):
                        Fill1D(histos['RecoStandardEle_CBVeto_Pt'], ch.GenPart_pt[i])
                        if(bUseVertex):
                            Fill1D(histos['RecoStandardEle_CBVeto_DxyAbs'], abs(GenPart_dxy))
                    if(abRecoAsStandardCutBased[1]):
                        Fill1D(histos['RecoStandardEle_CBLoose_Pt'], ch.GenPart_pt[i])
                        if(bUseVertex):
                            Fill1D(histos['RecoStandardEle_CBLoose_DxyAbs'], abs(GenPart_dxy))
                    if( abRecoAsStandardCutBased[2] ):
                        Fill1D(histos['RecoStandardEle_CBMedium_Pt'], ch.GenPart_pt[i])
                        if(bUseVertex):
                            Fill1D(histos['RecoStandardEle_CBMedium_DxyAbs'], abs(GenPart_dxy))
                    if( abRecoAsStandardCutBased[3] ):
                        Fill1D(histos['RecoStandardEle_CBTight_Pt'], ch.GenPart_pt[i])
                        if(bUseVertex):
                            Fill1D(histos['RecoStandardEle_CBTight_DxyAbs'], abs(GenPart_dxy))


                    if(abRecoAsStandardMVA[0]):
                        Fill1D(histos['RecoStandardEle_mvaNoIsoWPL_Pt'], ch.GenPart_pt[i])
                        Fill1D(histos['RecoStandardEle_mvaNoIsoWPL_Eta'], ch.GenPart_eta[i])
                        Fill1D(histos['RecoStandardEle_mvaNoIsoWPL_Phi'], ch.GenPart_phi[i])

                        if(bUseVertex):
                            Fill1D(histos['RecoStandardEle_mvaNoIsoWPL_DxyAbs'], abs(GenPart_dxy))
                            Fill1D(histos['RecoStandardEle_mvaNoIsoWPL_DzAbs'], abs(GenPart_dz))

                    if(abRecoAsStandardMVA[1]):
                        Fill1D(histos['RecoStandardEle_mvaNoIsoWP90_Pt'], ch.GenPart_pt[i])
                        Fill1D(histos['RecoStandardEle_mvaNoIsoWP90_Eta'], ch.GenPart_eta[i])
                        Fill1D(histos['RecoStandardEle_mvaNoIsoWP90_Phi'], ch.GenPart_phi[i])

                        if(bUseVertex):
                            Fill1D(histos['RecoStandardEle_mvaNoIsoWP90_DxyAbs'], abs(GenPart_dxy))
                            Fill1D(histos['RecoStandardEle_mvaNoIsoWP90_DzAbs'], abs(GenPart_dz))

                    if(abRecoAsStandardMVA[2]):
                        Fill1D(histos['RecoStandardEle_mvaNoIsoWP80_Pt'], ch.GenPart_pt[i])
                        Fill1D(histos['RecoStandardEle_mvaNoIsoWP80_Eta'], ch.GenPart_eta[i])
                        Fill1D(histos['RecoStandardEle_mvaNoIsoWP80_Phi'], ch.GenPart_phi[i])

                        if(bUseVertex):
                            Fill1D(histos['RecoStandardEle_mvaNoIsoWP80_DxyAbs'], abs(GenPart_dxy))
                            Fill1D(histos['RecoStandardEle_mvaNoIsoWP80_DzAbs'], abs(GenPart_dz))




                    for k in range(len(bdtcuts)):
                        if(abRecoAsLowPtBDT[k]):
                            Fill1D(histos['RecoLowPtEle_Pt_BDT'+str(bdtcuts[k])], ch.GenPart_pt[i])
                            if(bUseVertex):
                                Fill1D(histos['RecoLowPtEle_DxyAbs_BDT'+str(bdtcuts[k])], abs(GenPart_dxy))


                    if( bRecoAsLowPt):
                        # record lowptElectron reco efficiency, on its own 
                        Fill1D(histos['RecoLowPtEle_Pt'], ch.GenPart_pt[i])
                        Fill1D(histos['RecoLowPtEle_Eta'], ch.GenPart_eta[i])
                        Fill1D(histos['RecoLowPtEle_Phi'], ch.GenPart_phi[i])
                        if(bUseVertex):
                            Fill1D(histos['RecoLowPtEle_DxyAbs'], abs(GenPart_dxy))
                            Fill1D(histos['RecoLowPtEle_DzAbs'], abs(GenPart_dz))

                        Fill1D(histos['RecoLowPtEleBDTID'], ch.LowPtElectron_ID[eleIdx])
                        Fill1D(histos['RecoLowPtEleBDTembeddedID'], ch.LowPtElectron_embeddedID[eleIdx])
                        histos['RecoLowPtEleBDT2D'].Fill(ch.LowPtElectron_ID[eleIdx], ch.LowPtElectron_embeddedID[eleIdx])


                        histos['LowPtEle_Matched_BDT_RelIso_2D'].Fill(ch.LowPtElectron_ID[eleIdx], ch.LowPtElectron_miniPFRelIso_all[eleIdx])
                        histos['LowPtEle_Matched_BDT_AbsIso_2D'].Fill(ch.LowPtElectron_ID[eleIdx], ch.LowPtElectron_miniPFRelIso_all[eleIdx] * ch.LowPtElectron_pt[eleIdx])
                        if(ch.LowPtElectron_pt[eleIdx] < 3.5):
                            histos['LowPtEle_Matched_BDT_AbsIso_2D_pt0_3.5'].Fill(ch.LowPtElectron_ID[eleIdx], ch.LowPtElectron_miniPFRelIso_all[eleIdx] * ch.LowPtElectron_pt[eleIdx])
                        elif(ch.LowPtElectron_pt[eleIdx] < 10):
                            histos['LowPtEle_Matched_BDT_AbsIso_2D_pt3.5_10'].Fill(ch.LowPtElectron_ID[eleIdx], ch.LowPtElectron_miniPFRelIso_all[eleIdx] * ch.LowPtElectron_pt[eleIdx])
                        elif(ch.LowPtElectron_pt[eleIdx] < 25):
                            histos['LowPtEle_Matched_BDT_AbsIso_2D_pt10_25'].Fill(ch.LowPtElectron_ID[eleIdx], ch.LowPtElectron_miniPFRelIso_all[eleIdx] * ch.LowPtElectron_pt[eleIdx])

                        # FIXME! RENDES SELECTION!!!
                        miniIsoID = getMiniIsoType(ch.LowPtElectron_miniPFRelIso_all[eleIdx],ch.LowPtElectron_pt[eleIdx] )
                        if( miniIsoID == 0 or miniIsoID == 1):
                            histos['LowPtEle_Matched_BDT_RelIso_2D_Loose'].Fill(ch.LowPtElectron_ID[eleIdx], ch.LowPtElectron_miniPFRelIso_all[eleIdx])
                        if(miniIsoID == 1):
                            histos['LowPtEle_Matched_BDT_RelIso_2D_Tight'].Fill(ch.LowPtElectron_ID[eleIdx], ch.LowPtElectron_miniPFRelIso_all[eleIdx])


                        # also record to IPP||LowPt
                        Fill1D(histos['RecoIPP||LowPtEle_Pt'], ch.GenPart_pt[i])
                        Fill1D(histos['RecoIPP||LowPtEle_Eta'], ch.GenPart_eta[i])
                        Fill1D(histos['RecoIPP||LowPtEle_Phi'], ch.GenPart_phi[i])
                        if(bUseVertex):
                            Fill1D(histos['RecoIPP||LowPtEle_DxyAbs'], abs(GenPart_dxy))
                            Fill1D(histos['RecoIPP||LowPtEle_DzAbs'], abs(GenPart_dz))


                        # if there was no overlap, add it to the Standard || LowPt and Extended || LowPt sets
                        if(foundStandardElectron == -1):

                            # Fill Standard || LowPt -> moved to end

                            # also add it to Extended || LowPt. We recorded foundLowPtElectron, we don't add Extended into this later if already found here.
                            Fill1D(histos['RecoExtended||LowPtEle_Pt'], ch.GenPart_pt[i])
                            Fill1D(histos['RecoExtended||LowPtEle_Eta'], ch.GenPart_eta[i])
                            Fill1D(histos['RecoExtended||LowPtEle_Phi'], ch.GenPart_phi[i])

                            # extended+lowpt+extrapol doesnt exist yet
                            #Fill1D(histos['RecoExtended+ExtrapolEle_Pt'], ch.GenPart_pt[i])
                            #Fill1D(histos['RecoExtended+ExtrapolEle_Eta'], ch.GenPart_eta[i])
                            #Fill1D(histos['RecoExtended+ExtrapolEle_Phi'], ch.GenPart_phi[i])

                        
                            if(bUseVertex):
                                Fill1D(histos['RecoExtended||LowPtEle_Dxy'], GenPart_dxy)
                                Fill1D(histos['RecoExtended||LowPtEle_Dz'], GenPart_dz)
                                Fill1D(histos['RecoExtended||LowPtEle_DxyAbs'], abs(GenPart_dxy))
                                Fill1D(histos['RecoExtended||LowPtEle_DzAbs'], abs(GenPart_dz))

                                #Fill1D(histos['RecoExtended+ExtrapolEleVtxx'], ch.GenPart_vx[i])
                                #Fill1D(histos['RecoExtended+ExtrapolEleVtxy'], ch.GenPart_vy[i])
                                #Fill1D(histos['RecoExtended+ExtrapolEleVtxz'], ch.GenPart_vz[i])
                                #Fill1D(histos['RecoExtended+ExtrapolEle_Dxy'], GenPart_dxy)
                                #Fill1D(histos['RecoExtended+ExtrapolEle_Dz'], GenPart_dz)
                                #Fill1D(histos['RecoExtended+ExtrapolEle_DxyAbs'], abs(GenPart_dxy))
                                #Fill1D(histos['RecoExtended+ExtrapolEle_DzAbs'], abs(GenPart_dz))

                                #if( ovpt > 100):
                                #    ovpt = 99
                                #if( ovdxy > 15):
                                #    ovdxy = 14.8
                                #histos['ExtendedRecoExtrapolPtDxy2D'].Fill( ovpt, ovdxy )
                                #histos['ExtendedRecoPtDxy2D'].Fill( ovpt, ovdxy )
                        else:
                            # Also found as standard electron: add it to Overlaps: 
                            # Fill Standard Overlap (Standrd && LowPt):
                            Fill1D(histos['OverlapStandard&&LowPtEle_Pt'], ch.GenPart_pt[i])
                            Fill1D(histos['OverlapStandard&&LowPtEle_Eta'], ch.GenPart_eta[i])
                            Fill1D(histos['OverlapStandard&&LowPtEle_Phi'], ch.GenPart_phi[i])
                            if(bUseVertex):
                                Fill1D(histos['OverlapStandard&&LowPtEle_DxyAbs'], abs(GenPart_dxy))
                                Fill1D(histos['OverlapStandard&&LowPtEle_DzAbs'], abs(GenPart_dz))

                            # Fill Extended (IPP || standard) overlap with lowPt:
                            Fill1D(histos['OverlapExtended&&LowPtEle_Pt'], ch.GenPart_pt[i])
                            Fill1D(histos['OverlapExtended&&LowPtEle_Eta'], ch.GenPart_eta[i])
                            Fill1D(histos['OverlapExtended&&LowPtEle_Phi'], ch.GenPart_phi[i])
                            if(bUseVertex):
                                Fill1D(histos['OverlapExtended&&LowPtEle_DxyAbs'], abs(GenPart_dxy))
                                Fill1D(histos['OverlapExtended&&LowPtEle_DzAbs'], abs(GenPart_dz))



                    if( bRecoAsIPP):
                        Fill1D(histos['RecoIPP_Pt'], ch.GenPart_pt[i])
                        Fill1D(histos['RecoIPP_Eta'], ch.GenPart_eta[i])
                        Fill1D(histos['RecoIPP_Phi'], ch.GenPart_phi[i])
                        if(bUseVertex):
                            Fill1D(histos['RecoIPP_DxyAbs'], abs(GenPart_dxy))
                            Fill1D(histos['RecoIPP_DzAbs'], abs(GenPart_dz))


                        if(foundStandardElectron == -1):
                            # if not found in Standard, fill Extended
                            Fill1D(histos['RecoExtendedEle_Pt'], ch.GenPart_pt[i])
                            Fill1D(histos['RecoExtendedEle_Eta'], ch.GenPart_eta[i])
                            Fill1D(histos['RecoExtendedEle_Phi'], ch.GenPart_phi[i])

                            if(bUseVertex):
                                Fill1D(histos['RecoExtendedEle_Dxy'], GenPart_dxy)
                                Fill1D(histos['RecoExtendedEle_Dz'], GenPart_dz)
                                Fill1D(histos['RecoExtendedEle_DxyAbs'], abs(GenPart_dxy))
                                Fill1D(histos['RecoExtendedEle_DzAbs'], abs(GenPart_dz))

                                if( ovpt > 100):
                                    ovpt = 99
                                if( ovdxy > 15):
                                    ovdxy = 14.8
                                histos['ExtendedRecoPtDxy2D'].Fill( ovpt, ovdxy )
            

                            # if ALSO not already found as lowPT, fill Extended || LowPt
                            if(foundLowPtElectron == -1):
                                Fill1D(histos['RecoExtended||LowPtEle_Pt'], ch.GenPart_pt[i])
                                Fill1D(histos['RecoExtended||LowPtEle_Eta'], ch.GenPart_eta[i])
                                Fill1D(histos['RecoExtended||LowPtEle_Phi'], ch.GenPart_phi[i])

                                if(bUseVertex):
                                    Fill1D(histos['RecoExtended||LowPtEle_Dxy'], GenPart_dxy)
                                    Fill1D(histos['RecoExtended||LowPtEle_Dz'], GenPart_dz)
                                    Fill1D(histos['RecoExtended||LowPtEle_DxyAbs'], abs(GenPart_dxy))
                                    Fill1D(histos['RecoExtended||LowPtEle_DzAbs'], abs(GenPart_dz))
                            else:
                                # also found as lowPt and as track/pho pair: fill extended overlap
                                # Fill Extended (track/pho || standard) overlap with lowPt:
                                Fill1D(histos['OverlapExtended&&LowPtEle_Pt'], ch.GenPart_pt[i])
                                Fill1D(histos['OverlapExtended&&LowPtEle_Eta'], ch.GenPart_eta[i])
                                Fill1D(histos['OverlapExtended&&LowPtEle_Phi'], ch.GenPart_phi[i])
                                if(bUseVertex):
                                    Fill1D(histos['OverlapExtended&&LowPtEle_DxyAbs'], abs(GenPart_dxy))
                                    Fill1D(histos['OverlapExtended&&LowPtEle_DzAbs'], abs(GenPart_dz))

                            '''
                            if(bExtrapolate):
                                fPairDist = 100000
                                iPairIdx = -1
                                for j in range(len(aPhotonIsoTrackExtrapolPairs)):
                                    
                                    dist0 = dR(ch.GenPart_eta[i], aPhotonIsoTrackExtrapolPairs[j].GetEta(), ch.GenPart_phi[i], aPhotonIsoTrackExtrapolPairs[j].GetPhi()) 

                                    if( dist0 < fPairDist):
                                        fPairDist = dist0
                                        iPairIdx = j

                                
                                if( dRcut(fPairDist, 0.2) ):
                                    Fill1D(histos['RecoExtended+ExtrapolEle_Pt'], ch.GenPart_pt[i])
                                    Fill1D(histos['RecoExtended+ExtrapolEle_Eta'], ch.GenPart_eta[i])
                                    Fill1D(histos['RecoExtended+ExtrapolEle_Phi'], ch.GenPart_phi[i])

                                    if(bUseVertex):
                                        Fill1D(histos['RecoExtended+ExtrapolEle_Dxy'], GenPart_dxy)
                                        Fill1D(histos['RecoExtended+ExtrapolEle_Dz'], GenPart_dz)
                                        Fill1D(histos['RecoExtended+ExtrapolEle_DxyAbs'], abs(GenPart_dxy))
                                        Fill1D(histos['RecoExtended+ExtrapolEle_DzAbs'], abs(GenPart_dz))

                                        if( ovpt > 100):
                                            ovpt = 99
                                        if( ovdxy > 15):
                                            ovdxy = 14.8
                                        histos['ExtendedRecoExtrapolPtDxy2D'].Fill( ovpt, ovdxy )
                            '''
                        else:
                            # can be reconstructed as Standard and as IPP as well, fill Standard&&Extended overlap:
                            Fill1D(histos['OverlapIPP&&StandardEle_Pt'], ch.GenPart_pt[i])
                            Fill1D(histos['OverlapIPP&&StandardEle_Eta'], ch.GenPart_eta[i])
                            Fill1D(histos['OverlapIPP&&StandardEle_Phi'], ch.GenPart_phi[i])
                            if(bUseVertex):
                                Fill1D(histos['OverlapIPP&&StandardEle_DxyAbs'], abs(GenPart_dxy))
                                Fill1D(histos['OverlapIPP&&StandardEle_DzAbs'], abs(GenPart_dz))


                        if(foundLowPtElectron == -1):
                            # not found as LP, fill IPP||LP
                            Fill1D(histos['RecoIPP||LowPtEle_Pt'], ch.GenPart_pt[i])
                            Fill1D(histos['RecoIPP||LowPtEle_Eta'], ch.GenPart_eta[i])
                            Fill1D(histos['RecoIPP||LowPtEle_Phi'], ch.GenPart_phi[i])
                            if(bUseVertex):
                                Fill1D(histos['RecoIPP||LowPtEle_DxyAbs'], abs(GenPart_dxy))
                                Fill1D(histos['RecoIPP||LowPtEle_DzAbs'], abs(GenPart_dz))
                        else:
                            # reco as both LP and IPP, fill overlap
                            Fill1D(histos['OverlapIPP&&LowPtEle_Pt'], ch.GenPart_pt[i])
                            Fill1D(histos['OverlapIPP&&LowPtEle_Eta'], ch.GenPart_eta[i])
                            Fill1D(histos['OverlapIPP&&LowPtEle_Phi'], ch.GenPart_phi[i])
                            if(bUseVertex):
                                Fill1D(histos['OverlapIPP&&LowPtEle_DxyAbs'], abs(GenPart_dxy))
                                Fill1D(histos['OverlapIPP&&LowPtEle_DzAbs'], abs(GenPart_dz))

                    if(abRecoAsIPPCutBased[0]):
                        Fill1D(histos['RecoIPP_CBLoose_Pt'], ch.GenPart_pt[i])
                        if(bUseVertex):
                            Fill1D(histos['RecoIPP_CBLoose_DxyAbs'], abs(GenPart_dxy))
                    if(abRecoAsIPPCutBased[1]):
                        Fill1D(histos['RecoIPP_CBMedium_Pt'], ch.GenPart_pt[i])
                        if(bUseVertex):
                            Fill1D(histos['RecoIPP_CBMedium_DxyAbs'], abs(GenPart_dxy))
                    if(abRecoAsIPPCutBased[2]):
                        Fill1D(histos['RecoIPP_CBTight_Pt'], ch.GenPart_pt[i])
                        if(bUseVertex):
                            Fill1D(histos['RecoIPP_CBTight_DxyAbs'], abs(GenPart_dxy))
                    if(abRecoAsIPPMVA[0]):
                        Fill1D(histos['RecoIPP_MVA90_Pt'], ch.GenPart_pt[i])
                        if(bUseVertex):
                            Fill1D(histos['RecoIPP_MVA90_DxyAbs'], abs(GenPart_dxy))
                    if(abRecoAsIPPMVA[1]):
                        Fill1D(histos['RecoIPP_MVA80_Pt'], ch.GenPart_pt[i])
                        if(bUseVertex):
                            Fill1D(histos['RecoIPP_MVA80_DxyAbs'], abs(GenPart_dxy))


                    for k in range(len(bdtcuts)):
                        if(abRecoAsLowPtBDT[k] or bRecoAsIPP):
                            Fill1D(histos['RecoIPP||LowPtEle_Pt_BDT'+str(bdtcuts[k])], ch.GenPart_pt[i])
                            if(bUseVertex):
                                Fill1D(histos['RecoIPP||LowPtEle_DxyAbs_BDT'+str(bdtcuts[k])], abs(GenPart_dxy))

                        if(abRecoAsLowPtBDT[k] or abRecoAsIPPCutBased[0]):
                            Fill1D(histos['RecoIPP_CBLoose||LowPtEle_Pt_BDT'+str(bdtcuts[k])], ch.GenPart_pt[i])
                            if(bUseVertex):
                                Fill1D(histos['RecoIPP_CBLoose||LowPtEle_DxyAbs_BDT'+str(bdtcuts[k])], abs(GenPart_dxy))

                        if(abRecoAsLowPtBDT[k] or abRecoAsIPPCutBased[1]):
                            Fill1D(histos['RecoIPP_CBMedium||LowPtEle_Pt_BDT'+str(bdtcuts[k])], ch.GenPart_pt[i])
                            if(bUseVertex):
                                Fill1D(histos['RecoIPP_CBMedium||LowPtEle_DxyAbs_BDT'+str(bdtcuts[k])], abs(GenPart_dxy))

                        if(abRecoAsLowPtBDT[k] or abRecoAsIPPCutBased[2]):
                            Fill1D(histos['RecoIPP_CBTight||LowPtEle_Pt_BDT'+str(bdtcuts[k])], ch.GenPart_pt[i])
                            if(bUseVertex):
                                Fill1D(histos['RecoIPP_CBTight||LowPtEle_DxyAbs_BDT'+str(bdtcuts[k])], abs(GenPart_dxy))

                        if(abRecoAsLowPtBDT[k] or abRecoAsIPPMVA[0]):
                            Fill1D(histos['RecoIPP_MVA90||LowPtEle_Pt_BDT'+str(bdtcuts[k])], ch.GenPart_pt[i])
                            if(bUseVertex):
                                Fill1D(histos['RecoIPP_MVA90||LowPtEle_DxyAbs_BDT'+str(bdtcuts[k])], abs(GenPart_dxy))

                        if(abRecoAsLowPtBDT[k] or abRecoAsIPPMVA[1]):
                            Fill1D(histos['RecoIPP_MVA80||LowPtEle_Pt_BDT'+str(bdtcuts[k])], ch.GenPart_pt[i])
                            if(bUseVertex):
                                Fill1D(histos['RecoIPP_MVA80||LowPtEle_DxyAbs_BDT'+str(bdtcuts[k])], abs(GenPart_dxy))


                    if(abRecoAsStandardCutBased[0] or bRecoAsIPP):
                        Fill1D(histos['RecoExtendedEle_CBVeto_Pt'], ch.GenPart_pt[i])
                        if(bUseVertex):
                            Fill1D(histos['RecoExtendedEle_CBVeto_DxyAbs'], abs(GenPart_dxy))
                    if(abRecoAsStandardCutBased[1] or bRecoAsIPP):
                        Fill1D(histos['RecoExtendedEle_CBLoose_Pt'], ch.GenPart_pt[i])
                        if(bUseVertex):
                            Fill1D(histos['RecoExtendedEle_CBLoose_DxyAbs'], abs(GenPart_dxy))
                    if( abRecoAsStandardCutBased[2] or bRecoAsIPP ):
                        Fill1D(histos['RecoExtendedEle_CBMedium_Pt'], ch.GenPart_pt[i])
                        if(bUseVertex):
                            Fill1D(histos['RecoExtendedEle_CBMedium_DxyAbs'], abs(GenPart_dxy))
                    if( abRecoAsStandardCutBased[3] or bRecoAsIPP ):
                        Fill1D(histos['RecoExtendedEle_CBTight_Pt'], ch.GenPart_pt[i])
                        if(bUseVertex):
                            Fill1D(histos['RecoExtendedEle_CBTight_DxyAbs'], abs(GenPart_dxy))


                    if(abRecoAsStandardCutBased[0] or abRecoAsIPPCutBased[0]):
                        Fill1D(histos['RecoExtendedEle_CBVeto_PhoCBLoose_Pt'], ch.GenPart_pt[i])
                        if(bUseVertex):
                            Fill1D(histos['RecoExtendedEle_CBVeto_PhoCBLoose_DxyAbs'], abs(GenPart_dxy))
                    if(abRecoAsStandardCutBased[1] or abRecoAsIPPCutBased[0]):
                        Fill1D(histos['RecoExtendedEle_CBLoose_PhoCBLoose_Pt'], ch.GenPart_pt[i])
                        if(bUseVertex):
                            Fill1D(histos['RecoExtendedEle_CBLoose_PhoCBLoose_DxyAbs'], abs(GenPart_dxy))
                    if( abRecoAsStandardCutBased[2] or abRecoAsIPPCutBased[0] ):
                        Fill1D(histos['RecoExtendedEle_CBMedium_PhoCBLoose_Pt'], ch.GenPart_pt[i])
                        if(bUseVertex):
                            Fill1D(histos['RecoExtendedEle_CBMedium_PhoCBLoose_DxyAbs'], abs(GenPart_dxy))
                    if( abRecoAsStandardCutBased[3] or abRecoAsIPPCutBased[0] ):
                        Fill1D(histos['RecoExtendedEle_CBTight_PhoCBLoose_Pt'], ch.GenPart_pt[i])
                        if(bUseVertex):
                            Fill1D(histos['RecoExtendedEle_CBTight_PhoCBLoose_DxyAbs'], abs(GenPart_dxy))

                    if(abRecoAsStandardCutBased[0] or abRecoAsIPPCutBased[1]):
                        Fill1D(histos['RecoExtendedEle_CBVeto_PhoCBMedium_Pt'], ch.GenPart_pt[i])
                        if(bUseVertex):
                            Fill1D(histos['RecoExtendedEle_CBVeto_PhoCBMedium_DxyAbs'], abs(GenPart_dxy))
                    if(abRecoAsStandardCutBased[1] or abRecoAsIPPCutBased[1]):
                        Fill1D(histos['RecoExtendedEle_CBLoose_PhoCBMedium_Pt'], ch.GenPart_pt[i])
                        if(bUseVertex):
                            Fill1D(histos['RecoExtendedEle_CBLoose_PhoCBMedium_DxyAbs'], abs(GenPart_dxy))
                    if( abRecoAsStandardCutBased[2] or abRecoAsIPPCutBased[1] ):
                        Fill1D(histos['RecoExtendedEle_CBMedium_PhoCBMedium_Pt'], ch.GenPart_pt[i])
                        if(bUseVertex):
                            Fill1D(histos['RecoExtendedEle_CBMedium_PhoCBMedium_DxyAbs'], abs(GenPart_dxy))
                    if( abRecoAsStandardCutBased[3] or abRecoAsIPPCutBased[1] ):
                        Fill1D(histos['RecoExtendedEle_CBTight_PhoCBMedium_Pt'], ch.GenPart_pt[i])
                        if(bUseVertex):
                            Fill1D(histos['RecoExtendedEle_CBTight_PhoCBMedium_DxyAbs'], abs(GenPart_dxy))

                    if(abRecoAsStandardCutBased[0] or abRecoAsIPPCutBased[2]):
                        Fill1D(histos['RecoExtendedEle_CBVeto_PhoCBTight_Pt'], ch.GenPart_pt[i])
                        if(bUseVertex):
                            Fill1D(histos['RecoExtendedEle_CBVeto_PhoCBTight_DxyAbs'], abs(GenPart_dxy))
                    if(abRecoAsStandardCutBased[1] or abRecoAsIPPCutBased[2]):
                        Fill1D(histos['RecoExtendedEle_CBLoose_PhoCBTight_Pt'], ch.GenPart_pt[i])
                        if(bUseVertex):
                            Fill1D(histos['RecoExtendedEle_CBLoose_PhoCBTight_DxyAbs'], abs(GenPart_dxy))
                    if( abRecoAsStandardCutBased[2] or abRecoAsIPPCutBased[2] ):
                        Fill1D(histos['RecoExtendedEle_CBMedium_PhoCBTight_Pt'], ch.GenPart_pt[i])
                        if(bUseVertex):
                            Fill1D(histos['RecoExtendedEle_CBMedium_PhoCBTight_DxyAbs'], abs(GenPart_dxy))
                    if( abRecoAsStandardCutBased[3] or abRecoAsIPPCutBased[2] ):
                        Fill1D(histos['RecoExtendedEle_CBTight_PhoCBTight_Pt'], ch.GenPart_pt[i])
                        if(bUseVertex):
                            Fill1D(histos['RecoExtendedEle_CBTight_PhoCBTight_DxyAbs'], abs(GenPart_dxy))


                    if(abRecoAsStandardCutBased[0] or abRecoAsIPPMVA[0]):
                        Fill1D(histos['RecoExtendedEle_CBVeto_PhoMVA90_Pt'], ch.GenPart_pt[i])
                        if(bUseVertex):
                            Fill1D(histos['RecoExtendedEle_CBVeto_PhoMVA90_DxyAbs'], abs(GenPart_dxy))
                    if(abRecoAsStandardCutBased[1] or abRecoAsIPPMVA[0]):
                        Fill1D(histos['RecoExtendedEle_CBLoose_PhoMVA90_Pt'], ch.GenPart_pt[i])
                        if(bUseVertex):
                            Fill1D(histos['RecoExtendedEle_CBLoose_PhoMVA90_DxyAbs'], abs(GenPart_dxy))
                    if( abRecoAsStandardCutBased[2] or abRecoAsIPPMVA[0] ):
                        Fill1D(histos['RecoExtendedEle_CBMedium_PhoMVA90_Pt'], ch.GenPart_pt[i])
                        if(bUseVertex):
                            Fill1D(histos['RecoExtendedEle_CBMedium_PhoMVA90_DxyAbs'], abs(GenPart_dxy))
                    if( abRecoAsStandardCutBased[3] or abRecoAsIPPMVA[0] ):
                        Fill1D(histos['RecoExtendedEle_CBTight_PhoMVA90_Pt'], ch.GenPart_pt[i])
                        if(bUseVertex):
                            Fill1D(histos['RecoExtendedEle_CBTight_PhoMVA90_DxyAbs'], abs(GenPart_dxy))

                    if(abRecoAsStandardCutBased[0] or abRecoAsIPPMVA[1]):
                        Fill1D(histos['RecoExtendedEle_CBVeto_PhoMVA80_Pt'], ch.GenPart_pt[i])
                        if(bUseVertex):
                            Fill1D(histos['RecoExtendedEle_CBVeto_PhoMVA80_DxyAbs'], abs(GenPart_dxy))
                    if(abRecoAsStandardCutBased[1] or abRecoAsIPPMVA[1]):
                        Fill1D(histos['RecoExtendedEle_CBLoose_PhoMVA80_Pt'], ch.GenPart_pt[i])
                        if(bUseVertex):
                            Fill1D(histos['RecoExtendedEle_CBLoose_PhoMVA80_DxyAbs'], abs(GenPart_dxy))
                    if( abRecoAsStandardCutBased[2] or abRecoAsIPPMVA[1] ):
                        Fill1D(histos['RecoExtendedEle_CBMedium_PhoMVA80_Pt'], ch.GenPart_pt[i])
                        if(bUseVertex):
                            Fill1D(histos['RecoExtendedEle_CBMedium_PhoMVA80_DxyAbs'], abs(GenPart_dxy))
                    if( abRecoAsStandardCutBased[3] or abRecoAsIPPMVA[1] ):
                        Fill1D(histos['RecoExtendedEle_CBTight_PhoMVA80_Pt'], ch.GenPart_pt[i])
                        if(bUseVertex):
                            Fill1D(histos['RecoExtendedEle_CBTight_PhoMVA80_DxyAbs'], abs(GenPart_dxy))




                    # Standard or LowPt
                    if(bRecoAsStandard or bRecoAsLowPt):
                        Fill1D(histos['RecoStandardEle||LowPtEle_Pt'], ch.GenPart_pt[i])
                        Fill1D(histos['RecoStandardEle||LowPtEle_Eta'], ch.GenPart_eta[i])
                        Fill1D(histos['RecoStandardEle||LowPtEle_Phi'], ch.GenPart_phi[i])

                        if(bUseVertex):
                            Fill1D(histos['RecoStandardEle||LowPtEle_Dxy'], GenPart_dxy)
                            Fill1D(histos['RecoStandardEle||LowPtEle_Dz'], GenPart_dz)
                            Fill1D(histos['RecoStandardEle||LowPtEle_DxyAbs'], abs(GenPart_dxy))
                            Fill1D(histos['RecoStandardEle||LowPtEle_DzAbs'], abs(GenPart_dz))

                    if(abRecoAsStandardCutBased[0] or bRecoAsLowPt):
                        Fill1D(histos['RecoStandardEle_CBVeto||LowPtEle_Pt'], ch.GenPart_pt[i])
                        if(bUseVertex):
                            Fill1D(histos['RecoStandardEle_CBVeto||LowPtEle_DxyAbs'], abs(GenPart_dxy))
                    if(abRecoAsStandardCutBased[1] or bRecoAsLowPt):
                        Fill1D(histos['RecoStandardEle_CBLoose||LowPtEle_Pt'], ch.GenPart_pt[i])
                        if(bUseVertex):
                            Fill1D(histos['RecoStandardEle_CBLoose||LowPtEle_DxyAbs'], abs(GenPart_dxy))
                    if( abRecoAsStandardCutBased[2] or bRecoAsLowPt ):
                        Fill1D(histos['RecoStandardEle_CBMedium||LowPtEle_Pt'], ch.GenPart_pt[i])
                        if(bUseVertex):
                            Fill1D(histos['RecoStandardEle_CBMedium||LowPtEle_DxyAbs'], abs(GenPart_dxy))
                    if( abRecoAsStandardCutBased[3] or bRecoAsLowPt ):
                        Fill1D(histos['RecoStandardEle_CBTight||LowPtEle_Pt'], ch.GenPart_pt[i])
                        if(bUseVertex):
                            Fill1D(histos['RecoStandardEle_CBTight||LowPtEle_DxyAbs'], abs(GenPart_dxy))

                    if(abRecoAsStandardMVA[0] or bRecoAsLowPt):
                        Fill1D(histos['RecoStandardEle_mvaNoIsoWPL||LowPtEle_Pt'], ch.GenPart_pt[i])
                        if(bUseVertex):
                            Fill1D(histos['RecoStandardEle_mvaNoIsoWPL||LowPtEle_DxyAbs'], abs(GenPart_dxy))
                    if(abRecoAsStandardMVA[1] or bRecoAsLowPt):
                        Fill1D(histos['RecoStandardEle_mvaNoIsoWP90||LowPtEle_Pt'], ch.GenPart_pt[i])
                        if(bUseVertex):
                            Fill1D(histos['RecoStandardEle_mvaNoIsoWP90||LowPtEle_DxyAbs'], abs(GenPart_dxy))
                    if( abRecoAsStandardMVA[2] or bRecoAsLowPt ):
                        Fill1D(histos['RecoStandardEle_mvaNoIsoWP80||LowPtEle_Pt'], ch.GenPart_pt[i])
                        if(bUseVertex):
                            Fill1D(histos['RecoStandardEle_mvaNoIsoWP80||LowPtEle_DxyAbs'], abs(GenPart_dxy))

                    for k in range(len(bdtcuts)):
                        if(bRecoAsStandard or abRecoAsLowPtBDT[k]):
                            Fill1D(histos['RecoStandardEle||LowPtEle_Pt_BDT'+str(bdtcuts[k])], ch.GenPart_pt[i])
                            if(bUseVertex):
                                Fill1D(histos['RecoStandardEle||LowPtEle_DxyAbs_BDT'+str(bdtcuts[k])], abs(GenPart_dxy))
                        if(abRecoAsStandardCutBased[0] or abRecoAsLowPtBDT[k]):
                            Fill1D(histos['RecoStandardEle_CBVeto||LowPtEle_Pt_BDT'+str(bdtcuts[k])], ch.GenPart_pt[i])
                            if(bUseVertex):
                                Fill1D(histos['RecoStandardEle_CBVeto||LowPtEle_DxyAbs_BDT'+str(bdtcuts[k])], abs(GenPart_dxy))
                        if(abRecoAsStandardCutBased[1] or abRecoAsLowPtBDT[k]):
                            Fill1D(histos['RecoStandardEle_CBLoose||LowPtEle_Pt_BDT'+str(bdtcuts[k])], ch.GenPart_pt[i])
                            if(bUseVertex):
                                Fill1D(histos['RecoStandardEle_CBLoose||LowPtEle_DxyAbs_BDT'+str(bdtcuts[k])], abs(GenPart_dxy))
                        if( abRecoAsStandardCutBased[2] or abRecoAsLowPtBDT[k] ):
                            Fill1D(histos['RecoStandardEle_CBMedium||LowPtEle_Pt_BDT'+str(bdtcuts[k])], ch.GenPart_pt[i])
                            if(bUseVertex):
                                Fill1D(histos['RecoStandardEle_CBMedium||LowPtEle_DxyAbs_BDT'+str(bdtcuts[k])], abs(GenPart_dxy))
                        if( abRecoAsStandardCutBased[3] or abRecoAsLowPtBDT[k] ):
                            Fill1D(histos['RecoStandardEle_CBTight||LowPtEle_Pt_BDT'+str(bdtcuts[k])], ch.GenPart_pt[i])
                            if(bUseVertex):
                                Fill1D(histos['RecoStandardEle_CBTight||LowPtEle_DxyAbs_BDT'+str(bdtcuts[k])], abs(GenPart_dxy))

                        if(abRecoAsStandardMVA[0] or abRecoAsLowPtBDT[k]):
                            Fill1D(histos['RecoStandardEle_mvaNoIsoWPL||LowPtEle_Pt_BDT'+str(bdtcuts[k])], ch.GenPart_pt[i])
                            if(bUseVertex):
                                Fill1D(histos['RecoStandardEle_mvaNoIsoWPL||LowPtEle_DxyAbs_BDT'+str(bdtcuts[k])], abs(GenPart_dxy))
                        if(abRecoAsStandardMVA[1] or abRecoAsLowPtBDT[k]):
                            Fill1D(histos['RecoStandardEle_mvaNoIsoWP90||LowPtEle_Pt_BDT'+str(bdtcuts[k])], ch.GenPart_pt[i])
                            if(bUseVertex):
                                Fill1D(histos['RecoStandardEle_mvaNoIsoWP90||LowPtEle_DxyAbs_BDT'+str(bdtcuts[k])], abs(GenPart_dxy))
                        if( abRecoAsStandardMVA[2] or abRecoAsLowPtBDT[k] ):
                            Fill1D(histos['RecoStandardEle_mvaNoIsoWP80||LowPtEle_Pt_BDT'+str(bdtcuts[k])], ch.GenPart_pt[i])
                            if(bUseVertex):
                                Fill1D(histos['RecoStandardEle_mvaNoIsoWP80||LowPtEle_DxyAbs_BDT'+str(bdtcuts[k])], abs(GenPart_dxy))



                        # All 3:
                        if(abRecoAsStandardMVA[0] or abRecoAsIPPMVA[0] or abRecoAsLowPtBDT[k]):
                            Fill1D(histos['RecoExtendedEle_mvaNoIsoWPL_PhoMVA90||LowPtEle_Pt_BDT'+str(bdtcuts[k])], ch.GenPart_pt[i])
                            if(bUseVertex):
                                Fill1D(histos['RecoExtendedEle_mvaNoIsoWPL_PhoMVA90||LowPtEle_DxyAbs_BDT'+str(bdtcuts[k])], abs(GenPart_dxy))
                        if(abRecoAsStandardCutBased[0] or abRecoAsIPPCutBased[0] or abRecoAsLowPtBDT[k]):
                            Fill1D(histos['RecoExtendedEle_CBVeto_PhoCBLoose||LowPtEle_Pt_BDT'+str(bdtcuts[k])], ch.GenPart_pt[i])
                            if(bUseVertex):
                                Fill1D(histos['RecoExtendedEle_CBVeto_PhoCBLoose||LowPtEle_DxyAbs_BDT'+str(bdtcuts[k])], abs(GenPart_dxy))
                        if(abRecoAsStandardCutBased[1] or abRecoAsIPPCutBased[1] or abRecoAsLowPtBDT[k]):
                            Fill1D(histos['RecoExtendedEle_CBLoose_PhoCBMedium||LowPtEle_Pt_BDT'+str(bdtcuts[k])], ch.GenPart_pt[i])
                            if(bUseVertex):
                                Fill1D(histos['RecoExtendedEle_CBLoose_PhoCBMedium||LowPtEle_DxyAbs_BDT'+str(bdtcuts[k])], abs(GenPart_dxy))

                    # overlap study for standard vs lowpt at the given special workingpoints
                    if(abRecoAsLowPtBDT[0] and abRecoAsStandardCutBased[0]):
                        Fill1D(histos['OverlapStandard_CBVeto&&LowPtEle_BDT0_Pt'], ch.GenPart_pt[i])
                        if(bUseVertex):
                            Fill1D(histos['OverlapStandard_CBVeto&&LowPtEle_BDT0_DxyAbs'], abs(GenPart_dxy))





                    # final overlap study for ipp:
                    if(abRecoAsIPPCutBased[0] and ( abRecoAsStandardCutBased[0] or abRecoAsLowPtBDT[0] )):
                        Fill1D(histos['OverlapIPP_CBLoose&&(Standard_CBVeto||LowPtEle_BDT0)_Pt'], ch.GenPart_pt[i])
                        if(bUseVertex):
                            Fill1D(histos['OverlapIPP_CBLoose&&(Standard_CBVeto||LowPtEle_BDT0)_DxyAbs'], abs(GenPart_dxy))



                    # Other recos:

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
                        Fill1D(histos['IsoTrack_Pt'], ch.GenPart_pt[i])
                        Fill1D(histos['IsoTrack_Eta'], ch.GenPart_eta[i])
                        Fill1D(histos['IsoTrack_Phi'], ch.GenPart_phi[i])

                        if(bUseVertex):
                            Fill1D(histos['IsoTrack_Dxy'], GenPart_dxy)
                            Fill1D(histos['IsoTrack_Dz'], GenPart_dz)
                            Fill1D(histos['IsoTrack_DxyAbs'], abs(GenPart_dxy))
                            Fill1D(histos['IsoTrack_DzAbs'], abs(GenPart_dz))
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
                        Fill1D(histos['RecoPhoton_Pt'], ch.GenPart_pt[i])
                        Fill1D(histos['RecoPhoton_Eta'], ch.GenPart_eta[i])
                        Fill1D(histos['RecoPhoton_Phi'], ch.GenPart_phi[i])
                        if(bUseVertex):
                            Fill1D(histos['RecoPhoton_Dxy'], GenPart_dxy)
                            Fill1D(histos['RecoPhoton_Dz'], GenPart_dz)
                            Fill1D(histos['RecoPhoton_DxyAbs'], abs(GenPart_dxy))
                            Fill1D(histos['RecoPhoton_DzAbs'], abs(GenPart_dz))

                    histos['dRControlPlot2D'].Fill(isodist,phodist)



        ############################
        ## BACKGROUND CALCULATION ##
        ############################

        # create background
        if(bBackground):
            background = []
            for iTR in range(ch.nIsoTrack):
                min_dRele = 1000
                minIEle = -1
                #min_dRbackgnd = 1000
                #minIBackgnd = -1

                # eta and pt cut also here:
                if(abs(ch.IsoTrack_eta[iTR]) > 2.5 or ch.IsoTrack_pt[iTR] < ptcut):
                    continue

                for iMC in range(ch.nGenPart):
                    temp_dR = dR(ch.GenPart_eta[iMC], ch.IsoTrack_eta[iTR], ch.GenPart_phi[iMC], ch.IsoTrack_phi[iTR])
                    if abs(ch.GenPart_pdgId[iMC]) == 11:
                        if( temp_dR < min_dRele):
                            min_dRele = temp_dR
                            minIEle = iMC
                    # suggestion by Ivan: dont do this. 
                    #else:
                    #    if(temp_dR < min_dRbackgnd):
                    #        min_dRbackgnd = temp_dR
                    #        minIBackgnd = iMC

                #if( dRcut(min_dRbackgnd,0.2)):
                if(not dRcut(min_dRele,0.2) ):
                    #background.append(minIBackgnd)
                    background.append(iTR)

                    pt = ch.IsoTrack_pt[iTR]
                    # fixme: what is up with overflow
                    if(pt > 99):
                        pt = 99.0
                    # this is the denom:
                    histos['TrueBK_Pt'].Fill( pt)
                    histos['TrueBK_Eta'].Fill( ch.IsoTrack_eta[iTR])
                    histos['TrueBK_Phi'].Fill( ch.IsoTrack_phi[iTR])
                    if(bUseVertex):
                        dxy = abs(ch.IsoTrack_dxy[iTR])
                        dz = abs(ch.IsoTrack_dz[iTR])

                        # fixme: overflow?
                        if(dxy > 14 and binningSet == 0):
                            dxy = 14
                        elif(dxy > 2.5 and binningSet == 1):
                            dxy = 2.48

                        if(dz > 14 and binningSet == 0):
                            dz = 14
                        elif(dz > 5 and binningSet == 1):
                            dz = 4.9

                        histos['TrueBK_DxyAbs'].Fill(dxy)
                        histos['TrueBK_DzAbs'].Fill(dz)

            # To get background rejection: redo the matching, but work on the true background as input sample.

            for iBK in range(len(background)):

                i = background[iBK]
                if(bUseVertex):
                    IsoTrack_dxyAbs = abs(ch.IsoTrack_dxy[i])
                    IsoTrack_dzAbs = abs(ch.IsoTrack_dz[i])

                    # fixme: overflow?
                    if(IsoTrack_dxyAbs > 14 and binningSet == 0):
                        IsoTrack_dxyAbs = 14
                    elif(IsoTrack_dxyAbs > 2.5 and binningSet == 1):
                        IsoTrack_dxyAbs = 2.48

                    if(IsoTrack_dzAbs > 14 and binningSet == 0):
                        IsoTrack_dzAbs = 14
                    elif(IsoTrack_dzAbs > 5 and binningSet == 1):
                        IsoTrack_dzAbs = 4.9
                pt = ch.IsoTrack_pt[i]
                if(pt > 99.0):
                    pt = 99.0



                #############################
                ## TEST RECO METHODS ON BK ##
                #############################

                # Try as Standard:

                bRejectedAsStandard = False
                abRejectedAsStandardCutBased = [False, False, False, False]
                abRejectedAsStandardMVA = [False, False, False]
                eledist = 100000
                eleIdx = -1
                idx0 = -1
                for j in range(ch.nElectron):

                    # remove isolation
                    #eleVID(ch.Electron_vidNestedWPBitmap[j], 0, removedCuts=['pfRelIso03_all'])
                    dist0 = dR(ch.IsoTrack_eta[i], ch.Electron_eta[j], ch.IsoTrack_phi[i], ch.Electron_phi[j] ) 
                    if( dist0 < eledist):
                        eledist = dist0
                        idx0 = j
                # not dRcut: background REJECTION efficiency
                if( not dRcut(eledist) ):
                    bRejectedAsStandard = True
                else:
                    eleIdx = idx0 
                misidentifiedAsStandardElectron = eleIdx



                # standard with cut based id
                for iCut in range(4):
                    eledist = 100000
                    idx0 = -1
                    for j in range(ch.nElectron):
                        #if(ch.Electron_cutBased[j] >> iCut & 1):
                        #if(ch.Electron_cutBased[j] > iCut ):
                        if(eleVID( ch.Electron_vidNestedWPBitmap[j], iCut + 1, removedCuts=['pfRelIso03_all'])) :
                            dist0 = dR(ch.IsoTrack_eta[i], ch.Electron_eta[j], ch.IsoTrack_phi[i], ch.Electron_phi[j] ) 
                            if( dist0 < eledist):
                                eledist = dist0
                                idx0 = j
                    # not dRcut: background REJECTION efficiency
                    if( not dRcut(eledist) ):
                        abRejectedAsStandardCutBased[iCut] = True



                # with mva ID:
                eledist = 100000
                idx0 = -1
                for j in range(ch.nElectron):
                    if(ch.Electron_mvaFall17V2noIso_WPL[j]):
                        dist0 = dR(ch.IsoTrack_eta[i], ch.Electron_eta[j], ch.IsoTrack_phi[i], ch.Electron_phi[j] ) 
                        if( dist0 < eledist):
                            eledist = dist0
                            idx0 = j
                if( not dRcut(eledist) ):
                    abRejectedAsStandardMVA[0] = True

                eledist = 100000
                idx0 = -1
                for j in range(ch.nElectron):
                    if(ch.Electron_mvaFall17V2noIso_WP90[j]):
                        dist0 = dR(ch.IsoTrack_eta[i], ch.Electron_eta[j], ch.IsoTrack_phi[i], ch.Electron_phi[j] ) 
                        if( dist0 < eledist):
                            eledist = dist0
                            idx0 = j
                if( not dRcut(eledist) ):
                    abRejectedAsStandardMVA[1] = True

                eledist = 100000
                idx0 = -1
                for j in range(ch.nElectron):
                    if(ch.Electron_mvaFall17V2noIso_WP80[j]):
                        dist0 = dR(ch.IsoTrack_eta[i], ch.Electron_eta[j], ch.IsoTrack_phi[i], ch.Electron_phi[j] ) 
                        if( dist0 < eledist):
                            eledist = dist0
                            idx0 = j
                # not dRcut: background REJECTION efficiency
                if( not dRcut(eledist) ):
                    abRejectedAsStandardMVA[2] = True



                # V9: Try as LowPt Electron

                misidentifiedAsLowPtElectron = -1
                bRejectedAsLowPt = False
                abRejectedAsLowPtBDT = [False] * len(bdtcuts)
                bdtidx = [-1] * len(bdtcuts)
                bdtdist = [100000] * len(bdtcuts)
                if(bUseV9):
                    eledist = 100000
                    idx0 = -1
                    for j in range(ch.nLowPtElectron):

                        dist0 = dR(ch.IsoTrack_eta[i], ch.LowPtElectron_eta[j], ch.IsoTrack_phi[i], ch.LowPtElectron_phi[j] ) 
                        if( dist0 < eledist):
                            eledist = dist0
                            idx0 = j

                        for k in range(len(bdtcuts)):
                            if(bdtdist[k] > dist0 and ch.LowPtElectron_ID[j] > bdtcuts[k]):
                                bdtdist[k] = dist0
                                bdtidx[k] = j

                    if( not dRcut(eledist) ):
                        bRejectedAsLowPt = True
                    else:
                        misidentifiedAsLowPtElectron = idx0

                        Fill1D(histos['MatchedBKBDTID'], ch.LowPtElectron_ID[misidentifiedAsLowPtElectron])
                        Fill1D(histos['MatchedBKBDTembeddedID'], ch.LowPtElectron_embeddedID[misidentifiedAsLowPtElectron])
                        histos['MatchedBKBDT2D'].Fill(ch.LowPtElectron_ID[misidentifiedAsLowPtElectron], ch.LowPtElectron_embeddedID[misidentifiedAsLowPtElectron])

                    for k in range(len(bdtcuts)):
                        if(not dRcut(bdtdist[k])):
                            abRejectedAsLowPtBDT[k] = True

                    
                # Try as IPP
                bRejectedAsIPP = False
                abRejectedAsIPPCutBased = [False,False,False]
                abRejectedAsIPPMVA = [False,False]
                fPairDist = 100000
                iPairIdx = -1
                for j in range(len(aPhotonIsoTrackPairs)):
                    dist0 = dR(ch.IsoTrack_eta[i], aPhotonIsoTrackPairs[j].GetEta(), ch.IsoTrack_phi[i], aPhotonIsoTrackPairs[j].GetPhi() ) 
                    if( dist0 < fPairDist):
                        fPairDist = dist0
                        iPairIdx = j
                if( not dRcut(fPairDist, 0.2) ):
                    bRejectedAsIPP = True

                fPairDist = 100000
                iPairIdx = -1
                for j in range(len(aPhotonIsoTrackPairs_CBLoose)):
                    dist0 = dR(ch.IsoTrack_eta[i], aPhotonIsoTrackPairs_CBLoose[j].GetEta(), ch.IsoTrack_phi[i], aPhotonIsoTrackPairs_CBLoose[j].GetPhi() ) 
                    if( dist0 < fPairDist):
                        fPairDist = dist0
                        iPairIdx = j
                if( not dRcut(fPairDist, 0.2) ):
                    abRejectedAsIPPCutBased[0] = True

                fPairDist = 100000
                iPairIdx = -1
                for j in range(len(aPhotonIsoTrackPairs_CBMedium)):
                    dist0 = dR(ch.IsoTrack_eta[i], aPhotonIsoTrackPairs_CBMedium[j].GetEta(), ch.IsoTrack_phi[i], aPhotonIsoTrackPairs_CBMedium[j].GetPhi() ) 
                    if( dist0 < fPairDist):
                        fPairDist = dist0
                        iPairIdx = j
                if( not dRcut(fPairDist, 0.2) ):
                    abRejectedAsIPPCutBased[1] = True

                fPairDist = 100000
                iPairIdx = -1
                for j in range(len(aPhotonIsoTrackPairs_CBTight)):
                    dist0 = dR(ch.IsoTrack_eta[i], aPhotonIsoTrackPairs_CBTight[j].GetEta(), ch.IsoTrack_phi[i], aPhotonIsoTrackPairs_CBTight[j].GetPhi() ) 
                    if( dist0 < fPairDist):
                        fPairDist = dist0
                        iPairIdx = j
                if( not dRcut(fPairDist, 0.2) ):
                    abRejectedAsIPPCutBased[2] = True

                fPairDist = 100000
                iPairIdx = -1
                for j in range(len(aPhotonIsoTrackPairs_MVA90)):
                    dist0 = dR(ch.IsoTrack_eta[i], aPhotonIsoTrackPairs_MVA90[j].GetEta(), ch.IsoTrack_phi[i], aPhotonIsoTrackPairs_MVA90[j].GetPhi() ) 
                    if( dist0 < fPairDist):
                        fPairDist = dist0
                        iPairIdx = j
                if( not dRcut(fPairDist, 0.2) ):
                    abRejectedAsIPPMVA[0] = True

                fPairDist = 100000
                iPairIdx = -1
                for j in range(len(aPhotonIsoTrackPairs_MVA80)):
                    dist0 = dR(ch.IsoTrack_eta[i], aPhotonIsoTrackPairs_MVA80[j].GetEta(), ch.IsoTrack_phi[i], aPhotonIsoTrackPairs_MVA80[j].GetPhi() ) 
                    if( dist0 < fPairDist):
                        fPairDist = dist0
                        iPairIdx = j
                if( not dRcut(fPairDist, 0.2) ):
                    abRejectedAsIPPMVA[1] = True

                ###############################
                ## FILL BK REJECT HISTOGRAMS ##
                ###############################

                if(bRejectedAsStandard):
                    # idk what is up with overflow, but for some reason it needs to be added here:
                    histos['RejectBKasStandard_Pt'].Fill( pt)
                    histos['RejectBKasStandard_Eta'].Fill( ch.IsoTrack_eta[i])
                    histos['RejectBKasStandard_Phi'].Fill( ch.IsoTrack_phi[i])
                    if(bUseVertex):
                        histos['RejectBKasStandard_DxyAbs'].Fill( IsoTrack_dxyAbs)
                        histos['RejectBKasStandard_DzAbs'].Fill( IsoTrack_dzAbs)

                if(abRejectedAsStandardCutBased[0]):
                    # idk what is up with overflow, but for some reason it needs to be added here:
                    histos['RejectBKasStandard_CBVeto_Pt'].Fill( pt)
                    if(bUseVertex):
                        histos['RejectBKasStandard_CBVeto_DxyAbs'].Fill( IsoTrack_dxyAbs)
                if( abRejectedAsStandardCutBased[1] ):
                    histos['RejectBKasStandard_CBLoose_Pt'].Fill( pt)
                    if(bUseVertex):
                        histos['RejectBKasStandard_CBLoose_DxyAbs'].Fill( IsoTrack_dxyAbs)
                if( abRejectedAsStandardCutBased[2] ):
                    histos['RejectBKasStandard_CBMedium_Pt'].Fill( pt)
                    if(bUseVertex):
                        histos['RejectBKasStandard_CBMedium_DxyAbs'].Fill( IsoTrack_dxyAbs)
                if( abRejectedAsStandardCutBased[3] ):
                    histos['RejectBKasStandard_CBTight_Pt'].Fill( pt)
                    if(bUseVertex):
                        histos['RejectBKasStandard_CBTight_DxyAbs'].Fill( IsoTrack_dxyAbs)
  

                if(abRejectedAsStandardMVA[0]):
                    # idk what is up with overflow, but for some reason it needs to be added here:
                    histos['RejectBKasStandard_mvaNoIsoWPL_Pt'].Fill( pt)
                    histos['RejectBKasStandard_mvaNoIsoWPL_Eta'].Fill( ch.IsoTrack_eta[i])
                    histos['RejectBKasStandard_mvaNoIsoWPL_Phi'].Fill( ch.IsoTrack_phi[i])
                    if(bUseVertex):
                        histos['RejectBKasStandard_mvaNoIsoWPL_DxyAbs'].Fill( IsoTrack_dxyAbs)
                        histos['RejectBKasStandard_mvaNoIsoWPL_DzAbs'].Fill( IsoTrack_dzAbs)

                if(abRejectedAsStandardMVA[1]):
                    # idk what is up with overflow, but for some reason it needs to be added here:
                    histos['RejectBKasStandard_mvaNoIsoWP90_Pt'].Fill( pt)
                    histos['RejectBKasStandard_mvaNoIsoWP90_Eta'].Fill( ch.IsoTrack_eta[i])
                    histos['RejectBKasStandard_mvaNoIsoWP90_Phi'].Fill( ch.IsoTrack_phi[i])
                    if(bUseVertex):
                        histos['RejectBKasStandard_mvaNoIsoWP90_DxyAbs'].Fill( IsoTrack_dxyAbs)
                        histos['RejectBKasStandard_mvaNoIsoWP90_DzAbs'].Fill( IsoTrack_dzAbs)

                if(abRejectedAsStandardMVA[2]):
                    # idk what is up with overflow, but for some reason it needs to be added here:
                    histos['RejectBKasStandard_mvaNoIsoWP80_Pt'].Fill( pt)
                    histos['RejectBKasStandard_mvaNoIsoWP80_Eta'].Fill( ch.IsoTrack_eta[i])
                    histos['RejectBKasStandard_mvaNoIsoWP80_Phi'].Fill( ch.IsoTrack_phi[i])
                    if(bUseVertex):
                        histos['RejectBKasStandard_mvaNoIsoWP80_DxyAbs'].Fill( IsoTrack_dxyAbs)
                        histos['RejectBKasStandard_mvaNoIsoWP80_DzAbs'].Fill( IsoTrack_dzAbs)



                if( bRejectedAsLowPt ):
                    # record lowptElectron reco efficiency, on its own
                    Fill1D(histos['RejectBKasLowPt_Pt'], pt)
                    Fill1D(histos['RejectBKasLowPt_Eta'], ch.IsoTrack_eta[i])
                    Fill1D(histos['RejectBKasLowPt_Phi'], ch.IsoTrack_phi[i])
                    if(bUseVertex):
                        Fill1D(histos['RejectBKasLowPt_DxyAbs'], IsoTrack_dxyAbs)
                        Fill1D(histos['RejectBKasLowPt_DzAbs'], IsoTrack_dzAbs)

                    # if also no standard ele: rejected Standard+Lowpt ->> moved to end

                    
                for k in range(len(bdtcuts)):
                    if(abRejectedAsLowPtBDT[k]):
                        Fill1D(histos['RejectBKasLowPt_Pt_BDT'+str(bdtcuts[k])], pt)
                        if(bUseVertex):
                            Fill1D(histos['RejectBKasLowPt_DxyAbs_BDT'+str(bdtcuts[k])], IsoTrack_dxyAbs)



                if(bRejectedAsIPP):
                    Fill1D(histos['RejectBKasIPP_Pt'], pt)
                    Fill1D(histos['RejectBKasIPP_Eta'], ch.IsoTrack_eta[i])
                    Fill1D(histos['RejectBKasIPP_Phi'], ch.IsoTrack_phi[i])

                    if(bUseVertex):
                        Fill1D(histos['RejectBKasIPP_DxyAbs'], IsoTrack_dxyAbs)
                        Fill1D(histos['RejectBKasIPP_DzAbs'], IsoTrack_dzAbs)



                    if( misidentifiedAsStandardElectron == -1):
                        Fill1D(histos['RejectBKasExtended_Pt'], pt)
                        Fill1D(histos['RejectBKasExtended_Eta'], ch.IsoTrack_eta[i])
                        Fill1D(histos['RejectBKasExtended_Phi'], ch.IsoTrack_phi[i])

                        if(bUseVertex):
                            Fill1D(histos['RejectBKasExtended_DxyAbs'], IsoTrack_dxyAbs)
                            Fill1D(histos['RejectBKasExtended_DzAbs'], IsoTrack_dzAbs)

        
                        # if also no lowpt, then rejected extended+lowpt
                        if(misidentifiedAsLowPtElectron == -1):
                            Fill1D(histos['RejectBKasExtended&&LowPt_Pt'], pt)
                            Fill1D(histos['RejectBKasExtended&&LowPt_Eta'], ch.IsoTrack_eta[i])
                            Fill1D(histos['RejectBKasExtended&&LowPt_Phi'], ch.IsoTrack_phi[i])

                            if(bUseVertex):
                                Fill1D(histos['RejectBKasExtended&&LowPt_DxyAbs'], IsoTrack_dxyAbs)
                                Fill1D(histos['RejectBKasExtended&&LowPt_DzAbs'], IsoTrack_dzAbs)


                    if( misidentifiedAsLowPtElectron == -1):
                        Fill1D(histos['RejectBKasIPP&&LowPt_Pt'], pt)
                        Fill1D(histos['RejectBKasIPP&&LowPt_Eta'], ch.IsoTrack_eta[i])
                        Fill1D(histos['RejectBKasIPP&&LowPt_Phi'], ch.IsoTrack_phi[i])

                        if(bUseVertex):
                            Fill1D(histos['RejectBKasIPP&&LowPt_DxyAbs'], IsoTrack_dxyAbs)
                            Fill1D(histos['RejectBKasIPP&&LowPt_DzAbs'], IsoTrack_dzAbs)

    
                if(abRejectedAsIPPCutBased[0]):
                    Fill1D(histos['RejectBKasIPP_CBLoose_Pt'], pt)
                    if(bUseVertex):
                        Fill1D(histos['RejectBKasIPP_CBLoose_DxyAbs'], IsoTrack_dxyAbs)
                if(abRejectedAsIPPCutBased[1]):
                    Fill1D(histos['RejectBKasIPP_CBMedium_Pt'], pt)
                    if(bUseVertex):
                        Fill1D(histos['RejectBKasIPP_CBMedium_DxyAbs'], IsoTrack_dxyAbs)
                if(abRejectedAsIPPCutBased[2]):
                    Fill1D(histos['RejectBKasIPP_CBTight_Pt'], pt)
                    if(bUseVertex):
                        Fill1D(histos['RejectBKasIPP_CBTight_DxyAbs'], IsoTrack_dxyAbs)
                if(abRejectedAsIPPCutBased[0]):
                    Fill1D(histos['RejectBKasIPP_MVA90_Pt'], pt)
                    if(bUseVertex):
                        Fill1D(histos['RejectBKasIPP_MVA90_DxyAbs'], IsoTrack_dxyAbs)
                if(abRejectedAsIPPCutBased[1]):
                    Fill1D(histos['RejectBKasIPP_MVA80_Pt'], pt)
                    if(bUseVertex):
                        Fill1D(histos['RejectBKasIPP_MVA80_DxyAbs'], IsoTrack_dxyAbs)


                for k in range(len(bdtcuts)):
                    if(abRejectedAsLowPtBDT[k] and bRejectedAsIPP):
                        Fill1D(histos['RejectBKasIPP&&LowPt_Pt_BDT'+str(bdtcuts[k])], pt)
                        if(bUseVertex):
                            Fill1D(histos['RejectBKasIPP&&LowPt_DxyAbs_BDT'+str(bdtcuts[k])], IsoTrack_dxyAbs)

                    if(abRejectedAsLowPtBDT[k] and abRejectedAsIPPCutBased[0]):
                        Fill1D(histos['RejectBKasIPP_CBLoose&&LowPt_Pt_BDT'+str(bdtcuts[k])], pt)
                        if(bUseVertex):
                            Fill1D(histos['RejectBKasIPP_CBLoose&&LowPt_DxyAbs_BDT'+str(bdtcuts[k])], IsoTrack_dxyAbs)

                    if(abRejectedAsLowPtBDT[k] and abRejectedAsIPPCutBased[1]):
                        Fill1D(histos['RejectBKasIPP_CBMedium&&LowPt_Pt_BDT'+str(bdtcuts[k])], pt)
                        if(bUseVertex):
                            Fill1D(histos['RejectBKasIPP_CBMedium&&LowPt_DxyAbs_BDT'+str(bdtcuts[k])], IsoTrack_dxyAbs)

                    if(abRejectedAsLowPtBDT[k] and abRejectedAsIPPCutBased[2]):
                        Fill1D(histos['RejectBKasIPP_CBTight&&LowPt_Pt_BDT'+str(bdtcuts[k])], pt)
                        if(bUseVertex):
                            Fill1D(histos['RejectBKasIPP_CBTight&&LowPt_DxyAbs_BDT'+str(bdtcuts[k])], IsoTrack_dxyAbs)

                    if(abRejectedAsLowPtBDT[k] and abRejectedAsIPPMVA[0]):
                        Fill1D(histos['RejectBKasIPP_MVA90&&LowPt_Pt_BDT'+str(bdtcuts[k])], pt)
                        if(bUseVertex):
                            Fill1D(histos['RejectBKasIPP_MVA90&&LowPt_DxyAbs_BDT'+str(bdtcuts[k])], IsoTrack_dxyAbs)

                    if(abRejectedAsLowPtBDT[k] and abRejectedAsIPPMVA[1]):
                        Fill1D(histos['RejectBKasIPP_MVA80&&LowPt_Pt_BDT'+str(bdtcuts[k])], pt)
                        if(bUseVertex):
                            Fill1D(histos['RejectBKasIPP_MVA80&&LowPt_DxyAbs_BDT'+str(bdtcuts[k])], IsoTrack_dxyAbs)




                if(abRejectedAsStandardCutBased[0] and bRejectedAsIPP):
                    # idk what is up with overflow, but for some reason it needs to be added here:
                    histos['RejectBKasExtended_CBVeto_Pt'].Fill( pt)
                    if(bUseVertex):
                        histos['RejectBKasExtended_CBVeto_DxyAbs'].Fill( IsoTrack_dxyAbs)
                if( abRejectedAsStandardCutBased[1]  and bRejectedAsIPP):
                    histos['RejectBKasExtended_CBLoose_Pt'].Fill( pt)
                    if(bUseVertex):
                        histos['RejectBKasExtended_CBLoose_DxyAbs'].Fill( IsoTrack_dxyAbs)
                if( abRejectedAsStandardCutBased[2]  and bRejectedAsIPP):
                    histos['RejectBKasExtended_CBMedium_Pt'].Fill( pt)
                    if(bUseVertex):
                        histos['RejectBKasExtended_CBMedium_DxyAbs'].Fill( IsoTrack_dxyAbs)
                if( abRejectedAsStandardCutBased[3]  and bRejectedAsIPP):
                    histos['RejectBKasExtended_CBTight_Pt'].Fill( pt)
                    if(bUseVertex):
                        histos['RejectBKasExtended_CBTight_DxyAbs'].Fill( IsoTrack_dxyAbs)


                if(abRejectedAsStandardCutBased[0] and abRejectedAsIPPCutBased[0]):
                    histos['RejectBKasExtended_CBVeto_PhoCBLoose_Pt'].Fill( pt)
                    if(bUseVertex):
                        histos['RejectBKasExtended_CBVeto_PhoCBLoose_DxyAbs'].Fill( IsoTrack_dxyAbs)
                if( abRejectedAsStandardCutBased[1]  and abRejectedAsIPPCutBased[0]):
                    histos['RejectBKasExtended_CBLoose_PhoCBLoose_Pt'].Fill( pt)
                    if(bUseVertex):
                        histos['RejectBKasExtended_CBLoose_PhoCBLoose_DxyAbs'].Fill( IsoTrack_dxyAbs)
                if( abRejectedAsStandardCutBased[2]  and abRejectedAsIPPCutBased[0]):
                    histos['RejectBKasExtended_CBMedium_PhoCBLoose_Pt'].Fill( pt)
                    if(bUseVertex):
                        histos['RejectBKasExtended_CBMedium_PhoCBLoose_DxyAbs'].Fill( IsoTrack_dxyAbs)
                if( abRejectedAsStandardCutBased[3]  and abRejectedAsIPPCutBased[0]):
                    histos['RejectBKasExtended_CBTight_PhoCBLoose_Pt'].Fill( pt)
                    if(bUseVertex):
                        histos['RejectBKasExtended_CBTight_PhoCBLoose_DxyAbs'].Fill( IsoTrack_dxyAbs)


                if(abRejectedAsStandardCutBased[0] and abRejectedAsIPPCutBased[1]):
                    histos['RejectBKasExtended_CBVeto_PhoCBMedium_Pt'].Fill( pt)
                    if(bUseVertex):
                        histos['RejectBKasExtended_CBVeto_PhoCBMedium_DxyAbs'].Fill( IsoTrack_dxyAbs)
                if( abRejectedAsStandardCutBased[1]  and abRejectedAsIPPCutBased[1]):
                    histos['RejectBKasExtended_CBLoose_PhoCBMedium_Pt'].Fill( pt)
                    if(bUseVertex):
                        histos['RejectBKasExtended_CBLoose_PhoCBMedium_DxyAbs'].Fill( IsoTrack_dxyAbs)
                if( abRejectedAsStandardCutBased[2]  and abRejectedAsIPPCutBased[1]):
                    histos['RejectBKasExtended_CBMedium_PhoCBMedium_Pt'].Fill( pt)
                    if(bUseVertex):
                        histos['RejectBKasExtended_CBMedium_PhoCBMedium_DxyAbs'].Fill( IsoTrack_dxyAbs)
                if( abRejectedAsStandardCutBased[3]  and abRejectedAsIPPCutBased[1]):
                    histos['RejectBKasExtended_CBTight_PhoCBMedium_Pt'].Fill( pt)
                    if(bUseVertex):
                        histos['RejectBKasExtended_CBTight_PhoCBMedium_DxyAbs'].Fill( IsoTrack_dxyAbs)

                if(abRejectedAsStandardCutBased[0] and abRejectedAsIPPCutBased[2]):
                    histos['RejectBKasExtended_CBVeto_PhoCBTight_Pt'].Fill( pt)
                    if(bUseVertex):
                        histos['RejectBKasExtended_CBVeto_PhoCBTight_DxyAbs'].Fill( IsoTrack_dxyAbs)
                if( abRejectedAsStandardCutBased[1]  and abRejectedAsIPPCutBased[2]):
                    histos['RejectBKasExtended_CBLoose_PhoCBTight_Pt'].Fill( pt)
                    if(bUseVertex):
                        histos['RejectBKasExtended_CBLoose_PhoCBTight_DxyAbs'].Fill( IsoTrack_dxyAbs)
                if( abRejectedAsStandardCutBased[2]  and abRejectedAsIPPCutBased[2]):
                    histos['RejectBKasExtended_CBMedium_PhoCBTight_Pt'].Fill( pt)
                    if(bUseVertex):
                        histos['RejectBKasExtended_CBMedium_PhoCBTight_DxyAbs'].Fill( IsoTrack_dxyAbs)
                if( abRejectedAsStandardCutBased[3]  and abRejectedAsIPPCutBased[2]):
                    histos['RejectBKasExtended_CBTight_PhoCBTight_Pt'].Fill( pt)
                    if(bUseVertex):
                        histos['RejectBKasExtended_CBTight_PhoCBTight_DxyAbs'].Fill( IsoTrack_dxyAbs)

                if(abRejectedAsStandardCutBased[0] and abRejectedAsIPPMVA[0]):
                    histos['RejectBKasExtended_CBVeto_PhoMVA90_Pt'].Fill( pt)
                    if(bUseVertex):
                        histos['RejectBKasExtended_CBVeto_PhoMVA90_DxyAbs'].Fill( IsoTrack_dxyAbs)
                if( abRejectedAsStandardCutBased[1]  and abRejectedAsIPPMVA[0]):
                    histos['RejectBKasExtended_CBLoose_PhoMVA90_Pt'].Fill( pt)
                    if(bUseVertex):
                        histos['RejectBKasExtended_CBLoose_PhoMVA90_DxyAbs'].Fill( IsoTrack_dxyAbs)
                if( abRejectedAsStandardCutBased[2]  and abRejectedAsIPPMVA[0]):
                    histos['RejectBKasExtended_CBMedium_PhoMVA90_Pt'].Fill( pt)
                    if(bUseVertex):
                        histos['RejectBKasExtended_CBMedium_PhoMVA90_DxyAbs'].Fill( IsoTrack_dxyAbs)
                if( abRejectedAsStandardCutBased[3]  and abRejectedAsIPPMVA[0]):
                    histos['RejectBKasExtended_CBTight_PhoMVA90_Pt'].Fill( pt)
                    if(bUseVertex):
                        histos['RejectBKasExtended_CBTight_PhoMVA90_DxyAbs'].Fill( IsoTrack_dxyAbs)

                if(abRejectedAsStandardCutBased[0] and abRejectedAsIPPMVA[1]):
                    histos['RejectBKasExtended_CBVeto_PhoMVA80_Pt'].Fill( pt)
                    if(bUseVertex):
                        histos['RejectBKasExtended_CBVeto_PhoMVA80_DxyAbs'].Fill( IsoTrack_dxyAbs)
                if( abRejectedAsStandardCutBased[1]  and abRejectedAsIPPMVA[1]):
                    histos['RejectBKasExtended_CBLoose_PhoMVA80_Pt'].Fill( pt)
                    if(bUseVertex):
                        histos['RejectBKasExtended_CBLoose_PhoMVA80_DxyAbs'].Fill( IsoTrack_dxyAbs)
                if( abRejectedAsStandardCutBased[2]  and abRejectedAsIPPMVA[1]):
                    histos['RejectBKasExtended_CBMedium_PhoMVA80_Pt'].Fill( pt)
                    if(bUseVertex):
                        histos['RejectBKasExtended_CBMedium_PhoMVA80_DxyAbs'].Fill( IsoTrack_dxyAbs)
                if( abRejectedAsStandardCutBased[3]  and abRejectedAsIPPMVA[1]):
                    histos['RejectBKasExtended_CBTight_PhoMVA80_Pt'].Fill( pt)
                    if(bUseVertex):
                        histos['RejectBKasExtended_CBTight_PhoMVA80_DxyAbs'].Fill( IsoTrack_dxyAbs)





                # Standard and LowPt
                if(bRejectedAsStandard and bRejectedAsLowPt):
                    Fill1D(histos['RejectBKasStandard&&LowPt_Pt'], pt)
                    Fill1D(histos['RejectBKasStandard&&LowPt_Eta'], ch.IsoTrack_eta[i])
                    Fill1D(histos['RejectBKasStandard&&LowPt_Phi'], ch.IsoTrack_phi[i])
                    if(bUseVertex):
                        Fill1D(histos['RejectBKasStandard&&LowPt_DxyAbs'], IsoTrack_dxyAbs)
                        Fill1D(histos['RejectBKasStandard&&LowPt_DzAbs'], IsoTrack_dzAbs)

                if(abRejectedAsStandardCutBased[0] and bRejectedAsLowPt):
                    Fill1D(histos['RejectBKasStandard_CBVeto&&LowPt_Pt'], pt)
                    if(bUseVertex):
                        Fill1D(histos['RejectBKasStandard_CBVeto&&LowPt_DxyAbs'], IsoTrack_dxyAbs)
                if(abRejectedAsStandardCutBased[1] and bRejectedAsLowPt):
                    Fill1D(histos['RejectBKasStandard_CBLoose&&LowPt_Pt'], pt)
                    if(bUseVertex):
                        Fill1D(histos['RejectBKasStandard_CBLoose&&LowPt_DxyAbs'], IsoTrack_dxyAbs)
                if( abRejectedAsStandardCutBased[2] and bRejectedAsLowPt ):
                    Fill1D(histos['RejectBKasStandard_CBMedium&&LowPt_Pt'], pt)
                    if(bUseVertex):
                        Fill1D(histos['RejectBKasStandard_CBMedium&&LowPt_DxyAbs'], IsoTrack_dxyAbs)
                if( abRejectedAsStandardCutBased[3] and bRejectedAsLowPt ):
                    Fill1D(histos['RejectBKasStandard_CBTight&&LowPt_Pt'], pt)
                    if(bUseVertex):
                        Fill1D(histos['RejectBKasStandard_CBTight&&LowPt_DxyAbs'], IsoTrack_dxyAbs)

                if(abRejectedAsStandardMVA[0] and bRejectedAsLowPt):
                    Fill1D(histos['RejectBKasStandard_mvaNoIsoWPL&&LowPt_Pt'], pt)
                    if(bUseVertex):
                        Fill1D(histos['RejectBKasStandard_mvaNoIsoWPL&&LowPt_DxyAbs'], IsoTrack_dxyAbs)
                if(abRejectedAsStandardMVA[1] and bRejectedAsLowPt):
                    Fill1D(histos['RejectBKasStandard_mvaNoIsoWP90&&LowPt_Pt'], pt)
                    if(bUseVertex):
                        Fill1D(histos['RejectBKasStandard_mvaNoIsoWP90&&LowPt_DxyAbs'], IsoTrack_dxyAbs)
                if( abRejectedAsStandardMVA[2] and bRejectedAsLowPt ):
                    Fill1D(histos['RejectBKasStandard_mvaNoIsoWP80&&LowPt_Pt'], pt)
                    if(bUseVertex):
                        Fill1D(histos['RejectBKasStandard_mvaNoIsoWP80&&LowPt_DxyAbs'], IsoTrack_dxyAbs)


                        

                for k in range(len(bdtcuts)):
                    if(bRejectedAsStandard and abRejectedAsLowPtBDT[k]):
                        Fill1D(histos['RejectBKasStandard&&LowPt_Pt_BDT'+str(bdtcuts[k])], pt)
                        if(bUseVertex):
                            Fill1D(histos['RejectBKasStandard&&LowPt_DxyAbs_BDT'+str(bdtcuts[k])], IsoTrack_dxyAbs)
                    if(abRejectedAsStandardCutBased[0] and abRejectedAsLowPtBDT[k]):
                        Fill1D(histos['RejectBKasStandard_CBVeto&&LowPt_Pt_BDT'+str(bdtcuts[k])], pt)
                        if(bUseVertex):
                            Fill1D(histos['RejectBKasStandard_CBVeto&&LowPt_DxyAbs_BDT'+str(bdtcuts[k])], IsoTrack_dxyAbs)
                    if(abRejectedAsStandardCutBased[1] and abRejectedAsLowPtBDT[k]):
                        Fill1D(histos['RejectBKasStandard_CBLoose&&LowPt_Pt_BDT'+str(bdtcuts[k])], pt)
                        if(bUseVertex):
                            Fill1D(histos['RejectBKasStandard_CBLoose&&LowPt_DxyAbs_BDT'+str(bdtcuts[k])], IsoTrack_dxyAbs)
                    if( abRejectedAsStandardCutBased[2] and abRejectedAsLowPtBDT[k] ):
                        Fill1D(histos['RejectBKasStandard_CBMedium&&LowPt_Pt_BDT'+str(bdtcuts[k])], pt)
                        if(bUseVertex):
                            Fill1D(histos['RejectBKasStandard_CBMedium&&LowPt_DxyAbs_BDT'+str(bdtcuts[k])], IsoTrack_dxyAbs)
                    if( abRejectedAsStandardCutBased[3] and abRejectedAsLowPtBDT[k] ):
                        Fill1D(histos['RejectBKasStandard_CBTight&&LowPt_Pt_BDT'+str(bdtcuts[k])], pt)
                        if(bUseVertex):
                            Fill1D(histos['RejectBKasStandard_CBTight&&LowPt_DxyAbs_BDT'+str(bdtcuts[k])], IsoTrack_dxyAbs)

                    if(abRejectedAsStandardMVA[0] and abRejectedAsLowPtBDT[k]):
                        Fill1D(histos['RejectBKasStandard_mvaNoIsoWPL&&LowPt_Pt_BDT'+str(bdtcuts[k])], pt)
                        if(bUseVertex):
                            Fill1D(histos['RejectBKasStandard_mvaNoIsoWPL&&LowPt_DxyAbs_BDT'+str(bdtcuts[k])], IsoTrack_dxyAbs)
                    if(abRejectedAsStandardMVA[1] and abRejectedAsLowPtBDT[k]):
                        Fill1D(histos['RejectBKasStandard_mvaNoIsoWP90&&LowPt_Pt_BDT'+str(bdtcuts[k])], pt)
                        if(bUseVertex):
                            Fill1D(histos['RejectBKasStandard_mvaNoIsoWP90&&LowPt_DxyAbs_BDT'+str(bdtcuts[k])], IsoTrack_dxyAbs)
                    if( abRejectedAsStandardMVA[2] and abRejectedAsLowPtBDT[k] ):
                        Fill1D(histos['RejectBKasStandard_mvaNoIsoWP80&&LowPt_Pt_BDT'+str(bdtcuts[k])], pt)
                        if(bUseVertex):
                            Fill1D(histos['RejectBKasStandard_mvaNoIsoWP80&&LowPt_DxyAbs_BDT'+str(bdtcuts[k])], IsoTrack_dxyAbs)

                    # All 3:

                    if(abRejectedAsStandardMVA[0] and abRejectedAsIPPMVA[0] and abRejectedAsLowPtBDT[k]):
                        Fill1D(histos['RejectBKasExtended_mvaNoIsoWPL_PhoMVA90&&LowPt_Pt_BDT'+str(bdtcuts[k])], pt)
                        if(bUseVertex):
                            Fill1D(histos['RejectBKasExtended_mvaNoIsoWPL_PhoMVA90&&LowPt_DxyAbs_BDT'+str(bdtcuts[k])], IsoTrack_dxyAbs)

                    if(abRejectedAsStandardCutBased[0] and abRejectedAsIPPCutBased[0] and abRejectedAsLowPtBDT[k]):
                        Fill1D(histos['RejectBKasExtended_CBVeto_PhoCBLoose&&LowPt_Pt_BDT'+str(bdtcuts[k])], pt)
                        if(bUseVertex):
                            Fill1D(histos['RejectBKasExtended_CBVeto_PhoCBLoose&&LowPt_DxyAbs_BDT'+str(bdtcuts[k])], IsoTrack_dxyAbs)
                    if(abRejectedAsStandardCutBased[1] and abRejectedAsIPPCutBased[1] and abRejectedAsLowPtBDT[k]):
                        Fill1D(histos['RejectBKasExtended_CBLoose_PhoCBMedium&&LowPt_Pt_BDT'+str(bdtcuts[k])], pt)
                        if(bUseVertex):
                            Fill1D(histos['RejectBKasExtended_CBLoose_PhoCBMedium&&LowPt_DxyAbs_BDT'+str(bdtcuts[k])], IsoTrack_dxyAbs)





    # Save histos:
    extrapolate_log.close()
    hfile.Write()


# Make samplechain
ch = SampleChain(sample, options.startfile, options.nfiles, year).getchain()

n_entries = ch.GetEntries() #num entries in tree
nevtcut = n_entries -1 if nEvents == - 1 else nEvents - 1


if __name__ == '__main__':
    print "RecoFixThreaded calculating efficiencies with and without recofix on",Nthreads,"threads. Sample: "+sample+", ptcut: "+str(ptcut)+((", ptcutHigh: "+str(ptcutHigh)) if options.ptcutHigh != -1 else "")
    proc = {}
    for i in range(Nthreads):
        it0 = 0 + i * (nevtcut + 1) / Nthreads 
        it1 = 0 + (i+1) * (nevtcut + 1) / Nthreads 
        proc[i] = multiprocessing.Process(target=RecoAndRecoFix, args=(i, it0, it1, nevtcut, ch, ptcut, ptcutHigh, binningSet))
        proc[i].start()
        #print 'Started thread',i
else:
    print("Something is serously wrong.")
    quit()

for i in range(Nthreads):
        proc[i].join()

print("All threads finished. Stacking results...")
# Stack partial results

highptcutstr = '_ptcutHigh'+str(ptcutHigh) if options.ptcutHigh != -1 else ''
extrapolstr = '_IPPextrapol' if bExtrapolate else ''

hfile = ROOT.TFile( 'root_files/RecoFixThreadedSTACK_Sample'+sample+'_year'+year+'_pt'+str(ptcut)+highptcutstr+extrapolstr+'.root', 'RECREATE')
savehistos = Histograms.SpawnHistograms(binningSet)


for threadID in range(Nthreads):
    savedfiles = ROOT.TFile.Open('root_files/RecoFixThreaded_Sample'+sample+'_Thread-'+str(threadID)+'.root')

    for key in savehistos.keys():
        #print "Trying for key: "+str(key)
        savehistos[key].Add( savedfiles.Get( key+'_' ))


print("Saving stacked root file...")
# Save endresult:
hfile.Write()


print("Removing temporary files...")
for i in range(Nthreads):
    os.remove('root_files/RecoFixThreaded_Sample'+sample+'_Thread-'+str(i)+'.root')

print("All finished! Sample: "+sample+", ptcut: "+str(ptcut)+((", ptcutHigh: "+str(ptcutHigh)) if options.ptcutHigh != -1 else ""))