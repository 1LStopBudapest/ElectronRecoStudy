class isoTrackPhotonPair:
    eta = 0
    phi = 0
    etaSC = 0
    phiSC = 0
    pt = 0
    photon_cutBased = -1
    photon_cutBased_veto = False
    photon_cutBased_loose = False
    photon_cutBased_medium = False
    photon_cutBased_tight = False
    photon_mva80 = False
    photon_mva90 = False

    def __init__(self, ch, iTrack, iPhoton, bExtrapolate):
        self.eta = ch.IsoTrack_eta[iTrack]
        self.phi = ch.IsoTrack_phi[iTrack]
        if(bExtrapolate):
            self.etaSC = ch.IsoTrack_eta[iTrack] + ch.IsoTrack_deltaEta[iTrack]
            self.phiSC = ch.IsoTrack_phi[iTrack] + ch.IsoTrack_deltaPhi[iTrack]
        self.pt = ch.IsoTrack_pt[iTrack]

        if(ch.Photon_cutBased[iPhoton] >> 0 & 1):
            self.photon_cutBased_veto = True
            self.photon_cutBased = 0
        if(ch.Photon_cutBased[iPhoton] >> 1 & 1):
            self.photon_cutBased_loose = True
            self.photon_cutBased = 1
        if(ch.Photon_cutBased[iPhoton] >> 2 & 1):
            self.photon_cutBased_medium = True
            self.photon_cutBased = 2
        if(ch.Photon_cutBased[iPhoton] >> 3 & 1):
            self.photon_cutBased_tight = True
            self.photon_cutBased = 3
        


    def GetPt(self):
        return self.pt
    
    def GetEta(self):
        return self.eta
    
    def GetPhi(self):
        return self.phi

    def GetEtaSC(self):
        return self.etaSC

    def GetPhiSC(self):
        return self.phiSC

    def getPhotonCutBased(self):
        return self.photon_cutBased
    
    def getPhotonCutBasedBool(self, type):
        if(type == "veto"):
            return self.photon_cutBased >= 0
        elif(type == "loose"):
            return self.photon_cutBased >= 1
        elif(type == "medium"):
            return self.photon_cutBased >= 2
        elif(type == "tight"):
            return self.photon_cutBased >= 3
        return False
        
        



    
