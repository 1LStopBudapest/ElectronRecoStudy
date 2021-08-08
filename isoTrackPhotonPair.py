class isoTrackPhotonPair:
    eta = 0
    phi = 0
    etaSC = 0
    phiSC = 0
    pt = 0

    def __init__(self, ch, iTrack, bExtrapolate):
        self.eta = ch.IsoTrack_eta[iTrack]
        self.phi = ch.IsoTrack_phi[iTrack]
        if(bExtrapolate):
            self.etaSC = ch.IsoTrack_eta[iTrack] + ch.IsoTrack_deltaEta[iTrack]
            self.phiSC = ch.IsoTrack_phi[iTrack] + ch.IsoTrack_deltaPhi[iTrack]
        self.pt = ch.IsoTrack_pt[iTrack]