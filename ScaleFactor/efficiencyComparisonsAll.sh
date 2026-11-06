#!/bin/bash


python3 efficiencyCompare.py JPsi/UL2018miniAODv2/PUincl/tnpCombinedEle_StdPref_afIdCo_noisoCBLoose_IDs_HybridIso/passingBDT2.5/results_JPsi_chev3.yaml UL2018miniAODv2/tnpCombinedEle_StdPref_afIdCo_noisoCBLoose_IDs_HybridIso/passingBDT2.5/results_Z_chev3.yaml --fr1=False --fr2=False --fname1=JPsi-Combined --fname2=Z18-Combined-WindowFit

# JPsi fr is always false, it's just variable window minuit fit!
python3 efficiencyCompare.py JPsi/UL2018miniAODv2/PUincl/tnpCombinedEle_StdPref_afIdCo_noisoCBLoose_IDs_HybridIso/passingBDT2.5/results_JPsi_chev3.yaml UL2018miniAODv2/tnpCombinedEle_StdPref_afIdCo_noisoCBLoose_IDs_HybridIso/passingBDT2.5/results_Z_chev3.yaml --fr1=False --fr2=True --fname1=JPsi-Combined --fname2=Z18-Combined-FullRangeFit



python3 efficiencyCompare.py JPsi/UL2018miniAODv2/PUincl/RemLowPtEle_tnpCombinedEle_StdPref_afIdCo_noisoCBLoose_IDs_HybridIso/passingBDT2.5/results_JPsi_chev3.yaml UL2018miniAODv2/RemLowPt_HybridIso/passingBDT2.5/results_Z_chev3.yaml --fr1=False --fr2=False --fname1=JPsi-RemLowPt --fname2=Z18-RemLowPt-WindowFit


python3 efficiencyCompare.py JPsi/UL2018miniAODv2/PUincl/RemLowPtEle_tnpCombinedEle_StdPref_afIdCo_noisoCBLoose_IDs_HybridIso/passingBDT2.5/results_JPsi_chev3.yaml UL2018miniAODv2/RemLowPt_HybridIso/passingBDT2.5/results_Z_chev3.yaml --fr1=False --fr2=True --fname1=JPsi-RemLowPt --fname2=Z18-RemLowPt-FullRangeFit



python3 efficiencyCompare.py JPsi/UL2018miniAODv2/PUincl/tnpEleIDs_HybridIso/passingLoose/results_JPsi_chev3.yaml UL2018miniAODv2/tnpEleIDs_HybridIso/Loose/results_Z_chev3.yaml --fr1=False --fr2=False --fname1=JPsi-Std --fname2=Z18-Std-WindowFit

python3 efficiencyCompare.py JPsi/UL2018miniAODv2/PUincl/tnpEleIDs_HybridIso/passingLoose/results_JPsi_chev3.yaml UL2018miniAODv2/tnpEleIDs_HybridIso/Loose/results_Z_chev3.yaml --fr1=False --fr2=True --fname1=JPsi-Std --fname2=Z18-Std-FullRangeFit


# Jpsi vs Jpsi

python3 efficiencyCompare.py JPsi/UL2018miniAODv2/PUincl/tnpCombinedEle_StdPref_afIdCo_noisoCBLoose_IDs_HybridIso/passingBDT2.5/results_JPsi_chev3.yaml JPsi/UL2018miniAODv2/PUincl/RemLowPtEle_tnpCombinedEle_StdPref_afIdCo_noisoCBLoose_IDs_HybridIso/passingBDT2.5/results_JPsi_chev3.yaml --fr1=False --fr2=False --fname1=JPsi-Combined --fname2=JPsi-RemLowPt

python3 efficiencyCompare.py JPsi/UL2018miniAODv2/PUincl/tnpCombinedEle_StdPref_afIdCo_noisoCBLoose_IDs_HybridIso/passingBDT2.5/results_JPsi_chev3.yaml JPsi/UL2018miniAODv2/PUincl/tnpEleIDs_HybridIso/passingLoose/results_JPsi_chev3.yaml --fr1=False --fr2=False --fname1=JPsi-Combined --fname2=JPsi-Std



# Z vs Z

python3 efficiencyCompare.py UL2018miniAODv2/tnpCombinedEle_StdPref_afIdCo_noisoCBLoose_IDs_HybridIso/passingBDT2.5/results_Z_chev3.yaml UL2018miniAODv2/RemLowPt_HybridIso/passingBDT2.5/results_Z_chev3.yaml  --fr1=False --fr2=False --fname1=Z18-Combined-WindowFit --fname2=Z18-RemLowPt-WindowFit
python3 efficiencyCompare.py UL2018miniAODv2/tnpCombinedEle_StdPref_afIdCo_noisoCBLoose_IDs_HybridIso/passingBDT2.5/results_Z_chev3.yaml UL2018miniAODv2/RemLowPt_HybridIso/passingBDT2.5/results_Z_chev3.yaml  --fr1=True --fr2=True --fname1=Z18-Combined-FullRangeFit --fname2=Z18-RemLowPt-FullRangeFit


python3 efficiencyCompare.py UL2018miniAODv2/tnpCombinedEle_StdPref_afIdCo_noisoCBLoose_IDs_HybridIso/passingBDT2.5/results_Z_chev3.yaml UL2018miniAODv2/tnpEleIDs_HybridIso/Loose/results_Z_chev3.yaml  --fr1=False --fr2=False --fname1=Z18-Combined-WindowFit --fname2=Z18-Std-WindowFit
python3 efficiencyCompare.py UL2018miniAODv2/tnpCombinedEle_StdPref_afIdCo_noisoCBLoose_IDs_HybridIso/passingBDT2.5/results_Z_chev3.yaml UL2018miniAODv2/tnpEleIDs_HybridIso/Loose/results_Z_chev3.yaml  --fr1=True --fr2=True --fname1=Z18-Combined-FullRangeFit --fname2=Z18-Std-FullRangeFit