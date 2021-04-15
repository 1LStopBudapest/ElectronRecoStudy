# ElectronRecoStudy
Will be more detailed later.

RecoFix.py: 
Try to recover lost electrons by matching isotracks and photons.

RecoFixThreaded.py: 
Identical to RecoFix.py, except uses multitasking.Process() to gain speed. Default hardware concurrency is specific to Higgs with 8 cores! Use --Nthreads X to change number of threads. 
While useful to gain execution speed, always doublecheck the results with the simple code variant, at least once.

RecoFixPlot.py: 
Plot histos made by RecoFixThreaded from root_files/. 
Use --threading 0 to use simple code results instead.

IDStudy.py: Study electron identification workingpoints and isolation-cut effects (with/without)

IDStudyPlot.py Plot histos made by IDStudy from root_files/

MatchCheck.py: study trying to understand the shape of some efficiency(eta) plots result in IDStudy

MatchCheckPlot.py: plot them
