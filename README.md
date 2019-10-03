# TDR Plots and Documentation

For PhaseII TDR documentation

# Generating Input Ntuples

Follow latest PhaseII CMSSW recipe
https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideL1TPhase2Instructions#Recipe

May need to modify the recipe to get the latest Hybrid tracking code:
https://gitlab.cern.ch/cms-tracker-phase2-backend-development/BE_software/L1Tracking/tree/master

Use latest version of track ntuples with jet collections here (on branch TDRBranch): 
https://gitlab.cern.ch/L1TrackJetsPrimaryVtx/TrackletL1PVJets/tree/TDRBranch


# TDR Plotting Codes

This area now includes plotting macros and codes to run over ntuples and produce the plots for the L1 TDR. 

The files are contained in FinalVersionMacros/

The Folder Prompt4parameter/JetEfficiency contains the plot recipe for prompt track jet efficiency for TTBar( NoPU) 
and the Stop sample (PU200). To generate the input run the Event loop in Prompt4parameter/JetEfficiency.C

The files contained in Prompt4parameter/TurnOnRate contain the plot recipe for the prompt track jet turn on for 
Stop sample (PU200). To generate the input for the macro run the event loop in Prompt4parameter/TurnOnRate/L1RateDisplacedPU.C

For the Hidden Higgs samples in the menu section of the TDR, we produce a L1 HT Rate vs. signal effiency plot as well 
as the signal yields vs. lifetime for both prompt and extended tracks input to track jets. The macro for the input of the signal efficiency 
and L1Rate (from the minbias sample) is the same: ExtendedTracking/L1RateDisplacedPU.C . It includes the output for both prompt HT and 
HT+displaced tag. The plotting macro for the L1 HT Rate vs. signal effiency plot is ExtendedTracking/QuickRateEff.C . 

