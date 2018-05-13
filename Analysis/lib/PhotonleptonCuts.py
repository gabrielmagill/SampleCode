import numpy as np

##############################################################
##############################################################
#################  PHOTONLEPTON ANALYSIS!!!! ###################
##############################################################
##############################################################

# Notes:
# Everything grouped together with a flag is set by the flag

#set this to 1. at the end, it's used to relax constraints when developing code for small samples
scale = 1.
scale2= 1.

########## Reading Information #############
#max number of particles that on LHE event can have
maxParticlesPerEventLHE = 5
nEventsLHE = 50000

########## Loose Cuts ##########
#### Picks out individual particles from events ###
#Photon
ObjectIdPhotonFlag = True
etaGammaBad1_lb = 1.44
etaGammaBad1_ub = 10.0
etaGammaBad2_lb = 1.44
ETSubGammaBad_ub = 10.

#Lepton 
ObjectIdLeptonFlag = True
pTLeptonBad_ub = 25./scale
etaMuonBad_lb = 2.4*scale
etaElectronBad1_lb = 1.44
etaElectronBad1_ub = 1.56
etaElectronBad2_lb = 2.5*scale

########## Hard (Events Selection) Cuts ##########
#### Vetoes events from run ####

## Both
nLeptonFlag = True
nLepton_lb = 1
nLepton_up = 1

nPhotonFlag = True
nGamma_lb = 1
nGamma_up = 1000

ETLeadGammaFlag = True
ETLeadGamma_lb = 40./scale

METFlag = True
MET_lb = 120./scale

## MonoPhoton Only ##
DeltaPhiGammaMETFlag = False
DeltaPhiGammaMET_lb = 0.4/scale2

## Leptonphoton Only ##
DeltaRLeadGammaLepFlag = True
DeltaRLeadGammaLep_lb = 0.8/scale

GammaEInvMassFlag = True
GammaEInv_lb = 10.

MTFlag = True
MT_lb = 100./scale

########## Sensitivity Variables ##########

effElectron = 0.8
effMuon = 0.9

Luminosity = 19.7

rtransElectron = 0.2/1000. #m
rtransMuon = 2./1000. #m
zlongElectron = 1./1000. #m
zlongMuon = 5./1000. #m
zlongLeadGamma = 0.25 #m

L1gamma = 0.
L2gamma = 1.29 #m

nSR = 3
SRMETFlag = True
SRMET_lb = np.array([120.,200.,300.])
SRMET_up = np.array([200.,300.,140000.])
effGamma = np.array([0.9, 0.9, 0.9])
NCLSignal = np.array([19.7118, 8.16854, 5.07888])
# 95%CL NCLSignal were obtained using the Sensitivity.nb Notebook

####### Plotting and Output ######

dScale = np.array([0.2,1.,5.]) # dProd = d, dDecay = d * dscale

# Masses and indices to analyse
massLoop = np.array([0.001,0.003,0.01,0.03,0.1,0.3,1,3,10,30,100,300,1000,1500,2000,2500,3000])
maxIndex ={"0.001":5,
           "0.003":5,
           "0.01":5,
           "0.03":5,
           "0.1":5,
           "0.3":5,
           "1":5,
           "3":5,
           "10":5,
           "30":5,
           "100":1,
           "300":1,
           "1000":1,
           "1500":1,
           "2000":1,
           "2500":1,
           "3000":1}

#massLoop = np.array([2000, 2500, 3000])
#maxIndex={"2000":1,"2500":1,"3000":1}


# Masses to save plots to disk
massPlot = massLoop
