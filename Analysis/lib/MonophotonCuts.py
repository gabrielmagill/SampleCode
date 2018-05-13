import numpy as np

##############################################################
##############################################################
#################  MONOPHOTON ANALYSIS!!!! ###################
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
etaGammaBad1_lb = 1.37
etaGammaBad1_ub = 1.52
etaGammaBad2_lb = 2.37
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
nLepton_lb = 0
nLepton_up = 0

nPhotonFlag = True
nGamma_lb = 1
nGamma_up = 1000

ETLeadGammaFlag = True
ETLeadGamma_lb = 150./scale

METFlag = True
MET_lb = 150./scale

## MonoPhoton Only ##
DeltaPhiGammaMETFlag = True
DeltaPhiGammaMET_lb = 0.4/scale2


## LeptonPhoton Only ##
DeltaRLeadGammaLepFlag = False
DeltaRLeadGammaLep_lb = 0.8/scale

GammaEInvMassFlag = False
GammaEInv_lb = 10.

MTFlag = False
MT_lb = 100./scale

########## Sensitivity Variables ##########

effElectron = 0.8
effMuon = 0.9

Luminosity = 36.1

rtransElectron = 0.2/1000. #m
rtransMuon = 2./1000. #m
zlongElectron = 1./1000. #m
zlongMuon = 5./1000. #m
zlongLeadGamma = 0.25 #m

L1gamma = 0.
L2gamma = 1.5 #m

nSR = 5
SRMETFlag = True
SRMET_lb = np.array([150,   225,    300,    150,    225])
SRMET_up = np.array([140000,140000,140000,  225,    300])
effGamma = np.array([0.92,  0.92,   0.92,   0.92,   0.92])
NCLSignal = np.array([96.4397, 42.6703, 25.8965, 100.601, 37.4494])
# 95%CL NCLSignal were obtained using the Sensitivity.nb Notebook


####### Plotting and Output ######

# Masses and indices to analyse
dScale = np.array([1]) # dProd = d, dDecay = d * dscale

massLoop = np.array([0.001,0.003,0.01,0.03,0.1,0.3,1,3,10,30,100,300,1000,1500,2000,2500,3000,4000])
maxIndex ={"0.001":5,
           "0.003":5,
           "0.01":5,
           "0.03":2,
           "0.1":2,
           "0.3":1,
           "1":2,
           "3":1,
           "10":1,
           "30":1,
           "100":1,
           "300":1,
           "1000":1,
           "1500":1,
           "2000":1,
           "2500":1,
           "3000":1,
           "4000":1}


#massLoop = np.array([4000])
#maxIndex ={"4000":1}


# Masses to save plots to disk
massPlot = massLoop
