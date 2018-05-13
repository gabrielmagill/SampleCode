import numpy as np
import sys
sys.path.append('lib/')
import pandas as pd
import ParticlePhysicsLibrary as ppl
import ParticleConstants as pc
import matplotlib.pyplot as plt
import MyMethods as MyMeth
import time
from scipy.optimize import minimize
plt.style.use('seaborn-white')


HNLRepoDir = str(sys.argv[1])

outputColumns = ['HNLmass','dSR1','dSR2','dSR3','dSR4','dSR5','dSRBest','SRBestNumber','meanMasstruth','sigmaMCdequal1','sigmaEpsilonA','sigmaEpsilonAProbDecay','effMCToEventSelection','effEventSelectionToBestSR','NEventsPostSRBest','ProbDecaySRBest','NEventsdSRBest','NEventsdSRBestoverNEvents95CL','dScale','avgSRBestDeltaRPhoHNL']

# Write Columns to file
dt = open("outputColumns.csv",'w')
# Write data to file
MyMeth.WriteListToFile(dt,outputColumns)
dt.close()



if(len(sys.argv)!=3):
	print("\nRun as: python3 <program.py> <Directory> <StartOfFileName> ")
	print("Make sure all header files above in same directory\n")
	sys.exit()	


if("HNLMonoPhoton.13TeV.UFO." in str(sys.argv[2])):
	import MonophotonCuts as cuts
	AnalName = "HNLMonoPhoton.13TeV."

if("HNLPhotonLepton.8TeV.UFO." in str(sys.argv[2])):
	import PhotonleptonCuts as cuts
	AnalName = "HNLPhotonLepton.8TeV."


def DisplayEvents(hnlmass):
	
	# All this assumes 1 lepton, >=1 photons
	PhotonMask = events['pid']==pc.pid['ph']
	ManyPhotonMask = PhotonMask & events.groupby('evid')['pid'].transform(lambda x: (np.count_nonzero(x==pc.pid['ph']) >= 2)).astype(np.bool)
	
	LeadPhotons = events.groupby('evid').apply(lambda x: ppl.LeadParticle(x,pc.phoid))
	SubLeadPhotons = events[ManyPhotonMask].groupby('evid').apply(lambda x: x.iloc[1])
	LeadHNL = events.groupby('evid').apply(lambda x: ppl.LeadParticle(x,pc.hnlid))
	LeadMET = events.groupby('evid').apply(lambda x: ppl.evVecMET(x))
	
	#Make Plots
	# myHist.MyHistogram(pData,pRange,pBins,pSubTitle,pLabel,pLegends,pName)
	prename = '../Plots/'+AnalName+'Mass'+str(hnlmass)+'GeV.'

	# pT Plot
	if(len(SubLeadPhotons) >= 1):
		pData = [[LeadPhotons['pt'], SubLeadPhotons['pt']], [ppl.pt(LeadMET)]]
		pLegends = [["Lead Photon", "SubLead Photon"], ["MET"]]
		MyMeth.MyHistogram(pData,[0,800],80,["pT [GeV]","Counts"],pLegends,prename+'PhoMETpT.pdf')		
	else:
		pData = [[LeadPhotons['pt']], [ppl.pt(LeadMET)]]
		pLegends = [["Lead Photon"], ["MET"]]
		MyMeth.MyHistogram(pData,[0,800],80,["pT [GeV]","Counts"],pLegends,prename+'PhoMETpT.pdf')		

	# eta Plots
	if(len(SubLeadPhotons) > 0):
		pData = [[LeadPhotons['eta'], SubLeadPhotons['eta']]]
		pLegends = [["Lead Photon", "SubLead Photon"]]
		MyMeth.MyHistogram(pData,[-3,3],80,["eta","Counts"],pLegends,prename+'PhotonEta.pdf')		
	else:
		pData = [[LeadPhotons['eta']]]
		pLegends = [["Lead Photon"]]
		MyMeth.MyHistogram(pData,[-3,3],80,["eta","Counts"],pLegends,prename+'PhotonEta.pdf')		

	# HNL Energy Plot
	pData = [[LeadHNL['e']]]
	pLegends = [["HNL"]]
	MyMeth.MyHistogram(pData,[0,2000],80,["Energy [GeV]","Counts"],pLegends,prename+'HNLEnergy.pdf')		

	# HNL Mass Plot
	pData = [[ppl.mass(LeadHNL)]]
	pLegends = [["HNL Theoretical mass: "+str(hnlmass)+"GeV"]]
	MyMeth.MyHistogram(pData,[0,2000],50,["Mass [GeV]","Counts"],pLegends,prename+'HNLMass.pdf')		

	# deltaPhi(LeadPhoton,MET)		
	pData = [[ppl.deltaPhi(LeadPhotons,LeadMET)]]
	pLegends = [["deltaPhi(LeadPhoton,MET)"]]
	MyMeth.MyHistogram(pData,[-2*np.pi,2*np.pi],80,["deltaPhi(LeadPhoton,MET)","Counts"],pLegends,prename+'deltaPhiLeadPhotonMET.pdf')		

	# deltaR(LeadPhoton, HNL)
	pData = [[ppl.deltaR(LeadPhotons,LeadHNL)]]
	pLegends = [["deltaR(LeadPhoton,LeadHNL)"]]
	MyMeth.MyHistogram(pData,[0,2],80,["deltaR(LeadPhoton,LeadHNL)","Counts"],pLegends,prename+'deltaRLeadPhotonLeadHNL.pdf')		


	if(cuts.nLepton_lb >= 1):
		LeadLeptons = events.groupby('evid').apply(lambda x: ppl.LeadParticle(x, pc.chargedlepid 	))
		LeadMuons =   events.groupby('evid').apply(lambda x: ppl.LeadParticle(x, pc.muid			))
		LeadElectrons = events.groupby('evid').apply(lambda x: ppl.LeadParticle(x,pc.eid))

		# pT Plot
		pData = [[LeadLeptons['pt'], LeadMuons['pt'], LeadElectrons['pt']]]
		pLegends = [["Lead Lepton", "Lead Muon", "Lead Electron"]]
		MyMeth.MyHistogram(pData,[0,800],80,["pT [GeV]","Counts"],pLegends,prename+'LeptonpT.pdf')		

		#eta
		pData = [[LeadLeptons['eta'], LeadMuons['eta'], LeadElectrons['eta']]]
		pLegends = [["Lead Lepton", "Lead Muon", "Lead Electron"]]
		MyMeth.MyHistogram(pData,[-2.7,2.7],80,["eta","Counts"],pLegends,prename+'Leptoneta.pdf')		

	
		if(cuts.DeltaRLeadGammaLepFlag):
			pData = [[ppl.deltaR(LeadPhotons,LeadLeptons)]]
			pLegends = [["deltaR(LeadPhoton,LeadLepton)"]]
			MyMeth.MyHistogram(pData,[0,10],80,["deltaR(LeadPhoton,LeadLepton)","Counts"],pLegends,prename+'deltaRLeadPhotonLeadLepton.pdf')		
		if(cuts.MTFlag):
			pData = [[ppl.MT(LeadLeptons,LeadMET)]]
			pLegends = [["MT(LeadLepton,MET)"]]
			MyMeth.MyHistogram(pData,[0,1000],80,["MT(LeadLepton,MET)","Counts"],pLegends,prename+'MTLeadLeptonMET.pdf')		
			

massIndex = -1 

for lmass in cuts.massLoop:

	start_time = time.time()
	evId = -1
	partId = -1
	massIndex += 1
	
	#Used to fix a naming bug
	lmassStr = lmass
	if(lmass >=1):
		lmassStr = lmass.astype(np.int)

	print("Reading lmassStr=",lmassStr,"GeV. . .")
	
	# Record Cross Section
	with open(HNLRepoDir+"CrossSection."+str(sys.argv[2]) +str(lmassStr)+"GeV.1_unweighted_events.lhe") as CrossSection:
		sigma = pc.pbTofb*np.float64([line.split() for line in CrossSection])[0,0]
	
	#Setup container for all events
	NPanda = (cuts.nEventsLHE)*(cuts.maxParticlesPerEventLHE)*(cuts.maxIndex[str(lmassStr)])
	events = pd.DataFrame(np.full((NPanda,7), np.nan),columns=['evid','pid','im1','e','px','py','pz']) 
	
	for ind in range(1,1+cuts.maxIndex[str(lmassStr)]):
		fileName = HNLRepoDir+str(sys.argv[2]) +str(lmassStr)+"GeV."+str(ind)+"_unweighted_events.lhe"
		recording = False
		with open('%s' % fileName, 'r') as f:

			for line in f:
				if line == "<event>\n":
					recording = True
					evId += 1
					continue

				if recording == True:
					splitline = line.split()
					#find final state particles
					if len(splitline)==13:
						idtemp = np.int32(splitline[pc.ipid])
						if (splitline[pc.istatus] == '1') | (np.abs(idtemp)==pc.pid['N']):
							partId += 1
							events.iloc[partId] = {
								'evid':np.int64(evId),
								'pid':idtemp,
								'im1':np.int16(splitline[pc.im1]),
								'e':np.float64(splitline[pc.ie]),
								'px':np.float64(splitline[pc.ipx]),
								'py':np.float64(splitline[pc.ipy]), 
								'pz':np.float64(splitline[pc.ipz])
								}

				#End of event
				if line == "<mgrwt>\n":
					recording = False
	
	end_time = time.time()
	events = events.dropna()
	nEventsInitial = np.float64(evId+1.)

	print("Comment: Read m =",lmass,"GeV file in",np.round(end_time-start_time), "seconds")
	
	# File is done being read
	# Main analysis Code Here:

	#Add a pT and eta column to DataFrame
	events['pt'] = ppl.pt(events)
	events['eta'] = ppl.eta(events)

	#Sort events by pT
	events = events.sort_values(by=['evid','pt'],ascending=[1,0])

	#####   Loose Cuts   #####

	# Removes particles from events that don't satisfy requirements
	if(cuts.ObjectIdPhotonFlag):

		#All photon pT cuts
		PhotonMask = (events['pid']==pc.pid['ph'])
		cut1 = (events['pt'] < cuts.ETSubGammaBad_ub) & PhotonMask
		
		#All photon eta cuts
		cut2a =(
			(  cuts.etaGammaBad1_lb < np.abs(events['eta'])  ) & 
			(  np.abs(events['eta']) < cuts.etaGammaBad1_ub )
			)			
		cut2b = cuts.etaGammaBad2_lb < np.abs(events['eta'])
		cut2 = (cut2a | cut2b) & PhotonMask
			
			#Filter bad events
		events = events[~cut1 & ~cut2 ]		


	# Removes leptons from events that don't satisfy requirements
	if(cuts.ObjectIdLeptonFlag):		

		#Setup useful functions
		LeptonMask = (events['pid'].isin(pc.chargedlepid))
		ElectronMask = (events['pid'].isin(pc.eid))
		MuonMask = (events['pid'].isin(pc.muid))

		#Lepton pT cuts
		cut1 =  (events['pt'] < cuts.pTLeptonBad_ub) & LeptonMask
		
		# Eta Cuts
		cut2 = (np.abs(events['eta']) > cuts.etaMuonBad_lb) & MuonMask
		cut3a = (
			(	cuts.etaElectronBad1_lb < np.abs(events['eta'])	) & 
			(	np.abs(events['eta']) < cuts.etaElectronBad1_ub	)
			)
		cut3b = np.abs(events['eta']) > cuts.etaElectronBad2_lb
		cut3 = (cut3a | cut3b) & ElectronMask

		#Filter bad events
		events = events[ ~cut1 & ~cut2 & ~cut3]
	
	print("Comment: Loose Cuts Complete")


	####### Apply Event Selection Cuts #######
	# Number of Photons
	if(cuts.nPhotonFlag):
		events = events.groupby('evid').filter(lambda x: (ppl.CountParticles(x,pc.phoid)>=cuts.nGamma_lb) &
											(ppl.CountParticles(x,pc.phoid) <= cuts.nGamma_up) )
	# Number of Leptons
	if(cuts.nLeptonFlag):
		events = events.groupby('evid').filter(lambda x: (ppl.CountParticles(x,pc.chargedlepid)>=cuts.nLepton_lb) &
											(ppl.CountParticles(x,pc.chargedlepid) <= cuts.nLepton_up) )
	#ET Gamma Cut
	if(cuts.ETLeadGammaFlag):
		events = events.groupby('evid').filter(lambda x: ppl.LeadParticle(x,pc.phoid)['pt']>=cuts.ETLeadGamma_lb)	
	#Missing Energy Cuts
	if(cuts.METFlag):
		events = events.groupby('evid').filter(lambda x: ppl.evMET(x) >= cuts.MET_lb)		
	# Angle between MET and Lead Gamma	
	if(cuts.DeltaPhiGammaMETFlag):
		events = events.groupby('evid').filter(lambda x: np.abs(ppl.deltaPhi(ppl.LeadParticle(x,pc.phoid),ppl.evVecMET(x))) >= cuts.DeltaPhiGammaMET_lb)	
	# DeltaR(Gamma,ell)
	if(cuts.DeltaRLeadGammaLepFlag):
		events = events.groupby('evid').filter(lambda x: np.abs(ppl.deltaR(ppl.LeadParticle(x,pc.phoid),ppl.LeadParticle(x,pc.chargedlepid))) >= cuts.DeltaRLeadGammaLep_lb)	
	# MT = sqrt(2*MET*pT(1-cosDeltaPhi))
	if(cuts.MTFlag):
		events = events.groupby('evid').filter(lambda x: ppl.MT(ppl.LeadParticle(x,pc.chargedlepid),ppl.evVecMET(x)) >= cuts.MT_lb)	
	# |m(e,gamma)-m_Z| > lb
	if(cuts.GammaEInvMassFlag):
		events = events.groupby('evid').filter(lambda x: ppl.IsGammaEInvMassCut(x,cuts.GammaEInv_lb))	

	print("Comment: Event Selection Cuts Complete")
	
	#### Visualize the Event ####
	if(lmass in cuts.massPlot):
		DisplayEvents(lmassStr)

	
	for dscale in cuts.dScale:

		##### Do sensitivity #####
		bounds = np.full(10, np.nan)
		effSR = np.full(cuts.nSR, np.nan)
		nEventsPostSelection = len(np.unique(events['evid']))
		effEvent = nEventsPostSelection / nEventsInitial
		effLepton = np.full(cuts.nSR, 1.)
		ProbDecaySR = np.full(cuts.nSR,np.nan)
		NEventsdSol = np.full(cuts.nSR,np.nan)
		nEventsPostSR = np.full(cuts.nSR,np.nan)
		meanMassTruth = np.full(cuts.nSR,np.nan)
		meanDeltaRPhoHNL = np.full(cuts.nSR,np.nan)
	
		for ii in range(0,cuts.nSR):

			if(cuts.SRMETFlag):
				eventsSR = events.groupby('evid').filter(lambda x: (ppl.evMET(x) >= cuts.SRMET_lb[ii]) & (ppl.evMET(x) <= cuts.SRMET_up[ii]))
			LeadHNL = eventsSR.groupby('evid').apply(lambda x: ppl.LeadParticle(x,pc.hnlid))
			PhotonFromHNL = eventsSR[eventsSR['im1']>2.].groupby('evid').apply(lambda x: ppl.LeadParticle(x,pc.phoid))
			nEventsPostSR[ii] = len(np.unique(eventsSR['evid']))
			effSR[ii] = nEventsPostSR[ii] / nEventsPostSelection  
	
			meanMassTruth[ii] = np.mean(ppl.mass(LeadHNL))
		
			if(cuts.nLepton_lb >= 1):
				nElEv = np.sum(eventsSR.groupby('evid').apply(lambda x: ppl.CountParticles(x, pc.eid) > 0 ))
				nMuEv = np.sum(eventsSR.groupby('evid').apply(lambda x: ppl.CountParticles(x,pc.muid) > 0 ))
				assert (nElEv + nMuEv == nEventsPostSR[ii])
				effLepton[ii] = (cuts.effElectron*nElEv + cuts.effMuon*nMuEv)/nEventsPostSR[ii]

			EnMeans, EnWeights = ppl.RebinMeanWeights(LeadHNL['e'], 10, (np.min(LeadHNL['e'])-0.001, np.max(LeadHNL['e'])+0.001))
			SigmaEpsilonA = sigma * effSR[ii] * effEvent * cuts.effGamma[ii] * effLepton[ii] * cuts.Luminosity #d and ProbDecay are in SensitivityEquation
			guess = 1.e-5
			bounds[ii] = np.abs(minimize(ppl.SensitivityEquation, [guess], 
				args=(SigmaEpsilonA, cuts.NCLSignal[ii], lmass, EnMeans, EnWeights, cuts.L1gamma*(dscale**2), cuts.L2gamma*(dscale**2)), 
				tol=1e-10, method='Nelder-Mead').x[0])

		
			ProbDecaySR[ii] = ppl.averagedProbDecay(bounds[ii], lmass, EnMeans, EnWeights, cuts.L1gamma*(dscale**2), cuts.L2gamma*(dscale**2))
			NEventsdSol[ii] = SigmaEpsilonA * ProbDecaySR[ii] * (bounds[ii]**2)
		
			meanDeltaRPhoHNL[ii] = np.mean(ppl.deltaR(PhotonFromHNL,LeadHNL).dropna())
	
		argBest = np.argmin(bounds[:cuts.nSR])	
		Output = pd.DataFrame(columns=outputColumns)
		Output = Output.append({
			'HNLmass':lmass,
			'dSR1':bounds[0],
			'dSR2':bounds[1],
			'dSR3':bounds[2],
			'dSR4':bounds[3],
			'dSR5':bounds[4],
			'dSRBest':bounds[argBest],
			'SRBestNumber':argBest+1,
			'meanMasstruth':meanMassTruth[argBest],
			'sigmaMCdequal1':sigma,
			'sigmaEpsilonA':NEventsdSol[argBest]/(cuts.Luminosity*ProbDecaySR[argBest]),
			'sigmaEpsilonAProbDecay':NEventsdSol[argBest]/(cuts.Luminosity),
			'effMCToEventSelection':effEvent,
			'effEventSelectionToBestSR':effSR[argBest],
			'NEventsPostSRBest':nEventsPostSR[argBest],
			'ProbDecaySRBest':ProbDecaySR[argBest],
			'NEventsdSRBest':NEventsdSol[argBest],
			'NEventsdSRBestoverNEvents95CL': NEventsdSol[argBest]/cuts.NCLSignal[argBest],
			'dScale':dscale,
			'avgSRBestDeltaRPhoHNL':meanDeltaRPhoHNL[argBest]},ignore_index=True)

		print("Comment: Final Sensitivity:")
		print(np.array2string(Output.values[0], separator=', '))
		print("\n\n")
	
		# Write to file
		out = open("LHCSensitivity."+AnalName+"dscale"+str(dscale)+".csv",'a')
		MyMeth.WriteListToFile(out, Output.values[0])
		out.close()

print("### Don't forgot to copy paste output file into the Sensitivity folder! ###")
