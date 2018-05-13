import numpy as np
import pandas as pd
import ParticleConstants as pc
from scipy.stats import binned_statistic


def epsilon():
    return 0.0000000000001


#Assumes already ordered by pT
def LeadParticle(group,idArray):
    ParticleMask = group['pid'].isin(idArray)
    if np.any(ParticleMask):
        return group[ParticleMask].iloc[0]
def CountParticles(group,idArray):
    nParticle = np.count_nonzero(group['pid'].isin(idArray))
    return nParticle 

def IsGammaEInvMassCut(group,invMassCut):
    if(CountParticles(group,pc.eid) == 0):
        return True
    else:
        return np.abs(mass(LeadParticle(group,pc.eid)+LeadParticle(group,pc.phoid)) - pc.mZ) >= invMassCut  


class PJ:
    def __init__(self,e,px,py,pz):
        self.e = np.float64(e)
        self.px = np.float64(px)
        self.py = np.float64(py)
        self.pz = np.float64(pz)
        self.pvec = np.array([self.px,self.py,self.pz])
    def eta(self):
        return np.arctanh(self.pz/np.linalg.norm(self.pvec))
    def pt(self):
        return np.linalg.norm(np.array([self.px,self.py]))


def pvec(ev):
    return ev[['px','py','pz']]
def pnorm(ev):
    return np.sqrt(ev['px']**2+ev['py']**2+ev['pz']**2)
def eta(ev):
    return np.arctanh(ev['pz']/pnorm(ev))
def deltaEta(pj1,pj2):
    return eta(pj1)-eta(pj2)
def deltaR(pj1,pj2):
    return np.sqrt(deltaEta(pj1,pj2)**2+deltaPhi(pj1,pj2)**2)
def absEta(pj):
    return np.absolute(eta(pj))
def pt(pj):
    return np.sqrt(pj['px']**2+pj['py']**2)
def evVecMET(ev):
    visible = ~ev['pid'].isin(pc.invid) & ~ev['pid'].isin(pc.status2id)
    return -np.sum(ev[visible])[['px','py','pz']]
def evMET(ev):
    return pt(evVecMET(ev))


def phi(ev):
    phival = pt(ev)
    if(isinstance(phival,np.float64)):
        if( phival > 0. ):
            phival = np.arctan2(ev['py'],ev['px'])
        while( phival < 0. ):
            phival += 2*np.pi
        while( phival >= 2*np.pi):
            phival -= 2*np.pi
    if(isinstance(phival,pd.Series)):
        phival[phival > 0.] = np.arctan2(ev['py'],ev['px'])
        while( np.any(phival[phival < 0.]) ):
            phival[phival < 0.] += 2*np.pi
        while( np.any(phival[phival >= 2*np.pi]) ):
            phival[phival >= 2*np.pi] -= 2*np.pi
    return phival

def deltaPhi(ev1,ev2):
    dphi = phi(ev2)-phi(ev1)
    if(isinstance(dphi,np.float64)):
        while ( dphi > np.pi ): 
            dphi -= 2*np.pi;
        while ( dphi < -np.pi ): 
            dphi += 2*np.pi;

    if(isinstance(dphi,pd.Series)):
        while( np.any(dphi[dphi > np.pi])  ):
            dphi[dphi > np.pi] -= 2*np.pi
        while( np.any(dphi[dphi < -np.pi]) ):
            dphi[dphi < -np.pi] += 2*np.pi

    return dphi

def mass(ev):
    return np.sqrt(ev['e']**2 - pnorm(ev)** 2)
def MT(Lep4vec,METvec):
    a = 2*pt(METvec)*pt(Lep4vec)
    b = np.cos(deltaPhi(Lep4vec,METvec))
    c = np.sqrt(a*(1.-b))
    return c
def gamma(m,En):
    return En/m

def beta(m,En):
    g = gamma(m,En)
    return np.sqrt(g*g-1)/g



#def WolframGamma(a,x):
#    return scipy.special.gamma(a)*(1-scipy.special.gammainc(a,x))

#def CLEquation(Nsig, Nbg,alpha):
#    numerator = WolframGamma(1.+Nbg, Nbg + Nsig)
#    denom = WolframGamma(1. + Nbg, Nbg)
#    return alpha - (numerator / denom)

#def FindSensitivity(Nbg, alpha):        
#    TargetSignal = minimize(CLEquation, [2.], args=(Nbg,alpha), tol=1e-10, method='Nelder-Mead')    
#    return TargetSignal

def HNLGammaDecay(d, m):
    return (d**2)*(m**3)/(4*np.pi)  

def ProbDecay(d, m, En, L1, L2 ):
    conv = pc.OneEqualsXGeVs*pc.clight
    LDecay = conv*gamma(m, En) * beta(m, En) / HNLGammaDecay(d,m)
    return (np.exp(-L1/LDecay) - np.exp(-L2/LDecay))

def RebinMeanWeights(dat, nBins, Range):
    Means = binned_statistic(dat, dat, statistic='mean',bins=nBins, range=Range)[0]
    Weights = binned_statistic(dat, dat, statistic='count',bins=nBins, range=Range)[0].astype(np.float32)
    nonNan = ~np.isnan(Means)
    Means = Means[nonNan]
    Weights = Weights[nonNan]
    Weights_norm = Weights / np.sum(Weights)
    assert len(Means) == len(Weights_norm)
    return Means, Weights_norm

def averagedProbDecay(d, m, En, weights, L1, L2):
    return np.sum( weights * ProbDecay(d, m, En, L1, L2) ) 

def SensitivityEquation(d, SigmaEpsilonA, NCLSignalEvents, m, En, weights, L1, L2):
    averagedprobdecay = averagedProbDecay(d, m, En, weights, L1, L2) 
    return (xSquared(d)*SigmaEpsilonA*averagedprobdecay - NCLSignalEvents)**2
def xSquared(x):
    return x**2;
