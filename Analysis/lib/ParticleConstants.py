import pandas as pd
import numpy as np

mZ = 91.1876 #GeV

pid = pd.Series(
    {'d': 1, 'u':2, 's':3, 'c':4, 'b':5, 't':6,
    'e':11, 've':12, 'mu':13, 'vmu':14,'tau':15, 'vtau':16,
    'gl':21,'ph':22,'Z':23,'Wp':24,'h0':25,'Hp':37,
    'N':999999}, name='ParticleID')
chargedlepid = [pid['e'],-pid['e'],pid['mu'],-pid['mu'],pid['tau'],-pid['tau']]
muid = [pid['mu'],-pid['mu']]
eid = [pid['e'],-pid['e']]
phoid = [pid['ph']]
invid = [12,-12,14,-14,16,-16]
status2id = [999999,-999999]
hnlid = [999999,-999999]

# Set column names for LHE events
ipid = 0; istatus = 1; im1 = 2; im2 = 3; id1 = 4; id2 = 5; ipx = 6; ipy = 7; ipz = 8; ie = 9; im = 10; iend1 = 11; iend2 = 12

pbTofb = 1000
clight = np.float64(299792458.) #m/s
OneEqualsXGeVs=np.float64(6.58*(10**(-25)))
