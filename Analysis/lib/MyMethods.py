import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

kwargs = dict(histtype='stepfilled',alpha=0.3, normed=True ,edgecolor='none',fill=True)
col4 = ['black','blue','red','green']
col3 = ['black','blue','red']
col2 = ['black','blue']
col1 = ['black']
col=[col1,col2,col3,col4]


def MyHistogram(Data,Range,nBins,PlotLabel,Legends,name):

    nPlots = len(Data)    
    
    for i in range(0,nPlots):
        plt.subplots_adjust(hspace=.4)
        plt.subplot(str(nPlots)+str(1)+str(i+1))
        plt.hist(Data[i], label=Legends[i], range=Range, bins=nBins,   color=col[len(Data[i])-1], **kwargs)
        plt.xlabel(PlotLabel[0])
        plt.ylabel(PlotLabel[1])
        plt.legend(loc='best') 
    plt.savefig(name)
    plt.close()
    
def WriteListToFile(file,row):
    for r in row:
        file.write(str(r))
        if r == row[-1]:
            file.write("\n")
        else:
            file.write(",")
