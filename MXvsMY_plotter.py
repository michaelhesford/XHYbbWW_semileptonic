from argparse import Namespace
from glob import glob
from XHYbbWW_studies import XHYbbWW_studies
from TwoDHistStyle import style
from TIMBER.Tools.Common import DictStructureCopy, CompileCpp, ExecuteCmd, OpenJSON, StitchQCD
from TIMBER.Tools.Plot import CompareShapes, EasyPlots
from TIMBER.Analyzer import HistGroup
from collections import OrderedDict
import ROOT

# maybe think about adding args later to only do CERTAIN mass points (MX, MY), for now just do all
def plot():
    # we want to plot mX vs mY for the two score regions for QCD, V+Jets, ttbar, and Signal    
    #now do signal
    files = [f for f in glob('plots/MXvsMY_lepVetoSRorCR/*MXvsMY_lepVetoSRorCR.root')] #loop through all files
    #files = [f for f in glob('XHY*_MXvsMYstudies.root')]
    c = ROOT.TCanvas('c','c')
    c.cd()
    for f in files:
        inFile=ROOT.TFile.Open(f,'READ')
        setname=f.split('/')[2].split('_')[0]
        #setname = f.split('_')[0]
        histlist=[] #for 2D plots
        for key in inFile.GetListOfKeys():
            varname=key.GetName()
            inhist = inFile.Get(varname)  # get it from the file
            inhist.SetDirectory(0) # set the directory so hist is stored in memory and not as reference to TFile (this way it doesn't get tossed by python garbage collection when infile changes)      
            inhist.SetName('{}_'.format(setname)+varname)
            c.Clear
            style(inhist,xTitle='m_{X} (GeV)',yTitle='m_{Y} (GeV)',histTitle=inhist.GetTitle())
            inhist.Draw('COLZ')
            c.Print('plots/MXvsMY_lepVetoSRorCR/pdf/{}.pdf'.format(inhist.GetName()))
            #c.Print('{}.pdf'.format(inhist.GetName()))

if __name__=='__main__':
    plot()
