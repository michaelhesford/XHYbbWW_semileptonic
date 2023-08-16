from argparse import Namespace
from glob import glob
from TwoDHistStyle import style
from TIMBER.Tools.Common import DictStructureCopy, CompileCpp, ExecuteCmd, OpenJSON, StitchQCD
from TIMBER.Tools.Plot import CompareShapes, EasyPlots
from TIMBER.Analyzer import HistGroup
from collections import OrderedDict
import ROOT

def plot():
    #make plots of nominal input histograms used for the background fit
    #now do signal
    files = [f for f in glob('selection/XHYbbWWselection*.root')] #loop through all files
    c = ROOT.TCanvas('c','c')
    c.cd()
    for f in files:
        inFile=ROOT.TFile.Open(f,'READ')
        setname=f.split('/')[1].split('_')[1]
        histlist=[] #for 2D plots
        for key in inFile.GetListOfKeys():
            varname=key.GetName()
            if not varname.endswith('nominal'):
                continue
            inhist = inFile.Get(varname)  # get it from the file
            inhist.SetDirectory(0) # set the directory so hist is stored in memory and not as reference to TFile (this way it doesn't get tossed by python garbage collection when infile changes)      
            inhist.SetName('{}_'.format(setname)+varname)
            c.Clear
            style(inhist,xTitle='m_{X} (GeV)',yTitle='m_{Y} (GeV)',histTitle=inhist.GetTitle())
            inhist.Draw('COLZ')
            c.Print('selection/png/{}.png'.format(inhist.GetName()))

if __name__=='__main__':
    plot()               
