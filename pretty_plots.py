import ROOT
from TIMBER.Analyzer import analyzer, HistGroup
from TIMBER.Tools.Plot import *
#from CompareShapes import CompareShapes
from collections import OrderedDict
from argparse import ArgumentParser
ROOT.gROOT.SetBatch(True)

def PrettyPlots(dataset,varnames,mass_window,signals={},backgrounds={},colors={},new_bins={}):
    #mass_window = 'high' or 'low,' use to split plots for a given variable into 2 separate plots w/ different signal masses
    #new_bins used to manually change the binning on a histogram. I haven't been able to get this to work so far :(
    SB=OrderedDict()
    SB.update(signals)
    SB.update(backgrounds)
    histgroups = {}
    for s in SB.keys():
        inFile = ROOT.TFile.Open('plots/{}/{}_{}.root'.format(dataset,s,dataset),'READ')
        print(inFile)
        # Put histograms into the HistGroups
        histgroups[s] = HistGroup(s)
        for key in inFile.GetListOfKeys():      # loop over histograms in the file
            varname = key.GetName()             # gets the histogram name 
            inhist = inFile.Get(key.GetName())  # get it from the file
            inhist.SetDirectory(0) # set the directory so hist is stored in memory and not as reference to TFile (this way it doesn't get tossed by python garbage collection when infile changes)
            if varname in new_bins.keys():
                newhist=inhist.Clone(varname)
                newhist.SetDirectory(0)
                newhist.GetXaxis().SetRangeUser(new_bins[varname][0],new_bins[varname][1])
                histgroups[s].Add(varname,newhist) # add to our group
                print('newhist added for {},{}'.format(varname,s))
            else:
                histgroups[s].Add(varname,inhist)             
            print('Added {} distribution for {}'.format(varname, s))
        inFile.Close()
    
    # Now plot the variables
    for varname in varnames.keys():
        sig_hists = OrderedDict() #dict of signal hists 
        bkg_hists = OrderedDict()
        for sig in signals.keys():
            sig_hists[sig] = histgroups[sig][varname]
            sig_hists[sig].SetTitle(varname)
        for bkg in backgrounds.keys():
            bkg_hists[bkg] = histgroups[bkg][varname]
            bkg_hists[bkg].SetTitle(varname)

        plot_filename = 'plots/{}/png/{}.png'.format(dataset,varname+'_'+mass_window)#dataset
        CompareShapes(
            outfilename = plot_filename,
            year = 2,
            prettyvarname = varnames[varname],
            bkgs = bkg_hists,
            signals = sig_hists,
            colors = colors,
            names = SB,
            scale = True,
            stackBkg = True
        )

