import ROOT
from TIMBER.Analyzer import analyzer, HistGroup
from TIMBER.Tools.Plot import *
#from CompareShapes import CompareShapes
from collections import OrderedDict
from argparse import ArgumentParser
ROOT.gROOT.SetBatch(True)

####################################
# Store some global variables here #
####################################

def PrettyPlots(dataset,varnames,mass_window,signals={},backgrounds={},colors={},new_bins={}):
    #mass_window = 'high' or 'low,' use to split plots for a given variable into 2 separate plots w/ different signal masses
    #new_bins used to manually change the binning on a histogram. I haven't been able to get this to work so far :(
    SB=OrderedDict()
    SB.update(signals)
    SB.update(backgrounds)
    histgroups = {}
    for s in SB.keys():
        inFile = ROOT.TFile.Open('plots/{}/{}.root'.format(dataset,s),'READ')
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
    '''
    for b in backgrounds.keys():
        inFile = ROOT.TFile.Open('plots/{}/{}.root'.format(dataset,b),'READ')
        print(inFile)
        # Put histograms into the HistGroups
        histgroups[b] = HistGroup(b)
        for key in inFile.GetListOfKeys():      # loop over histograms in the file
            varname = key.GetName()             # gets the histogram name 
            inhist = inFile.Get(key.GetName())  # get it from the file
            inhist.SetDirectory(0) # set the directory so hist is stored in memory and not as reference to TFile (this way it doesn't get tossed by python garbage collection when infile changes)
            if varname in new_bins.keys():
                newhist=inhist.Clone(varname)
                newhist.SetDirectory(0)
                newhist.GetXaxis().SetRangeUser(new_bins[varname][0],new_bins[varname][1])
                histgroups[b].Add(varname,newhist) # add to our group
            else:
                histgroups[b].Add(varname,inhist)
            print('Added {} distribution for background {}'.format(varname, b))
        inFile.Close()
    '''
    # Now plot the variables up in the global definitions above
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

if __name__=="__main__": #eventually clean this up to do something
    # store the histograms we want to track
    parser = ArgumentParser()
    #parser.add_argument('-y',type=str,dest='year',
    #                    action='store',help='year to plot')
    args=parser.parse_args()
    #year = args.year

    #variable names to plot
    varnames = {
        'FatJet_msoftdrop':'FatJet m_{SD} (GeV)',
        'IsoLepton_pt':'Lepton p_{T} (GeV)',
        'HighPtElectron_miniPFRelIso_all':'e^{-} Rel_Iso',
        'HighPtMuon_miniPFRelIso_all':'#mu^{-} mini rel_iso',
        'Wqq_msoftdrop':'W_{qq} m_{SD} (GeV)',
        'Higgs_msoftdrop':'H_{bb} m_{SD} (GeV)',
        'DeltaPhi_Higgs_Wqq':'#Delta #phi H_{bb},W_{qq}',
        'DeltaEta_Higgs_Wqq':'#Delta #eta H_{bb},W_{qq}',
        'DeltaPhi_IsoLepton_MET':'#Delta #phi lepton,p_{T}^{miss}',
        'DeltaPhi_IsoLepton_Higgs':'#Delta #phi Lepton,H_{bb}',
        'DeltaPhi_IsoLepton_Wqq':'#Delta #phi Lepton,W_{qq}',
        'DeltaPhi_FatJets':'#Delta #phi j_{L},j_{S}',
        'DeltaEta_FatJets':'#Delta #eta j_{L},j_{S}',
        'Wqq_particleNetMD_WqqvsQCD':'W_{qq} particleNetMD_WqqvsQCD',
        'Wqq_particleNetMD_HbbvsQCD':'W_{qq} particleNetMD_HbbvsQCD',
        'Higgs_particleNetMD_HbbvsQCD':'H_{bb} particleNetMD_HbbvsQCD',
        'Higgs_particleNetMD_WqqvsQCD':'H_{bb} particleNetMD_WqqvsQCD',
        'Y_mass_allCuts':'Y mass_{inv} (GeV)',
        'W_massInv_allCuts':'W mass_{inv} (GeV)',
        'W_massTran_allCuts':'W mass_{tran} (GeV)',
        'W_massTran_genMET_allCuts':'W mass_{tran} (generator MET) (GeV)',
        'Wqq_particleNetMD_cut_nminus1':'W_{qq} particleNetMD_WqqvsQCD (N-1)',
        'Higgs_particleNetMD_cut_nminus1':'H_{bb} particleNetMD_HbbvsQCD (N-1)',
        'DeltaPhi_Lepton_Higgs_cut_nminus1':'#Delta #phi Lepton,H_{bb} (N-1)',
        'DeltaPhi_Lepton_Wqq_cut_nminus1':'#Delta #phi Lepton,W_{qq} (N-1)',
        'Isolation_cut_nminus1':'Lepton mini rel_iso (N-1)',
        'mWqq_particleNetMD_cut_nminus1':'W_{qq} m_{SD} (N-1) (GeV)',
        'mHiggs_particleNetMD_cut_nminus1':'H_{bb} m_{SD} (N-1) (GeV)'
    }



    signals = {

#        'XHY-4000-2500' : 'm_{X}=4000, m_{Y}=2500',
#        'XHY-3500-1800' : 'm_{X}=3500, m_{Y}=1800',
#        'XHY-3000-1400' : 'm_{X}=3000, m_{Y}=1400' #add comma!!
        'XHY-2000-1000' : 'm_{X}=2000, m_{Y}=1000',
        'XHY-1600-700' : 'm_{X}=1600, m_{Y}=700',
        'XHY-1000-500' : 'm_{X}=1000, m_{Y}=500'
    }
#        'XHY-2500-1300' : 'm_d{X}=2500, m_{Y}=1300 (GeV)',
    backgrounds = {
        'QCD' : 'QCD',
        'WJets': 'WJets',
        'ZJets' : 'ZJets',
#        'ttbar-semilep' : 'ttbar_{semilep}',
#        'ttbar-allhad' : 'ttbar_{allhad}'
        'ttbar' : 'ttbar'
    }
    
    SB = { # to use for names (not the most efficient thing but whatever)
#        'XHY-4000-2500' : 'm_{X}=4000, m_{Y}=2500',
#        'XHY-3500-1800' : 'm_{X}=3500, m_{Y}=1800',
#        'XHY-3000-1400' : 'm_{X}=3000, m_{Y}=1400',
        'XHY-2000-1000' : 'm_{X}=2000, m_{Y}=1000',
        'XHY-1600-700' : 'm_{X}=1600, m_{Y}=700',
        'XHY-1000-500' : 'm_{X}=1000, m_{Y}=500',        
        'QCD' : 'QCD',
        'WJets' : 'WJets',
        'ZJets' : 'ZJets',
#        'ttbar-semilep' : 'ttbar_{semilep}',
#        'ttbar-allhad' : 'ttbar_{allhad}'
        'ttbar' : 'ttbar'
    }
    
    # store the colors to plot them
    colors = {
#        'XHY-4000-2500' : ROOT.kRed,
#        'XHY-3500-1800' : ROOT.kBlue,
#        'XHY-3000-1400' : ROOT.kBlack,
#        'XHY-2500-1300' : ROOT.kGreen,
        'XHY-2000-1000' : ROOT.kBlue,
        'XHY-1600-700' : ROOT.kBlack,
        'XHY-1000-500' : ROOT.kRed,
        'QCD' : ROOT.kGreen,
        'WJets' : ROOT.kCyan,
        'ZJets' : ROOT.kOrange,
#        'ttbar-semilep' : ROOT.kYellow,
#        'ttbar-allhad' : ROOT.kGray
        'ttbar' : ROOT.kYellow
    }
    
