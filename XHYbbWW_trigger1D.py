import sys, time, ROOT
from collections import OrderedDict
from array import array
from TIMBER.Analyzer import HistGroup
from TIMBER.Tools.Common import CompileCpp, ExecuteCmd
from XHYbbWW_class import XHYbbWW

def MakeEfficiency(year,ijob,njobs):
    '''
    year (str) : 16, 17, 17B, 17noB, 18
    '''
    if year == '17B':
        fName = 'snapshots/DataB_17_snapshot.txt'
    elif year == '17noB':
        fName = 'snapshots/DataNoB_17_snapshot.txt'
    else:
        fName = 'snapshots/Data_{}_snapshot.txt'.format(year)

    selection = XHYbbWW(fName,ijob,njobs)
    #selection.ApplyJMECorrections('None')

    #Want to replicate full selection so that trigger efficiencies are "meaningful"
    selection.KinematicLepton()
    #selection.Dijets() # I think I can leave out jet-related selection - will increase statistics without impacting efficiency measurement for leptonic triggers
    #selection.ApplyMassCuts()
    selection.a.Define('Lepton_abseta','abs(Lepton_eta)') # bin efficiencies in |eta|

    '''
    if '16' in year:
        refTrigger = ["HLT_IsoMu24","HLT_IsoTkMu24","HLT_Mu50"] #2016
    elif '17' in year:
        refTrigger = ["HLT_IsoMu27","HLT_Mu50"] #2017
    elif '18' in year:
        refTrigger = ["HLT_IsoMu24","HLT_Mu50"] #2018
    '''

    if '16' in year:
        refTrigger = ["HLT_PFHT800","HLT_PFHT900","HLT_PFJet450"] #2016
    elif '17' in year:
        refTrigger = ["HLT_PFHT1050","HLT_AK8PFJet400_TrimMass30"] #2017
    elif '18' in year:
        refTrigger = ["HLT_PFHT1050","HLT_AK8PFJet400_TrimMass30"] #2018

    all_trigs = selection.trigs[selection.year if '16' not in selection.year else '16']
    if selection.year == '18':
        trigs = all_trigs['EGAMMA'] + all_trigs['MUON']
    else:
        trigs = all_trigs['ELECTRON'] + all_trigs['MUON'] + all_trigs['PHOTON']

    #### IMPORTANT: check if reference/target triggers are correlated ####

    #Reference
    start = selection.a.GetActiveNode()
    nstart = selection.a.DataFrame.Count().GetValue()
    selection.a.Cut('n_ref',selection.a.GetTriggerString(refTrigger))
    nref = selection.a.DataFrame.Count().GetValue()
    e_ref = nref/nstart

    #Target
    selection.a.SetActiveNode(start)
    selection.a.Cut('n_trig',selection.a.GetTriggerString(trigs))
    ntest = selection.a.DataFrame.Count().GetValue()
    e_test = ntest/nstart

    #Reference + target
    selection.a.SetActiveNode(start)
    selection.a.Cut('n_ref_comb',selection.a.GetTriggerString(refTrigger))
    selection.a.Cut('n_trig_comb',selection.a.GetTriggerString(trigs))
    n_combined = selection.a.DataFrame.Count().GetValue()
    e_combined = n_combined/nstart
    
    a = e_test*e_ref / e_combined #should be close to 1 for uncorrelated triggers

    selection.a.SetActiveNode(start)

    #Print result in a text file
    out = open('out_{}.txt'.format(year),'w')
    out.write('Year: 20'+year+'\n')
    out.write('Reference triggers: '+', '.join(refTrigger)+'\n')
    out.write('Target triggers: '+', '.join(trigs)+'\n')
    out.write('alpha = e_ref*e_target/e_(ref&target) = {}'.format(a))

    ######################################################################

    noTag = selection.a.Cut('pretrig',selection.a.GetTriggerString(refTrigger))

    # Include tagging
    hists = HistGroup('out') 

    pt_bins = [25,50,75,100,125,150,175,200,225,250,275,300,325,350,400,500,1000]  #[25,50,75,100,150,300,600,1000]
    eta_bins = [0,0.3,0.6,0.9,1.2,1.5,1.8,2.5]  #[0,0.9,1.2,2.1,2.5]

    all_lep = selection.a.GetActiveNode() # measure efficiency separately for electron/muon events
    eleEvents=selection.a.Cut('eleEvents','LeptonType == 0')

    hists.Add('lepPtDenominator_ele',selection.a.DataFrame.Histo1D(('lepPtDenominator_ele','',len(pt_bins)-1, array('d',pt_bins)),'Lepton_pt'))
    hists.Add('lepEtaDenominator_ele',selection.a.DataFrame.Histo1D(('lepEtaDenominator_ele','',len(eta_bins)-1, array('d',eta_bins)),'Lepton_abseta'))
    selection.a.Cut('trigger_ele',selection.a.GetTriggerString(trigs))
    hists.Add('lepPtNumerator_ele',selection.a.DataFrame.Histo1D(('lepPtNumerator_ele','', len(pt_bins)-1, array('d',pt_bins)),'Lepton_pt'))
    hists.Add('lepEtaNumerator_ele',selection.a.DataFrame.Histo1D(('lepEtaNumerator_ele','', len(eta_bins)-1, array('d',eta_bins)),'Lepton_abseta'))

    selection.a.SetActiveNode(all_lep)
    muEvents = selection.a.Cut('muEvents','LeptonType == 1')

    hists.Add('lepPtDenominator_mu',selection.a.DataFrame.Histo1D(('lepPtDenominator_mu','',len(pt_bins)-1, array('d',pt_bins)),'Lepton_pt'))
    hists.Add('lepEtaDenominator_mu',selection.a.DataFrame.Histo1D(('lepEtaDenominator_mu','',len(eta_bins)-1, array('d',eta_bins)),'Lepton_abseta'))
    selection.a.Cut('trigger_mu',selection.a.GetTriggerString(trigs))
    hists.Add('lepPtNumerator_mu',selection.a.DataFrame.Histo1D(('lepPtNumerator_mu','', len(pt_bins)-1, array('d',pt_bins)),'Lepton_pt'))
    hists.Add('lepEtaNumerator_mu',selection.a.DataFrame.Histo1D(('lepEtaNumerator_mu','', len(eta_bins)-1, array('d',eta_bins)),'Lepton_abseta'))

    effs = {
        "lepPtEfficiency_ele": ROOT.TEfficiency(hists['lepPtNumerator_ele'], hists['lepPtDenominator_ele']),
        "lepEtaEfficiency_ele": ROOT.TEfficiency(hists['lepEtaNumerator_ele'], hists['lepEtaDenominator_ele']),
        "lepPtEfficiency_mu": ROOT.TEfficiency(hists['lepPtNumerator_mu'], hists['lepPtDenominator_mu']),
        "lepEtaEfficiency_mu": ROOT.TEfficiency(hists['lepEtaNumerator_mu'], hists['lepEtaDenominator_mu'])
    }

    out = ROOT.TFile.Open('XHYbbWWtrigger1D_{}.root'.format(year),'RECREATE')
    #out = ROOT.TFile.Open('XHYbbWWtrigger1D_{}_{}_{}.root'.format(year,ijob,njobs),'RECREATE')
    out.cd()

    hists['lepPtDenominator_ele'].Write()
    hists['lepPtNumerator_ele'].Write()
    hists['lepEtaDenominator_ele'].Write()
    hists['lepEtaNumerator_ele'].Write()

    hists['lepPtDenominator_mu'].Write()
    hists['lepPtNumerator_mu'].Write()
    hists['lepEtaDenominator_mu'].Write()
    hists['lepEtaNumerator_mu'].Write()

    for name, eff in effs.items():
        eff.SetName(name)
        eff.Write()

    out.Close()

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-y', type=str, dest='year',
                        action='store', required=False,
                        help='Year of set (16, 17, 18).')
    parser.add_argument('-j', type=int, dest='ijob',required=False,
                        action='store', default=1, help='current job')
    parser.add_argument('-n', type=int, dest='njobs',required=False,
                        action='store', default=1, help='number of jobs')
    parser.add_argument('--recycle', dest='recycle',
                        action='store_true', default=False,
                        help='Recycle existing files and just plot.')
    args = parser.parse_args()
    start = time.time()
    if not args.recycle:
        MakeEfficiency(args.year,args.ijob,args.njobs)

    else:
        files = {
            #'17noB': ROOT.TFile.Open('plots/trigger1D/XHYbbWWtrigger1D_17noB.root'),
            #'17B': ROOT.TFile.Open('plots/trigger1D/XHYbbWWtrigger1D_17B.root'),
            '18': ROOT.TFile.Open('plots/trigger1D/XHYbbWWtrigger1D_18.root')
        }

        #Make file for full 2016 data
        ExecuteCmd('hadd -f plots/trigger1D/XHYbbWWtrigger1D_16all.root plots/trigger1D/XHYbbWWtrigger1D_16.root plots/trigger1D/XHYbbWWtrigger1D_16APV.root')
        files['16'] = ROOT.TFile.Open('plots/trigger1D/XHYbbWWtrigger1D_16all.root')

        #Do the same for 2017
        ExecuteCmd('hadd -f plots/trigger1D/XHYbbWWtrigger1D_17.root plots/trigger1D/XHYbbWWtrigger1D_17B.root plots/trigger1D/XHYbbWWtrigger1D_17noB.root')
        files['17'] = ROOT.TFile.Open('plots/trigger1D/XHYbbWWtrigger1D_17.root')

        #Make nicely formatted histograms
        effs = {effname.GetName():[files[y].Get(effname.GetName()) for y in ['16','17','18']] for effname in files['16'].GetListOfKeys()}
        colors = [ROOT.kBlack, ROOT.kRed, ROOT.kOrange-3, ROOT.kGreen+1]
        legendNames = ['2016','2017','2018']
        for effname in effs.keys():
            c = ROOT.TCanvas('c','c',2000,500)
            c.Divide(3,1)
            obj = 'Electron' if 'ele' in effname else 'Muon'
            if 'Pt' in effname:
                varname = f'{obj} p_{{T}} (GeV)'
            elif 'Eta' in effname:
                varname = f'{obj} |#eta|'
            for i,h in enumerate(effs[effname]):
                c.cd(i+1)
                if 'Efficiency' in effname:
                    h = h.CreateGraph()
                    h.SetName(effname+'_graph')
                    h.SetTitle(effname)
                    h.GetXaxis().SetTitle(varname)
                    h.GetYaxis().SetTitle('Efficiency')
                    h.SetMinimum(0.0)
                    h.SetMaximum(1.0)
             
                h.GetYaxis().SetTitleOffset(1.7)
                h.SetLineColor(colors[i])
                h.SetTitle(legendNames[i])
                ROOT.gPad.SetLeftMargin(0.13)
                ROOT.gPad.SetRightMargin(0.16)
 
                if 'Efficiency' in effname:
                    h.Draw('AP')

                elif 'Numerator' in effname or 'Denominator' in effname:
                    h.GetXaxis().SetTitle(varname)
                    h.GetYaxis().SetTitle('Events')
                    h.Draw()

            c.Print('plots/trigger1D/pdf/Trigger1D_{}.pdf'.format(effname),'pdf')

    print('%s sec'%(time.time()-start))
