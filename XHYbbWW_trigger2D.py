import sys, time, ROOT
from collections import OrderedDict
from array import array
from TIMBER.Analyzer import HistGroup
from TIMBER.Tools.Common import CompileCpp, ExecuteCmd
from XHYbbWW_class import XHYbbWW

def MakeEfficiency(year,ijob,njobs):
    '''
    year (str) : 16, 17B, 17noB, 18
    '''
    if year == '17B':
        fName = 'snapshots/DataB_17_snapshot.txt'
    elif year == '17noB':
        fName = 'snapshots/DataNoB_17_snapshot.txt'
    else:
        fName = 'snapshots/Data_{}_snapshot.txt'.format(year)

    selection = XHYbbWW(fName,ijob,njobs)

    #Want to replicate full selection so that trigger efficiencies are "meaningful"
    selection.KinematicLepton()

    if '16' in year:
        refTrigger = ["HLT_PFHT800","HLT_PFHT900","HLT_PFJet450"] #2016
    elif '17' in year:
        refTrigger = ["HLT_PFHT1050","HLT_AK8PFJet500","HLT_AK8PFHT750_TrimMass50","HLT_AK8PFHT800_TrimMass50","HLT_AK8PFJet400_TrimMass30"] #2017
    elif '18' in year:
        refTrigger = ["HLT_PFHT1050","HLT_AK8PFJet400_TrimMass30","HLT_AK8PFHT850_TrimMass50"] #2018

    all_trigs = selection.trigs[selection.year if '16' not in selection.year else '16']
    if selection.year == '18':
        trigs = all_trigs['EGAMMA'] + all_trigs['MUON']
    else:
        trigs = all_trigs['ELECTRON'] + all_trigs['MUON'] + all_trigs['PHOTON']

    noTag = selection.a.Cut('pretrig',selection.a.GetTriggerString(refTrigger))

    # Include tagging
    hists = HistGroup('out') 

    #pt_bins = [25,50,75,100,125,150,175,200,225,250,275,300,325,350,400,500,1000]  #[25,50,75,100,150,300,600,1000]
    pt_bins = [35,56,77,98,119,140,161,182,203,224,245,266,287,308,329,350,400,500,1000]
    eta_bins = [0,0.3,0.6,0.9,1.2,1.5,1.8,2.5]

    # Calculate overall efficiency
    hists.Add('Denominator',selection.a.DataFrame.Histo2D(('Denominator','', len(eta_bins)-1, array('d',eta_bins), len(pt_bins)-1, array('d',pt_bins)),'Lepton_abseta','Lepton_pt'))
    selection.a.Cut('trigger',selection.a.GetTriggerString(trigs))
    hists.Add('Numerator',selection.a.DataFrame.Histo2D(('Numerator','', len(eta_bins)-1, array('d',eta_bins), len(pt_bins)-1, array('d',pt_bins)),'Lepton_abseta','Lepton_pt'))

    # Electron events
    selection.a.SetActiveNode(noTag)
    eleEvents = selection.a.Cut('eleEvents','LeptonType == 0')

    hists.Add('Denominator_ele',selection.a.DataFrame.Histo2D(('Denominator_ele','', len(eta_bins)-1, array('d',eta_bins), len(pt_bins)-1, array('d',pt_bins)),'Lepton_abseta','Lepton_pt'))
    selection.a.Cut('trigger_ele',selection.a.GetTriggerString(trigs))
    hists.Add('Numerator_ele',selection.a.DataFrame.Histo2D(('Numerator_ele','', len(eta_bins)-1, array('d',eta_bins), len(pt_bins)-1, array('d',pt_bins)),'Lepton_abseta','Lepton_pt'))

    # Muon events
    selection.a.SetActiveNode(noTag)
    muEvents = selection.a.Cut('muEvents','LeptonType == 1')

    hists.Add('Denominator_mu',selection.a.DataFrame.Histo2D(('Denominator_mu','', len(eta_bins)-1, array('d',eta_bins), len(pt_bins)-1, array('d',pt_bins)),'Lepton_abseta','Lepton_pt'))
    selection.a.Cut('trigger_mu',selection.a.GetTriggerString(trigs))
    hists.Add('Numerator_mu',selection.a.DataFrame.Histo2D(('Numerator_mu','', len(eta_bins)-1, array('d',eta_bins), len(pt_bins)-1, array('d',pt_bins)),'Lepton_abseta','Lepton_pt'))

    effs = {
        "Efficiency_ele": ROOT.TEfficiency(hists['Numerator_ele'], hists['Denominator_ele']),
        "Efficiency_mu": ROOT.TEfficiency(hists['Numerator_mu'], hists['Denominator_mu']),
        "Efficiency": ROOT.TEfficiency(hists['Numerator'], hists['Denominator'])
    }

    out = ROOT.TFile.Open('XHYbbWWtrigger2D_{}.root'.format(year),'RECREATE')
    #out = ROOT.TFile.Open('XHYbbWWtrigger2D_{}_{}_{}.root'.format(year,ijob,njobs),'RECREATE')
    out.cd()

    hists['Denominator_ele'].Write()
    hists['Numerator_ele'].Write()
    hists['Denominator_mu'].Write()
    hists['Numerator_mu'].Write()
    hists['Denominator'].Write()
    hists['Numerator'].Write()

    for name, eff in effs.items():
        eff.SetName(name)
        eff.Write()

    out.Close()

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-y', type=str, dest='year',
                        action='store', required=False,
                        help='Year of set (16, 17, 17B, 18).')
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
            '18': ROOT.TFile.Open('plots/trigger2D/XHYbbWWtrigger2D_18.root')
        }

        #Make file for full 2016 data
        ExecuteCmd('hadd -f plots/trigger2D/XHYbbWWtrigger2D_16all.root plots/trigger2D/XHYbbWWtrigger2D_16.root plots/trigger2D/XHYbbWWtrigger2D_16APV.root')
        files['16'] = ROOT.TFile.Open('plots/trigger2D/XHYbbWWtrigger2D_16all.root')

        #Do the same for 2017
        ExecuteCmd('hadd -f plots/trigger2D/XHYbbWWtrigger2D_17.root plots/trigger2D/XHYbbWWtrigger2D_17B.root plots/trigger2D/XHYbbWWtrigger2D_17noB.root')
        files['17'] = ROOT.TFile.Open('plots/trigger2D/XHYbbWWtrigger2D_17.root')

        #Make nicely formatted histograms
        effs = {effname.GetName():[files[y].Get(effname.GetName()) for y in ['16','17','18']] for effname in files['16'].GetListOfKeys()}
        colors = [ROOT.kBlack, ROOT.kRed, ROOT.kOrange-3, ROOT.kGreen+1]
        legendNames = ['2016','2017','2018']
        c = ROOT.TCanvas('c','c',2000,2000)
        for effname in effs.keys():
            for i,h in enumerate(effs[effname]):
                c.Clear()
                
                if 'ele' in effname:
                    obj = 'Electron'
                elif 'mu' in effname:
                    obj = 'Muon'
                else:
                    obj = 'Lepton'
                
                if 'Efficiency' in effname:
                    h = h.CreateHistogram()
                    h.SetName(effname+'_hist')
                    h.SetTitle(effname)
                    h.GetZaxis().SetTitle('Efficiency')
                    h.SetMinimum(0.0)
                    h.SetMaximum(1.0)
                else:
                    h.GetZaxis().SetTitle('Events') 

                ROOT.gPad.SetLeftMargin(0.13)
                ROOT.gPad.SetRightMargin(0.16) 

                h.GetYaxis().SetTitle(f'{obj} p_{{T}} (GeV)')
                h.GetXaxis().SetTitle(f'{obj} |#eta|')           
                h.GetZaxis().SetTitleOffset(1.7)
                h.SetLineColor(colors[i])
                h.SetTitle(legendNames[i])
                h.Draw('COLZ TEXT')

                '''
                # Manually draw outlines around each bin
                boxes = []  # Store references to TBox objects to avoid garbage collection
                for a in range(h.GetNbinsX()):
                    for b in range(h.GetNbinsY()):
                        x1 = h.GetXaxis().GetBinLowEdge(a)
                        x2 = h.GetXaxis().GetBinUpEdge(a)
                        y1 = h.GetYaxis().GetBinLowEdge(b)
                        y2 = h.GetYaxis().GetBinUpEdge(b)
                        box = ROOT.TBox(x1, y1, x2, y2)
                        box.SetLineColor(ROOT.kBlack)
                        box.SetFillStyle(0)  # Transparent fill
                        box.Draw("L")
                        boxes.append(box)
        
                # Update the canvas to show the drawn histogram
                c.Update()
                '''
                c.Print('plots/trigger2D/pdf/Trigger2D_{}_{}.pdf'.format(effname,legendNames[i]),'pdf')

    print('%s sec'%(time.time()-start))
