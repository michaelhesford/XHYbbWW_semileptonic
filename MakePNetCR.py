import ROOT, time
ROOT.gROOT.SetBatch(True)
from argparse import ArgumentParser
from collections import OrderedDict
from TIMBER.Tools.Common import CompileCpp
from XHYbbWW_class import XHYbbWW

def PNetSnapshot(self, node=None):
    startNode = self.a.GetActiveNode()
    if node == None:
        node = self.a.GetActiveNode()

    columns = [
        'nJet','BJet_pt','BJet_eta','BJet_phi','BJet_mass','BJet_jetId',
        'nFatJet','ProbeJet_eta','ProbeJet_msoftdrop','ProbeJet_pt','ProbeJet_phi','ProbeJet_JES_nom','ProbeJet_particleNetMD*',
        'ProbeJet_jetId','ProbeJet_jetFlavor','nMuon','Muon_*','nElectron','Electron_*','LeptonType','Lepton_*',
        'ChsMET_phi','ChsMET_pt','ChsMET_sumEt',
        'genWeight','event','eventWeight','luminosityBlock','run','NSTART','NFLAGS','nkinLep','nBtag','nMET','nProbe'
    ]
    columns.extend(['nGenJet','nSoftActivityJet','nSubJet'])

    #Triggers
    all_trigs = self.trigs[self.year if 'APV' not in self.year else '16']
    for group in ['HADRONIC','ELECTRON','MUON']:
        columns.extend(all_trigs[group])

    if not self.a.isData:
        columns.extend(['nGenPart','GenPart_*','GenMET_*','genWeight'])
        columns.extend(['ProbeJet_JES_up','ProbeJet_JES_down',
                        'ProbeJet_JER_nom','ProbeJet_JER_up','ProbeJet_JER_down',
                        'ProbeJet_JMS_nom','ProbeJet_JMS_up','ProbeJet_JMS_down',
                        'ProbeJet_JMR_nom','ProbeJet_JMR_up','ProbeJet_JMR_down'])
        columns.extend(['Pileup__nom','Pileup__up','Pileup__down','Pdfweight__up','Pdfweight__down'])
        if self.year == '16' or self.year == '17' or 'APV' in self.year:
            columns.extend(['L1PreFiringWeight_Nom', 'L1PreFiringWeight_Up', 'L1PreFiringWeight_Dn'])       # these are the default columns in NanoAODv9
            columns.extend(['L1PreFiringWeight__nom','L1PreFiringWeight__up','L1PreFiringWeight__down'])    # these are the weight columns created by the BranchCorrection module
        elif self.year == '18':
            columns.append('HEM_drop__nom')

    # get ready to send out snapshot
    self.a.SetActiveNode(node)
    self.a.Snapshot(columns, 'PNetCRSnapshot_{}_{}_{}_{}.root'.format(self.setname,self.year,self.ijob,self.njobs),'Events', openOption='RECREATE',saveRunChain=True)
    self.a.SetActiveNode(startNode)

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-s', type=str, dest='setname',
		    action='store', required=True,
		    help='Setname to process.')
    parser.add_argument('-y', type=str, dest='year',
		    action='store', required=True,
		    help='Year of set (16, 17, 18).')
    parser.add_argument('-j', type=int, dest='ijob',
		    action='store', default=1,
		    help='Job number')
    parser.add_argument('-n', type=int, dest='njobs',
		    action='store', default=1,
		    help='Number of jobs')
    args = parser.parse_args()

    start = time.time()

    setname=args.setname
    year=args.year
    ijob=args.ijob
    njobs=args.njobs

    filename='raw_nano/{}_{}.txt'.format(setname,year)
    selection = XHYbbWW(filename,ijob,njobs)
    
    #Get initial event count
    selection.NSTART = selection.getNweighted()
    selection.AddCutflowColumn(selection.NSTART,'NSTART')
    
    #Apply MET filters
    flags = [
        'Flag_goodVertices',
        'Flag_globalSuperTightHalo2016Filter',
        'Flag_HBHENoiseFilter',
        'Flag_HBHENoiseIsoFilter',
        'Flag_EcalDeadCellTriggerPrimitiveFilter',
        'Flag_BadPFMuonFilter',
        'Flag_BadPFMuonDzFilter',
        'Flag_eeBadScFilter'
     ]
    if selection.year == '17' or selection.year == '18':
        flags.append('Flag_ecalBadCalibFilter')
    MET_filters = selection.a.GetFlagString(flags)       # string valid (existing in RDataFrame node) flags together w logical and
    selection.a.Cut('flags', MET_filters)
    selection.NFLAGS = selection.getNweighted()
    selection.AddCutflowColumn(selection.NFLAGS, "NFLAGS")
    
    #Apply some corrections
    selection.ApplyStandardCorrections(snapshot=True)

    '''
    #Find true MC flavor for fat jet
    if not selection.a.isData:
        selection.a.Define('FatJet_jetFlavor', 'JetFlavor_ttbar(FatJet_eta, FatJet_phi, GenPart_eta, GenPart_phi, GenPart_pdgId, GenPart_genPartIdxMother)')
    '''

    #Make ttbar-enriched control region
    selection.MakePNetCR()
  
    #create cut flow table
    outFile = ROOT.TFile.Open('{}_{}_{}_PNetSnapCutFlow.root'.format(setname,year,ijob),'RECREATE')
    outFile.cd()
        
    cutflowInfo = OrderedDict([
        ('start',selection.NSTART),
        ('flags',selection.NFLAGS),
        ('tag lepton',selection.nkinLep),
        ('tag bjet',selection.nBtag),
        ('MET>30',selection.nMET),
        ('probe ak8',selection.nProbe)	
    ])

    nLabels = len(cutflowInfo)
    hCutflow = ROOT.TH1F('cutflow', 'Number of events after each cut', nLabels, 0.5, nLabels+0.5)
    nBin = 1
    for label, value in cutflowInfo.items():
        hCutflow.GetXaxis().SetBinLabel(nBin, label)
        hCutflow.AddBinContent(nBin, value)
        nBin += 1
    hCutflow.Write()
    
    outFile.Close()

    PNetSnapshot(selection)

    print('%s sec'%(time.time()-start))
