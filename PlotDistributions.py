import ROOT
from XHYbbWW_class import XHYbbWW
from TIMBER.Analyzer import HistGroup
from TIMBER.Tools.Common import CompileCpp
from TIMBER.Tools.Plot import *
from collections import OrderedDict
from argparse import ArgumentParser

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('-s', type=str, dest='setname',
                        action='store',help='name of data set to run on')
    parser.add_argument('-y', type=str, dest='year',
			action='store', help='year',required=False)
    parser.add_argument('-j', type=int, dest='ijob',required=False,
                        action='store', help='current job')
    parser.add_argument('-n', type=int, dest='njobs',required=False,
                        action='store', help='number of jobs')
    args = parser.parse_args()

    setname=args.setname
    year=args.year

    variables = [
        'CloseJet_btag',
        'Y_massInv',
        'CloseLeptonMET_massInv',
        'W_massTran',
        'W_massTran_METfid',
        'FatJet_pt',
        'FatJet_msoftdrop',
        'FatJet_eta',
        'Electron_pt',
        'Electron_eta',
        'Electron_miniPFRelIso_all',
        'Muon_pt',
        'Muon_eta',
        'Muon_miniPFRelIso_all',
        'MET_pt',
        'MET_phi',
        'MET_fiducialGenPt',
        'METfid_mass',
        'Wqq_msoftdrop',
        'Higgs_msoftdrop',
        'DeltaPhi_Jets',
        'DeltaPhi_Electron_Lead',
        'DeltaPhi_Muon_Lead',
        'DeltaPhi_Electron_Sublead',
        'DeltaPhi_Muon_Sublead',
        'DeltaPhi_Lead_MET',
        'DeltaPhi_Sublead_MET',
        'DeltaPhi_Electron_MET',
        'DeltaPhi_Muon_MET',
        'DeltaPhi_TagJets',
        'DeltaPhi_Lepton_Higgs',
        'DeltaPhi_Lepton_Wqq',
        'DeltaPhi_Higgs_MET',
        'DeltaPhi_Wqq_MET',
        'DeltaPhi_Lepton_MET'
    ]


    if 'XHY' in setname:
        signal=True
        filename = 'raw_nano/Signal/{}_{}.txt'.format(setname,year)
    else:
        signal=False
        if 'QCD' in setname:
            filename = 'raw_nano/Background/QCD/{}_{}.txt'.format(setname,year)
        elif 'ttbar' in setname:
            filename = 'raw_nano/Background/ttbar/{}_{}.txt'.format(setname,year)
        elif 'Jets' in setname:
            filename = 'raw_nano/Background/V+Jets/{}_{}.txt'.format(setname,year)

    ana = XHYbbWW(filename,args.ijob,args.njobs)

    CompileCpp('delta_phi.cc')   

    #function to use for later
    def GetDeltaPhi(name,objects,idxs): #objects: [obj1, obj2] indexes: [idx1, idx2]
        title = 'DeltaPhi_'+name
        if idxs[0] == -1 and idxs[1] == -1:
            ana.a.Define(title,'DeltaPhiSingular2(%s,%s)' %(objects[0],objects[1]))
        elif -1 in idxs:
            idx=idxs[0]
            ana.a.Define(title,'DeltaPhiSingular(%s,%s,%i)' %(objects[0],objects[1],idx))
        else:
            ana.a.Define(title,'DeltaPhi(%s,%s,{%i,%i})' %(objects[0],objects[1],idxs[0],idxs[1]))

        return ana.a.GetActiveNode() 

    #first ensure enough leptons/jets exist to avoid referencing indexes which don't exist
    ana.a.Cut('nLepton','nElectron > 1 && nMuon > 1')
    ana.NLEPTON = ana.getNweighted()
    ana.AddCutflowColumn(ana.NLEPTON,'NLEPTON')

    ana.a.Cut('nFatJet','nFatJet > 1')
    ana.NJET = ana.getNweighted()
    ana.AddCutflowColumn(ana.NJET,'NJET')

    #calculate delta phi between different objects
    DeltaPhi_objs = {}
    DeltaPhi_idxs = {}

    DeltaPhi_objs['Electron_MET'] = ['Electron_phi','MET_phi']
    DeltaPhi_objs['Muon_MET']= ['Muon_phi','MET_phi']
    DeltaPhi_objs['Jets']= ['FatJet_phi','FatJet_phi']
    DeltaPhi_objs['Electron_Lead']= ['Electron_phi','FatJet_phi']
    DeltaPhi_objs['Muon_Lead']= ['Muon_phi','FatJet_phi']
    DeltaPhi_objs['Electron_Sublead']= ['Electron_phi','FatJet_phi']
    DeltaPhi_objs['Muon_Sublead']= ['Muon_phi','FatJet_phi']
    DeltaPhi_objs['Lead_MET'] = ['FatJet_phi','MET_phi']
    DeltaPhi_objs['Sublead_MET'] = ['FatJet_phi','MET_phi']   
 
    DeltaPhi_idxs['Electron_MET'] = [0,-1]
    DeltaPhi_idxs['Muon_MET'] = [0,-1]
    DeltaPhi_idxs['Jets']= [0,1]
    DeltaPhi_idxs['Electron_Lead']= [0,0]
    DeltaPhi_idxs['Muon_Lead']= [0,0]
    DeltaPhi_idxs['Electron_Sublead']= [0,1]
    DeltaPhi_idxs['Muon_Sublead']= [0,1]
    DeltaPhi_idxs['Lead_MET'] = [0,-1]
    DeltaPhi_idxs['Sublead_MET'] = [1,-1]

    for name in DeltaPhi_objs.keys():
        objects = DeltaPhi_objs[name]
        idxs = DeltaPhi_idxs[name]
        GetDeltaPhi(name,objects,idxs)

    PreSelection_node = ana.a.GetActiveNode() # use later when making delta phi plots

    #Now after getting the delta phi's lets do a more specific jet selection
    ana.SignalLepton()
    ana.Dijets()
    ana.JetSelection(mW=[60,100], mH=[100,150], Wtag=0.8, Htag=0.8)

    #calculate b-tag score of jet closeset to lepton
    ana.a.Define('CloseIdx','GetClosestJet(LeptonIdxs,Electron_phi,Muon_phi,Jet_phi)')
    ana.a.Define('CloseJet_btag','Jet_btagDeepFlavB[CloseIdx]') 

    #Calculate invariant mass of (lepton+MET+W) and (lepton+MET+closest jet)
    #First merge some lepton columns
    ana.a.Define('Lepton_pt','LeptonIdxs[0] == -1 ? Muon_pt[LeptonIdxs[1]] : Electron_pt[LeptonIdxs[0]]')
    ana.a.Define('Lepton_eta','LeptonIdxs[0] == -1 ? Muon_eta[LeptonIdxs[1]] : Electron_eta[LeptonIdxs[0]]')
    ana.a.Define('Lepton_phi','LeptonIdxs[0] == -1 ? Muon_phi[LeptonIdxs[1]] : Electron_phi[LeptonIdxs[0]]')
    ana.a.Define('Lepton_mass','LeptonIdxs[0] == -1 ? Muon_mass[LeptonIdxs[1]] : Electron_mass[LeptonIdxs[0]]')

    #for H/Wqq mass plots
    ana.a.Define('Higgs_msoftdrop','Dijet_msoftdrop[1]')
    ana.a.Define('Wqq_msoftdrop','Dijet_msoftdrop[0]')

    #calculate more delta phi's (because this is my life now)
    DeltaPhi_objs = {}
    DeltaPhi_idxs = {}

    DeltaPhi_objs['TagJets'] = ['Dijet_pt','Dijet_pt']
    DeltaPhi_objs['Lepton_Higgs'] = ['Dijet_phi','Lepton_phi']
    DeltaPhi_objs['Lepton_Wqq'] = ['Dijet_phi','Lepton_phi']
    DeltaPhi_objs['Higgs_MET'] = ['Dijet_phi','MET_phi']
    DeltaPhi_objs['Wqq_MET'] = ['Dijet_phi','MET_phi']
    DeltaPhi_objs['Lepton_MET'] = ['Lepton_phi','MET_phi']

    DeltaPhi_idxs['TagJets'] = [0,1]
    DeltaPhi_idxs['Lepton_Higgs'] = [1,-1]
    DeltaPhi_idxs['Lepton_Wqq'] = [0,-1]
    DeltaPhi_idxs['Higgs_MET'] = [1,-1]
    DeltaPhi_idxs['Wqq_MET'] = [0,-1]
    DeltaPhi_idxs['Lepton_MET'] = [-1,-1]

    for name in DeltaPhi_objs.keys():
        objects = DeltaPhi_objs[name]
        idxs = DeltaPhi_idxs[name]
        GetDeltaPhi(name,objects,idxs)

    #Use m^2 = E^2 - |p|^2 to calculate mass contribution of MET
    ana.a.Define('MET_mass','GetMETmass(MET_pt,MET_sumEt)')
    ana.a.Define('METfid_mass','GetMETmass(MET_fiducialGenPt,MET_sumEt)')
    ana.a.Define('MET_eta','GetMETeta()') #return eta = 0

    ana.a.Define('W_vect','hardware::TLvector(Dijet_pt[0],Dijet_eta[0],Dijet_phi[0],Dijet_msoftdrop[0])')
    ana.a.Define('CloseJet_vect','hardware::TLvector(Jet_pt[CloseIdx],Jet_eta[CloseIdx],Jet_phi[CloseIdx],Jet_mass[CloseIdx])')
    ana.a.Define('Lepton_vect','hardware::TLvector(Lepton_pt,Lepton_eta,Lepton_phi,Lepton_mass)') 
    ana.a.Define('MET_vect','hardware::TLvector(MET_pt,MET_eta,MET_phi,MET_mass)') #figure out how to incorporate MET into this calculation

    ana.a.Define('W_massTran','TransverseMass(Lepton_pt,MET_pt)')
    ana.a.Define('W_massTran_METfid','TransverseMass(Lepton_pt,MET_fiducialGenPt)') #see if this differerent MET pt variable makes a big difference
 
    ana.a.Define('Y_massInv','hardware::InvariantMass({W_vect,Lepton_vect,MET_vect})')
    ana.a.Define('CloseLeptonMET_massInv','hardware::InvariantMass({CloseJet_vect,Lepton_vect,MET_vect})')

    PostSelection_node = ana.a.GetActiveNode() #use for plots involving jet/lepton selection

    # create a histgroup to store all the histograms of variables easily
    histgroup = HistGroup(setname)
    for varname in variables: # Define binning for different plotting variables
        ana.a.SetActiveNode(PreSelection_node) 
	if 'msoftdrop' in varname:
            if not 'Jet' in varname:
                ana.a.SetActiveNode(PostSelection_node)
	    binning = [20,0,200]
	elif 'DeltaPhi' in varname:
            if 'Tag' or 'Lepton' or 'Higgs' or 'Wqq' in varname:
                ana.a.SetActiveNode(PostSelection_node)
	    binning = [30,-3.2,3.2]
        elif 'on_pt' in varname:
            binning = [30,0,400]
        elif 'et_pt' in varname:
            binning = [40,0,2000]
        elif 'eta' in varname:
            binning = [40,-3,3]
        elif 'dxy' in varname or 'dz' in varname: 
            binning = [30,-0.2,0.2]
        elif 'hoe' in varname:
            binning = [50,0,0.5]
        elif 'sieie' in varname:
            binning = [30,0,0.05]
        elif 'PInv' in varname:
            binning = [30,-0.15,0.15]
        elif 'Comp' in varname or 'mini' in varname:
            binning = [10,0,1]
        elif 'massInv' in varname:
            ana.a.SetActiveNode(PostSelection_node)
            binning = [50,0,1000]
        elif 'massTran' in varname:
            ana.a.SetActiveNode(PostSelection_node)
            binning = [40,0,400]
        elif 'id_mass' in varname:
            ana.a.SetActiveNode(PostSelection_node)
            binning = [20,0,2*10**(-9)]
        elif 'btag' in varname:
            ana.a.SetActiveNode(PostSelection_node)
            binning = [10,0,1] 

        #plot the histograms
        histname = '{}_{}'.format(setname,varname)
        hist_tuple = (histname,histname,binning[0],binning[1],binning[2])
        print('creating {} histogram'.format(histname))
        hist = ana.a.GetActiveNode().DataFrame.Histo1D(hist_tuple,varname)
        hist.GetValue() # This gets the actual TH1 instead of a pointer to the TH1
        histgroup.Add(varname,hist)

    # save the raw histos to a file
    if signal: #AHHH CHANGE THIS AFTER THE TEST WORKS CHANGE THE DEST AFTER THE TEST AHHHH
        action='RECREATE'
    elif not signal:
        action='RECREATE'
    #file_path='plots/distributions/root/{}_TEST.root'.format(setname) #GET RID OF THE TEST AHHH
    outFile = ROOT.TFile.Open('{}_{}_{}.root'.format(setname,year,args.ijob),'RECREATE')
    outFile.cd()
    histgroup.Do('Write') # This will call TH1.Write() for all of the histograms in the group

    #create cut flow table
    ana.a.SetActiveNode(PostSelection_node)
    cutflowInfo = OrderedDict([
        ('NLEPTON',ana.NLEPTON),
        ('NJET',ana.NJET),
        ('SIGLEP',ana.SIGLEP),
        ('NDIJETS',ana.NDIJETS),
        ('JETSELECTION',ana.JETSELECTION)
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

'''
#Variables

        'FatJet_msoftdrop',
        'DeltaPhi_Jets',
        'DeltaPhi_Lead_MET',
        'DeltaPhi_Sublead_MET',
        'DeltaPhi_Electron_Lead',
        'DeltaPhi_Electron_Sublead',
        'DeltaPhi_Muon_Lead',
        'DeltaPhi_Muon_Sublead',
        'FatJet_eta',
        'Electron_eta',
        'Electron_pt',
        'Electron_hoe',
        'Electron_sieie',
        'Electron_eInvMinusPInv',
        'Electron_dxy',
        'Electron_dz',
        'Electron_miniPFRelIso_all',
        'Muon_eta',
        'Muon_pt',
        'Muon_segmentComp',
        'Muon_dxy',
        'Muon_dz',
        'Muon_miniPFRelIso_all'

    def HistPlotter(varname,binning): #binning [nbins,val_min,val_max]
        histname = '{}_{}'.format(setname,varname)
        hist_tuple = (histname,histname,binning[0],binning[1],binning[2])
        print('creating {} histogram'.format(histname))
        hist = ana.a.GetActiveNode().DataFrame.Histo1D(hist_tuple,varname)
        return hist.GetValue() # This gets the actual TH1 instead of a pointer to the TH1

        'CloseLepton_mass',
        'Y_mass',

'''
