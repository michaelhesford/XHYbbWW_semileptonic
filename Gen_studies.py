from XHYbbWW_class import XHYbbWW
from XHYbbWW_studies import XHYbbWW_studies
from TIMBER.Analyzer import HistGroup, VarGroup
from TIMBER.Tools.Common import CompileCpp
from TIMBER.Tools.Plot import *
from collections import OrderedDict
from argparse import ArgumentParser

def Plotter(self,variables):
    if not 'histgroup' in globals(): #no histgroup created yet, use for first instance of function
        global histgroup
        histgroup = HistGroup(self.setname)
    for varname in variables:
        if 'msoftdrop' in varname:
            binning = [20,0,200]
        elif '_pt' in varname:
            binning = [40,150,2000]
        elif 'particleNet' in varname:
            binning = [20,0,1]
        elif 'DeltaPhi' in varname:
            binning = [30,-3.2,3.2]
        elif 'DeltaEta' in varname:
            binning = [20,0,4.8]
        elif varname.startswith('n'):
            binning = [10,0,10]
        elif 'RelIso' in varname:
            binning = [20,0,1]
        elif 'W_mass' in varname:
            if 'XHY' in self.setname: #looking at signal
                if int(self.setname.split('_')[0].split('-')[2]) > 1000: #expected y mass > 1000   
                    binning = [25,0,1000]
                else:
                    binning = [20,0,800]
            else: #looking at background
                binning = [20,0,800] 
        elif 'Y_mass' in varname:
            if 'XHY' in self.setname:
                if int(self.setname.split('_')[0].split('-')[2]) > 1000: #expected y mass > 1000   
                    binning = [60,100,3000]
                else:
                    binning = [30,100,1600]
            else:
                binning = [30,100,1600]

        histname = '{}'.format(varname)
        hist_tuple = (histname,histname,binning[0],binning[1],binning[2])
        print('creating {} histogram'.format(histname))
        hist = self.a.GetActiveNode().DataFrame.Histo1D(hist_tuple,varname,'weight__nominal')
        hist.GetValue() # This gets the actual TH1 instead of a pointer to the TH1
        histgroup.Add(varname,hist)

def GenPartMatching(self):

    self.a.SubCollection('GenLep','GenPart','GenPart_pdgId == 11 || GenPart_pdgId == 13 || GenPart_pdgId == -11 || GenPart_pdgId == -13')
    self.a.Define("GenLep_pdgIdMother","FindMothersPdgId(GenPart_pdgId,GenLep_genPartIdxMother)")
    self.a.SubCollection('GenLepfromW','GenLep','GenLep_pdgIdMother == 24 || GenLep_pdgIdMother == -24',skip=['idx'])

    self.a.Define('GenLepIdx','HighestPtIdx(GenLepfromW_pt)') #get index of highest pt generator electron/muon w/ parent W

    self.a.Cut('GenLepExists','GenLepIdx != -1')
    self.GENLEP = self.getNweighted()
    self.AddCutflowColumn(self.GENLEP,'GENLEP')

    self.a.ObjectFromCollection('GenLepLead','GenLepfromW','GenLepIdx') #pull out highest pt generator electron/muon w/ parent W
    
    self.a.Cut('GenLep_pt','GenLepLead_pt > 25')

    self.a.Define('DeltaPhi_GenLep_Higgs','hardware::DeltaPhi(Higgs_phi,GenLepLead_phi)')
    self.a.Define('DeltaPhi_GenLep_Wqq','hardware::DeltaPhi(Wqq_phi,GenLepLead_phi)')

    self.a.Cut('DeltaPhi_GenLep_Jets','abs(DeltaPhi_GenLep_Higgs) > 1.7 && abs(DeltaPhi_GenLep_Wqq) > 1 && abs(DeltaPhi_GenLep_Wqq) < 2.5')
    self.GENLEPPHI = self.getNweighted()
    self.AddCutflowColumn(self.GENLEPPHI,'GENLEPPHI')
    
    return self.a.GetActiveNode()

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
    ijob=args.ijob
    njobs=args.njobs

    #####NOTE:THIS BIT HAS BEEN MOVED TO THE CLASS FILE, NOW JUST INPUT SETNAME/YEAR INTO XHYBBWW_CLASS()#####

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

    ###########################################################################################################

    ana = XHYbbWW(filename,ijob,njobs)

    #reweight events based on lumi*xsec/geneventweight
    ana.a.MakeWeightCols(extraNominal='' if ana.a.isData else 'genWeight*%s'%ana.GetXsecScale())
    
    ana.BasicKinematics()
    ana.Dijets()
    mW=[60,100]
    mH=[100,150]
    Wtag=0.8
    Htag=0.8
    ana.JetSelection(mW,mH,Wtag,Htag)

    GenPartMatching(ana)

    #Lorentz vectors to use for invariant mass calculations
    ana.a.Define('W_vect','hardware::TLvector(Wqq_pt,Wqq_eta,Wqq_phi,Wqq_msoftdrop)')
    ana.a.Define('Lepton_vect','hardware::TLvector(GenLepLead_pt,GenLepLead_eta,GenLepLead_phi,GenLepLead_mass)')
    ana.a.Define('MET_vect','hardware::TLvector(MET_fiducialGenPt,0,MET_fiducialGenPhi,0)') #neutrino mass negligable - for now, assuming pt_miss has no z component 

    # Transverse/invariant masses of W -> l/v (should come out around 80 GeV)
    ana.a.Define('W_massTran','hardware::TransverseMass(MET_fiducialGenPt,GenLepLead_pt,MET_fiducialGenPhi,GenLepLead_phi)') #using generator-level MET variables
    ana.a.Define('W_massInv','hardware::InvariantMass({Lepton_vect,MET_vect})')

    #Finally, invariant mass of Y!
    ana.a.Define('Y_massInv','hardware::InvariantMass({W_vect,Lepton_vect,MET_vect})')

    Plotter(ana,['W_massTran','W_massInv','Y_massInv'])

    # save the raw histos to a file
    outFile = ROOT.TFile.Open('{}_{}_{}_GenMatching.root'.format(setname,year,ijob),'RECREATE')
    outFile.cd()
    histgroup.Do('Write') # This will call TH1.Write() for all of the histograms in the group
    
    #create cut flow table
    cutflowInfo = OrderedDict([
        ('start',ana.NSTART),
        ('flags',ana.NFLAGS),
        ('1 lepton',ana.NLEPTON),
        ('lepton pt>25',ana.LEPPT),
        ('lepton iso<1',ana.LEPISO),
        ('2 fat jets',ana.NJET),
        ('2 jets mass>50',ana.NFATJETMASS),
        ('lead/sublead pt',ana.JETPT),
        ('b-to-b dijets',ana.NDIJETS),
        ('H/Wqq masses',ana.JETMASS),
        ('H/Wqq tagged',ana.JETTAG),
        ('genLep exists',ana.GENLEP),
        ('genLep phi',ana.GENLEPPHI)
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




