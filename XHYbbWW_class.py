import ROOT
from TIMBER.Analyzer import Correction, CutGroup, VarGroup, ModuleWorker, analyzer
from TIMBER.Tools.Common import CompileCpp, OpenJSON
from TIMBER.Tools.AutoPU import ApplyPU
from JMEvalsOnly import JMEvalsOnly
import TIMBER.Tools.AutoJME as AutoJME
from collections import OrderedDict

# Helper file for dealing with .txt files containing NanoAOD file locs
def SplitUp(filename,npieces,nFiles=False):
    '''Take in a txt file name where the contents are root
    file names separated by new lines. Split up the files
    into N lists where N is `npieces` in the case that `nFiles == False`.
    In the case that `nFiles == True`, `npieces` is treated as the
    number of files to have per list.
    '''
    files = open(filename,'r').readlines()
    nfiles = len(files)

    if npieces > nfiles:
        npieces = nfiles
    
    if not nFiles: files_per_piece = float(nfiles)/float(npieces)
    else: files_per_piece = npieces
    
    out = []
    iend = 0
    for ipiece in range(1,npieces+1):
        piece = []
        for ifile in range(iend,min(nfiles,int(ipiece*files_per_piece))):
            piece.append(files[ifile].strip())

        iend = int(ipiece*files_per_piece)
        out.append(piece)
    return out

class XHYbbWW:
    def __init__(self, inputfile, ijob=1, njobs=1):
        #inputfile format is '(raw_nano or snapshots)/SETNAME_YEAR(maybe _+something else).txt'
        infiles = SplitUp(inputfile,njobs)[ijob-1]
        self.setname = inputfile.split('/')[1].split('_')[0]
        self.year = inputfile.split('/')[1].split('_')[1].split('.')[0]
        self.ijob = ijob
        self.njobs = njobs

        self.a = analyzer(infiles)

        if 'Data' in inputfile:
          # self.a = analyzer(infiles)
            self.a.isData = True
        else:
          # MY=self.setname.split('-')[2]
          # self.a = analyzer(infiles,multiSampleStr=MY)
          # self.a.Cut('CorrectMass', 'GenModel_YMass_{} == 1'.format(MY)) 
            self.a.isData = False
 
	# Cuts from config file
        self.config = OpenJSON('XHYbbWWconfig.json')
        self.cuts = self.config['CUTS']

	# triggers
        self.trigs = self.config['TRIGS']

    CompileCpp('HWWmodules.cc') 

    def AddCutflowColumn(self, var, varName):
        '''
        for future reference:
        https://root-forum.cern.ch/t/rdataframe-define-column-of-same-constant-value/34851
        '''
        print('Adding cutflow information...\n\t{}\t{}'.format(varName, var))
        self.a.Define('{}'.format(varName), str(var))
 
    def getNweighted(self):
        if not self.a.isData:
            return self.a.DataFrame.Sum("genWeight").GetValue()
        else:
            return self.a.DataFrame.Count().GetValue()

    def BasicKinematics(self,make_ttCR=False): #nothing super selective, but should cut out some obviously uninteresting events
        #make separate snapshots which will be used to create ttbar control region for measurement of systematics
        ttCR = 'true' if make_ttCR else 'false' #C++ booleans represented with lowercase
        self.NSTART = self.getNweighted()
        self.AddCutflowColumn(self.NSTART,'NSTART')
       
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
        if self.year == '17' or self.year == '18':
            flags.append('Flag_ecalBadCalibFilter')
        MET_filters = self.a.GetFlagString(flags)       # string valid (existing in RDataFrame node) flags together w logical and
        self.a.Cut('flags', MET_filters)
        self.NFLAGS = self.getNweighted()
        self.AddCutflowColumn(self.NFLAGS, "NFLAGS")
        
        #at least one lepton passing pre-selection criteria
        self.a.Define('isLeptonPre','isLeptonPreselected(nElectron, Electron_pt, Electron_eta, Electron_miniPFRelIso_all, nMuon, Muon_pt, Muon_eta, Muon_miniPFRelIso_all)')
        self.a.Cut('isLeptonPreselected','isLeptonPre')
        self.LEPPRE=self.getNweighted()
        self.AddCutflowColumn(self.LEPPRE,'LEPPRE')

        #events pass jet preselection
        self.a.Define('isJetPre','isJetPreselected(nFatJet,FatJet_pt,FatJet_eta,FatJet_msoftdrop)')    
        self.a.Cut('isJetPreselected','isJetPre')
        self.JETPRE = self.getNweighted()
        self.AddCutflowColumn(self.JETPRE,'JETPRE') 
        
        #MET cut
        self.a.Cut('MET_cut','MET_pt > 25')
        self.METPT = self.getNweighted()
        self.AddCutflowColumn(self.METPT, "METPT")
        
        return self.a.GetActiveNode()

    def ApplyStandardCorrections(self, snapshot=False):
        # first apply corrections for snapshot phase
        if snapshot:
            if self.a.isData:
                # NOTE: LumiFilter requires the year as an integer 
                lumiFilter = ModuleWorker('LumiFilter','TIMBER/Framework/include/LumiFilter.h',[int(self.year) if 'APV' not in self.year else 16])    # defaults to perform "eval" method 
                self.a.Cut('lumiFilter',lumiFilter.GetCall(evalArgs={"lumi":"luminosityBlock"}))	       # replace lumi with luminosityBlock
                #if self.year == '18':
                    #HEM_worker = ModuleWorker('HEM_drop','TIMBER/Framework/include/HEM_drop.h',[self.setname if 'Muon' not in self.setname else self.setname[10:]])
                    #self.a.Cut('HEM','%s[0] > 0'%(HEM_worker.GetCall(evalArgs={"FatJet_eta":"FatJet_eta","FatJet_phi":"FatJet_phi"})))
            else:
                self.a = ApplyPU(self.a, 'XHYbbWWpileup.root', '20{}'.format(self.year), ULflag=True, histname='{}_{}'.format(self.setname,self.year))
                self.a.AddCorrection(Correction('Pdfweight','TIMBER/Framework/include/PDFweight_uncert.h',[self.a.lhaid],corrtype='uncert'))
                if self.year == '16' or self.year == '17' or 'APV' in self.year:
                    L1PreFiringWeight = Correction("L1PreFiringWeight","TIMBER/Framework/TopPhi_modules/BranchCorrection.cc",constructor=[],mainFunc='evalWeight',corrtype='weight',columnList=['L1PreFiringWeight_Nom','L1PreFiringWeight_Up','L1PreFiringWeight_Dn'])
                    self.a.AddCorrection(L1PreFiringWeight, evalArgs={'val':'L1PreFiringWeight_Nom','valUp':'L1PreFiringWeight_Up','valDown':'L1PreFiringWeight_Dn'})
                #elif self.year == '18':
                    #self.a.AddCorrection(Correction('HEM_drop','TIMBER/Framework/include/HEM_drop.h',[self.setname],corrtype='corr'))
            self.a = AutoJME.AutoJME(self.a, 'FatJet', '20{}'.format(self.year), self.setname if 'Muon' not in self.setname else self.setname[10:])
        # now for selection
        else:
            if not self.a.isData:
                self.a.AddCorrection(Correction('Pileup',corrtype='weight'))
                self.a.AddCorrection(Correction('Pdfweight',corrtype='uncert'))                
                if 'ttbar' in self.setname:
                    self.a.Define('GenParticle_vect','hardware::TLvector(GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)')
                    self.a.AddCorrection(
                        Correction('TptReweight','TIMBER/Framework/include/TopPt_weight.h',corrtype='weight'),
                        evalArgs={
                            'jet0':'hardware::TLvector(Wqq_pt_corr,Wqq_eta,Wqq_phi,Wqq_msoftdrop_corr)',
                            'jet1':'hardware::TLvector(Higgs_pt_corr,Higgs_eta,Higgs_phi,Higgs_msoftdrop_corr)',
                            'GenPart_vect':'GenParticle_vect',
                            'scale':'1'
                        }
                    )
                if self.year == '16' or self.year == '17' or 'APV' in self.year:
                    self.a.AddCorrection(Correction('L1PreFiringWeight',corrtype='weight'))
                elif self.year == '18':
                    self.a.AddCorrection(Correction('HEM_drop',corrtype='corr'))
        return self.a.GetActiveNode()

    
    def ApplyTopPtReweight(self):
        if 'ttbar' in self.setname:
            self.a.Define('GenParticle_vect','hardware::TLvector(GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)')
            self.a.AddCorrection(
                Correction('TptReweight','TIMBER/Framework/include/TopPt_weight.h',corrtype='weight'),
                evalArgs={
                    #'jet0':'hardware::TLvector(Wqq_pt_corr,Wqq_eta,Wqq_phi,Wqq_msoftdrop_corr)',
                    #'jet1':'hardware::TLvector(Higgs_pt_corr,Higgs_eta,Higgs_phi,Higgs_msoftdrop_corr)',
                    'jet0':'hardware::TLvector(bJet_pt,bJet_eta,bJet_phi,bJet_mass)',
                    'jet1':'hardware::TLvector(bqqJet_pt,bqqJet_eta,bqqJet_phi,bqqJet_msoftdrop)',
                    'GenPart_vect':'GenParticle_vect'
                }
            )
    
    def ApplyJMECorrections(self, variation):
	# for trigger effs
	#self.a.Define('Trijet_vect_trig','hardware::TLvector(Trijet_pt, Trijet_eta, Trijet_phi, Trijet_msoftdrop)')
	#self.a.Define('mhww_trig','hardware::InvariantMass(Trijet_vect_trig)')
	
        # JME variations - we only do this for MC
	if not self.a.isData:
	    # since H, W close enough in mass, we can treat them the same. 
	    # Higgs, W will have same pt and mass calibrations
	    # WARNING --------------------------------------------------------------------------------------------------------------
	    # IS THIS ACTUALLY TRUE?
            # IS IT, AMITAV? IS IT???
	    # ----------------------------------------------------------------------------------------------------------------------
	    pt_calibs, mass_calibs = JMEvariationStr('Higgs',variation)
	    self.a.Define('FatJet_pt_corr','hardware::MultiHadamardProduct(FatJet_pt,%s)'%pt_calibs)
	    self.a.Define('FatJet_msoftdrop_corr','hardware::MultiHadamardProduct(FatJet_msoftdrop,{})'.format(mass_calibs))
	else:
	    self.a.Define('FatJet_pt_corr','hardware::MultiHadamardProduct(FatJet_pt,{FatJet_JES_nom})')
	    self.a.Define('FatJet_msoftdrop_corr','hardware::MultiHadamardProduct(FatJet_msoftdrop,{FatJet_JES_nom})')
	return self.a.GetActiveNode()

    def ApplyPNetHbbReweight(self,wp=0.98):
        #Definition of taggers - moved here from Dijets() fuction so I can run the scale factors beforehand
        if 'XHY' in self.setname:   # only run on signal for now 
            # Get the efficiency map from the rootfile
            #effmap = 'root://cmseos.fnal.gov//store/user/mhesford/XHYbbWW_semileptonic/plots/HbbEfficiencies/{}_{}_HbbEfficiencies.root'.format(self.setname,self.year)
            #effmap = '/uscms/home/mhesford/nobackup/XHYbbWW/CMSSW_11_1_4/src/semileptonic/plots/HbbEfficiencies/{}_{}_HbbEfficiencies.root'.format(self.setname,self.year)
            effmap = '../semileptonic/plots/HbbEfficiencies/{}_{}_HbbEfficiencies.root'.format(self.setname,self.year)
            # choose the working point to delineate Fail and Pass Hbb tagging
            # the column names that will be used for the `eval()` function
            cols = [ 'FatJet_pt_corr', 'FatJet_eta', 'FatJet_particleNetMD_HbbvsQCD']
            PNetHbbWeight = Correction('PNetHbbWeight', 'TIMBER/Framework/include/PNetHbb_weight.h', constructor=[str(self.year), effmap, wp], mainFunc='eval', corrtype='weight', columnList=cols)
            self.a.AddCorrection(PNetHbbWeight, evalArgs={'FatJet_pT':'FatJet_pt_corr', 'FatJet_eta':'FatJet_eta','FatJet_PNetHbbScore':'FatJet_particleNetMD_HbbvsQCD'})

    # for trigger efficiencies
    def ApplyTrigs(self):
        #For now, just apply triggers straight up
        #Eventually, once trigger efficiencies in data measured, can just apply them as a weight correction to MC
        all_trigs = self.trigs[self.year if 'APV' not in self.year else '16']
        trigs = []
        for group in ['HADRONIC','ELECTRON','MUON']:
            for trig in all_trigs[group]:
                trigs.append(trig.encode('ascii','ignore')) #something about unicode	
        print('Trigger list: {}'.format(trigs))
        self.a.Cut('trigger',self.a.GetTriggerString(trigs))
        self.nTrigs = self.getNweighted()
        self.AddCutflowColumn(self.nTrigs, "nTrigs")
	return self.a.GetActiveNode()

    def Snapshot(self, node=None):
        startNode = self.a.GetActiveNode()
        if node == None:
            node = self.a.GetActiveNode()
        
        columns = [       
            'nJet','Jet_btagDeepFlavB','Jet_pt','Jet_eta','Jet_phi','Jet_mass','Jet_btagDeepB','Jet_jetId',
            'nFatJet','FatJet_eta','FatJet_msoftdrop','FatJet_pt','FatJet_phi','FatJet_JES_nom','FatJet_particleNetMD*', 'FatJet_rawFactor',
            'FatJet_subJetIdx1','FatJet_subJetIdx2','FatJet_hadronFlavour','FatJet_nBHadrons','FatJet_nCHadrons','SubJet*'
            'FatJet_jetId','nCorrT1METJet','CorrT1METJet_*','nMuon','Muon_*','nElectron','Electron_*','HLT_*',
            'RawMET_phi','RawMET_pt','RawMET_sumEt','ChsMET_phi','ChsMET_pt','ChsMET_sumEt','MET_phi','MET_pt','MET_sumEt',
            'genWeight','event','eventWeight','luminosityBlock','run','NSTART','NFLAGS','LEPPRE','JETPRE','METPT'
        ]
        columns.extend(['nGenJet','nSoftActivityJet','nSubJet'])
        '''
        columns = [
            'nElectron','Electron_*','nMuon','Muon_*','nCorrT1METJet','CorrT1METJet_*',
            'nJet','nJet','Jet_btagDeepFlavB','Jet_pt','Jet_eta','Jet_phi','Jet_mass','Jet_btagDeepB','Jet_jetId',
            'nFatJet','FatJet_pt','FatJet_eta','FatJet_phi','FatJet_msoftdrop','FatJet_jetId','FatJet_JES_nom','FatJet_particleNetMD*',
            'HLT_*','RawMET_phi','RawMET_pt','RawMET_sumEt','ChsMET_phi','ChsMET_pt','ChsMET_sumEt','MET_phi','MET_pt','MET_sumEt',
            'genWeight','event','eventWeight','luminosityBlock','run','NSTART','NFLAGS','LEPPRE','JETPRE','METPT'
        ]

        columns.extend(['nSoftActivityJet','nSubJet']) #stupid columns I have to add to keep this stupid snapshot function happy
        ''' 
           
        if not self.a.isData:
            columns.extend(['nGenPart','GenPart_*','GenMET_*','genWeight']) 
            columns.extend(['FatJet_JES_up','FatJet_JES_down',
                            'FatJet_JER_nom','FatJet_JER_up','FatJet_JER_down',
                            'FatJet_JMS_nom','FatJet_JMS_up','FatJet_JMS_down',
                            'FatJet_JMR_nom','FatJet_JMR_up','FatJet_JMR_down'])
            columns.extend(['Pileup__nom','Pileup__up','Pileup__down','Pdfweight__up','Pdfweight__down'])
            if self.year == '16' or self.year == '17' or 'APV' in self.year:
                columns.extend(['L1PreFiringWeight_Nom', 'L1PreFiringWeight_Up', 'L1PreFiringWeight_Dn'])       # these are the default columns in NanoAODv9
                columns.extend(['L1PreFiringWeight__nom','L1PreFiringWeight__up','L1PreFiringWeight__down'])    # these are the weight columns created by the BranchCorrection module
            elif self.year == '18':
                columns.append('HEM_drop__nom')
        
        # get ready to send out snapshot
        self.a.SetActiveNode(node)
        self.a.Snapshot(columns, 'XHYbbWWsnapshot_{}_{}_{}_{}.root'.format(self.setname,self.year,self.ijob,self.njobs),'Events', openOption='RECREATE',saveRunChain=True)
        self.a.SetActiveNode(startNode)

    def Dijets(self):

        #Pick highest pt back to back fat jets
        self.a.Define('DijetIdxs','PickDijets(FatJet_pt,FatJet_eta,FatJet_phi,FatJet_msoftdrop)')
        self.a.Cut('Dijets_exist','DijetIdxs[0] != -1 && DijetIdxs[1] != -1')

        self.nDijets = self.getNweighted()
        self.AddCutflowColumn(self.nDijets, "nDijets")

        self.a.SubCollection('Dijet','FatJet','DijetIdxs',useTake=True)

        #Define mass-decorrelated Higgs/Wqq taggers
        #Mass-decorrelated W tagger discriminant is defined by inclusion of X->cc
        #See slide 16: https://indico.cern.ch/event/809820/contributions/3632617/attachments/1970786/3278138/MassDecorrelation_ML4Jets_H_Qu.pdf
        #Right now I'm testing some older MC samples for WZ/ZZ-->llqq which do not have particleNet tagger info - run instead using deepTagMD, will treat exactly the same as particleNetMD for these initial tests
        
        #Testing new method of initial jet IDing 
        #mark one jet as Higgs and other as Wqq
        self.a.Define('WqqIdx','Dijet_particleNetMD_WqqvsQCD[0] > Dijet_particleNetMD_WqqvsQCD[1] ? 0 : 1')
        self.a.Define('HiggsIdx','WqqIdx == 0 ? 1 : 0')

        self.a.ObjectFromCollection('Higgs','Dijet','HiggsIdx')
        self.a.ObjectFromCollection('Wqq','Dijet','WqqIdx')

        return self.a.GetActiveNode()

    def KinematicLepton(self):

        self.a.Define('kinEleIdx','kinElectron(Electron_pt,Electron_eta,Electron_miniPFRelIso_all)')
        self.a.Define('kinMuIdx','kinMuon(Muon_pt,Muon_eta,Muon_miniPFRelIso_all)')  
        self.a.Cut('kinLepton_cut','kinEleIdx != -1 || kinMuIdx != -1') #at least one good lepton
        self.a.Define('LeptonType','LeptonIdx(kinEleIdx,kinMuIdx,Electron_pt,Muon_pt)') #picks higher pt signal lepton - output = 0 (lepton is electron) or 1 (lepton is muon)

        self.a.Cut('leptonQuality','LeptonType == 0? Electron_mvaFall17V2noIso_WP80[kinEleIdx] : Muon_mediumId[kinMuIdx]')

        self.nkinLep = self.getNweighted()
        self.AddCutflowColumn(self.nkinLep,'nkinLep')

        #For ease, merge some lepton columns that will be useful later (for lepton-type specific variables, use LeptonType to determine if electron or muon)
        self.a.Define('Lepton_pt','LeptonType == 1 ? Muon_pt[kinMuIdx] : Electron_pt[kinEleIdx]')
        self.a.Define('Lepton_eta','LeptonType == 1 ? Muon_eta[kinMuIdx] : Electron_eta[kinEleIdx]')
        self.a.Define('Lepton_phi','LeptonType == 1 ? Muon_phi[kinMuIdx] : Electron_phi[kinEleIdx]')
        self.a.Define('Lepton_mass','LeptonType == 1 ? Muon_mass[kinMuIdx] : Electron_mass[kinEleIdx]')
        self.a.Define('Lepton_vect','hardware::TLvector(Lepton_pt,Lepton_eta,Lepton_phi,Lepton_mass)')

        return self.a.GetActiveNode()

    def JetSelection(self,mW,mH,Wtag,Htag): #ex: mW = [60,100], Wtag = 0.8

        #mass cuts
        self.a.Cut('WqqMass','Wqq_msoftdrop > {} && Wqq_msoftdrop < {}'.format(mW[0],mW[1]))
        self.a.Cut('HiggsMass','Higgs_msoftdrop > {} && Higgs_msoftdrop < {}'.format(mH[0],mH[1]))

        self.JETMASS = self.getNweighted()
        self.AddCutflowColumn(self.JETMASS,'JETMASS')

        #tagging cuts
        self.a.Cut('WqqTag','Wqq_particleNetMD_WqqvsQCD > {}'.format(Wtag))
        self.a.Cut('HiggsTag','Higgs_particleNetMD_HbbvsQCD > {}'.format(Htag))

        self.JETTAG = self.getNweighted()
        self.AddCutflowColumn(self.JETTAG,'JETTAG')

        return self.a.GetActiveNode()

    def GetXsecScale(self):
        lumi = self.config['lumi{}'.format(self.year)]
        xsec = self.config['XSECS'][self.setname]
        if self.a.genEventSumw == 0:
            raise ValueError('%s %s: genEventSumw is 0'%(self.setname, self.year))
        print('Normalizing by lumi*xsec/genEventSumw:\n\t{} * {} / {} = {}'.format(lumi,xsec,self.a.genEventSumw,lumi*xsec/self.a.genEventSumw))
        return lumi*xsec/self.a.genEventSumw

    def ApplyMassCuts(self):
        # perform Higgs mass window cut, save cutflow info
	self.a.Cut('mH_{}_cut'.format('window'),'Higgs_msoftdrop > {0} && Higgs_msoftdrop < {1}'.format(*self.cuts['JET']['mh']))
	self.nHiggs = self.getNweighted()
	self.AddCutflowColumn(self.nHiggs, 'nHiggsMassCut')
	# Lead W mass window cut, cutflow
	self.a.Cut('mW_{}_cut'.format('window'),'Wqq_msoftdrop > {0} && Wqq_msoftdrop < {1}'.format(*self.cuts['JET']['mw']))
	self.nWqq = self.getNweighted()
	self.AddCutflowColumn(self.nWqq, 'nWqqMassCut')
	
        return self.a.GetActiveNode()
    
    def ApplySRorCR(self, SRorCR, tagger): #This is getting frequently updated as I play w/ different definitions
	'''
	SRorCR [str] = 'SR' or 'CR', used to generate cutflow information
	tagger [str] = name of tagger to be used 
        
        Work with two control regions to measure ttbar - one defined by Hbb tag fail and other defined by inverting delta phi cut between lepton/Higgs
	Selection criteria:
		SR: Hbb > 0.94
                    |deltaPhi| lepton,Higgs > 0.5
		ttCR_HF: Hbb < 0.94
                    b-tagged AK4 jet exists w/ pt > 25, |eta| < 2.4, |deltaR| lepton/jet < 2, pass medium WP of deepJet algorithm
                ttCR_dPhi: |deltaPhi| lepton/Higgs < 0.5
                    b-tagged AK4 jet exists w/ pt > 25, |eta| < 2.4, |deltaR| lepton/jet < 2, pass medium WP of deepJet algorithm
                Note: We want the control region to be void of signal/high in background, but we don't want it to be so different from the signal region. 
                      It should be similar enough to the signal region, but statistically disjunct.
	'''
	assert(SRorCR=='SR' or SRorCR=='ttCR_HF' or SRorCR=='ttCR_dPhi' or SRorCR=='ttCR_combined')
	# tagger values for each 
	if SRorCR == 'SR':
	    # Signal region 
	    self.a.Cut('Higgs_{}_cut_{}'.format(tagger,SRorCR), 'Higgs_{}_HbbvsQCD > {}'.format(tagger,self.cuts['JET']['particleNetMD_HbbvsQCD']))
            self.nHtag_SR = self.getNweighted()
            self.AddCutflowColumn(self.nHtag_SR, 'nHtag_SR')
            
            self.a.Cut('DeltaPhi_Lepton_Higgs_{}'.format(SRorCR),'abs(hardware::DeltaPhi(Lepton_phi,Higgs_phi)) > 0.5')
            self.nDPhiLH_SR = self.getNweighted()
            self.AddCutflowColumn(self.nDPhiLH_SR,'nDPhiLH_SR')
 
        elif SRorCR == 'ttCR_HF':
            #ttbar control region w/ Hbb tag fail
            self.a.Cut('Hbb_{}_cut_{}'.format(tagger,SRorCR), 'Higgs_{}_HbbvsQCD < {}'.format(tagger,self.cuts['JET']['particleNetMD_HbbvsQCD']))
            self.nHtag_ttCR = self.getNweighted()
            self.AddCutflowColumn(self.nHtag_ttCR, 'nHtag_ttCR')

            #self.a.Define('DeltaR_Jet_Lepton','hardware::DeltaR(Jet_vect,Lepton_vect)') #Can't talke delta R of single vector w/ list of vectors
            #self.a.SubCollection('JetB','Jet','Jet_btagDeepFlavB > {}'.format(self.cuts['JET']['btagM_{}'.format(self.year)]))
            self.a.Define('Jet_vect','hardware::TLvector(Jet_pt,Jet_eta,Jet_phi,Jet_mass)')
            self.a.Define('btag_idx','PickJetB(Jet_vect,Lepton_vect,Jet_btagDeepFlavB,{})'.format(self.cuts['JET']['btagM_{}'.format(self.year)]))
            self.a.Cut('nbtaggedJet_HF','btag_idx != -1')
            self.nJetB_ttCR_HF = self.getNweighted()
            self.AddCutflowColumn(self.nJetB_ttCR_HF, 'nJetB_ttCR_HF')
            
        elif SRorCR == 'ttCR_dPhi':
            self.a.Cut('DeltaPhi_Lepton_Higgs_{}'.format(SRorCR),'abs(hardware::DeltaPhi(Lepton_phi,Higgs_phi)) < 0.5')
            self.nDPhiLH_ttCR = self.getNweighted()
            self.AddCutflowColumn(self.nDPhiLH_ttCR,'nDPhiLH_ttCR')        

            self.a.Define('Jet_vect','hardware::TLvector(Jet_pt,Jet_eta,Jet_phi,Jet_mass)')
            self.a.Define('btag_idx','PickJetB(Jet_vect,Lepton_vect,Jet_btagDeepFlavB,{})'.format(self.cuts['JET']['btagM_{}'.format(self.year)]))
            self.a.Cut('nbtaggedJet_dPhi','btag_idx != -1')
            self.nJetB_ttCR_dPhi = self.getNweighted()
            self.AddCutflowColumn(self.nJetB_ttCR_dPhi, 'nJetB_ttCR_dPhi')

        elif SRorCR == 'ttCR_combined':
            self.a.Cut('ttCR_combined_cut','abs(hardware::DeltaPhi(Lepton_phi,Higgs_phi)) < 0.5 || Higgs_{}_HbbvsQCD < {}'.format(tagger,self.cuts['JET']['particleNetMD_HbbvsQCD']))
            self.nttCR_combined = self.getNweighted()
            self.AddCutflowColumn(self.nttCR_combined,'nttCR_combined')        

            self.a.Define('Jet_vect','hardware::TLvector(Jet_pt,Jet_eta,Jet_phi,Jet_mass)')
            self.a.Define('btag_idx','PickJetB(Jet_vect,Lepton_vect,Jet_btagDeepFlavB,{})'.format(self.cuts['JET']['btagM_{}'.format(self.year)]))
            self.a.Cut('nbtaggeJet_dPhi','btag_idx != -1')
            self.nJetB_ttCR_combined = self.getNweighted()
            self.AddCutflowColumn(self.nJetB_ttCR_combined, 'nJetB_ttCR_combined')


        return self.a.GetActiveNode()
    
    def ApplyPassFail(self, SRorCR, tagger):
        
        assert(SRorCR=='SR' or SRorCR=='ttCR')
        checkpoint = self.a.GetActiveNode()
        FP = OrderedDict()
        if SRorCR == 'SR':
            # Wqq fail + cutflow info
            FP['fail_'+SRorCR] = self.a.Cut('WqqTag_fail','Wqq_{0}_WqqvsQCD < {1}'.format(tagger, self.cuts['JET']['particleNetMD_WqqvsQCD']))
            self.nWF_SR = self.getNweighted()
            self.AddCutflowColumn(self.nWF_SR, 'WqqF_SR')

            # Wqq pass + cutflow
            self.a.SetActiveNode(checkpoint)
            FP['pass_'+SRorCR] = self.a.Cut('WqqTag_pass','Wqq_{0}_WqqvsQCD > {1}'.format(tagger, self.cuts['JET']['particleNetMD_WqqvsQCD']))
            self.nWP_SR = self.getNweighted()
            self.AddCutflowColumn(self.nWP_SR, 'WqqP_SR')
        elif SRorCR == 'ttCR':
            FP['fail_'+SRorCR] = self.a.Cut('WqqTag_fail','Wqq_{0}_WqqvsQCD < {1}'.format(tagger, self.cuts['JET']['particleNetMD_WqqvsQCD']))
            self.nWF_ttCR = self.getNweighted()
            self.AddCutflowColumn(self.nWF_ttCR, 'WqqF_ttCR')

            self.a.SetActiveNode(checkpoint)
            FP['pass_'+SRorCR] = self.a.Cut('WqqTag_pass','Wqq_{0}_WqqvsQCD > {1}'.format(tagger, self.cuts['JET']['particleNetMD_WqqvsQCD']))
            self.nWP_ttCR = self.getNweighted()
            self.AddCutflowColumn(self.nWP_ttCR, 'WqqP_ttCR')
        
        # reset state, return dict
        self.a.SetActiveNode(checkpoint)
        return FP

    def make_ttCR(self):
        #Built from ttCR snapshots - fuction to create ttbar control region for use in fitting ttbar to data
        #Presumably at this point we've already picked a signal-like lepton, so all that's left is to select jets 
        #But first, impose lepton id
        self.a.Cut('leptonQuality_ttCR','LeptonType == 0? {}[kinEleIdx] : {}[kinMuIdx]'.format(self.cuts['ELECTRON']['mva80'],self.cuts['MUON']['mediumId']))        
        self.a.Define('ttbar_jetIdxs','PickJets_ttCR(FatJet_pt,FatJet_eta,FatJet_phi,Jet_pt,Jet_eta,Jet_phi,Jet_btagDeepFlavB,{},Lepton_phi)'.format(self.cuts['JET']['btag_{}'.format(self.year)])) #use medium WP of deepJet algorithm

        self.a.Cut('ttCRJets_exist','ttbar_jetIdxs[0] != -1 && ttbar_jetIdxs[1] != -1')
        self.nJetsTTCR = self.getNweighted()
        self.AddCutflowColumn(self.nJetsTTCR, "nJetsTTCR")

        self.a.ObjectFromCollection('bqqJet','FatJet','tbar_jetIdxs[0]')
        self.a.ObjectFromCollection('bJet','Jet','tbar_jetIdxs[1]')

        #define Lorentz 4-vectors for both jets
        self.a.Define('bqqJet_vect','hardware::TLvector(bqqJet_pt_corr,bqqJet_eta,bqqJet_phi,bqqJet_msoftdrop_corr)')
        self.a.Define('bJet_vect','hardware::TLvector(bJet_pt,bJet_eta,bJet_phi,bJet_mass)')

        self.a.Define('mT','hardware::InvariantMass({Lepton_vect,MET_vect,bJet_vect})')
        self.a.Define('mTT','hardware::InvariantMass({Lepton_vect,MET_vect,bJet_vect,bqqJet_vect})')

    def get_standard_int(self,region,binsX,binsY,node):
        #The idea behind this function is to perform all the "standard cuts" for a given region (pass/fail for SR, CR, or ttCR) and produce a histogram
        #Then, we count the total # events in the histogram so that when we apply comparatively looser cuts to the MX vs MY phase space we can scale by (# events passing loose selection)/(# events passing standard selection)
        self.a.SetActiveNode(node)
        #first perform all universal cuts
        self.a.Define('isLeptonPre_{}'.format(region),'isLeptonPreselected(nElectron, Electron_pt, Electron_eta, Electron_miniPFRelIso_all, nMuon, Muon_pt, Muon_eta, Muon_miniPFRelIso_all)')
        self.a.Cut('isLeptonPreselected_{}'.format(region),'isLeptonPre_{}'.format(region))
        self.a.Define('FatJetHM_msoftdrop_{}'.format(region),'FatJet_msoftdrop > 50')
        self.a.Cut('nFatJet_HighMass_{}'.format(region),'FatJetHM_msoftdrop_{}.size() > 1'.format(region))
        self.a.Cut('Jet_pt_{}'.format(region),'FatJet_pt[0] > 300 && FatJet_pt[1] > 200')
        self.a.Cut('MET_cut_{}'.format(region),'MET_pt > 25')
        self.a.Cut('Lepton_pt_{}'.format(region),'Lepton_pt > 25')
        self.a.Cut('mH_{}'.format(region),'Higgs_msoftdrop > 100 && Higgs_msoftdrop < 145')
        self.a.Cut('mW_{}'.format(region),'Wqq_msoftdrop > 65 && Wqq_msoftdrop < 100')

        if region == 'pass_SR':
            self.a.Cut('Wtag_passSR','Wqq_particleNetMD_WqqvsQCD > 0.8')
            self.a.Cut('lepIso_passSR','LeptonType == 0? Electron_miniPFRelIso_all[kinEleIdx] < 0.2 : Muon_miniPFRelIso_all[kinMuIdx] < 0.2')
            self.a.Cut('lepQ_passSR','LeptonType == 0? Electron_mvaFall17V2noIso_WP80[kinEleIdx] : Muon_mediumId[kinMuIdx]')
            self.a.Cut('Htag_passSR','Higgs_particleNetMD_HbbvsQCD > 0.94')
        elif region == 'pass_CR':
            self.a.Cut('Wtag_passCR','Wqq_particleNetMD_WqqvsQCD > 0.4 && Wqq_particleNetMD_WqqvsQCD < 0.8')
            self.a.Cut('lepIso_passCR','LeptonType == 0? Electron_miniPFRelIso_all[kinEleIdx] > 0.2 && Electron_miniPFRelIso_all[kinEleIdx] < 0.5 : Muon_miniPFRelIso_all[kinMuIdx] > 0.2 && Muon_miniPFRelIso_all[kinMuIdx] < 0.5')
            self.a.Cut('lepQ_passCR','LeptonType == 0? Electron_mvaFall17V2noIso_WP80[kinEleIdx] : Muon_mediumId[kinMuIdx]')
            self.a.Cut('Htag_passCR','Higgs_particleNetMD_HbbvsQCD > 0.94')
        elif region == 'fail_SR':
            self.a.Cut('Wtag_failSR','Wqq_particleNetMD_WqqvsQCD > 0.8')
            self.a.Cut('lepIso_failSR','LeptonType == 0? Electron_miniPFRelIso_all[kinEleIdx] < 0.2 : Muon_miniPFRelIso_all[kinMuIdx] < 0.2')
            self.a.Cut('lepQ_failSR','LeptonType == 0? Electron_mvaFall17V2noIso_WP80[kinEleIdx] : Muon_mediumId[kinMuIdx]')
            self.a.Cut('Htag_failSR','Higgs_particleNetMD_HbbvsQCD < 0.94')
        elif region == 'fail_CR':
            self.a.Cut('Wtag_failCR','Wqq_particleNetMD_WqqvsQCD > 0.4 && Wqq_particleNetMD_WqqvsQCD < 0.8')
            self.a.Cut('lepIso_failCR','LeptonType == 0? Electron_miniPFRelIso_all[kinEleIdx] > 0.2 && Electron_miniPFRelIso_all[kinEleIdx] < 0.5 : Muon_miniPFRelIso_all[kinMuIdx] > 0.2 && Muon_miniPFRelIso_all[kinMuIdx] < 0.5')
            self.a.Cut('lepQ_failCR','LeptonType == 0? Electron_mvaFall17V2noIso_WP80[kinEleIdx] : Muon_mediumId[kinMuIdx]')
            self.a.Cut('Htag_failCR','Higgs_particleNetMD_HbbvsQCD < 0.94')
        elif region == 'pass_ttCR':
            self.a.Cut('Wtag_passttCR','Wqq_particleNetMD_WqqvsQCD < 0.8')
            self.a.Cut('lepIso_passttCR','LeptonType == 0? Electron_miniPFRelIso_all[kinEleIdx] < 0.2 : Muon_miniPFRelIso_all[kinMuIdx] < 0.2')
            self.a.Cut('lepQ_passttCR','LeptonType == 0? Electron_mvaFall17V2noIso_WP80[kinEleIdx] : Muon_mediumId[kinMuIdx]')
            self.a.Define('Jet_btag_passttCR','Jet_btagDeepFlavB > 0.8')
            self.a.Cut('nbtaggedJet_passttCR','Jet_btag_passttCR.size() > 0')
            self.a.Cut('Htag_passttCR','Higgs_particleNetMD_HbbvsQCD > 0.94')
        else: #region == 'fail_ttCR'
            self.a.Cut('Wtag_failttCR','Wqq_particleNetMD_WqqvsQCD < 0.8')
            self.a.Cut('lepIso_failttCR','LeptonType == 0? Electron_miniPFRelIso_all[kinEleIdx] < 0.2 : Muon_miniPFRelIso_all[kinMuIdx] < 0.2')
            self.a.Cut('lepQ_failttCR','LeptonType == 0? Electron_mvaFall17V2noIso_WP80[kinEleIdx] : Muon_mediumId[kinMuIdx]')
            self.a.Define('Jet_btag_failttCR','Jet_btagDeepFlavB > 0.8')
            self.a.Cut('nbtaggedJet_failttCR','Jet_btag_failttCR.size() > 0')
            self.a.Cut('Htag_failttCR','Higgs_particleNetMD_HbbvsQCD < 0.94')

        syst_hists = self.a.MakeTemplateHistos(ROOT.TH2F('MXvMY_{}'.format(region),'temp',binsX[0],binsX[1],binsX[2],binsY[0],binsY[1],binsY[2]),['mhwlv','mwlv'])
        stdInts = {}
        for name in syst_hists.keys():
            hist=syst_hists[name]
            hist.SetDirectory(0)
            stdInts[name] = hist.Integral()
        
        self.a.SetActiveNode(node)
        return stdInts

# for use in selection - essentially just creates combinations of all the JME variations
def JMEvariationStr(p, variation):
    base_calibs = ['FatJet_JES_nom','FatJet_JER_nom','FatJet_JMS_nom','FatJet_JMR_nom']
    variationType = variation.split('_')[0]
    pt_calib_vect = '{'
    mass_calib_vect = '{'
    for c in base_calibs:
	if 'JM' in c and p != 'Top':	# WARNING - might need to change this if we treat W, H differently for mass and pt calibrations  
	    mass_calib_vect += '%s,'%('FatJet_'+variation if variationType in c else c)
	elif 'JE' in c:
	    pt_calib_vect += '%s,'%('FatJet_'+variation if variationType in c else c)
	    mass_calib_vect += '%s,'%('FatJet_'+variation if variationType in c else c)
    pt_calib_vect = pt_calib_vect[:-1]+'}'
    mass_calib_vect = mass_calib_vect[:-1]+'}'
    return pt_calib_vect, mass_calib_vect
