import ROOT
from TIMBER.Analyzer import Correction, CutGroup, VarGroup, ModuleWorker, analyzer
from TIMBER.Tools.Common import CompileCpp, OpenJSON
from TIMBER.Tools.AutoPU import ApplyPU
from JMEvalsOnly import JMEvalsOnly
import TIMBER.Tools.AutoJME as AutoJME

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
    def __init__(self, inputfile, ijob, njobs):
        # format: (raw_nano/data_type/setname_year.txt
        infiles = SplitUp(inputfile,njobs)[ijob-1]
        self.setname=inputfile.split('/')[-1].split('_')[0]
        self.year=inputfile.split('/')[-1].split('_')[-1].split('.')[0]
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
 
        self.ijob = ijob
        self.njobs = njobs
     
	# Cuts from config file
        self.config = OpenJSON('XHYbbWWconfig.json')
        self.cuts = self.config['CUTS']

	# triggers
        self.trigs = {
            16:['HLT_PFHT800','HLT_PFHT900'],
            17:["HLT_PFHT1050","HLT_AK8PFJet500","HLT_AK8PFHT750_TrimMass50","HLT_AK8PFHT800_TrimMass50","HLT_AK8PFJet400_TrimMass30"],
            18:['HLT_AK8PFJet400_TrimMass30','HLT_AK8PFHT850_TrimMass50','HLT_PFHT1050']
        }

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

    def BasicKinematics(self):

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

        #first ensure enough leptons/jets exist to avoid referencing indexes which don't exist
        self.a.Cut('nLepton','nElectron > 0 || nMuon > 0')
        self.NLEPTON = self.getNweighted()
        self.AddCutflowColumn(self.NLEPTON,'NLEPTON')

        #at least one high pt lepton
        self.a.Cut('Lepton_pt','nElectron == 0 ? Muon_pt[0] > 25 : nMuon == 0 ? Electron_pt[0] > 25 : (Electron_pt[0] > 25 || Muon_pt[0] > 25)')
        self.LEPPT=self.getNweighted()
        self.AddCutflowColumn(self.LEPPT,"LEPPT")

        #at least one isolated lepton
        self.a.Define('IsoElectron','Electron_miniPFRelIso_all < 1')
        self.a.Define('IsoMuon','Muon_miniPFRelIso_all < 1') 
        self.a.Cut('nIsoLepton','IsoElectron.size() > 0 || IsoMuon.size() > 0')

        self.LEPISO=self.getNweighted()
        self.AddCutflowColumn(self.LEPISO,"LEPISO")

        #at least two fat jets
        self.a.Cut('nFatJet','nFatJet > 1')
        self.NJET = self.getNweighted()
        self.AddCutflowColumn(self.NJET,'NJET')

        #at least two fat jet w/ mass > 50
        self.a.Define('FatJetHM_msoftdrop','FatJet_msoftdrop > 50')
        self.a.Cut('nFatJet_HighMass','FatJetHM_msoftdrop.size() > 1')
        self.JETMASS = self.getNweighted()
        self.AddCutflowColumn(self.JETMASS,'JETMASS')

        #high pt fat jets
        self.a.Cut('Jet_pt','FatJet_pt[0] > 400 && FatJet_pt[1] > 200')
        self.JETPT=self.getNweighted()
        self.AddCutflowColumn(self.JETPT,"JETPT")

        #jet eta
        self.a.Cut('eta_cut', 'abs(FatJet_eta[0]) < 2.4 && abs(FatJet_eta[1]) < 2.4') # no super forward jets
	self.JETETA = self.getNweighted()
	self.AddCutflowColumn(self.JETETA, "JETETA")

 # corrections - used in both snapshots and selection
    def ApplyStandardCorrections(self, snapshot=False):
	# first apply corrections for snapshot phase
	if snapshot:
	    if self.a.isData:
		# NOTE: LumiFilter requires the year as an integer 
		lumiFilter = ModuleWorker('LumiFilter','TIMBER/Framework/include/LumiFilter.h',[int(self.year) if 'APV' not in self.year else 16])    # defaults to perform "eval" method 
		self.a.Cut('lumiFilter',lumiFilter.GetCall(evalArgs={"lumi":"luminosityBlock"}))	       # replace lumi with luminosityBlock
		if self.year == '18':
		    HEM_worker = ModuleWorker('HEM_drop','TIMBER/Framework/include/HEM_drop.h',[self.setname if 'Muon' not in self.setname else self.setname[10:]])
		    self.a.Cut('HEM','%s[0] > 0'%(HEM_worker.GetCall(evalArgs={"FatJet_eta":"Trijet_eta","FatJet_phi":"Trijet_phi"})))
	    else:
		self.a = ApplyPU(self.a, 'XHYbbWWpileup.root', '20{}'.format(self.year), ULflag=True, histname='{}_{}'.format(self.setname,self.year))
		if self.year == '16' or self.year == '17' or 'APV' in self.year:
		    #self.a.AddCorrection(Correction("Prefire","TIMBER/Framework/include/Prefire_weight.h",[self.year],corrtype='weight'))
		    L1PreFiringWeight = Correction("L1PreFiringWeight","TIMBER/Framework/TopPhi_modules/BranchCorrection.cc",constructor=[],mainFunc='evalWeight',corrtype='weight',columnList=['L1PreFiringWeight_Nom','L1PreFiringWeight_Up','L1PreFiringWeight_Dn'])
		    self.a.AddCorrection(L1PreFiringWeight, evalArgs={'val':'L1PreFiringWeight_Nom','valUp':'L1PreFiringWeight_Up','valDown':'L1PreFiringWeight_Dn'})
		elif self.year == '18':
		    self.a.AddCorrection(Correction('HEM_drop','TIMBER/Framework/include/HEM_drop.h',[self.setname],corrtype='corr'))

	    self.a = AutoJME.AutoJME(self.a,'FatJet', '20{}'.format(self.year), self.setname if 'Muon' not in self.setname else self.setname[10:])
	    self.a.MakeWeightCols(extraNominal='genWeight' if not self.a.isData else '')

	# now for selection
	else:
	    if not self.a.isData:
		self.a.AddCorrection(Correction('Pileup',corrtype='weight'))
		self.a.AddCorrection(Correction('Pdfweight',corrtype='uncert'))
                if self.year == '16' or self.year == '17' or 'APV' in self.year:
                    #self.a.AddCorrection(Correction('Prefire',corrtype='weight'))
		    # Instead, instantiate ModuleWorker to handle the C++ code via clang. This uses the branches already existing in NanoAODv9
		    self.a.AddCorrection(Correction('L1PreFiringWeight',corrtype='weight'))
                elif self.year == '18':
                    self.a.AddCorrection(Correction('HEM_drop',corrtype='corr'))
                #if 'ttbar' in self.setname:
                    #self.a.AddCorrection(Correction('TptReweight',corrtype='weight')
	return self.a.GetActiveNode()

    def Snapshot(self, node=None):
        startNode = self.a.GetActiveNode()
        if node == None:
            node = self.a.GetActiveNode()
        
        
        columns = [
        'nJet','nFatJet','FatJet_JES_nom','FatJet_eta','FatJet_msoftdrop','FatJet_pt','FatJet_phi',
        'FatJet_deepTagMD_HbbvsQCD', 'FatJet_deepTagMD_ZHbbvsQCD',
        'FatJet_deepTagMD_WvsQCD', 'FatJet_deepTag_TvsQCD', 'FatJet_particleNet_HbbvsQCD',
        'FatJet_particleNet_TvsQCD', 'FatJet_particleNetMD.*', 'FatJet_rawFactor', 'FatJet_tau*',
        'FatJet_jetId', 'nFatJet','FatJet_particleNetMD_Xqq', 'FatJet_particleNetMD_Xcc', 'FatJet_particleNet_QCD',
        'FatJet_particleNet_WvsQCD','nElectron','Electron_eta','Electron_pt','Electron_mass','Electron_phi','Electron_miniPFRelIso_all',
        'Electron_dxy','Electron_dz','Electron_sip3d','Electron_hoe','Electron_eInvMinusPInv','Electron_deltaEtaSC','Electron_sieie',
        'Electron_scEtOverPt','nMuon','Muon_eta','Muon_pt','Muon_mass','Muon_phi','Muon_miniPFRelIso_all','Muon_dxy','Muon_dz', 'Muon_sip3d',
        'Muon_mediumPromptId','Muon_segmentComp','Muon_isGlobal','HLT_PFHT.*', 'HLT_PFJet.*', 'HLT_AK8.*', 'HLT_Mu50',
        'HLT_TkMu50','HLT_TkMu100','HLT_IsoMu24','HLT_IsoMu27','HLT_Ele27_WPTight_Gsf','HLT_Ele35_WPTight_Gsf',
        'HLT_Ele32_WPTight_Gsf','HLT_Photon175','HLT_Photon200','HLT_Ele115_CaloIdVT_GsfTrkIdT','HLT_OldMu100',
        'event', 'eventWeight', 'luminosityBlock', 'run','NSTART','NFLAGS','NLEPTON','LEPPT','LEPISO','NJET','JETMASS','JETPT','JETETA']
        

        # append to columns list if not Data
        if not self.a.isData:
            columns.extend(['GenPart_.*', 'nGenPart', 'genWeight', 'GenModel*'])
            columns.extend(['FatJet_JES_up','FatJet_JES_down',
                            'FatJet_JER_nom','FatJet_JER_up','FatJet_JER_down',
                            'FatJet_JMS_nom','FatJet_JMS_up','FatJet_JMS_down',
                            'FatJet_JMR_nom','FatJet_JMR_up','FatJet_JMR_down'])
            columns.extend(['Pileup__nom','Pileup__up','Pileup__down','Pdfweight__nom','Pdfweight__up','Pdfweight__down'])
            if self.year == '16' or self.year == '17' or 'APV' in self.year:
		#columns.extend(['Prefire__nom','Prefire__up','Prefire__down'])
                columns.extend(['L1PreFiringWeight_Nom', 'L1PreFiringWeight_Up', 'L1PreFiringWeight_Dn'])	# these are the default columns in NanoAODv9
                columns.extend(['L1PreFiringWeight__nom','L1PreFiringWeight__up','L1PreFiringWeight__down'])    # these are the weight columns created by the BranchCorrection module
            elif self.year == '18':
                columns.append('HEM_drop__nom')
            
        # get ready to send out snapshot
        self.a.SetActiveNode(node)
        self.a.Snapshot(columns, 'XHYbbWWsnapshot_{}_{}_{}_{}.root'.format(self.setname,self.year,self.ijob,self.njobs),'Events', openOption='RECREATE',saveRunChain=True)
        self.a.SetActiveNode(startNode)

    def Dijets(self):

        #Pick highest pt back to back fat jets
        self.a.Define('DijetIdxs','PickDijets(FatJet_eta,FatJet_phi,FatJet_msoftdrop)')
        self.a.Cut('Dijets_exist','DijetIdxs[0] != -1 && DijetIdxs[1] != -1')

        self.NDIJETS = self.getNweighted()
        self.AddCutflowColumn(self.NDIJETS, "NDIJETS")

        self.a.SubCollection('Dijet','FatJet','DijetIdxs',useTake=True)

        #Define mass-decorrelated Higgs/Wqq taggers
        self.a.Define('Dijet_particleNetMD_HbbvsQCD','Dijet_particleNetMD_Xbb/(Dijet_particleNetMD_Xbb+Dijet_particleNetMD_QCD)')
        self.a.Define('Dijet_particleNetMD_WqqvsQCD','Dijet_particleNetMD_Xqq/(Dijet_particleNetMD_Xqq+Dijet_particleNetMD_QCD)')

        #mark one jet as Higgs and other as Wqq
        self.a.Define('HiggsIdx','Dijet_particleNetMD_HbbvsQCD[0] > Dijet_particleNetMD_HbbvsQCD[1] ? 0 : 1')
        self.a.Define('WqqIdx','HiggsIdx == 0 ? 1 : 0')

        self.a.ObjectFromCollection('Higgs','Dijet','HiggsIdx')
        self.a.ObjectFromCollection('Wqq','Dijet','WqqIdx')

        return self.a.GetActiveNode()

    def SignalLepton(self):

        self.a.Define('ElectronIdx','SignalLepton(Electron_pt,Electron_eta,Electron_phi,Electron_miniPFRelIso_all,Higgs_phi,Wqq_phi)')
        self.a.Define('MuonIdx','SignalLepton(Muon_pt,Muon_eta,Muon_phi,Muon_miniPFRelIso_all,Higgs_phi,Wqq_phi)')  
        self.a.Define('LeptonIdxs','LeptonIdx(ElectronIdx,MuonIdx,Electron_pt,Muon_pt)') #output = {electronIdx,muonIdx}
        self.a.Cut('Lepton_cut','LeptonIdxs[0] != -1 || LeptonIdxs[1] != -1') #at least one good lepton

        self.SIGLEP = self.getNweighted()
        self.AddCutflowColumn(self.SIGLEP,'SIGLEP')

        #For ease, merge some lepton columns that will be useful later (for lepton-type specific variables, use LeptonIdxs to determine lepton type)
        self.a.Define('Lepton_pt','LeptonIdxs[0] == -1 ? Muon_pt[LeptonIdxs[1]] : Electron_pt[LeptonIdxs[0]]')
        self.a.Define('Lepton_eta','LeptonIdxs[0] == -1 ? Muon_eta[LeptonIdxs[1]] : Electron_eta[LeptonIdxs[0]]')
        self.a.Define('Lepton_phi','LeptonIdxs[0] == -1 ? Muon_phi[LeptonIdxs[1]] : Electron_phi[LeptonIdxs[0]]')
        self.a.Define('Lepton_mass','LeptonIdxs[0] == -1 ? Muon_mass[LeptonIdxs[1]] : Electron_mass[LeptonIdxs[0]]')

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
        lumi = self.config['lumi{}'.format(self.year if 'APV' not in self.year else 16)]
        xsec = self.config['XSECS'][self.setname]
        if self.a.genEventSumw == 0:
            raise ValueError('%s %s: genEventSumw is 0'%(self.setname, self.year))
        return lumi*xsec/self.a.genEventSumw

    def ApplyMassCuts(self):
	# perform Higgs mass window cut, save cutflow info
	self.a.Cut('mH_{}_cut'.format('window'),'Higgs_msoftdrop > {0} && Higgs_msoftdrop < {1}'.format(*self.cuts['mh']))
	self.nHiggs = self.getNweighted()
	self.AddCutflowColumn(self.nHiggs, 'nHiggsMassCut')
	# Lead W mass window cut, cutflow
	self.a.Cut('mW_{}_cut'.format('window'),'Wqq_msoftdrop > {0} && Wqq_msoftdrop < {1}'.format(*self.cuts['mw']))
	self.nLeadW = self.getNweighted()
	self.AddCutflowColumn(self.nLeadW, 'nWqqMassCut')
	
        return self.a.GetActiveNode()

    def ApplyWTag(self, SRorCR, tagger):
	'''
	SRorCR [str] = 'SR' or 'CR', used to generate cutflow information
	tagger [str] = name of tagger to be used 
	W tagging criteria:
		SR: W > 0.8
		CR: W < 0.8
	'''
	assert(SRorCR=='SR' or SRorCR=='CR')
	# tagger values for each 
        Wtag = 0.8
	if SRorCR == 'SR':
	    # Signal region
	    self.a.Cut('Wqq_{}_cut_{}'.format(tagger,SRorCR), 'Wqq_{0}_WqqvsQCD > {1}'.format(tagger,Wtag))
	else:
	    # Control Region
	    self.a.Cut('Wqq_{}_cut_{}'.format(tagger,SRorCR), 'Wqq_{0}_WqqvsQCD < {1}'.format(tagger,Wtag))
	# save cutflow info
	self.nWTag = self.getNweighted()
	self.AddCutflowColumn(self.nWTag, 'nWTag_{}'.format(SRorCR))
	return self.a.GetActiveNode()

    def ApplyHiggsTag(self, SRorCR, tagger):
	'''
	Fail:	H < 0.94
	Pass: 	H > 0.94
	'''
	assert(SRorCR=='SR' or SRorCR=='CR')
	checkpoint = self.a.GetActiveNode()
	Htag = 0.94
	FLP = {}
	# Higgs fail + cutflow info
	FLP['fail'] = self.a.Cut('HbbTag_fail','Hbb_{0}_HbbvsQCD < {1}'.format(tagger, Htag))
	if SRorCR=='SR':
	    self.nHF_SR = self.getNweighted()
	    self.AddCutflowColumn(self.nHF_SR, 'higgsF_SR')
	else:
	    self.nHF_CR = self.getNweighted()
	    self.AddCutflowColumn(self.nHF_CR, 'higgsF_CR')

	# Higgs Pass + cutflow
	self.a.SetActiveNode(checkpoint)
	FLP['pass'] = self.a.Cut('HbbTag_pass','Hbb_{0}_HbbvsQCD > {1}'.format(tagger, Htag))
	if SRorCR == 'SR':
	    self.nHP_SR = self.getNweighted()
	    self.AddCutflowColumn(self.nHP_SR, 'higgsP_SR')
	else:
	    self.nHP_CR = self.getNweighted()
	    self.AddCutflowColumn(self.nHP_CR, 'higgsP_CR')

	# reset state, return dict
	self.a.SetActiveNode(checkpoint)
	return FLP
