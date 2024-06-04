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
        print(infiles)
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

    CompileCpp('HWWmodules.cc') # Compile c++ helper functions (for use throughout analysis)

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

    def MakePNetCR(self):
        '''
        Creates control region used for calculating particleNet mistagging SF's in ttbar (ttbar enriched region)
        Requirements: 
            Leptonic top decay signatures:
            -1 electron(muon) w/ pT > 25 GeV, |eta| < 2.5(2.4), mini_iso < 0.1, pass tight ID
            -1 b-tagged AK4 jet w/ pT > 25, |eta| < 2.4, pass jetId 
            -MET > 30
            Search for (what is hopefully) a top jet
            -AK8 jet w/ pT > 200, |eta| < 2.4
                -DeltaR(AK8, lepton) > 0.8
                -DeltaR(AK8, AK4) > 0.8
        '''
        self.KinematicLepton() #should cover all the lepton criteria

        #Now pick AK4 jet            
        self.a.Define('Jet_vect','hardware::TLvector(Jet_pt,Jet_eta,Jet_phi,Jet_mass)')
        self.a.Define('btag_idx','PickJetB(Jet_vect,Lepton_vect,Jet_btagDeepFlavB,{})'.format(self.cuts['JET']['btagM_{}'.format(self.year)]))
        self.a.Cut('nbtaggedJet','btag_idx != -1')
        self.nBtag = self.getNweighted()
        self.AddCutflowColumn(self.nBtag, 'nBtag')
        self.a.ObjectFromCollection('BJet','Jet','btag_idx')

        #MET
        self.a.Cut('ChsMET','ChsMET_pt > 30')
        self.nMET = self.getNweighted()
        self.AddCutflowColumn(self.nMET, 'nMET')
                        
        #Now look for AK8 jet
        self.a.Define('FatJet_vect','hardware::TLvector(FatJet_pt,FatJet_eta,FatJet_phi,FatJet_msoftdrop)')
        self.a.Define('ak8jet_idx','PickAK8(FatJet_vect,BJet_vect,Lepton_vect)')
        self.a.Cut('nak8jet','ak8jet_idx != -1')
        self.nProbe = self.getNweighted()
        self.AddCutflowColumn(self.nProbe, 'nProbe')
        self.a.ObjectFromCollection('ProbeJet','FatJet','ak8jet_idx')  

        #Define taggers used
        self.a.Define('ProbeJet_particleNetMD_HbbvsQCD','ProbeJet_particleNetMD_Xbb/(ProbeJet_particleNetMD_Xbb+ProbeJet_particleNetMD_QCD)')
        self.a.Define('ProbeJet_particleNetMD_WqqvsQCD','(ProbeJet_particleNetMD_Xqq+ProbeJet_particleNetMD_Xcc)/(ProbeJet_particleNetMD_Xqq+ProbeJet_particleNetMD_Xcc+ProbeJet_particleNetMD_QCD)')

        return self.a.GetActiveNode()

    def ApplyStandardCorrections(self, snapshot=False):
        # first apply corrections for snapshot phase
        if snapshot:
            if self.a.isData:
                # NOTE: LumiFilter requires the year as an integer 
                lumiFilter = ModuleWorker('LumiFilter','TIMBER/Framework/include/LumiFilter.h',[int(self.year) if 'APV' not in self.year else 16])    # defaults to perform "eval" method 
                self.a.Cut('lumiFilter',lumiFilter.GetCall(evalArgs={"lumi":"luminosityBlock"}))	       # replace lumi with luminosityBlock
                if self.year == '18':
                    HEM_worker = ModuleWorker('HEM_drop','TIMBER/Framework/include/HEM_drop.h',[self.setname if 'Data' not in self.setname else self.setname[-5:]])
                    self.a.Cut('HEM','%s[0] > 0'%(HEM_worker.GetCall(evalArgs={"FatJet_eta":"FatJet_eta","FatJet_phi":"FatJet_phi"})))
            else:
                self.a = ApplyPU(self.a, 'XHYbbWWpileup.root', '20{}'.format(self.year), ULflag=True, histname='{}_{}'.format(self.setname,self.year))
                if self.a.lhaid != -1: #some backgrounds (diboson) don't have LHA ID in MC file
                    self.a.AddCorrection(Correction('Pdfweight','TIMBER/Framework/include/PDFweight_uncert.h',[self.a.lhaid],corrtype='uncert'))
                if self.year == '16' or self.year == '17' or 'APV' in self.year:
                    L1PreFiringWeight = Correction("L1PreFiringWeight","TIMBER/Framework/Zbb_modules/BranchCorrection.cc",constructor=[],mainFunc='evalWeight',corrtype='weight',columnList=['L1PreFiringWeight_Nom','L1PreFiringWeight_Up','L1PreFiringWeight_Dn'])
                    self.a.AddCorrection(L1PreFiringWeight, evalArgs={'val':'L1PreFiringWeight_Nom','valUp':'L1PreFiringWeight_Up','valDown':'L1PreFiringWeight_Dn'})
                elif self.year == '18':
                    self.a.AddCorrection(Correction('HEM_drop','TIMBER/Framework/include/HEM_drop.h',[self.setname],corrtype='corr'))

                # Parton shower weights 
                #	- https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopSystematics#Parton_shower_uncertainties
                #	- "Default" variation: https://twiki.cern.ch/twiki/bin/view/CMS/HowToPDF#Which_set_of_weights_to_use
                #	- https://github.com/mroguljic/Hgamma/blob/409622121e8ab28bc1072c6d8981162baf46aebc/templateMaker.py#L210
                self.a.Define("ISR__up","PSWeight[2]")
                self.a.Define("ISR__down","PSWeight[0]")
                self.a.Define("FSR__up","PSWeight[3]")
                self.a.Define("FSR__down","PSWeight[1]")
                genWCorr = Correction('genW','TIMBER/Framework/Zbb_modules/BranchCorrection.cc',corrtype='corr',mainFunc='evalCorrection') # workaround so we can have multiple BCs
                self.a.AddCorrection(genWCorr, evalArgs={'val':'genWeight'})
                #ISRcorr = Correction('ISRunc', 'TIMBER/Framework/TopPhi_modules/BranchCorrection.cc', mainFunc='evalUncert', corrtype='uncert')
                #FSRcorr = Correction('FSRunc', 'TIMBER/Framework/TopPhi_modules/BranchCorrection.cc', mainFunc='evalUncert', corrtype='uncert')
                ISRcorr = genWCorr.Clone("ISRunc",newMainFunc="evalUncert",newType="uncert")
                FSRcorr = genWCorr.Clone("FSRunc",newMainFunc="evalUncert",newType="uncert")
                self.a.AddCorrection(ISRcorr, evalArgs={'valUp':'ISR__up','valDown':'ISR__down'})
                self.a.AddCorrection(FSRcorr, evalArgs={'valUp':'FSR__up','valDown':'FSR__down'})
                # QCD factorization and renormalization corrections (only to non-signal MC)
                if (('WW' not in self.setname) and ('WZ' not in self.setname) and ('ZZ' not in self.setname)):
                    # First instatiate a correction module for the factorization correction
                    facCorr = Correction('QCDscale_factorization','LHEScaleWeights.cc',corrtype='weight',mainFunc='evalFactorization')
                    self.a.AddCorrection(facCorr, evalArgs={'LHEScaleWeights':'LHEScaleWeight'})
                    # Now clone it and call evalRenormalization for the renormalization correction
                    renormCorr = facCorr.Clone('QCDscale_renormalization',newMainFunc='evalRenormalization',newType='weight')
                    self.a.AddCorrection(renormCorr, evalArgs={'LHEScaleWeights':'LHEScaleWeight'})
                    # Now do one for the combined correction
                    combCorr = facCorr.Clone('QCDscale_combined',newMainFunc='evalCombined',newType='weight')
                    self.a.AddCorrection(combCorr, evalArgs={'LHEScaleWeights':'LHEScaleWeight'})
                    # And finally, do one for the uncertainty
                    # See: https://indico.cern.ch/event/938672/contributions/3943718/attachments/2073936/3482265/MC_ContactReport_v3.pdf (slide 27)
                    QCDScaleUncert = facCorr.Clone('QCDscale_uncert',newMainFunc='evalUncert',newType='uncert')
                    self.a.AddCorrection(QCDScaleUncert, evalArgs={'LHEScaleWeights':'LHEScaleWeight'})

            if 'Data' not in self.setname:
                JMEname = self.setname
            elif 'DataB' in self.setname: #account for extra v1/v2 at the end of the setname
                JMEname = 'DataB'
            else:
                JMEname = self.setname[-5:]
            self.a = AutoJME.AutoJME(self.a, 'FatJet', '20{}'.format(self.year), JMEname)

        # now for selection
        else:
            if not self.a.isData:
                # I forgot to add the `genW` branch in snapshots, so just redo it here...
                # In the end it doesn't really matter, since the correction just uses genWeight.
                # One could also opt to add genWeight*GetXsecScale() in the MakeWeightCols() call as well..
                # This is EXTREMELY IMPORTANT for getting the event weighting correct
                genWCorr = Correction('genW','TIMBER/Framework/Zbb_modules/BranchCorrection.cc',corrtype='corr',mainFunc='evalCorrection')
                #self.a.AddCorrection(Correction('genW',corrtype='corr'))
                self.a.AddCorrection(genWCorr, evalArgs={'val':'genWeight'})

                self.a.AddCorrection(Correction('Pileup',corrtype='weight'))
                self.a.AddCorrection(Correction('Pdfweight',corrtype='uncert'))                
                self.a.AddCorrection(Correction('ISRunc',corrtype='uncert'))
                self.a.AddCorrection(Correction('FSRunc',corrtype='uncert'))

                # perhaps the first three should be uncert types, but because nominal = 1.0, it's functionally equivalent
                self.a.AddCorrection(Correction('QCDscale_factorization',corrtype='weight'))
                self.a.AddCorrection(Correction('QCDscale_renormalization',corrtype='weight'))
                self.a.AddCorrection(Correction('QCDscale_combined',corrtype='weight'))
                self.a.AddCorrection(Correction('QCDscale_uncert',corrtype='uncert'))

                if self.year == '16' or self.year == '17' or 'APV' in self.year:
                    self.a.AddCorrection(Correction('L1PreFiringWeight',corrtype='weight'))
                elif self.year == '18':
                    self.a.AddCorrection(Correction('HEM_drop',corrtype='corr'))
        return self.a.GetActiveNode()

    
    def ApplyTopPtReweight(self,jet0,jet1,scale = 1, isJet1AK4 = False):
        #isJet1AK4 asks if the second jet is an AK4 jet, in which case definitions of mass/momentum are different
        if 'ttbar' in self.setname:
            if isJet1AK4:
                jet1_pt = 'pt'
                jet1_mass = 'mass'
            else:
                jet1_pt = 'pt_corr'
                jet1_mass = 'msoftdrop_corr'
            self.a.Define('GenParticle_vect','hardware::TLvector(GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)')
            self.a.AddCorrection(
                Correction('TptReweight','TIMBER/Framework/include/TopPt_weight.h',corrtype='weight'),
                    evalArgs={
                        'jet0':'hardware::TLvector({0}_pt_corr,{0}_eta,{0}_phi,{0}_msoftdrop_corr)'.format(jet0),
                        'jet1':'hardware::TLvector({0}_{1},{0}_eta,{0}_phi,{0}_{2})'.format(jet1,jet1_pt,jet1_mass),
                        'GenPart_vect':'GenParticle_vect',
                        'scale':scale
                    }
             )  
   
 
    def ApplyJMECorrections(self, variation, jets = ['FatJet']):
	# for trigger effs
	#self.a.Define('Trijet_vect_trig','hardware::TLvector(Trijet_pt, Trijet_eta, Trijet_phi, Trijet_msoftdrop)')
	#self.a.Define('mhww_trig','hardware::InvariantMass(Trijet_vect_trig)')
	
        # JME variations - we only do this for MC
        for jet in jets:
            if not self.a.isData:
    	        # since H, W close enough in mass, we can treat them the same. 
    	        # Higgs, W will have same pt and mass calibrations
    	        # WARNING --------------------------------------------------------------------------------------------------------------
	        # IS THIS ACTUALLY TRUE?
                # IS IT, AMITAV? IS IT???
	        # ----------------------------------------------------------------------------------------------------------------------
                pt_calibs, mass_calibs = JMEvariationStr('Higgs',jet,variation)
                if jet == 'FatJet': #dealing with RVec quantities, need to do hadamard product
                    self.a.Define('{}_pt_corr'.format(jet),'hardware::MultiHadamardProduct({}_pt,{})'.format(jet,pt_calibs))
                    self.a.Define('{}_msoftdrop_corr'.format(jet),'hardware::MultiHadamardProduct({}_msoftdrop,{})'.format(jet,mass_calibs))
                else: #Assume dealing with a single jet and not a collection
                    self.a.Define('{}_pt_corr'.format(jet),'%s_pt%s'%(jet,pt_calibs))
                    self.a.Define('{}_msoftdrop_corr'.format(jet),'%s_msoftdrop%s'%(jet,mass_calibs))
            else:
                if jet == 'FatJet':
                    self.a.Define('{}_pt_corr'.format(jet),'hardware::MultiHadamardProduct(%s_pt,{%s_JES_nom})'%(jet,jet))
                    self.a.Define('{}_msoftdrop_corr'.format(jet),'hardware::MultiHadamardProduct(%s_msoftdrop,{%s_JES_nom})'%(jet,jet))
                else: 
                    self.a.Define('{}_pt_corr'.format(jet),'%s_pt * %s_JES_nom'%(jet,jet))
                    self.a.Define('{}_msoftdrop_corr'.format(jet),'%s_msoftdrop * %s_JES_nom'%(jet,jet))
        return self.a.GetActiveNode()

    def make_other_efficiency_map(self, tagger, cats):
        #Helper function to make efficiency map for "other" jet category - combines efficiency maps for multiple jet flavors
        print('Making combined efficiency map for tagger {} and jet flavors {}'.format(tagger,cats))
        all_hist = ROOT.TH2F('Jets all','Jets all'+';pt;eta',15,0,3000,12,-2.4,2.4)
        tagged_hist = ROOT.TH2F('Jets tagged','Jets tagged'+';pt;eta',15,0,3000,12,-2.4,2.4)        
        outfile_name = 'plots/ParticleNet_Efficiencies/{}_{}_otherJets_particleNet{}_Efficiencies.root'.format(self.setname,self.year,tagger)
        outfile = ROOT.TFile.Open(outfile_name,'RECREATE')
        for cat in cats:
            infile = ROOT.TFile.Open('plots/ParticleNet_Efficiencies/{}_{}_{}Jets_particleNet{}_Efficiencies.root'.format(self.setname,self.year,cat,tagger),'READ')
            inhist_all = infile.Get('Jets all')
            inhist_tagged = infile.Get('Jets tagged')
            all_hist.Add(inhist_all)
            tagged_hist.Add(inhist_tagged)
        ratio_hist = make_ratio_hist(tagged_hist,all_hist)
        ratio_hist.SetName('ratio')
        ratio_hist.SetTitle('ratio')
   
        outfile.cd()
        all_hist.Write()
        tagged_hist.Write()
        ratio_hist.Write() 

        return outfile_name

    def ApplyPNetReweight(self,Hbb_wp=0.98, Wqq_wp=0.80): # wp chosen to delineate Hbb pass/fail
        if 'XHY' in self.setname:
            # Signal, apply "true" Higgs/W scale factors
            # Get the efficiency map from the rootfile
            Hbb_effmap = 'plots/ParticleNet_Efficiencies/{}_{}_HiggsJets_particleNetHbb_Efficiencies.root'.format(self.setname,self.year)
            Wqq_effmap = 'plots/ParticleNet_Efficiencies/{}_{}_WJets_particleNetWqq_Efficiencies.root'.format(self.setname,self.year)
            #column names for eval function
            Hbb_cols = ['Higgs_pt_corr', 'Higgs_eta', 'Higgs_particleNetMD_HbbvsQCD', 'Higgs_jetFlavor']
            Wqq_cols = ['Wqq_pt_corr', 'Wqq_eta', 'Wqq_particleNetMD_WqqvsQCD', 'Wqq_jetFlavor']

            #Make scale factors/uncertainties
            PNetHbb_HJ = Correction('PNetHbb_HJ','TIMBER/Framework/include/PNetHbb_weight.h',constructor=[str(self.year), Hbb_effmap, Hbb_wp, 4], mainFunc='eval', corrtype='weight', columnList=Hbb_cols)
            self.a.AddCorrection(PNetHbb_HJ, evalArgs={'Higgs_pt':'Higgs_pt_corr', 'Higgs_eta':'Higgs_eta','Higgs_PNetHbbScore':'Higgs_particleNetMD_HbbvsQCD','Higgs_jetFlavor':'Higgs_jetFlavor'})
            PNetWqq_WJ = Correction('PNetWqq_WJ','TIMBER/Framework/include/PNetWqq_weight.h',constructor=[str(self.year), Wqq_effmap, Wqq_wp, 1], mainFunc='eval', corrtype='weight', columnList=Wqq_cols)
            self.a.AddCorrection(PNetWqq_WJ, evalArgs={'Wqq_pt':'Wqq_pt_corr', 'Wqq_eta':'Wqq_eta','Wqq_PNetWqqScore':'Wqq_particleNetMD_WqqvsQCD','Wqq_jetFlavor':'Wqq_jetFlavor'})

        elif 'ttbar' in self.setname: 
            # Apply top/bq scale factors for Hbb tagger, top/w/other(bq+unmerged) scale factors for w tagger
            # Note: -1 jet index in the constructor denotes "other" (a combination of multiple flavors/indicies) 
                # Indicies to consider as "other" can be changed by adding other = {idx1,idx2,etc...}
                # Default indices: {1,3} (w, unmerged) for Hbb tagger, {2,3} (bq, unmerged for Wqq tagger)  

            # Hbb scale factors
            Hbb_cols = ['Higgs_pt_corr', 'Higgs_eta', 'Higgs_particleNetMD_HbbvsQCD', 'Higgs_jetFlavor']
            Hbb_tJ_effmap = 'plots/ParticleNet_Efficiencies/{}_{}_TopJets_particleNetHbb_Efficiencies.root'.format(self.setname,self.year)
            PNetHbb_tJ = Correction('PNetHbb_tJ', 'TIMBER/Framework/include/PNetHbb_weight.h', constructor=[str(self.year),Hbb_tJ_effmap,Hbb_wp,0],mainFunc='eval', corrtype='weight', columnList=Hbb_cols)
            self.a.AddCorrection(PNetHbb_tJ, evalArgs={'Higgs_pt':'Higgs_pt_corr', 'Higgs_eta':'Higgs_eta','Higgs_PNetHbbScore':'Higgs_particleNetMD_HbbvsQCD','Higgs_jetFlavor':'Higgs_jetFlavor'})

            Hbb_bqJ_effmap = 'plots/ParticleNet_Efficiencies/{}_{}_bqJets_particleNetHbb_Efficiencies.root'.format(self.setname,self.year)
            PNetHbb_bqJ = Correction('PNetHbb_bqJ', 'TIMBER/Framework/include/PNetHbb_weight.h', constructor=[str(self.year),Hbb_bqJ_effmap,Hbb_wp,2],mainFunc='eval',corrtype='weight',columnList=Hbb_cols)
            self.a.AddCorrection(PNetHbb_bqJ, evalArgs={'Higgs_pt':'Higgs_pt_corr', 'Higgs_eta':'Higgs_eta','Higgs_PNetHbbScore':'Higgs_particleNetMD_HbbvsQCD','Higgs_jetFlavor':'Higgs_jetFlavor'})
            
            # Now do Wqq scale factors
            Wqq_cols = ['Wqq_pt_corr', 'Wqq_eta', 'Wqq_particleNetMD_WqqvsQCD', 'Wqq_jetFlavor']
            Wqq_tJ_effmap = 'plots/ParticleNet_Efficiencies/{}_{}_TopJets_particleNetWqq_Efficiencies.root'.format(self.setname,self.year)
            PNetWqq_tJ = Correction('PNetWqq_tJ', 'TIMBER/Framework/include/PNetWqq_weight.h', constructor=[str(self.year),Wqq_tJ_effmap,Wqq_wp,0], mainFunc='eval', corrtype='weight', columnList=Wqq_cols)
            self.a.AddCorrection(PNetWqq_tJ, evalArgs={'Wqq_pt':'Wqq_pt_corr', 'Wqq_eta':'Wqq_eta','Wqq_PNetWqqScore':'Wqq_particleNetMD_WqqvsQCD','Wqq_jetFlavor':'Wqq_jetFlavor'})

            Wqq_WJ_effmap = 'plots/ParticleNet_Efficiencies/{}_{}_WJets_particleNetWqq_Efficiencies.root'.format(self.setname,self.year)
            PNetWqq_WJ = Correction('PNetWqq_WJ', 'TIMBER/Framework/include/PNetWqq_weight.h', constructor=[str(self.year),Wqq_WJ_effmap,Wqq_wp,1], mainFunc='eval', corrtype='weight', columnList=Wqq_cols)
            self.a.AddCorrection(PNetWqq_WJ, evalArgs={'Wqq_pt':'Wqq_pt_corr', 'Wqq_eta':'Wqq_eta','Wqq_PNetWqqScore':'Wqq_particleNetMD_WqqvsQCD','Wqq_jetFlavor':'Wqq_jetFlavor'})

            Wqq_oJ_effmap = self.make_other_efficiency_map(tagger = 'Wqq', cats = ['bq','UM'])
            PNetWqq_oJ = Correction('PNetWqq_oJ', 'TIMBER/Framework/include/PNetWqq_weight.h', constructor=[str(self.year),Wqq_oJ_effmap,Wqq_wp,-1],mainFunc='eval', corrtype='weight', columnList=Wqq_cols)
            self.a.AddCorrection(PNetWqq_oJ, evalArgs={'Wqq_pt':'Wqq_pt_corr', 'Wqq_eta':'Wqq_eta','Wqq_PNetWqqScore':'Wqq_particleNetMD_WqqvsQCD','Wqq_jetFlavor':'Wqq_jetFlavor'})
            

    def ApplyLeptonCorrections(self):
        if not self.a.isData:
            #Corrections for electron mva ID (wp80 no isolation)
            cols = ['Lepton_eta','Lepton_pt','LeptonType']
            ele_filepath = 'corrections/electron_{}.json'.format(self.year) #path to json file - for constructor
            ElectronIDWeight = Correction('ElectronIDWeight','TIMBER/Framework/include/ElectronID_weight.h',
                constructor=[str(self.year), "wp80noiso", ele_filepath], #proper pt, eta bounds included as default argument in c++ code
                mainFunc='eval', 
                corrtype='weight', 
                columnList=cols
            )
            self.a.AddCorrection(ElectronIDWeight, evalArgs={'Lepton_eta':'Lepton_eta','Lepton_pt':'Lepton_pt','LeptonType':'LeptonType'})

            #Can use same TIMBER correction module to obtain electron reconstruction corrections
            ElectronRecoWeight = Correction('ElectronRecoWeight','TIMBER/Framework/include/ElectronID_weight.h',
                constructor=[str(self.year), "RecoAbove20", ele_filepath], 
                mainFunc='eval', 
                corrtype='weight', 
                columnList=cols
            )
            self.a.AddCorrection(ElectronRecoWeight, evalArgs={'Lepton_eta':'Lepton_eta','Lepton_pt':'Lepton_pt','LeptonType':'LeptonType'})
            
            #Apply muon ID corrections
            muoID_filepath = 'corrections/ScaleFactors_Muon_highPt_IDISO_20{}_schemaV2.json'.format(self.year)
            self.a.Define('Lepton_abseta','abs(Lepton_eta)') # Muon scale factors organized by |eta|
            muoID_cols = ['Lepton_abseta','Lepton_pt','LeptonType']
            MuonIDWeight = Correction('MuonIDWeight','TIMBER/Framework/include/MuonID_weight.h',
                constructor=["NUM_MediumID_DEN_GlobalMuonProbes", muoID_filepath], 
                mainFunc='eval', 
                corrtype='weight', 
                columnList=muoID_cols #same column list
            )
            self.a.AddCorrection(MuonIDWeight, evalArgs={'Lepton_eta':'Lepton_abseta','Lepton_pt':'Lepton_pt','LeptonType':'LeptonType'})

            #Muon reconstruction scale factors
            self.a.Define('Lepton_p','Lepton_pt*cosh(Lepton_eta)') #Muon reco scale factors indexed by p (not pt)
            muoRECO_filepath = 'corrections/ScaleFactors_Muon_highPt_RECO_20{}_schemaV2.json'.format(self.year)
            muoRECO_cols = ['Lepton_abseta','Lepton_p','LeptonType']
            MuonRecoWeight = Correction('MuonRecoWeight','TIMBER/Framework/include/MuonID_weight.h',
                constructor=["NUM_GlobalMuons_DEN_TrackerMuonProbes",muoRECO_filepath], 
                mainFunc='eval', 
                corrtype='weight', 
                columnList=muoRECO_cols
            )
            self.a.AddCorrection(MuonRecoWeight, evalArgs={'Lepton_eta':'Lepton_abseta','Lepton_pt':'Lepton_p','LeptonType':'LeptonType'})

    def GetXsecScale(self):
        lumi = self.config['lumi{}'.format(self.year)]
        xsec = self.config['XSECS'][self.setname]
        if self.a.genEventSumw == 0:
            raise ValueError('%s %s: genEventSumw is 0'%(self.setname, self.year))
        print('Normalizing by lumi*xsec/genEventSumw:\n\t{} * {} / {} = {}'.format(lumi,xsec,self.a.genEventSumw,lumi*xsec/self.a.genEventSumw))
        return lumi*xsec/self.a.genEventSumw

    # for trigger efficiencies
    def ApplyTrigs(self, corr=None):
        if self.a.isData:
            #The basic logic is that for a given primary dataset, we need to accept events which pass the corresponding triggers for that dataset, but reject events which also pass the triggers of other datasets (in order to avoid double counting). We start by accepting all events which pass SingleMuon triggers and then go down the ladder. Note that for 2018, the electron/photon triggers were combined into a single PD (EGamma), but are separate (Single Electron + Single Photon) for the other years.
            all_trigs = self.trigs[self.year if '16' not in self.year else '16']
            h_trigs = all_trigs['JETHT']
            e_trigs = all_trigs['ELECTRON'] if self.year != '18' else all_trigs['EGAMMA']
            m_trigs = all_trigs['MUON']
            if self.year != '18':
                p_trigs = all_trigs['PHOTON']

            if 'Muon' in self.setname:
                print('Muon trigger list: {}'.format(m_trigs))
                self.a.Cut('m_trigger',self.a.GetTriggerString(m_trigs))

            elif 'Electron' in self.setname or 'EGamma' in self.setname:
                print('Electron/EGamma trigger list: {}'.format(e_trigs))
                self.a.Cut('e_trigger',self.a.GetTriggerString(e_trigs))
                self.a.Cut('e_muOrthogonal','!'+self.a.GetTriggerString(m_trigs)) #remove events which pass electron AND muon triggers to avoid double-counting

            elif 'Photon' in self.setname: #for year != 2018
                print('Photon trigger list: {}'.format(p_trigs))
                self.a.Cut('p_trigger',self.a.GetTriggerString(p_trigs))
                self.a.Cut('p_muOrthogonal','!'+self.a.GetTriggerString(m_trigs))
                self.a.Cut('p_eleOrthogonal','!'+self.a.GetTriggerString(e_trigs))

            else: #JetHT data
                print('Hadronic trigger list: {}'.format(h_trigs))
                self.a.Cut('h_trigger',self.a.GetTriggerString(h_trigs))
                self.a.Cut('h_muOrthogonal','!'+self.a.GetTriggerString(m_trigs))
                self.a.Cut('h_eleOrthogonal','!'+self.a.GetTriggerString(e_trigs))
                if self.year != '18':
                    self.a.Cut('h_gamOrthogonal','!'+self.a.GetTriggerString(p_trigs))
        else:
            self.a.AddCorrection(corr, evalArgs={'xval':'mhww_msoftdrop','yval':'Lepton_abseta'})
        
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
            'FatJet_jetId','nCorrT1METJet','CorrT1METJet_*','nMuon','Muon_*','nElectron','Electron_*',
            'RawMET_phi','RawMET_pt','RawMET_sumEt','ChsMET_phi','ChsMET_pt','ChsMET_sumEt','MET_phi','MET_pt','MET_sumEt',
            'genWeight','event','eventWeight','luminosityBlock','run','NSTART','NFLAGS','LEPPRE','JETPRE','METPT'
        ]
        columns.extend(['nGenJet','nSoftActivityJet','nSubJet'])
        columns.extend(['HLT_AK8*','HLT_PFHT*','HLT_PFJet*','HLT_Iso*','HLT_Mu*','HLT_TkMu*','HLT_OldMu*','HLT_Ele*','HLT_Photon*','HLT_MET*','HLT_PFMET*'])
          
        if not self.a.isData:
            columns.extend(['nGenPart','GenPart_*','GenMET_*','genWeight']) 
            columns.extend(['FatJet_JES_up','FatJet_JES_down',
                            'FatJet_JER_nom','FatJet_JER_up','FatJet_JER_down',
                            'FatJet_JMS_nom','FatJet_JMS_up','FatJet_JMS_down',
                            'FatJet_JMR_nom','FatJet_JMR_up','FatJet_JMR_down'])
            columns.extend(['Pileup__nom','Pileup__up','Pileup__down','Pdfweight__up','Pdfweight__down','ISR__up','ISR__down','FSR__up','FSR__down'])
            if self.year == '16' or self.year == '17' or 'APV' in self.year:
                columns.extend(['L1PreFiringWeight_Nom', 'L1PreFiringWeight_Up', 'L1PreFiringWeight_Dn'])       # these are the default columns in NanoAODv9
                columns.extend(['L1PreFiringWeight__nom','L1PreFiringWeight__up','L1PreFiringWeight__down'])    # these are the weight columns created by the BranchCorrection module
            elif self.year == '18':
                columns.append('HEM_drop__nom')

            if (('WW' not in self.setname) and ('WZ' not in self.setname) and ('ZZ' not in self.setname)):
                columns.extend(['QCDscale_factorization__nom','QCDscale_factorization__up','QCDscale_factorization__down'])
                columns.extend(['QCDscale_renormalization__nom','QCDscale_renormalization__up','QCDscale_renormalization__down'])
                columns.extend(['QCDscale_combined__nom','QCDscale_combined__up','QCDscale_combined__down'])
                columns.extend(['QCDscale_uncert__up','QCDscale_uncert__down'])
        
        # get ready to send out snapshot
        self.a.SetActiveNode(node)
        self.a.Snapshot(columns, 'XHYbbWWsnapshot_{}_{}_{}_{}.root'.format(self.setname,self.year,self.ijob,self.njobs),'Events', openOption='RECREATE',saveRunChain=True)
        self.a.SetActiveNode(startNode)

    def Dijets(self):
        #Define mass-decorrelated Higgs/Wqq taggers
        #Mass-decorrelated W tagger discriminant is defined by inclusion of X->cc
        #See slide 16: https://indico.cern.ch/event/809820/contributions/3632617/attachments/1970786/3278138/MassDecorrelation_ML4Jets_H_Qu.pdf
        self.a.Define('FatJet_particleNetMD_HbbvsQCD','FatJet_particleNetMD_Xbb/(FatJet_particleNetMD_Xbb+FatJet_particleNetMD_QCD)')
        self.a.Define('FatJet_particleNetMD_WqqvsQCD','(FatJet_particleNetMD_Xqq+FatJet_particleNetMD_Xcc)/(FatJet_particleNetMD_Xqq+FatJet_particleNetMD_Xcc+FatJet_particleNetMD_QCD)')

        #Pick two highest pt back to back fat jets
        self.a.Define('DijetIdxs','PickDijets(FatJet_pt_corr,FatJet_eta,FatJet_phi,FatJet_msoftdrop_corr)')
        self.a.Cut('Dijets_exist','DijetIdxs[0] != -1 && DijetIdxs[1] != -1')

        self.nDijets = self.getNweighted()
        self.AddCutflowColumn(self.nDijets, "nDijets")

        self.a.SubCollection('Dijet','FatJet','DijetIdxs',useTake=True)

        # Apply an ID tag to MC jets from truth matching:
            # 0: Merged top jet
            # 1: Merged W jet
            # 2: Merged Higgs jet
            # 3: Unmerged jet (other)
        # Setup to work for signal and ttbar only
        if not self.a.isData:
            self.a.Define('Dijet_jetFlavor', 'JetFlavor_{}_vec(Dijet_eta, Dijet_phi, GenPart_eta, GenPart_phi, GenPart_pdgId, GenPart_genPartIdxMother)'.format('signal' if 'XHY' in self.setname else 'ttbar'))

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

        self.a.Cut('leptonQuality','LeptonType == 0? Electron_mvaFall17V2noIso_WP90[kinEleIdx] : Muon_mediumId[kinMuIdx]')

        self.nkinLep = self.getNweighted()
        self.AddCutflowColumn(self.nkinLep,'nkinLep')

        #For ease, merge some lepton columns that will be useful later (for lepton-type specific variables, use LeptonType to determine if electron or muon)
        self.a.Define('Lepton_pt','LeptonType == 1 ? Muon_pt[kinMuIdx] : Electron_pt[kinEleIdx]')
        self.a.Define('Lepton_eta','LeptonType == 1 ? Muon_eta[kinMuIdx] : Electron_eta[kinEleIdx]')
        self.a.Define('Lepton_phi','LeptonType == 1 ? Muon_phi[kinMuIdx] : Electron_phi[kinEleIdx]')
        self.a.Define('Lepton_mass','LeptonType == 1 ? Muon_mass[kinMuIdx] : Electron_mass[kinEleIdx]')
        self.a.Define('Lepton_vect','hardware::TLvector(Lepton_pt,Lepton_eta,Lepton_phi,Lepton_mass)')

        return self.a.GetActiveNode()

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
                SR: Hbb > 0.98
                    |deltaPhi| lepton,Higgs > 0.5
                ttCR_HF: Hbb < 0.98
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
        # SRorCR = "SR" or "CR" (signal or control - current region in which to create pass/fail subregions)
        # Return dictionary of 2 nodes corresponding to pass/fail regions within either signal or control region
        # Pass: Wqq > 0.8
        # Fail: Wqq < 0.8
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

    def get_standard_int(self,region,binsX,binsY,node):
        # For use with VV MC samples for which there appear bins w/ negative event count in MX vs MY histograms
        # Method: loosen standard cuts in order to increase systematics in histogram (should help eliminate negative bins), then rescale to maintain proper normalization
        # The idea behind this function is to perform all the "standard cuts" for a given region (pass/fail for SR, CR, or ttCR) and produce a histogram
        # Then, we count the total # events in the histogram so that when we apply comparatively looser cuts to the MX vs MY phase space we can scale by (# events passing loose selection)/(# events passing standard selection)
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
def JMEvariationStr(p, jet, variation):
    base_calibs = ['{}_JES_nom'.format(jet),'{}_JER_nom'.format(jet),'{}_JMS_nom'.format(jet),'{}_JMR_nom'.format(jet)]
    variationType = variation.split('_')[0]
    if jet == 'FatJet': #RVec quantities, set up the vector for a multi-hadamard product
        pt_calib_str = '{'
        mass_calib_str = '{'
        for c in base_calibs:
            if 'JM' in c and p != 'Top':	# WARNING - might need to change this if we treat W, H differently for mass and pt calibrations 
                mass_calib_str += '%s,'%(jet+'_'+variation if variationType in c else c)
            elif 'JE' in c:
                pt_calib_str += '%s,'%(jet+'_'+variation if variationType in c else c)
                mass_calib_str += '%s,'%(jet+'_'+variation if variationType in c else c)
        pt_calib_str += '}'
        mass_calib_str += '}'
    else: #Assuming input is a single jet and not a collection - create a string to multiply all the relevant calibrations
        pt_calib_str = ''
        mass_calib_str = ''
        for c in base_calibs:
            if 'JM' in c and p != 'Top':        # WARNING - might need to change this if we treat W, H differently for mass and pt calibrations 
                mass_calib_str += '*%s'%(jet+'_'+variation if variationType in c else c)
            elif 'JE' in c:
                pt_calib_str += '*%s'%(jet+'_'+variation if variationType in c else c)
                mass_calib_str += '*%s'%(jet+'_'+variation if variationType in c else c)

    return pt_calib_str, mass_calib_str

# helper function to make the ratio of two 2D histograms
def make_ratio_hist(histN,histD):
    ratio_hist = histN.Clone('ratio_hist') #New histogram - ratio of (tagged jets)/(all jets) per bin
    #Loop over all bins, manually calculate ratios + set the bin contents of new histogram (I don't think there is a better way to do this for TH2?)
    for x in range(1,ratio_hist.GetNbinsX()+1):
        for y in range(1,ratio_hist.GetNbinsY()+1):
            val_tagged = histN.GetBinContent(x,y)
            val_all = histD.GetBinContent(x,y)
            if val_all == 0:
                ratio = 0
            else:
                ratio = val_tagged/val_all
            ratio_hist.SetBinContent(x,y,ratio)
    return ratio_hist
