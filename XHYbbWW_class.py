import ROOT
from TIMBER.Analyzer import Correction, CutGroup, VarGroup, ModuleWorker, analyzer, HistGroup
from TIMBER.Tools.Common import CompileCpp, OpenJSON
from TIMBER.Tools.AutoPU import ApplyPU
from JMEvalsOnly import JMEvalsOnly
import TIMBER.Tools.AutoJME as AutoJME
import TIMBER.Tools.AutoJME_correctionlib as AutoJME_correctionlib
from collections import OrderedDict
from array import array

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
    def __init__(self, inputfile, ijob=1, njobs=1, use_xrd_global = False):
        #inputfile format is '(raw_nano or snapshots)/SETNAME_YEAR(maybe _+something else).txt'
        infiles = SplitUp(inputfile,njobs)[ijob-1]
        print(infiles)
        if use_xrd_global:
            for i in range(len(infiles)):
                infiles[i] = infiles[i].replace('cmsxrootd.fnal.gov','cms-xrd-global.cern.ch')
            print('SEARCHING FOR FILES WITH GLOBAL REDIRECTOR')
            print(infiles)

        self.setname = inputfile.split('/')[1].split('_')[0]
        print(self.setname)
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
 
    def getNweighted(self, MCWeight='genWeight'):
        if not self.a.isData:
            return self.a.DataFrame.Sum(f'{MCWeight}').GetValue()
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

    def ApplyStandardCorrections(self, snapshot=False, do_ak4JEC = True, JMC = 'all'):
        # str JMC: specify which mass collection(s) to which JEC's should be applied (defaults to 'all', additional options are 'softdrop' or 'regressed')
        # first apply corrections for snapshot phase
        if snapshot:
            if self.a.isData:
                # DATA - only filter valid lumi blocks and drop HEM-affected events
                # NOTE: LumiFilter requires the year as an integer 
                lumiFilter = ModuleWorker('LumiFilter','TIMBER/Framework/include/LumiFilter.h',[int(self.year) if 'APV' not in self.year else 16]) # defaults to perform "eval" method 
                self.a.Cut('lumiFilter',lumiFilter.GetCall(evalArgs={"lumi":"luminosityBlock"})) # replace lumi with luminosityBlock
                '''
                if self.year == '18':
                    print('pre HEM: {} Events'.format(self.a.DataFrame.Count().GetValue()))
                    HEM_worker = ModuleWorker('HEM_drop','TIMBER/Framework/include/HEM_drop.h',[self.setname if 'Data' not in self.setname else "data"+self.setname[-1]])
                    self.a.Cut('HEM','%s[0] > 0'%HEM_worker.GetCall(evalArgs={"FatJet_eta":"FatJet_eta","FatJet_phi":"FatJet_phi"}))
                    print('post HEM: {} Events'.format(self.a.DataFrame.Count().GetValue()))
                '''
            else:
                # MC - apply corrections
                # Parton shower weights 
                #	- https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopSystematics#Parton_shower_uncertainties
                #	- "Default" variation: https://twiki.cern.ch/twiki/bin/view/CMS/HowToPDF#Which_set_of_weights_to_use
                #	- https://github.com/mroguljic/Hgamma/blob/409622121e8ab28bc1072c6d8981162baf46aebc/templateMaker.py#L210
                self.a.Define("ISR__up","PSWeight[0]")
                self.a.Define("ISR__down","PSWeight[2]")
                self.a.Define("FSR__up","PSWeight[1]")
                self.a.Define("FSR__down","PSWeight[3]")
                genWCorr = Correction('genW','TIMBER/Framework/TopPhi_modules/BranchCorrection.cc',corrtype='corr',mainFunc='evalCorrection') # workaround so we can have multiple BCs
                self.a.AddCorrection(genWCorr, evalArgs={'val':'genWeight'})
                ISRcorr = genWCorr.Clone("ISRunc",newMainFunc="evalUncert",newType="uncert")
                FSRcorr = genWCorr.Clone("FSRunc",newMainFunc="evalUncert",newType="uncert")
                self.a.AddCorrection(ISRcorr, evalArgs={'valUp':'ISR__up','valDown':'ISR__down'})
                self.a.AddCorrection(FSRcorr, evalArgs={'valUp':'FSR__up','valDown':'FSR__down'})
                # QCD factorization and renormalization corrections (only to non-signal MC, but generate signal w this variation just in case..)
                # For some reason, diboson processes don't have the LHEScaleWeight branch, so don't apply to those either.
                # Reminder: nominal values of these corrections are 1 (effectively a "corr" type correction)
                if (('WW' not in self.setname) and ('WZ' not in self.setname) and ('ZZ' not in self.setname) and ('XHY' not in self.setname)):
                    # First instatiate a correction module for the factorization correction
                    facCorr = Correction('QCDscale_factorization','modules/LHEScaleWeights.cc',corrtype='weight',mainFunc='evalFactorization')
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

                # Pileup correction
                self.a = ApplyPU(self.a, 'XHYbbWWpileup.root', '20{}'.format(self.year), ULflag=True, histname='{}_{}'.format(self.setname,self.year))

                # PDF weight correction - https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopSystematics#PDF
                if self.a.lhaid != -1: #some backgrounds (diboson) don't have LHA ID in MC file
                    PdfWeight = Correction('Pdfweight','TIMBER/Framework/include/PDFweight_uncert.h',[self.a.lhaid],corrtype='uncert')
                    self.a.AddCorrection(PdfWeight)

                # Level-1 prefire corrections
                if self.year == '16' or self.year == '17' or 'APV' in self.year:
                    L1PreFiringWeight = genWCorr.Clone('L1PreFiringWeight',newMainFunc='evalWeight',newType='weight')
                    self.a.AddCorrection(L1PreFiringWeight, evalArgs={'val':'L1PreFiringWeight_Nom','valUp':'L1PreFiringWeight_Up','valDown':'L1PreFiringWeight_Dn'})

                # HEM drop to 2018 MC
                elif self.year == '18':
                    HEM = Correction('HEM_drop','TIMBER/Framework/include/HEM_drop.h',[self.setname],corrtype='corr')
                    self.a.AddCorrection(HEM)

            if 'Data' not in self.setname:
                JMEname = self.setname
            elif 'DataB' in self.setname: #account for extra v1/v2 at the end of the setname
                JMEname = 'DataB'
            else:
                JMEname = self.setname[-5:]
            self.a = AutoJME.AutoJME(self.a, 'FatJet', '20{}'.format(self.year), JMEname, mass = JMC)
            #self.a = AutoJME_correctionlib.AutoJME(self.a, 'FatJet', '20{}_UL'.format(self.year), dataEra=JMEname, calibrate=True) # just using this to get the L1 corrections

            if do_ak4JEC:
                self.a = AutoJME.AutoJME(self.a, 'Jet', '20{}'.format(self.year), JMEname)
                #self.a = AutoJME_correctionlib.AutoJME(self.a, 'Jet', '20{}_UL'.format(self.year), dataEra=JMEname, calibrate=True)

            #print('CORRECTIONS:')
            #for correction in self.a.GetCorrectionNames():
            #    print(' - '+correction)
            #self.a.MakeWeightCols(correctionNames=list(self.a.GetCorrectionNames()),extraNominal='' if self.a.isData else str(self.GetXsecScale()))

        # now for selection
        else:
            if self.a.isData and self.year == '18':
                # I messed up the HEM drop cut for the snapshots, so I will apply it again here
                print('pre HEM: {}'.format(self.a.DataFrame.Count().GetValue()))
                HEM_worker = ModuleWorker('HEM_drop','TIMBER/Framework/include/HEM_drop.h',[self.setname if 'Data' not in self.setname else "data"+self.setname[-1]])
                self.a.Cut('HEM','%s[0] > 0'%HEM_worker.GetCall(evalArgs={"FatJet_eta":"FatJet_eta","FatJet_phi":"FatJet_phi"}))
                print('post HEM: {}'.format(self.a.DataFrame.Count().GetValue()))
          
            elif not self.a.isData:
                self.a.AddCorrection(Correction('genW',corrtype='corr'))
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

    
    def ApplyTopPtReweight(self,jet0,jet1,scale = 0.5, isJet1AK4 = False):
        #isJet1AK4 asks if the second jet is an AK4 jet, in which case definitions of mass/momentum are different
        r0 = 0.8
        r1 = 0.8 #0.4 if isJet1AK4 else 0.8
        if 'ttbar' in self.setname:
            if isJet1AK4:
                jet1_pt = 'pt'
                jet1_mass = 'mass'
            else:
                jet1_pt = 'pt'
                jet1_mass = 'msoftdrop'
            self.a.Define('GenParticle_vect','hardware::TLvector(GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)')
            self.a.AddCorrection(
                Correction('TptReweight','TIMBER/Framework/include/TopPt_weight.h',corrtype='weight'),
                    evalArgs={
                        'jet0':'hardware::TLvector({0}_pt,{0}_eta,{0}_phi,{0}_msoftdrop)'.format(jet0),
                        'jet1':'hardware::TLvector({0}_{1},{0}_eta,{0}_phi,{0}_{2})'.format(jet1,jet1_pt,jet1_mass),
                        'GenPart_vect':'GenParticle_vect',
                        'r0':r0,
                        'r1':r1,
                        'scale':scale
                    }
             )  
   
    def ApplyJMECorrections(self, variation, jets = ['FatJet','Jet']):
        # for trigger effs
        #self.a.Define('Trijet_vect_trig','hardware::TLvector(Trijet_pt, Trijet_eta, Trijet_phi, Trijet_msoftdrop)')
        #self.a.Define('mhww_trig','hardware::InvariantMass(Trijet_vect_trig)')

        if 'Jet' in jets:
            # Get uncorrected AK4 jet pt, so that we can apply the JME variations
            if '2' in self.setname and '12' not in self.setname:
                #self.a.Define('Jet_JES_nom_inv','invert_vector(Jet_JES_nom)')
                #self.a.Define('Jet_pt_raw','hardware::HadamardProduct(Jet_pt, Jet_JES_nom_inv)')
                self.a.Define('Jet_rawFactor','zeros_like(Jet_pt)') # For later
            
            self.a.Define('Jet_pt_raw','Jet_pt - hardware::HadamardProduct(Jet_pt,Jet_rawFactor)') # pt_raw = pt * (1 - rawFactor)

        # JME variations - we only do this for MC
        for jet in jets:
            if not self.a.isData:
                # since H, W close enough in mass, we can treat them the same. 
                # Higgs, W will have same pt and mass calibrations
                # WARNING --------------------------------------------------------------------------------------------------------------
                # IS THIS ACTUALLY TRUE?
                # IS IT, AMITAV? IS IT???
                # ----------------------------------------------------------------------------------------------------------------------
                pt_calibs, sd_mass_calibs, reg_mass_calibs = JMEvariationStr('Higgs',jet,variation)
                if jet == 'FatJet': #dealing with RVec quantities, need to do hadamard product
                    self.a.Define(f'{jet}_pt_corr',f'hardware::MultiHadamardProduct({jet}_pt,{pt_calibs})')
                    self.a.Define(f'{jet}_msoftdrop_corr',f'hardware::MultiHadamardProduct({jet}_msoftdrop,{sd_mass_calibs})')
                    self.a.Define(f'{jet}_particleNet_mass_corr',f'hardware::MultiHadamardProduct({jet}_particleNet_mass,{reg_mass_calibs})')
                elif jet == 'Jet':
                    self.a.Define(f'{jet}_pt_corr',f'hardware::MultiHadamardProduct({jet}_pt_raw,{pt_calibs})')
                else: #Assume dealing with a single jet and not a collection
                    self.a.Define(f'{jet}_pt_corr'.format,'%s_pt%s'%(jet,pt_calibs))
                    self.a.Define(f'{jet}_msoftdrop_corr','%s_msoftdrop%s'%(jet,sd_mass_calibs))
                    self.a.Define(f'{jet}_particleNet_mass_corr','%s_particleNet_mass%s'%(jet,reg_mass_calibs))
            else:
                if jet == 'FatJet':
                    self.a.Define(f'{jet}_pt_corr','hardware::MultiHadamardProduct(%s_pt,{%s_JES_nom})'%(jet,jet))
                    self.a.Define(f'{jet}_msoftdrop_corr','hardware::MultiHadamardProduct(%s_msoftdrop,{%s_JES_nom})'%(jet,jet))
                    self.a.Define(f'{jet}_particleNet_mass_corr','hardware::MultiHadamardProduct(%s_particleNet_mass,{%s_JES_nom})'%(jet,jet))
                elif jet == 'Jet':
                    self.a.Define(f'{jet}_pt_corr','hardware::MultiHadamardProduct(%s_pt_raw,{%s_JES_nom})'%(jet,jet))
                else:
                    self.a.Define(f'{jet}_pt_corr','%s_pt * %s_JES_nom'%(jet,jet))
                    self.a.Define(f'{jet}_msoftdrop_corr','%s_msoftdrop * %s_JES_nom'%(jet,jet))
                    self.a.Define(f'{jet}_particleNet_mass_corr','%s_particleNet_mass * %s_JES_nom'%(jet,jet))
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

    def ApplyPNetReweight(self,Hbb_wp=0.98, Wqq_wp=0.80, doWtag=True): # wp chosen to delineate Hbb pass/fail
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

            if doWtag: #might not use W tag in pass/fail definition, in which case don't apply W scale factors 
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
            print(Hbb_wp)
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

    def ApplyLeptonCorrections(self,eleId='medium',muId='medium'):
        if not self.a.isData:
            eleId_tags = {'medium':'wp90noiso','tight':'wp80noiso'}
            muId_tags = {'medium':'NUM_MediumID_DEN_GlobalMuonProbes','tight':'NUM_TightID_DEN_GlobalMuonProbes'}
        
            ele_tag = eleId_tags[eleId]
            mu_tag = muId_tags[muId]

            #Corrections for electron mva ID (wp80 no isolation)
            cols = ['Lepton_eta','Lepton_pt','LeptonType']
            ele_filepath = 'corrections/electron_{}.json'.format(self.year) #path to json file - for constructor
            ElectronIDWeight = Correction('ElectronIDWeight','TIMBER/Framework/include/ElectronID_weight.h',
                constructor=[str(self.year), ele_tag, ele_filepath], #proper pt, eta bounds included as default argument in c++ code
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
            muoID_cols = ['Lepton_abseta','Lepton_pt','LeptonType']
            MuonIDWeight = Correction('MuonIDWeight','TIMBER/Framework/include/MuonID_weight.h',
                constructor=[mu_tag, muoID_filepath], 
                mainFunc='eval', 
                corrtype='corr', 
                columnList=muoID_cols 
            )
            MuonIDWeight_uncert_stat = MuonIDWeight.Clone('MuonIDWeight_uncert_stat',newType='uncert')
            MuonIDWeight_uncert_syst = MuonIDWeight.Clone('MuonIDWeight_uncert_syst',newType='uncert')
            self.a.AddCorrection(MuonIDWeight, evalArgs={'Lepton_eta':'Lepton_abseta','Lepton_pt':'Lepton_pt','LeptonType':'LeptonType','uncert':'"nominal"'})
            self.a.AddCorrection(MuonIDWeight_uncert_stat, evalArgs={'Lepton_eta':'Lepton_abseta','Lepton_pt':'Lepton_pt','LeptonType':'LeptonType','uncert':'"stat"'})
            self.a.AddCorrection(MuonIDWeight_uncert_syst, evalArgs={'Lepton_eta':'Lepton_abseta','Lepton_pt':'Lepton_pt','LeptonType':'LeptonType','uncert':'"syst"'})

            #Muon reconstruction scale factors
            self.a.Define('Lepton_p','Lepton_pt*cosh(Lepton_eta)') #Muon reco scale factors indexed by p (not pt)
            muoRECO_filepath = 'corrections/ScaleFactors_Muon_highPt_RECO_20{}_schemaV2.json'.format(self.year)
            muoRECO_cols = ['Lepton_abseta','Lepton_p','LeptonType']
            MuonRecoWeight = Correction('MuonRecoWeight','TIMBER/Framework/include/MuonID_weight.h',
                constructor=["NUM_GlobalMuons_DEN_TrackerMuonProbes",muoRECO_filepath], 
                mainFunc='eval', 
                corrtype='corr', 
                columnList=muoRECO_cols
            )
            MuonRecoWeight_uncert_stat = MuonRecoWeight.Clone('MuonRecoWeight_uncert_stat',newType='uncert')
            MuonRecoWeight_uncert_syst = MuonRecoWeight.Clone('MuonRecoWeight_uncert_syst',newType='uncert')
            self.a.AddCorrection(MuonRecoWeight, evalArgs={'Lepton_eta':'Lepton_abseta','Lepton_pt':'Lepton_p','LeptonType':'LeptonType','uncert':'"nominal"'})
            self.a.AddCorrection(MuonRecoWeight_uncert_stat, evalArgs={'Lepton_eta':'Lepton_abseta','Lepton_pt':'Lepton_p','LeptonType':'LeptonType','uncert':'"stat"'})
            self.a.AddCorrection(MuonRecoWeight_uncert_syst, evalArgs={'Lepton_eta':'Lepton_abseta','Lepton_pt':'Lepton_p','LeptonType':'LeptonType','uncert':'"syst"'})

    def GetXsecScale(self):
        lumi = self.config['lumi{}'.format(self.year)]
        xsec = self.config['XSECS'][self.setname]
        if self.a.genEventSumw == 0:
            raise ValueError('%s %s: genEventSumw is 0'%(self.setname, self.year))
        print('Normalizing by lumi*xsec/genEventSumw:\n\t{} * {} / {} = {}'.format(lumi,xsec,self.a.genEventSumw,lumi*xsec/self.a.genEventSumw))
        return lumi*xsec/self.a.genEventSumw

    # for trigger efficiencies
    def ApplyTrigs(self):

        from modules.TriggerSF_weight import MuonSF, ElectronSF

        all_trigs = self.trigs[self.year if '16' not in self.year else '16']
        h_trigs = all_trigs['JETHT']
        e_trigs = all_trigs['ELECTRON'] if self.year != '18' else all_trigs['EGAMMA']

        m_trigs_lowPt = all_trigs['MUON']['lowPt']
        m_trigs_lowPt_str = self.a.GetTriggerString(m_trigs_lowPt)

        m_trigs_highPt = all_trigs['MUON']['highPt']
        m_trigs_highPt_str = self.a.GetTriggerString(m_trigs_highPt)

        m_trigs = m_trigs_lowPt + m_trigs_highPt

        if self.year != '18':
            p_trigs = all_trigs['PHOTON']
        else:
            p_trigs = [] # Quick fix for applying trigger bits in MC

        if self.a.isData:
            if 'Muon' in self.setname:
                # First veto electron events
                self.a.Cut('ele_eventVeto','LeptonType == 1')
                print('Muon trigger list: {}'.format(m_trigs))
                # Now apply the muon triggers
                self.a.Cut('m_trigger',f'Lepton_pt < 55 ? {m_trigs_lowPt_str} : {m_trigs_highPt_str}')

            elif 'Electron' in self.setname or 'EGamma' in self.setname:
                # Veto muon events
                self.a.Cut('mu_eventVeto','LeptonType == 0')
                print('Electron/EGamma trigger list: {}'.format(e_trigs))
                self.a.Cut('e_trigger',self.a.GetTriggerString(e_trigs))

            elif 'Photon' in self.setname: #for year != 2018
                # Veto muon events
                self.a.Cut('ele_eventVeto','LeptonType == 0')
                print('Photon trigger list: {}'.format(p_trigs))
                self.a.Cut('p_eleOrthogonal','!'+self.a.GetTriggerString(e_trigs)) # must veto electron triggers to ensure orthogonality
                self.a.Cut('p_trigger',self.a.GetTriggerString(p_trigs))

        else: # Apply simulated trigger bits and corrections
            egm_string = self.a.GetTriggerString(e_trigs+p_trigs)
            muon_string_low = '{'+(',').join(m_trigs_lowPt)+'}'
            muon_string_high = '{'+(',').join(m_trigs_highPt)+'}'

            self.a.Define('muon_trigger_status',f'check_muon_trigger_status(LeptonType, Lepton_pt, {muon_string_low}, {muon_string_high})')
            self.a.Cut('m_trigger','muon_trigger_status')
            self.a.Cut('egm_trigger',f'LeptonType == 0 ? {egm_string} : true')
            
            #Now apply corrections
            MuonSF(self)
            ElectronSF(self)

        self.nTrigs = self.getNweighted()
        self.AddCutflowColumn(self.nTrigs, "nTrigs")
        return self.a.GetActiveNode()
    
    '''
    # for trigger efficiencies
    def ApplyTrigs(self, corr=None, useMCtrigs = False, channel = None):
        all_trigs = self.trigs[self.year if '16' not in self.year else '16']
        h_trigs = all_trigs['JETHT']
        e_trigs = all_trigs['ELECTRON'] if self.year != '18' else all_trigs['EGAMMA']

        m_trigs_lowPt = all_trigs['MUON']['lowPt']
        m_trigs_highPt = all_trigs['MUON']['highPt']
        m_trigs = m_trigs_lowPt + m_trigs_highPt

        if self.year != '18':
            p_trigs = all_trigs['PHOTON']
        else:
            p_trigs = [] # Quick fix for applying trigger bits in MC

        if self.a.isData:
            #The basic logic is that for a given primary dataset, we need to accept events which pass the corresponding triggers for that dataset, but reject events which also pass the triggers of other datasets (in order to avoid double counting). We start by accepting all events which pass SingleMuon triggers and then go down the ladder. Note that for 2018, the electron/photon triggers were combined into a single PD (EGamma), but are separate (Single Electron + Single Photon) for the other years.
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

        elif (self.a.isData == False and useMCtrigs == False):
            #self.a.AddCorrection(corr, evalArgs={'xval':'Lepton_abseta','yval':'Lepton_pt','zval':'Lepton_miniPFRelIso_all'})
            self.a.AddCorrection(corr, evalArgs={'xval':'Lepton_abseta','yval':'Lepton_pt','zval':'LeptonType'})
        
        else: # Apply simulated trigger bits
            trig_list = e_trigs+m_trigs+p_trigs
            if channel == None:
                print('Applying all triggers to MC: {}'.format(trig_list))
                self.a.Cut('all_trigger',self.a.GetTriggerString(trig_list))
            elif channel == 'Electron':
                print('Applying EGamma triggers to MC: {}'.format(e_trigs+p_trigs))
                self.a.Cut('egamma_trigger',self.a.GetTriggerString(e_trigs+p_trigs))
            elif channel == 'Muon':
                print('Applying Muon triggers to MC: {}'.format(m_trigs))
                self.a.Cut('muon_trigger',self.a.GetTriggerString(m_trigs))

        self.nTrigs = self.getNweighted()
        self.AddCutflowColumn(self.nTrigs, "nTrigs")
        return self.a.GetActiveNode()
    '''

    def Snapshot(self, node=None):
        startNode = self.a.GetActiveNode()
        if node == None:
            node = self.a.GetActiveNode()
        columns = [       
            'nJet','Jet_btagDeepFlavB','Jet_pt','Jet_eta','Jet_phi','Jet_mass','Jet_btagDeepB','Jet_jetId','Jet_hadronFlavour','Jet_partonFlavour','Jet_JES_nom',
            'Jet_rawFactor','Jet_puId','Jet_puIdDisc','nFatJet','FatJet_eta','FatJet_msoftdrop','FatJet_pt','FatJet_phi','FatJet_JES_nom','FatJet_particleNetMD*', 
            'FatJet_rawFactor','FatJet_subJetIdx1','FatJet_subJetIdx2','FatJet_hadronFlavour','FatJet_nBHadrons','FatJet_nCHadrons','FatJet_btagDeepB','nSubJet','SubJet_*'
            'FatJet_btagCMVA','FatJet_btagCSVV2','FatJet_jetId','FatJet_particleNet_mass','nCorrT1METJet','CorrT1METJet_*','nMuon','Muon_*','nElectron','Electron_*',
            'fixedGridRhoFastjetAll','RawMET_*','ChsMET_*','MET_*','PV_*','nOtherPV','OtherPV_z','fixedGridRhoFastjetAll',
            'genWeight','event','eventWeight','luminosityBlock','run','NSTART','NFLAGS','LEPPRE','JETPRE','METPT'
        ]
        #'Jet_puId','Jet_puIdDisc',
        columns.extend(['nSoftActivityJet','nGenJet','nboostedTau','nGenVisTau','nIsoTrack','nTau']) # Columns I'm forced to add because why
        columns.extend(['HLT_AK8*','HLT_PFHT*','HLT_PFJet*','HLT_Iso*','HLT_Mu*','HLT_TkMu*','HLT_OldMu*','HLT_Ele*','HLT_Photon*','HLT_MET*','HLT_PFMET*'])
        
        # Columns to save potentially for MET re-correction 
        columns.extend([
            'Jet_area',
            'Jet_chEmEF',
            'Jet_chFPV0EF',
            'Jet_chFPV1EF',
            'Jet_chFPV2EF',
            'Jet_chFPV3EF',
            'Jet_chHEF',
            'Jet_cleanmask',
            'Jet_electronIdx1',
            'Jet_electronIdx2',
            'Jet_genJetIdx',
            'Jet_muEF',
            'Jet_muonIdx1',
            'Jet_muonIdx2',
            'Jet_muonSubtrFactor',
            'Jet_nConstituents',
            'Jet_nElectrons',
            'Jet_nMuons',
            'Jet_neEmEF',
            'Jet_neHEF'
        ])
         
        if not self.a.isData:
            
            columns.extend(['nGenPart','GenPart_*','GenMET_*','genWeight']) 
            
            columns.extend(['Jet_JES_up','Jet_JES_down','Jet_JER_nom','Jet_JER_up','Jet_JER_down'])
            columns.extend(['FatJet_JES_up','FatJet_JES_down',
                            'FatJet_JER_nom','FatJet_JER_up','FatJet_JER_down',
                            'FatJet_JMS_softdrop_nom','FatJet_JMS_softdrop_up','FatJet_JMS_softdrop_down',
                            'FatJet_JMR_softdrop_nom','FatJet_JMR_softdrop_up','FatJet_JMR_softdrop_down',
                            'FatJet_JMS_regressed_nom','FatJet_JMS_regressed_up','FatJet_JMS_regressed_down',
                            'FatJet_JMR_regressed_nom','FatJet_JMR_regressed_up','FatJet_JMR_regressed_down'])
            columns.extend(['genW__nom','Pileup__nom','Pileup__up','Pileup__down','Pdfweight__up','Pdfweight__down','ISRunc__up','ISRunc__down','FSRunc__up','FSRunc__down'])
        
            if self.year == '16' or self.year == '17' or 'APV' in self.year:
                columns.extend(['L1PreFiringWeight_Nom', 'L1PreFiringWeight_Up', 'L1PreFiringWeight_Dn'])       # these are the default columns in NanoAODv9
                columns.extend(['L1PreFiringWeight__nom','L1PreFiringWeight__up','L1PreFiringWeight__down'])    # these are the weight columns created by the BranchCorrection module
            elif self.year == '18':
                columns.append('HEM_drop__nom')
            if (('WW' not in self.setname) and ('WZ' not in self.setname) and ('ZZ' not in self.setname) and ('XHY' not in self.setname)):
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
        self.a.Define('DijetIdxs','PickDijets(FatJet_pt_corr,FatJet_eta,FatJet_particleNet_mass_corr,FatJet_jetId)')
        self.a.Cut('Dijets_exist','DijetIdxs[0] != -1 && DijetIdxs[1] != -1')

        self.nDijets = self.getNweighted()
        self.AddCutflowColumn(self.nDijets, "nDijets")

        self.a.SubCollection('Dijet','FatJet','DijetIdxs',useTake=True)
        print(self.a.DataFrame.Sum('nFatJet').GetValue())
        print(self.a.DataFrame.Sum('nDijet').GetValue())

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

        self.a.Define('kinEleIdx','kinElectron(Electron_pt,Electron_eta,Electron_miniPFRelIso_all,Electron_dxy,Electron_dz,Electron_sip3d)')
        self.a.Define('kinMuIdx','kinMuon(Muon_pt,Muon_eta,Muon_miniPFRelIso_all,Muon_dxy,Muon_dz,Muon_sip3d)')  
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
        self.a.Define('Lepton_miniPFRelIso_all','LeptonType == 1 ? Muon_miniPFRelIso_all[kinMuIdx] : Electron_miniPFRelIso_all[kinEleIdx]')
        self.a.Define('Lepton_vect','hardware::TLvector(Lepton_pt,Lepton_eta,Lepton_phi,Lepton_mass)')
        self.a.Define('Lepton_abseta','abs(Lepton_eta)') #used for trigger efficiencies/lepton corrections

        return self.a.GetActiveNode()

    def ApplyHiggsMass(self,sidebands = False):
        if not sidebands:
            # perform Higgs mass window cut, save cutflow info
            self.a.Cut('mH_window_cut','Higgs_particleNet_mass_corr >= {0} && Higgs_particleNet_mass_corr <= {1}'.format(*self.cuts['JET']['mh']))
        else:
            self.a.Cut('mH_sideband_cut','Higgs_particleNet_mass_corr < {0} || Higgs_particleNet_mass_corr > {1}'.format(*self.cuts['JET']['mh']))

        return self.a.GetActiveNode()

    def ApplyWMass(self): 
        # W mass window cut, cutflow
        self.a.Cut('mW_{}_cut'.format('window'),'Wqq_particleNet_mass_corr >= {0} && Wqq_particleNet_mass_corr <= {1}'.format(self.cuts['JET']['mw'][0], self.cuts['JET']['mw'][1]))
        self.nWqq = self.getNweighted()
        self.AddCutflowColumn(self.nWqq, 'nWqqMassCut')
	
        return self.a.GetActiveNode()

    def ApplySRorCR(self, SRorCR, tagger): #This is getting frequently updated as I play w/ different definitions
        '''
        SRorCR [str] = 'SR' or 'CR' or 'ttCR', used to generate cutflow information
        tagger [str] = name of tagger to be used

        Work with two control regions to measure ttbar - one defined by Hbb tag fail and other defined by inverting delta phi cut between lepton/Higgs
        Selection criteria:
                SR: Hbb > 0.98, Wqq > 0.8, b-tag veto
                ttCR: Hbb > 0.8, Wqq > 0.8, invert b-tag veto
                b-tag veto involves removing events with an AK4 jet satisfying the following:
                    -pT > 25 
                    -|eta| < 2.4
                    -DeltaR(jet,Higgs) > 1.2
                    -DeepJet score > tight wp
        '''
        assert(SRorCR=='SR' or SRorCR=='CR' or SRorCR=='ttCR')
        # tagger values for each 
        HMP = self.cuts['JET']['particleNetMD_HbbvsQCD'][0]
        HHP = self.cuts['JET']['particleNetMD_HbbvsQCD'][1]

        if SRorCR == 'SR':
            # Signal region
            self.a.Cut('Higgs_{}_cut_{}'.format(tagger,SRorCR), 'Higgs_{0}_HbbvsQCD >= {1}'.format(tagger,HHP))
            self.nHtag_SR = self.getNweighted()
            self.AddCutflowColumn(self.nHtag_SR, 'nHtag_SR')

            #b-tag AK4 veto
            self.a.Cut('ak4_veto_SR','!bjet_exists')
            self.nAK4VetoSR = self.getNweighted()
            self.AddCutflowColumn(self.nAK4VetoSR, 'nAK4VetoSR')

        elif SRorCR == 'ttCR':
            #ttbar-enriched CR
            self.a.Cut('Higgs_{}_cut_{}'.format(tagger,SRorCR), 'Higgs_{0}_HbbvsQCD >= {1}'.format(tagger,HMP))
            self.nHtag_ttCR = self.getNweighted()
            self.AddCutflowColumn(self.nHtag_ttCR, 'nHtag_ttCR')

            #Invert b-tag veto
            self.a.Cut('ak4_antiveto_ttCR','bjet_exists')
            self.nAK4AntiVetottCR = self.getNweighted()
            self.AddCutflowColumn(self.nAK4AntiVetottCR, 'nAK4AntiVetottCR')

        return self.a.GetActiveNode()

    def ApplyPassFail(self, SRorCR, tagger):
        # Generates isolation pass/fail subregions inside the SR/ttCRT
        # Creates a histogram of the isolation score in the Higgs mass sidebands
        assert(SRorCR=='SR' or SRorCR=='CR' or SRorCR=='ttCR')

        pre_massCut = self.a.GetActiveNode()
        FP = OrderedDict()
        W = self.cuts['JET']['particleNetMD_WqqvsQCD'][self.year]

        if SRorCR == 'SR':
            checkpoint = self.ApplyHiggsMass(sidebands=False) #First make all regions inside Higgs mass window
            self.a.SetActiveNode(checkpoint)
            self.nHiggs_SR = self.getNweighted()
            self.AddCutflowColumn(self.nHiggs_SR, 'nHiggs_SR')
 
            # W tag fail + cutflow info
            FP['fail_'+SRorCR] = self.a.Cut('WTag_fail_SR',f'Wqq_{tagger}_WqqvsQCD < {W}')
            self.nF_SR = self.getNweighted()
            self.AddCutflowColumn(self.nF_SR, 'nF_SR')            

            # W tag pass + cutflow
            self.a.SetActiveNode(checkpoint)
            FP['pass_'+SRorCR] = self.a.Cut('WTag_pass_SR',f'Wqq_{tagger}_WqqvsQCD >= {W}')
            self.nP_SR = self.getNweighted()
            self.AddCutflowColumn(self.nP_SR, 'nP_SR')

        elif SRorCR == 'CR':
            print('There is no life in the control region. Only death.')

        elif SRorCR == 'ttCR':
            checkpoint = self.ApplyHiggsMass(sidebands=False) #First make all regions inside Higgs mass window
            self.a.SetActiveNode(checkpoint)
            self.nHiggs_ttCR = self.getNweighted()
            self.AddCutflowColumn(self.nHiggs_ttCR, 'nHiggs_ttCR')

            # W tag fail + cutflow info
            FP['fail_'+SRorCR] = self.a.Cut('WTag_fail_ttCR',f'Wqq_{tagger}_WqqvsQCD < {W}')
            self.nF_ttCR = self.getNweighted()
            self.AddCutflowColumn(self.nF_ttCR, 'nF_ttCR') 

            # W tag pass + cutflow
            self.a.SetActiveNode(checkpoint)
            FP['pass_'+SRorCR] = self.a.Cut('WTag_pass_ttCR',f'Wqq_{tagger}_WqqvsQCD >= {W}')
            self.nP_ttCR = self.getNweighted()
            self.AddCutflowColumn(self.nP_ttCR, 'nP_ttCR')       

        #Now, return to pre-isolation cut stage and go into Higgs sideband in order to plot isolation score
        self.a.SetActiveNode(pre_massCut)
        tagGroup = HistGroup('tagGroup')
        
        sideband = self.ApplyHiggsMass(sidebands=True)
        self.a.SetActiveNode(sideband)
        WtagHist = self.a.GetActiveNode().DataFrame.Histo1D(('Wtag_{}'.format(SRorCR),'Wtag_{}'.format(SRorCR),50,0,1),'Wqq_particleNetMD_WqqvsQCD','weight__nominal')
        tagGroup.Add('Wtag_{}'.format(SRorCR),WtagHist)

        self.a.SetActiveNode(pre_massCut)
        
        return FP, tagGroup

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
    base_calibs = [f'{jet}_JES_nom',f'{jet}_JER_nom',f'{jet}_JMS_softdrop_nom',f'{jet}_JMR_softdrop_nom',f'{jet}_JMS_regressed_nom',f'{jet}_JMR_regressed_nom']
    variationType = variation.split('_')[0]
    if variationType == 'None' or variationType == 'UE':
        ud = 'nom'
    else:
        ud = variation.split('_')[1] #up or down
    if jet == 'FatJet' or jet == 'Jet': #RVec quantities, set up the vector for a multi-hadamard product
        pt_calib_str = '{'
        sd_mass_calib_str = '{'
        reg_mass_calib_str = '{'
        for c in base_calibs:
            if 'JM' in c and p != 'Top':	# WARNING - might need to change this if we treat W, H differently for mass and pt calibrations 
                if 'softdrop' in c:
                    sd_mass_calib_str += '%s,'%(jet+'_'+variationType+'_softdrop_'+ud if variationType in c else c)
                elif 'regressed' in c:
                    reg_mass_calib_str += '%s,'%(jet+'_'+variationType+'_regressed_'+ud if variationType in c else c)
            elif 'JE' in c:
                pt_calib_str += '%s,'%(jet+'_'+variation if variationType in c else c)
                sd_mass_calib_str += '%s,'%(jet+'_'+variation if variationType in c else c)
                reg_mass_calib_str += '%s,'%(jet+'_'+variation if variationType in c else c)
        pt_calib_str += '}'
        sd_mass_calib_str += '}'
        reg_mass_calib_str += '}'
    else: #Assuming input is a single jet and not a collection - create a string to multiply all the relevant calibrations
        pt_calib_str = ''
        sd_mass_calib_str = ''
        reg_mass_calib_str = ''
        for c in base_calibs:
            if 'JM' in c and p != 'Top':        # WARNING - might need to change this if we treat W, H differently for mass and pt calibrations 
                if 'softdrop' in c:
                    sd_mass_calib_str += '*%s'%(jet+'_'+variationType+'_softdrop_'+ud if variationType in c else c)
                elif 'regressed' in c:
                    reg_mass_calib_str += '*%s'%(jet+'_'+variationType+'_regressed_'+ud if variationType in c else c)
            elif 'JE' in c:
                pt_calib_str += '*%s'%(jet+'_'+variation if variationType in c else c)
                mass_calib_str += '*%s'%(jet+'_'+variation if variationType in c else c)

    return pt_calib_str, sd_mass_calib_str, reg_mass_calib_str

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
