from argparse import Namespace
from array import array
from glob import glob
from TIMBER.Tools.Common import DictStructureCopy, CompileCpp, ExecuteCmd, OpenJSON, StitchQCD
from TIMBER.Tools.Plot import CompareShapes, EasyPlots
from TIMBER.Analyzer import HistGroup
from collections import OrderedDict
from make_QCD import Rebin2D
import ROOT
import sys
sys.path.append('./modules')
from TwoDHistStyle import setTDRStyle

def outline(hist, canvas):
    #Draw outlines around bins of an input 2D histogram
    for i in range(1,hist.GetNbinsX()+1):
        for j in range(1,hist.GetNbinsY()+1):
            #Get the edges of the bin
            x1 = hist.GetXaxis().GetBinLowEdge(i)
            x2 = hist.GetXaxis().GetBinUpEdge(i)
            y1 = hist.GetYaxis().GetBinLowEdge(j)
            y2 = hist.GetYaxis().GetBinUpEdge(j)

            #Draw lines around the bin
            lineTop = ROOT.TLine(x1, y2, x2, y2)    #Top line
            lineBottom = ROOT.TLine(x1, y1, x2, y1) #Bottom line
            lineLeft = ROOT.TLine(x1, y1, x1, y2)   #Left line
            lineRight = ROOT.TLine(x2, y1, x2, y2)  #Right line

            lineTop.SetLineColor(ROOT.kBlack)
            lineBottom.SetLineColor(ROOT.kBlack)
            lineLeft.SetLineColor(ROOT.kBlack)
            lineRight.SetLineColor(ROOT.kBlack)

            lineTop.SetLineWidth(5)
            lineBottom.SetLineWidth(5)
            lineLeft.SetLineWidth(5)
            lineRight.SetLineWidth(5)
           
            lineTop.Draw("SAME")
            lineBottom.Draw("SAME")
            lineLeft.Draw("SAME")
            lineRight.Draw("SAME")
            #canvas.Update()

    ROOT.gPad.Modified()
    ROOT.gPad.Update()

def plot(proj=False,asimov=False):
    #make plots of nominal input histograms used for the background fit
    #now do signal
    files = [f for f in glob('selection/XHYbbWWselection*.root')] #loop through all files
    if asimov:
        files.append('/uscms/home/mhesford/nobackup/XHYbbWW/CMSSW_14_1_0_pre4/src/XHYbbWW_BackgroundEstimate/PseudoDataToy_Poisson_method_seed12345_unblindedFail.root')
    tagHists = {'tag_SR':HistGroup('tag_SR'),'tag_ttCR':HistGroup('tag_ttCR')}

    style = setTDRStyle()
    style.cd() # style for 2D histograms

    cms_text = ROOT.TLatex()
    cms_text.SetTextSize(0.04)
    cms_text.SetTextFont(42)
    cms_text.SetTextColor(ROOT.kBlack)
    cms_text.SetTextAlign(13)

    lumi_text = ROOT.TLatex()
    lumi_text.SetTextSize(0.03)
    lumi_text.SetTextFont(42)
    lumi_text.SetTextAlign(31)

    for f in files:
        if 'JES' in f or 'JER' in f or 'JMS' in f or 'JMR' in f: continue
        print(f)
        can = ROOT.TCanvas('c','c',2800,2100)
        inFile=ROOT.TFile.Open(f,'READ')
        setname=f.split('/')[-1].split('_')[1].split('.')[0]
        histlist=[] #for 2D plots
        for key in inFile.GetListOfKeys():
            varname=key.GetName()
            if not varname.endswith('nominal') and 'Asimov' not in varname:
                continue
            inhist = inFile.Get(varname)  # get it from the file
            binsX = [600,1000,1200,1400,2000,4500]
            binsY = [100,400,600,750,1000,2500]
            outhist = Rebin2D(inhist,binsX,binsY)
            outhist.SetDirectory(0) # set the directory so hist is stored in memory and not as reference to TFile (this way it doesn't get tossed by python garbage collection when infile changes)      
            if 'Poisson' in setname:
                proc_title = 'Asimov'
            elif 'QCD' in setname:
                proc_title = 'QCD'
            else:
                proc_title = setname
            region_txt = outhist.GetTitle().split('_')
            outhist.SetTitle(f'{proc_title} ({region_txt[2]} {region_txt[1]})')
            outhist.SetName('{}_'.format(setname)+varname)
            outhist.GetXaxis().SetTitle('m_{X} (GeV)')
            outhist.GetYaxis().SetTitle('m_{Y} (GeV)')
            outhist.GetZaxis().SetTitle('Events / bin')
            outhist.GetXaxis().SetTitleOffset(1.3)
            can.cd()
            can.Clear
            
            outhist.Draw('COLZ TEXT') #draw histogram to canvas
            cms_text.DrawLatexNDC(0.18,0.87,'#splitline{#bf{CMS}}{#it{Work in Progress}}')
            lumi_text.DrawLatexNDC(0.83,0.91,'137 fb^{-1} (13 TeV)')

            can.Update()
            can.Print('selection/png/{}.png'.format(outhist.GetName()))
            
            if proj:
                #Now make projection on to M X,Y planes
                x_bins = [600,900,1100,1300,1500,2000,4500]
                pX = outhist.ProjectionX()
                pX_clone = ROOT.TH1D('x','x',len(x_bins)-1, array('d',x_bins))

                for b in range(1, pX.GetNbinsX() + 1):
                    bin_content = pX.GetBinContent(b)
                    bin_center = pX.GetBinCenter(b)
                    new_bin = pX_clone.FindBin(bin_center)
                    pX_clone.AddBinContent(new_bin, bin_content)

                CompareShapes(
                    outfilename = 'selection/png/{}_ProjX.pdf'.format(inhist.GetName()),
                    year = 1,
                    prettyvarname = 'm_{X}',
                    bkgs = {setname:pX_clone},
                    signals = {},
                    colors = {setname:ROOT.kRed},
                    names = {setname:setname},
                    scale = False,
                    stackBkg = False
                )

                y_bins = [100,300,500,700,1000,2500]
                pY = outhist.ProjectionX()
                pY_clone = ROOT.TH1D('y','y',len(y_bins)-1, array('d',y_bins))

                for b in range(1, pY.GetNbinsX() + 1):
                    bin_content = pY.GetBinContent(b)
                    bin_center = pY.GetBinCenter(b)
                    new_bin = pY_clone.FindBin(bin_center)
                    pY_clone.AddBinContent(new_bin, bin_content)
                
                CompareShapes(
                    outfilename = 'selection/png/{}_ProjY.pdf'.format(inhist.GetName()),
                    year = 1,
                    prettyvarname = 'm_{Y}',
                    bkgs = {setname:pY_clone},
                    signals = {},
                    colors = {setname:ROOT.kRed},
                    names = {setname:setname},
                    scale = False,
                    stackBkg = False
                )
    
        if not 'Data' and not 'QCD' and not 'PseudoDataToy' in f:
            cutflow = inFile.Get('cutflow_'+setname)
            cutflow.SetDirectory(0)
            c1 = ROOT.TCanvas('c1','c1')
            c1.cd()
            cutflow.Draw()
            c1.Print('selection/png/cutflow_{}.png'.format(setname))
    
        #plot isolation histogram for data-driven QCD
        print(f)
        if 'ttbar' in f or 'WJets' in f or 'ZJets' in f or 'DYJets' in f or 'ST' in f or 'Data' in f or 'QCD' in f:
            if 'PseudoDataToy' in f: continue
            process = f.split('/')[-1].split('_')[1].split('.')[0]
            for region in ['SR','ttCR']:
                grab_hist = inFile.Get('Wtag_{}'.format(region))
                if 'QCD' not in process:
                    grab_hist.Rebin(2)
                grab_hist.SetDirectory(0)
                tagHists['tag'+'_'+region][process] = grab_hist
     
    #Make not plots of isolation histograms
    processNames = OrderedDict([
       ('DataRun2','Data'), 
       ('QCD','QCD'),
       ('ttbar','ttbar'),
       ('ST','Single-top'),
       ('WJets','W+Jets'),
       ('ZJets','Z+Jets'),
       ('DYJets','Drell-Yan')
    ])
    processColors = OrderedDict([
       ('DataRun2',ROOT.kBlack),
       ('QCD',ROOT.kRed),
       ('ttbar',ROOT.kGreen-5),
       ('WJets',ROOT.kAzure+5),
       ('DYJets',ROOT.kViolet+6),
       ('ST',ROOT.kOrange-3),
       ('ZJets',ROOT.kGray)
    ])
    for var in tagHists.keys(): #iso<LEP>_<REGION>
        all_hists = tagHists[var]
        signals = OrderedDict([
            ('DataRun2',all_hists['DataRun2']),
            ('QCD',all_hists['QCD'])
        ])
        backgrounds = {process:all_hists[process] for process in ['ttbar','ST','WJets','ZJets','DYJets']}

        CompareShapes(
            outfilename = 'selection/png/{}.pdf'.format(var),
            year = 1,
            prettyvarname = 'W_{qq} D_{Wqq}',
            bkgs = backgrounds,
            signals = signals,
            colors = processColors,
            names = processNames,
            scale = False,
            stackBkg = True                    
        )    
    
if __name__=='__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--asimov',action='store_true',help='tell the program to plot the Asimov data, stored in a separate directory')
    args = parser.parse_args()
    plot(asimov = args.asimov)               


