from glob import glob
import subprocess, os
from TIMBER.Tools.Common import ExecuteCmd
from TIMBER.Analyzer import HistGroup
import ROOT
from collections import OrderedDict
from array import array
import math
import sympy as sp

def Rebin2D(oldHist, binsX, binsY):
    nbinsX = len(binsX) - 1
    nbinsY = len(binsY) - 1
    binsX_array = array('d',binsX)
    binsY_array = array('d',binsY)
    name = oldHist.GetName()
    title = oldHist.GetTitle()
    oldHist.SetName('old')
    oldHist.SetTitle('old')
    newHist = ROOT.TH2D(name,title,nbinsX,binsX_array,nbinsY,binsY_array)

    for i in range(1,oldHist.GetNbinsX()+1):
        for j in range(1,oldHist.GetNbinsY()+1):
            content = oldHist.GetBinContent(i,j)
            x = oldHist.GetXaxis().GetBinCenter(i)
            y = oldHist.GetYaxis().GetBinCenter(j)
            newHist.Fill(x,y,content)

    return newHist

def make_QCD(backgrounds):

    binsX = [600,1000,1200,1400,2000,4500] #Values for re-binning X axis 
    binsY = [100,400,600,750,1000,2500] #Values for re-binning Y axis
    data_file = ROOT.TFile.Open('selection/XHYbbWWselection_DataRun2.root','READ')
    QCD_2d = HistGroup('QCD_2d')
    QCD_tag = HistGroup('QCD_tag')
    for region in ['SR','ttCR']:
        data_2d_in = data_file.Get(f'MXvMY_fail_{region}__nominal')
        print(f'Rebinning data histogram for {region} events')
        data_2d = Rebin2D(data_2d_in,binsX,binsY) #Rebin data histogram
        data_2d.SetDirectory(0)
        data_2d.SetTitle(f'MXvMY_fail_{region}__nominal')
        data_tag = data_file.Get(f'Wtag_{region}')
        data_tag = data_tag.Rebin(2)
        data_tag.SetDirectory(0)
        data_tag.SetTitle(f'Wtag_{region}')
        for bkg in backgrounds:
            bkg_file = ROOT.TFile.Open(f'selection/XHYbbWWselection_{bkg}.root')
            bkg_2d_in = bkg_file.Get(f'MXvMY_fail_{region}__nominal')
            if bkg == 'ST':
                for x in range(1,bkg_2d_in.GetNbinsX()+1):
                    for y in range(1,bkg_2d_in.GetNbinsY()+1):
                        print(f'{region} ST fail bin {x},{y}: ',bkg_2d_in.GetBinContent(x,y))
                        if math.isnan(bkg_2d_in.GetBinContent(x,y)):
                            bkg_2d_in.SetBinContent(x,y,0)

            print(f'Rebinning {bkg} histogram for {region} events')
            bkg_2d = Rebin2D(bkg_2d_in,binsX,binsY) #Rebin background histogram
            bkg_tag = bkg_file.Get(f'Wtag_{region}')
            bkg_tag = bkg_tag.Rebin(2)
            if bkg == 'ST': 
                for x in range(1,bkg_tag.GetNbinsX()+1):
                    print(f'{region} ST tag bin {x}: ',bkg_tag.GetBinContent(x))
                    if math.isnan(bkg_tag.GetBinContent(x)):
                        bkg_tag.SetBinContent(x,0)


            data_2d.Add(bkg_2d,-1)
            data_tag.Add(bkg_tag,-1)

        #Set any negative bins to 0
        for x in range(1,data_2d.GetNbinsX()+1):
            for y in range(1,data_2d.GetNbinsY()+1):
                if data_2d.GetBinContent(x,y) < 0:
                    data_2d.SetBinContent(x,y,0)
                    print(f'Setting negative QCD bin ({str(x)},{str(y)}) to 0 in {region} events')
            
            
        QCD_2d[f'{region}'] = data_2d
        QCD_tag[f'{region}'] = data_tag
    '''
    #Now extract pass/fail ratio for isolation, with the updated method
    tag_ratios = OrderedDict()
    QCD_pass = HistGroup('QCD_pass')
    for region in ['SR','ttCR']:
        #Get the parameters from the tf
        #param_file = open(f'/uscms/home/mhesford/nobackup/XHYbbWW/CMSSW_14_1_0_pre4/src/XHYbbWW_BackgroundEstimate/QCDrpf_fits/{region}_fits/{region}-1x0_area/rpf_params_Background_pf_1x0_fitb.txt','r').readlines()
        param_file = open(f'/uscms/home/mhesford/nobackup/XHYbbWW/CMSSW_14_1_0_pre4/src/XHYbbWW_BackgroundEstimate/QCDrpf_fits/QCDrpf_fits/QCDrpf_area/rpf_params_Background_{region}_pf_1x0_fitb.txt','r').readlines()
        p0 = 0.1*float(param_file[0].split()[1])
        p0_err = 0.1*float(param_file[0].split()[3])
        p1 = 0.1*float(param_file[1].split()[1])
        p1_err = 0.1*float(param_file[1].split()[3])
        params = [p0,p0_err,p1,p1_err]

        yearlyW = {
            '16APV':0.637,
            '16':0.642,
            '17':0.579,
            '18':0.59
        }
        lumis = {
            '16APV':19.51,
            '16':16.81,
            '17':41.48,
            '18':59.83
        }
        lumi_total = 137.63

        x = sp.symbols('x')
        y = (p1-p1_err)*x + p0
        nfail = 0
        npass = 0
        for year in lumi:
            wp = yearlyW[year]
            lumi = lumis[year]
            dfail = sp.integrate(y,(x,0,wp))
            dpass = sp.integrate(y,(x,wp,1))
            
            nfail += dfail*(lumi/lumi_total)
            npass += dpass*(lumi/lumi_total)

        ratio = npass/nfail
        tag_ratios[f'{region}'] = ratio

        pass_hist = QCD_2d[f'{region}'].Clone('QCD_pass')
        pass_hist.Scale(ratio)
        pass_hist.SetName(f'MXvMY_pass_{region}__nominal')
        pass_hist.SetTitle(f'MXvMY_pass_{region}__nominal')
        pass_hist.SetDirectory(0)
        QCD_pass[region] = pass_hist
    '''
    out1 = ROOT.TFile.Open(f'selection/XHYbbWWselection_QCD.root','RECREATE')
    QCD_2d.Do('Write')
    #QCD_pass.Do('Write')
    #QCD_tag.Do('Write')
    out1.Close()
    '''
    out2 = open(f'selection/tag_ratios.txt','w')
    for r in tag_ratios.keys():
        out2.write('Pass/fail for {}: {}\n'.format(r,tag_ratios[r]))
    out2.close()
    '''


if __name__=='__main__':

    groupnames = ['ttbar','WJets','ZJets','DYJets','ST']
    make_QCD(groupnames)



