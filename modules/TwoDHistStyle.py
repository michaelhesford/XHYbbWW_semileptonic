#Style function for 2D MX vs MY histograms 
#Modeled from: https://github.com/lcorcodilos/TIMBER/blob/183db7261ae93f9cdd68bf31859ebdb9b83b1630/TIMBER/Tools/Plot.py

import ROOT as rt
#from TIMBER.Tools.CMS import tdrstyle

def setTDRStyle(): #From: https://github.com/lcorcodilos/TIMBER/blob/183db7261ae93f9cdd68bf31859ebdb9b83b1630/TIMBER/Tools/CMS/tdrstyle.py
  tdrStyle =  rt.TStyle("tdrStyle","Style for P-TDR")

  #for the canvas:
  tdrStyle.SetCanvasBorderMode(0)
  tdrStyle.SetCanvasColor(rt.kWhite)
  tdrStyle.SetCanvasDefH(1200) #Height of canvas
  tdrStyle.SetCanvasDefW(1200) #Width of canvas
  tdrStyle.SetCanvasDefX(0)   #Position on screen
  tdrStyle.SetCanvasDefY(0)


  tdrStyle.SetPadBorderMode(0)#0
  #tdrStyle.SetPadBorderSize(Width_t size = 1)
  tdrStyle.SetPadColor(rt.kWhite)
  tdrStyle.SetPadGridX(False)
  tdrStyle.SetPadGridY(False)
  tdrStyle.SetGridColor(0)
  tdrStyle.SetGridStyle(3)
  tdrStyle.SetGridWidth(1)

#For the frame:
  tdrStyle.SetFrameBorderMode(0) #0
  tdrStyle.SetFrameBorderSize(1)
  tdrStyle.SetFrameFillColor(0)
  tdrStyle.SetFrameFillStyle(0)
  tdrStyle.SetFrameLineColor(1)
  tdrStyle.SetFrameLineStyle(1)
  tdrStyle.SetFrameLineWidth(1)
 
#For the histo:
  #tdrStyle.SetHistFillColor(1)
  #tdrStyle.SetHistFillStyle(0)
  tdrStyle.SetHistLineColor(1)
  tdrStyle.SetHistLineStyle(0)
  tdrStyle.SetHistLineWidth(1)
  #tdrStyle.SetLegoInnerR(Float_t rad = 0.5)
  #tdrStyle.SetNumberContours(Int_t number = 20)

  tdrStyle.SetEndErrorSize(2)
  #tdrStyle.SetErrorMarker(20)
  #tdrStyle.SetErrorX(0.)
  
  tdrStyle.SetMarkerStyle(20)
  
#For the fit/function:
  tdrStyle.SetOptFit(1)
  tdrStyle.SetFitFormat("5.4g")
  tdrStyle.SetFuncColor(2)
  tdrStyle.SetFuncStyle(1)
  tdrStyle.SetFuncWidth(1)

#For the date:
  tdrStyle.SetOptDate(0)
  # tdrStyle.SetDateX(Float_t x = 0.01)
  # tdrStyle.SetDateY(Float_t y = 0.01)

# For the statistics box:
  tdrStyle.SetOptFile(0)
  tdrStyle.SetOptStat(0) # To display the mean and RMS:   SetOptStat("mr")
  tdrStyle.SetStatColor(rt.kWhite)
  tdrStyle.SetStatFont(42)
  tdrStyle.SetStatFontSize(0.025)
  tdrStyle.SetStatTextColor(1)
  tdrStyle.SetStatFormat("6.4g")
  tdrStyle.SetStatBorderSize(1)
  tdrStyle.SetStatH(0.1)
  tdrStyle.SetStatW(0.15)
  # tdrStyle.SetStatStyle(Style_t style = 1001)
  # tdrStyle.SetStatX(Float_t x = 0)
  # tdrStyle.SetStatY(Float_t y = 0)

# Margins:
  tdrStyle.SetPadTopMargin(0.1) #0.05
  #rt.gPad.SetTopMargin(0.05)
  tdrStyle.SetPadBottomMargin(0.13) #0.13
  #rt.gPad.SetBottomMargin(0.12)
  tdrStyle.SetPadLeftMargin(0.16) #0.16
  #rt.gPad.SetLeftMargin(0.10)
  tdrStyle.SetPadRightMargin(0.17) #0.16
  #rt.gPad.SetRightMargin(0.17) #(0.17 most things), for some reason this is the only way that adjusting the margins will work

# For the Global title:

  tdrStyle.SetOptTitle(1)
  tdrStyle.SetTitleFont(42)
  tdrStyle.SetTitleColor(1)
  tdrStyle.SetTitleTextColor(1)
  tdrStyle.SetTitleFillColor(0) #10
  tdrStyle.SetTitleFontSize(0.04)
  tdrStyle.SetTitleH(0) # Set the height of the title box
  tdrStyle.SetTitleW(0) # Set the width of the title box
  #tdrStyle.SetTitleAlign(13)
  #tdrStyle.SetTitleOffset(0.2,'X')
  tdrStyle.SetTitleAlign(11)
  tdrStyle.SetTitleX(0.16) #PNetEfficiencies (0.1), PNetTemplates (0.17)
  #tdrStyle.SetTitleX(0.15) # Set the position of the title box
  tdrStyle.SetTitleY(0.90)#0.93 # Set the position of the title box
  tdrStyle.SetTitleStyle(1001)
  tdrStyle.SetTitleBorderSize(0) #1

# For the axis titles:

  tdrStyle.SetTitleColor(1, "XYZ")
  tdrStyle.SetTitleFont(42, "XYZ")
  tdrStyle.SetTitleSize(0.04, "XYZ") #0.06
  # tdrStyle.SetTitleXSize(Float_t size = 0.02) # Another way to set the size?
  # tdrStyle.SetTitleYSize(Float_t size = 0.02)
  #tdrStyle.SetTitleXOffset(1.5) #0.9
  #tdrStyle.SetTitleYOffset(0.7)
  tdrStyle.SetTitleOffset(1.1, "X")
  tdrStyle.SetTitleOffset(1.2, "Y")
  tdrStyle.SetTitleOffset(1.3, "Z")

# For the axis labels:

  tdrStyle.SetLabelColor(1, "XYZ")
  tdrStyle.SetLabelFont(42, "XYZ")
  tdrStyle.SetLabelOffset(0.007, "XYZ")
  tdrStyle.SetLabelSize(0.03, "XYZ")

# For the axis:

  tdrStyle.SetAxisColor(1, "XYZ")
  tdrStyle.SetStripDecimals(True)
  tdrStyle.SetTickLength(0.02, "XYZ")
  tdrStyle.SetNdivisions(510, "XYZ")
  tdrStyle.SetPadTickX(True)  # To get tick marks on the opposite side of the frame
  tdrStyle.SetPadTickY(True)

# Change for log plots:
  tdrStyle.SetOptLogx(0)
  tdrStyle.SetOptLogy(0)
  tdrStyle.SetOptLogz(0)

# Postscript options:
  tdrStyle.SetPaperSize(20.,20.)
  # tdrStyle.SetLineScalePS(Float_t scale = 3)
  # tdrStyle.SetLineStyleString(Int_t i, const char* text)
  # tdrStyle.SetHeaderPS(const char* header)
  # tdrStyle.SetTitlePS(const char* pstitle)

  # tdrStyle.SetBarOffset(Float_t baroff = 0.5)
  #tdrStyle.SetBarWidth(0.5)
  tdrStyle.SetPaintTextFormat("1.1f") #(const char* format = "g")
  # tdrStyle.SetPalette(Int_t ncolors = 0, Int_t* colors = 0)
  # tdrStyle.SetTimeOffset(Double_t toffset)
  # tdrStyle.SetHistMinimumZero(kTRUE)

  tdrStyle.SetHatchesLineWidth(5)
  tdrStyle.SetHatchesSpacing(0.05)

  return tdrStyle
