# XHYbbWW semileptonic

Adapted from [XHYbbWW](https://github.com/ammitra/XHYbbWW) (Amitav Mitra).

Description: Semileptonic $X \to HY \to bbWW$ analysis, with $W^+W^- \to \ell \nu q\bar{q}$. Analysis is performed predominantly using [TIMBER](https://lucascorcodilos.com/TIMBER/index.html), and the background (QCD) estimate + calibration of Monte Carlo samples is done with [2DAlphabet](https://lucascorcodilos.com/2DAlphabet/) (see [XHYbbWW_BackgroundEstimate](https://github.com/michaelhesford/XHYbbWW_BackgroundEstimate) for this work).

## 1) Set up analysis workspace

### **Condor**
i) Use `condor/tar_env.sh` to create a tarball of the current environment and store it in the EOS.

ii) Create a symlink to the TIMBER CondorHelper script: 
```
ln -s $TIMBERPATH/TIMBER/Utilities/Condor/CondorHelper.py
```

iii) To submit jobs using condor: 
```
python CondorHelper.py -r <RUN SCRIPT>.sh -a <ARGUMENTS FILE>.txt -i "<LOCAL SCRIPTS>"
```

The `-i` option is used to incorporate any scripts which have been created/modified since the last time the environment tarball was updated.


## 2) Grab raw data/MC files using CERN's Data Aggregation System (DAS)

Perform

```
python raw_nano/get_Data.py
python raw_nano/get_MC.py 
```

to get locations of NanoAOD data/MC and store in `.txt` files in the `raw_nano` directory

## 3) Create pileup distributions for pileup weights

Individual jobs: `python XHYbbWW_pileup.py -s <SETNAME> -y <YEAR>`

Batch: run `python condor/pileup_args.py` to produce the arguments file and then

```
python CondorHelper.py -r condor/run_pileup.sh -a condor/pileup_args.txt -i "XHYbbWWpileup.py raw_nano/"
```

to submit the pileup jobs using condor. The outputs are then collected into one local ROOT file called `XHybbWWpileup.root` using

```
python scripts/get_pileup_file.py
```

## 3) Perform snapshots

Snapshots are skims performed on the raw NanoAOD files in order to isolate a smaller set of events/variables which are interesting for the analysis. To do this we impose a very loose set of cuts on the data, and save any desired columns of the resulting RDataFrame in a separate `.root` file. Any future scripts can be run on these snapshot files, saving an enormous amount of time.

Individual job: `python XHYbbWW_snapshot.py -s <SETNAME> -y <YEAR> -j <IJOB> -n <NJOBS>`

The `-n` argument is used to break the initial `.txt` file of ROOT files into a specified number of smaller chunks. The specific chunk to run on is given by the `-i` argument. 

Batch: 

```
python condor/snapshot_args.py
python CondorHelper.py -r condor/run_snapshot.sh -a condor/snapshot_args.txt -i "XHYbbWW_snapshot.py XHYbbWW_class.py HWWmodules.cc"
``` 

Locations of the otuputs can be collected into `.txt` files using `python snapshots/get_snapshots.py`.

## 4) Run basic studies

This step produces basic kinematic plots and N-1 plots for key variables used in the analysis. Given an ensemble of cuts on different variables, we can use the `Nminus1()` function of TIMBER's`analyzer` class to create a set of nodes in which all but one cut is applied. We then plot each variable, for signal and background MC, in the node where cuts on all variables besides itself are applied. These so-called N-1 plots provide insight into the unique contributions of specific cuts on individual variables; we use them to help select working points which maximize the ratio $s/\sqrt{b}$, where $s$ and $b$ are the signal and background yields, respectively.

Individual job: `python XHYbbWW_studies.py -s <SETNAME> -y <YEAR> -j <IJOB> -n <NJOBS>`
Batch: 

```
python condor/studies_args.py
python CondorHelper.py -r condor/run_studies.sh -a condor/studies_args.txt -i "XHYbbWW_studies.py XHYbbWW_class.py HWWmodules.cc"
```

Histograms are stored in `.root` files, which can be grabbed from the EOS and stored in `plots/studies/` using
```
python plots/get_all.py -f studies -g <LIST OF SETNAMES> --combine_years
```

The `--combine_years` option uses the `hadd` function to combine outputs from different years for each setname into full Run2 plots.

## 5) Make trigger efficiencies

The set of triggers used in the analysis is listed below and is also stored in `XHYbbWWconfig.json`.

| Year | Triggers |
| ---- | -------- |
| 2016 | SingleMuon: `HLT_IsoMu24` `HLT_IsoTkMu24` `HLT_Mu50` `HLT_TkMu50` \\ SingleElectron: `HLT_Ele27_WPTight_Gsf` `HLT_Ele45_WPLoose_Gsf` `HLT_Ele115_CaloIdVT_GsfTrkIdT` \\ SinglePhoton: `HLT_Photon175` |
| 2017 | SingleMuon: `HLT_IsoMu27`HLT_Mu50` SingleElectron: `HLT_Ele35_WPTight_Gsf``HLT_Ele32_WPTight_Gsf_L1DoubleEG``HLT_Ele115_CaloIdVT_GsfTrkIdT` SinglePhoton: `HLT_Photon200` |
| 2018 | SingleMuon: 'HLT_IsoMu24''HLT_Mu50' EGamma: `HLT_Ele32_WPTight_Gsf``HLT_Ele115_CaloIdVT_GsfTrkIdT``HLT_Photon200` |

## 5) Run selection + make template histograms




