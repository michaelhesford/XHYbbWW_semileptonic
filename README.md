# XHYbbWW semileptonic

Adapted from [XHYbbWW](https://github.com/ammitra/XHYbbWW) (Amitav Mitra).

Description: Semileptonic $X \to HY \to bbWW$ analysis, with $W^+W^- \to \ell \nu q\bar{q}$. Analysis is performed predominantly using [TIMBER](https://lucascorcodilos.com/TIMBER/index.html), and the background (QCD) estimate + calibration of Monte Carlo samples is done with [2DAlphabet](https://lucascorcodilos.com/2DAlphabet/) (see [XHYbbWW_BackgroundEstimate](https://github.com/michaelhesford/XHYbbWW_BackgroundEstimate) for this work).

## 1) Set up analysis workspace

### **Condor**
i) Use `condor/tar_env.sh` to create a tarball of the current environment and store it in the EOS.

ii) Create a symlink to the TIMBER condorHelper script: 
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



