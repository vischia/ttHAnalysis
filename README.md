# ttH Analyzer

   * Note: must be installed after the [Framework](https://github.com/cp3-llbb/Framework)

```
cd ${CMSSW_BASE}/src
git clone -o upstream git@github.com:vischia/ttHAnalysis.git cp3_llbb/ttHAnalysis
cd ${CMSSW_BASE}/src/cp3_llbb/ttHAnalysis
source setup.sh
cd ${CMSSW_BASE}/src/
scram b -j 20
```


   * TODO:
      - Update JEC, electron corrections, muon corrections
      - Implement baseline event selection
      - Implement datacards creator
      - Import unfolding software
