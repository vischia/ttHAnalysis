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


   * STATUS:
      - at the moment, it just runs on 2017 data and MC with old corrections, and locally (not on grid)
      - Run with ```cmsRun test/ttHConfiguration.py``` (eventually change the value of runOnData to switch between running on data and on MC, and uncomment different MCs)

   * NOTE: all temporary synchronization numbers and details on procedures (CMS private numbers) in private repo ( https://gitlab.cern.ch/vischia/uclouvain_tthsync/ )
   
   * TODO:
      - Update JEC, electron corrections, muon corrections -> Hesham
      - Synchronization -> Hesham.
      - Update datasets database (RunIIFall17MiniAODv2 branch) -> Pietro
      - Implement baseline event selection
      - Implement datacards creator ->
      - Import unfolding software -> Pietro
