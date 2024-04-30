# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line args: myNanoProdMc2018 -s NANO --mc --eventcontent NANOAOD --datatier NANOAOD --no_exec --conditions 106X_upgrade2018_realistic_v16_L1v1 --era Run2_2018,run2_nanoAOD_106Xv2 --customise_commands=process.add_(cms.Service("InitRootHandlers", EnableIMT = cms.untracked.bool(False)))
import importlib

import FWCore.ParameterSet.Config as cms

from LLStaus_Run2.Production.arg_config import *
args = get_args()

d_procConfig = {
    "Data": {
        "2016": {
            "condition": "auto:run2_data",
            "era": "Run2_2016",
            "eramodifier": "run2_nanoAOD_106Xv2",
        },
        "2017":{
            "condition": "auto:run2_data",
            "era": "Run2_2017",
            "eramodifier": "run2_nanoAOD_106Xv2",
        },
        "2018":{
            "condition": "auto:run2_data",
            "era": "Run2_2018",
            "eramodifier": "run2_nanoAOD_106Xv2",
        },
    },
    
    "MC": {
        # 2016 conditions not checked yet; just a placeholder for now
        "2016": {
            "condition": "auto:phase1_2016_realistic",
            "era": "Run2_2016",
            "eramodifier": "run2_nanoAOD_106Xv2",
        },
        "2017":{
            "condition": "auto:phase1_2017_realistic",
            "era": "Run2_2017",
            "eramodifier": "run2_nanoAOD_106Xv2",
        },
        "2018":{
            "condition": "auto:phase1_2018_realistic",
            "era": "Run2_2018",
            "eramodifier": "run2_nanoAOD_106Xv2",
        },
        "2022EE":{
            "condition": "auto:phase1_2022_realistic_postEE",
            "era": "Run3",
            "eramodifier": "run3_nanoAOD_122"
        },
        "2022":{
            "condition": "auto:phase1_2022_realistic",
            "era": "Run3",
            "eramodifier": "run3_nanoAOD_122"
        }
    }
}

#isMC = (args.sampleType == "MC")
isMC = True
condition_str = d_procConfig[args.sampleType][args.era]["condition"]
era_str = d_procConfig[args.sampleType][args.era]["era"]
eramodifier_str = d_procConfig[args.sampleType][args.era]["eramodifier"]

era_cff = importlib.import_module(f"Configuration.Eras.Era_{era_str}_cff")
era = getattr(era_cff, era_str)

eramodifier_cff = importlib.import_module(f"Configuration.Eras.Modifier_{eramodifier_str}_cff")
eramodifier = getattr(eramodifier_cff, eramodifier_str)

process = cms.Process("NANO", era, eramodifier)

# import of standard configurations
process.load("Configuration.StandardSequences.Services_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.EventContent.EventContent_cff")
process.load("SimGeneral.MixingModule.mixNoPU_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("PhysicsTools.NanoAOD.nano_cff")
process.load("Configuration.StandardSequences.EndOfProcess_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, condition_str, "")

process.MessageLogger.cerr.enableStatistics = True

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(args.maxEvents)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(),
    secondaryFileNames = cms.untracked.vstring(),
    #lumisToProcess = cms.untracked.VLuminosityBlockRange('1:16-1:16', '1:21-1:21')
    )
from LLStaus_Run2.Production.readFileList import *
if len(args.inputFiles) > 0:
    addList(process.source.fileNames, args.inputFiles, fileNamePrefix=args.fileNamePrefix)
    args.outFile = args.inputFiles[0].split("/")[-5]+".root"
elif len(args.sourceFile) > 0:
    addList(process.source.fileNames, args.inputFiles, fileNamePrefix=None)
    args.outFile = args.inputFiles[0].split("/")[-5]+".root"

if len(args.lumiFile) > 0:
    import FWCore.PythonUtilities.LumiList as LumiList
    process.source.lumisToProcess = LumiList.LumiList(filename = args.lumiFile).getVLuminosityBlockRange()

if args.eventRange != "":
    process.source.eventsToProcess = cms.untracked.VEventRange(re.split(",", args.eventRange))

if args.maxEvents > 0:
    process.maxEvents.input = args.maxEvents

process.options = cms.untracked.PSet()

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string("myNanoProd{args.sampleType}{args.era}"),
    name = cms.untracked.string("Applications"),
    version = cms.untracked.string("$Revision: 1.19 $")
)

## let's try sara
if (args.era == "2022"): 
    eramodifier.toModify(
        process.linkedObjects, jets="finalJets"
    )

    # Replace AK4 Puppi with AK4 CHS for Run-2
    _nanoTableTaskCommonRun2 = process.nanoTableTaskCommon.copy()
    _nanoTableTaskCommonRun2.replace(process.jetPuppiTask, process.jetTask)
    _nanoTableTaskCommonRun2.replace(process.jetPuppiForMETTask, process.jetForMETTask)
    _nanoTableTaskCommonRun2.replace(process.jetPuppiTablesTask, process.jetTablesTask)
    eramodifier.toReplaceWith(
        process.nanoTableTaskCommon, _nanoTableTaskCommonRun2
    )

    eramodifier.toModify(
        process.ptRatioRelForEle, srcJet="updatedJets"
    )   
    eramodifier.toModify(
        process.ptRatioRelForMu, srcJet="updatedJets"
    )

process.trigoutput = cms.EDFilter("TriggerResultsFilter",
      
    ### Filter events on trigger
    #SelectEvents = cms.untracked.PSet(
    triggerConditions = cms.vstring("HLT_PFMET120_PFMHT120_IDTight_v20",
                                      "HLT_PFMET130_PFMHT130_IDTight_v20",
                                      "HLT_PFMET140_PFMHT140_IDTight_v20",
                                      "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v20",
                                      "HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v19", 
                                      "HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v19",
                                      "HLT_PFMET120_PFMHT120_IDTight_PFHT60_v9",
                                      "HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight_v20",
                                      "HLT_PFMETTypeOne140_PFMHT140_IDTight_v11",
                                      "HLT_MET105_IsoTrk50_v9",
                                      "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF_v1",
                                      "HLT_MET120_IsoTrk50_v9",
                                      "HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS35_L2NN_eta2p1_CrossL1_v1",
                                      "HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_CrossL1_v1",
                                      "HLT_Ele30_WPTight_Gsf_v1",
                                      'HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1_v1',
                                      'HLT_DoubleMediumChargedIsoDisplacedPFTauHPS32_Trk1_eta2p1_v1',
                                      'HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_eta2p1_v1'
            ),
    hltResults = cms.InputTag("TriggerResults", "", "HLT"),
    l1tResults = cms.InputTag("gtStage2Digis"),
    throw = cms.bool(True)
    #),
)

# Output definition
assert(args.disTauTagOutputOpt in [0, 1, 2])

if isMC :
    outputCommands = process.NANOAODSIMEventContent.outputCommands
else :
    outputCommands = process.NANOAODEventContent.outputCommands

if args.disTauTagOutputOpt == 1 :
    outputCommands += cms.untracked.vstring("keep *_TriggerResults_*_*")
    args.outFile = args.outFile.replace(".root", "_with-disTauTagScore.root")

elif args.disTauTagOutputOpt == 2 :
    
    outputCommands = cms.untracked.vstring(
        "drop *",
        #"keep *_*_*disTauTag*_*",
        #"keep nanoaodFlatTable_*Table_*_*",
        #"keep nanoaodFlatTable_*Table_*_*",
        "keep nanoaodFlatTable_jetTable_*_*",
        "keep *_TriggerResults_*_*",
    )
    
    args.outFile = args.outFile.replace(".root", "_only-disTauTagScore.root")

process.NANOAODoutput = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm = cms.untracked.string("LZMA"),
    compressionLevel = cms.untracked.int32(9),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string("NANOAODSIM") if isMC else cms.untracked.string("NANOAOD"),
        filterName = cms.untracked.string("")
    ),
    fileName = cms.untracked.string("With_trigselec.root"),
    outputCommands = outputCommands,
)



# Additional output definition

# Other statements

# Path and EndPath definitions
if isMC :
    process.nanoAOD_step = cms.Path(process.nanoSequenceMC)
else :
    process.nanoAOD_step = cms.Path(process.nanoSequence)

process.endjob_step = cms.EndPath(process.endOfProcess)
process.NANOAODoutput_step = cms.EndPath(process.trigoutput*process.NANOAODoutput)

# Schedule definition
process.schedule = cms.Schedule(process.nanoAOD_step,process.endjob_step,process.NANOAODoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# customisation of the process.

if isMC :
    from PhysicsTools.NanoAOD.nano_cff import nanoAOD_customizeCommon 
    process = nanoAOD_customizeCommon(process)

else :
    from PhysicsTools.NanoAOD.nano_cff import nanoAOD_customizeCommon
    process = nanoAOD_customizeCommon(process)

from LLStaus_Run2.Production.customize_nanoaod_eventcontent_cff import *
customize_process_and_associate(process, isMC = isMC, disTauTagOutputOpt = args.disTauTagOutputOpt)


# End of customisation functions

# Customisation from command line

process.add_(cms.Service("InitRootHandlers", EnableIMT = cms.untracked.bool(False)))
# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion


# Debug EDM
if (args.debugEDM) :
    
    process.out = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string("debugEDM.root")
    )
    
    process.output_step = cms.EndPath(process.out)
    process.schedule.extend([process.output_step])
