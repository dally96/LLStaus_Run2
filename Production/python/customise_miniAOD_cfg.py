import importlib

import FWCore.ParameterSet.Config as cms

from LLStaus_Run2.Production.arg_config import *
args = get_args()

def addTauTagger(process):

  process.disTauTag = cms.EDProducer(
        "DisTauTag",
        graphPath = cms.string("data/particlenet_v1_a27159734e304ea4b7f9e0042baa9e22.pb"),
        jets = cms.InputTag("slimmedJets"),
        pfCandidates = cms.InputTag('packedPFCandidates'),
        save_inputs  = cms.bool(False)
  )
  process.distau_path = cms.Path(process.disTauTag)
  
  
process = cms.Process("TAUTAGGER")
process.load("Configuration.StandardSequences.MagneticField_cff") # for CH reco
process.load("Configuration.Geometry.GeometryRecoDB_cff")
addTauTagger(process)


process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
if args.era == 2022:
  process.GlobalTag = GlobalTag(process.GlobalTag, '130X_mcRun3_2022_realistic_v5', '')
if args.era == 2018:
  process.GlobalTag = GlobalTag(process.GlobalTag, '106X_upgrade2018_realistic_v16_L1v1', '')

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(),
    secondaryFileNames = cms.untracked.vstring()
    )
from LLStaus_Run2.Production.readFileList import *
if len(args.inputFiles) > 0:
    addList(process.source.fileNames, args.inputFiles, fileNamePrefix=args.fileNamePrefix)
    args.outFile = args.inputFiles[0].split("/")[-5]+".root"
elif len(args.sourceFile) > 0:
    addList(process.source.fileNames, args.inputFiles, fileNamePrefix=None)
    args.outFile = args.inputFiles[0].split("/")[-5]+".root"

## input files
#readFiles = cms.untracked.vstring()
#secFiles = cms.untracked.vstring()
#process.source = cms.Source(
#    "PoolSource", fileNames=readFiles, secondaryFileNames=secFiles)

#readFiles.extend([
#    '/store/group/lpcdisptau/Staus_M_100_100mm_13p6TeV_Run3Summer22EE/MINIAODSIM/230503_144828/0000/TSG-Run3Summer22EEMiniAOD_inMINIAODSIM_1.root',
#    ])

# limit the number of events to be processed
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32( 20000 )
)

## output
process.load('Configuration.EventContent.EventContent_cff')
process.output = cms.OutputModule('PoolOutputModule',
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('MINIAODSIM'),
        filterName = cms.untracked.string('')
    ),
    dropMetaData = cms.untracked.string('ALL'),
    eventAutoFlushCompressedSize = cms.untracked.int32(-900),
    fastCloning = cms.untracked.bool(False),
    fileName = cms.untracked.string('file:TSG-Run3Summer22EEMiniAOD_inMINIAODSIM.root'),
    outputCommands = process.MINIAODSIMEventContent.outputCommands,
    overrideInputFileSplitLevels = cms.untracked.bool(True),
    splitLevel = cms.untracked.int32(0)
)

process.out = cms.EndPath(process.output)
process.schedule = cms.Schedule(process.distau_path, process.out)

process.output.fileName = cms.untracked.string("/eos/user/d/dally/DisplacedTauAnalysis/"+args.inputFiles[0].split("/")[-5]+"_"+args.inputFiles[0].split(".")[0].split("/")[-1].split("_")[-1]+".root")
process.output.outputCommands.append('keep *_genParticlePlusGeant_*_*')
process.output.outputCommands.append('keep *_disTauTag_*_*')
