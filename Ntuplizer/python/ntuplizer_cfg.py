import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.dump=cms.EDAnalyzer('EventContentAnalyzer')
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.EventContent.EventContent_cff')

process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate 
fileFormat = 'AOD'

if fileFormat == 'AOD' :
    dataFormat = DataFormat.AOD
else :
    dataFormat = DataFormat.MiniAOD

switchOnVIDElectronIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = [
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_Trig_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_50ns_Trig_V1_cff',
                ]

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)



process.ntuple = cms.EDAnalyzer('Ntuplizer',
                                   inputFileFormat = cms.untracked.string(fileFormat),
                                   beamSpot = cms.InputTag('offlineBeamSpot'),
                                   # input collection names AOD
                                   electronsAOD = cms.InputTag('gedGsfElectrons'),
                                   verticesAOD = cms.InputTag('offlinePrimaryVertices'),
                                   conversionsAOD = cms.InputTag('allConversions'),
                                   genParticlesAOD = cms.InputTag('genParticles'), 
                                   genEventInfoProductAOD = cms.InputTag('generator'),
                                   PFMETAOD = cms.InputTag('pfMet'),                                
    
                                   electronsMiniAOD = cms.InputTag('slimmedElectrons'),
                                   verticesMiniAOD = cms.InputTag('offlineSlimmedPrimaryVertices'),
                                   conversionsMiniAOD = cms.InputTag('reducedEgamma:reducedConversions'),
                                   genParticlesMiniAOD = cms.InputTag('prunedGenParticles'), 
                                   PFMETMiniAOD = cms.InputTag('slimmedMETs'),                                

                                   HLTTag          = cms.InputTag('TriggerResults','','HLT'),
                                   isMC = cms.bool(True),
                                   MVAId  = cms.VInputTag(),
                            
                                   electronID1 = cms.InputTag('electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig25nsV1Values'),
                                   electronID2 = cms.InputTag('electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig50nsV1Values'),

)
fileNameForSample = 'ntuple'


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/home/llr/cms/pigard/data/testFiles/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v3_AODSIM_file1.root')

#fileNames = cms.untracked.vstring('file:/home/llr/cms/pigard/GluGluHToZZTo4L_M125_13TeV_powheg_JHUgen_pythia8_MINIAODSIM.root')

)
 
#process.GlobalTag.globaltag = ' MCRUN2_74_V8::All'
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag

process.GlobalTag = GlobalTag(process.GlobalTag, 'MCRUN2_74_V9', '') # MCRUN2_74_V8

process.TFileService = cms.Service("TFileService", fileName = cms.string(fileNameForSample + '.root') )

process.MessageLogger.cerr.threshold = 'WARNING'
process.MessageLogger.categories.append('Demo')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(-1)
)
process.MessageLogger = cms.Service(
    "MessageLogger",
    destinations = cms.untracked.vstring(
        'LOG',
        'critical'
    ),
    LOG = cms.untracked.PSet(
        threshold  = cms.untracked.string('WARNING'), # DEBUG 
        filename  = cms.untracked.string(fileNameForSample  + '.log')
    ),
    debugModules = cms.untracked.vstring('*'),
       
    statistics     = cms.untracked.vstring('STAT'),
        STAT = cms.untracked.PSet(
            threshold = cms.untracked.string('WARNING'),
            filename  = cms.untracked.string(fileNameForSample  + '_stats.log')
    )
)

process.p = cms.Path(process.ntuple)#*process.dump)
