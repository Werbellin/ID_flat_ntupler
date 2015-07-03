import FWCore.ParameterSet.Config as cms

#process.load("RecoTracker.Configuration.python.RecoTracker_cff")
process = cms.Process("Demo")

process.dump=cms.EDAnalyzer('EventContentAnalyzer')
process.load("Configuration.StandardSequences.MagneticField_cff")

#process.load("Configuration.StandardSequences.Geometry_cff")
#process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.EventContent.EventContent_cff')

process.load('Configuration.StandardSequences.Reconstruction_cff')
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")




process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )


#process.load("RecoTracker.TrackProducer.KfTrackFromGsfTrack_cfi")
#process.KfTrackFromGsfTrack.src = cms.InputTag("electronGsfTracks")
#theCKFTrajectoryInput = cms.string('KfTrackFromGsfTrack')
#process.KFFittingSmootherWithOutliersRejectionAndRK.EstimateCut = 10000

#process.load("RecoTracker.TrackProducer.TrackProducerOwnLast5_cfi")
#process.TrackProducerOwnLast5.GSFTrajectoryInput  = cms.string('GsfTrackRefitter')

#process.load("RecoTracker.TrackProducer.TrackRefittersMod_cff")
#process.TrackRefitter.src = cms.InputTag("generalTracks")

#process.load("RecoTracker.TrackProducer.GsfTrackRefitterMod_cff")
#process.GsfTrackRefitter.src = cms.InputTag("electronGsfTracks")
#process.GsfTrackRefitter.TrajectoryInEvent = cms.bool(True)

#theGSFTrajectoryInput = cms.string('GsfTrackRefitter')



#therunGsfRefitter      = cms.bool(True)
#therunKfWithGsfRefitter = cms.bool(True)

process.ntuple = cms.EDAnalyzer('Ntuplizer',
                                   inputFileFormat = cms.untracked.string('AOD'),
                                   beamSpot = cms.InputTag('offlineBeamSpot'),
                                   # input collection names AOD
                                   electronsAOD = cms.InputTag('gedGsfElectrons'),
                                   verticesAOD = cms.InputTag('offlinePrimaryVertices'),
                                   conversionsAOD = cms.InputTag('allConversions'),
                                   genParticlesAOD = cms.InputTag('genParticles'), 
                                   PFMETAOD = cms.InputTag('pfMet'),                                
    
                                   electronsMiniAOD = cms.InputTag('slimmedElectrons'),
                                   verticesMiniAOD = cms.InputTag('offlineSlimmedPrimaryVertices'),
                                   conversionsMiniAOD = cms.InputTag('reducedEgamma:reducedConversions'),
                                   genParticlesMiniAOD = cms.InputTag('prunedGenParticles'), 
                                   PFMETMiniAOD = cms.InputTag('slimmedMETs'),                                



                                   HLTTag          = cms.InputTag('TriggerResults','','HLT'),
                                   isMC = cms.bool(True),
                                   ispythia6 = cms.bool(True),
                                   #runGsfRefitter      = therunGsfRefitter,
                                   #GSFTrajectoryInput  = theGSFTrajectoryInput,
                                   #runKfWithGsfRefitter = therunKfWithGsfRefitter,
                                   #CKFTrajectoryInput  = theCKFTrajectoryInput,

                                   MVAId  = cms.VInputTag()#"mvaNonTrigV025nsPHYS14", "mvaNonTrigV025nsPHYS14FIX")
)
fileNameForSample = 'ntuple'


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.source = cms.Source("PoolSource",
#    fileNames = cms.untracked.vstring(
        #'root://xrootd.unl.edu//store/user/ppigard/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6_Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1_MODMINIAOD_All/ff0f7a48dd492a8600ee2d8b0c10e377/step2_RAW2DIGI_L1Reco_RECO_100_1_Ngh.root')
#fileNames = cms.untracked.vstring('root://xrootd.unl.edu//store/user/ppigard/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6_Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1_MODMINIAOD_All/ff0f7a48dd492a8600ee2d8b0c10e377/step2_RAW2DIGI_L1Reco_RECO_103_1_WWl.root')

fileNames = cms.untracked.vstring('file:/home/llr/cms/pigard/CMSSW_7_4_1_patch1/src/Analyzer/Ntuplizer/python/0033A97B-8707-E511-9D3B-008CFA1980B8.root')
#fileNames = cms.untracked.vstring('file:/home/llr/cms/pigard/GluGluHToZZTo4L_M125_13TeV_powheg_JHUgen_pythia8_MINIAODSIM.root')

)
 
#process.GlobalTag.globaltag = ' MCRUN2_74_V8::All'
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag

#process.GlobalTag.globaltag = 'PHYS14_25_V1::All'
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
        'LOG', #output//' + fileNameForSample + '.log',
        'critical'
    ),
    LOG = cms.untracked.PSet(
        threshold  = cms.untracked.string('WARNING'), # DEBUG 
        filename  = cms.untracked.string(fileNameForSample  + '.log')
    ),
    debugModules = cms.untracked.vstring('*'), # *
       
    statistics     = cms.untracked.vstring('STAT'#'output//' + fileNameForSample  + '_stats.log',
    ),
        STAT = cms.untracked.PSet(
            threshold = cms.untracked.string('WARNING'),
            filename  = cms.untracked.string(fileNameForSample  + '_stats.log')
    )
)

process.p = cms.Path(process.ntuple)#*process.dump)
