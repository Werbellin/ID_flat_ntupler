import FWCore.ParameterSet.Config as cms

#process.load("RecoTracker.Configuration.python.RecoTracker_cff")
process = cms.Process("Demo")

process.dump=cms.EDAnalyzer('EventContentAnalyzer')
process.load("Configuration.StandardSequences.MagneticField_cff")

#process.load("Configuration.StandardSequences.Geometry_cff")
process.load('Configuration.Geometry.GeometryIdeal_cff')

process.load('Configuration.EventContent.EventContent_cff')

process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")




process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )


process.load("RecoTracker.TrackProducer.KfTrackFromGsfTrack_cfi")
process.KfTrackFromGsfTrack.src = cms.InputTag("electronGsfTracks")
theCKFTrajectoryInput = cms.string('KfTrackFromGsfTrack')
process.KFFittingSmootherWithOutliersRejectionAndRK.EstimateCut = 10000

#process.load("RecoTracker.TrackProducer.TrackProducerOwnLast5_cfi")
#process.TrackProducerOwnLast5.GSFTrajectoryInput  = cms.string('GsfTrackRefitter')

process.load("RecoTracker.TrackProducer.TrackRefittersMod_cff")
process.TrackRefitter.src = cms.InputTag("generalTracks")

process.load("RecoTracker.TrackProducer.GsfTrackRefitterMod_cff")
process.GsfTrackRefitter.src = cms.InputTag("electronGsfTracks")
process.GsfTrackRefitter.TrajectoryInEvent = cms.bool(True)

theGSFTrajectoryInput = cms.string('GsfTrackRefitter')



therunGsfRefitter      = cms.bool(True)
therunKfWithGsfRefitter = cms.bool(True)

process.ntuple = cms.EDAnalyzer('Ntuplizer',
                                   EleTag      = cms.InputTag('gedGsfElectrons'),
                                   VerticesTag = cms.InputTag('offlinePrimaryVertices'),
                                   HLTTag          = cms.InputTag('TriggerResults','','HLT'),
                                   isMC = cms.bool(True),
                                   ispythia6 = cms.bool(True),
                                   runGsfRefitter      = therunGsfRefitter,
                                   GSFTrajectoryInput  = theGSFTrajectoryInput,
                                   runKfWithGsfRefitter = therunKfWithGsfRefitter,
                                   CKFTrajectoryInput  = theCKFTrajectoryInput,

                                   MVAId  = cms.VInputTag()#"mvaNonTrigV025nsPHYS14", "mvaNonTrigV025nsPHYS14FIX")
)
fileNameForSample = 'ntuple'


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",
#    fileNames = cms.untracked.vstring(
        #'root://xrootd.unl.edu//store/user/ppigard/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6_Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1_MODMINIAOD_All/ff0f7a48dd492a8600ee2d8b0c10e377/step2_RAW2DIGI_L1Reco_RECO_100_1_Ngh.root')
#fileNames = cms.untracked.vstring('root://xrootd.unl.edu//store/user/ppigard/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6_Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1_MODMINIAOD_All/ff0f7a48dd492a8600ee2d8b0c10e377/step2_RAW2DIGI_L1Reco_RECO_103_1_WWl.root')
fileNames = cms.untracked.vstring('file:/home/llr/cms/pigard/CMSSW_7_2_3/res/ALCARECO-test/02405CE0-8E6B-E411-A37C-20CF3019DF03.root')
)
 
#if process.demo.vCMSSW == '7' :
process.GlobalTag.globaltag = 'PHYS14_25_V1::All'
#else :
 #   process.GlobalTag.globaltag = 'START53_V29A::All'


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

process.p = cms.Path(process.MeasurementTrackerEvent * process.KfTrackFromGsfTrack * process.GsfTrackRefitter * process.ntuple)#*process.dump)
#

#if process.demo.runGsfRefitter == cms.bool(False) and process.demo.runKfRefitter == cms.bool(False) and process.demo.runKfWithGsfRefitter == cms.bool(True):
#    print 'Running custom refitter'
#    process.p = cms.Path(process.TrackRefitterWithGsfHits * process.demo)
#if process.demo.runGsfRefitter == cms.bool(False) and process.demo.runKfRefitter == cms.bool(False) and process.demo.runKfWithGsfRefitter == cms.bool(False):
#    process.p = cms.Path(process.demo)


#if process.demo.runGsfRefitter == cms.bool(True) and process.demo.runKfRefitter == cms.bool(True) and process.demo.runKfWithGsfRefitter == cms.bool(True):
#    print 'This is gonna take some time' 
#    process.p = cms.Path(process.TrackRefitterWithGsfHits * process.TrackRefitter * process.GsfTrackRefitter *  process.demo)#*process.dump)
#* process.TrackProducerOwnLast5 
#process.p = cms.Path(process.GsfTrackRefitter * process.demo)#*process.dump)
#process.p = cms.Path(process.MeasurementTrackerEvent * process.TrackRefitter * process.GsfTrackRefitter * process.demo * process.ntuple)#*process.dump)

#if process.demo.runKfWithGsfRefitter == cms.bool(False) :
#    if process.demo.runGsfRefitter == cms.bool(False) and process.demo.runKfRefitter == cms.bool(False) :
#        process.p = cms.Path(process.demo)#*process.dump)
#    if process.demo.runGsfRefitter == cms.bool(True) and process.demo.runKfRefitter == cms.bool(True) :
#        process.p = cms.Path(process.dump*process.TrackRefitter * process.GsfTrackRefitter * process.demo)#*process.dump)
#    if process.demo.runGsfRefitter == cms.bool(True) and process.demo.runKfRefitter == cms.bool(False) :
#        process.p = cms.Path(process.GsfTrackRefitter * process.demo)#*process.dump)
#    if process.demo.runGsfRefitter == cms.bool(False) and process.demo.runKfRefitter == cms.bool(True) :
#        process.p = cms.Path(process.TrackRefitter * process.demo)#*process.dump)
