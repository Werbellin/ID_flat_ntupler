from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.requestName = 'Signal-run2b'
config.General.transferOutputs = True

config.section_('JobType')
config.JobType.psetName = 'ntuplizer_cfg.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['ntuple.root', 'ntuple.log', 'ntuple_stats.log']

config.section_('Data')
config.Data.inputDataset = '/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6/ppigard-GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6_Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1_MODMINIAOD_All-ff0f7a48dd492a8600ee2d8b0c10e377/USER'
#config.Data.inputDataset = '/DYJetsToLL_M-50_13TeV-madgraph-pythia8/ppigard-DYJetsToLL_M-50_13TeV-madgraph-pythia8_Phys14DR-PU20bx25_PHYS14_25_V1-v1_MODAODSIM_All-ff0f7a48dd492a8600ee2d8b0c10e377/USER'
config.Data.inputDBS = 'phys03'
config.Data.outLFNDirBase = '/store/user/ppigard/eID/MC/HToZZTo4L/run2b'
config.Data.totalUnits = -1
config.Data.unitsPerJob = 1
config.Data.splitting = 'FileBased'
config.Data.publication = False

config.section_('User')

config.section_('Site')
#config.Site.whitelist = ['T2_FR_GRIF_LLR']
config.Site.storageSite = 'T2_FR_GRIF_LLR'
