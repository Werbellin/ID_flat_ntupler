from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.requestName = 'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v3_AODSIM_v4'
config.General.transferOutputs = True

config.section_('JobType')
config.JobType.psetName = 'ntuplizer_cfg.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['ntuple.root', 'ntuple.log', 'ntuple_stats.log']

config.section_('Data')

config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v3/AODSIM'


config.Data.inputDBS = 'global'
config.Data.outLFNDirBase = '/store/user/ppigard/eID/MC/Spring15/25ns/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v3_AODSIM/v4'
config.Data.totalUnits = -1
config.Data.unitsPerJob = 10
config.Data.splitting = 'FileBased'
config.Data.publication = False

config.section_('User')

config.section_('Site')
#config.Site.whitelist = ['T2_FR_GRIF_LLR']
config.Site.storageSite = 'T2_FR_GRIF_LLR'
