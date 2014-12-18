import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing ('standard')
options.register ('partonStatusList',
                  [],
                  VarParsing.VarParsing.multiplicity.list, # singleton or list
                  VarParsing.VarParsing.varType.int,         # string, int, or float
                  "Status to consider as parton (default: status 3 only)")
options.register ('writeAllGenParticles',
                  False,
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.bool,         # string, int, or float
                  "write all gen particles on parton collection")
options.register ('GlobalTag',
                  'START53_LV4::All',
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.string,         # string, int, or float
                  "Global tag (if unspecified, use START53_LV4)")
options.maxEvents=-1
options.partonStatusList.append(3)
options.files='file:bla.root'
options.parseArguments()
print 'PARTONS ARE GENPARTICLES WITH STATUS = '+str(options.partonStatusList)

process = cms.Process("GenLevelAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string(options.GlobalTag)

#process.load("Configuration.StandardSequences.GeometryDB_cff")
#process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load('Configuration.StandardSequences.Reconstruction_cff')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.files),
#    duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)

process.analyze = cms.EDAnalyzer("GenLevelAnalyzer",
                                 genParticlesTag = cms.InputTag("genParticles"),
                                 genJetsTag = cms.InputTag("ak5GenJets"),
                                 minPtPhotons = cms.double(25),
                                 maxEtaPhotons = cms.double(2.5),
                                 minPtJets = cms.double(25),
                                 maxEtaJets = cms.double(4.7),
                                 partonStatusList = cms.vint32(options.partonStatusList),
                                 writeAllGenParticles = cms.bool(options.writeAllGenParticles),
                                 OutputFile = cms.string('outfile.root')
                                 )

process.p = cms.Path(process.analyze)
process.schedule = cms.Schedule(process.p)
