import FWCore.ParameterSet.Config as cms

generator = cms.EDFilter("Pythia8HadronizerFilter",
    maxEventsToPrint = cms.untracked.int32(1),
    pythiaPylistVerbosity = cms.untracked.int32(1),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(7000.),
    jetMatching = cms.untracked.PSet(
      scheme = cms.string("MadgraphPy8Internal"),
    ),
    PythiaParameters = cms.PSet(
        processParameters = cms.vstring(
            'Main:timesAllowErrors = 10000',
            'Check:epTolErr = 0.01',
            'SLHA:keepSM = on',
            'SLHA:minMassSM = 1000.',
            'ParticleDecays:limitTau0 = on',
            'ParticleDecays:tau0Max = 10',
            'ParticleDecays:allowPhotonRadiation = on',
            'Tune:pp 14',
            'Tune:ee 7',
            'MultipartonInteractions:pT0Ref=2.4024',
            'MultipartonInteractions:ecmPow=0.25208',
            'MultipartonInteractions:expPow=1.6',
        'SpaceShower:pTmaxMatch = 1',
        'SpaceShower:pTmaxFudge = 1',
        'SpaceShower:MEcorrections = off',
        'TimeShower:pTmaxMatch = 1',
        'TimeShower:pTmaxFudge = 1',
        'TimeShower:MEcorrections = off',
        'TimeShower:globalRecoil = on',
        'TimeShower:limitPTmaxGlobal = on',
        'TimeShower:nMaxGlobalRecoil = 1',
        'TimeShower:globalRecoilMode = 2',
        'TimeShower:nMaxGlobalBranch = 1',
        'JetMatching:setMad = off',
        'JetMatching:scheme = 1',
        'JetMatching:merge = on',
        'JetMatching:jetAlgorithm = 2',
        'JetMatching:etaJetMax = 999.',
        'JetMatching:coneRadius = 1.',
        'JetMatching:slowJetPower = 1',
        'JetMatching:qCut = 30.', #this is the actual merging scale
        'JetMatching:doFxFx = on',
        'JetMatching:qCutME = 10.',#this must match the ptj cut in the lhe generation step
        'JetMatching:nQmatch = 4', #4 corresponds to 4-flavour scheme (no matching of b-quarks), 5 for 5-flavour scheme
        'JetMatching:nPartonsNow = 0', #number of partons in born matrix element
        'TimeShower:nPartonsInBorn = 0', #number of partons in born matrix element
        'JetMatching:nJetMax = 1', #number of partons in born matrix element for highest multiplicity
        ),
        parameterSets = cms.vstring('processParameters')
    )
)
