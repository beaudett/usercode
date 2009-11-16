import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")
process.load("DQMServices.Core.DQM_cfg")
process.DQM.collectorHost = ''

# Source : general definition
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(),                            
                            secondaryFileNames = cms.untracked.vstring(),
                            noEventSort = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
                            )

# Update input files
#print dbs_discovery.search()
#process.source.fileNames=['/store/relval/CMSSW_3_3_0/RelValSingleElectronPt35/GEN-SIM-RECO/MC_31X_V9-v1/0009/98C1C609-75B7-DE11-941B-001D09F295FB.root','/store/relval/CMSSW_3_3_0/RelValSingleElectronPt35/GEN-SIM-RECO/MC_31X_V9-v1/0008/100559BC-EDB6-DE11-BA68-000423D6006E.root']
#process.source.fileNames=['file:98C1C609-75B7-DE11-941B-001D09F295FB.root']
process.source.fileNames = [
#        '/store/relval/CMSSW_3_3_0/RelValQCD_Pt_80_120/GEN-SIM-RECO/MC_31X_V9-v1/0009/880AA097-75B7-DE11-B848-001D09F23C73.root',
        '/store/relval/CMSSW_3_3_0/RelValQCD_Pt_80_120/GEN-SIM-RECO/MC_31X_V9-v1/0008/CC9C9FAC-86B6-DE11-8A1B-001D09F24FEC.root', 
        '/store/relval/CMSSW_3_3_0/RelValQCD_Pt_80_120/GEN-SIM-RECO/MC_31X_V9-v1/0008/A44075CB-7FB6-DE11-8010-000423D98868.root', 
        '/store/relval/CMSSW_3_3_0/RelValQCD_Pt_80_120/GEN-SIM-RECO/MC_31X_V9-v1/0008/6E433585-84B6-DE11-A3E1-001D09F2437B.root', 
        '/store/relval/CMSSW_3_3_0/RelValQCD_Pt_80_120/GEN-SIM-RECO/MC_31X_V9-v1/0008/285A7AC1-82B6-DE11-BA29-001D09F2525D.root', 
        '/store/relval/CMSSW_3_3_0/RelValQCD_Pt_80_120/GEN-SIM-RECO/MC_31X_V9-v1/0008/12755919-88B6-DE11-8FF2-000423D996C8.root', 
        '/store/relval/CMSSW_3_3_0/RelValQCD_Pt_80_120/GEN-SIM-RECO/MC_31X_V9-v1/0008/04B93222-89B6-DE11-93E9-001D09F29524.root']

process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(-1)
        )

# Calo geometry service model
process.load("Configuration.StandardSequences.Geometry_cff")

from FastSimulation.Calorimetry.Calorimetry_cff import *
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.pfAllElectrons = cms.EDProducer("PdgIdPFCandidateSelector",
    pdgId = cms.vint32(11, -11),
    src = cms.InputTag("pfNoPileUp")
)


process.gensource = cms.EDProducer("GenParticlePruner",
    src = cms.InputTag("genParticles"),
    select = cms.vstring('drop *', 
        'keep pdgId = 211', 
        'keep pdgId = -211')
)


process.pfPileUp = cms.EDProducer("PFPileUp",
    PFCandidates = cms.InputTag("particleFlow"),
    verbose = cms.untracked.bool(False),
    Vertices = cms.InputTag("offlinePrimaryVerticesWithBS")
)


process.pfNoPileUp = cms.EDProducer("TPPileUpPFCandidatesOnPFCandidates",
    bottomCollection = cms.InputTag("particleFlow"),
    topCollection = cms.InputTag("pfPileUp"),
    name = cms.untracked.string('pileUpOnPFCandidates'),
    verbose = cms.untracked.bool(False)
)
process.pfNoPileUpSequence = cms.Sequence(process.pfPileUp+process.pfNoPileUp)


process.trackExtra = cms.EDAnalyzer("TrackExtra",FamosCalorimetryBlock,
                           GSFTracks=cms.InputTag("electronGsfTracks"),
                           BarrelRecHits=cms.InputTag("ecalRecHit:EcalRecHitsEB"),
                           EndcapRecHits=cms.InputTag("ecalRecHit:EcalRecHitsEE"),
                           PFCandidate = cms.InputTag("particleFlow:electrons"),
                           InputTruthLabel = cms.InputTag("gensource"),
                           DeltaRMatch = cms.double(0.05),
                           MatchWithMC = cms.bool(True),
                                    
                           TestParticleFilter = cms.PSet(
    # Particles with |eta| > etaMax (momentum direction at primary vertex)
    # are not simulated
    etaMax = cms.double(5.0),
    # Charged particles with pT < pTMin (GeV/c) are not simulated
    pTMin = cms.double(0.0),
    # Particles with energy smaller than EMin (GeV) are not simulated
    EMin = cms.double(0.0),
    # Protons with energy in excess of this value (GeV) will kept no matter what
    EProton = cms.double(99999.0)
    ))

process.p =cms.Path(process.pfNoPileUpSequence+process.pfAllElectrons+
                    process.gensource+
                    process.trackExtra)


#process.out = cms.OutputModule("PoolOutputModule",
#                               outputCommands = cms.untracked.vstring('keep *'),
#                               outputFile = cms.string(os.environ['TEST_OUTPUT_FILE'])
#                               )
#process.outpath = cms.EndPath(process.out)

process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.MessageLogger.cerr.FwkReport.reportEvery = 100



