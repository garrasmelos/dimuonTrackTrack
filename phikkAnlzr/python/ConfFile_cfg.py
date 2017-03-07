import FWCore.ParameterSet.Config as cms

process = cms.Process("Dimuonkk")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data','')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        '/store/group/phys_bphys/asanchez/MuOnia/bphskim/170215_050232/0000/BPHSkim_1.root'
    )
)

#process.dimuonkk = cms.EDAnalyzer('phikkAnlzr'
#)
process.load('dimuonTrackTrack.phikkAnlzr.CfiFile_cfi')
process.dimuonkk.HLTLastFilters = cms.vstring (
      'hltDisplacedmumuFilterDimuon0PhiBarrel',
      'hltDiMuonGlb16Trk0DzFiltered0p2Phi'
)
process.TFileService = cms.Service("TFileService",
                   fileName = cms.string("hTest.root")
                   )

process.p = cms.Path(process.dimuonkk)
