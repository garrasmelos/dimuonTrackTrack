import FWCore.ParameterSet.Config as cms

dimuonkk = cms.EDAnalyzer('phikkAnlzr',
								oniaLabel = cms.InputTag("onia2MuMuPAT"),
								muonLabel = cms.InputTag("replaceme"),
								trackLabel = cms.InputTag("oniaSelectedTracks"),
								pvsLabel = cms.InputTag("offlinePrimaryVertices"),
								triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                     	filterNames = cms.vstring(
                           'HLT_Dimuon0_Phi_Barrel',
                           'HLT_Mu16_TkMu0_dEta18_Phi'
                     	),
                     	diMuMassRange = cms.vdouble(0.85,1.2),
                     	TrackQuality = cms.string("highPurity")
                     	
                     	
)
