// -*- C++ -*-
//
// Package:    dimuonTrackTrack/phikkAnlzr
// Class:      phikkAnlzr
// 
/**\class phikkAnlzr phikkAnlzr.cc dimuonTrackTrack/phikkAnlzr/plugins/phikkAnlzr.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Gabriel Ramirez Sanchez
//         Created:  Tue, 21 Feb 2017 20:10:05 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateClosestToPoint.h"
#include "FWCore/Utilities/interface/GCC11Compatibility.h"


#include <vector>
#include "TMath.h"
#include "Math/VectorUtil.h"
#include "TVector3.h"
#include "TTree.h"
#include "TLorentzVector.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class phikkAnlzr : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit phikkAnlzr(const edm::ParameterSet&);
      ~phikkAnlzr();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
   	bool isInPV(const edm::Event& , const edm::EventSetup&, reco::Track ,reco::Vertex );
   	UInt_t getTriggerBits(const edm::Event &, std::vector<std::string>);
   	UInt_t isTriggerMatched(const pat::CompositeCandidate *);
   	
   	
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      TTree* oniaTree;
      TTree* dikaonTree;
      UInt_t nPVs=0;
      UInt_t tksInPV=0;
      UInt_t isMatched=0;
      UInt_t trig=0;
      double kaonMass=0.493677;
      TLorentzVector dimuon_p4;
      TLorentzVector muonP_p4;
      TLorentzVector muonN_p4;
      
      TLorentzVector dikaon_p4;
      TLorentzVector kaonP_p4;
      TLorentzVector kaonN_p4;   
      
		UInt_t run=0;
		UInt_t nevent=0;
      Double_t vProb_diMu=0.;
      
      const reco::TrackBase::TrackQuality trackQuality_;
      edm::EDGetTokenT<pat::CompositeCandidateCollection> oniaToken_;
      edm::EDGetTokenT<reco::TrackCollection> trackToken_;
      edm::EDGetTokenT<pat::MuonCollection> muonToken_;
      edm::EDGetTokenT<reco::VertexCollection> pvsToken_;
      edm::EDGetTokenT<edm::TriggerResults> triggerResults_;
      
      std::vector<std::string> HLTLastFilters;
      std::vector<std::string> FilterNames_;
      std::vector<double> diMuMassRange_;
};

//
// constants, enums and typedefs
//
//Test from mac
//
// static data member definitions
//

//
// constructors and destructor
//
phikkAnlzr::phikkAnlzr(const edm::ParameterSet& iConfig)
: 	trackQuality_(reco::TrackBase::qualityByName(iConfig.getParameter<std::string>("TrackQuality"))),
	oniaToken_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("oniaLabel"))),
	trackToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("trackLabel"))),
	muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter< edm::InputTag>("muonLabel"))),
	pvsToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("pvsLabel"))),
	triggerResults_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
	HLTLastFilters(iConfig.getParameter<std::vector<std::string>>("HLTLastFilters")),
	FilterNames_(iConfig.getParameter<std::vector<std::string>>("filterNames")),
	diMuMassRange_(iConfig.getParameter<std::vector<double>>("diMuMassRange"))
	
{
   //now do what ever initialization is needed
   //usesResource("TFileService");

}


phikkAnlzr::~phikkAnlzr()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

bool phikkAnlzr::isInPV(const edm::Event& iEvent, const edm::EventSetup& iSetup, reco::Track trk,reco::Vertex thepv)
{
	using namespace reco;
	using namespace edm;
	using namespace std;
	edm::ESHandle<TransientTrackBuilder> theB;
	iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", theB);
	Handle<reco::VertexCollection> pvts;
	iEvent.getByToken(pvsToken_, pvts);
		
	TransientTrack t_trk = (*theB).build(trk);
	bool accepted=false;
	Vertex vtx_min;
	for(reco::VertexCollection::const_iterator ipv = pvts->begin(); ipv != pvts->end(); ipv++)
	{
		double dz_min=10000.;
		GlobalPoint vert(ipv->x(), ipv->y(), ipv->z());
		TrajectoryStateClosestToPoint  traj = t_trk.trajectoryStateClosestToPoint(vert);
		double dz = traj.perigeeParameters().longitudinalImpactParameter();
 				
		if(ipv == pvts->begin()) 
		{
			dz_min = TMath::Abs(dz);
			vtx_min = Vertex(*ipv);
		}
		else 
		{
			if(dz_min > TMath::Abs(dz))
			{
				vtx_min = Vertex(*ipv); 
			}
		}
	}
	if (vtx_min.x() == thepv.x() && vtx_min.y() == thepv.y() && vtx_min.z() == thepv.z()) accepted=true;
	return accepted;
}


UInt_t phikkAnlzr::getTriggerBits(const edm::Event& iEvent, std::vector<std::string> TestFilterNames_)
{
	UInt_t trigger=0;
	edm::Handle<edm::TriggerResults> triggerResults_handle;
	iEvent.getByToken(triggerResults_,triggerResults_handle);
	if(triggerResults_handle.isValid())
	{
		const edm::TriggerNames & triggerNames = iEvent.triggerNames(*triggerResults_handle);
		for(unsigned i=0 ; i< TestFilterNames_.size();i++)
		{
			for(int ver=1 ; ver < 9 ; ver++)
			{
				std::stringstream ss;
				ss << TestFilterNames_[i] << "_v" << ver;
				unsigned int bit = triggerNames.triggerIndex(edm::InputTag(ss.str()).label().c_str());
				if(bit < triggerResults_handle->size() && triggerResults_handle->accept(bit) && !triggerResults_handle->error(bit))
				{
					//trigger += 1 << bit;
					trigger += 1 << (i*9 + ver);
					break;
				}
			}
		}
	} else std::cout << "miniAnalyzer::getTriggerBits: **** No triggers found ****" << std::endl;
	return trigger;
}


UInt_t phikkAnlzr::isTriggerMatched(const pat::CompositeCandidate *diMuon_cand) 
{
	const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(diMuon_cand->daughter("muon1"));
	const pat::Muon* muon2 = dynamic_cast<const pat::Muon*>(diMuon_cand->daughter("muon2"));
	UInt_t matched = 0;  // if no list is given, is not matched 
	// if matched a given trigger, set the bit, in the same order as listed
	for (unsigned int iTr = 0; iTr<HLTLastFilters.size(); iTr++ ) 
	{
		const pat::TriggerObjectStandAloneCollection mu1HLTMatches = muon1->triggerObjectMatchesByFilter(HLTLastFilters[iTr]);
		const pat::TriggerObjectStandAloneCollection mu2HLTMatches = muon2->triggerObjectMatchesByFilter(HLTLastFilters[iTr]);
		if (mu1HLTMatches.size() > 0 && mu2HLTMatches.size() > 0) matched += (1<<iTr); 
	}
	return matched;
}


// ------------ method called for each event  ------------
void
phikkAnlzr::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	using namespace std;
	using namespace reco;
	
	Handle<pat::CompositeCandidateCollection> onias;
	iEvent.getByToken(oniaToken_, onias);
	
	Handle<reco::TrackCollection> tracks;
	iEvent.getByToken(trackToken_,tracks);
	
	Handle<reco::VertexCollection> pvs;
	iEvent.getByToken(pvsToken_, pvs);
	
	dimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
	muonP_p4.SetPtEtaPhiM(0.,0.,0.,0.);
	muonN_p4.SetPtEtaPhiM(0.,0.,0.,0.);
   
	nPVs = pvs->size();
   
   //cout << "First track charge: " << tracks->at(0).charge() << endl;
	run = iEvent.id().run();
	nevent = iEvent.id().event();
	trig       = getTriggerBits(iEvent,FilterNames_);
	//cout << "trig: "<< trig << endl;
	if(trig)
	{
		Vertex thePV;
		UInt_t diMutight = 0;
		double diMuMinMass = diMuMassRange_[0];
		double diMuMaxMass = diMuMassRange_[1];
		//cout << "Number of onias: " << onias->size() << endl;
		if(onias->size()>0)
		{	
			for(pat::CompositeCandidateCollection::const_iterator idimuon = onias->begin(); idimuon != onias->end(); idimuon++)
			{
				const pat::Muon* mu1 = dynamic_cast<const pat::Muon*>(idimuon->daughter("muon1"));
				const pat::Muon* mu2 = dynamic_cast<const pat::Muon*>(idimuon->daughter("muon2"));
 	    		
   			
				vProb_diMu = idimuon->userFloat("vProb");
				if(idimuon->mass()>diMuMinMass && idimuon->mass()<diMuMaxMass)
				{
					isMatched = isTriggerMatched(&*idimuon);
					dimuon_p4.SetPtEtaPhiM(idimuon->pt(),idimuon->eta(),idimuon->phi(),idimuon->mass());
					reco::Candidate::LorentzVector vP = idimuon->daughter("muon1")->p4();
					reco::Candidate::LorentzVector vM = idimuon->daughter("muon2")->p4();
					if ( idimuon->daughter("muon1")->charge() < 0 )
					{
						vP = idimuon->daughter("muon2")->p4();
						vM = idimuon->daughter("muon1")->p4();
					}
					muonP_p4.SetPtEtaPhiM(vP.pt(),vP.eta(),vP.phi(),vP.mass());
					muonN_p4.SetPtEtaPhiM(vM.pt(),vM.eta(),vM.phi(),vM.mass());
					if(pvs->size()>0)
					{
						if(idimuon->userData<reco::Vertex>("commonVertex"))
						{
							const reco::Vertex* dimuvertex = idimuon->userData<reco::Vertex>("commonVertex");
							//Taking the first vertex of the list as the PV.
							//Choosing the closest PV to the dimuon vertex (in z).
							for(reco::VertexCollection::const_iterator ipv = pvs->begin(); ipv != pvs->end(); ipv++)
							{                              
								if(ipv!=pvs->begin())
								{
									if(TMath::Abs(dimuvertex->z()-thePV.z()) > TMath::Abs(dimuvertex->z() - ipv->z())) thePV = Vertex(*ipv);//thePV = ipv;
								}	
								else
								{
									thePV = Vertex(*ipv);
								}                   
							}
							diMutight = 0; 
							if (mu1->isTightMuon(thePV) && mu2->isTightMuon(thePV)) diMutight=1;
							if(!diMutight) break;
							oniaTree->Fill();
							break;
						}
					}
				}
			}   		
		}else
		{
			cout << "There was not dimuon in this event." <<endl;
		}
		vector<reco::Track> tracksPInPV, tracksNInPV;
		cout << "tracks size " << tracks->size() << " " << diMutight << endl;
		if(tracks->size()>0 && diMutight)
		{	   		   	
			for(auto& itk : *tracks)
			{
				bool inPV = isInPV(iEvent,iSetup,itk, thePV);
				if(inPV)
				{  
					if(itk.charge() > 0.) tracksPInPV.push_back(itk);
					else tracksNInPV.push_back(itk);
				}
			}
			cout << "%% Here 1 %%" << endl;   		
			for(unsigned int i =0; i < tracksPInPV.size(); i++)
			{   			
				if(tracksPInPV[i].quality(trackQuality_))
				{
					for(unsigned int j = 0 ; j < tracksNInPV.size(); j++)
					{
						if(tracksNInPV[j].quality(trackQuality_))
						{
							kaonP_p4.SetPtEtaPhiM(tracksPInPV[i].pt(),tracksPInPV[i].eta(),tracksPInPV[i].phi(),0.493677);
							kaonN_p4.SetPtEtaPhiM(tracksNInPV[j].pt(),tracksNInPV[j].eta(),tracksNInPV[j].phi(),0.493677);
							dikaon_p4 = kaonP_p4 + kaonN_p4;	
							//phikk_q = tracksPInPV[i].charge() + tracksNInPV[j].charge();
							if(dikaon_p4.M()>0.5 && dikaon_p4.M() <1.15) 
							{
								dikaonTree->Fill();
								cout << "%% Here 2 %%" << endl;
							}
						}
					}
				}
			}
		}
	}
}
	// ------------ method called once each job just before starting event loop  ------------
	
void 
phikkAnlzr::beginJob()
{
	edm::Service<TFileService> fs;
	oniaTree = fs->make<TTree>("oniaTree","Tree of phi and tracks");
	oniaTree->Branch("nPVs",				&nPVs,						"nPVs/i");
	oniaTree->Branch("tksInPV",			&tksInPV,					"tksInPV/i");
	oniaTree->Branch("isMatched",			&isMatched,					"isMatched/i");
	oniaTree->Branch("trig",				&trig,						"trig/i");
	oniaTree->Branch("vProb_diMu",		&vProb_diMu,				"vProb_diMu/d");
	oniaTree->Branch("dimuon_p4",			"TLorentzVector",			&dimuon_p4);
	oniaTree->Branch("muonP_p4",			"TLorentzVector",			&muonP_p4);
	oniaTree->Branch("muonN_p4",			"TLorentzVector",			&muonN_p4);

	dikaonTree = fs->make<TTree>("dikaonTree","Tree of phi to kk");
	dikaonTree->Branch("dikaon_p4",		"TLorentzVector",			&dikaon_p4);
	dikaonTree->Branch("kaonP_p4",		"TLorentzVector",			&kaonP_p4);
	dikaonTree->Branch("kaonN_p4",		"TLorentzVector",			&kaonN_p4);
	dikaonTree->Branch("run",				&run,							"run/i");
	dikaonTree->Branch("nevent",			&nevent,						"nevent/i");
	return;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
phikkAnlzr::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
phikkAnlzr::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(phikkAnlzr);
