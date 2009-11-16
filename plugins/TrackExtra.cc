#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "UserCode/FBeaudette/plugins/TrackExtra.h"
#include "FastSimulation/CaloGeometryTools/interface/CaloGeometryHelper.h"
#include "FastSimulation/CaloHitMakers/interface/EcalHitMaker.h"
#include "FastSimulation/CaloGeometryTools/interface/CaloGeometryHelper.h"
#include "FastSimulation/ShowerDevelopment/interface/EMECALShowerParametrization.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"


#include "FastSimulation/Event/interface/FSimEvent.h"


TrackExtra::TrackExtra(const edm::ParameterSet& iConfig)
{
  inputTagGSFTracks_  = iConfig.getParameter<edm::InputTag>("GSFTracks");  
  inputTagBarrelRecHits_ = iConfig.getParameter<edm::InputTag>("BarrelRecHits");  
  inputTagEndcapRecHits_ = iConfig.getParameter<edm::InputTag>("EndcapRecHits");  
  inputTagPFCandidates_ = iConfig.getParameter<edm::InputTag>("PFCandidate");
  inputTagTruth_ = iConfig.getParameter<edm::InputTag>("InputTruthLabel");
  deltarMatch_ = iConfig.getParameter<double>("DeltaRMatch");
  matchWithMC_ = iConfig.getParameter<bool>("MatchWithMC");
  myGeometry = new CaloGeometryHelper(iConfig.getParameter<edm::ParameterSet>("Calorimetry"));
  mySimEvent = new FSimEvent(iConfig.getParameter<edm::ParameterSet>( "TestParticleFilter" ));

}

TrackExtra::~TrackExtra()
{
  dbe->save("track.root");
}

void 
TrackExtra::beginRun(edm::Run const&, edm::EventSetup const& iSetup)
{
  dbe = edm::Service<DQMStore>().operator->();
  h0 = dbe->book2D("h0","Profile all ",100,0.,35.,100,0,1.5);
  h2 = dbe->book2D("h2","Profile barrel ",100,0.,35.,100,0,1.5);
  h22 = dbe->book2D("h22","Profile barrel N=2",100,0.,35.,100,0,1.5);
  h23 = dbe->book2D("h23","Profile barrel N=3",100,0.,35.,100,0,1.5);
  h4 = dbe->book2D("h4","Profile endcap ",100,0.,35.,100,0,1.5);
  h42 = dbe->book2D("h42","Profile endcap ",100,0.,35.,100,0,1.5);
  h43 = dbe->book2D("h43","Profile endcap ",100,0.,35.,100,0,1.5);

  h100 = dbe->book2D("h100","Profile all ",100,0.,35.,100,0,1.5);
  h102 = dbe->book2D("h102","Profile barrel ",100,0.,35.,100,0,1.5);
  h122 = dbe->book2D("h122","Profile barrel N=2",100,0.,35.,100,0,1.5);
  h123 = dbe->book2D("h123","Profile barrel N=3",100,0.,35.,100,0,1.5);
  h104 = dbe->book2D("h104","Profile endcap ",100,0.,35.,100,0,1.5);
  h142 = dbe->book2D("h142","Profile endcap ",100,0.,35.,100,0,1.5);
  h143 = dbe->book2D("h143","Profile endcap ",100,0.,35.,100,0,1.5);
  
  h200 = dbe->book2D("h200","expected profile all ",100,0.,35.,100,0,2.5);
  h202 = dbe->book2D("h202","expected profile all N=2",100,0.,35.,100,0,2.5);
  h203 = dbe->book2D("h203","expected profile all N=3",100,0.,35.,100,0,2.5);

  h502 = dbe->book2D("h502","measured/expected profile all N=2",100,0.,35.,100,0,2.5);
  h503 = dbe->book2D("h503","measured/expected profile all N=3",100,0.,35.,100,0,2.5);
  h602 = dbe->book2D("h602","measured/expected profile all N=2",100,0.,35.,100,0,2.5);
  h603 = dbe->book2D("h603","measured/expected profile all N=3",100,0.,35.,100,0,2.5);

  h10 = dbe->book1D("h10","N crystals",10,-0.5,9.5);
  h12 = dbe->book1D("h12","N crystals barrel",10,-0.5,9.5);
  h14 = dbe->book1D("h14","N crystals endcaps",10,-0.5,9.5);
  h700 = dbe->book1D("h700","Delta R",100,0.,1.);
  h30 = dbe->book1D("h30","eta E=0",140,-3.5,3.5);
  h32 = dbe->book2D("h32","ECAL entrance",350,-350,350,150,0,150);
  edm::ESHandle<CaloTopology> theCaloTopology;
  iSetup.get<CaloTopologyRecord>().get(theCaloTopology);     
  
  edm::ESHandle<CaloGeometry> pG;
  iSetup.get<CaloGeometryRecord>().get(pG);     
  
  
  // Setup the tools
  double bField000 = 4.;
  myGeometry->setupGeometry(*pG);
  myGeometry->setupTopology(*theCaloTopology);
  myGeometry->initialize(bField000);
  
  // init Particle data table (from Pythia)
  edm::ESHandle < HepPDT::ParticleDataTable > pdt;
  iSetup.getData(pdt);
  mySimEvent->initializePdt(&(*pdt));
  
}

void TrackExtra::analyze(const edm::Event & iEvent,const edm::EventSetup & c)
{
  edm::Handle<reco::GsfTrackCollection> gsfTracksH;
  bool found=iEvent.getByLabel(inputTagGSFTracks_,gsfTracksH);
  if(!found ) {
    std::ostringstream  err;
    err<<" cannot get GsfTracks: "
       <<inputTagGSFTracks_<<std::endl;
    edm::LogError("TrackExtra")<<err.str();
    throw cms::Exception( "MissingProduct", err.str());
  }  

  edm::Handle<reco::PFCandidateCollection> pfCandidatesH;
  iEvent.getByLabel(inputTagPFCandidates_, pfCandidatesH);


  edm::Handle<EBRecHitCollection> hrechit_EB_col;
  iEvent.getByLabel(inputTagBarrelRecHits_,hrechit_EB_col);
  edm::Handle<EERecHitCollection> hrechit_EE_col;
  iEvent.getByLabel(inputTagEndcapRecHits_,hrechit_EE_col);

  typedef edm::View<reco::Candidate> candidateCollection ;
  edm::Handle<candidateCollection> truth_hnd;
  bool isGen = iEvent.getByLabel(inputTagTruth_, truth_hnd);

  if ( !isGen ) { 
    std::cout << "Warning : no Gen !" << std::endl;
    matchWithMC_=false;
  }
  const candidateCollection *truth_candidates=0;
  if(isGen&& matchWithMC_)
    truth_candidates = truth_hnd.product();

  const EcalRecHitCollection& barrelRecHits= *hrechit_EB_col.product();
  const EcalRecHitCollection&  endcapRecHits= *hrechit_EE_col.product();

  bool runFromCandidates=true;
  
  unsigned ncandidates=(runFromCandidates)?pfCandidatesH->size():gsfTracksH->size();
  
  for(unsigned igsf=0;igsf<ncandidates;++igsf)
    {      
      reco::GsfTrackRef theTrackRef;
      std::cout << " Track " << igsf << std::endl;
      if(runFromCandidates)
	{
	  theTrackRef=(*pfCandidatesH)[igsf].gsfTrackRef();
	  if(abs((*pfCandidatesH)[igsf].pdgId()!=11)) continue;
	}
      else
	{
	  theTrackRef=reco::GsfTrackRef(gsfTracksH, igsf);      
	}

      const reco::Candidate * matchCand=0;
      float deltar=0;
      if(matchWithMC_)
	matchCand = mcMatch(truth_candidates,theTrackRef,deltar);
      h700->Fill(deltar);
      
      if(matchWithMC_ && (deltar > deltarMatch_ || matchCand==0)) continue;
      const math::XYZVector & outMom = theTrackRef->outerMomentum() ;
      const math::XYZPoint & outPos = theTrackRef->outerPosition();
      math::XYZVectorD pos(outPos.x(),outPos.y(),outPos.z());
      math::XYZTLorentzVectorD mom;
      mom.SetCoordinates(outMom.x(),outMom.y(),outMom.z(),std::sqrt(outMom.Mag2()));
      std::vector<CaloSegment> segments;
      doExtrapolation(pos,mom,theTrackRef->charge(),segments);
      std::vector<float> energyDeposits;
      std::vector<float> xtalsmiddle;
      int nxtals = doEnergyDeposits(segments,barrelRecHits,endcapRecHits,energyDeposits,xtalsmiddle);
      //      std::cout << "nxtals " << nxtals << std::endl;
      bool barrel=(fabs(outPos.eta())<1.5);
      
      h10->Fill(nxtals);
      if(barrel)
	{
	  h12->Fill(nxtals);
	}
      else
	{
	  h14->Fill(nxtals);
	}

      float etot=0;
      float x0tot=0.;
      double p=std::sqrt(outMom.Mag2());
      for(unsigned in=0;in<energyDeposits.size();++in)
	{
	  //	  std::cout << xtalsmiddle[in] << " " << energyDeposits[in] << std::endl;
	  etot+=energyDeposits[in];
	  if(energyDeposits[in]==-1.) continue;
	  h0->Fill(xtalsmiddle[in],energyDeposits[in]/p);	  
	  if(barrel)
	    {
	      
	      h2->Fill(xtalsmiddle[in],energyDeposits[in]/p);
	      if(nxtals==2)
		{
		  std::cout << in << " " << xtalsmiddle[in] << std::endl; 
		  h22->Fill(xtalsmiddle[in],energyDeposits[in]/p);
		}
	      if(nxtals==3)
		h23->Fill(xtalsmiddle[in],energyDeposits[in]/p);
	    }
	  else
	    {
	      h4->Fill(xtalsmiddle[in],energyDeposits[in]/p);
	      if(nxtals==2)
		h42->Fill(xtalsmiddle[in],energyDeposits[in]/p);
	      if(nxtals==3)
		h43->Fill(xtalsmiddle[in],energyDeposits[in]/p);
	    }
	}
      if(etot==0)
	h30->Fill(outPos.eta());
      // now that etot is known let's start over again
      for(unsigned in=0;in<energyDeposits.size();++in)
	{
	  if(energyDeposits[in]==-1.) continue;
	  h100->Fill(xtalsmiddle[in],energyDeposits[in]/etot);	  
	  if(barrel)
	    {
	      h102->Fill(xtalsmiddle[in],energyDeposits[in]/etot);
	      if(nxtals==2)
		{
		  std::cout << in << " " << xtalsmiddle[in] << std::endl; 
		  h122->Fill(xtalsmiddle[in],energyDeposits[in]/etot);
		}
	      if(nxtals==3)
		h123->Fill(xtalsmiddle[in],energyDeposits[in]/etot);
	    }
	  else
	    {
	      h104->Fill(xtalsmiddle[in],energyDeposits[in]/etot);
	      if(nxtals==2)
		h142->Fill(xtalsmiddle[in],energyDeposits[in]/etot);
	      if(nxtals==3)
		h143->Fill(xtalsmiddle[in],energyDeposits[in]/etot);
	    }
	}    
      // Now try to compare with the expected energy deposit
      double a,b;
      if(!mySimEvent->nTracks()) continue;
      bool status = showerParam(mySimEvent->track(0),a,b);
      if(!status) continue;
      std::cout << " Back from showerParam " << std::endl;
      unsigned nsegments=segments.size();
      x0tot=0.;
      for(unsigned iseg=0; iseg<nsegments; ++iseg)
	{
	  if(segments[iseg].material()==CaloSegment::PbWO4)
	    {
	      if(segments[iseg].entrance().getDetId().subdetId()==EcalEndcap) continue;
	      if(segments[iseg].X0length()==0) continue;
	      std::cout << segments[iseg] << std::endl;
	      std::cout << " Calling expected E " << a << " " << b << " " << p << std::endl;
	      double expectedE=deposit(segments[iseg].sX0Exit(),a,b,segments[iseg].X0length());
	      float x0entrance=x0tot;
	      x0tot += segments[iseg].X0length();
	      std::cout << " Segment " << segments[iseg].sX0Exit()-segments[iseg].X0length() << " " ;
	      std::cout << segments[iseg].sX0Exit() << " calc " << x0entrance << std::endl;
	      float xtalmiddle=x0entrance+0.5*segments[iseg].X0length();
	      
	      float measuredE=getEnergy(segments[iseg],barrelRecHits,endcapRecHits);
	      h200->Fill(xtalmiddle,expectedE);
	      if(nxtals==2)
		{
		  h202->Fill(xtalmiddle,expectedE);
		  h502->Fill(xtalmiddle,measuredE/expectedE/p);
		  h602->Fill(xtalmiddle,measuredE/expectedE/etot);
		}
	      if(nxtals==3)
		{
		  h203->Fill(xtalmiddle,expectedE);
		  h503->Fill(xtalmiddle,measuredE/expectedE/p);
		  h603->Fill(xtalmiddle,measuredE/expectedE/etot);
		}
	    }
	} 
    }
}

void TrackExtra::doExtrapolation(const math::XYZVectorD& pos,const math::XYZTLorentzVectorD & mom, int charge,std::vector<CaloSegment>& segments )
{
  std::vector<SimTrack> mySimTracks;
  int pid=(charge<0) ? 11 : -11 ;
  SimTrack myTrack(pid,mom,0,-1,pos,mom);
  mySimTracks.push_back(myTrack);
  std::vector<SimVertex> mySimVertices;
  mySimVertices.push_back(SimVertex(pos,0.));
  std::cout << " Creating vertex with pos " << pos << " and mom " << mom <<  " eta " << mom.eta() << std::endl;
  mySimEvent->fill(mySimTracks,mySimVertices);
  if(!mySimEvent->nTracks()) return;
  
  FSimTrack  & mySimTrack (mySimEvent->track(0));

  std::cout <<  " Getting vertex " << mySimTrack.vertex().position().Vect() << std::endl;
  
  RawParticle myPart = mySimTrack.ecalEntrance(); 
  std::vector<const RawParticle*> thePart;
  thePart.push_back(&myPart);
  // no preshower
  XYZPoint ecalentrance = myPart.vertex().Vect();
  h32->Fill(ecalentrance.z(),std::sqrt(ecalentrance.x()*ecalentrance.x()+ecalentrance.y()*ecalentrance.y()));
  
  DetId pivot(myGeometry->getClosestCell(ecalentrance, true, mySimTrack.onEcal()==1));
  
  const RandomEngine* random =0 ; 
  EcalHitMaker myGrid(myGeometry,ecalentrance,pivot,mySimTrack.onEcal(),7,0,random);
  myGrid.setTrackParameters(myPart.Vect().Unit(),0.,mySimTrack);

  segments= myGrid.getSegments();
}

float TrackExtra::getEnergy(const CaloSegment& seg, const EcalRecHitCollection & barrel, const EcalRecHitCollection &endcap)
{
  DetId gDetId=seg.entrance().getDetId();
  const EcalRecHitCollection * coll=0;
  if(seg.material()==CaloSegment::PbWO4)
    {
        if(gDetId.subdetId()==EcalBarrel)
	  coll = & barrel;
	if(gDetId.subdetId()==EcalEndcap)
	  coll = & endcap;
	
	EcalRecHitCollection::const_iterator  itcheck=coll->find(gDetId);
	if(itcheck==coll->end())
	  {
	    return 0;
	    } 
	else
	  {
	    return itcheck->energy();
	  }
    }
  return 0.;
}

int TrackExtra::doEnergyDeposits(const std::vector<CaloSegment>& segments, const EcalRecHitCollection & barrel, const EcalRecHitCollection &endcap,std::vector  <float> & deposits,std::vector<float> &xtalpos)
{
  int result=0;
  unsigned nsegments=segments.size();
  float x0tot=0.;

  for(unsigned iseg=0; iseg<nsegments; ++iseg)
    {

      if(segments[iseg].material()==CaloSegment::PbWO4)
	{
	  float x0entrance = x0tot;
	  x0tot += segments[iseg].X0length();	  
	  float xtalmiddle = x0entrance + 0.5*segments[iseg].X0length();

	  xtalmiddle=(segments[iseg].sX0Entrance()+.5*segments[iseg].X0length());

	  DetId gDetId=segments[iseg].entrance().getDetId();
	  const EcalRecHitCollection * coll=0;
	  if(gDetId.subdetId()==EcalBarrel)
	    {	      
	      coll = & barrel;
	      //	      std::cout << " Looking for " << EBDetId(gDetId) << std::endl;
	      ++result;
	      std::cout << " Seg " << iseg << " Mid " << xtalmiddle << " X0tot " << x0tot << std::endl;
	    }
	  else if(gDetId.subdetId()==EcalEndcap)
	    {


	      coll = & endcap;
	      //	      std::cout << " Looking for " << EEDetId(gDetId) << std::endl;
	      ++result;
	    }
	  EcalRecHitCollection::const_iterator  itcheck=coll->find(gDetId);
	  if(itcheck==coll->end())
	    {
	      deposits.push_back(0.);
	      xtalpos.push_back(xtalmiddle);
	      std::cout << " Not found " << std::endl;
	    } 
	  else
	    {
	      deposits.push_back(itcheck->energy());
	      xtalpos.push_back(xtalmiddle);
	      //	      std::cout << " Found " << itcheck->energy() << std::endl;
	    }
	}
      else
	{
	  deposits.push_back(-1.);
	  xtalpos.push_back(0.);	  
	}
    }  
  return result;
}

void TrackExtra::endRun()
{
  ;
}

bool TrackExtra::showerParam(const FSimTrack & mySimTrack,double &a, double &b)
{  
  //  if(!mySimTrack.onEcal() || !mySimTrack.onHcal()) 
  if(!mySimTrack.onEcal() ) 
    {
      std::cout << " Unexpected track " << mySimTrack << std::endl;
      std::cout << " onEcal " <<  mySimTrack.onEcal()<< " onHcal " << mySimTrack.onHcal() <<std::endl;
      return false;
    }
  static std::vector<double> coreParams;
  if(coreParams.size()==0)
    {
      coreParams.push_back(100);
      coreParams.push_back(0.1);
    }
  static std::vector<double> tailParams;
  if(tailParams.size()==0)
    {
      tailParams.push_back(1);
      tailParams.push_back(0.1);
      tailParams.push_back(100);
      tailParams.push_back(1);
    }
  
  // define the calorimeter properties 
  EMECALShowerParametrization 
    showerparam(myGeometry->ecalProperties(mySimTrack.onEcal()), 
		myGeometry->hcalProperties(mySimTrack.onHcal()), 
		myGeometry->layer1Properties(mySimTrack.onLayer1()), 
		myGeometry->layer2Properties(mySimTrack.onLayer2()),
		coreParams,
		tailParams);
  
  double lny = std::log ( mySimTrack.momentum().E() / showerparam.ecalProperties()->criticalEnergy() );
  double theMeanT = showerparam.meanT(lny);
  a = std::exp(showerparam.meanLnT(lny));
  b = (a-1.)/theMeanT;
  std::cout << " Shower param : lny" << lny << " meanT " << theMeanT << " meanLnT " << showerparam.meanLnT(lny) << " a " << a << " b " << b <<std::endl;
  return true;
}

double TrackExtra::deposit(double t, double a, double b, double dt)
{
  std::cout << " Incomplete Gamma " << a << std::endl;
  myIncompleteGamma.a().setValue(a);
  double b1=b*(t-dt);
  double b2=b*t;
  std::cout << " Incomplete Gamma b1 b2 " << b1 << " " << b2 << std::endl;
  return (myIncompleteGamma(b2)-myIncompleteGamma(b1));
}

const reco::Candidate * TrackExtra::mcMatch(const edm::View<reco::Candidate> * genColl,const reco::GsfTrackRef theTrackRef, float& deltarMin)
{
  const reco::Candidate * result=0;
  const math::XYZVector & inMom = theTrackRef->innerMomentum() ;
  deltarMin=9999.;
  for (unsigned int i = 0; i < genColl->size(); i++) {
    const reco::Candidate *particle = &(*genColl)[i];
    float deta=particle->eta()-inMom.eta();
    float dphi=particle->phi()-inMom.phi();
    if(dphi>M_PI)
      dphi-=2.*M_PI;
    if(dphi<-M_PI)
      dphi+=2.*M_PI;
    float dr=std::sqrt(dphi*dphi+deta*deta);
    if(dr<deltarMin)
      {
	deltarMin=dr;
	result=particle;
      }
  } 
  return result;
}

DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(TrackExtra);




