#ifndef TRACKEXTRA_H
#define TRACKEXTRA_H
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "FastSimulation/CaloGeometryTools/interface/CaloSegment.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "CLHEP/GenericFunctions/IncompleteGamma.hh"
#include "FastSimulation/Event/interface/FSimTrack.h"
#include <iostream>
#include <string>
#include <map>

class CaloGeometryHelper;
class FSimEvent;
class reco::Candidate;
typedef math::XYZVector XYZVector;
typedef math::XYZVector XYZPoint;

class TrackExtra : public edm::EDAnalyzer
{
 public:
  explicit TrackExtra(const edm::ParameterSet&);
  ~TrackExtra();
  virtual void beginRun(edm::Run const&, edm::EventSetup const& );
  virtual void analyze(const edm::Event & iEvent,const edm::EventSetup & c);
  virtual void endRun();

 private:
  void doExtrapolation(const math::XYZVectorD& pos,const math::XYZTLorentzVectorD & mom,int charge,std::vector<CaloSegment>& segments  );
  float getEnergy(const CaloSegment& seg, const EcalRecHitCollection & barrel, const EcalRecHitCollection &endcap);
  int doEnergyDeposits(const std::vector<CaloSegment>& segments, const EcalRecHitCollection & barrel, const EcalRecHitCollection &endcap,std::vector<float> & deposits,std::vector<float>&);
  bool  showerParam(const FSimTrack & mySimTrack, double &a, double &b);
  double deposit(double t, double a, double b, double dt) ;
  const reco::Candidate * mcMatch(const edm::View<reco::Candidate> *,const reco::GsfTrackRef,float &);

 private:
  edm::InputTag inputTagGSFTracks_;
  edm::InputTag inputTagBarrelRecHits_;
  edm::InputTag inputTagEndcapRecHits_;
  edm::InputTag inputTagPFCandidates_;
  edm::InputTag inputTagTruth_;
  Genfun::IncompleteGamma myIncompleteGamma;
  bool matchWithMC_;
  double deltarMatch_ ;

  CaloGeometryHelper * myGeometry;
  FSimEvent * mySimEvent;
  DQMStore * dbe;
  MonitorElement* h0,*h2,*h4,*h10,*h12,*h14;
  MonitorElement* h42,*h43;
  MonitorElement* h22,*h23;
  MonitorElement* h100,*h102,*h104;
  MonitorElement* h142,*h143;
  MonitorElement* h122,*h123;
  MonitorElement* h200, *h202, *h203;
  MonitorElement* h30, *h32;
  MonitorElement* h502, *h503;
  MonitorElement* h602, *h603;
  MonitorElement* h700;
};
#endif
