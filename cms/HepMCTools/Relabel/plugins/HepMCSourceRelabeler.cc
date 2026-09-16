// Re-labels the products of MCFileSource ("source","generator") to the
// standard generator-module layout expected by the CMSSW GEN/SIM/NANO chain:
//   edm::HepMCProduct   generator:unsmeared
//   GenEventInfoProduct generator
//   GenRunInfoProduct   generator   (run)   -- needed by NanoAOD Rivet/HTXS modules
//   GenLumiInfoHeader   generator   (lumi)  -- needed by NanoAOD genWeightsTable
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoHeader.h"

class HepMCSourceRelabeler : public edm::one::EDProducer<edm::EndRunProducer, edm::BeginLuminosityBlockProducer> {
public:
  explicit HepMCSourceRelabeler(const edm::ParameterSet& ps)
      : hepmcToken_(consumes<edm::HepMCProduct>(ps.getParameter<edm::InputTag>("src"))),
        infoToken_(consumes<GenEventInfoProduct>(ps.getParameter<edm::InputTag>("src"))) {
    produces<edm::HepMCProduct>("unsmeared");
    produces<GenEventInfoProduct>();
    produces<GenRunInfoProduct, edm::Transition::EndRun>();
    produces<GenLumiInfoHeader, edm::Transition::BeginLuminosityBlock>();
  }
  void produce(edm::Event& e, const edm::EventSetup&) override {
    e.put(std::make_unique<edm::HepMCProduct>(e.get(hepmcToken_)), "unsmeared");
    e.put(std::make_unique<GenEventInfoProduct>(e.get(infoToken_)));
  }
  void beginLuminosityBlockProduce(edm::LuminosityBlock& lb, const edm::EventSetup&) override {
    lb.put(std::make_unique<GenLumiInfoHeader>());
  }
  void endRunProduce(edm::Run& r, const edm::EventSetup&) override {
    r.put(std::make_unique<GenRunInfoProduct>());
  }
  static void fillDescriptions(edm::ConfigurationDescriptions& d) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src", edm::InputTag("source", "generator"));
    d.addDefault(desc);
  }
private:
  const edm::EDGetTokenT<edm::HepMCProduct> hepmcToken_;
  const edm::EDGetTokenT<GenEventInfoProduct> infoToken_;
};
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(HepMCSourceRelabeler);
