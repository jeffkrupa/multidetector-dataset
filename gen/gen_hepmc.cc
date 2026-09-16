// gen_hepmc.cc -- Pythia8 -> HepMC3 (and HepMC2-ASCII) event generator.
//
// One program, one card file, one seed => one byte-identical event file that is
// fed to BOTH the ATLAS (Athena) and CMS (CMSSW) full simulation chains.
//
// Usage:
//   gen_hepmc <card.cmnd> <seed> <nevents> <output_basename>
// Produces:
//   <output_basename>.hepmc3   HepMC3 ASCII  (Athena >= 22, CMSSW >= 14 HepMC3 reader)
//   <output_basename>.hepmc2   HepMC2 IO_GenEvent ASCII (legacy readers, e.g. CMSSW MCFileSource)
//   <output_basename>.json     provenance: Pythia version, card, seed, xsec, N
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC3.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/WriterAsciiHepMC2.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/Version.h"
#include <fstream>
#include <iostream>
#include <string>

int main(int argc, char** argv) {
  if (argc != 5) {
    std::cerr << "usage: " << argv[0] << " <card.cmnd> <seed> <nevents> <output_basename>\n";
    return 2;
  }
  const std::string card = argv[1];
  const int seed = std::stoi(argv[2]);
  const int nev = std::stoi(argv[3]);
  const std::string base = argv[4];

  Pythia8::Pythia pythia;
  pythia.readFile(card);
  // Seed is a command-line argument so the card is process-only and the seed is explicit provenance.
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = " + std::to_string(seed));
  pythia.readString("Main:numberOfEvents = " + std::to_string(nev));
  if (!pythia.init()) { std::cerr << "Pythia init failed\n"; return 1; }

  Pythia8::Pythia8ToHepMC toHepMC;
  toHepMC.set_print_inconsistency(true);
  toHepMC.set_store_pdf(true);
  toHepMC.set_store_proc(true);
  toHepMC.set_store_xsec(true);

  HepMC3::WriterAscii w3(base + ".hepmc3");
  HepMC3::WriterAsciiHepMC2 w2(base + ".hepmc2");

  int nGood = 0;
  for (int i = 0; i < nev; ++i) {
    if (!pythia.next()) { std::cerr << "event " << i << " failed, retrying\n"; --i; continue; }
    HepMC3::GenEvent ev(HepMC3::Units::GEV, HepMC3::Units::MM);
    toHepMC.fill_next_event(pythia, &ev);
    ev.set_event_number(i + 1);            // 1-based, same number in both experiments' outputs
    w3.write_event(ev);
    w2.write_event(ev);
    ++nGood;
    if ((i + 1) % 100 == 0) std::cerr << "generated " << i + 1 << " events\n";
  }
  w3.close(); w2.close();
  pythia.stat();

  std::ofstream js(base + ".json");
  js << "{\n"
     << "  \"generator\": \"Pythia8\",\n"
     << "  \"pythia_version\": \"" << PYTHIA_VERSION << "\",\n"
     << "  \"hepmc3_version\": \"" << HEPMC3_VERSION << "\",\n"
     << "  \"card\": \"" << card << "\",\n"
     << "  \"seed\": " << seed << ",\n"
     << "  \"nevents\": " << nGood << ",\n"
     << "  \"sigma_gen_mb\": " << pythia.info.sigmaGen() << ",\n"
     << "  \"sigma_err_mb\": " << pythia.info.sigmaErr() << ",\n"
     << "  \"weight_sum\": " << pythia.info.weightSum() << "\n"
     << "}\n";
  std::cout << "wrote " << nGood << " events to " << base << ".hepmc3 / .hepmc2\n";
  return 0;
}
