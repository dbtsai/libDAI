#include<iostream>
#include<fstream>
#include<string>

#include<dai/factorgraph.h>
#include<dai/evidence.h>
#include<dai/alldai.h>

using namespace std;
using namespace dai;

void usage(const string& msg) {
  cerr << msg << endl;
  cerr << "Usage:" << endl;
  cerr << " testem factorgraph.fg evidence.tab emconfig.em" << endl;
  exit(1);
}

int main(int argc, char** argv) {
  if (argc != 4) {
    usage("Incorrect number of arguments.");
  }
  
  FactorGraph fg;
  ifstream fgstream(argv[1]);
  fgstream >> fg;

  PropertySet infprops;
  infprops.Set("verbose", (size_t)1);
  infprops.Set("updates", string("HUGIN"));
  InfAlg* inf = newInfAlg("JTREE", fg, infprops);
  inf->init();

  Evidence e;
  ifstream estream(argv[2]);
  e.addEvidenceTabFile(estream, fg);

  cout << "Number of samples: " << e.nrSamples() << endl;
  Evidence::iterator ps = e.begin();
  for (; ps != e.end(); ps++) {
    cout << "Sample " << ps->first << " has " 
	 << ps->second.observations().size() << " observations." << endl;
  }

  ifstream emstream(argv[3]);
  EMAlg em(e, *inf, emstream);

  for (size_t i = 0; i < 10; ++i) {
    Real l = em.iterate();
    cout << "Iteration " << i << " likelihood: " << l <<endl;
  }

  return 0;
}
