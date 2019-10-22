#include "matching.hh"
#include "fastjet/FunctionOfPseudoJet.hh"

#include <iomanip> // for std::setw

using namespace TemplateOverlap;

//----------------------------------------------------------------------
// forward declaration for printing out info about a jet
//----------------------------------------------------------------------
std::ostream & operator<<(std::ostream &, const fastjet::PseudoJet &);


int main(int argc, char* argv[]) {
  
  if (argc !=2) {
    std::cout << "Usage:  "<< argv[0] << " <template.dat>" << std::endl;
    return 1;
  }
   /// Set Run Parameters 
  double R = 1.2; // anti-kt parameter
  double R2 = 0.40; // template subcone radius
  double sigma = 0.333; // template Gaussian width

  /// An event with three-particles
  std::vector<fastjet::PseudoJet> particles;
  particles.push_back( fastjet::PseudoJet( 112.0, -19.8, -56.1, 126.9 ));
  particles.push_back( fastjet::PseudoJet( 110.6,  13.9,  25.3, 114.4 ));
  particles.push_back( fastjet::PseudoJet( 102.3,   5.9,  30.8, 107.1 ));

  /// Create an instance of MatchingMethod for the analysis.  
  /// The ctor also loads the templates
  MatchingMethod myCone(argv[1], TemplateDefinition(R2, sigma)); 
  
  /// Jet definitions 
  fastjet::JetDefinition         *jetDef = NULL;
  jetDef = new fastjet::JetDefinition(fastjet::antikt_algorithm, R);
  
  /// Find antikt jets
  vector <fastjet::PseudoJet> jets;
  fastjet::ClusterSequence clustSeq(particles, *jetDef);    
  jets = sorted_by_pt(clustSeq.inclusive_jets(25.));
    
  /// Find the best matched template
  templ_t result = myCone.getOv(jets[0]);
  std::vector<fastjet::PseudoJet> maxTempl = result.second;
  double maxOv = result.first;
  
  /// Show me what you got
  std::cout << "Hardest jet: " << jets[0] << std::endl << std::endl;
  std::cout << " The best-matched templates are: "<< "(Ov2 = " << maxOv <<")" << std::endl <<std::endl;     
  for (unsigned i = 0; i < maxTempl.size(); ++i)
    std::cout << maxTempl[i] << std::endl << std::endl;
  
  return 0;
}

//----------------------------------------------------------------------
// does the actual work for printing out a jet
//----------------------------------------------------------------------
std::ostream & operator<<(std::ostream & ostr, const fastjet::PseudoJet & jet) {
  ostr << "pt, y, phi =" 
       << " " << std::setw(10) << jet.perp() 
       << " " << std::setw(6) <<  jet.rap()  
       << " " << std::setw(6) <<  jet.phi()  
       << ", mass = " << std::setw(10) << jet.m()
       << ", btag = " << jet.user_index();
  return ostr;	
}
