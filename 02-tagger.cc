#include "matching.hh"
#include "TemplateTagger.hh"
#include "fastjet/FunctionOfPseudoJet.hh"

#include <iomanip> 


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
  const double R = 1.2;
  const double R2 = 0.40;
  const double sigma = 0.333;
  const double minMass = 110., maxMass = 140.;  
  const double ovcut = 0.6; 

  /// Create event containers
  std::vector<fastjet::PseudoJet> particles;
  
  /// An event with three-particles
  particles.push_back( fastjet::PseudoJet( 112.0, -19.8, -56.1, 126.9 ));
  particles.push_back( fastjet::PseudoJet( 110.6,  13.9,  25.3, 114.4 ));
  particles.push_back( fastjet::PseudoJet( 102.3,   5.9,  30.8, 107.1 ));

  /// Set up the template tagger
  TemplateOverlap::TemplateDefinition templDef(R2,sigma,FIXED,CONE,ovcut);
  fastjet::TemplateTagger tagger(argv[1], templDef, fastjet::SelectorMassRange(minMass,maxMass));
  
  /// Jet definitions 
  fastjet::JetDefinition         *jetDef = NULL;
  jetDef = new fastjet::JetDefinition(fastjet::antikt_algorithm, R);
  
  /// Find antikt jets
  vector <fastjet::PseudoJet> jets;
  fastjet::ClusterSequence clustSeq(particles, *jetDef);    
  jets = sorted_by_pt(clustSeq.inclusive_jets(25.));
  
  std::cout << "Ran: " << jetDef->description() << std::endl << std::endl;
  std::cout << "Hardest jet: " << jets[0] << std::endl << std::endl;
  
  /// Now tag the leading jet using template tagger
  fastjet::PseudoJet tagged = tagger(jets[0]);

  /// Show me what you got
  if (tagged == 0){
    std::cout << "No substructure found" << std::endl;
  }  
  else{
  
  fastjet::PseudoJet parent1 = tagged.pieces()[0];
  fastjet::PseudoJet parent2 = tagged.pieces()[1];
  
  double maxOv = tagged.structure_of<fastjet::TemplateTagger>().ov();
  std::vector<fastjet::PseudoJet> maxTempl = tagged.structure_of<fastjet::TemplateTagger>().maxTempl();
  
  std::cout << "Found suitable pair of subjets: " << std::endl;
  std::cout << " " << parent1 << std::endl;
  std::cout << " " << parent2 << std::endl;
  std::cout << "Total = " << std::endl;
  std::cout << " " << tagged << std::endl << std::endl;
  
  std::cout << " The best-matched templates are: "<< "(Ov2 = " << maxOv <<")" << std::endl <<std::endl;     
  for (unsigned i = 0; i < maxTempl.size(); ++i)
    std::cout << maxTempl[i] << std::endl << std::endl;
  }

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
