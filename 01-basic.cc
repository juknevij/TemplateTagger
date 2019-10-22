#include "matching.hh"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>



using namespace std;
using namespace TemplateOverlap;


///  Helper function for printing jets
void PrintJets(const vector <fastjet::PseudoJet>& jets) {
   printf("%5s %15s %15s %15s %15s %15s\n","jet #", "rapidity", "phi", "pt","m","e"); // label the columns
   for (unsigned int i = 0; i < jets.size(); i++) {
      printf("%5u %15.8f %15.8f %15.8f %15.8f %15.8f\n",i, jets[i].rap(),jets[i].phi(),jets[i].perp(),jets[i].m(),jets[i].e());
   }
}



/// Read input to fastjet
bool ReadEvent (ifstream& is,  
		vector<fastjet::PseudoJet> & input_particles, 
		vector<fastjet::PseudoJet> & hard_partons)
{

  string line;
  double pup1, pup2, pup3, pup4;
  int nhep, i, isthep,idhep,evtnumber;
  
  input_particles.resize(0);
  hard_partons.resize(0);
  
  if (!getline(is,line)) return false;	
  istringstream getpro(line);
  getpro  >> evtnumber >> nhep;

  for (int j=0; j<nhep;++j){
    if (!getline(is,line)) return false;	
    istringstream getall(line);
    getall  >> i >> isthep >>idhep>> pup1 >> pup2 >> pup3 >> pup4;
	
    if ( isthep == 21 )
    {
      fastjet::PseudoJet pJet = fastjet::PseudoJet(pup1,pup2,pup3,pup4);
      pJet.set_user_index(idhep);  
      hard_partons.push_back(pJet);
    } 

    if (fabs(idhep) == 12 || fabs(idhep) == 14 || fabs(idhep) == 16) continue;
    
     if ( isthep == 1 )
    {
      fastjet::PseudoJet pJet = fastjet::PseudoJet(pup1,pup2,pup3,pup4);
      pJet.set_user_index(idhep);  
      input_particles.push_back(pJet);
    } 
  } 

return evtnumber;

}



int main(int argc, char* argv[]) {
  
  if (argc !=4) {
    cout << "Usage:  "<< argv[0] << " <input.dat> <template.dat> <output.dat>" << endl;
    return 1;
  }
   /// Set Run Parameters 
  double subConeRadius = 0.20;
  double coneRadius = 1.0;
  double sigma = 0.333 ;
  TemplateDefinition myParms(subConeRadius, sigma);

 /// Create event containers
  std::vector<fastjet::PseudoJet> input_particles, partons;
 
 /// Create an instance of Cone for the analysis.  
  /// The ctor also loads the templates
  MatchingMethod myCone(argv[2], myParms); 

 /// Load event
  ifstream file(argv[1]);
  ReadEvent(file, input_particles, partons);
  
   /// Jet definitions 
  fastjet::JetDefinition         *jetDef = NULL;
  jetDef = new fastjet::JetDefinition(fastjet::antikt_algorithm, coneRadius,
                                      fastjet::E_scheme, fastjet::Best);
 /// Find antikt jets
  vector <fastjet::PseudoJet> fatjets;
  fastjet::ClusterSequence clustSeq(input_particles, *jetDef);    
  fatjets = sorted_by_pt(clustSeq.inclusive_jets());
 
 ///  Do the matching and localize the best match
 templ_t maxTempl = myCone.getOv(fatjets[0]);
  
  
  /// Show me what you got
  
  printf("------------------------ Template Tagger v0.3.0 --------------------------------------"); printf("\n");
  cout << "pT (GeV) = "<< fatjets[0].perp() << ", m (GeV) = " << fatjets[0].m() << endl;
  cout << "subConeRadius = " << subConeRadius << ", Sigma = "<< sigma << endl;
  printf("-------------------------------------------------------------------------------------"); printf("\n");
  std::cout << "Max. overlap (Ov"<<maxTempl.second.size()<<"):" << maxTempl.first << std::endl;

  cout << "Max. overlap template:" << endl;
  PrintJets(maxTempl.second);
  
  
  return 0;
}



