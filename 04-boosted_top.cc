#include "matching.hh"
#include "TemplateTagger.hh"
#include "Settings.hh"

#include "fastjet/FunctionOfPseudoJet.hh"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <time.h> // For timing
#include <iomanip> 

using namespace std;

using namespace TemplateOverlap;

unsigned int WhichPtBin (double PtJet, double PtMin, double PtBinWidth, int nptbins)
{
if (PtJet < PtMin) return 0;
int i =  int(floor  ((PtJet - PtMin)/PtBinWidth));
if (i > nptbins) i = nptbins;
if (i<0) return 0;
return i;
}

//==========================================================================

// Method to pick a number according to a Poissonian distribution.

int poisson(double nAvg, double rndm) {
	
	if (nAvg < 0.0001) return 0;

  // Set maximum to avoid overflow.
  const int NMAX = 100;

  // Random number.
  double rPoisson = rndm * exp(nAvg);
  // Initialize.
  double rSum  = 0.;
  double rTerm = 1.;
  
  // Add to sum and check whether done.
  for (int i = 0; i < NMAX; ) {
    rSum += rTerm;
    if (rSum > rPoisson) return i;

    // Evaluate next term. 
    ++i;
    rTerm *= nAvg / i;
  }

  // Emergency return.
  return NMAX; 
}

//==========================================================================

// Generate a random number between 0 and 1
// return a uniform number in [0,1].
double unifRand()
{
    return rand() / double(RAND_MAX);
}
//
// Generate a random number in a real interval.
// param a one end point of the interval
// param b the other end of the interval
// return a inform rand numberin [a,b].
double unifRand(double a, double b)
{
    return (b-a)*unifRand() + a;
}




///  Helper function for printing jets
void PrintJets(const vector <fastjet::PseudoJet>& jets) {
   printf("%5s %15s %15s %15s %15s %15s\n","jet #", "rapidity", "phi", "pt","m","e"); // label the columns
   for (unsigned int i = 0; i < jets.size(); i++) {
      printf("%5u %15.8f %15.8f %15.8f %15.8f %15.8f\n",i, jets[i].rap(),jets[i].phi(),jets[i].perp(),jets[i].m(),jets[i].e());
   }
}


//----------------------------------------------------------------------
// forward declaration for printing out info about a jet
//----------------------------------------------------------------------
ostream & operator<<(ostream &, const fastjet::PseudoJet &);


//----------------------------------------------------------------------
// class for reading in a jet
//----------------------------------------------------------------------

bool ReadEvent (ifstream& is, vector<fastjet::PseudoJet> & input_particles, vector<fastjet::PseudoJet> & partons)
{

  string line;
  double pup1, pup2, pup3, pup4;
  int nhep, i, isthep,idhep,evtnumber;
  
  input_particles.resize(0);
  partons.resize(0);
    
  if (!getline(is,line)) return false;	
  istringstream getpro(line);
  getpro  >> evtnumber >> nhep;

  for (int j=0; j<nhep;++j){
    if (!getline(is,line)) return false;	
    istringstream getall(line);
    getall  >> i >> isthep >>idhep>> pup1 >> pup2 >> pup3 >> pup4;
	
    fastjet::PseudoJet pJet = fastjet::PseudoJet(pup1,pup2,pup3,pup4);
    pJet.set_user_index(idhep);      

    if ( isthep == 21 && j > 6 ) {partons.push_back(pJet);continue;}
    if ( isthep == 21) continue;

    if (isthep ==1 )  input_particles.push_back(pJet);
  } 
 
  
  return evtnumber;
}

void splitPartons(const std::vector<fastjet::PseudoJet> & partons, std::vector<fastjet::PseudoJet> & top1, std::vector<fastjet::PseudoJet> & top2)
{
  top1.resize(0);
  top2.resize(0);
  for (unsigned int i=0;i<partons.size();++i)
  {
    int id = partons[i].user_index();
    if (id==5|| id==2 || id==4 || id==-1 ||id==-3|| id==-11 || id==-13) top1.push_back(partons[i]);
    if (id==-5|| id==-2 || id==-4|| id==1 ||id==3 || id==11 || id==13) top2.push_back(partons[i]);
  } 
}

//---------------------------------------------------------
/////////////////Add pileup
//---------------------------------------------------------

 /// Overlay pileup 
 /// Only for demonstration, we can overlay a fixed number of pp interactions
void addPileup(std::vector<fastjet::PseudoJet> & particles, std::ifstream & pileupFile, const int nPileup)
{
  std::vector <fastjet::PseudoJet> pileups, partons, leptons;
  fastjet::PseudoJet met;
  for (int iPileup = 0; iPileup < nPileup; ++iPileup) 
  {
    if(!ReadEvent(pileupFile,pileups, partons)) 
    {
      pileupFile.clear() ;
      pileupFile.seekg(0, ios::beg) ;
      ReadEvent(pileupFile,pileups, partons);
      
    }
    for (unsigned int i = 0; i < pileups.size(); ++i)
      particles.push_back(pileups[i]);
  }
}

/// an example program showing how to use TemplateTagger in FastJet
/// in a semi-realistic Top identification analysis
int main(int argc, char* argv[]) {
  
  if (argc != 3  ){
    cout << "Usage: " << argv[0] << " <param.dat> <output.txt>" << endl;
    return 0;}
  
  ///-----------------------------------------------------------------
  /// Read parameters from param.dat file
  ///----------------------------------------------------------------
 
  Settings settings;
  settings.readFile(argv[1]);
  
  const int width = 12;

  double coneRadius =  settings.parm("coneRadius"), 
      subConeRadius3 = settings.parm("subConeRadius3"),
      sigma = settings.parm("sigma");
  double pTbinWidth = settings.parm ("ptbinwidth");
	
  double minPt = settings.parm("minPt"), 
      maxPt = settings.parm("maxPt"), 
      minMass = settings.parm("minmass"), 
      maxMass = settings.parm("maxmass"); 
    

  double ptjetmin = settings.parm("ptjetmin");
  int nptbins = settings.mode("nptbins");
  
  string fileName3[] = {settings.word("file3b1"),
      settings.word("file3b2"),
      settings.word("file3b3"),
      settings.word("file3b4"),
      settings.word("file3b5"),
      settings.word("file3b6"),
      settings.word("file3b7")};

  
   /// Set Run Parameters 
  double ovcut = 0.; /// For TemplateTagger
  int nPileup = 0; /// For pileup overlay
 
 /// Create best matched template containers
  std::vector<fastjet::PseudoJet> maxTempl;

 /// Create event containers
  std::vector<fastjet::PseudoJet> input_particles, partons;

 /// I/O files
  string fileName = settings.word("eventFile");
  const char*  cstring = fileName.c_str();
  ifstream eventFile(cstring);
   
  fileName = settings.word("pileupFile");
  cstring = fileName.c_str();
  ifstream pileupFile(cstring);
  
  ofstream mOutput(argv[2]);
 
  ///----------------------------------------------------------
  /// Create 7 top taggers for different pT's
  ///----------------------------------------------------------
  std::vector<fastjet::SharedPtr<fastjet::Transformer> > HadTopTaggers;
  //-------------------------------------------------------------
  // Here we define a template tagger using variable size cones
  // This allows us to correct for out-of-cone QCD radiation
  // ------------------------------------------------------------
  TemplateOverlap::TemplateDefinition HadTopDef(subConeRadius3, sigma, VARIABLE, CONE, ovcut);
  fastjet::Selector sel = fastjet::SelectorMassRange(minMass,maxMass) && fastjet::SelectorPtRange(minPt, maxPt);
  
  
  for( int i=0; i<nptbins;++i)
  {
    fastjet::SharedPtr<fastjet::Transformer> tagger(new fastjet::TemplateTagger(fileName3[i], HadTopDef, sel));
    HadTopTaggers.push_back(tagger);
  }
  int Nevent = 0;
  
  clock_t start = clock(); // start clock
  
  // read in input particles
  //----------------------------------------------------------
   while (ReadEvent(eventFile, input_particles, partons)){
    ++Nevent;
  
    /// Jet definitions 
    fastjet::JetDefinition         *jetDef = NULL;
    jetDef = new fastjet::JetDefinition(fastjet::antikt_algorithm, coneRadius,
                                      fastjet::E_scheme, fastjet::Best);

    /// Overlay pileup 
    /// Only for demonstration, we can overlay a fixed number of pp interactions
    if (nPileup > 0) addPileup(input_particles, pileupFile, nPileup);
  
    /// get the resulting jets ordered in pt
    vector <fastjet::PseudoJet> fatjets;
    fastjet::ClusterSequence clustSeq(input_particles, *jetDef);    
    fatjets = sorted_by_pt(clustSeq.inclusive_jets(ptjetmin));
  
    /// Choose a tagger (and a filter) based on jet pT
    int iTemp = WhichPtBin(fatjets[0].perp(), minPt, pTbinWidth, nptbins-1); 
    const fastjet::Transformer & fhad = *HadTopTaggers[iTemp];
  
    /// We will tag the hardest jet in the event
    /// this is our hadronic top jet candidate
    
   // fastjet::PseudoJet HadTopJet = fatjets[0];  
    /// Optionally, we can precluster the original hadrons to speed up the analysis
    fastjet::PseudoJet HadTopJet = fastjet::join(precluster(fatjets[0].constituents(), 0.05, 1.0));
  
    /// Print jet information 
    cout << "Run: " << jetDef->description() << endl << endl;
    cout << "Hadronic Top: " << HadTopJet << endl;
  
    /// Now tag the leading jet using template tagger
    fastjet::PseudoJet tagged = fhad(HadTopJet);
  
    std::vector<fastjet::PseudoJet> top1, top2;
    splitPartons(partons, top1, top2);
    
    /// Several variables
    double JetPt(HadTopJet.perp()), 
      JetMass(HadTopJet.m()), 
      JetRap(HadTopJet.rapidity()),
      ov3(0.),
      JetPf(0.),
      tPf(0.),
      ThetaBar(0.);   
       
    /// Show me what you got
    if (tagged == 0)
    {
      std::cout << "No 3-body substructure found in Hadronic Top" << std::endl;
    }  
    else
    {
  
    /// We extract substructure information from filtered jet  
    /// Tagged jets contain two kind of substructure information
    /// 1) A maximum overlap template
    /// 2) The value of the maximum overlap
    const fastjet::TemplateTagger::StructureType & templ3b = tagged.structure_of<fastjet::TemplateTagger>(); 
    std::vector<fastjet::PseudoJet> my3bTemplate = templ3b.maxTempl();
    fastjet::PseudoJet HadTopTempl = fastjet::join(my3bTemplate);
    ov3 = templ3b.ov();


    /// Compute the actual jet planar flow
    fastjet::planarflow pf;
    JetPf = pf(HadTopJet);
    
    /// Compute several template jet shapes
    /// These variables are computed from the template momenta
    /// The definition of planar flow below smears the energy according to the template subcones
    /// We use it to compute template planar flow
    std::vector<double> radii = templateRadii(my3bTemplate, subConeRadius3);
    fastjet::planarflow tpf(radii);
    tPf = tpf(HadTopTempl);
    
    /// Compute thetaBar
    fastjet::thetaBar tb;
    ThetaBar = tb(HadTopTempl);
    }
    
   ///Write results to output file
    mOutput << fixed;
    mOutput << setprecision (3) <<setw(width) << JetRap
          << setprecision (3) <<setw(width) << JetMass
          << setprecision (3) <<setw(width) << JetPt
          << setprecision (3) <<setw(width) << ov3
	  << setprecision (3) <<setw(width) << JetPf
	  << setprecision (3) <<setw(width) << tPf
	  << setprecision (3) <<setw(width) << ThetaBar

	  << std::endl;
    }
    /// Print execution time
  std::cout << "Time Elapsed: " << ((double)clock() - start) / CLOCKS_PER_SEC << std::endl;  
  
  return 0;
}


//----------------------------------------------------------------------
// does the actual work for printing out a jet
//----------------------------------------------------------------------
ostream & operator<<(ostream & ostr, const fastjet::PseudoJet & jet) {
  ostr << "pt, y, phi =" 
       << " " << setw(10) << jet.perp() 
       << " " << setw(6) <<  jet.rap()  
       << " " << setw(6) <<  jet.phi()  
       << ", mass = " << setw(10) << jet.m()
       << ", btag = " << jet.user_index();
  return ostr;
	
	
}
