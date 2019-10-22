//  Template Tagger Package
//  Version 1.0.3 (January 1, 2013)
//  Questions/Comments?  jose.juknevich@weizmann.ac.il

// Copyright (c) 2011-13, Mihailo Backovic, Jose Juknevich,.
//
//----------------------------------------------------------------------
// This file is part of the Template Tagger package ("Template Overlap").

#ifndef MATCHING_H
#define MATCHING_H

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

using namespace std;

namespace TemplateOverlap{

enum Jet_shape_scheme {CONE, GAUSSIAN, TOP};
enum Resolution_scale_scheme {FIXED, VARIABLE};
enum Particle_id {HADRON, LEPTON, MET};

///Additional info for templates
class MyInfo: public fastjet::PseudoJet::UserInfoBase {
public:  
  MyInfo(const double & radius, const double & sigma) : _radius(radius), _sigma(sigma){}
  double radius() const {return _radius;}
  double sigma() const {return _sigma;}
protected:
  double _radius;
  double _sigma;  
};

/// Parameters that define TemplateOverlap
class TemplateDefinition {
private:
   double _subConeRadius;     // subCone radius
   double _sigma;  // Gaussian energy resolution relative to parton pT
   Resolution_scale_scheme _subConeMode; //Fixed or varying subcones
   Jet_shape_scheme _mode; // Functional measure
   double _ovcut; // minimum Score for the template to be considered a good match
   double _minPtParton; //For infrared safety, partons are not too soft

public:
   TemplateDefinition(const double subConeRadiusIn=0.20, 
	    const double sigmaIn=0.33, 
	    Resolution_scale_scheme variableConeIn=FIXED, Jet_shape_scheme modeIn = CONE, const double ovcut = 0.0, 
	     const double minPtPartonIn = 10.) :
   _subConeRadius(subConeRadiusIn),  _sigma(sigmaIn),
   _subConeMode(variableConeIn), _mode(modeIn), _ovcut(ovcut) ,_minPtParton(minPtPartonIn) {}
    // Returns the value of the template sub cone radius
   double r() const {return _subConeRadius;}
   //Returns the fraction of template parton p_T used as \sigma_a.
   double sigma() const {return _sigma;}
    //Returns the status value for the overlap calculation mode, 
   //i.e. fixed or varying cone
   Resolution_scale_scheme variableCone() const {return _subConeMode;}
   //Returns the kernel function mode 
   //i.e. Gaussian or Cone
   Jet_shape_scheme mode() const {return _mode;}
   //Returns the minimum p_T of the template parton with smaller transverse momentum. 
   double minPtParton() const {return _minPtParton;}
   
   double ovcut() const {return _ovcut;}
};

/// Energy profile inside cones

double eShape (double x)
{
 // double aux =  0.422258 - 0.00377161* x + 0.0000174186 * x*x -   3.50639e-8 * x *x*x + 2.53302e-11 * x*x*x*x;
   double aux =  0.399 - 0.00535*x + 0.0000352*x*x - 1.16318E-7*x*x*x + 1.85503E-10*x*x*x*x - 1.12941E-13*x*x*x*x*x;

  
  return aux;
}
/// Calculate template radii from template pTs
std::vector<double> templateRadii (const std::vector<fastjet::PseudoJet> & particles, const double & r0)
{ 
  std::vector<double> radii;
  for (unsigned int i = 0; i < particles.size(); ++i){
    double r = eShape(particles[i].perp()) - eShape(100.) + r0;
  if (r > 0.4) r=0.4;
  if (r < 0.05) r =0.05;
    radii.push_back( r);
    
  }
  return radii;
}

/// A helper class for MatchingMethod
class SingleParticle {

public:

  // Constructors.
  SingleParticle(  double pTIn = 0., double yIn = 0., double phiIn = 0.)
  :  pT(pTIn), y(yIn), phi(phiIn), mult(1), isUsed(false), radius(0.), sigma(0.), id(0) { }
  SingleParticle(const SingleParticle& ssj) :  pT(ssj.pT),
    y(ssj.y), phi(ssj.phi), mult(ssj.mult), isUsed(ssj.isUsed), radius(ssj.radius), sigma(ssj.sigma), id(ssj.id) { }
  SingleParticle& operator=(const SingleParticle& ssj) { if (this != &ssj)
    {  pT = ssj.pT; y = ssj.y; phi = ssj.phi; 
    mult = ssj.mult; isUsed = ssj.isUsed; radius =ssj.radius; sigma=ssj.sigma; id=ssj.id;} return *this; }

  // Properties of particle.
  double pT, y, phi;
  int    mult; 
  bool   isUsed;
  double radius; //For templates
  double sigma; // For templates
  int id; // particle id


  double deltaR2(const SingleParticle & other) const {
    double  dPhi = abs(phi - other.phi );
    if (dPhi > M_PI) dPhi = 2. * M_PI - dPhi;
    double  dEta =  y - other.y;
    return (dEta * dEta + dPhi * dPhi );
  }
 
  double deltaPhi(const SingleParticle & other) const {
    double  dPhi = abs(phi - other.phi );
    return (dPhi > M_PI) ?  2. * M_PI - dPhi : dPhi;
  } 
  
};

typedef std::vector<SingleParticle> jet_t;
typedef std::pair<double, std::vector<fastjet::PseudoJet> > templ_t;


std::vector<fastjet::PseudoJet> precluster (const std::vector<fastjet::PseudoJet> & constituents, const double & r_cell, const double & pt_threshold)
{
  std::vector<fastjet::PseudoJet> my_cells;
  std::vector<fastjet::PseudoJet> particles = sorted_by_pt(constituents);
  my_cells.push_back(sorted_by_pt(particles)[0]);		
  for (unsigned int i = 1; i < particles.size(); i++)
  {
    bool found = false;
    
    std::vector<fastjet::PseudoJet >::iterator it = my_cells.begin(); 
    for (; it != my_cells.end(); ++it)
    {
      if(particles[i].plain_distance(*it) < r_cell * r_cell)
      {
	found = true;
	double eNow  = particles[i].e() +(*it).e();
	double pxNow = particles[i].px()+(*it).px();
	double pyNow = particles[i].py()+(*it).py();
	double pzNow = particles[i].pz()+(*it).pz();
	double pAbs = sqrt(pxNow*pxNow+pyNow*pyNow+pzNow*pzNow); 

	double pxNew = eNow*pxNow/pAbs;
	double pyNew = eNow*pyNow/pAbs;
	double pzNew = eNow*pzNow/pAbs;
	
	fastjet::PseudoJet jetNew(pxNew,pyNew,pzNew,eNow);
	my_cells.erase(it);
	my_cells.push_back(jetNew);
	break;
      }
    }
    if(!found)
    {
      my_cells.push_back(particles[i]);
      my_cells = sorted_by_pt(my_cells);
    }
    
  }
  
  std::vector<fastjet::PseudoJet >::iterator it = my_cells.begin(); 
  for (; it != my_cells.end();){
    if(((*it).perp() < pt_threshold))
      it = my_cells.erase(it);
    else
      ++it;
  }
  return(my_cells);	
}




class MeasureFunctor {

   protected:
      MeasureFunctor() {}

   public:
      virtual  ~MeasureFunctor() {}

      virtual double distance(const SingleParticle& particle, const SingleParticle& axis) = 0;
      virtual double numerator(const SingleParticle& particle, const SingleParticle& axis) = 0;

      double overlap(jet_t & jet, const jet_t& axes, const double minNow);

};


double MeasureFunctor::overlap(jet_t& jet, const jet_t& axes, const double minNow) {
  // Calculates Ov
   double ov = 0.0;
   
    for(jet_t::const_iterator t_particle = axes.begin(); t_particle != axes.end(); ++t_particle) 
    {
      double eTjetNow = 0.;
      //  Sum up unused particles within required distance of parton.
      for ( jet_t::iterator particle = jet.begin(); particle != jet.end(); ++particle) 
      {
	if ((*particle).isUsed) continue;
	if ((*particle).id == MET && (*particle).deltaPhi(*t_particle) > (*t_particle).radius) continue;

	if ((*particle).id != MET && (*particle).deltaR2(*t_particle) > (*t_particle).radius * (*t_particle).radius ) continue;
	eTjetNow += numerator((*particle),(*t_particle));
	(*particle).isUsed = true;
      }
      double eTcenterNow = (*t_particle).pT;

      ov += 0.5 * 1/ ( ((*t_particle).sigma)*((*t_particle).sigma) ) * (eTjetNow - eTcenterNow) * (eTjetNow - eTcenterNow)/(eTcenterNow * eTcenterNow); 
      
      if (ov > 4.6 ) break;     // if partial sum too large, no chance this will give large overlap 
      
      if (ov > minNow) break;  // if partial sum larger than current minimum, exit
      
    } 
   
   return ov;
}

///Default Measure 
class DefaultMeasure : public MeasureFunctor {

   private:
      TemplateDefinition _templ_def;

   public:
      DefaultMeasure(){}
      virtual double distance(const SingleParticle& particle, const SingleParticle& axis) {
         return std::sqrt(particle.deltaR2(axis));
      }

      virtual double numerator(const SingleParticle& particle, const SingleParticle& axis) {
         double deltaR = std::sqrt(particle.deltaR2(axis));
         if (deltaR > axis.radius) return 0.0;
         return particle.pT;
      }

};

///Modified Measure 
class ModifiedMeasure : public MeasureFunctor {

   private:
      TemplateDefinition _templ_def;

   public:
      ModifiedMeasure(){}
      virtual double distance(const SingleParticle& particle, const SingleParticle& axis) {
         return std::sqrt(particle.deltaR2(axis));
      }

      virtual double numerator(const SingleParticle& particle, const SingleParticle& axis) {
	if (particle.id != axis.id) return 0.0;
	switch (particle.id){
	  case HADRON:
	  {
	    double deltaR = std::sqrt(particle.deltaR2(axis));
	    if (deltaR > axis.radius) return 0.0;
	  } return particle.pT; 
	  case LEPTON:
	  {
	    double deltaR = std::sqrt(particle.deltaR2(axis));
	    if (deltaR > axis.radius) return 0.0;
	    return particle.pT; 
	  }
	  case MET:
	  {
	    double deltaPhi = particle.deltaPhi(axis);
	    if (deltaPhi > axis.radius) return 0.0;
	    return particle.pT;  
	  }
	  default:  return 0;
	}
      }

};

/// Gaussian Measure 
class GaussianMeasure : public MeasureFunctor {

   private:
      TemplateDefinition _templ_def;

   public:
      GaussianMeasure(){}
      virtual double distance(const SingleParticle& particle, const SingleParticle& axis) {
         return std::sqrt(particle.deltaR2(axis));
      }

      virtual double numerator(const SingleParticle& particle, const SingleParticle& axis) {
        double etaNow = particle.y;
	double phiNow = particle.phi;
	double omegaNow = particle.radius/sqrt(2);
        return particle.pT * exp (-( etaNow * etaNow + phiNow * phiNow)/(2 * omegaNow * omegaNow)   );
      }

};

/// Reset jet flags.
void reset(jet_t & jet)
{
  for (jet_t::iterator j = jet.begin(); j != jet.end(); ++j) 
    if ((*j).isUsed) (*j).isUsed = false;	        
}

/// Sum pT around template momenta and return overlap.                  
double overlapDistance (jet_t & jet1, jet_t & jet2, double R0) 
{
  double sum = 0.;
    for(jet_t::const_iterator k = jet1.begin(); k != jet1.end(); ++k) {
            double eTjetNow = 0.;              
       //  Sum up unused particles within required distance of parton.
    for (jet_t::iterator j = jet2.begin(); j !=jet2.end(); ++j) {
      if ((*j).isUsed) continue;
    if ((*j).deltaR2(*k) > R0 * R0 ) continue;
         eTjetNow += (*j).pT;
        (*j).isUsed = true;
    }              
    double eTcenterNow = (*k).pT;
    //change sigma to E/3
    sum += 4.5 * (eTjetNow - eTcenterNow) * (eTjetNow - eTcenterNow)/(eTcenterNow * eTcenterNow);

}

reset(jet2);
return (exp(-sum));
    }
    
inline double fun (const double ei, const double et, const double s)
{
  return (ei - et)*(ei - et)/ (2*s*s);
}

/// Old overlap function 
void matchTemplate(jet_t & jet, const  std::vector<jet_t> & templates, double & maxVal, int & maxLoc, const double & pT_correction)
{ 
  double minSum = 100;
  for (unsigned int i = 0; i < templates.size(); ++i)
  {
    double sum = 0.;
    for(jet_t::const_iterator t_particle = templates[i].begin(); t_particle != templates[i].end(); ++t_particle) 
    {
      double eTjetNow = 0.;
      //  Sum up unused particles within required distance of parton.
      for ( jet_t::iterator particle = jet.begin(); particle != jet.end(); ++particle) 
      {
	if ((*particle).isUsed) continue;
	if ((*particle).deltaR2(*t_particle) > (*t_particle).radius * (*t_particle).radius ) continue;
        eTjetNow += (*particle).pT;
	(*particle).isUsed = true;
      }
      // Note we rescale the template pT by a factor pT_correction to account for out-of-cone QCD radiation
      // This is valid for VARIABLE cone
      double eTcenterNow = pT_correction * (*t_particle).pT;
      double sigmaNow = (*t_particle).sigma * eTcenterNow;
   //   sum += 0.5 * (eTjetNow - eTcenterNow) * (eTjetNow - eTcenterNow)/( sigmaNow * sigmaNow ); 
      sum += fun (eTjetNow, eTcenterNow, sigmaNow);
      
      if (sum > 4.6 ) break;     // if partial overlap is too low, we skip the rest. 
      //if (sum > 0.7 ) break;   // Can be selected based on an a priori Ov cut?
      if (sum > minSum) break;   // if overlap already lower than current maximum, exit. 
      
    } 
// Clean flags. 
  reset(jet);
  if (sum < minSum) {minSum = sum; maxLoc = i;}
  //if (minSum < 0.2) break; // Current minimum is already acceptable 
    
  }
 
  maxVal = exp(-minSum);

}


/// New overlap function 
/// We allow for some flexibility in defining new functional measures
void matchTemplate2(jet_t & jet, const  std::vector<jet_t> & templates, double & maxVal, int & maxLoc, MeasureFunctor* measure)
{ 
  double minSum = 1000;
  maxLoc = -1;
  for (int i = 0; i < int (templates.size()); ++i)
  { 
    double sum = measure->overlap(jet,templates[i], minSum);
    if (sum < minSum) {maxLoc = i; minSum = sum;}
  //  if (minSum < 0.2) break;
    reset(jet); //Clean flags. 
  }
  maxVal = exp(-minSum);

}

/// New overlap function 
/// We allow for some flexibility in defining new functional measures
void matchTemplate2(std::vector <jet_t> & jets, const  std::vector<jet_t> & templates, double & maxVal, int & maxLoc, MeasureFunctor* measure)
{ 
  double minSum = 1000;
  maxLoc = -1;
  
  for (int i = 0; i < int (templates.size()); ++i)
  { //Loop over the pieces of the event or jet,
    double sum = 0;
    for (unsigned int j = 0; j < jets.size(); ++j )
    {
     sum += measure->overlap(jets[j], jet_t (1,templates[i][j]), minSum);    
     reset(jets[j]); //Clean flags. 

    }
    if (sum < minSum) {maxLoc = i; minSum = sum;}
  //  if (minSum < 0.2) break;
  }
  maxVal = exp(-minSum);

}

void maximize(const vector<double>  & result, double & maxVal, int & maxLoc)
{
  maxVal = 0;
  maxLoc = 0;

  for (int i = 0; i < int (result.size()); ++i)
  {
    double thisVal = result[i];
    if (thisVal > maxVal){ maxLoc = i; maxVal = thisVal;}
  }
}


//--------------------------------------------------------------------------


std::vector<fastjet::PseudoJet> ConvertToPseudoJet(const jet_t& particles) 
{
   std::vector<fastjet::PseudoJet> Jets;
   for (jet_t::const_iterator part = particles.begin(); part != particles.end(); part++) {    
      fastjet::PseudoJet temp = fastjet::PtYPhiM((*part).pT, (*part).y, (*part).phi);
      Jets.push_back(temp);
   }
   return Jets;
}

std::vector<fastjet::PseudoJet> ConvertToPseudoJet(const jet_t& particles, const fastjet::PseudoJet & tilt_axis) 
{
   std::vector<fastjet::PseudoJet> Jets;
   double etaJet = tilt_axis.rapidity();
   double phiJet = tilt_axis.phi_std();
   for (jet_t::const_iterator part = particles.begin(); part != particles.end(); part++) 
   {     
      fastjet::PseudoJet temp = fastjet::PtYPhiM((*part).pT, (*part).y + etaJet, (*part).phi + phiJet); 
      
      temp.set_user_info(new MyInfo((*part).radius, (*part).sigma));
      Jets.push_back(temp);
   }
   return Jets;
}

jet_t ConvertToMat(const fastjet::PseudoJet& jet)
{
  std::vector<fastjet::PseudoJet> particles = jet.constituents();
  jet_t jet_mat; 
  
  double jetRapidity = jet.rapidity();	  
 // double jetPhi =      jet.phi();

  for (unsigned int i = 0; i < particles.size(); ++i) {	
    double etaNow = particles[i].rapidity() - jetRapidity;
 //   double phiNow = particles[i].phi() - jetPhi;
    double phiNow = jet.delta_phi_to(particles[i]);
    double pTnow = particles[i].perp();
    
    SingleParticle particleNow( pTnow, etaNow, phiNow);
   // particleNow.id = particles[i].user_index();
    
    jet_mat.push_back( particleNow);
  }
  return jet_mat;
}

jet_t ConvertToMat(const fastjet::PseudoJet& jet, const fastjet::PseudoJet & axis)
{
  std::vector<fastjet::PseudoJet> particles;
  if( jet.has_constituents()) 
    particles = jet.constituents();
  else
    particles.push_back(jet);
  jet_t jet_mat; 
  
  double jetRapidity = axis.rapidity();	  
//  double jetPhi =      axis.phi();

  for (unsigned int i = 0; i < particles.size(); ++i) {	
    double etaNow = particles[i].rapidity() - jetRapidity;
   // double phiNow = particles[i].phi();
    double phiNow = axis.delta_phi_to(particles[i]);
    double pTnow = particles[i].perp();
    
    SingleParticle particleNow( pTnow, etaNow, phiNow);
   // particleNow.id = particles[i].user_index();
    
    jet_mat.push_back( particleNow);
  }
  return jet_mat;
}

//--------------------------------------------------------------------------
double DeltaPhi(double phi1, double phi2)
{
  double dphi = phi1-phi2;
  
  if (dphi >  M_PI) dphi -= 2* M_PI;
  if (dphi < -M_PI) dphi += 2* M_PI;
 // if (abs (temp) > M_PI) temp = 2 * M_PI - abs(temp);
return dphi;  
}

/// The main analysis class
class MatchingMethod {
private:

   ifstream _mystream;
   vector<jet_t> _templates; 
   double _subConeRadius;
   double _sigma;
   double _minPtParton;
   double _maxVal;
   int _maxLoc; 
   Resolution_scale_scheme _subConeMode;
   MeasureFunctor * _functor;
   int _number_of_templates;
   
   int templateRead(ifstream& is);
   bool checkTemplate (const jet_t & probe_template);
   bool _top_decay_mode; 
   
   

public:
   MatchingMethod(const string & templateFile, const TemplateDefinition & myTemplDef) : 
   _mystream(templateFile.c_str()), _subConeRadius(myTemplDef.r()),
   _sigma(myTemplDef.sigma()), _minPtParton(myTemplDef.minPtParton()),
   _subConeMode(myTemplDef.variableCone()) {
     _number_of_templates = templateRead(_mystream);
     switch (myTemplDef.mode()){
       case CONE :
	  _functor = new DefaultMeasure;
	  _top_decay_mode = false;
	  break;
       case GAUSSIAN:	
	  _functor = new GaussianMeasure;
	  _top_decay_mode = false;
	  break;
       case TOP:
	  _functor = new DefaultMeasure;
	 _top_decay_mode = true;
	 break;
       default :
	  _functor = 0;
	  break;
	  
     }
  }
   
   ~MatchingMethod() { delete _functor;}

   templ_t getOv (const fastjet::PseudoJet & jet);
   int get_number_of_templates() {return _number_of_templates;}
  
};

 /// Main method to perform template matching
templ_t  MatchingMethod::getOv(const fastjet::PseudoJet& jet)
{ 
  /// Defining TemplateOverlap parameters
  /// Using normalized step function around template momenta
  /// DefaultMeasure StepFunctionMeasure; 
  
   /// Source jet to matrix
  jet_t jet_mat; 
  jet_mat = ConvertToMat(jet); 
  
  /// Correct for QCD radiation outside subcones 
  /// when using VARIABLE subcone mode
  /// Note in VARIABLE subcone mode, radii are chosen to collect 80% of parton energy
  double out_of_cone_correction = (_subConeMode == VARIABLE) ? 0.8 : 1.0;
  
  /// Do the matching and localize the best match
  matchTemplate(jet_mat, _templates, _maxVal, _maxLoc, out_of_cone_correction);
 //  matchTemplate(jet_mat, _templates, _maxVal, _maxLoc, _functor);

  /// Convert to PseudoJet and return result
  std::vector<fastjet::PseudoJet> peak_template;
  peak_template = ConvertToPseudoJet(_templates[_maxLoc], jet);  
 
  return std::make_pair<double, std::vector<fastjet::PseudoJet> >(_maxVal, peak_template);
} 

/// Load templates from file.
int MatchingMethod::templateRead(ifstream& is) {

  string line, tag;
 
 //Event particlesSave;
 
jet_t particlesSave; 
_templates.clear();

  while (getline(is,line)){
  istringstream getpro(line);
  int nupSave;
  getpro >> nupSave;
  if (!getpro) return false;

  // Reset particlesSave vector.
  particlesSave.clear(); 

  // Read in particle info one by one, and store it.
  double pup1, pup2, pup3, pup4;
  
   if (!getline(is, line)) return false;
    istringstream getall(line);
    getall  >> pup1 >> pup2 >> pup3 >> pup4;
  
  for (int ip = 0; ip < nupSave; ++ip) { 
    if (!getline(is, line)) return false;
    istringstream getall(line);
    getall  >> pup1 >> pup2 >> pup3 >> pup4;
    if (!getall) return false;   
    
   fastjet::PseudoJet pTemp(pup1,pup2,pup3,pup4);
   SingleParticle pAux= SingleParticle( pTemp.perp(), pTemp.rapidity(), pTemp.phi_std());
   pAux.radius = (_subConeMode==VARIABLE) ? 
	_subConeRadius + (eShape(pTemp.perp())-eShape(100.) )
	:  _subConeRadius ;
   if (pAux.radius< 0.05) pAux.radius = 0.05; // Minimum subcone radius according to Hcal granularity 
   if (pAux.radius > 0.4)  pAux.radius = 0.4;

   pAux.sigma = _sigma;
   particlesSave.push_back(pAux);
  }
  
 
  if (checkTemplate (particlesSave) ) 
  {
    _templates.push_back(particlesSave);
  }
  
  // Reading worked.
  
}
  return (_templates.size());

}

/// Make sure that template subcones do not overlap.
bool MatchingMethod::checkTemplate(const jet_t& probe_template)
{  
  bool isOK = true; 
  for (jet_t::const_iterator i = probe_template.begin() ; i != probe_template.end() && isOK ; ++i)  {  
    if ((*i).pT < _minPtParton) { isOK = false; break;}
    for (jet_t::const_iterator j = i + 1 ; j != probe_template.end(); ++j){	
      if((*i).deltaR2((*j)) <   ((*i).radius+(*j).radius) *  ((*i).radius+(*j).radius))
	{ isOK = false; break ;}
    }
  } 
  return (isOK);
}




} //end namespace


#endif // end MATCHING_H
