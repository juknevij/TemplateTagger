//  Template Tagger Package
//  Version 1.0.3 (January 1, 2013)
//  Questions/Comments?  jose.juknevich@weizmann.ac.il

// Copyright (c) 2011-13, Mihailo Backovic, Jose Juknevich,.
//
//----------------------------------------------------------------------
// This file is part of the Template Tagger package ("Template Overlap").
//
//  Template Tagger is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 3 of the License, or
//  (at your option) any later version.
//----------------------------------------------------------------------


#ifndef __TEMPLATE_TAGGER_HH__
#define __TEMPLATE_TAGGER_HH__

#include "matching.hh"

#include "fastjet/FunctionOfPseudoJet.hh"
#include <fastjet/Selector.hh>
#include <fastjet/CompositeJetStructure.hh> // to derive the FilterStructure from CompositeJetStructure
#include <fastjet/tools/Transformer.hh>     // to derive Filter from Transformer



#include <string>
#include <climits>
#include <utility>

static const double TINY = 1e-20;


using namespace TemplateOverlap;


FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh


/// Forward declarations for TemplateTagger and TemplateTaggerStructure
class TemplateTagger;
class TemplateTaggerStructure;

/// This class wraps the core Template Tagger code to provide the fastjet::FunctionOfPseudoJet
/// interface for convenience in larger analyses.  See matching.hh for
/// definitions of Ov_N and the constructor options.

class Noverlap : public FunctionOfPseudoJet<templ_t> {
public:

  Noverlap(double subConeRadius, double sigma,  const string & templateFile, Resolution_scale_scheme variableCone, Jet_shape_scheme mode)
  :  _templateFinder(templateFile, TemplateDefinition(subConeRadius, sigma,variableCone, mode))
{_number_of_templates=_templateFinder.get_number_of_templates();
  std::cout << "Using " << _number_of_templates <<" templates"<<std::endl;
}
   
   /// returns Ov_N, measured on the constituents of this jet
  templ_t result(const PseudoJet& jet) const {return _templateFinder.getOv(jet);}

private:

   mutable TemplateOverlap::MatchingMethod _templateFinder; 
   int _number_of_templates;
};


class planarflow : public FunctionOfPseudoJet<double> {
public:

  planarflow(double r1, double r2, double r3) {_r.push_back(r1);_r.push_back(r2);_r.push_back(r3);}
  planarflow(std::vector<double> r) : _r(r) {}
  planarflow() {}
   
   /// returns Pf, measured on the constituents of this jet
  double result(const PseudoJet& jet) const {
    double Ixy = 0.,Ixx = 0.,Iyy = 0.;
    double etaJet = jet.rapidity();
    double phiJet = jet.phi();
 
    std::vector<fastjet::PseudoJet> particles = jet.constituents();
    _r.resize(particles.size(), 0.);

    for ( unsigned int i = 0; i < particles.size(); ++i) {
      //Get x and y components in that basis.
      double etaNow = particles[i].rapidity() - etaJet;
      double phiNow = particles[i].phi() - phiJet;
      double ee = particles[i].perp();
      if (abs(phiNow) > M_PI) phiNow = 2 * M_PI - abs(phiNow);
      // Calculate matrix to be used in the computation of planar flow
      Ixx += ee * (etaNow  * etaNow + _r[i] * _r[i]/4) ;  
      Iyy += ee * (phiNow * phiNow + _r[i] * _r[i]/4) ;
      Ixy += ee * etaNow * phiNow; 
    }         

    //Return planar flow.
    return(4 * (Ixx * Iyy -Ixy * Ixy)/((Ixx + Iyy )*(Ixx + Iyy )) );  
  }

private:
 mutable std::vector<double> _r;   
};

class angularity : public FunctionOfPseudoJet<double> {
public:

  angularity(int N): _N(N) {}
   
  double result(const PseudoJet& jet) const {
     std::vector<PseudoJet> particles = jet.constituents();
    double sum = 0.;
    for (unsigned int i = 0; i < particles.size(); ++i)
    {
      double Theta = max (sqrt(particles[i].plain_distance(jet)), TINY);
      double cosTheta = cos(Theta);
      double sinTheta = sqrt(1 - cosTheta * cosTheta);
      sum+= particles[i].perp() * pow(sinTheta, _N) * pow(1-cosTheta, 1-_N);
    }
    return (1/jet.m() * sum);
  }

private:
   int _N;
   
};

class thetaBar : public FunctionOfPseudoJet<double> {
public:

  thetaBar() {}
   
  double result(const PseudoJet& jet) const {
    std::vector<PseudoJet> particles = jet.constituents();
    double sumTheta = 0.;
    for (unsigned int i = 0; i < particles.size();++i) 
      sumTheta += sin( max (sqrt(particles[i].plain_distance(jet)), TINY));
	
    return sumTheta;
  }

private:
   int _N;
   
};

class thetaS : public FunctionOfPseudoJet<double> {
public:

  thetaS() {}
   
  double result(const PseudoJet& jet) const {
    std::vector<PseudoJet> particles = jet.constituents();
    int iMin  = 0;
    double pTmin = particles[0].perp();
    for (unsigned int i = 1; i < particles.size();++i) 
    { double pTnow = particles[i].perp();
      if (pTnow < pTmin) {iMin = i; pTmin = pTnow;}
    }
    return sqrt(particles[iMin].plain_distance(jet));
  }
   
};



//---------------------------------------------------------------------------------
// A new implementation of the Template Tagger that does not make use of Noverlap
//---------------------------------------------------------------------------------

class TemplateTagger : public Transformer {
public:
  /// ctor with arguments (see the class description above)
   TemplateTagger(const string & templateFile, const TemplateDefinition & templDef , Selector selector ) : 
    _templateFinder(templateFile, templDef), _selector(selector), 
    _subConeRadius(templDef.r()), _ovcut(templDef.ovcut()), _mode(templDef.mode())
    {_number_of_templates=_templateFinder.get_number_of_templates();
  std::cout << "Using " << _number_of_templates <<" templates"<<std::endl;}
    
  /// returns a textual description of the tagger 
  virtual std::string description() const;

  /// runs the tagger on the given jet and
  /// returns the tagged PseudoJet if successful, a PseudoJet==0 otherwise
  /// (standard access is through operator()).
  virtual PseudoJet result(const PseudoJet & jet) const;

  /// the type of Structure returned
  typedef TemplateTaggerStructure StructureType;

protected:
  
  mutable TemplateOverlap::MatchingMethod _templateFinder; 
  mutable Selector _selector;  /// the subjet selection criterium
  double _subConeRadius; /// The size of the subcone around each template parton
  double _ovcut;                 /// the minimum overlap value
  double _sigma;
  int _number_of_templates;
  Jet_shape_scheme _mode; 
  
   
};

class TemplateTaggerStructure : public CompositeJetStructure{
public:
  /// ctor with pieces initialisation
  TemplateTaggerStructure(const std::vector<PseudoJet> & pieces_in) :
    CompositeJetStructure(pieces_in), _ov(0.0){}

  /// returns the associated N-overlap
  inline double ov() const{return _ov;}
  inline std::vector<PseudoJet> maxTempl() const {return _peak_template;}
  
  virtual std::string description() const { return "Template Tagger"; }


protected:
  double _ov;      ///< the value of the N-overlap
  std::vector<PseudoJet> _peak_template; ///< the best matched template

  // allow the tagger to set these
  friend class TemplateTagger;
};

//----------------------------------------------------------------------
// TemplateTagger class implementation
//----------------------------------------------------------------------

// class description
string TemplateTagger::description() const {
  ostringstream ostr;
  ostr << "Filter with subjet_def = ";
  if (_number_of_templates>0)
    ostr << "Template Tagger with cone";
  ostr<< ", selection " << _selector.description();
  
  return ostr.str();
}


//----------------------------------------------------------------------

// return a vector of subjets, which are the ones that would be kept
// by the filtering
PseudoJet TemplateTagger::result(const PseudoJet &jet) const {
  templ_t my_result =  _templateFinder.getOv(jet);
  if (my_result.first < _ovcut || ! _selector.pass(jet) ) return PseudoJet();
  
    std::vector<fastjet::PseudoJet> subjets;
    fastjet::Selector sel = fastjet::SelectorCircle(_subConeRadius);


  if (_mode == TOP)
  {
    for (unsigned int i=0; i< my_result.second.size(); ++i)
    {
      sel.set_reference(my_result.second[i]);
      if (jet.pieces()[i].has_constituents())
      {
	std::vector<fastjet::PseudoJet> particles =  jet.pieces()[i].constituents();
	subjets.push_back(fastjet::join(sel(particles)));
      }
      else
	subjets.push_back(jet.pieces()[i]);
    }
 
  }
  else
  {
    std::vector<fastjet::PseudoJet> particles =  jet.constituents();
    for (unsigned int i=0; i< my_result.second.size(); ++i)
    {
      sel.set_reference(my_result.second[i]);
      fastjet::PseudoJet subjet = fastjet::join(sel(particles));
      subjets.push_back(subjet);
    }
  }
  
  PseudoJet result = join<StructureType>(subjets);
  
  StructureType * s = (StructureType *) result.structure_non_const_ptr();
  s->_ov = my_result.first;

  s->_peak_template = my_result.second;
  

  return result;
}


FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh

#endif  // __TEMPLATE_TAGGER_HH__
