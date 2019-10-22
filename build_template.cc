//  Template Tagger Package
//  Version 1.0.3 (January 1, 2013)
//  Questions/Comments?  jose.juknevich@weizmann.ac.il

// Copyright (c) 2011-13, Mihailo Backovic, Jose Juknevich,.
//
//----------------------------------------------------------------------
// This file is part of the Template Tagger package ("Template Overlap").

#include "TemplateBuilder.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include <iomanip> // for std::setw

TemplateModel to_enum(int n)
{
  switch( n )
  {
    case 0 :  return TOP;
    case 1 : return HIGGS2;
    case 2 : return HIGGS3;
    default: return HIGGS2; //Make 2-body Higgs template default
  }
}


int main(int argc, char* argv[]) {
  
  if (argc < 7) {
    std::cout << std::endl;
    std::cout << "Usage:  "<< argv[0] << " Model <template.dat> R Mass pT nEta [nPt]" << std::endl;
    std::cout << "Model (=0,1,2) is one of the followings (0=TOP, 1=HIGGS2, 2=HIGGS3)" << std::endl;
    std::cout << std::endl;

    return 1;
  }
   /// Set Run Parameters 
  TemplateModel mode = to_enum(atoi(argv[1]));
  std::ofstream oFile(argv[2]);
  double R = atof(argv[3]); // anti-kt parameter
  double Mass = atof(argv[4]);
  double pT = atof(argv[5]);
  double px = pT, py = 0, pz = 0.;
  double ee = sqrt (px * px + py * py + pz * pz + Mass * Mass); 			
  fastjet::PseudoJet axis(px,py,pz,ee);	  
  
  double etaMax = R;
  double phiMax = R;
  double minPt = 10.; 
  int nEta = atoi(argv[6]);
  int nPhi = nEta;
  int nPt = 25; 
  
  if (atoi(argv[1]) == 2 && argc==8) nPt = atoi(argv[7]);
  
  
  TemplateBuilder(oFile, axis, etaMax, phiMax, minPt, nEta, nPhi, nPt, R, mode);

 
  return 0;
}

