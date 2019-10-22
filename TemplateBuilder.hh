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

#include <iomanip>

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

enum TemplateModel {TOP, HIGGS2, HIGGS3};

std::ostream & operator<<(std::ostream & ostr, const fastjet::PseudoJet & jet) {
  
  ostr << " " << std::setw(7) <<  jet.px() 
       << " " << std::setw(7) <<  jet.py()  
       << " " << std::setw(7) <<  jet.pz()  
       << " " << std::setw(7) <<  jet.e()  ;
  return ostr;
	
	
}

double  operator *(const fastjet::PseudoJet & jet1, const fastjet::PseudoJet & jet2) {
  
  return dot_product(jet1,jet2);
	
	
}


void TemplateBuilder (std::ofstream & File,
		      const fastjet::PseudoJet & axis,
		      const double etaMax, 
		      const double phiMax, 
		      const double minPt, 
		      const int nEta, 
		      const int nPhi,
		      const int nPt, 
		      const double R, 
		      TemplateModel mode
 		    )
{
  
  fastjet::Strategy               strategy = fastjet::Best;
  fastjet::RecombinationScheme    recombScheme = fastjet::E_scheme;
  fastjet::JetDefinition         *jetDef = NULL;
  jetDef = new fastjet::JetDefinition(fastjet::antikt_algorithm, R,
                                      recombScheme, strategy);
  
    std::vector <fastjet::PseudoJet> fjInputs;

  
  int nJets = 0;
  
  const double mW = 80.; 
  const double Mass = axis.m();
  
  const double  pTmin = 0., pTmax = axis.perp();
  
  
  switch(mode){
    case HIGGS2:
     for( int iEtaNow2 = 0; iEtaNow2 < nEta; ++iEtaNow2 )	{
      for( int iPhiNow2 = 0; iPhiNow2 < nPhi; ++iPhiNow2 )	{
	double etaCell2 = (etaMax / nEta) * (2 * iEtaNow2 - 1 - nEta);
	double phiCell2 = (phiMax / nPhi) * (2 * iPhiNow2 - 1 - nPhi);
	double pTcell2 = Mass * Mass / ( 2 * (axis.e() * cosh(etaCell2) - axis.perp() * cos(phiCell2)) );
	fastjet::PseudoJet p2 = fastjet::PtYPhiM(pTcell2, etaCell2, phiCell2);
	fastjet::PseudoJet p1 = axis -  p2;
	fjInputs.resize(0);
	if (p1.perp() < minPt) continue;
	if (p2.perp() < minPt) continue;
	fjInputs.push_back(p1); 
	fjInputs.push_back(p2); 
	std::vector <fastjet::PseudoJet> inclusiveJets, sortedJets;
	fastjet::ClusterSequence clustSeq(fjInputs, *jetDef);
	inclusiveJets = clustSeq.inclusive_jets();
	sortedJets    = sorted_by_pt(inclusiveJets);  
	if (sortedJets.size()>1) continue;
	File << int(2) <<std::endl;
	File << axis << std::endl;
	File << p1 << std::endl;
	File << p2 << std::endl;
	++nJets;}}
	std::cout << "Found "<<nJets << std::endl;
	break;
    case HIGGS3:
      for( int ipTnow = 0; ipTnow < nPt; ++ipTnow )	{
  for( int iEtaNow1 = 0; iEtaNow1 < nEta; ++iEtaNow1 )	{
    for( int iPhiNow1 = 0; iPhiNow1 < nPhi; ++iPhiNow1 )	{
	for( int iEtaNow2 = 0; iEtaNow2 < nEta; ++iEtaNow2 )	{
	  for( int iPhiNow2 = 0; iPhiNow2 < nPhi; ++iPhiNow2 )	{	 
	    double etaCell1 = (etaMax / nEta) * (2 * iEtaNow1 - 1 - nEta);
	    double phiCell1 = (phiMax / nPhi) * (2 * iPhiNow1 - 1 - nPhi);
	    double etaCell2 = (etaMax / nEta) * (2 * iEtaNow2 - 1 - nEta);
	    double phiCell2 = (phiMax / nPhi) * (2 * iPhiNow2 - 1 - nPhi);
	    double pTcell1 = (  (pTmax - pTmin)/nPt)  *  ipTnow + pTmin;
	    fastjet::PseudoJet p1 = fastjet::PtYPhiM( pTcell1, etaCell1, phiCell1);
	    double p1P =axis * p1;
	    double pTcell2 = (Mass * Mass - 2 * p1P) / ( 2 * pTcell1 * (cos(phiCell1 -phiCell2) - cosh(etaCell1 - etaCell2) ) 
							+ 2 * (axis.e() * cosh(etaCell2) - axis.perp() * cos(phiCell2)) );
	    fastjet::PseudoJet p2 = fastjet::PtYPhiM(pTcell2, etaCell2, phiCell2);
	    fastjet::PseudoJet p3 = axis - p1 - p2;
							
	    // Reset Fastjet input
	    fjInputs.resize(0);
		//min pt cut on momenta					
	    if (p1.perp() < minPt) continue;
	    if (p2.perp() < minPt) continue;
	    if (p3.perp() < minPt) continue;
		  
	    if (p1.m() < -0.0001) continue;
	    if (p2.m() < -0.0001) continue;
	    if (p3.m() < -0.0001) continue;
	    
	    if (p1.e() < 0.) continue;
	    if (p2.e() < 0.) continue;
	    if (p3.e() < 0.) continue;
	    
	    if (p1.perp() <p2.perp() || p1.perp() < p3.perp() ) continue;

	    fjInputs.push_back(p1); 
	    fjInputs.push_back(p2); 
	    fjInputs.push_back(p3); 

	// Run Fastjet algorithm
	    std::vector <fastjet::PseudoJet> inclusiveJets, sortedJets;
	    fastjet::ClusterSequence clustSeq(fjInputs, *jetDef);
	    inclusiveJets = clustSeq.inclusive_jets();
	    sortedJets    = sorted_by_pt(inclusiveJets);  
    
	    if (sortedJets.size()>1) continue;
	   File << int(3) <<std::endl;
	
	    File << axis <<std::endl;
	    File << p1 << std::endl;
	    File << p2 << std::endl;
	    File << p3 << std::endl;	
	     ++nJets;
      }}}}}
      std::cout << "Found "<<nJets << std::endl;

      break;
    case TOP:
   for( int iEtaNow1 = 0; iEtaNow1 < nEta; ++iEtaNow1 )	{
    for( int iPhiNow1 = 0; iPhiNow1 < nPhi; ++iPhiNow1 )	{
      for( int iEtaNow2 = 0; iEtaNow2 < nEta; ++iEtaNow2 )	{
	for( int iPhiNow2 = 0; iPhiNow2 < nPhi; ++iPhiNow2 )	{
	  double etaCell1 = (etaMax / nEta) * (2 * iEtaNow1 - 1 - nEta);
	  double phiCell1 = (phiMax / nPhi) * (2 * iPhiNow1 - 1 - nPhi);
	  double etaCell2 = (etaMax / nEta) * (2 * iEtaNow2 - 1 - nEta);
	  double phiCell2 = (phiMax / nPhi) * (2 * iPhiNow2 - 1 - nPhi);
	  fastjet::PseudoJet e1 = fastjet::PtYPhiM(1., etaCell1, phiCell1);
	  fastjet::PseudoJet e2 = fastjet::PtYPhiM(1., etaCell2, phiCell2);
	  double p1P = axis * e1;
	  double p2P = axis * e2;
	  double p12 = e1 * e2;
	  double Delta = 1- (8* mW*mW*p1P*p2P)/((Mass*Mass+mW*mW)*(Mass*Mass+mW*mW) *p12);
	  if (Delta <0.) continue;
	  double pT1 = ( (Mass * Mass + mW * mW) / (4 * p1P ) ) *( 1+ sqrt(Delta ));		
	  double pT2 = mW*mW/(2 * pT1 *p12); 
	  fastjet::PseudoJet p1 = pT1 * e1;
	  fastjet::PseudoJet p2 = pT2 * e2;
	  fastjet::PseudoJet p3 = axis - p1 - p2;
	// Reset Fastjet input
	  fjInputs.resize(0);
	  if (p1.m() < -0.0001) continue;
	  if (p2.m() < -0.0001) continue;
	  if (p3.m() < -0.0001) continue;  
	  if (p1.e() < 0.) continue;
	  if (p2.e() < 0.) continue;
	  if (p3.e() < 0.) continue;
	  if (p1.perp() < minPt) continue;
	  if (p2.perp() < minPt) continue;
	  if (p3.perp() < minPt) continue;
	  fjInputs.push_back(p1); 
	  fjInputs.push_back(p2); 
	  fjInputs.push_back(p3); 
	  // Run Fastjet algorithm
	  std::vector <fastjet::PseudoJet> inclusiveJets, sortedJets;
	  fastjet::ClusterSequence clustSeq(fjInputs, *jetDef);
	  inclusiveJets = clustSeq.inclusive_jets();
	  sortedJets    = sorted_by_pt(inclusiveJets);  
	  if (sortedJets.size()>1) continue;
	  File << int(3) <<std::endl;
	  File << axis <<std::endl;
	  File << p3 << std::endl;
	  File << p2 << std::endl;
	  File << p1 << std::endl;	
	  ++nJets;
	 }}}}
	 std::cout << "Found "<<nJets << std::endl;

	 break;
  default: std::cout << "This does not correspond to a valid mode"<<std::endl;
  break;

}
  
}




#endif // end MATCHING_H
