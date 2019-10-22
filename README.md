# TemplateTagger
A Template Matching Tool for Jet Substructure

The TemplateTagger Package
Version 1.0.3 
Questions/Comments?  juknevich@gmail.com


Files included:
   README:  This file
   matching.hh:  Core code for template matching
   TemplateTagger.hh:  A simple FunctionOfPseudoJet<templ_t> interface to measure TemplateOverlap
   short_example.cc:  Example program
   example.cc:
   example7.cc:  A longer example program
   Makefile:  Basic makefile to compile example program
   WHsample.dat:  sample WH(->bb) events in text format for example program
   template2bf.dat: sample H->bb templates 

Short Example Program:

   The example program assumes that you have FastJet 3+ installed.
   You can run the sample program:

   make short
   ./short_example template2bf.dat

Long Example Program:

   Works similarly, using a full Pythia data sample and TemplateTagger
   class for overlap measurement.  FastJet 3+ is required.

   make ex3
   ./example5 WHsample.dat 
