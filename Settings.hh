// This is a modified version of Pythia 8 Setting class
#ifndef Settings_H
#define Settings_H

#include <iostream>
#include <map>
#include <utility> 
#include <string>
#include <sstream> 
#include <fstream>

// Allow string and character manipulation.
#include <cctype>

using namespace std;

 
class Settings {

public:
// Constructor.
  Settings() :isInit(false){init();}
  bool init() ;
// Read in one update from a single line.
  bool readFile(string fileName) ;
  bool readFile(istream& is) ;
  bool readString(string line) ; 

  // Query existence of an entry.
  bool isFlag(string keyIn) {
    return (flags.find(toLower(keyIn)) != flags.end()); }
  bool isMode(string keyIn) { 
    return (modes.find(toLower(keyIn)) != modes.end()); }
  bool isParm(string keyIn) {
    return (params.find(toLower(keyIn)) != params.end()); }
  bool isWord(string keyIn) {
    return (words.find(toLower(keyIn)) != words.end()); }
 
  // Add new entry.
  void addFlag(string keyIn, bool defaultIn) {flags[toLower(keyIn)] = defaultIn; }  
  void addMode(string keyIn, int defaultIn) { modes[toLower(keyIn)] = defaultIn; }      
  void addParm(string keyIn, double defaultIn) { params[toLower(keyIn)] = defaultIn; }  
  void addWord(string keyIn, string defaultIn) {words[toLower(keyIn)] =  defaultIn; }  

  // Give back current value, with check that key exists. 
  bool   flag(string keyIn) {  if (isFlag(keyIn)) return flags[toLower(keyIn)]; return 0; }
  int   mode(string keyIn) {  if (isMode(keyIn)) return modes[toLower(keyIn)]; return 0;}
  double   parm(string keyIn) {  if (isParm(keyIn)) return params[toLower(keyIn)]; return 0.;}
  string  word(string keyIn) {  if (isWord(keyIn)) return words[toLower(keyIn)];  return " ";}

  // Change current value, respecting limits.
  void flag(string keyIn, bool nowIn); 
  void mode(string keyIn, int nowIn);
  void parm(string keyIn, double nowIn); 
  void word(string keyIn, string nowIn); 

private:

  // Maps
  map<string, bool> flags;
  map<string, int> modes;
  map<string, double> params;
  map<string, string> words;

  // Flag that initialization has been performed.
  bool isInit;

  string toLower(const string& name);
};

bool Settings::init() {
// I/O files
  words["eventfile"] =  "WH.txt";
  words["pileupfile"] = "pileups.txt";
  words["file2b1"]  =       "template2bf.txt";
  words["file3b1"]  =       "template3bf.txt";
  words["file2b2"]  =       "template2bf.txt";
  words["file3b2"]  =       "template3bf.txt";
  words["file2b3"]  =       "template2bf.txt";	
  words["file3b3"]  =       "template3bf.txt";
  words["file2b4"]  =       "template2bf.txt";
  words["file3b4"]  =       "template3bf.txt";
  words["file2b5"]  =       "template2bf.txt";
  words["file3b5"]  =       "template3bf.txt";
  words["file2b6"]  =       "template2bf.txt";
  words["file3b6"]  =       "template3bf.txt";
  words["file2b7"]  =       "template2bf.txt";
  words["file3b7"]  =       "template3bf.txt";
  words["outputfile"] =      "overlap.txt";
  words["jetfile"] =      "jets.txt";
  words["minioutputfile"] =	   "output.txt";
  
  modes["nptbins"] = 7;
  params["ptbinwidth"] = 30.;

  //overlap
  params["subconeradius2"] =      0.4;
  params["subconeradius3"] =      0.2  ;
  params["subconecomparison2"] =    0.4; //  ! 0. =dissabled
  flags["mode"] = 0 ;  //   ! 0: Step function 1: Gaussian
  params["sigma"] = 0.333;

  //pileup overlay
  params["npileupavg"] =   8.8;
  params["minptclus"]  =   0.8;
  params["ov2min"] = 0.;

  //fastjet
  params["coneradius"] =     1.4;
  params["coneradius2"] =     0.4;
  params["minpt"]=          200.;
  params["maxpt"] =       7000.;
  params["minmass"] =      0.;
  params["maxmass"] =     7000.;
  modes["nmaxjets"] =     1000;
  modes["btags"] = 0; 
  modes["njetoutside"] = 10 ; 
  params["minptw"] = 150.;
  params["minmet"] = 40.;	
  params["ptjetmin"] =  15.;
  isInit = true;

  flags["debug"] = false;

return true;
}

bool Settings::readFile(string fileName) {

  // Open file for reading.
  const char* cstring = fileName.c_str();
  ifstream is(cstring);  
  if (!is.good()) {return false;}
  
  
  // Hand over real work to next method.
  return readFile( is);

}

//read in setting file.	
bool Settings::readFile(istream& is) {

  // Read in one line at a time.
  string line;
  while ( getline(is, line) ) {

    readString(line) ;
  
  }
  return true;

}

// Read in updates from a character string, like a line of a file. 

bool Settings::readString(string line) {

  // If empty line then done.
  if (line.find_first_not_of(" \n\t\v\b\r\f\a") == string::npos) return true;

  // If first character is not a letter, then taken to be a comment line.
  string lineNow = line;
  int firstChar = lineNow.find_first_not_of(" \n\t\v\b\r\f\a");
  if (!isalpha(lineNow[firstChar])) return true; 

  // Replace an equal sign by a blank to make parsing simpler.
  while (lineNow.find("=") != string::npos) {
    int firstEqual = lineNow.find_first_of("=");
    lineNow.replace(firstEqual, 1, " ");   
  }

  // Get first word of a line.
  istringstream splitLine(lineNow);
  string name;
  splitLine >> name;
     
  // Check whether this is in the database. Done if not.
  int inDataBase = 0;
  if      (isFlag(name)) inDataBase = 1;   
  else if (isMode(name)) inDataBase = 2;   
  else if (isParm(name)) inDataBase = 3; 
  else if (isWord(name)) inDataBase = 4; 
  if (inDataBase == 0) return false;  
    

  // Find value. Warn if none found.
  string valueString;
  splitLine >> valueString;
  if (!splitLine) {return false;}
    

  // Update flag map; allow many ways to say yes.
  if (inDataBase == 1) {
   istringstream flagData(valueString);

    bool value ;
	     flagData >> value;
    if (!flagData) return false;  
	  
    flag(name, value);

  // Update mode map.
  } else if (inDataBase == 2) {
    istringstream modeData(valueString);
    int value;
    modeData >> value;
    if (!modeData) return false;  
    
    mode(name, value);
        
  // Update parm map.
  } else if (inDataBase == 3) {
    istringstream parmData(valueString);
    double value;
    parmData >> value;
    if (!parmData) return false;  
      
    parm(name, value);
        
  // Update word map.
  } else {
    word(name, valueString);
  }

  // Done.
  return true;
}

// Change current value.

void Settings::flag(string keyIn, bool nowIn) { 
    if (isFlag(keyIn)) flags[toLower(keyIn)] = nowIn; 
}

void Settings::mode(string keyIn, int nowIn) { 
  if (isMode(keyIn)) modes[toLower(keyIn)] = nowIn; 
}

void Settings::parm(string keyIn, double nowIn) { 
  if (isParm(keyIn)) params[toLower(keyIn)] = nowIn;
}  

void Settings::word(string keyIn, string nowIn) { 
    if (isWord(keyIn)) words[toLower(keyIn)] = nowIn; 
}

// Convert string to lowercase for case-insensitive comparisons.
// Also remove initial and trailing blanks, if any.

string Settings::toLower(const string& name) { 

  // Copy string without initial and trailing blanks.
  if (name.find_first_not_of(" \n\t\v\b\r\f\a") == string::npos) return "";
  int firstChar = name.find_first_not_of(" \n\t\v\b\r\f\a");
  int lastChar  = name.find_last_not_of(" \n\t\v\b\r\f\a");
  string temp   = name.substr( firstChar, lastChar + 1 - firstChar);

  // Convert to lowercase letter by letter.
  for (int i = 0; i < int(temp.length()); ++i) 
    temp[i] = std::tolower(temp[i]); 
  return temp; 

}


#endif // Settings_H
