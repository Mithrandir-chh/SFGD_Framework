// ROOT libraries
#include "TROOT.h"
#include "TRint.h"
#include "TApplication.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TNtuple.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TObject.h"
#include "TRandom3.h"
#include <TStyle.h>
#include <TTree.h>
#include <TBranch.h>
#include <TObject.h>
#include <TClonesArray.h>
#include <TClass.h>
#include <TChain.h>
#include <TLegend.h>
#include "TDirectory.h"
#include "TGraph.h"



//C++ libraries
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include <array>
#include <time.h>
#include <cmath>
#include <iterator>
#include <algorithm>
#include <sstream>

//Simplify things using std namespace
using namespace std;


//Classes
#include "../../analysis/src/classes/Event.cc"
#include "../../analysis/src/classes/Hit.cc"
