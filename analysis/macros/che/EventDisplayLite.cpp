// EventDisplayLite.cpp - Simple event display for SFGD detector with strict filtering
// Compatible with the SFGD framework build system

#define THIS_NAME EventDisplayLite
#define OVERRIDE_OPTIONS

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TClonesArray.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TPaveText.h"
#include "TApplication.h"

// Include the framework headers FIRST (before defining any globals)
#include "../../src/tools/global_header.hh"

using namespace std;

// Global parameters for the detector geometry
const int DETECTOR_X = 24;  // X dimension in units
const int DETECTOR_Y = 8;   // Y dimension in units  
const int DETECTOR_Z = 48;  // Z dimension in units

// Strict filtering parameters
const double MIN_CHARGE_THRESHOLD = 10.0;  // Minimum PE per hit
const int MIN_HIGH_CHARGE_HITS = 15;       // Minimum number of hits above threshold

// Main function following framework convention
void EventDisplayLite() {
    
    // Default values
    TString fileIn = "";
    TString fileOut = "event_display.root";
    int evtIni = 0;
    int evtFin = -1;
    int minHits = 15;  // This is the old simple filter
    bool batch = true;
    bool saveImages = true;
    TString imageDir = "./";
    
    // Get command line arguments from the global application
    int argc = gApplication->Argc();
    char** argv = gApplication->Argv();
    
    // Parse command line arguments
    for (int iarg = 0; iarg < argc; iarg++) {
        string arg = argv[iarg];
        
        if (arg == "-h" || arg == "--help") {
            cout << "Simple Event Display for SFGD Detector with Strict Filtering\n";
            cout << "Usage: " << argv[0] << " -i input.root [options]\n";
            cout << "Options:\n";
            cout << "  -i, --input FILE     Input ROOT file (required)\n";
            cout << "  -o, --output FILE    Output ROOT file (default: event_display.root)\n";
            cout << "  -a, --evtIni N       Start event number (default: 0)\n";
            cout << "  -z, --evtFin N       End event number (default: all events)\n";
            cout << "  -m, --minhits N      Minimum hits to display event (default: 15)\n";
            cout << "  -b, --batch          Run in batch mode (no display)\n";
            cout << "  --save-images        Save PNG images of events\n";
            cout << "  --image-dir DIR      Directory for saving images (default: ./)\n";
            cout << "  -h, --help           Show this help message\n";
            cout << "\nStrict Filtering Criteria:\n";
            cout << "  - Minimum charge per hit: " << MIN_CHARGE_THRESHOLD << " PE\n";
            cout << "  - Minimum high-charge hits: " << MIN_HIGH_CHARGE_HITS << " hits\n";
            cout << "  - Only events passing BOTH criteria will generate PNG output\n";
            return;
        }
        else if ((arg == "-i" || arg == "--input") && iarg + 1 < argc) {
            iarg++;
            fileIn = argv[iarg];
        }
        else if ((arg == "-o" || arg == "--output") && iarg + 1 < argc) {
            iarg++;
            fileOut = argv[iarg];
        }
        else if ((arg == "-a" || arg == "--evtIni") && iarg + 1 < argc) {
            iarg++;
            evtIni = atoi(argv[iarg]);
        }
        else if ((arg == "-z" || arg == "--evtFin") && iarg + 1 < argc) {
            iarg++;
            evtFin = atoi(argv[iarg]);
        }
        else if ((arg == "-m" || arg == "--minhits") && iarg + 1 < argc) {
            iarg++;
            minHits = atoi(argv[iarg]);
        }
        else if (arg == "-b" || arg == "--batch") {
            batch = true;
            gROOT->SetBatch(kTRUE);
        }
        else if (arg == "--save-images") {
            saveImages = true;
        }
        else if (arg == "--image-dir" && iarg + 1 < argc) {
            iarg++;
            imageDir = argv[iarg];
            // Ensure directory ends with /
            if (!imageDir.EndsWith("/")) imageDir += "/";
        }
    }
    
    // Check if input file was provided
    if (fileIn == "") {
        cerr << "Error: Input file is required. Use -h for help.\n";
        return;
    }
    
    // Set up output paths based on input file location
    TString inputDir = gSystem->DirName(fileIn);
    TString inputBaseName = gSystem->BaseName(fileIn);
    inputBaseName.ReplaceAll(".root", "");
    
    // If no output file specified, create it in the same directory as input
    if (fileOut == "event_display.root") {
        fileOut = inputDir + "/" + inputBaseName + "_display.root";
    }
    
    // Set up image directory with specific naming for filtered events
    if (saveImages && imageDir == "./") {
        imageDir = inputDir + "/EventDisplay_PNG_Filtered/";
    }
    
    // Create image directory if saving images
    if (saveImages) {
        gSystem->mkdir(imageDir, kTRUE);  // kTRUE = create parent directories if needed
        cout << "PNG files for filtered events will be saved to: " << imageDir << endl;
    }
    
    // Set style
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kRainBow);
    
    cout << "Input file: " << fileIn << endl;
    cout << "Output file: " << fileOut << endl;
    cout << "\n=== STRICT FILTERING CRITERIA ===" << endl;
    cout << "Minimum charge per hit: " << MIN_CHARGE_THRESHOLD << " PE" << endl;
    cout << "Minimum high-charge hits required: " << MIN_HIGH_CHARGE_HITS << " hits" << endl;
    cout << "=================================" << endl;
    
    // Open input file
    TFile* fInput = TFile::Open(fileIn, "READ");
    if (!fInput || fInput->IsZombie()) {
        cerr << "Error: Cannot open input file: " << fileIn << endl;
        return;
    }
    
    // Get the event tree - try both possible tree names
    TTree* eventTree = (TTree*)fInput->Get("TimeGroupedEvents");
    if (!eventTree) {
        eventTree = (TTree*)fInput->Get("AllEvents");
    }
    
    if (!eventTree) {
        cerr << "Error: Cannot find event tree in input file" << endl;
        fInput->Close();
        return;
    }
    
    // Create event object and set branch address
    Event* event = new Event();
    eventTree->SetBranchAddress("Event", &event);
    
    // Open output file
    TFile* fOutput = new TFile(fileOut, "RECREATE");
    if (!fOutput || fOutput->IsZombie()) {
        cerr << "Error: Cannot create output file: " << fileOut << endl;
        fInput->Close();
        return;
    }
    
    // Create histograms
    TH2F* hXY = new TH2F("hXY", "XY View;X [channels];Y [channels];Charge [PE]", 
                         DETECTOR_X, 0, DETECTOR_X, DETECTOR_Y, 0, DETECTOR_Y);
    TH2F* hXZ = new TH2F("hXZ", "XZ View;X [channels];Z [channels];Charge [PE]", 
                         DETECTOR_X, 0, DETECTOR_X, DETECTOR_Z, 0, DETECTOR_Z);
    TH2F* hZY = new TH2F("hZY", "ZY View;Z [channels];Y [channels];Charge [PE]", 
                         DETECTOR_Z, 0, DETECTOR_Z, DETECTOR_Y, 0, DETECTOR_Y);
    
    // Summary histograms
    TH1F* hHitMultiplicity = new TH1F("hHitMultiplicity", 
                                      "Hit Multiplicity;Number of Hits;Events", 
                                      100, 0, 100);
    TH1F* hTotalCharge = new TH1F("hTotalCharge", 
                                   "Total Event Charge;Total Charge [PE];Events", 
                                   200, 0, 2000);
    TH1F* hMaxCharge = new TH1F("hMaxCharge", 
                                "Maximum Hit Charge;Max Charge [PE];Events", 
                                200, 0, 200);
    
    // New histograms for strict filtering analysis
    TH1F* hHighChargeHits = new TH1F("hHighChargeHits", 
                                     Form("High Charge Hits (>%.1f PE);Number of High Charge Hits;Events", MIN_CHARGE_THRESHOLD), 
                                     100, 0, 100);
    TH1F* hChargeDistribution = new TH1F("hChargeDistribution", 
                                         "Hit Charge Distribution;Charge [PE];Hits", 
                                         200, 0, 200);
    TH1F* hFilteredHitMultiplicity = new TH1F("hFilteredHitMultiplicity", 
                                              "Hit Multiplicity (Filtered Events);Number of Hits;Events", 
                                              100, 0, 100);
    
    hXY->SetDirectory(nullptr);
    hXZ->SetDirectory(nullptr);
    hZY->SetDirectory(nullptr);
    hHitMultiplicity->SetDirectory(nullptr);
    hTotalCharge->SetDirectory(nullptr);
    hMaxCharge->SetDirectory(nullptr);
    hHighChargeHits->SetDirectory(nullptr);
    hChargeDistribution->SetDirectory(nullptr);
    hFilteredHitMultiplicity->SetDirectory(nullptr);

    // Get number of events
    Long64_t nEvents = eventTree->GetEntries();
    cout << "Total events in file: " << nEvents << endl;
    
    // Set event range
    if (evtFin < 0 || evtFin > nEvents) evtFin = nEvents;
    
    cout << "Processing events " << evtIni << " to " << evtFin - 1 << endl;
    cout << "Event selection criteria:" << endl;
    cout << "  - Old minimum hits filter: " << minHits << " (for histogram filling)" << endl;
    cout << "  - NEW strict filter: " << MIN_HIGH_CHARGE_HITS << " hits with >" << MIN_CHARGE_THRESHOLD << " PE (for PNG generation)" << endl;
    
    // Create canvas for interactive display
    TCanvas* canvas = new TCanvas("EventDisplay", "SFGD Event Display", 1200, 900);
    canvas->Divide(2, 2);
    
    // Counters
    int nProcessedEvents = 0;
    int nDisplayedEvents = 0;
    int nStrictFilteredEvents = 0;  // Events passing the strict filter
    int nOldFilteredEvents = 0;     // Events passing the old filter
    
    // Event loop
    for (int iEvent = evtIni; iEvent < evtFin; iEvent++) {
        eventTree->GetEntry(iEvent);
        
        if (iEvent % 100 == 0) {
            cout << "Processing event " << iEvent << " (" 
                 << (100.0 * (iEvent - evtIni) / (evtFin - evtIni)) 
                 << "% done) - Old filter: " << nOldFilteredEvents 
                 << ", Strict filter: " << nStrictFilteredEvents << "\r" << flush;
        }
        
        nProcessedEvents++;
        
        // Clear histograms
        hXY->Reset();
        hXZ->Reset();
        hZY->Reset();
        
        // Process hits and apply strict filtering
        int nHits = 0;
        int nHighChargeHits = 0;  // NEW: count hits above charge threshold
        double totalCharge = 0;
        double maxCharge = 0;
        
        TClonesArray* hits = event->GetHits();
        if (!hits) continue;
        
        // First pass: analyze all hits and count high-charge hits
        for (int iHit = 0; iHit < hits->GetEntries(); iHit++) {
            Hit* hit = (Hit*)hits->At(iHit);
            if (!hit) continue;
            
            double charge = hit->GetPE();
            
            nHits++;
            totalCharge += charge;
            if (charge > maxCharge) maxCharge = charge;
            
            // Fill charge distribution for all hits
            hChargeDistribution->Fill(charge);
            
            // NEW: Count hits above charge threshold
            if (charge >= MIN_CHARGE_THRESHOLD) {
                nHighChargeHits++;
            }
        }
        
        // Fill summary histogram for high charge hits (for all events)
        hHighChargeHits->Fill(nHighChargeHits);
        
        // Apply OLD filter (for histogram filling and compatibility)
        bool passesOldFilter = (nHits >= minHits);
        
        // Apply NEW strict filter (for PNG generation)
        bool passesStrictFilter = (nHighChargeHits >= MIN_HIGH_CHARGE_HITS);
        
        if (passesOldFilter) {
            nOldFilteredEvents++;
            
            // Fill summary histograms only for events passing the old filter
            hHitMultiplicity->Fill(nHits);
            hTotalCharge->Fill(totalCharge);
            hMaxCharge->Fill(maxCharge);
        }
        
        if (passesStrictFilter) {
            nStrictFilteredEvents++;
            
            // Fill filtered event histograms
            hFilteredHitMultiplicity->Fill(nHits);
        }
        
        // Only display and save images for events passing STRICT filter
        if (passesStrictFilter) {
            nDisplayedEvents++;
            
            // Second pass: fill histograms for display (only for strict filtered events)
            for (int iHit = 0; iHit < hits->GetEntries(); iHit++) {
                Hit* hit = (Hit*)hits->At(iHit);
                if (!hit) continue;
                
                double charge = hit->GetPE();
                
                // Fill appropriate histogram based on view
                if (hit->GetView() == 0) {  // XY view
                    hXY->Fill(hit->GetX(), hit->GetY(), charge);
                }
                else if (hit->GetView() == 1) {  // XZ view
                    hXZ->Fill(hit->GetX(), hit->GetZ(), charge);
                }
                else if (hit->GetView() == 2) {  // ZY view
                    hZY->Fill(hit->GetZ(), hit->GetY(), charge);
                }
            }
            
            // Display event
            if (canvas) {
                canvas->cd(1);
                hXY->Draw("COLZ");
                gPad->SetRightMargin(0.15);
                
                canvas->cd(2);
                hXZ->Draw("COLZ");
                gPad->SetRightMargin(0.15);
                
                canvas->cd(3);
                hZY->Draw("COLZ");
                gPad->SetRightMargin(0.15);
                
                canvas->cd(4);
                TPaveText* info = new TPaveText(0.1, 0.1, 0.9, 0.9, "NDC");
                info->AddText(Form("Event: %d", iEvent));
                info->AddText(Form("Total hits: %d", nHits));
                info->AddText(Form("High charge hits (>%.1f PE): %d", MIN_CHARGE_THRESHOLD, nHighChargeHits));
                info->AddText(Form("Max charge: %.1f PE", maxCharge));
                info->AddText(Form("Total charge: %.1f PE", totalCharge));
                info->AddText(Form("STRICT FILTERED EVENT #%d", nDisplayedEvents));
                info->AddText("*** PASSES STRICT FILTER ***");
                info->SetTextColor(kRed);
                info->Draw();
                
                canvas->Update();
                
                // Save image (ONLY for events passing strict filter)
                if (saveImages) {
                    TString imageName = Form("%sevent_%06d.png", imageDir.Data(), iEvent);
                    canvas->SaveAs(imageName);
                }
                
                delete info;
            }
            
            // Save event histograms to file (ONLY for strict filtered events)
            fOutput->cd();
            hXY->Write(Form("Event%06d_XY_Filtered", iEvent));
            hXZ->Write(Form("Event%06d_XZ_Filtered", iEvent));
            hZY->Write(Form("Event%06d_ZY_Filtered", iEvent));
        }
    }
    
    cout << "\nProcessing complete!" << endl;
    cout << "======================" << endl;
    cout << "Total events processed: " << nProcessedEvents << endl;
    cout << "Events passing OLD filter (>=" << minHits << " hits): " << nOldFilteredEvents << endl;
    cout << "Events passing STRICT filter (>=" << MIN_HIGH_CHARGE_HITS << " hits with >=" << MIN_CHARGE_THRESHOLD << " PE): " << nStrictFilteredEvents << endl;
    cout << "PNG images generated: " << nStrictFilteredEvents << endl;
    cout << "Strict filter efficiency: " << (100.0 * nStrictFilteredEvents / nProcessedEvents) << "%" << endl;
    cout << "Old filter efficiency: " << (100.0 * nOldFilteredEvents / nProcessedEvents) << "%" << endl;
    
    // Save summary histograms
    fOutput->cd();
    hHitMultiplicity->Write();
    hTotalCharge->Write();
    hMaxCharge->Write();
    hHighChargeHits->Write();
    hChargeDistribution->Write();
    hFilteredHitMultiplicity->Write();
    
    // Save filtering statistics
    TH1F* hFilterStats = new TH1F("hFilterStats", "Filtering Statistics", 4, 0, 4);
    hFilterStats->SetDirectory(nullptr);
    hFilterStats->SetBinContent(1, nProcessedEvents);
    hFilterStats->SetBinContent(2, nOldFilteredEvents);
    hFilterStats->SetBinContent(3, nStrictFilteredEvents);
    hFilterStats->SetBinContent(4, nStrictFilteredEvents); // PNG count
    hFilterStats->GetXaxis()->SetBinLabel(1, "Total Events");
    hFilterStats->GetXaxis()->SetBinLabel(2, "Old Filter");
    hFilterStats->GetXaxis()->SetBinLabel(3, "Strict Filter");
    hFilterStats->GetXaxis()->SetBinLabel(4, "PNG Generated");
    hFilterStats->Write();
    
    // Cleanup
    cout << "Output saved to: " << fileOut << endl;
    cout << "PNG images saved to: " << imageDir << endl;

    // Close files first
    fOutput->Close();
    delete fOutput;

    fInput->Close();
    delete fInput;

    // Delete histograms
    delete hXY;
    delete hXZ;
    delete hZY;
    delete hHitMultiplicity;
    delete hTotalCharge;
    delete hMaxCharge;
    delete hHighChargeHits;
    delete hChargeDistribution;
    delete hFilteredHitMultiplicity;
    delete hFilterStats;

    // Delete canvas
    delete canvas;

    // Delete event object last
    delete event;
    exit(0);

    return;
}