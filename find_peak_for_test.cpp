// This example is developed by Meng Lyu based on PeakFinder class
// GitHub link: https://github.com/claydergc/find-peaks
// Editted the min_v max_v part by Beining Rao on 2024/11/23
//added the peak range function by Beining Rao on 2024/11/27

#include <iostream>
#include <algorithm>
#include "PeakFinder.cpp"
#include "TChain.h"
#include "TH1I.h"
#include "../ChannelReader.h"


void loader(int rate);

void find_peak_for_test()
{
    Double_t offset1 = 150;  // The offset of A B channels (set in picoscope software)  AB is LabR3 channel
    Double_t offset2 = 0;  // The offset of C D channels (set in picoscope software) CD is sipm channel
    Double_t threshold1 = -10;   // The threshold, below which it won't find peak for A B
    Double_t threshold2 = 50;   // The threshold, below which it won't find peak for C D

    TString infpath = "../../../Test1030/20241030/RootData/20241030_merged.root";
    TString outfpath = "20241030_testwithaerogel_peaks.root";


    TChain* tt = new TChain("wfm");
    tt->Add(infpath);

    // Initialize ChannelReader Class AFTER reading the tree and BEFORE creating any instance
    ChannelReader::Initialize(tt);

    // Input parameter is the channel name that determines which channel the reader will read
    ChannelReader ChA1("ChA1"), ChB1("ChB1"), ChC1("ChC1"), ChD1("ChD1"),
            ChA2("ChA2"), ChB2("ChB2"), ChC2("ChC2"), ChD2("ChD2");

    TFile* fo = TFile::Open(outfpath, "RECREATE");
    TTree* peakTree = new TTree("PeakTree", "Tree storing peaks data");

    int eventNumber;
    std::vector<float> peakTimesA, peakValuesA;
    std::vector<float> peakTimesB, peakValuesB;
    std::vector<float> peakTimesC, peakValuesC;
    std::vector<float> peakTimesD, peakValuesD;
    std::vector<float> peak_A1_rangeleft,peak_B1_rangeleft,peak_C1_rangeleft,peak_D1_rangeleft;
    std::vector<float> peak_A1_rangeright,peak_B1_rangeright,peak_C1_rangeright,peak_D1_rangeright;
   
    peakTree->Branch("EventNumber", &eventNumber);
    peakTree->Branch("PeakTimesA", &peakTimesA);
    peakTree->Branch("PeakValuesA", &peakValuesA);
    peakTree->Branch("PeakrangeAleft",&peak_A1_rangeleft);
    peakTree->Branch("PeakrangeAright",&peak_A1_rangeright);
    peakTree->Branch("PeakTimesB", &peakTimesB);
    peakTree->Branch("PeakrangeBleft",&peak_B1_rangeleft);
    peakTree->Branch("PeakrangeBright",&peak_B1_rangeright);
    peakTree->Branch("PeakValuesB", &peakValuesB);
    peakTree->Branch("PeakTimesC", &peakTimesC);
    peakTree->Branch("PeakrangeCleft",&peak_C1_rangeleft);
    peakTree->Branch("PeakrangeCright",&peak_C1_rangeright);
    peakTree->Branch("PeakValuesC", &peakValuesC);
    peakTree->Branch("PeakTimesD", &peakTimesD);
    peakTree->Branch("PeakValuesD", &peakValuesD);
    peakTree->Branch("PeakrangeDleft",&peak_D1_rangeleft);
    peakTree->Branch("PeakrangeDright",&peak_D1_rangeright);

    int nentries = tt->GetEntries();
    std::cout << "Processing " << nentries << " events..." << std::endl;
     
    std::vector<int> peakn_A1, peakn_B1, peakn_C1, peakn_D1,
            peakn_A2, peakn_B2, peakn_C2, peakn_D2;
    std::vector<int> validEventsA; 
    std::vector<int> validEventsB; 
    // std::vector<std::pair<int, int>> peakn_A1_range,peakn_B1_range,peakn_C1_range,peakn_D1_range;
   
    for (int ientry = 0; ientry < nentries; ientry++) {
        if (nentries >= 100 && ientry % (nentries/100) == 0) loader(ientry/(nentries/100));

        tt->GetEntry(ientry);

        peakn_A1.clear();
        peakn_A1.shrink_to_fit();
        peakn_B1.clear();
        peakn_B1.shrink_to_fit();
        
         // PeakFinder only supports float
        peakTimesA.clear();
        peakValuesA.clear();
        peak_A1_rangeleft.clear();
        peak_A1_rangeright.clear();
        // peakn_A1_range.clear();
        peakn_A1.shrink_to_fit();
        // peakn_A1_range.shrink_to_fit();
        peakTimesB.clear();
        peakValuesB.clear();
        peak_B1_rangeleft.clear();
        peak_B1_rangeright.clear();
        // peakn_B1_range.clear();
        peakn_B1.shrink_to_fit();
        // peakn_B1_range.shrink_to_fit();
        std::vector<float> ChA1_vec_float, ChB1_vec_float, ChA2_vec_float, ChB2_vec_float;

        for (int i=0; i<ChA1.V->size(); i++){
            ChA1_vec_float.push_back(ChA1.V->at(i) - offset1);
            ChB1_vec_float.push_back(ChB1.V->at(i) - offset1);
        }
        ChA1.CalculateMinVoltage();
        ChB1.CalculateMinVoltage();
        ChA1.CalculateMaxVoltage();
        ChB1.CalculateMaxVoltage();
        // std::cout << std::endl << "C1max:" << ChC1.max_v <<std::endl;
        // std::cout << std::endl << "D1max:" << ChD1.max_v <<std::endl;
        // std::cout << std::endl << "A1min:" << ChA1.min_v <<std::endl;
        // std::cout << std::endl << "B1min:" << ChB1.min_v <<std::endl;
        if(ChA1.min_v <= threshold1 + offset1) PeakFinder::findPeaks(ChA1_vec_float, peakn_A1,false, -1); // -1 is find minimum, 1 is find maximum
        if(ChB1.min_v <= threshold1 + offset1) PeakFinder::findPeaks(ChB1_vec_float, peakn_B1, false, -1);
        if (peakn_A1.size() == 0 && peakn_B1.size() == 0) continue;

        // ONLY DO NECESSARY THINGS BEFORE VALID DATA FOUND
        
        peakn_C1.clear();
        peakn_C1.shrink_to_fit();
        peakn_D1.clear();
        peakn_D1.shrink_to_fit();
        
        peakTimesC.clear();
        peakValuesC.clear();
        peak_C1_rangeleft.clear();
        peak_C1_rangeright.clear();
        // peakn_C1_range.clear();
        // peakn_C1_range.shrink_to_fit();
        peakTimesD.clear();
        peakValuesD.clear();
        peak_D1_rangeleft.clear();
        peak_D1_rangeright.clear();
        // peakn_D1_range.clear();
        // peakn_D1_range.shrink_to_fit();

        std::vector<float> ChC1_vec_float, ChD1_vec_float, ChC2_vec_float, ChD2_vec_float;

        for (int i=0; i<ChA1.V->size(); i++){
            ChC1_vec_float.push_back(ChC1.V->at(i) - offset2);
            ChD1_vec_float.push_back(ChD1.V->at(i) - offset2);
        }
        ChC1.CalculateMaxVoltage();
        ChD1.CalculateMaxVoltage();
        if(ChC1.max_v >= threshold2 + offset2 && !std::isinf(ChC1.max_v)) PeakFinder::findPeaks(ChC1_vec_float, peakn_C1, false, 1);
        if(ChD1.max_v >= threshold2 + offset2 && !std::isinf(ChD1.max_v)) PeakFinder::findPeaks(ChD1_vec_float, peakn_D1, false, 1);
        if (peakn_C1.size() == 0 && peakn_D1.size() == 0) continue;

        std::cout << std::endl << "Potential Valid Event: " << ientry << std::endl;

        Double_t timeLimit = 999;

        std::cout << "ChA1: ";
        // if (peakn_A1_range.size()==0) continue;

        // size_t counter =0;
   
        for (auto i :peakn_A1) {
            std::cout << "\t" << i;
            std::cout << "\t" << ChA1.T->at(i);
            peakTimesA.push_back(ChA1.T->at(i));
            peakValuesA.push_back(ChA1.V->at(i));
            // std::cout << "peak_A1_range" << peakn_A1_range.size() << std::endl;
       
            
            int leftBoundary = i;
            while (leftBoundary >0 && std::abs(ChA1_vec_float[leftBoundary]) >= 0.1) {
                leftBoundary--;
            }

            int rightBoundary = i;
            while (rightBoundary < ChA1_vec_float.size()-1 && std::abs(ChA1_vec_float[rightBoundary]) >= 0.1) {
                rightBoundary++;
            }

            
        

            float leftBoundaryTime = ChA1.T->at(leftBoundary);
            float rightBoundaryTime = ChA1.T->at(rightBoundary);
            std::cout << "Left Boundary: " << leftBoundaryTime << ", Right Boundary: " << rightBoundaryTime << std::endl;
            // if (peakn_A1_range[i].first >= 0 && peakn_A1_range[i].first < ChA1.T->size() &&
            //     peakn_A1_range[i].second >= 0 && peakn_A1_range[i].second < ChA1.T->size()) {
            //     int leftBoundaryIndex = ChA1.T->at(peakn_A1_range[i].first);
            //     int rightBoundaryIndex = ChA1.T->at(peakn_A1_range[i].second);
            peak_A1_rangeleft.push_back(leftBoundaryTime);
            peak_A1_rangeright.push_back(rightBoundaryTime);
        
            // } else {
            //     std::cerr << "Warning: Index out of range while accessing peakn_A1_range at index " << i << std::endl;
            // }
            // counter++;
            
            timeLimit = ChA1.T->at(i);
            break;  
            }
        
        std::cout << std::endl;


        std::cout << "ChB1: ";
       
        for (auto i : peakn_B1) {
            std::cout << "\t" << i;
            std::cout << "\t" << ChB1.T->at(i);
            peakTimesB.push_back(ChB1.T->at(i));
            peakValuesB.push_back(ChB1.V->at(i));
 
            int leftBoundary = i;
            while (leftBoundary > 0 && std::abs(ChB1_vec_float[leftBoundary]) >= 0.1) {
                leftBoundary--;
            }

            int rightBoundary = i;
            while (rightBoundary <ChB1_vec_float.size()-1 && std::abs(ChB1_vec_float[rightBoundary]) >= 0.1) {
                rightBoundary++;
            }

            
        

            float leftBoundaryTime = ChB1.T->at(leftBoundary);
            float rightBoundaryTime = ChB1.T->at(rightBoundary);
            std::cout << "Left Boundary: " << leftBoundaryTime << ", Right Boundary: " << rightBoundaryTime << std::endl;
            // if (peakn_A1_range[i].first >= 0 && peakn_A1_range[i].first < ChA1.T->size() &&
            //     peakn_A1_range[i].second >= 0 && peakn_A1_range[i].second < ChA1.T->size()) {
            //     int leftBoundaryIndex = ChA1.T->at(peakn_A1_range[i].first);
            //     int rightBoundaryIndex = ChA1.T->at(peakn_A1_range[i].second);
            peak_B1_rangeleft.push_back(leftBoundaryTime);
            peak_B1_rangeright.push_back(rightBoundaryTime);
            
            // if (i >= 0 && i < peakn_B1_range.size()) {
            //     if (peakn_B1_range[i].first >= 0 && peakn_B1_range[i].first < ChB1.T->size() &&
            //         peakn_B1_range[i].second >= 0 && peakn_B1_range[i].second < ChB1.T->size()) {
            //         int leftBoundaryIndex = ChB1.T->at(peakn_B1_range[i].first);
            //         int rightBoundaryIndex = ChB1.T->at(peakn_B1_range[i].second);
            //         peak_B1_range.push_back({leftBoundaryIndex, rightBoundaryIndex});
            //     } else {
            //         std::cerr << "Warning: Boundary indices are out of range for peakn_B1_range at index " << i << std::endl;
            //     }
            // } else {
            //     std::cerr << "Warning: Index out of range while accessing peakn_B1_range at index " << i << std::endl;
            // }

            // update timeLimit
            timeLimit = ChB1.T->at(i) < timeLimit ? ChB1.T->at(i) : timeLimit;

            break;  // stop at the first peak
            }
        std::cout << std::endl;

        std::cout << "ChC1: ";
        for (auto i : peakn_C1) {
            if (ChC1.T->at(i) > timeLimit || ChC1.V->at(i) < threshold2) continue;

            std::cout << "\t" << ChC1.T->at(i);
            std::cout << "\t" << ChC1.V->at(i);
            peakTimesC.push_back(ChC1.T->at(i));
            peakValuesC.push_back(ChC1.V->at(i));
            int leftBoundary = i;
            while (leftBoundary > 0 && std::abs(ChC1_vec_float[leftBoundary]) >= 0.1) {
                leftBoundary--;
            }

            int rightBoundary = i;
            while (rightBoundary <ChC1_vec_float.size()-1 && std::abs(ChC1_vec_float[rightBoundary]) >= 0.1) {
                rightBoundary++;
            }

            
        

            float leftBoundaryTime = ChC1.T->at(leftBoundary);
            float rightBoundaryTime = ChC1.T->at(rightBoundary);
            std::cout << "Left Boundary: " << leftBoundaryTime << ", Right Boundary: " << rightBoundaryTime << std::endl;
            // if (peakn_A1_range[i].first >= 0 && peakn_A1_range[i].first < ChA1.T->size() &&
            //     peakn_A1_range[i].second >= 0 && peakn_A1_range[i].second < ChA1.T->size()) {
            //     int leftBoundaryIndex = ChA1.T->at(peakn_A1_range[i].first);
            //     int rightBoundaryIndex = ChA1.T->at(peakn_A1_range[i].second);
            peak_C1_rangeleft.push_back(leftBoundaryTime);
            peak_C1_rangeright.push_back(rightBoundaryTime);
            break;
        }
        std::cout << std::endl;

        std::cout << "ChD1: ";
        for (auto i : peakn_D1) {
            if (ChD1.T->at(i) > timeLimit || ChD1.V->at(i) < threshold2) continue;

            std::cout << "\t" << ChD1.T->at(i);
            std::cout << "\t" << ChD1.V->at(i);
            peakTimesD.push_back(ChD1.T->at(i));
            peakValuesD.push_back(ChD1.V->at(i));
            int leftBoundary = i;
            while (leftBoundary > 0 && std::abs(ChD1_vec_float[leftBoundary]) >= 0.1) {
                leftBoundary--;
            }

            int rightBoundary = i;
            while (rightBoundary <ChD1_vec_float.size()-1 && std::abs(ChD1_vec_float[rightBoundary]) >= 0.1) {
                rightBoundary++;
            }

            
        

            float leftBoundaryTime = ChD1.T->at(leftBoundary);
            float rightBoundaryTime = ChD1.T->at(rightBoundary);
            std::cout << "Left Boundary: " << leftBoundaryTime << ", Right Boundary: " << rightBoundaryTime << std::endl;
            // if (peakn_A1_range[i].first >= 0 && peakn_A1_range[i].first < ChA1.T->size() &&
            //     peakn_A1_range[i].second >= 0 && peakn_A1_range[i].second < ChA1.T->size()) {
            //     int leftBoundaryIndex = ChA1.T->at(peakn_A1_range[i].first);
            //     int rightBoundaryIndex = ChA1.T->at(peakn_A1_range[i].second);
            peak_D1_rangeleft.push_back(leftBoundaryTime);
            peak_D1_rangeright.push_back(rightBoundaryTime);
            break;
        }
        std::cout << std::endl;

        if (peakTimesC.size() == 0 && peakTimesD.size() == 0) continue;  // Skip recording if valid signal CD doesn't exist

        if (peakn_A1.size() != 0) {
            validEventsA.push_back(ientry); // if condition is satisfied , store the entry
        }
        if (peakn_B1.size() != 0) {
            validEventsB.push_back(ientry);
        }

        std::cout << std::endl;
        if (std::find(validEventsA.begin(), validEventsA.end(), ientry) != validEventsA.end() ||
            std::find(validEventsB.begin(), validEventsB.end(), ientry) != validEventsB.end()) {
            eventNumber = ientry;
            peakTree->Fill();
        }
    }
    std::cout << "Events with peaks in ChA1:" << std::endl;
    for (auto event : validEventsA) {
    std::cout << event << " ";
    }
    std::cout << std::endl;
    std::cout << "Events with peaks in ChB1:" << std::endl;
    for (auto event : validEventsB) {
    std::cout << event << " ";
    }

	fo->Write();
	fo->Close();
	std::cout << "\nSave file as: " << outfpath << std::endl;

}

void loader(int rate)
{
    char proc[22];
    memset(proc, '\0', sizeof(proc));

    for (int i = 0; i < rate/5; i++)
    {
        proc[i] = '#';
    }

    printf("\r[%-20s] [%d%%]", proc, rate);        
    fflush(stdout);                                 
}
