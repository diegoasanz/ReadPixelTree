//============================================================================
// Name        : ReadPixelTree.cpp
// Author      : Diego Alejandro
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
using namespace std;
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include "TString.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2C.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TMacro.h"
#include "TF1.h"
#include "TNtuple.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLinearFitter.h"
#include "TSpectrum.h"
#include "TPolyMarker.h"
#include "TVirtualFFT.h"
#include "TRandom3.h"
#include "TRandom.h"
#include "TChain.h"

void Prueba(char *rootFilePath, char *outputPath);

TChain * ReadPixRootFile(const char *fileName);

void PhHistogramExtraction(int numCLROC, vector<float> *phROC, vector<float> *clust_ROC_X, vector<float> *clust_ROC_Y, TH1F *phDUT, TH2F *avPHDUT, TH2F *hitMapPHDUT);

void PlotHistogram1D(TH1F *histogram, Bool_t logY, Bool_t logX, const char *title, const char *outputPath, const char *suffix);

void PlotHistogram2D(TH2F *histogram, const char *title, const char *outputPath, const char *suffix);

void FindAverageHistogramDUT(TH2F *histToAvDUT1, TH2F *histHitsDUT1, TH2F *histToAvDUT2, TH2F *histHitsDUT2, TH2F *histToAvDUT3, TH2F *histHitsDUT3, Int_t xbins, Int_t ybins);

void FindAverageHistogramTPlanes(TH2F *histToAvTPlane0, TH2F *histHitsTPlane0, TH2F *histToAvTPlane1, TH2F *histHitsTPlane1, TH2F *histToAvTPlane2, TH2F *histHitsTPlane2, TH2F *histToAvTPlane3, TH2F *histHitsTPlane3, Int_t xbins, Int_t ybins);

void AverageBinHistogram(TH2F *histToAv, TH2F *histHits, Int_t binx, Int_t biny);

const char *rootFilePathComplete = "~/Data/psi_2015_10/root/pixel/TrackedRun312.root";
const char *rootFilePathSmall = "~/DataSmall/psi_2015_10/root/pixel/TrackedRun312.root";

const char *outputPathComplete = "~/Dropbox/201601/DiamondPixels/ReadPixelTree/";
const char *outputPathSmall = "~/Dropbox/201601/DiamondPixels/ReadPixelTree/Small/";

const char *histogramOutputFile = "Histograms.root";

Int_t contour2D = 1024;
Float_t xminR1 = -6.075, xmaxR1 = 6.075, yminR1 = -6.05, ymaxR1 = 6.05;
Int_t divXR1 = 81, divYR1 = 121;
Float_t deltaXR1 = (Float_t)((xmaxR1-xminR1)/divXR1), deltaYR1 = (Float_t)((ymaxR1-yminR1)/divYR1);

int main() {
	gROOT->ProcessLine("#include <vector>");
	gSystem->Load("libMathCore");
	gSystem->Load("libPhysics");
	const char *rootFilePath = rootFilePathComplete;
	const char *outputPath = outputPathComplete;
	Prueba((char *)rootFilePath, (char *)outputPath);
	return 0;
}

TChain * ReadPixRootFile(const char *fileName){
	TChain *chain1 = new TChain("tree");
	chain1->Add(fileName);
	return chain1;
}

void Prueba(char *rootFilePath,char *outputPath){
	//Create Chain to read tree
	TChain *chain1 = ReadPixRootFile(rootFilePath);
	//Make pointers to branches
	TBranch *braPlane = chain1->GetBranch("plane");
	//Create variables to read the tree
	vector<int> *plane=0, *col=0, *row=0, *adc=0;
	vector<unsigned char> *clust_per_plane=0;
	vector<float> *clust_ROC0_X = 0, *clust_ROC0_Y = 0, *ph_ROC0_1cl = 0, *clust_ROC1_X = 0, *clust_ROC1_Y = 0, *ph_ROC1_1cl = 0;
	vector<float> *clust_ROC2_X = 0, *clust_ROC2_Y = 0, *ph_ROC2_1cl = 0, *clust_ROC3_X = 0, *clust_ROC3_Y = 0, *ph_ROC3_1cl = 0;
	vector<float> *clust_ROC4_X = 0, *clust_ROC4_Y = 0, *ph_ROC4_1cl = 0, *clust_ROC5_X = 0, *clust_ROC5_Y = 0, *ph_ROC5_1cl = 0, *clust_ROC6_X = 0, *clust_ROC6_Y = 0, *ph_ROC6_1cl = 0;
	Float_t time=0, dia1X=0, dia1Y=0, dia2X=0, dia2Y=0, chi2_tracks=0, chi2_x=0, chi2_y=0, slope_x=0, slope_y=0;
	Int_t event_number=0;
	//Int_t n_tracks = 0, n_cluster = 0;
	//dimension variables
	vector<int> numCLROC;
	numCLROC.resize(7);
	Int_t numPlane = 0, numCol = 0, numRow = 0, numADC = 0;
	Long64_t numEntries = braPlane->GetEntries();

	chain1->SetBranchAddress("event_number",&event_number);
	chain1->SetBranchAddress("time",&time);
	chain1->SetBranchAddress("plane",&plane);
	chain1->SetBranchAddress("col",&col);
	chain1->SetBranchAddress("row",&row);
	chain1->SetBranchAddress("adc",&adc);
	chain1->SetBranchAddress("diam1_track_x",&dia1X);
	chain1->SetBranchAddress("diam1_track_y",&dia1Y);
	chain1->SetBranchAddress("diam2_track_x",&dia2X);
	chain1->SetBranchAddress("diam2_track_y",&dia2Y);
	chain1->SetBranchAddress("chi2_tracks",&chi2_tracks);
	chain1->SetBranchAddress("chi2_x",&chi2_x);
	chain1->SetBranchAddress("chi2_y",&chi2_y);
	chain1->SetBranchAddress("slope_x",&slope_x);
	chain1->SetBranchAddress("slope_y",&slope_y);
	//chain1->SetBranchAddress("n_tracks",&n_tracks);
	//chain1->SetBranchAddress("n_cluster",&n_cluster);
	chain1->SetBranchAddress("clusters_per_plane",&clust_per_plane);
	chain1->SetBranchAddress("cluster_pos_ROC0_X",&clust_ROC0_X);
	chain1->SetBranchAddress("cluster_pos_ROC0_Y",&clust_ROC0_Y);
	chain1->SetBranchAddress("pulse_height_ROC0_1_cluster",&ph_ROC0_1cl);
	chain1->SetBranchAddress("cluster_pos_ROC1_X",&clust_ROC1_X);
	chain1->SetBranchAddress("cluster_pos_ROC1_Y",&clust_ROC1_Y);
	chain1->SetBranchAddress("pulse_height_ROC1_1_cluster",&ph_ROC1_1cl);
	chain1->SetBranchAddress("cluster_pos_ROC2_X",&clust_ROC2_X);
	chain1->SetBranchAddress("cluster_pos_ROC2_Y",&clust_ROC2_Y);
	chain1->SetBranchAddress("pulse_height_ROC2_1_cluster",&ph_ROC2_1cl);
	chain1->SetBranchAddress("cluster_pos_ROC3_X",&clust_ROC3_X);
	chain1->SetBranchAddress("cluster_pos_ROC3_Y",&clust_ROC3_Y);
	chain1->SetBranchAddress("pulse_height_ROC3_1_cluster",&ph_ROC3_1cl);
	chain1->SetBranchAddress("cluster_pos_ROC4_X",&clust_ROC4_X);
	chain1->SetBranchAddress("cluster_pos_ROC4_Y",&clust_ROC4_Y);
	chain1->SetBranchAddress("pulse_height_ROC4_1_cluster",&ph_ROC4_1cl);
	chain1->SetBranchAddress("cluster_pos_ROC5_X",&clust_ROC5_X);
	chain1->SetBranchAddress("cluster_pos_ROC5_Y",&clust_ROC5_Y);
	chain1->SetBranchAddress("pulse_height_ROC5_1_cluster",&ph_ROC5_1cl);
	chain1->SetBranchAddress("cluster_pos_ROC6_X",&clust_ROC6_X);
	chain1->SetBranchAddress("cluster_pos_ROC6_Y",&clust_ROC6_Y);
	chain1->SetBranchAddress("pulse_height_ROC6_1_cluster",&ph_ROC6_1cl);

	TH2F *hitPlane0 = new TH2F("hitPlane0","Hit map Plane 0",52,-0.5,51.5,80,-0.5,79.5);
	TH2F *hitPlane1 = new TH2F("hitPlane1","Hit map Plane 1",52,-0.5,51.5,80,-0.5,79.5);
	TH2F *hitPlane2 = new TH2F("hitPlane2","Hit map Plane 2",52,-0.5,51.5,80,-0.5,79.5);
	TH2F *hitPlane3 = new TH2F("hitPlane3","Hit map Plane 3",52,-0.5,51.5,80,-0.5,79.5);
	TH2F *hitPlane4 = new TH2F("hitPlane4","Hit map Plane 4",52,-0.5,51.5,80,-0.5,79.5);
	TH2F *hitPlane5 = new TH2F("hitPlane5","Hit map Plane 5",52,-0.5,51.5,80,-0.5,79.5);
	TH2F *hitPlane6 = new TH2F("hitPlane6","Hit map Plane 6",52,-0.5,51.5,80,-0.5,79.5);
	TH2F *avADCDia1 = new TH2F("avADCDia1","Average ADC Diamond 1",52,-0.5,51.5,80,-0.5,79.5);
	TH2F *avADCDia2 = new TH2F("avADCDia2","Average ADC Diamond 2",52,-0.5,51.5,80,-0.5,79.5);
	TH2F *avADCSi = new TH2F("avADCSi","Average ADC Si",52,-0.5,51.5,80,-0.5,79.5);
	TH1F *adcDia1 = new TH1F("adcDia1","ADC Diamond 1",541,-270.5,270.5);
	TH1F *adcDia2 = new TH1F("adcDia2","ADC Diamond 2",541,-270.5,270.5);
	TH1F *adcSi = new TH1F("adcSi","ADC Si",541,-270.5,270.5);

	TH2F *avPHDia1 = new TH2F("avPHDia1","Average PH Diamond 1",divXR1,xminR1,xmaxR1,divYR1,yminR1,ymaxR1);
	TH1F *phDia1 = new TH1F("phDia1","Pulse Height Diamond 1",501,-250,250250);
	TH2F *avPHDia2 = new TH2F("avPHDia2","Average PH Diamond 2",divXR1,xminR1,xmaxR1,divYR1,yminR1,ymaxR1);
	TH1F *phDia2 = new TH1F("phDia2","Pulse Height Diamond 2",501,-2052,50052);
	TH2F *avPHSi = new TH2F("avPHSi","Average PH Silicon",divXR1,xminR1,xmaxR1,divYR1,yminR1,ymaxR1);
	TH1F *phSi = new TH1F("phSi","Pulse Height Silicon",501,-10210,200210);

	TH2F *avPHPlane0 = new TH2F("avPHPlane0","Average PH Plane 0",divXR1,xminR1,xmaxR1,divYR1,yminR1,ymaxR1);
	TH1F *phPlane0 = new TH1F("phPlane0","Pulse Height Plane 0",501,-250,250250);
	TH2F *avPHPlane1 = new TH2F("avPHPlane1","Average PH Plane 1",divXR1,xminR1,xmaxR1,divYR1,yminR1,ymaxR1);
	TH1F *phPlane1 = new TH1F("phPlane1","Pulse Height Plane 1",501,-250,250250);
	TH2F *avPHPlane2 = new TH2F("avPHPlane2","Average PH Plane 2",divXR1,xminR1,xmaxR1,divYR1,yminR1,ymaxR1);
	TH1F *phPlane2 = new TH1F("phPlane2","Pulse Height Plane 2",501,-250,250250);
	TH2F *avPHPlane3 = new TH2F("avPHPlane3","Average PH Plane 3",divXR1,xminR1,xmaxR1,divYR1,yminR1,ymaxR1);
	TH1F *phPlane3 = new TH1F("phPlane3","Pulse Height Plane 3",501,-250,250250);

	TH2F *hitMapPHDia1 = new TH2F("hitMapPHDia1","Hit Map Diamond 1 for PH",divXR1,xminR1,xmaxR1,divYR1,yminR1,ymaxR1);
	TH2F *hitMapPHDia2 = new TH2F("hitMapPHDia2","Hit Map Diamond 2 for PH",divXR1,xminR1,xmaxR1,divYR1,yminR1,ymaxR1);
	TH2F *hitMapPHSi = new TH2F("hitMapPHSi","Hit Map Si for PH",divXR1,xminR1,xmaxR1,divYR1,yminR1,ymaxR1);

	TH2F *hitMapPHPlane0 = new TH2F("hitMapPHPlane0","Hit Map Plane 0 for PH",divXR1,xminR1,xmaxR1,divYR1,yminR1,ymaxR1);
	TH2F *hitMapPHPlane1 = new TH2F("hitMapPHPlane1","Hit Map Plane 1 for PH",divXR1,xminR1,xmaxR1,divYR1,yminR1,ymaxR1);
	TH2F *hitMapPHPlane2 = new TH2F("hitMapPHPlane2","Hit Map Plane 2 for PH",divXR1,xminR1,xmaxR1,divYR1,yminR1,ymaxR1);
	TH2F *hitMapPHPlane3 = new TH2F("hitMapPHPlane3","Hit Map Plane 3 for PH",divXR1,xminR1,xmaxR1,divYR1,yminR1,ymaxR1);

	Double_t tempADC = 0;
	avADCDia1->SetContour(contour2D);
	avADCDia2->SetContour(contour2D);
	avADCSi->SetContour(contour2D);
	hitPlane0->SetContour(contour2D);
	hitPlane1->SetContour(contour2D);
	hitPlane2->SetContour(contour2D);
	hitPlane3->SetContour(contour2D);
	hitPlane4->SetContour(contour2D);
	hitPlane5->SetContour(contour2D);
	hitPlane6->SetContour(contour2D);

	avPHDia1->SetContour(contour2D);
	avPHDia2->SetContour(contour2D);
	avPHSi->SetContour(contour2D);

	avPHPlane0->SetContour(contour2D);
	avPHPlane1->SetContour(contour2D);
	avPHPlane2->SetContour(contour2D);
	avPHPlane3->SetContour(contour2D);

	avADCDia1->SetMinimum(0);
	avADCDia2->SetMinimum(0);
	avADCSi->SetMinimum(0);
	hitPlane0->SetMinimum(0);
	hitPlane1->SetMinimum(0);
	hitPlane2->SetMinimum(0);
	hitPlane3->SetMinimum(0);
	hitPlane4->SetMinimum(0);
	hitPlane5->SetMinimum(0);
	hitPlane6->SetMinimum(0);

	avPHDia1->SetMinimum(0);
	avPHDia2->SetMinimum(0);
	avPHSi->SetMinimum(0);

	avPHPlane0->SetMinimum(0);
	avPHPlane1->SetMinimum(0);
	avPHPlane2->SetMinimum(0);
	avPHPlane3->SetMinimum(0);

	TH1F *colP3 = new TH1F("colP3","colP3",61,-0.5,60.5);
	TH1F *rowP3 = new TH1F("rowP3","rowP3",91,-0.5,90.5);

	for(Long64_t i = 0; i < numEntries; i++){
		chain1->GetEntry(i);
		numPlane = plane->size();
		numCol = col->size();
		numRow = row->size();
		numADC = adc->size();

		if(numPlane & numCol & numRow == numADC){
			for(Int_t j = 0; j < numPlane; j++){
				colP3->Fill(col->at(j));
				rowP3->Fill(row->at(j));
				if(plane->at(j)==0){
					hitPlane0->Fill(col->at(j),row->at(j));
				}
				else if(plane->at(j) == 1){
					hitPlane1->Fill(col->at(j),row->at(j));
				}
				else if(plane->at(j) == 2){
					hitPlane2->Fill(col->at(j),row->at(j));
				}
				else if(plane->at(j) == 3){
					hitPlane3->Fill(col->at(j),row->at(j));
				}
				else if(plane->at(j) == 4){
					hitPlane4->Fill(col->at(j),row->at(j));
					tempADC = (Double_t)avADCDia1->GetBinContent(col->at(j)+1,row->at(j)+1)+adc->at(j);
					avADCDia1->SetBinContent(col->at(j)+1,row->at(j)+1,tempADC);
					adcDia1->Fill(adc->at(j));
				}
				else if(plane->at(j) == 5){
					hitPlane5->Fill(col->at(j),row->at(j));
					tempADC = (Double_t)avADCDia2->GetBinContent(col->at(j)+1,row->at(j)+1)+adc->at(j);
					avADCDia2->SetBinContent(col->at(j)+1,row->at(j)+1,tempADC);
					adcDia2->Fill(adc->at(j));
				}
				else if(plane->at(j) == 6){
					hitPlane6->Fill(col->at(j),row->at(j));
					tempADC = (Double_t)avADCSi->GetBinContent(col->at(j)+1,row->at(j)+1)+adc->at(j);
					avADCSi->SetBinContent(col->at(j)+1,row->at(j)+1,tempADC);
					adcSi->Fill(adc->at(j));
				}
				else{
					cout << "ERROR. PLANE DOES NOT EXIST"<<endl;
					break;
				}
			}
		}

		// Get the number of cluster for each DUT roc
		for(Int_t i = 0; i < 7; i++){
			numCLROC[i] = (Int_t)clust_per_plane->at(i);
		}
		// For Diamond 1
		PhHistogramExtraction(numCLROC[4], ph_ROC4_1cl, clust_ROC4_X, clust_ROC4_Y, phDia1, avPHDia1, hitMapPHDia1);
		// For Diamond 2
		PhHistogramExtraction(numCLROC[5], ph_ROC5_1cl, clust_ROC5_X, clust_ROC5_Y, phDia2, avPHDia2, hitMapPHDia2);
		// For Silicon
		PhHistogramExtraction(numCLROC[6], ph_ROC6_1cl, clust_ROC6_X, clust_ROC6_Y, phSi, avPHSi, hitMapPHSi);
		// For Plane 0
		PhHistogramExtraction(numCLROC[0], ph_ROC0_1cl, clust_ROC0_X, clust_ROC0_Y, phPlane0, avPHPlane0, hitMapPHPlane0);
		// For Plane 1
		PhHistogramExtraction(numCLROC[1], ph_ROC1_1cl, clust_ROC1_X, clust_ROC1_Y, phPlane1, avPHPlane1, hitMapPHPlane1);
		// For Plane 2
		PhHistogramExtraction(numCLROC[2], ph_ROC2_1cl, clust_ROC2_X, clust_ROC2_Y, phPlane2, avPHPlane2, hitMapPHPlane2);
		// For Plane 3
		PhHistogramExtraction(numCLROC[3], ph_ROC3_1cl, clust_ROC3_X, clust_ROC3_Y, phPlane3, avPHPlane3, hitMapPHPlane3);
	}

	// average for ADC DUT
	FindAverageHistogramDUT(avADCDia1, hitPlane4, avADCDia2, hitPlane5, avADCSi, hitPlane6, 52, 80);

	// average for PH DUT
	FindAverageHistogramDUT(avPHDia1,hitMapPHDia1,avPHDia2,hitMapPHDia2,avPHSi,hitMapPHSi,divXR1,divYR1);

	//average for PH Telescope Planes
	FindAverageHistogramTPlanes(avPHPlane0,hitMapPHPlane0,avPHPlane1,hitMapPHPlane1,avPHPlane2,hitMapPHPlane2,avPHPlane3,hitMapPHPlane3,divXR1,divYR1);

	// Plot histograms
	PlotHistogram2D(hitPlane0,"Hit map Plane 0",outputPath,"HitMapPlane0.png");
	PlotHistogram2D(hitPlane1,"Hit map Plane 1",outputPath,"HitMapPlane1.png");
	PlotHistogram2D(hitPlane2,"Hit map Plane 2",outputPath,"HitMapPlane2.png");
	PlotHistogram2D(hitPlane3,"Hit map Plane 3",outputPath,"HitMapPlane3.png");
	PlotHistogram2D(hitPlane4,"Hit map Plane 4",outputPath,"HitMapPlane4.png");
	PlotHistogram2D(hitPlane5,"Hit map Plane 5",outputPath,"HitMapPlane5.png");
	PlotHistogram2D(hitPlane6,"Hit map Plane 6",outputPath,"HitMapPlane6.png");

	PlotHistogram2D(avADCDia1,"ADC Average Diamond 1",outputPath,"AvADCD1.png");
	PlotHistogram2D(avADCDia2,"ADC Average Diamond 2",outputPath,"AvADCD2.png");
	PlotHistogram2D(avADCSi,"ADC Average Silicon",outputPath,"AvADCSi.png");

	PlotHistogram1D(adcDia1, kTRUE, kFALSE, "ADC Diamond 1", outputPath, "ADCD1.png");
	PlotHistogram1D(adcDia2, kTRUE, kFALSE, "ADC Diamond 2", outputPath, "ADCD2.png");
	PlotHistogram1D(adcSi, kTRUE, kFALSE, "ADC Silicon", outputPath, "ADCSi.png");

	PlotHistogram2D(avPHDia1,"PH Average Diamond 1",outputPath,"AvPHD1.png");
	PlotHistogram2D(avPHDia2,"PH Average Diamond 2",outputPath,"AvPHD2.png");
	PlotHistogram2D(avPHSi,"PH Average Silicon",outputPath,"AvPHSi.png");

	PlotHistogram2D(avPHPlane0,"PH Average Plane 0",outputPath,"AvPHTPlane0.png");
	PlotHistogram2D(avPHPlane1,"PH Average Plane 1",outputPath,"AvPHTPlane1.png");
	PlotHistogram2D(avPHPlane2,"PH Average Plane 2",outputPath,"AvPHTPlane2.png");
	PlotHistogram2D(avPHPlane3,"PH Average Plane 3",outputPath,"AvPHTPlane3.png");

	PlotHistogram2D(hitMapPHDia1,"Hit Map Diamond 1 for PH",outputPath,"HitMapPHDia1.png");
	PlotHistogram2D(hitMapPHDia2,"Hit Map Diamond 2 for PH",outputPath,"HitMapPHDia2.png");
	PlotHistogram2D(hitMapPHSi,"Hit Map Silicon for PH",outputPath,"HitMapPHSi.png");

	PlotHistogram2D(hitMapPHPlane0,"Hit Map Plane 0 for PH",outputPath,"HitMapPHTPlane0.png");
	PlotHistogram2D(hitMapPHPlane1,"Hit Map Plane 1 for PH",outputPath,"HitMapPHTPlane1.png");
	PlotHistogram2D(hitMapPHPlane2,"Hit Map Plane 2 for PH",outputPath,"HitMapPHTPlane2.png");
	PlotHistogram2D(hitMapPHPlane3,"Hit Map Plane 3 for PH",outputPath,"HitMapPHTPlane3.png");

	PlotHistogram1D(phDia1, kTRUE, kFALSE, "PH Diamond 1", outputPath, "PHD1.png");
	PlotHistogram1D(phDia2, kTRUE, kFALSE, "PH Diamond 2", outputPath, "PHD2.png");
	PlotHistogram1D(phSi, kTRUE, kFALSE, "PH Silicon", outputPath, "PHSi.png");

	PlotHistogram1D(phPlane0, kTRUE, kFALSE, "PH Plane 0", outputPath, "PHTPlane0.png");
	PlotHistogram1D(phPlane1, kTRUE, kFALSE, "PH Plane 1", outputPath, "PHTPlane1.png");
	PlotHistogram1D(phPlane2, kTRUE, kFALSE, "PH Plane 2", outputPath, "PHTPlane2.png");
	PlotHistogram1D(phPlane3, kTRUE, kFALSE, "PH Plane 3", outputPath, "PHTPlane3.png");

	// Save histograms
	TFile f(Form("%s%s",outputPath,histogramOutputFile),"RECREATE");
	adcDia1->Write();
	adcDia2->Write();
	adcSi->Write();
	phDia1->Write();
	phDia2->Write();
	phSi->Write();
	hitPlane0->Write();
	hitPlane1->Write();
	hitPlane2->Write();
	hitPlane3->Write();
	hitPlane4->Write();
	hitPlane5->Write();
	hitPlane6->Write();
	avADCDia1->Write();
	avADCDia2->Write();
	avADCSi->Write();

	avPHDia1->Write();
	avPHDia2->Write();
	avPHSi->Write();
	avPHPlane0->Write();
	avPHPlane1->Write();
	avPHPlane2->Write();
	avPHPlane3->Write();

	phDia1->Write();
	phDia2->Write();
	phSi->Write();
	phPlane0->Write();
	phPlane1->Write();
	phPlane2->Write();
	phPlane3->Write();

	hitMapPHDia1->Write();
	hitMapPHDia2->Write();
	hitMapPHSi->Write();
	hitMapPHPlane0->Write();
	hitMapPHPlane1->Write();
	hitMapPHPlane2->Write();
	hitMapPHPlane3->Write();

	f.Close();


}

void PhHistogramExtraction(int numCLROC, vector<float> *phROC, vector<float> *clust_ROC_X, vector<float> *clust_ROC_Y, TH1F *phDUT, TH2F *avPHDUT, TH2F *hitMapPHDUT){
	if(numCLROC != 0){
		int numPhROC = phROC->size();
		Int_t tempBinX = 0, tempBinY = 0;
		Double_t tempPH = 0;
		for(Int_t i = 0; i < numPhROC; i++){
			if(clust_ROC_X->size() >= numPhROC && clust_ROC_Y->size() >= numPhROC){
				tempBinX = (Int_t)(TMath::Ceil((Float_t)((clust_ROC_X->at(i)*10-xminR1)/deltaXR1)));
				tempBinY = (Int_t)(TMath::Ceil((Float_t)((clust_ROC_Y->at(i)*10-yminR1)/deltaYR1)));
				phDUT->Fill(phROC->at(i));
				if( tempBinX >= 1 && tempBinY >= 1 && tempBinX <= divXR1 && tempBinY <= divYR1){
					tempPH = (Double_t)(avPHDUT->GetBinContent(tempBinX,tempBinY)+phROC->at(i));
					avPHDUT->SetBinContent(tempBinX,tempBinY,tempPH);
					hitMapPHDUT->Fill(tempBinX,tempBinY);
				}
			}
		}
	}
}

void PlotHistogram1D(TH1F *histogram, Bool_t logY, Bool_t logX, const char *title, const char *outputPath, const char *suffix){
	TCanvas *c1 = new TCanvas("c1",title,1);
	c1->cd();
	if(logX)
		c1->SetLogx();
	if(logY)
		c1->SetLogy();
	histogram->Draw("E1 SAME HIST");
	c1->SaveAs(Form("%s%s",outputPath,suffix));
	delete c1;
}

void PlotHistogram2D(TH2F *histogram, const char *title, const char *outputPath, const char *suffix){
	TCanvas *c1 = new TCanvas("c1",title,1);
	c1->cd();
	histogram->Draw("COLZ");
	c1->SaveAs(Form("%s%s",outputPath,suffix));
	delete c1;
}

void FindAverageHistogramDUT(TH2F *histToAvDUT1, TH2F *histHitsDUT1, TH2F *histToAvDUT2, TH2F *histHitsDUT2, TH2F *histToAvDUT3, TH2F *histHitsDUT3, Int_t xbins, Int_t ybins){
	for (Int_t i = 1; i <= xbins; i++){
		for (Int_t j = 1; j <= ybins; j++){
			if(histHitsDUT1->GetBinContent(i,j) >= 1){
				AverageBinHistogram(histToAvDUT1,histHitsDUT1,i,j);
			}
			if(histHitsDUT2->GetBinContent(i,j) >= 1){
				AverageBinHistogram(histToAvDUT2,histHitsDUT2,i,j);
			}
			if(histHitsDUT3->GetBinContent(i,j) >= 1){
				AverageBinHistogram(histToAvDUT3,histHitsDUT3,i,j);
			}
		}
	}
}

void FindAverageHistogramTPlanes(TH2F *histToAvTPlane0, TH2F *histHitsTPlane0, TH2F *histToAvTPlane1, TH2F *histHitsTPlane1, TH2F *histToAvTPlane2, TH2F *histHitsTPlane2, TH2F *histToAvTPlane3, TH2F *histHitsTPlane3, Int_t xbins, Int_t ybins){
	for (Int_t i = 1; i <= xbins; i++){
		for (Int_t j = 1; j <= ybins; j++){
			if(histHitsTPlane0->GetBinContent(i,j) >= 1){
				AverageBinHistogram(histToAvTPlane0,histHitsTPlane0,i,j);
			}
			if(histHitsTPlane1->GetBinContent(i,j) >= 1){
				AverageBinHistogram(histToAvTPlane1,histHitsTPlane1,i,j);
			}
			if(histHitsTPlane2->GetBinContent(i,j) >= 1){
				AverageBinHistogram(histToAvTPlane2,histHitsTPlane2,i,j);
			}
			if(histHitsTPlane3->GetBinContent(i,j) >= 1){
				AverageBinHistogram(histToAvTPlane3,histHitsTPlane3,i,j);
			}
		}
	}
}

void AverageBinHistogram(TH2F *histToAv, TH2F *histHits, Int_t binx, Int_t biny){
	Double_t temp = (Double_t)(histToAv->GetBinContent(binx,biny)/(Double_t)histHits->GetBinContent(binx,biny));
	histToAv->SetBinContent(binx,biny,temp);
}