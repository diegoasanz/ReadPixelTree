//============================================================================
// Name        : ReadPixelTree.cpp
// Author      : Diego Alejandro
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <stdlib.h>
#include <stdint.h>

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
#include "TGraphErrors.h"
#include "TStyle.h"

void Prueba(char *rootFilePath, char *outputPath);

TChain * ReadPixRootFile(const char *fileName);

void Ph1DHistogramsExtraction(int numCLROC, vector<float> *ph1, vector<float> *ph2, vector<float> *ph3, vector<float> *phM4, TH1F *phHall, TH1F *phH1, TH1F *phH2, TH1F *phH3, TH1F *phHM4);

void Ph2DHistogramExtraction(int numCLROC, vector<float> *phROC, vector<float> *clust_ROC_X, vector<float> *clust_ROC_Y, TH2F *avPHDUT, TH2F *hitMapPHDUT, float deltaX, float deltaY, Bool_t isRowCol);

void PlotPulseHeightsOverlay(TH1F *phall, TH1F *ph1, TH1F *ph2, TH1F *ph3, TH1F *phM4, Bool_t logY, Bool_t logX, const char *title, const char *outputPath, const char *suffix);

void PlotHistogram2D(TH2F *histogram, const char *title, const char *outputPath, const char *suffix);

void FindAverageHistogramDUT(TH2F *histToAvDUT1, TH2F *histHitsDUT1, TH2F *histToAvDUT2, TH2F *histHitsDUT2, TH2F *histToAvDUT3, TH2F *histHitsDUT3, Int_t xbins, Int_t ybins);

void FindAverageHistogramTPlanes(TH2F *histToAvTPlane0, TH2F *histHitsTPlane0, TH2F *histToAvTPlane1, TH2F *histHitsTPlane1, TH2F *histToAvTPlane2, TH2F *histHitsTPlane2, TH2F *histToAvTPlane3, TH2F *histHitsTPlane3, Int_t xbins, Int_t ybins);

void AverageBinHistogram(TH2F *histToAv, TH2F *histHits, Int_t binx, Int_t biny);

void MaskPixelsFromFile(const string maskFile);

Bool_t IsMasked(Int_t roc, Int_t col, Int_t row);

const char *rootFilePathComplete = "~/Data/psi_2015_10/root/pixel/TrackedRun312.root";
const char *rootFilePathSmall = "~/DataSmall/psi_2015_10/root/pixel/TrackedRun312.root";

const char *outputPathComplete = "~/Dropbox/201601/DiamondPixels/ReadPixelTree/";
const char *outputPathSmall = "~/Dropbox/201601/DiamondPixels/ReadPixelTree/Small/";

const char *histogramOutputFile = "Histograms.root";

const string maskFile = "~/Dropbox/201601/DiamondPixels/ReadPixelTree/maskFile.msk";

vector<int> *mskdroc = 0, *mskdcol = 0, *mskdrow = 0;

Int_t contour2D = 1024, mean1Dbins = 100, numBinCol = 52, numBinRow = 80, minCol = 0, maxCol = 51, minRow = 0, maxRow = 79, binWidthTGraph=2000;
Float_t xminR1 = -6.075, xmaxR1 = 6.075, yminR1 = -6.05, ymaxR1 = 6.05, ph1Dmin=0, ph1Dmax=50000;
Int_t divXR1 = 81, divYR1 = 121, ph1Dbins=100;
Float_t deltaXR1 = (Float_t)((xmaxR1-xminR1)/divXR1), deltaYR1 = (Float_t)((ymaxR1-yminR1)/divYR1);

int main() {
	gROOT->ProcessLine("#include <vector>");
	gSystem->Load("libMathCore");
	gSystem->Load("libPhysics");
	const char *rootFilePath = rootFilePathSmall;
	const char *outputPath = outputPathSmall;
	cout << "Starting Analysis" << endl;
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
	cout << "Reading Root file: " << rootFilePath << endl;
	TChain *chain1 = ReadPixRootFile(rootFilePath);
	//Make pointers to branches
	TBranch *braPlane = chain1->GetBranch("plane");
	//Create variables to read the tree
	vector<int> *plane=0, *col=0, *row=0, *adc=0;
	vector<vector<float> *> clust_row, clust_col;
	vector<unsigned char> *clust_per_plane=0;
	vector<vector<float> *> charge_all, clust_Telescope_X, clust_Telescope_Y, clust_Local_X, clust_Local_Y, ph_1cl, ph_2cl, ph_3cl, ph_M4cl;
	Float_t time=0, dia1X=0, dia1Y=0, dia2X=0, dia2Y=0, chi2_tracks=0, chi2_x=0, chi2_y=0, slope_x=0, slope_y=0, coincidence_map=0;
	Int_t event_number=0;
	UChar_t hit_plane_bits=0, n_tracks=0, n_clusters=0;
	//Int_t n_tracks = 0, n_cluster = 0;
	//dimension variables
	vector<int> numCLROC;
	numCLROC.resize(7);
	Int_t numPlane = 0, numCol = 0, numRow = 0, numADC = 0;
	Long64_t numEntries = braPlane->GetEntries();
	Long64_t maxphplots=(Long64_t)8*numEntries/100;

	//MaskPixelsFromFile(maskFile);

	clust_row.resize(7);
	clust_col.resize(7);
	charge_all.resize(7);
	clust_Telescope_X.resize(7);
	clust_Telescope_Y.resize(7);
	clust_Local_X.resize(7);
	clust_Local_Y.resize(7);
	ph_1cl.resize(7);
	ph_2cl.resize(7);
	ph_3cl.resize(7);
	ph_M4cl.resize(7);

	chain1->SetBranchAddress("event_number",&event_number);
	chain1->SetBranchAddress("time",&time);
	chain1->SetBranchAddress("plane",&plane);
	chain1->SetBranchAddress("col",&col);
	chain1->SetBranchAddress("row",&row);
	chain1->SetBranchAddress("adc",&adc);
	chain1->SetBranchAddress("hit_plane_bits",&hit_plane_bits);
	chain1->SetBranchAddress("diam1_track_x",&dia1X);
	chain1->SetBranchAddress("diam1_track_y",&dia1Y);
	chain1->SetBranchAddress("diam2_track_x",&dia2X);
	chain1->SetBranchAddress("diam2_track_y",&dia2Y);
	chain1->SetBranchAddress("chi2_tracks",&chi2_tracks);
	chain1->SetBranchAddress("chi2_x",&chi2_x);
	chain1->SetBranchAddress("chi2_y",&chi2_y);
	chain1->SetBranchAddress("slope_x",&slope_x);
	chain1->SetBranchAddress("slope_y",&slope_y);
	chain1->SetBranchAddress("n_tracks",&n_tracks);
	chain1->SetBranchAddress("n_clusters",&n_clusters);
	chain1->SetBranchAddress("clusters_per_plane",&clust_per_plane);
	chain1->SetBranchAddress("coincidence_map",&coincidence_map);
	// ROCs
	for(Int_t iROC = 0; iROC<7; iROC++){
		TString branch_name_charge = TString::Format("charge_all_ROC%d",iROC);
		chain1->SetBranchAddress(branch_name_charge,&(charge_all[iROC]));
		TString branch_name_cluster_Telescope_X = TString::Format("cluster_pos_ROC%d_Telescope_X",iROC);
		chain1->SetBranchAddress(branch_name_cluster_Telescope_X,&(clust_Telescope_X[iROC]));
		TString branch_name_cluster_Telescope_Y = TString::Format("cluster_pos_ROC%d_Telescope_Y",iROC);
		chain1->SetBranchAddress(branch_name_cluster_Telescope_Y,&(clust_Telescope_Y[iROC]));
		TString branch_name_cluster_Local_X = TString::Format("cluster_pos_ROC%d_Local_X",iROC);
		chain1->SetBranchAddress(branch_name_cluster_Local_X,&(clust_Local_X[iROC]));
		TString branch_name_cluster_Local_Y = TString::Format("cluster_pos_ROC%d_Local_Y",iROC);
		chain1->SetBranchAddress(branch_name_cluster_Local_Y,&(clust_Local_Y[iROC]));
		TString branch_name_cluster_row = TString::Format("cluster_row_ROC%d",iROC);
		chain1->SetBranchAddress(branch_name_cluster_row,&(clust_row[iROC]));
		TString branch_name_cluster_col = TString::Format("cluster_col_ROC%d",iROC);
		chain1->SetBranchAddress(branch_name_cluster_col,&(clust_col[iROC]));
		TString branch_name_pulse_height_1pix_cluster = TString::Format("pulse_height_ROC%d_1_cluster",iROC);
		chain1->SetBranchAddress(branch_name_pulse_height_1pix_cluster,&(ph_1cl[iROC]));
		TString branch_name_pulse_height_2pix_cluster = TString::Format("pulse_height_ROC%d_2_cluster",iROC);
		chain1->SetBranchAddress(branch_name_pulse_height_2pix_cluster,&(ph_2cl[iROC]));
		TString branch_name_pulse_height_3pix_cluster = TString::Format("pulse_height_ROC%d_3_cluster",iROC);
		chain1->SetBranchAddress(branch_name_pulse_height_3pix_cluster,&(ph_3cl[iROC]));
		TString branch_name_pulse_height_M4pix_cluster = TString::Format("pulse_height_ROC%d_More4_cluster",iROC);
		chain1->SetBranchAddress(branch_name_pulse_height_M4pix_cluster,&(ph_M4cl[iROC]));
	}


	// 1D Histograms
	TH1F *coincidenceMap = new TH1F("coincidenceMap","Coincidence Map",141,-0.5,140.5);

	// ROCs
	TH1F *phROC_all[7], *phROC_1cl[7], *phROC_2cl[7], *phROC_3cl[7], *phROC_M4cl[7];
	TGraphErrors *meanPhROC_all[7], *meanPhROC_1cl[7], *meanPhROC_2cl[7], *meanPhROC_3cl[7], *meanPhROC_M4cl[7];
	TH2F *hitMap[7], *avPhROC_local_1cl[7], *phROC_hitMap_local_1cl[7], *avPhROC_pix_1cl[7], *phROC_hitMap_pix_1cl[7], *avPhROC_telescope_1cl[7], *phROC_hitMap_telescope_1cl[7];
	for(Int_t iROC = 0; iROC<7; iROC++){
		// 1D Histograms Draw with ape and same
		TString hist_name_pulse_height_all = TString::Format("phROC%d_all",iROC);
		TString hist_title_pulse_height_all = TString::Format("Pulse Height ROC%d all cluster",iROC);
		phROC_all[iROC] = new TH1F(hist_name_pulse_height_all,hist_title_pulse_height_all,ph1Dbins+1,ph1Dmin-(ph1Dmax-ph1Dmin)/(2*ph1Dbins),ph1Dmax+(ph1Dmax-ph1Dmin)/(2*ph1Dbins));
		phROC_all[iROC]->GetXaxis()->SetTitle("Charge (e)");
		phROC_all[iROC]->GetYaxis()->SetTitle("Num Clusters");
		phROC_all[iROC]->SetMaximum(maxphplots);
		phROC_all[iROC]->SetMinimum(0);
		phROC_all[iROC]->SetLineColor(1); // black
		TString hist_name_pulse_height_1cl = TString::Format("phROC%d_1cl",iROC);
		TString hist_title_pulse_height_1cl = TString::Format("Pulse Height ROC%d 1pix cluster",iROC);
		phROC_1cl[iROC] = new TH1F(hist_name_pulse_height_1cl,hist_title_pulse_height_1cl,ph1Dbins+1,ph1Dmin-(ph1Dmax-ph1Dmin)/(2*ph1Dbins),ph1Dmax+(ph1Dmax-ph1Dmin)/(2*ph1Dbins));
		phROC_1cl[iROC]->GetXaxis()->SetTitle("Charge (e)");
		phROC_1cl[iROC]->GetYaxis()->SetTitle("Num Clusters");
		phROC_1cl[iROC]->SetMaximum(maxphplots);
		phROC_1cl[iROC]->SetMinimum(0);
		phROC_1cl[iROC]->SetLineColor(4); // blue
		TString hist_name_pulse_height_2cl = TString::Format("phROC%d_2cl",iROC);
		TString hist_title_pulse_height_2cl = TString::Format("Pulse Height ROC%d 2pix cluster",iROC);
		phROC_2cl[iROC] = new TH1F(hist_name_pulse_height_2cl,hist_title_pulse_height_2cl,ph1Dbins+1,ph1Dmin-(ph1Dmax-ph1Dmin)/(2*ph1Dbins),ph1Dmax+(ph1Dmax-ph1Dmin)/(2*ph1Dbins));
		phROC_2cl[iROC]->GetXaxis()->SetTitle("Charge (e)");
		phROC_2cl[iROC]->GetYaxis()->SetTitle("Num Clusters");
		phROC_2cl[iROC]->SetMaximum(maxphplots);
		phROC_2cl[iROC]->SetMinimum(0);
		phROC_2cl[iROC]->SetLineColor(3); // green
		TString hist_name_pulse_height_3cl = TString::Format("phROC%d_3cl",iROC);
		TString hist_title_pulse_height_3cl = TString::Format("Pulse Height ROC%d 3pix cluster",iROC);
		phROC_3cl[iROC] = new TH1F(hist_name_pulse_height_3cl,hist_title_pulse_height_3cl,ph1Dbins+1,ph1Dmin-(ph1Dmax-ph1Dmin)/(2*ph1Dbins),ph1Dmax+(ph1Dmax-ph1Dmin)/(2*ph1Dbins));
		phROC_3cl[iROC]->GetXaxis()->SetTitle("Charge (e)");
		phROC_3cl[iROC]->GetYaxis()->SetTitle("Num Clusters");
		phROC_3cl[iROC]->SetMaximum(maxphplots);
		phROC_3cl[iROC]->SetMinimum(0);
		phROC_3cl[iROC]->SetLineColor(2); // red
		TString hist_name_pulse_height_M4cl = TString::Format("phROC%d_M4cl",iROC);
		TString hist_title_pulse_height_M4cl = TString::Format("Pulse Height ROC%d M4pix cluster",iROC);
		phROC_M4cl[iROC] = new TH1F(hist_name_pulse_height_M4cl,hist_title_pulse_height_M4cl,ph1Dbins+1,ph1Dmin-(ph1Dmax-ph1Dmin)/(2*ph1Dbins),ph1Dmax+(ph1Dmax-ph1Dmin)/(2*ph1Dbins));
		phROC_M4cl[iROC]->GetXaxis()->SetTitle("Charge (e)");
		phROC_M4cl[iROC]->GetYaxis()->SetTitle("Num Clusters");
		phROC_M4cl[iROC]->SetMaximum(maxphplots);
		phROC_M4cl[iROC]->SetMinimum(0);
		phROC_M4cl[iROC]->SetLineColor(6); // magenta
		// 1D TGraphErrors Draw with hist and samehist
		meanPhROC_all[iROC] = new TGraphErrors(mean1Dbins);
		TString tgraph_name_mean_pulse_height_all = TString::Format("meanPHROC%d_all",iROC);
		TString tgraph_title_mean_pulse_height_all = TString::Format("Mean Pulse Height ROC%d all cluster",iROC);
		meanPhROC_all[iROC]->SetName(tgraph_name_mean_pulse_height_all);
		meanPhROC_all[iROC]->SetTitle(tgraph_title_mean_pulse_height_all);
		meanPhROC_all[iROC]->GetXaxis()->SetTitle("Event Number");
		meanPhROC_all[iROC]->GetYaxis()->SetTitle("Average Pulse Height");
		meanPhROC_all[iROC]->SetMinimum(0);
		meanPhROC_all[iROC]->SetMaximum(numEntries);
		meanPhROC_all[iROC]->SetLineColor(1); // black
		meanPhROC_all[iROC]->SetMarkerColor(1); // black

		meanPhROC_1cl[iROC] = new TGraphErrors(mean1Dbins);
		TString tgraph_name_mean_pulse_height_1cl = TString::Format("meanPHROC%d_1cl",iROC);
		TString tgraph_title_mean_pulse_height_1cl = TString::Format("Mean Pulse Height ROC%d 1pix cluster",iROC);
		meanPhROC_1cl[iROC]->SetName(tgraph_name_mean_pulse_height_1cl);
		meanPhROC_1cl[iROC]->SetTitle(tgraph_title_mean_pulse_height_1cl);
		meanPhROC_1cl[iROC]->GetXaxis()->SetTitle("Event Number");
		meanPhROC_1cl[iROC]->GetYaxis()->SetTitle("Average Pulse Height");
		meanPhROC_1cl[iROC]->SetMinimum(0);
		meanPhROC_1cl[iROC]->SetMaximum(numEntries);
		meanPhROC_1cl[iROC]->SetLineColor(4); // blue
		meanPhROC_1cl[iROC]->SetMarkerColor(4); // blue

		meanPhROC_2cl[iROC] = new TGraphErrors(mean1Dbins);
		TString tgraph_name_mean_pulse_height_2cl = TString::Format("meanPHROC%d_2cl",iROC);
		TString tgraph_title_mean_pulse_height_2cl = TString::Format("Mean Pulse Height ROC%d 2pix cluster",iROC);
		meanPhROC_2cl[iROC]->SetName(tgraph_name_mean_pulse_height_2cl);
		meanPhROC_2cl[iROC]->SetTitle(tgraph_title_mean_pulse_height_2cl);
		meanPhROC_2cl[iROC]->GetXaxis()->SetTitle("Event Number");
		meanPhROC_2cl[iROC]->GetYaxis()->SetTitle("Average Pulse Height");
		meanPhROC_2cl[iROC]->SetMinimum(0);
		meanPhROC_2cl[iROC]->SetMaximum(numEntries);
		meanPhROC_2cl[iROC]->SetLineColor(3); // green
		meanPhROC_2cl[iROC]->SetMarkerColor(3); // green

		meanPhROC_3cl[iROC] = new TGraphErrors(mean1Dbins);
		TString tgraph_name_mean_pulse_height_3cl = TString::Format("meanPHROC%d_3cl",iROC);
		TString tgraph_title_mean_pulse_height_3cl = TString::Format("Mean Pulse Height ROC%d 3pix cluster",iROC);
		meanPhROC_3cl[iROC]->SetName(tgraph_name_mean_pulse_height_3cl);
		meanPhROC_3cl[iROC]->SetTitle(tgraph_title_mean_pulse_height_3cl);
		meanPhROC_3cl[iROC]->GetXaxis()->SetTitle("Event Number");
		meanPhROC_3cl[iROC]->GetYaxis()->SetTitle("Average Pulse Height");
		meanPhROC_3cl[iROC]->SetMinimum(0);
		meanPhROC_3cl[iROC]->SetMaximum(numEntries);
		meanPhROC_3cl[iROC]->SetLineColor(2); // red
		meanPhROC_3cl[iROC]->SetMarkerColor(2); // red

		meanPhROC_M4cl[iROC] = new TGraphErrors(mean1Dbins);
		TString tgraph_name_mean_pulse_height_M4cl = TString::Format("meanPHROC%d_M4cl",iROC);
		TString tgraph_title_mean_pulse_height_M4cl = TString::Format("Mean Pulse Height ROC%d M4pix cluster",iROC);
		meanPhROC_M4cl[iROC]->SetName(tgraph_name_mean_pulse_height_M4cl);
		meanPhROC_M4cl[iROC]->SetTitle(tgraph_title_mean_pulse_height_M4cl);
		meanPhROC_M4cl[iROC]->GetXaxis()->SetTitle("Event Number");
		meanPhROC_M4cl[iROC]->GetYaxis()->SetTitle("Average Pulse Height");
		meanPhROC_M4cl[iROC]->SetMinimum(0);
		meanPhROC_M4cl[iROC]->SetMaximum(numEntries);
		meanPhROC_M4cl[iROC]->SetLineColor(6); // magenta
		meanPhROC_M4cl[iROC]->SetMarkerColor(6); // magenta

		// 2D Histograms
		// ROC
		TString hist_name_hitMap = TString::Format("hitMapROC%d",iROC);
		TString hist_title_histMap = TString::Format("Hit Map ROC%d",iROC);
		hitMap[iROC] = new TH2F(hist_name_hitMap,hist_title_histMap,numBinCol,minCol-0.5,maxCol+0.5,numBinRow,minRow-0.5,maxRow+0.5);
		hitMap[iROC]->GetXaxis()->SetTitle("Column");
		hitMap[iROC]->GetYaxis()->SetTitle("Row");
		hitMap[iROC]->SetContour(contour2D);
		hitMap[iROC]->SetMinimum(0);
		// LOCAL
		TString hist_name_pulse_height_local_1cl = TString::Format("avPh_ROC%d_local_1cl",iROC);
		TString hist_title_pulse_height_local_1cl = TString::Format("Average Pulse Height ROC%d Local Coord. 1pix cluster",iROC);
		avPhROC_local_1cl[iROC] = new TH2F(hist_name_pulse_height_local_1cl,hist_title_pulse_height_local_1cl,divXR1,xminR1,xmaxR1,divYR1,yminR1,ymaxR1);
		avPhROC_local_1cl[iROC]->GetXaxis()->SetTitle("x (mm)");
		avPhROC_local_1cl[iROC]->GetYaxis()->SetTitle("y (mm)");
		avPhROC_local_1cl[iROC]->SetContour(contour2D);
		avPhROC_local_1cl[iROC]->SetMinimum(0);

		TString hist_name_pulse_height_hitMap_local_1cl = TString::Format("ph_ROC%d_hitMap_local_1cl",iROC);
		TString hist_title_pulse_height_hitMap_local_1cl = TString::Format("Pulse Height ROC%d Hit Map Local Coord. 1pix cluster",iROC);
		phROC_hitMap_local_1cl[iROC] = new TH2F(hist_name_pulse_height_hitMap_local_1cl,hist_title_pulse_height_hitMap_local_1cl,divXR1,xminR1,xmaxR1,divYR1,yminR1,ymaxR1);
		phROC_hitMap_local_1cl[iROC]->GetXaxis()->SetTitle("x (mm)");
		phROC_hitMap_local_1cl[iROC]->GetYaxis()->SetTitle("y (mm)");
		phROC_hitMap_local_1cl[iROC]->SetContour(contour2D);
		phROC_hitMap_local_1cl[iROC]->SetMinimum(0);
		// TELESCOPE
		TString hist_name_pulse_height_telescope_1cl = TString::Format("avPh_ROC%d_telescope_1cl",iROC);
		TString hist_title_pulse_height_telescope_1cl = TString::Format("Average Pulse Height ROC%d telescope Coord. 1pix cluster",iROC);
		avPhROC_telescope_1cl[iROC] = new TH2F(hist_name_pulse_height_telescope_1cl,hist_title_pulse_height_telescope_1cl,divXR1,xminR1,xmaxR1,divYR1,yminR1,ymaxR1);
		avPhROC_telescope_1cl[iROC]->GetXaxis()->SetTitle("x (mm)");
		avPhROC_telescope_1cl[iROC]->GetYaxis()->SetTitle("y (mm)");
		avPhROC_telescope_1cl[iROC]->SetContour(contour2D);
		avPhROC_telescope_1cl[iROC]->SetMinimum(0);

		TString hist_name_pulse_height_hitMap_telescope_1cl = TString::Format("ph_ROC%d_hitMap_telescope_1cl",iROC);
		TString hist_title_pulse_height_hitMap_telescope_1cl = TString::Format("Pulse Height ROC%d Hit Map telescope Coord. 1pix cluster",iROC);
		phROC_hitMap_telescope_1cl[iROC] = new TH2F(hist_name_pulse_height_hitMap_telescope_1cl,hist_title_pulse_height_hitMap_telescope_1cl,divXR1,xminR1,xmaxR1,divYR1,yminR1,ymaxR1);
		phROC_hitMap_telescope_1cl[iROC]->GetXaxis()->SetTitle("x (mm)");
		phROC_hitMap_telescope_1cl[iROC]->GetYaxis()->SetTitle("y (mm)");
		phROC_hitMap_telescope_1cl[iROC]->SetContour(contour2D);
		phROC_hitMap_telescope_1cl[iROC]->SetMinimum(0);
		// ROC
		TString hist_name_pulse_height_pix_1cl = TString::Format("avPh_ROC%d_pix_1cl",iROC);
		TString hist_title_pulse_height_pix_1cl = TString::Format("Average Pulse Height ROC%d pixelated Coord. 1pix cluster",iROC);
		avPhROC_pix_1cl[iROC] = new TH2F(hist_name_pulse_height_pix_1cl,hist_title_pulse_height_pix_1cl,numBinCol,minCol-0.5,maxCol+0.5,numBinRow,minRow-0.5,maxRow+0.5);
		avPhROC_pix_1cl[iROC]->GetXaxis()->SetTitle("Column");
		avPhROC_pix_1cl[iROC]->GetYaxis()->SetTitle("Row");
		avPhROC_pix_1cl[iROC]->SetContour(contour2D);
		avPhROC_pix_1cl[iROC]->SetMinimum(0);

		TString hist_name_pulse_height_hitMap_pix_1cl = TString::Format("ph_ROC%d_hitMap_pix_1cl",iROC);
		TString hist_title_pulse_height_hitMap_pix_1cl = TString::Format("Pulse Height ROC%d Hit Map pixelated Coord. 1pix cluster",iROC);
		phROC_hitMap_pix_1cl[iROC] = new TH2F(hist_name_pulse_height_hitMap_pix_1cl,hist_title_pulse_height_hitMap_pix_1cl,numBinCol,minCol-0.5,maxCol+0.5,numBinRow,minRow-0.5,maxRow+0.5);
		phROC_hitMap_pix_1cl[iROC]->GetXaxis()->SetTitle("Column");
		phROC_hitMap_pix_1cl[iROC]->GetYaxis()->SetTitle("Row");
		phROC_hitMap_pix_1cl[iROC]->SetContour(contour2D);
		phROC_hitMap_pix_1cl[iROC]->SetMinimum(0);
	}

	cout << "Reading first entry." << endl;

	for(Long64_t i = 0; i < numEntries; i++){
		chain1->GetEntry(i);
		numPlane = plane->size();
		numCol = col->size();
		numRow = row->size();
		numADC = adc->size();
		coincidenceMap->Fill(coincidence_map);
		if(numPlane & numCol & numRow == numADC){

			for(Int_t j = 0; j < numPlane; j++){
				if(plane->at(j)==0){
					hitMap[0]->Fill(col->at(j),row->at(j));
				}
				else if(plane->at(j) == 1){
					hitMap[1]->Fill(col->at(j),row->at(j));
				}
				else if(plane->at(j) == 2){
					hitMap[2]->Fill(col->at(j),row->at(j));
				}
				else if(plane->at(j) == 3){
					hitMap[3]->Fill(col->at(j),row->at(j));
				}
				else if(plane->at(j) == 4){
					hitMap[4]->Fill(col->at(j),row->at(j));
				}
				else if(plane->at(j) == 5){
					hitMap[5]->Fill(col->at(j),row->at(j));
				}
				else if(plane->at(j) == 6){
					hitMap[6]->Fill(col->at(j),row->at(j));
				}
				else{
					cout << "ERROR. PLANE DOES NOT EXIST"<<endl;
					break;
				}
			}
		}

		for(Int_t iROC = 0; iROC < 7; iROC++){
			Int_t numClusters = (Int_t)clust_per_plane->at(iROC);
			Ph1DHistogramsExtraction(numClusters, ph_1cl[iROC], ph_2cl[iROC], ph_3cl[iROC], ph_M4cl[iROC], phROC_all[iROC], phROC_1cl[iROC], phROC_2cl[iROC], phROC_3cl[iROC], phROC_M4cl[iROC]);
			Ph2DHistogramExtraction(numClusters, ph_1cl[iROC], clust_Local_X[iROC],clust_Local_Y[iROC],avPhROC_local_1cl[iROC],phROC_hitMap_local_1cl[iROC],deltaXR1,deltaYR1,kFALSE);
			Ph2DHistogramExtraction(numClusters, ph_1cl[iROC], clust_Telescope_X[iROC], clust_Telescope_Y[iROC], avPhROC_telescope_1cl[iROC], phROC_hitMap_telescope_1cl[iROC],deltaXR1, deltaYR1, kFALSE);
			Ph2DHistogramExtraction(numClusters, ph_1cl[iROC], clust_col[iROC],clust_row[iROC],avPhROC_pix_1cl[iROC],phROC_hitMap_pix_1cl[iROC],1,1,kTRUE);
		}
		Float_t progress;
		if(i == 0) cout << endl;
		else if(i%100 == 0){
			cout << "\r";
			//cout << "\e[A\r";
			progress = i*100/(Float_t)numEntries;
			cout << "Progress: " << setprecision(2) << setw(4) << setfill('0') << fixed << progress;
			cout << setw(1) << "%";
		}
	}

	// average for PH DUT
	FindAverageHistogramDUT(avPhROC_local_1cl[4],phROC_hitMap_local_1cl[4],avPhROC_local_1cl[5],phROC_hitMap_local_1cl[5],avPhROC_local_1cl[6],phROC_hitMap_local_1cl[6],divXR1,divYR1);
	FindAverageHistogramDUT(avPhROC_telescope_1cl[4],phROC_hitMap_telescope_1cl[4],avPhROC_telescope_1cl[5],phROC_hitMap_telescope_1cl[5],avPhROC_telescope_1cl[6],phROC_hitMap_telescope_1cl[6],divXR1,divYR1);
	FindAverageHistogramDUT(avPhROC_pix_1cl[4],phROC_hitMap_pix_1cl[4],avPhROC_pix_1cl[5],phROC_hitMap_pix_1cl[5],avPhROC_pix_1cl[6],phROC_hitMap_pix_1cl[6],numBinCol,numBinRow);

	//average for PH Telescope Planes
	FindAverageHistogramTPlanes(avPhROC_local_1cl[0],phROC_hitMap_local_1cl[0],avPhROC_local_1cl[1],phROC_hitMap_local_1cl[1],avPhROC_local_1cl[2],phROC_hitMap_local_1cl[2],avPhROC_local_1cl[3],phROC_hitMap_local_1cl[3],divXR1,divYR1);
	FindAverageHistogramTPlanes(avPhROC_telescope_1cl[0],phROC_hitMap_telescope_1cl[0],avPhROC_telescope_1cl[1],phROC_hitMap_telescope_1cl[1],avPhROC_telescope_1cl[2],phROC_hitMap_telescope_1cl[2],avPhROC_telescope_1cl[3],phROC_hitMap_telescope_1cl[3],divXR1,divYR1);
	FindAverageHistogramTPlanes(avPhROC_pix_1cl[0],phROC_hitMap_pix_1cl[0],avPhROC_pix_1cl[1],phROC_hitMap_pix_1cl[1],avPhROC_pix_1cl[2],phROC_hitMap_pix_1cl[2],avPhROC_pix_1cl[3],phROC_hitMap_pix_1cl[3],numBinCol,numBinRow);

	// Plot
	TCanvas *c0 = new TCanvas("c0","Coincidence Map",1);
	c0->cd();
	coincidenceMap->Draw("hist");
	c0->SaveAs(Form("%sCoincidenceMap.png",outputPath));
	delete c0;

	for(Int_t iROC=0;iROC<7;iROC++){
		if(iROC==4){
			// Histograms
			// 2D
			TString canvas_name_HitMapD1("Hit Map ROC4 Diamond 1");
			TString image_name_HitMapD1("HitMapROC4D1.png");
			PlotHistogram2D(hitMap[iROC],canvas_name_HitMapD1,outputPath,image_name_HitMapD1);
			TString canvas_name_avPHD1_local("Average PH ROC4 Diamond 1 Local Coord.");
			TString image_name_avPHD1_local("AvPHROC4D1Local.png");
			PlotHistogram2D(avPhROC_local_1cl[iROC],canvas_name_avPHD1_local,outputPath,image_name_avPHD1_local);
			TString canvas_name_avPHD1_telescope("Average PH ROC4 Diamond 1 telescope Coord.");
			TString image_name_avPHD1_telescope("AvPHROC4D1telescope.png");
			PlotHistogram2D(avPhROC_telescope_1cl[iROC],canvas_name_avPHD1_telescope,outputPath,image_name_avPHD1_telescope);
			TString canvas_name_avPHD1_pix("Average PH ROC4 Diamond 1 Pixelated Coord.");
			TString image_name_avPHD1_pix("AvPHROC4D1Pix.png");
			PlotHistogram2D(avPhROC_pix_1cl[iROC],canvas_name_avPHD1_pix,outputPath,image_name_avPHD1_pix);
			// 1D
			TString canvas_name_phD1("Pulse Height ROC4 Diamond 1");
			TString image_name_phD1("PhROC4D1.png");
			PlotPulseHeightsOverlay(phROC_all[iROC],phROC_1cl[iROC],phROC_2cl[iROC],phROC_3cl[iROC],phROC_M4cl[iROC],kFALSE,kFALSE,canvas_name_phD1,outputPath,image_name_phD1);
		}
		else if(iROC == 5){
			// 2D
			TString canvas_name_HitMapD2("Hit Map ROC5 Diamond 2");
			TString image_name_HitMapD2("HitMapROC5D2.png");
			PlotHistogram2D(hitMap[iROC],canvas_name_HitMapD2,outputPath,image_name_HitMapD2);
			TString canvas_name_avPHD2_local("Average PH ROC5 Diamond 2 Local Coord.");
			TString image_name_avPHD2_local("AvPHROC5D2Local.png");
			PlotHistogram2D(avPhROC_local_1cl[iROC],canvas_name_avPHD2_local,outputPath,image_name_avPHD2_local);
			TString canvas_name_avPHD2_telescope("Average PH ROC5 Diamond 2 telescope Coord.");
			TString image_name_avPHD2_telescope("AvPHROC5D2telescope.png");
			PlotHistogram2D(avPhROC_telescope_1cl[iROC],canvas_name_avPHD2_telescope,outputPath,image_name_avPHD2_telescope);
			TString canvas_name_avPHD2_pix("Average PH ROC5 Diamond 2 Pixelated Coord.");
			TString image_name_avPHD2_pix("AvPHROC5D2Pix.png");
			PlotHistogram2D(avPhROC_pix_1cl[iROC],canvas_name_avPHD2_pix,outputPath,image_name_avPHD2_pix);
			// 1D
			TString canvas_name_phD2("Pulse Height ROC5 Diamond 2");
			TString image_name_phD2("PhROC5D2.png");
			PlotPulseHeightsOverlay(phROC_all[iROC],phROC_1cl[iROC],phROC_2cl[iROC],phROC_3cl[iROC],phROC_M4cl[iROC],kFALSE,kFALSE,canvas_name_phD2,outputPath,image_name_phD2);
		}
		else if(iROC == 6){
			// 2D
			TString canvas_name_HitMapSi("Hit Map ROC6 Silicon");
			TString image_name_HitMapSi("HitMapROC6Si.png");
			PlotHistogram2D(hitMap[iROC],canvas_name_HitMapSi,outputPath,image_name_HitMapSi);
			TString canvas_name_avPHSi_local("Average PH ROC6 Silicon Local Coord.");
			TString image_name_avPHSi_local("AvPHROC6SiLocal.png");
			PlotHistogram2D(avPhROC_local_1cl[iROC],canvas_name_avPHSi_local,outputPath,image_name_avPHSi_local);
			TString canvas_name_avPHSi_telescope("Average PH ROC6 Silicon telescope Coord.");
			TString image_name_avPHSi_telescope("AvPHROC6Sitelescope.png");
			PlotHistogram2D(avPhROC_telescope_1cl[iROC],canvas_name_avPHSi_telescope,outputPath,image_name_avPHSi_telescope);
			TString canvas_name_avPHSi_pix("Average PH ROC6 Silicon Pixelated Coord.");
			TString image_name_avPHSi_pix("AvPHROC6SiPix.png");
			PlotHistogram2D(avPhROC_pix_1cl[iROC],canvas_name_avPHSi_pix,outputPath,image_name_avPHSi_pix);
			// 1D
			TString canvas_name_phSi("Pulse Height ROC6 Silicon");
			TString image_name_phSi("PhROC6Si.png");
			PlotPulseHeightsOverlay(phROC_all[iROC],phROC_1cl[iROC],phROC_2cl[iROC],phROC_3cl[iROC],phROC_M4cl[iROC],kFALSE,kFALSE,canvas_name_phSi,outputPath,image_name_phSi);
		}
		else{
			// 2D
			TString canvas_name_HitMap = TString::Format("Hit Map ROC%d",iROC);
			TString image_name_HitMap = TString::Format("HitMapROC%d.png",iROC);
			PlotHistogram2D(hitMap[iROC],canvas_name_HitMap,outputPath,image_name_HitMap);
			TString canvas_name_avPH_local = TString::Format("Average PH ROC%d Local Coord.",iROC);
			TString image_name_avPH_local = TString::Format("AvPHROC%dLocal.png",iROC);
			PlotHistogram2D(avPhROC_local_1cl[iROC],canvas_name_avPH_local,outputPath,image_name_avPH_local);
			TString canvas_name_avPH_telescope = TString::Format("Average PH ROC%d telescope Coord.",iROC);
			TString image_name_avPH_telescope = TString::Format("AvPHROC%dtelescope.png",iROC);
			PlotHistogram2D(avPhROC_telescope_1cl[iROC],canvas_name_avPH_telescope,outputPath,image_name_avPH_telescope);
			TString canvas_name_avPH_pix = TString::Format("Average PH ROC%d Pixelated Coord.",iROC);
			TString image_name_avPH_pix = TString::Format("AvPHROC%dPix.png",iROC);
			PlotHistogram2D(avPhROC_pix_1cl[iROC],canvas_name_avPH_pix,outputPath,image_name_avPH_pix);
			// 1D
			TString canvas_name_ph(Form("Pulse Height ROC%d",iROC));
			TString image_name_ph(Form("PhROC%d.png",iROC));
			PlotPulseHeightsOverlay(phROC_all[iROC],phROC_1cl[iROC],phROC_2cl[iROC],phROC_3cl[iROC],phROC_M4cl[iROC],kFALSE,kFALSE,canvas_name_ph,outputPath,image_name_ph);
		}
	}

	// Save histograms
	TFile f(Form("%s%s",outputPath,histogramOutputFile),"RECREATE");
	coincidenceMap->Write();
	for(Int_t iROC = 0; iROC<7;iROC++){
		hitMap[iROC]->Write();
		avPhROC_local_1cl[iROC]->Write();
		avPhROC_pix_1cl[iROC]->Write();
		phROC_all[iROC]->Write();
		phROC_1cl[iROC]->Write();
		phROC_2cl[iROC]->Write();
		phROC_3cl[iROC]->Write();
		phROC_M4cl[iROC]->Write();
	}

	f.Close();


}

void DoAveragePulseHeight(int numCl, vector<float> *ph1, vector<float> *ph2, vector<float> *ph3, vector<float> *phM4, TGraphErrors *phGall, TGraphErrors *phG1, TGraphErrors *phG2, TGraphErrors *phG3, TGraphErrors *phGM4){

}

void Ph1DHistogramsExtraction(int numCLROC, vector<float> *ph1, vector<float> *ph2, vector<float> *ph3, vector<float> *phM4, TH1F *phHall, TH1F *phH1, TH1F *phH2, TH1F *phH3, TH1F *phHM4){
	if(numCLROC != 0){
		for(Int_t i = 0; i < ph1->size(); i++){
			phH1->Fill(ph1->at(i));
			phHall->Fill(ph1->at(i));
		}
		for(Int_t i = 0; i < ph2->size(); i++){
			phH2->Fill(ph2->at(i));
			phHall->Fill(ph2->at(i));
		}
		for(Int_t i = 0; i < ph3->size(); i++){
			phH3->Fill(ph3->at(i));
			phHall->Fill(ph3->at(i));
		}
		for (Int_t i = 0; i < phM4->size(); i++) {
			phHM4->Fill(phM4->at(i));
			phHall->Fill(phM4->at(i));
		}
	}
}

void Ph2DHistogramExtraction(int numCLROC, vector<float> *phROC, vector<float> *clust_ROC_X, vector<float> *clust_ROC_Y, TH2F *avPHDUT, TH2F *hitMapPHDUT, float deltaX, float deltaY ,Bool_t isRowCol){
	if(numCLROC != 0){
		int numPhROC = phROC->size();
		Int_t tempBinX = 0, tempBinY = 0;
		Double_t tempPH = 0;
		for(Int_t i = 0; i < numPhROC; i++){
			if(clust_ROC_X->size() >= numPhROC && clust_ROC_Y->size() >= numPhROC){
				if(isRowCol){
					tempBinX = (Int_t)clust_ROC_X->at(i)+1;
					tempBinY = (Int_t)clust_ROC_Y->at(i)+1;
				}
				else{
					tempBinX = (Int_t)(TMath::CeilNint((Float_t)((clust_ROC_X->at(i)*10-xminR1)/(Float_t)deltaX)));
					tempBinY = (Int_t)(TMath::CeilNint((Float_t)((clust_ROC_Y->at(i)*10-yminR1)/(Float_t)deltaY)));
				}
				//if(!IsMasked(roc,tempBinX,tempBinY)){
					if( tempBinX >= 1 && tempBinY >= 1 && tempBinX <= divXR1 && tempBinY <= divYR1){
						tempPH = (Double_t)(avPHDUT->GetBinContent(tempBinX,tempBinY)+phROC->at(i));
						avPHDUT->SetBinContent(tempBinX,tempBinY,tempPH);
						if(isRowCol)
							hitMapPHDUT->Fill(clust_ROC_X->at(i),clust_ROC_Y->at(i));
						else
							hitMapPHDUT->Fill(10*clust_ROC_X->at(i),10*clust_ROC_Y->at(i));
					}
				//}
			}
		}
	}
}

void PlotPulseHeightsOverlay(TH1F *phall, TH1F *ph1, TH1F *ph2, TH1F *ph3, TH1F *phM4, Bool_t logY, Bool_t logX, const char *title, const char *outputPath, const char *suffix){
	gStyle->SetOptStat(0);
	TCanvas *c1 = new TCanvas("c1",title,1);
	c1->cd();
	if(logX) {c1->SetLogx();}
	if(logY) {c1->SetLogy();}
	phall->Draw("p e0");
	phall->Draw("hist same");
	ph1->Draw("p e0 same");
	ph1->Draw("hist same");
	ph2->Draw("p e0 same");
	ph2->Draw("hist same");
	ph3->Draw("p e0 same");
	ph3->Draw("hist same");
	phM4->Draw("p e0 same");
	phM4->Draw("hist same");
	c1->SaveAs(Form("%s%s",outputPath,suffix));
}

void PlotHistogram2D(TH2F *histogram, const char *title, const char *outputPath, const char *suffix){
	gStyle->SetOptStat(0);
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

void MaskPixelsFromFile(const string maskFile){
	cout << "Reading the mask file: " << maskFile << endl;
	vector<int> temproc, tempcol, temprow;
	std::ifstream inFile(maskFile.c_str());
	if(!inFile.is_open()){
		cerr << "ERROR: Cannot open mask file: " << maskFile << endl;
		throw;
	}
	for(std::string line; getline(inFile,line);){
		int row = 0, col = 0, roc=0;
		if(line[0] == '#' || line[0] == ' ' || line[0] == '/' || line[0] == '%') continue;
		std::istringstream linestream;
		linestream.str(line);
		linestream >> roc >> col >> row;
		temproc.push_back(roc);
		temprow.push_back(row);
		tempcol.push_back(col);
	}
	mskdroc = &temproc;
	mskdcol = &tempcol;
	mskdrow = &temprow;
}

Bool_t IsMasked(Int_t roc, Int_t col, Int_t row){
	if(mskdcol->empty()) return kFALSE;
	for(Int_t i = 0; i < mskdcol->size(); i++){
		if(mskdrow->at(i) == row && mskdcol->at(i) == col && mskdroc->at(i) == roc)
			return kTRUE;
	}
	return kFALSE;
}