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

const char *rootFilePathComplete = "~/Data/psi_2015_10/root/pixel/test151000312_withTracks.root";
const char *rootFilePathSmall = "~/DataSmall/psi_2015_10/root/pixel/TrackedRun312.root";

const char *outputPathComplete = "~/Dropbox/201601/DiamondPixels/ReadPixelTree/";
const char *outputPathSmall = "~/Dropbox/201601/DiamondPixels/ReadPixelTree/Small/";

Int_t contour2D = 256;
Float_t xminR1 = -6.075, xminR2 = -6.0415, xmaxR1 = 6.075, xmaxR2 = 6.0415, yminR1 = -6.05, yminR2 = -6.0175, ymaxR1 = 6.05, ymaxR2 = 6.0175;
Int_t divXR1 = 81, divXR2 = 281, divYR1 = 121, divYR2 = 415;
Float_t deltaXR1 = (Float_t)((xmaxR1-xminR1)/divXR1), deltaXR2 = (Float_t)((xmaxR2-xminR2)/divXR2), deltaYR1 = (Float_t)((ymaxR1-yminR1)/divYR1), deltaYR2 = (Float_t)((ymaxR2-yminR2)/divYR2);

int main() {
	gROOT->ProcessLine("#include <vector>");
	gSystem->Load("libMathCore");
	gSystem->Load("libPhysics");
	const char *rootFilePath = rootFilePathSmall;
	const char *outputPath = outputPathSmall;
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
	TBranch *braTime = chain1->GetBranch("time");
	TBranch *braPlane = chain1->GetBranch("plane");
	TBranch *braCol = chain1->GetBranch("col");
	TBranch *braRow = chain1->GetBranch("row");
	TBranch *braADC = chain1->GetBranch("adc");
	TBranch *braDia1X = chain1->GetBranch("diam1_track_x");
	TBranch *braDia1Y = chain1->GetBranch("diam1_track_y");
	TBranch *braDia2X = chain1->GetBranch("diam2_track_x");
	TBranch *braDia2Y = chain1->GetBranch("diam2_track_y");
	//Create variables to read the tree
	vector<int> *plane=0, *col=0, *row=0, *adc=0;
	vector<unsigned char> *clust_per_plane=0;
	vector<float> *clust_ROC0_X = 0, *clust_ROC0_Y = 0, *ph_ROC0_1cl = 0, *clust_ROC1_X = 0, *clust_ROC1_Y = 0, *ph_ROC1_1cl = 0;
	vector<float> *clust_ROC2_X = 0, *clust_ROC2_Y = 0, *ph_ROC2_1cl = 0, *clust_ROC3_X = 0, *clust_ROC3_Y = 0, *ph_ROC3_1cl = 0;
	vector<float> *clust_ROC4_X = 0, *clust_ROC4_Y = 0, *ph_ROC4_1cl = 0, *clust_ROC5_X = 0, *clust_ROC5_Y = 0, *ph_ROC5_1cl = 0, *clust_ROC6_X = 0, *clust_ROC6_Y = 0, *ph_ROC6_1cl = 0;
	Float_t time=0, dia1X=0, dia1Y=0, dia2X=0, dia2Y=0, chi2_tracks=0, chi2_x=0, chi2_y=0, slope_x=0, slope_y=0;
	Long64_t event_number=0;
	//Int_t n_tracks = 0, n_cluster = 0;
	//dimension variables
	vector<int> numCLROC;
	numCLROC.resize(7);
	Int_t numPlane = 0, numCol = 0, numRow = 0, numADC = 0, numPhR0_1cl = 0, numCLR1X = 0, numCLR1Y = 0, numPhR1_1cl = 0, numCLR2X = 0, numCLR2Y = 0, numPhR2_1cl = 0;
	Int_t numCLR3X = 0, numCLR3Y = 0, numPhR3_1cl = 0, numCLR4X = 0, numCLR4Y = 0, numPhR4_1cl = 0, numCLR5X = 0, numCLR5Y = 0, numPhR5_1cl = 0, numCLR6X = 0, numCLR6Y = 0, numPhR6_1cl = 0;
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

	TH2F *avPHDia1Res1 = new TH2F("avPHDia1Res1","Average PH Diamond 1 R1",divXR1,xminR1,xmaxR1,divYR1,yminR1,ymaxR1);
	TH2F *avPHDia1Res2 = new TH2F("avPHDia1Res2","Average PH Diamond 1 R2",divXR2,xminR2,xmaxR2,divYR2,yminR2,ymaxR2);
	TH1F *phDia1 = new TH1F("phDia1","Pulse Height Diamond 1",501,-250,250250);
	TH2F *avPHDia2Res1 = new TH2F("avPHDia2Res1","Average PH Diamond 2 R1",divXR1,xminR1,xmaxR1,divYR1,yminR1,ymaxR1);
	TH2F *avPHDia2Res2 = new TH2F("avPHDia2Res2","Average PH Diamond 2 R2",divXR2,xminR2,xmaxR2,divYR2,yminR2,ymaxR2);
	TH1F *phDia2 = new TH1F("phDia2","Pulse Height Diamond 2",501,-2052,50052);
	TH2F *avPHSiRes1 = new TH2F("avPHSiRes1","Average PH Silicon R1",divXR1,xminR1,xmaxR1,divYR1,yminR1,ymaxR1);
	TH2F *avPHSiRes2 = new TH2F("avPHSiRes2","Average PH Silicon R2",divXR2,xminR2,xmaxR2,divYR2,yminR2,ymaxR2);
	TH1F *phSi = new TH1F("phSi","Pulse Height Silicon",501,-10210,200210);

	TH2F *hitMapPHDia1Res1 = new TH2F("hitMapPHDia1Res1","Hit Map Diamond 1 for PH R1",divXR1,xminR1,xmaxR1,divYR1,yminR1,ymaxR1);
	TH2F *hitMapPHDia1Res2 = new TH2F("hitMapPHDia1Res2","Hit Map Diamond 1 for PH R2",divXR2,xminR2,xmaxR2,divYR2,yminR2,ymaxR2);
	TH2F *hitMapPHDia2Res1 = new TH2F("hitMapPHDia2Res1","Hit Map Diamond 2 for PH R1",divXR1,xminR1,xmaxR1,divYR1,yminR1,ymaxR1);
	TH2F *hitMapPHDia2Res2 = new TH2F("hitMapPHDia2Res2","Hit Map Diamond 2 for PH R2",divXR2,xminR2,xmaxR2,divYR2,yminR2,ymaxR2);
	TH2F *hitMapPHSiRes1 = new TH2F("hitMapPHSiRes1","Hit Map Si for PH R1",divXR1,xminR1,xmaxR1,divYR1,yminR1,ymaxR1);
	TH2F *hitMapPHSiRes2 = new TH2F("hitMapPHSiRes2","Hit Map Si for PH R2",divXR2,xminR2,xmaxR2,divYR2,yminR2,ymaxR2);


	avADCDia1->SetContour(contour2D);
	avADCDia2->SetContour(contour2D);
	avADCSi->SetContour(contour2D);
	Double_t tempADC = 0;
	hitPlane0->SetContour(contour2D);
	hitPlane1->SetContour(contour2D);
	hitPlane2->SetContour(contour2D);
	hitPlane3->SetContour(contour2D);
	hitPlane4->SetContour(contour2D);
	hitPlane5->SetContour(contour2D);
	hitPlane6->SetContour(contour2D);

	Double_t tempPH = 0;
	Int_t tempBinX = 0, tempBinY = 0;
	avPHDia1Res1->SetContour(contour2D);
	avPHDia1Res2->SetContour(contour2D);
	avPHDia2Res1->SetContour(contour2D);
	avPHDia2Res2->SetContour(contour2D);
	avPHSiRes1->SetContour(contour2D);
	avPHSiRes2->SetContour(contour2D);

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
		for(Int_t i = 4; i < 7; i++){
			numCLROC[i] = (Int_t)clust_per_plane->at(i);
		}
		// For Diamond 1
		if(numCLROC[4] != 0){
			for(Int_t i = 0; i < numPhR4_1cl; i++){
				// For resolution 1
				if(clust_ROC4_X->size() >= numPhR4_1cl && clust_ROC4_Y->size() >= numPhR4_1cl){
					tempBinX = (Int_t)(TMath::Ceil((Float_t)((clust_ROC4_X->at(i)*10-xminR1)/deltaXR1)));
					tempBinY = (Int_t)(TMath::Ceil((Float_t)((clust_ROC4_Y->at(i)*10-yminR1)/deltaYR1)));
					phDia1->Fill(ph_ROC4_1cl->at(i));
					if( tempBinX >= 1 && tempBinY >= 1 && tempBinX <= divXR1 && tempBinY <= divYR1){
						tempPH = (Double_t)(avPHDia1Res1->GetBinContent(tempBinX,tempBinY)+ph_ROC4_1cl->at(i));
						avPHDia1Res1->SetBinContent(tempBinX,tempBinY,tempPH);
						hitMapPHDia1Res1->Fill(tempBinX,tempBinY);
					}
					// For resolution 2
					tempBinX = (Int_t)(TMath::Ceil((Float_t)((clust_ROC4_X->at(i)*10-xminR2)/deltaXR2)));
					tempBinY = (Int_t)(TMath::Ceil((Float_t)((clust_ROC4_Y->at(i)*10-yminR2)/deltaYR2)));
					if( tempBinX >= 1 && tempBinY >= 1 && tempBinX <= divXR2 && tempBinY <= divYR2){
						tempPH = (Double_t)(avPHDia1Res2->GetBinContent(tempBinX,tempBinY)+ph_ROC4_1cl->at(i));
						avPHDia1Res2->SetBinContent(tempBinX,tempBinY,tempPH);
						hitMapPHDia1Res2->Fill(tempBinX,tempBinY);
					}
				}
			}

		}
		// For Diamond 2
		if(numCLROC[5] != 0){
			numPhR5_1cl = ph_ROC5_1cl->size();
			for(Int_t i = 0; i < numPhR5_1cl; i++){
				// For resolution 1
				if(clust_ROC5_X->size() >= numPhR5_1cl && clust_ROC5_Y->size() >= numPhR5_1cl){
					tempBinX = (Int_t)(TMath::Ceil((Float_t)((clust_ROC5_X->at(i)*10-xminR1)/deltaXR1)));
					tempBinY = (Int_t)(TMath::Ceil((Float_t)((clust_ROC5_Y->at(i)*10-yminR1)/deltaYR1)));
					phDia2->Fill(ph_ROC5_1cl->at(i));
					if( tempBinX >= 1 && tempBinY >= 1 && tempBinX <= divXR1 && tempBinY <= divYR1){
						tempPH = (Double_t)(avPHDia2Res1->GetBinContent(tempBinX,tempBinY)+ph_ROC5_1cl->at(i));
						avPHDia2Res1->SetBinContent(tempBinX,tempBinY,tempPH);
						hitMapPHDia2Res1->Fill(tempBinX,tempBinY);
					}
					// For resolution 2
					tempBinX = (Int_t)(TMath::Ceil((Float_t)((clust_ROC5_X->at(i)*10-xminR2)/deltaXR2)));
					tempBinY = (Int_t)(TMath::Ceil((Float_t)((clust_ROC5_Y->at(i)*10-yminR2)/deltaYR2)));
					if( tempBinX >= 1 && tempBinY >= 1 && tempBinX <= divXR2 && tempBinY <= divYR2){
						tempPH = (Double_t)(avPHDia2Res2->GetBinContent(tempBinX,tempBinY)+ph_ROC5_1cl->at(i));
						avPHDia2Res2->SetBinContent(tempBinX,tempBinY,tempPH);
						hitMapPHDia2Res2->Fill(tempBinX,tempBinY);
					}
				}
			}

		}
		// For Silicon
		if(numCLROC[6] != 0){
			numPhR6_1cl = ph_ROC6_1cl->size();
			for(Int_t i = 0; i < numPhR6_1cl; i++){
				// For resolution 1
				if(clust_ROC6_X->size() >= numPhR6_1cl && clust_ROC6_Y->size() >= numPhR6_1cl){
					tempBinX = (Int_t)(TMath::Ceil((Float_t)((clust_ROC6_X->at(i)*10-xminR1)/deltaXR1)));
					tempBinY = (Int_t)(TMath::Ceil((Float_t)((clust_ROC6_Y->at(i)*10-yminR1)/deltaYR1)));
					phSi->Fill(ph_ROC6_1cl->at(i));
					if( tempBinX >= 1 && tempBinY >= 1 && tempBinX <= divXR1 && tempBinY <= divYR1){
						tempPH = (Double_t)(avPHSiRes1->GetBinContent(tempBinX,tempBinY)+ph_ROC6_1cl->at(i));
						avPHSiRes1->SetBinContent(tempBinX,tempBinY,tempPH);
						hitMapPHSiRes1->Fill(tempBinX,tempBinY);
					}
					// For resolution 2
					tempBinX = (Int_t)(TMath::Ceil((Float_t)((clust_ROC6_X->at(i)*10-xminR2)/deltaXR2)));
					tempBinY = (Int_t)(TMath::Ceil((Float_t)((clust_ROC6_Y->at(i)*10-yminR2)/deltaYR2)));
					if( tempBinX >= 1 && tempBinY >= 1 && tempBinX <= divXR2 && tempBinY <= divYR2){
						tempPH = (Double_t)(avPHSiRes2->GetBinContent(tempBinX,tempBinY)+ph_ROC6_1cl->at(i));
						avPHSiRes2->SetBinContent(tempBinX,tempBinY,tempPH);
						hitMapPHSiRes2->Fill(tempBinX,tempBinY);
					}
				}
			}
		}
	}

	// average for ADC
	for(Int_t ii = 1; ii <= 52; ii++){
		for(Int_t jj = 1; jj <= 80; jj++){
			if(hitPlane4->GetBinContent(ii,jj) >= 1){
				tempADC = (Double_t)(avADCDia1->GetBinContent(ii,jj)/(Double_t)hitPlane4->GetBinContent(ii,jj));
				avADCDia1->SetBinContent(ii,jj,tempADC);

			}
			if(hitPlane5->GetBinContent(ii,jj) >= 1){
				tempADC = (Double_t)(avADCDia2->GetBinContent(ii,jj)/(Double_t)hitPlane5->GetBinContent(ii,jj));
				avADCDia2->SetBinContent(ii,jj,tempADC);
			}
			if(hitPlane6->GetBinContent(ii,jj) >= 1){
				tempADC = (Double_t)(avADCSi->GetBinContent(ii,jj)/(Double_t)hitPlane6->GetBinContent(ii,jj));
				avADCSi->SetBinContent(ii,jj,tempADC);
			}
		}
	}
	// average for PH
	// Res1
	for(Int_t i = 1; i <= divXR1; i++){
		for(Int_t j = 1; j <= divYR1; j++){
			if(hitMapPHDia1Res1->GetBinContent(i,j) >= 1){
				tempPH = (Double_t)(avPHDia1Res1->GetBinContent(i,j)/(Double_t)hitMapPHDia1Res1->GetBinContent(i,j));
				avPHDia1Res1->SetBinContent(i,j,tempPH);
			}
			if(hitMapPHDia2Res1->GetBinContent(i,j) >= 1){
				tempPH = (Double_t)(avPHDia2Res1->GetBinContent(i,j)/(Double_t)hitMapPHDia2Res1->GetBinContent(i,j));
				avPHDia2Res1->SetBinContent(i,j,tempPH);
			}
			if(hitMapPHSiRes1->GetBinContent(i,j) >= 1){
				tempPH = (Double_t)(avPHSiRes1->GetBinContent(i,j)/(Double_t)hitMapPHSiRes1->GetBinContent(i,j));
				avPHSiRes1->SetBinContent(i,j,tempPH);
			}
		}
	}
	// Res2
	for(Int_t i = 1; i <= divXR2; i++){
		for(Int_t j = 1; j <= divYR2; j++){
			if(hitMapPHDia1Res2->GetBinContent(i,j) >= 1){
				tempPH = (Double_t)(avPHDia1Res2->GetBinContent(i,j)/(Double_t)hitMapPHDia1Res2->GetBinContent(i,j));
				avPHDia1Res2->SetBinContent(i,j,tempPH);
			}
			if(hitMapPHDia2Res2->GetBinContent(i,j) >= 1){
				tempPH = (Double_t)(avPHDia2Res2->GetBinContent(i,j)/(Double_t)hitMapPHDia2Res2->GetBinContent(i,j));
				avPHDia2Res2->SetBinContent(i,j,tempPH);
			}
			if(hitMapPHSiRes2->GetBinContent(i,j) >= 1){
				tempPH = (Double_t)(avPHSiRes2->GetBinContent(i,j)/(Double_t)hitMapPHSiRes2->GetBinContent(i,j));
				avPHSiRes2->SetBinContent(i,j,tempPH);
			}
		}
	}
	
	TCanvas *c1 = new TCanvas("c1","Hit map Plane 0",1);
	c1->cd();
	hitPlane0->Draw("colz");
	c1->SaveAs(Form("%sPlane0.png",outputPath));
	TCanvas *c2 = new TCanvas("c2","Hit map Plane 1",1);
	c2->cd();
	hitPlane1->Draw("colz");
	c2->SaveAs(Form("%sPlane1.png",outputPath));
	TCanvas *c3 = new TCanvas("c3","Hit map Plane 2",1);
	c3->cd();
	hitPlane2->Draw("colz");
	c3->SaveAs(Form("%sPlane2.png",outputPath));
	TCanvas *c4 = new TCanvas("c4","Hit map Plane 3",1);
	c4->cd();
	hitPlane3->Draw("colz");
	c4->SaveAs(Form("%sPlane3.png",outputPath));
	TCanvas *c5 = new TCanvas("c5","Hit map Plane 4",1);
	c5->cd();
	hitPlane4->Draw("colz");
	c5->SaveAs(Form("%sPlane4.png",outputPath));
	TCanvas *c6 = new TCanvas("c6","Hit map Plane 5",1);
	c6->cd();
	hitPlane5->Draw("colz");
	c6->SaveAs(Form("%sPlane5.png",outputPath));
	TCanvas *c7 = new TCanvas("c7","Hit map Plane 6",1);
	c7->cd();
	hitPlane6->Draw("colz");
	c7->SaveAs(Form("%sPlane6.png",outputPath));
	TCanvas *c8 = new TCanvas("c8","ADC Average Diamond 1",1);
	c8->cd();
	avADCDia1->Draw("colz");
	c8->SaveAs(Form("%sAvADCD1.png",outputPath));
	TCanvas *c9 = new TCanvas("c9","ADC Average Diamond 2",1);
	c9->cd();
	avADCDia2->Draw("colz");
	c9->SaveAs(Form("%sAvADCD2.png",outputPath));
	TCanvas *c12 = new TCanvas("c12","ADC Average Si",1);
	c12->cd();
	avADCSi->Draw("colz");
	c12->SaveAs(Form("%sAvADCSi.png",outputPath));
	TCanvas *c10 = new TCanvas("c10","ADC Diamond 1",1);
	c10->cd();
	c10->SetLogy();
	adcDia1->Draw();
	c10->SaveAs(Form("%sADCD1.png",outputPath));
	TCanvas *c11 = new TCanvas("c11","ADC Diamond 2",1);
	c11->cd();
	c11->SetLogy();
	adcDia1->Draw();
	c11->SaveAs(Form("%sADCD2.png",outputPath));
	TCanvas *c13 = new TCanvas("c13","ADC Si",1);
	c13->cd();
	c13->SetLogy();
	adcSi->Draw();
	c13->SaveAs(Form("%sADCSi.png",outputPath));
	TCanvas *c14 = new TCanvas("c14","PH Average Diamond 1 R1",1);
	c14->cd();
	avPHDia1Res1->Draw("colz");
	c14->SaveAs(Form("%sAvPHD1R1.png",outputPath));
	TCanvas *c15 = new TCanvas("c15","PH Average Diamond 1 R2",1);
	c15->cd();
	avPHDia1Res2->Draw("colz");
	c15->SaveAs(Form("%sAvPHD1R2.png",outputPath));
	TCanvas *c16 = new TCanvas("c16","PH Average Diamond 2 R1",1);
	c16->cd();
	avPHDia2Res1->Draw("colz");
	c16->SaveAs(Form("%sAvPHD2R1.png",outputPath));
	TCanvas *c17 = new TCanvas("c17","PH Average Diamond 2 R2",1);
	c17->cd();
	avPHDia2Res2->Draw("colz");
	c17->SaveAs(Form("%sAvPHD2R2.png",outputPath));
	TCanvas *c18 = new TCanvas("c18","PH Average Silicon R1",1);
	c18->cd();
	avPHSiRes1->Draw("colz");
	c18->SaveAs(Form("%sAvPHSiR1.png",outputPath));
	TCanvas *c19 = new TCanvas("c19","PH Average Silicon R2",1);
	c19->cd();
	avPHSiRes2->Draw("colz");
	c19->SaveAs(Form("%sAvPHSiR2.png",outputPath));
	TCanvas *c20 = new TCanvas("c20","PH Diamond 1",1);
	c20->cd();
	c20->SetLogy();
	phDia1->Draw("B E1");
	c20->SaveAs(Form("%sPHD1.png",outputPath));
	TCanvas *c21 = new TCanvas("c21","PH Diamond 2",1);
	c21->cd();
	c21->SetLogy();
	phDia2->Draw("B E1");
	c21->SaveAs(Form("%sPHD2.png",outputPath));
	TCanvas *c22 = new TCanvas("c22","PH Silicon",1);
	c22->cd();
	c22->SetLogy();
	phSi->Draw("B E1");
	c22->SaveAs(Form("%sPHSi.png",outputPath));
}
