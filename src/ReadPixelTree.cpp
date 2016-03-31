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

int main() {
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
	TBranch *braCharge = chain1->GetBranch("charge");
	TBranch *braDia1X = chain1->GetBranch("diam1_track_x");
	TBranch *braDia1Y = chain1->GetBranch("diam1_track_y");
	TBranch *braDia2X = chain1->GetBranch("diam2_track_x");
	TBranch *braDia2Y = chain1->GetBranch("diam2_track_y");
	//Create variables to read the tree
	vector<int> *plane=0, *col=0, *row=0, *adc=0, *charge=0;
	Float_t time=0, dia1X=0, dia1Y=0, dia2X=0, dia2Y=0;
	//dimension variables
	Int_t numPlane = 0, numCol = 0, numRow = 0, numADC = 0, numCharge = 0;
	Long64_t numEntries = braPlane->GetEntries();

	chain1->SetBranchAddress("time",&time);
	chain1->SetBranchAddress("plane",&plane);
	chain1->SetBranchAddress("col",&col);
	chain1->SetBranchAddress("row",&row);
	chain1->SetBranchAddress("adc",&adc);
	chain1->SetBranchAddress("charge",&charge);
	chain1->SetBranchAddress("diam1_track_x",&dia1X);
	chain1->SetBranchAddress("diam1_track_y",&dia1Y);
	chain1->SetBranchAddress("diam2_track_x",&dia2X);
	chain1->SetBranchAddress("diam2_track_y",&dia2Y);

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
	avADCDia1->SetContour(256);
	avADCDia2->SetContour(256);
	avADCSi->SetContour(256);
	Double_t tempADC = 0;
	hitPlane0->SetContour(256);
	hitPlane1->SetContour(256);
	hitPlane2->SetContour(256);
	hitPlane3->SetContour(256);
	hitPlane4->SetContour(256);
	hitPlane5->SetContour(256);
	hitPlane6->SetContour(256);

	TH1F *colP3 = new TH1F("colP3","colP3",61,-0.5,60.5);
	TH1F *rowP3 = new TH1F("rowP3","rowP3",91,-0.5,90.5);

	vector<int> columnas, filas, planos, adeces, cargas;

	for(Long64_t i = 0; i < numEntries; i++){
		chain1->GetEntry(i);
		numPlane = plane->size();
		numCol = col->size();
		numRow = row->size();
		numADC = adc->size();
		numCharge = charge->size();

		if(numPlane & numCol & numRow & numADC == numCharge){
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
	}
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
}
