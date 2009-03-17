/*
// This reads in the output file from simpleGeantSim and calculates the point
// of closest approach for those muons that hit the top and bottom planes
// of the scintillator
//
// This is based on StepThrough
//
// Ryan Nichol <rjn@hep.ucl.ac.uk>
*/


// Headers from simpleGeantSim
#include "DetectorDefs.hh"

// root headers
#include <TH3F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TStyle.h>
#include <TSystem.h>

// std headers
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>


//The steering file
#define SOURCE_F "sourcePca.dat"
#define SMALL_NUM  0.00000001

void CloseApproach(double xz[2], double yz[2], double inxz[2], double inyz[2], double ca[3])
{ 
  Double_t a=xz[0]*xz[0]+yz[0]*yz[0]+1;
  Double_t b=xz[0]*xz[1]+yz[0]*yz[1]+1;
  Double_t c=xz[1]*xz[1]+yz[1]*yz[1]+1;
  Double_t d=xz[0]*(inxz[0]-inxz[1]) + yz[0]*(inyz[0]-inyz[1]);
  Double_t e=xz[1]*(inxz[0]-inxz[1]) + yz[1]*(inyz[0]-inyz[1]);

  Double_t denom=a*c-b*b;
  Double_t sc=(a*e - c*d)/denom;
  Double_t tc=(a*e - c*d)/denom;

  //  std::cout << "CA:\t" << denom << "\t" << sc << "\t" << tc << "\n";
  ca[0]=0.5*((xz[0]*sc +inxz[0])+(xz[1]*tc+inxz[1]));
  ca[1]=0.5*((yz[0]*sc +inyz[0])+(yz[1]*tc+inyz[1]));
  ca[2]=0.5*(sc+tc);
}



void fitxyz(int nPoints, float X[], float Y[], float Z[], double gradXZ[2], double gradYZ[2])
{

  double sumZ2 = 0.0, sumX = 0.0, sumXZ = 0.0, sumZ = 0.0;

  double sumY2 = 0.0, sumY = 0.0, sumYZ = 0.0;
  
  double sumX2 = 0.0;

  double deltaX = 0.0, deltaY = 0.0;


  for (int i = 0; i < nPoints; i++)
    { 
      sumZ2 += double( Z[i]*Z[i] );       
      sumXZ += double( X[i]*Z[i] );
      sumYZ += double( Y[i]*Z[i] );
      sumX  += double( X[i] );       
      sumY  += double( Y[i] );       
      sumZ  += double( Z[i] );       
    }

  sumZ2/=nPoints;
  sumXZ/=nPoints;
  sumYZ/=nPoints;
  sumX/=nPoints;
  sumY/=nPoints;
  sumZ/=nPoints;
  
  double covxz=sumXZ-sumX*sumZ;
  double covyz=sumYZ-sumY*sumZ;
  double sigmaz2=sumZ2-sumZ*sumZ;

  gradXZ[0] = covxz/sigmaz2;
  gradXZ[1] = sumX - gradXZ[0]*sumZ;

  gradYZ[0] = covyz/sigmaz2;
  gradYZ[1] = sumY - gradYZ[0]*sumZ;


}



#define MAX_SCINT_HITS PLANES_PER_SIDE*2*10 //Just guess for now


int main(int argc, char**argv)
{
  float histoRange = SIDELENGTH*1000;

  int noBins = 100; //default value
  char inputFile[256];
  char outputFile[256];

  ifstream in0;
  in0.open(SOURCE_F);
  in0>>noBins>>inputFile;

  //  std::string outName = "pca_";
  //  outName.append (inputFile);

  //Note the below code is lazy and causes a memory leak
  sprintf(outputFile,"%s/pca_%s",gSystem->DirName(inputFile),gSystem->BaseName(inputFile));

  TFile* newfile = new TFile(outputFile,"RECREATE");
  


  double xPos = 0, yPos = 0, zPos = 0;

  TTree* PCA = new TTree ("PCA","PCA");
  
  PCA->Branch ("xPos", &xPos,"xPos/D");
  PCA->Branch ("yPos", &yPos,"yPos/D");
  PCA->Branch ("zPos", &zPos,"zPos/D");
 
  TH3F* PCAh = new TH3F("PCAh","PCAh",noBins,-histoRange/2,histoRange/2,  noBins,-histoRange/2,histoRange/2, noBins,-histoRange/2,histoRange/2); 
  
  double AGradx = 0, AGrady = 0, ACutx = 0, ACuty = 0; //cuts and gradients for y as a function of x or z  (ie y = a*x + b)

  TTree* Absorbed = new TTree("Absorbed","Absorbed");

  Absorbed->Branch("xGrad", &AGradx, "gradX/D" );
  Absorbed->Branch("xCut",  &ACutx,  "cutX/D" );
  Absorbed->Branch("yGrad", &AGrady, "gradY/D" );
  Absorbed->Branch("yCut",  &ACuty,  "cutY/D" );

  TH2F* AbsorbedHisto = new TH2F("AbsorbedH","AbsorbedH", noBins, -histoRange/2, histoRange/2, noBins, -histoRange/2, histoRange/2);

  TH1F* energy = new TH1F("energy","energy",100,0,100);
 
  std::cout<<inputFile<<std::endl;

  TFile* fp = new TFile(inputFile,"READ");
  in0.close();

  Int_t fRun,fEvent;
  Int_t numScintHits;
  Int_t side[MAX_SCINT_HITS];
  Int_t plane[MAX_SCINT_HITS];
  Int_t strip[MAX_SCINT_HITS];
  Double_t truePos[MAX_SCINT_HITS][3];
  Double_t energyDep[MAX_SCINT_HITS];

  TTree* scintTree = (TTree*) fp->Get("scintTree");
  scintTree->SetMakeClass(1);
  scintTree->SetBranchAddress("run",&fRun);
  scintTree->SetBranchAddress("event",&fEvent);
  scintTree->SetBranchAddress("ScintHitInfo",&numScintHits); 
  scintTree->SetBranchAddress("ScintHitInfo.side",side); 
  scintTree->SetBranchAddress("ScintHitInfo.plane",plane); 
  scintTree->SetBranchAddress("ScintHitInfo.strip",strip); 
  scintTree->SetBranchAddress("ScintHitInfo.truePos[3]",truePos); 
  scintTree->SetBranchAddress("ScintHitInfo.energyDep",energyDep); 

  float xTop[PLANES_PER_SIDE],yTop[PLANES_PER_SIDE],zTop[PLANES_PER_SIDE];
  float xBot[PLANES_PER_SIDE],yBot[PLANES_PER_SIDE],zBot[PLANES_PER_SIDE];

  double xzGrad[2] = {0.0}; // gradient for z = ax + b
  double yzGrad[2] = {0.0}; // as above in yz plane
  double xzCut[2]  = {0.0};  // axis intercept xz plane (b)
  double yzCut[2]  = {0.0};  // as above yz plane

  int nLines = 0;
  double gradxz[2];
  double gradyz[2];
  

  int nEntries = scintTree->GetEntries();
  //  nEntries=100;
  std::cout << "There are " << nEntries << " entries\n";
 
  for(int entry=0;  entry < nEntries; entry++)
    {
      if( !(entry%1000) ) std::cout<< entry << " of " << nEntries <<std::endl;
      scintTree->GetEntry(entry);


      if(numScintHits>3) {
	//	std::cout << "Num hits:\t" << numScintHits << "\n";
	//Minimum requirement
	int gotPlane[PLANES_PER_SIDE*2]={0};
	Double_t xTemp[PLANES_PER_SIDE*2]={0},yTemp[PLANES_PER_SIDE*2]={0},zTemp[PLANES_PER_SIDE*2]={0};
	Double_t weights[PLANES_PER_SIDE*2]={0};
	for(int hit=0;hit<numScintHits;hit++) {
	  //	  std::cout << plane[hit] << "\t" << strip[hit] << "\t" << truePos[hit][2] << "\t" << energyDep[hit] << "\n";
	  gotPlane[plane[hit]]=1;
	  xTemp[plane[hit]]+=truePos[hit][0]*energyDep[hit];
	  yTemp[plane[hit]]+=truePos[hit][1]*energyDep[hit];
	  zTemp[plane[hit]]+=truePos[hit][2]*energyDep[hit];
	  weights[plane[hit]]+=energyDep[hit];
	}	
	
	int topPlanes=0;
	int botPlanes=0;
	for(int plane=0;plane<PLANES_PER_SIDE;plane++) {
	  
	  if(gotPlane[plane]) {
	    xTemp[plane]/=weights[plane];
	    yTemp[plane]/=weights[plane];
	    zTemp[plane]/=weights[plane];
	    
	    //	    std::cout << plane << "\t" << xTemp[plane] << "\t" << yTemp[plane] << "\t" << zTemp[plane] << "\n";

	    xTop[topPlanes]=xTemp[plane];
	    yTop[topPlanes]=yTemp[plane];
	    zTop[topPlanes]=zTemp[plane];
	    topPlanes++;

	  }
	  if(gotPlane[plane+PLANES_PER_SIDE]) {
	    xTemp[plane+PLANES_PER_SIDE]/=weights[plane+PLANES_PER_SIDE];
	    yTemp[plane+PLANES_PER_SIDE]/=weights[plane+PLANES_PER_SIDE];
	    zTemp[plane+PLANES_PER_SIDE]/=weights[plane+PLANES_PER_SIDE];

	    xBot[botPlanes]=xTemp[plane+PLANES_PER_SIDE];
	    yBot[botPlanes]=yTemp[plane+PLANES_PER_SIDE];
	    zBot[botPlanes]=zTemp[plane+PLANES_PER_SIDE];
	    botPlanes++;
	  }

	}

	//	std::cout << "Planes:\t" << topPlanes << "\t" << botPlanes << "\n";
	if(topPlanes<3) continue; //Didn't hit the top
	fitxyz(topPlanes,xTop,yTop,zTop,gradxz,gradyz);
	xzGrad[0]=gradxz[0];
	yzGrad[0]=gradyz[0];
	xzCut[0]=gradxz[1];
	yzCut[0]=gradyz[1];
	
	//	std::cout << gradxz[0] << "\t" << gradxz[1] << "\t" << gradxz[1] + 1000*gradxz[0] << "\n";;
	//	std::cout << gradyz[0] << "\t" << gradyz[1] << "\t" << gradyz[1] + 1000*gradyz[0] << "\n";;
	 

	if(botPlanes>2) {
	  //Potential scatter
	  fitxyz(botPlanes,xBot,yBot,zBot,gradxz,gradyz);
	  xzGrad[1]=gradxz[0];
	  yzGrad[1]=gradyz[0];
	  xzCut[1]=gradxz[1];
	  yzCut[1]=gradyz[1];

	  Double_t pca[3];
	  CloseApproach(xzGrad,yzGrad,xzCut,yzCut,pca);
// 	  for(double z=-2000;z<2000;z+=100) {
// 	    double x1=xzGrad[0]*z + xzCut[0];
// 	    double x2=xzGrad[1]*z + xzCut[1];
// 	    double y1=yzGrad[0]*z + yzCut[0];
// 	    double y2=yzGrad[1]*z + yzCut[1];
// 	    double diff=(x2-x1)*(x2-x1) + (y2-y1)*(y2-y1);
// 	    std::cout << "Scan: " << z << "\t" << diff << "\n";
// 	  }

//	  std::cout << pca[0] << "\t" << pca[1] << "\t" << pca[2] << "\n";
	  xPos = pca[0];
	  yPos = pca[1];
	  zPos = pca[2];

	  PCA->Fill ();

	  PCAh->Fill(pca[0],pca[1],pca[2]);
	  
	
	}
	else {
	  //Potential Absorbtion ... need to tweak them
 	  AbsorbedHisto->Fill(xTop[0],yTop[0]);

	  AGradx=gradxz[0];
	  AGrady=gradyz[0];
	  ACutx=gradxz[1];
	  ACuty=gradyz[1];
 	  Absorbed ->Fill();
	}

      }

    }

  newfile->cd();

  newfile->Write();

}


