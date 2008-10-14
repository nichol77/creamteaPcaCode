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

void CloseApproach(double xy[2], double zy[2], double inxy[2], double inzy[2], double ca[3])
{ 
  double A1 = xy[0], B1 = inxy[0], C1 = zy[0], D1 = inzy[0];

  double A2 = xy[1], B2 = inxy[1], C2 = zy[1], D2 = inzy[1];

  double M11 = 1.0 + (A1*A1) + (A1*A1)/(C1*C1);
  double M12 =-1.0 - (A1*A2) - (A2*A1)/(C1*C2);
  double M22 = 1.0 + (A2*A2) + (A2*A2)/(C2*C2);

  double K1 = (A1*B1) - (A1*B2) + (A1*B1)/(C1*C1) - (A1*D1)/(C1*C1) - (A1*B2)/(C1*C2) + (D2*A1)/(C1*C2);  
  double K2 = (A2*B2) - (A2*B1) + (A2*D1)/(C1*C2) - (A2*B1)/(C1*C2) - (A2*D2)/(C2*C2) + (A2*B2)/(C2*C2);
 
  double xinter1 = (K1*M22 - K2*M12)/(M12*M12 - M11*M22);
  double xinter2 = (K1*M12 - K2*M11)/(M11*M22 - M12*M12);

  double yinter2 = (A2*xinter2) + B2;
  double yinter1 = (A1*xinter1) + B1;

  double zinter2 = (A2*xinter2 + B2 - D2)/C2;
  double zinter1 = (A1*xinter1 + B1 - D1)/C1;

  ca[0]= (xinter1 + xinter2)/2.;
  ca[1]= (yinter1 + yinter2)/2.;
  ca[2]= (zinter1 + zinter2)/2.;

  // double R = std::sqrt((xinter2-xinter1)*(xinter2-xinter1) + (yinter2-yinter1)*(yinter2-yinter1) + (zinter2-zinter1)*(zinter2-zinter1));
}



void fitxyz(float X[], float Y[], float Z[], double gradXY[2], double gradZY[2])
{

  double sumX2 = 0.0, sumX = 0.0, sumXY = 0.0, sumY = 0.0;

  double sumZ2 = 0.0, sumZ = 0.0, sumZY = 0.0;

  double deltaX = 0.0, deltaZ = 0.0;

  int nPoints=3;

  for (int i = 0; i < nPoints; i++)
    { 
      sumX2 += double( X[i]*X[i] );       
      sumXY += double( X[i]*Y[i] );
      sumX  += double( X[i] );       
      sumY  += double( Y[i] );       
    }

  deltaX    = (nPoints*sumX2) - (sumX*sumX);
  gradXY[1] = ( (sumX2*sumY) - (sumX*sumXY) )/deltaX;
  gradXY[0] = ( (nPoints*sumXY) - (sumX*sumY) )/deltaX;

  sumY=0.0;

  for (int i = 0; i < nPoints; i++)
    { 
      sumZ2 += double( Z[i]*Z[i] );       
      sumZY += double( Z[i]*Y[i] );
      sumZ  += double( Z[i] );       
      sumY  += double( Y[i] );       
    }



  deltaZ    = (nPoints*sumZ2) - (sumZ*sumZ);
  gradZY[1] = ( (sumZ2*sumY) - (sumZ*sumZY) )/ deltaZ;
  gradZY[0] = ( (nPoints*sumZY) - (sumZ*sumY) )/deltaZ;
}



#define MAX_SCINT_HITS PLANES_PER_SIDE*2*10 //Just guess for now


int main(int argc, char**argv)
{
  float histoRange = SIDELENGTH*1000;

  int noBins = 100; //default value
  char inputFile[256];

  ifstream in0;
  in0.open(SOURCE_F);
  in0>>noBins>>inputFile;

  std::string outName = "pca_";
  outName.append (inputFile);

  TFile* newfile = new TFile(outName.c_str(),"RECREATE");
  
  double xPos = 0, yPos = 0, zPos = 0;

  TTree* PCA = new TTree ("PCA","PCA");
  
  PCA->Branch ("xPos", &xPos,"xPos/D");
  PCA->Branch ("yPos", &yPos,"yPos/D");
  PCA->Branch ("zPos", &zPos,"zPos/D");
 
  TH3F* PCAh = new TH3F("PCAh","PCAh",noBins,-histoRange/2,histoRange/2,  noBins,-6,histoRange, noBins,-histoRange/2,histoRange/2); 
  
  double AGradx = 0, AGradz = 0, ACutx = 0, ACutz = 0; //cuts and gradients for y as a function of x or z  (ie y = a*x + b)

  TTree* Absorbed = new TTree("Absorbed","Absorbed");

  Absorbed->Branch("xGrad", &AGradx, "gradX/D" );
  Absorbed->Branch("xCut",  &ACutx,  "cutX/D" );
  Absorbed->Branch("zGrad", &AGradz, "gradZ/D" );
  Absorbed->Branch("zCut",  &ACutz,  "cutZ/D" );

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

  float x1[PLANES_PER_SIDE],y1[PLANES_PER_SIDE],z1[PLANES_PER_SIDE];

  double xyGrad[5] = {0.0,0.0,0.0,0.0,0.0}; // gradient for y = ax + b
  double zyGrad[5] = {0.0,0.0,0.0,0.0,0.0}; // as above in yz plane
  double xyCut[5]  = {0.0,0.0,0.0,0.0,0.0};  // axis intercept xy plane (b)
  double zyCut[5]  = {0.0,0.0,0.0,0.0,0.0};  // as above yz plane

  int nLines = 0;
  double gradx[2];
  double gradz[2];
  
  int nEntries = scintTree->GetEntries();

 
  for(int entry=0;  entry < nEntries; entry++)
    {
      if( !(entry%1000) ) std::cout<< entry << " of " << nEntries <<std::endl;
      scintTree->GetEntry(entry);

      std::cout << numScintHits << "\n";

    //   for (int i = 0; i<6; i++)
// 	{  
// 	  int hits=0; // assume every third hit is a new set of events ie either entrance or exit

// 	  for (int j=0; j<3; j++)
// 	    {   
// 	      if (gdata->boxy[j].nhits[i]>=1) //if there has been an impact on this layer + face
// 		{
// 		  int k = 0; 
// 		  x1[hits] = gdata->boxy[j].scintx[i][k];
// 		  y1[hits] = gdata->boxy[j].scinty[i][k];
// 		  z1[hits] = gdata->boxy[j].scintz[i][k];

// 		  energy->Fill(gdata->boxy[j].xdedx[i][k]);

// 		  previousEvtNo = currentEvtNo;
// 		  currentEvtNo = gdata->evtno;
		 
// 		  hits++;

// 		  impactNumber [i]++;
		  
// 		}
// 	    }

// 	  if (hits==3)
// 	    { 
// 	      fitxyz(x1,y1,z1, gradx,gradz); 

// 	      xyGrad[nLines] = gradx[0];
// 	      zyGrad[nLines] = gradz[0];

// 	      xyCut[nLines]  = gradx[1];
// 	      zyCut[nLines]  = gradz[1];

// 	      nLines++;
// 	    }
// 	}

//       double pca[3]; // point of closest approach x = 0 etc;

//       if(nLines==2 && previousEvtNo == currentEvtNo)
// 	{ 
// 	  CloseApproach(xyGrad,zyGrad,xyCut,zyCut,pca);

// 	  xPos = pca[0];
// 	  yPos = pca[1];
// 	  zPos = pca[2];

// 	  PCA->Fill ();

// 	  PCAh->Fill(pca[0],pca[1],pca[2]);
// 	}
//       else if (nLines == 1 && previousEvtNo == currentEvtNo) // assume all muons enter through top
// 	{
// 	  AbsorbedHisto->Fill(x1[2],z1[2]);

// 	  AGradx = gradx[0];
// 	  ACutx  = gradx[1];

// 	  AGradz = gradz[0];
// 	  ACutz  = gradz[1];

// 	  Absorbed ->Fill();
// 	}
    }

  newfile->cd();

  newfile->Write();

}


