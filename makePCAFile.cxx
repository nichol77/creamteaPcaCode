/*
// This reads in the output file from simpleGeantSim and calculates the point
// of closest approach for those muons that hit the top and bottom planes
// of the scintillator
//
// This is based on StepThrough
//
// Ryan Nichol <rjn@hep.ucl.ac.uk>
*/


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
#include <TVector3.h>

// std headers
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>


//Some global variables read from the world tree
Double_t scintLength;
Double_t planeWidth;
Double_t gapBetweenPlanes;
Int_t planesPerSide;
Int_t stripsPerPlane;

//The steering file
#define SOURCE_F "sourcePca.dat"
#define SMALL_NUM  0.00000001

namespace PlaneView {
  typedef enum EPlaneView{
    kXView=0,
    kYView=1
  } PlaneView_t;
};


PlaneView::PlaneView_t getStripCoords(Int_t planeNum, Int_t stripNum, Double_t xyz[3])
{
  //  std::cout << scintLength << "\t" << planeWidth << "\t" << gapBetweenPlanes;
  PlaneView::PlaneView_t view=PlaneView::kXView;
  static Double_t stripWidth=scintLength/stripsPerPlane;
  Int_t onTop=1;
  if(planeNum>=planesPerSide) {
    onTop=0;
    planeNum-=planesPerSide;
  }

  Double_t zPos=(scintLength/2)+gapBetweenPlanes*planeNum + planeWidth*0.5;
  Double_t xPos=0;
  Double_t yPos=0;
  if(!onTop) {
    zPos=(-1*zPos)+planeWidth;
  }
  if(planeNum%2==0) {
    //These planes are segmented along the x direction
    xPos=(-1*scintLength/2) + (stripNum+0.5)*stripWidth;   
  }
  else {
    view=PlaneView::kYView;
    yPos=(-1*scintLength/2) + (stripNum+0.5)*stripWidth;   
  }
  xyz[0]=xPos*1000; //mm
  xyz[1]=yPos*1000; //mm
  xyz[2]=zPos*1000; //mm
  return view;
}

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


void PointOfCloseApproach(double gradXz[2], double gradYz[2], double intXz[2], double intYz[2], double pca[3])
{ 
//   //Here we switch from z=mx+c to x=az+b
//   Double_t a1=1./gradXz[0];
//   Double_t b1=-1*intXz[0]/gradXz[0];
//   Double_t a2=1./gradXz[1];
//   Double_t b2=-1*intXz[1]/gradXz[1];
//   //Here we switch from z=my+c to y=cz+d
//   Double_t c1=1./gradYz[0];
//   Double_t d1=-1*intYz[0]/gradYz[0];
//   Double_t c2=1./gradYz[1];
//   Double_t d2=-1*intYz[1]/gradYz[1];
  Double_t a1=gradXz[0];
  Double_t b1=intXz[0];
  Double_t a2=gradXz[1];
  Double_t b2=intXz[1];
  Double_t c1=gradYz[0];
  Double_t d1=intYz[0];
  Double_t c2=gradYz[1];
  Double_t d2=intYz[1];
  
  //  std::cout << a1 << "\t" << b1 << "\t" << a2 << "\t" << b2 << "\n";
  //  std::cout << c1 << "\t" << d1 << "\t" << c2 << "\t" << d2 << "\n";


  Double_t bestZ=-1*((a1-a2)*(b1-b2)+(c1-c2)*(d1-d2))/((a1-a2)*(a1-a2)+(c1-c2)*(c1-c2));
  Double_t bestX=0.5*((a1*bestZ+b1)+(a2*bestZ+b2));
  Double_t bestY=0.5*((c1*bestZ+d1)+(c2*bestZ+d2));
  
  //  std::cout << "Vals: " << bestX << "\t" << bestY << "\t" << bestZ << "\n";

  pca[0]=bestX;
  pca[1]=bestY;
  pca[2]=bestZ;
  
}


void fitxyz(int nPoints, double X[], double Y[], double Z[], double gradXZ[2], double gradYZ[2])
{

  double sumZ2 = 0.0, sumX = 0.0, sumXZ = 0.0, sumZ = 0.0;
  //  double sumY2 = 0.0, 
  double sumY = 0.0, sumYZ = 0.0;
  //  double sumX2 = 0.0;

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

void findGradients(int numPoints, double X[], double Z[], double gradXZ[])
{
 
  double sumZ2 = 0.0, sumX = 0.0, sumXZ = 0.0, sumZ = 0.0;
  //  double sumX2 = 0.0;


  for (int i = 0; i < numPoints; i++)
    { 
      sumZ2 += double( Z[i]*Z[i] );       
      sumXZ += double( X[i]*Z[i] );
      sumX  += double( X[i] );       
      sumZ  += double( Z[i] );       
    }

  sumZ2/=numPoints;
  sumXZ/=numPoints;
  sumX/=numPoints;
  sumZ/=numPoints;
  
  double covxz=sumXZ-sumX*sumZ;
  double sigmaz2=sumZ2-sumZ*sumZ;
  //  std::cout << sumX << "\t" << sumZ << "\t" << covxz << "\t" << sigmaz2 << "\n";

  gradXZ[0] = covxz/sigmaz2;
  gradXZ[1] = sumX - gradXZ[0]*sumZ;
}

#define MAX_SCINT_HITS 200 //Just guess for now
#define MAX_PLANES_PER_SIDE 100

int main(int argc, char**argv)
{


  int noBins = 100; //default value
  char inputFile[256];
  char outputFile[256];

  ifstream in0;
  in0.open(SOURCE_F);
  in0>>noBins>>inputFile;
  in0.close();
  std::cout<<inputFile<<std::endl;

  TFile* fp = new TFile(inputFile,"READ");
  TTree *worldTree = (TTree*) fp->Get("worldTree");
  worldTree->SetMakeClass(1);
  worldTree->SetBranchAddress("scintLength",&scintLength);
  worldTree->SetBranchAddress("planeWidth",&planeWidth);
  worldTree->SetBranchAddress("gapBetweenPlanes",&gapBetweenPlanes);
  worldTree->SetBranchAddress("stripsPerPlane",&stripsPerPlane);
  worldTree->SetBranchAddress("planesPerSide",&planesPerSide);
  worldTree->GetEntry(0);

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

  //Note the below code is lazy and causes a memory leak
  sprintf(outputFile,"%s/pca_%s",gSystem->DirName(inputFile),gSystem->BaseName(inputFile));
  double histoRange = scintLength*1000;
  TFile* newfile = new TFile(outputFile,"RECREATE");

  double xPosTrue = 0, yPosTrue = 0, zPosTrue = 0;
  double xPosReco = 0, yPosReco = 0, zPosReco = 0;
  Double_t thetaTrue=0;
  Double_t thetaReco=0;
  TTree* pcaTree = new TTree ("pcaTree","pcaTree");
  pcaTree->Branch ("xPosTrue", &xPosTrue,"xPosTrue/D");
  pcaTree->Branch ("yPosTrue", &yPosTrue,"yPosTrue/D");
  pcaTree->Branch ("zPosTrue", &zPosTrue,"zPosTrue/D");
  pcaTree->Branch("thetaTrue",&thetaTrue,"thetaTrue/D");

  pcaTree->Branch ("xPosReco", &xPosReco,"xPosReco/D");
  pcaTree->Branch ("yPosReco", &yPosReco,"yPosReco/D");
  pcaTree->Branch ("zPosReco", &zPosReco,"zPosReco/D");
  pcaTree->Branch("thetaReco",&thetaReco,"thetaReco/D");

  TH3F* PCAh = new TH3F("PCAh","PCAh",noBins,-histoRange/2,histoRange/2,  noBins,-histoRange/2,histoRange/2, noBins,-histoRange/2,histoRange/2); 
  
  double AGradx = 0, AGrady = 0, ACutx = 0, ACuty = 0; //cuts and gradients for y as a function of x or z  (ie y = a*x + b)

  TTree* Absorbed = new TTree("Absorbed","Absorbed");
  Absorbed->Branch("xGrad", &AGradx, "gradX/D" );
  Absorbed->Branch("xCut",  &ACutx,  "cutX/D" );
  Absorbed->Branch("yGrad", &AGrady, "gradY/D" );
  Absorbed->Branch("yCut",  &ACuty,  "cutY/D" );

  TH2F* AbsorbedHisto = new TH2F("AbsorbedH","AbsorbedH", noBins, -histoRange/2, histoRange/2, noBins, -histoRange/2, histoRange/2);


  double xTop[MAX_PLANES_PER_SIDE],yTop[MAX_PLANES_PER_SIDE],zTop[MAX_PLANES_PER_SIDE];
  double xBot[MAX_PLANES_PER_SIDE],yBot[MAX_PLANES_PER_SIDE],zBot[MAX_PLANES_PER_SIDE];

  double xTopReco[MAX_PLANES_PER_SIDE],yTopReco[MAX_PLANES_PER_SIDE],xzTopReco[MAX_PLANES_PER_SIDE],yzTopReco[MAX_PLANES_PER_SIDE];
  double xBotReco[MAX_PLANES_PER_SIDE],yBotReco[MAX_PLANES_PER_SIDE],xzBotReco[MAX_PLANES_PER_SIDE],yzBotReco[MAX_PLANES_PER_SIDE];

  //The [2] are for the top and bottom
  double xzGradTrue[2] = {0.0}; // gradient for z = ax + b
  double yzGradTrue[2] = {0.0}; // as above in yz plane
  double xzCutTrue[2]  = {0.0};  // axis intercept xz plane (b)
  double yzCutTrue[2]  = {0.0};  // as above yz plane
  double xzGradReco[2] = {0.0}; // gradient for z = ax + b
  double yzGradReco[2] = {0.0}; // as above in yz plane
  double xzCutReco[2]  = {0.0};  // axis intercept xz plane (b)
  double yzCutReco[2]  = {0.0};  // as above yz plane


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
	int gotPlane[MAX_PLANES_PER_SIDE*2]={0};
	Double_t xTrue[MAX_PLANES_PER_SIDE*2]={0},yTrue[MAX_PLANES_PER_SIDE*2]={0},zTrue[MAX_PLANES_PER_SIDE*2]={0};
	Double_t xReco[MAX_PLANES_PER_SIDE*2]={0},yReco[MAX_PLANES_PER_SIDE*2]={0},zReco[MAX_PLANES_PER_SIDE*2]={0};
	Int_t whichView[MAX_PLANES_PER_SIDE*2]={0};
	Double_t weights[MAX_PLANES_PER_SIDE*2]={0};
	for(int hit=0;hit<numScintHits;hit++) {
	  //	  std::cout << plane[hit] << "\t" << strip[hit] << "\t" << truePos[hit][2] << "\t" << energyDep[hit] << "\n";
	  Double_t recoXYZ[3];
	  whichView[plane[hit]]=getStripCoords(plane[hit],strip[hit],recoXYZ);
	  //	  std::cout << plane[hit] << "\t" << strip[hit] << "\n";
	  //	  std::cout << "Reco: " << recoXYZ[0] << "\t" << recoXYZ[1] << "\t"
	  //		    << recoXYZ[2] << "\n";
	  //	  std::cout << "True: " << truePos[hit][0] << "\t" << truePos[hit][1] << "\t"
	  //		    << truePos[hit][2] << "\n";

	  gotPlane[plane[hit]]=1;
	  xTrue[plane[hit]]+=truePos[hit][0]*energyDep[hit];
	  yTrue[plane[hit]]+=truePos[hit][1]*energyDep[hit];
	  zTrue[plane[hit]]+=truePos[hit][2]*energyDep[hit];
	  xReco[plane[hit]]+=recoXYZ[0]*energyDep[hit];
	  yReco[plane[hit]]+=recoXYZ[1]*energyDep[hit];
	  zReco[plane[hit]]+=recoXYZ[2]*energyDep[hit];
	  weights[plane[hit]]+=energyDep[hit];
	}	
	
	int topPlanes=0;
	int botPlanes=0;
	int topPlanesXView=0;
	int topPlanesYView=0;
	int botPlanesXView=0;
	int botPlanesYView=0;
	

	//Normalise to find the actual values
	for(int plane=0;plane<planesPerSide*2;plane++) {	  
	  if(gotPlane[plane]) {
	    xTrue[plane]/=weights[plane];
	    yTrue[plane]/=weights[plane];
	    zTrue[plane]/=weights[plane];
	    xReco[plane]/=weights[plane];
	    yReco[plane]/=weights[plane];
	    zReco[plane]/=weights[plane];
	  }
	}
	

	//Now determine the point of closest approach using the true values
	for(int plane=0;plane<planesPerSide;plane++) {	  
	    //	    std::cout << plane << "\t" << xTrue[plane] << "\t" << yTrue[plane] << "\t" << zTrue[plane] << "\n";
	  if(gotPlane[plane]) {
	    xTop[topPlanes]=xTrue[plane];
	    yTop[topPlanes]=yTrue[plane];
	    zTop[topPlanes]=zTrue[plane];
	    topPlanes++;

	    if(whichView[plane]==PlaneView::kXView) {
	      xTopReco[topPlanesXView]=xReco[plane];
	      xzTopReco[topPlanesXView]=zReco[plane];
	      topPlanesXView++;
	    }
	    else {
	      yTopReco[topPlanesYView]=yReco[plane];
	      yzTopReco[topPlanesYView]=zReco[plane];
	      topPlanesYView++;
	    }
	  }
	  if(gotPlane[plane+ planesPerSide]) {
	    xBot[botPlanes]=xTrue[plane+ planesPerSide];
	    yBot[botPlanes]=yTrue[plane+ planesPerSide];
	    zBot[botPlanes]=zTrue[plane+ planesPerSide];
	    botPlanes++;

	    if(whichView[plane]==PlaneView::kXView) {
	      xBotReco[botPlanesXView]=xReco[plane+planesPerSide];
	      xzBotReco[botPlanesXView]=zReco[plane+planesPerSide];
	      botPlanesXView++;
	    }
	    else {
	      yBotReco[botPlanesYView]=yReco[plane+planesPerSide];
	      yzBotReco[botPlanesYView]=zReco[plane+planesPerSide];
	      botPlanesYView++;
	    }

	  }
	}

	//	std::cout << "Planes:\t" << topPlanes << "\t" << botPlanes << "\n";
	if(topPlanes<4) continue; //Didn't hit the top

	//Now find the gradients using the true positions
	Double_t gradxz[2];
	Double_t gradyz[2];
	fitxyz(topPlanes,xTop,yTop,zTop,gradxz,gradyz);
	xzGradTrue[0]=gradxz[0];
	yzGradTrue[0]=gradyz[0];
	xzCutTrue[0]=gradxz[1];
	yzCutTrue[0]=gradyz[1];
       
	//Now find the gradients using the reco positions
	Double_t recoGradXZ[2];
	Double_t recoGradYZ[2];
	findGradients(topPlanesXView,xTopReco,xzTopReco,recoGradXZ);
	findGradients(topPlanesYView,yTopReco,yzTopReco,recoGradYZ);
	xzGradReco[0]=recoGradXZ[0];
	yzGradReco[0]=recoGradYZ[0];
	xzCutReco[0]=recoGradXZ[1];
	yzCutReco[0]=recoGradYZ[1];

	//	std::cout << topPlanesXView << "\t" << recoGradXZ[0] << "\t" << recoGradXZ[1] << "\n";
	//	std::cout << topPlanes << "\t" << gradxz[0] << "\t" << gradxz[1] << "\n";


	if(botPlanes>3) {
	  //Potential scatter
	  fitxyz(botPlanes,xBot,yBot,zBot,gradxz,gradyz);
	  xzGradTrue[1]=gradxz[0];
	  yzGradTrue[1]=gradyz[0];
	  xzCutTrue[1]=gradxz[1];
	  yzCutTrue[1]=gradyz[1];


	  findGradients(botPlanesXView,xBotReco,xzBotReco,recoGradXZ);
	  findGradients(botPlanesYView,yBotReco,yzBotReco,recoGradYZ);
	  xzGradReco[1]=recoGradXZ[0];
	  yzGradReco[1]=recoGradYZ[0];
	  xzCutReco[1]=recoGradXZ[1];
	  yzCutReco[1]=recoGradYZ[1];

	  //	  std::cout << "RecoXZ: " << xzGradReco[0] << "\t" << xzGradReco[1] << "\n";

	  Double_t pcaTrue[3];
	  PointOfCloseApproach(xzGradTrue,yzGradTrue,xzCutTrue,yzCutTrue,pcaTrue);
	  xPosTrue = pcaTrue[0];
	  yPosTrue = pcaTrue[1];
	  zPosTrue = pcaTrue[2];

	  TVector3 intialDirTrue(xzGradTrue[0],yzGradTrue[0],1);
	  TVector3 finalDirTrue(xzGradTrue[1],yzGradTrue[1],1);
	  thetaTrue=intialDirTrue.Angle(finalDirTrue);

	  //Now do reco PCA
	  Double_t pcaReco[3];
	  PointOfCloseApproach(xzGradReco,yzGradReco,xzCutReco,yzCutReco,pcaReco);
	  xPosReco = pcaReco[0];
	  yPosReco = pcaReco[1];
	  zPosReco = pcaReco[2];
	  //	  std::cout << xPosReco << "\t" << yPosReco << "\t" << zPosReco << "\n";

	  TVector3 intialDirReco(xzGradReco[0],yzGradReco[0],1);
	  TVector3 finalDirReco(xzGradReco[1],yzGradReco[1],1);
	  thetaReco=intialDirReco.Angle(finalDirReco);


	 
	  pcaTree->Fill();

	  PCAh->Fill(pcaTrue[0],pcaTrue[1],pcaTrue[2]);
	  
	
	}
	else {
	  //Potential Absorbtion ... need to tweak them
 	  AbsorbedHisto->Fill(xTop[0],yTop[0]);

	  AGradx=xzGradTrue[0];
	  AGrady=yzGradTrue[0];
	  ACutx=xzCutTrue[0];
	  ACuty=yzCutTrue[0];
 	  Absorbed ->Fill();
	}

      }

    }

  newfile->cd();

  newfile->Write();

}


