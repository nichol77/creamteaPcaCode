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
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TStyle.h>
#include <TMath.h>
#include <TSystem.h>
#include <TVector3.h>
#include <TMarker.h>

// std headers
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>



//Some global variables read from the world tree
Int_t minervaStrips;
Double_t scintTriBase;
Double_t scintTriHeight;
Double_t scintTriLengthX;
Double_t scintTriLengthY;
Int_t numScintTriX;
Int_t numScintTriY;
Double_t scintLength;
Double_t planeWidth;
Double_t gapBetweenPlanes;
Double_t verticalSeparation;
Int_t planesPerSide=8;
Int_t stripsPerPlane=650;

//The steering file
#define SMALL_NUM  0.00000001
#define MIN_ENERGY_THRESHOLD 0
#define MAX_SCINT_HITS 400 //Just guess for now
#define MAX_PLANES_PER_SIDE 100

namespace PlaneView {
  typedef enum EPlaneView{
    kXView=0,
    kYView=1
  } PlaneView_t;
};


PlaneView::PlaneView_t getStripCoords(Int_t planeNum, Int_t stripNum, Double_t xyz[3])
{
  //  std::cout << scintLength << "\t" << planeWidth << "\t" << gapBetweenPlanes;
  static Int_t firstTime=1;
  static Double_t stripZeroX=0;  
  static Double_t stripZeroY=0;  
  if(firstTime) {
    if(minervaStrips) {
      stripZeroX=-1*((numScintTriX/2)*(scintTriBase/2.));
      if(numScintTriX%2==0) {
	stripZeroX+=scintTriBase/4.;
      }
      stripZeroY=-1*((numScintTriY/2)*(scintTriBase/2.));
      if(numScintTriY%2==0) {
	stripZeroY+=scintTriBase/4.;
      }
    }
    firstTime=0;
  }

  if(!minervaStrips) {
    PlaneView::PlaneView_t view=PlaneView::kXView;
    static Double_t stripWidth=scintLength/stripsPerPlane;

    Int_t onTop=1;
    if(planeNum>=planesPerSide) {
      onTop=0;
      planeNum-=planesPerSide;
    }
    
    Double_t zPos=(verticalSeparation/2)+gapBetweenPlanes*planeNum ;//+ planeWidth*0.5;
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
  else {
    PlaneView::PlaneView_t view=PlaneView::kXView;
    Int_t onTop=1;
    if(planeNum>=planesPerSide) {
      onTop=0;
      planeNum-=planesPerSide;
    }
    Double_t zPos=(verticalSeparation/2)+gapBetweenPlanes*planeNum;// + scintTriHeight*0.5;
    Double_t xPos=0;
    Double_t yPos=0;
    if(!onTop) {
      zPos=(-1*zPos)+planeWidth;
    }
    if(planeNum%2==0) {
      //These planes are segmented along the x direction
      xPos=stripZeroX+(stripNum*scintTriBase/2.);
    }
    else {
      view=PlaneView::kYView;
      yPos=stripZeroY+(stripNum*scintTriBase/2.);
    }
    xyz[0]=xPos*1000; //mm
    xyz[1]=yPos*1000; //mm
    xyz[2]=zPos*1000; //mm
    return view;



  }
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
 
  //If S=dx^2+dy^2 = az^2 + bz +c
  //  Double_t a=(a1-a2)*(a1-a2) + (c1-c2)*(c1-c2);



  Double_t bestZ=-1*(((a1-a2)*(b1-b2))+((c1-c2)*(d1-d2)))/(((a1-a2)*(a1-a2))+((c1-c2)*(c1-c2)));
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


void fitxyzw(int nPoints, double X[], double Y[], double Z[], double W[], double gradXZ[2], double gradYZ[2])
{

 //  Double_t sumW=0;

//   double sumZ2 = 0.0, sumX = 0.0, sumXZ = 0.0, sumZ = 0.0;
//   double sumY = 0.0, sumYZ = 0.0;
//   double sumW2 =0.0;


//   for (int i = 0; i < nPoints; i++)
//     { 
//       sumW  += double ( W[i] );
//       sumW2 += double( W[i]*W[i] );
//       sumZ2 += double( W[i]*Z[i]*Z[i] );       
//       sumXZ += double( W[i]*X[i]*Z[i] );
//       sumYZ += double( W[i]*Y[i]*Z[i] );
//       sumX  += double( W[i]*X[i] );       
//       sumY  += double( W[i]*Y[i] );       
//       sumZ  += double( W[i]*Z[i] );       
//     }
//   sumW2 /= sumW;
//   sumZ2/=sumW;
//   sumXZ/=sumW;
//   sumYZ/=sumW;
//   sumX/=sumW;
//   sumY/=sumW;
//   sumZ/=sumW;
  
//   double covxz=(sumXZ-sumX*sumZ)/(1-sumW2);
//   double covyz=(sumYZ-sumY*sumZ)/(1-sumW2);
//   double sigmaz2=(sumZ2-sumZ*sumZ)/(1-sumW2);

//   gradXZ[0] = covxz/sigmaz2;
//   gradXZ[1] = sumX - gradXZ[0]*sumZ;

//   gradYZ[0] = covyz/sigmaz2;
//   gradYZ[1] = sumY - gradYZ[0]*sumZ;
  TF1 fpol1("fpol1","pol1",-8000,8000);
  TGraphErrors fredX(nPoints,Z,X,0,W);
  fredX.Fit("fpol1","Q");
  gradXZ[0]=fpol1.GetParameter(1);
  gradXZ[1]=fpol1.GetParameter(0);
  TGraphErrors fredY(nPoints,Z,Y,0,W);
  fredY.Fit("fpol1","Q");
  gradYZ[0]=fpol1.GetParameter(1);
  gradYZ[1]=fpol1.GetParameter(0);

}


void fitxzw(int nPoints, double X[], double Z[], double W[], double gradXZ[2])
{
  TF1 fpol1("fpol1","pol1",-8000,8000);
  TGraphErrors fredX(nPoints,Z,X,0,W);
  fredX.Fit("fpol1","Q");
  gradXZ[0]=fpol1.GetParameter(1);
  gradXZ[1]=fpol1.GetParameter(0);
  
}

Double_t iterativelyFitXYZWithWeights(int nPoints, double X[], double Y[], double Z[], double W[], double gradXZ[2], double gradYZ[2])
{
  fitxyzw(nPoints,X,Y,Z,W,gradXZ,gradYZ);
  Double_t newX[MAX_SCINT_HITS];
  Double_t newY[MAX_SCINT_HITS];
  Double_t newZ[MAX_SCINT_HITS];
  Double_t newW[MAX_SCINT_HITS];
  Int_t newPoints=0;
  Int_t excludePoint[MAX_SCINT_HITS]={0};
  
  Double_t newGradXZ[2]={gradXZ[0],gradXZ[1]};
  Double_t newGradYZ[2]={gradYZ[0],gradYZ[1]};

  Double_t meanDiff=0;
  for(int it=0;it<nPoints-1;it++) {
    Double_t maxDiff=0;
    Int_t maxPoint=0;   
    meanDiff=0;
    Int_t countDiffs=0;
    for(int i=0;i<nPoints;i++) {    
      if(excludePoint[i]) continue;
      Double_t diff=((X[i]-(newGradXZ[0]*Z[i]+newGradXZ[1]))*
		     (X[i]-(newGradXZ[0]*Z[i]+newGradXZ[1])));
      diff+=((Y[i]-(newGradYZ[0]*Z[i]+newGradYZ[1]))*
	     (Y[i]-(newGradYZ[0]*Z[i]+newGradYZ[1])));
      meanDiff+=diff;
      countDiffs++;
      //      std::cout << "Point " << i << "\t" << diff << "\n";
      if(diff>maxDiff) {
	maxDiff=diff;
	maxPoint=i;
      }
    }
    meanDiff/=countDiffs;
    //    std::cout << "Max Point " << maxPoint << "\t" << maxDiff << "\t" << meanDiff <<  "\n";
    if(maxDiff<1) break;

    excludePoint[maxPoint]=1;
    newPoints=0;
    for(int i=0;i<nPoints;i++) {    
      if(excludePoint[i]) continue;
      newX[newPoints]=X[i];
      newY[newPoints]=Y[i];
      newZ[newPoints]=Z[i];
      newW[newPoints]=W[i];
      newPoints++;
    }
    fitxyzw(newPoints,newX,newY,newZ,newW,newGradXZ,newGradYZ);
    

  //   std::cout << TMath::Abs(gradXZ[0]-newGradXZ[0]) << "\t"
// 	      << TMath::Abs(gradXZ[1]-newGradXZ[1]) << "\t"
// 	      << TMath::Abs(gradYZ[0]-newGradYZ[0]) << "\t"
// 	      << TMath::Abs(gradYZ[1]-newGradYZ[1]) << "\n";

    gradXZ[0]=newGradXZ[0];
    gradXZ[1]=newGradXZ[1];
    gradYZ[0]=newGradYZ[0];
    gradYZ[1]=newGradYZ[1];

  }

  if(newPoints>0 && newPoints<3 && nPoints>6) {
    std::cout << "Redoing after using: " << newPoints << " of " << nPoints << "\n";
    //Try another tack and just use all the planes
   fitxyzw(nPoints,X,Y,Z,W,gradXZ,gradYZ);
   meanDiff=0;
   for(int i=0;i<nPoints;i++) {    
     Double_t diff=((X[i]-(gradXZ[0]*Z[i]+gradXZ[1]))*
		    (X[i]-(gradXZ[0]*Z[i]+gradXZ[1])));
     diff+=((Y[i]-(gradYZ[0]*Z[i]+gradYZ[1]))*
	    (Y[i]-(gradYZ[0]*Z[i]+gradYZ[1])));
     meanDiff+=diff;
   }
   meanDiff/=nPoints;
  }
  return meanDiff;
}


Double_t iterativelyFitXZWithWeights(int nPoints, double X[], double Z[], double W[], double gradXZ[2])
{
  fitxzw(nPoints,X,Z,W,gradXZ);
//   Double_t newX[MAX_SCINT_HITS];
//   Double_t newZ[MAX_SCINT_HITS];
//   Double_t newW[MAX_SCINT_HITS];
//   Int_t newPoints=0;
//   Int_t excludePoint[MAX_SCINT_HITS]={0};
  
  Double_t newGradXZ[2]={gradXZ[0],gradXZ[1]};

  Double_t meanDiff=0;
//   for(int it=0;it<nPoints-1;it++) {
//     Double_t maxDiff=0;
//     Int_t maxPoint=0;   
//     meanDiff=0;
//     Int_t countDiffs=0;
//     for(int i=0;i<nPoints;i++) {    
//       if(excludePoint[i]) continue;
//       Double_t diff=((X[i]-(newGradXZ[0]*Z[i]+newGradXZ[1]))*
// 		     (X[i]-(newGradXZ[0]*Z[i]+newGradXZ[1])));
//       meanDiff+=diff;
//       countDiffs++;
//       std::cout << "Point " << i << "\t" << X[i] << "\t" << Z[i] 
// 		<< "\t" << diff << "\n";
//       if(diff>maxDiff) {
// 	maxDiff=diff;
// 	maxPoint=i;
//       }
//     }
//     meanDiff/=countDiffs;
//     std::cout << "Max Point " << maxPoint << "\t" << maxDiff << "\t" << meanDiff <<  "\n";
//     if(maxDiff<1) break;

//     excludePoint[maxPoint]=1;
//     newPoints=0;
//     for(int i=0;i<nPoints;i++) {    
//       if(excludePoint[i]) continue;
//       newX[newPoints]=X[i];
//       newZ[newPoints]=Z[i];
//       newW[newPoints]=W[i];
//       newPoints++;
//     }
//     fitxzw(newPoints,newX,newZ,newW,newGradXZ);
    

//   //   std::cout << TMath::Abs(gradXZ[0]-newGradXZ[0]) << "\t"
// // 	      << TMath::Abs(gradXZ[1]-newGradXZ[1]) << "\n";


//     gradXZ[0]=newGradXZ[0];
//     gradXZ[1]=newGradXZ[1];
//   }

//   if(newPoints>0 && newPoints<3 && nPoints>2) {
//     std::cout << "Redoing after using: " << newPoints << " of " << nPoints << "\n";
//     //Try another tack and just use all the planes
//    fitxzw(nPoints,X,Z,W,gradXZ);
//    meanDiff=0;
//    for(int i=0;i<nPoints;i++) {    
//      Double_t diff=((X[i]-(newGradXZ[0]*Z[i]+newGradXZ[1]))*
// 		    (X[i]-(newGradXZ[0]*Z[i]+newGradXZ[1])));
//      meanDiff+=diff;
//    }
//    meanDiff/=nPoints;
//   }
  

  meanDiff=0;
  for(int i=0;i<nPoints;i++) {    
    Double_t diff=((X[i]-(gradXZ[0]*Z[i]+gradXZ[1]))*
		   (X[i]-(gradXZ[0]*Z[i]+gradXZ[1])));
    meanDiff+=diff;
  }
  meanDiff/=nPoints;
  return meanDiff;
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


void iterativelyFindGradients(Int_t numPoints, Double_t xVals[], Double_t zVals[], Double_t gradXZ[]) 
{
  Double_t newX[MAX_PLANES_PER_SIDE];
  Double_t newZ[MAX_PLANES_PER_SIDE];
  
  
  for(int i=0;i< numPoints;i++) {
    newX[i]=xVals[i];
    newZ[i]=zVals[i];
  }
  findGradients(numPoints,newX,newZ,gradXZ);
    
}


Double_t getSSquared(Double_t *x, Double_t *par) {
  Double_t m[4]={par[0],par[2],par[4],par[6]};
  Double_t c[4]={par[1],par[3],par[5],par[7]};
  
  Double_t z=x[0];
  Double_t x1=m[0]*z+c[0];
  Double_t x2=m[1]*z+c[1];
  Double_t y1=m[2]*z+c[2];
  Double_t y2=m[3]*z+c[3];
  Double_t sSq=((x1-x2)*(x1-x2))+((y1-y2)*(y1-y2));
  return sSq;
}




TFile *fp=0;
TTree *worldTree=0;
TTree *scintChain=0;
Int_t fRun,fEvent;
Int_t numScintHits;
Int_t side[MAX_SCINT_HITS];
Int_t plane[MAX_SCINT_HITS];
Int_t strip[MAX_SCINT_HITS];
Double_t truePos[MAX_SCINT_HITS][3];
Double_t energyDep[MAX_SCINT_HITS];

void loadFiles(char *inputDir, char *rootName, int fileNum)
{  
  std::cout << inputDir << "\t" << rootName << "\t" << fileNum << "\n";

  //  int noBins = 100; //default value

  char inputFile[256];
  //  Int_t countFiles=0;
  
  sprintf(inputFile,"%s/%s_%d.root",inputDir,rootName,fileNum);
  fp = new TFile(inputFile,"READ");
  //    std::cout << inputFile << "\t" << fp << "\n";
  if(!worldTree) {
    worldTree = (TTree*) fp->Get("worldTree");
    worldTree->SetMakeClass(1);
    worldTree->SetBranchAddress("scintLength",&scintLength);
    worldTree->SetBranchAddress("planeWidth",&planeWidth);
    worldTree->SetBranchAddress("gapBetweenPlanes",&gapBetweenPlanes);
    worldTree->SetBranchAddress("stripsPerPlane",&stripsPerPlane);
    worldTree->SetBranchAddress("planesPerSide",&planesPerSide);
    worldTree->SetBranchAddress("verticalSeparation",&verticalSeparation);
    worldTree->SetBranchAddress("minervaStrips",&minervaStrips);
    worldTree->SetBranchAddress("scintTriBase",&scintTriBase);
    worldTree->SetBranchAddress("scintTriHeight",&scintTriHeight);
    worldTree->SetBranchAddress("scintTriLengthX",&scintTriLengthX);
    worldTree->SetBranchAddress("scintTriLengthY",&scintTriLengthY);
    worldTree->SetBranchAddress("numScintTriX",&numScintTriX);
    worldTree->SetBranchAddress("numScintTriY",&numScintTriY);
    worldTree->GetEntry(0);
    //RJN hack for now
    //      verticalSeparation=3;
  }
  
  scintChain = (TTree*) fp->Get("scintTree");
  

  
  scintChain->SetMakeClass(1);
  scintChain->SetBranchAddress("run",&fRun);
  scintChain->SetBranchAddress("event",&fEvent);
  scintChain->SetBranchAddress("ScintHitInfo",&numScintHits); 
  scintChain->SetBranchAddress("ScintHitInfo.side",side); 
  scintChain->SetBranchAddress("ScintHitInfo.plane",plane); 
  scintChain->SetBranchAddress("ScintHitInfo.strip",strip); 
  scintChain->SetBranchAddress("ScintHitInfo.truePos[3]",truePos); 
  scintChain->SetBranchAddress("ScintHitInfo.energyDep",energyDep); 
}



void plotEvent(int entryNum) {
  //Note the below code is lazy and causes a memory leak
  //  double histoRange = scintLength*1000;

  double xPosTrue = 0, yPosTrue = 0, zPosTrue = 0;
  double xPosReco = 0, yPosReco = 0, zPosReco = 0;
  Double_t thetaTrue=0;
  Double_t thetaReco=0;
  Double_t thetaxzReco=0;
  Double_t thetayzReco=0;
  Double_t thetaxzTrue=0;
  Double_t thetayzTrue=0;

  //The [2] are for the top and bottom
  double xzGradTrue[2] = {0.0}; // gradient for z = ax + b
  double yzGradTrue[2] = {0.0}; // as above in yz plane
  double xzCutTrue[2]  = {0.0};  // axis intercept xz plane (b)
  double yzCutTrue[2]  = {0.0};  // as above yz plane
  double xzGradReco[2] = {0.0}; // gradient for z = ax + b
  double yzGradReco[2] = {0.0}; // as above in yz plane
  double xzCutReco[2]  = {0.0};  // axis intercept xz plane (b)
  double yzCutReco[2]  = {0.0};  // as above yz plane
  int gotRecoPCA=0;
  PlaneView::PlaneView_t whichView[MAX_SCINT_HITS];

 
  //   TTree* pcaTree = new TTree ("pcaTree","pcaTree");
  //   pcaTree->Branch ("xPosTrue", &xPosTrue,"xPosTrue/D");
  //   pcaTree->Branch ("yPosTrue", &yPosTrue,"yPosTrue/D");
  //   pcaTree->Branch ("zPosTrue", &zPosTrue,"zPosTrue/D");
  //   pcaTree->Branch("thetaTrue",&thetaTrue,"thetaTrue/D");
  //   pcaTree->Branch("thetaxzTrue",&thetaxzTrue,"thetaxzTrue/D");
  //   pcaTree->Branch("thetayzTrue",&thetayzTrue,"thetayzTrue/D");
  //   pcaTree->Branch("xzGradTrue",&xzGradTrue,"xzGradTrue[2]/D");
  //   pcaTree->Branch("yzGradTrue",&yzGradTrue,"yzGradTrue[2]/D");


  //   pcaTree->Branch("gotRecoPCA", &gotRecoPCA,"gotRecoPCA/I");
  //   pcaTree->Branch ("xPosReco", &xPosReco,"xPosReco/D");
  //   pcaTree->Branch ("yPosReco", &yPosReco,"yPosReco/D");
  //   pcaTree->Branch ("zPosReco", &zPosReco,"zPosReco/D");
  //   pcaTree->Branch("thetaReco",&thetaReco,"thetaReco/D");

  //   pcaTree->Branch("thetaxzReco",&thetaxzReco,"thetaxzReco/D");
  //   pcaTree->Branch("thetayzReco",&thetayzReco,"thetayzReco/D");
  //   pcaTree->Branch("xzGradReco",&xzGradReco,"xzGradReco[2]/D");
  //   pcaTree->Branch("yzGradReco",&yzGradReco,"yzGradReco[2]/D");

  //   TH3F* PCAh = new TH3F("PCAh","PCAh",noBins,-histoRange/2,histoRange/2,  noBins,-histoRange/2,histoRange/2, noBins,-histoRange/2,histoRange/2); 
  
  double AGradx = 0, AGrady = 0, ACutx = 0, ACuty = 0; //cuts and gradients for y as a function of x or z  (ie y = a*x + b)

  //   TTree* Absorbed = new TTree("Absorbed","Absorbed");
  //   Absorbed->Branch("xGrad", &AGradx, "gradX/D" );
  //   Absorbed->Branch("xCut",  &ACutx,  "cutX/D" );
  //   Absorbed->Branch("yGrad", &AGrady, "gradY/D" );
  //   Absorbed->Branch("yCut",  &ACuty,  "cutY/D" );

  //   TH2F* AbsorbedHisto = new TH2F("AbsorbedH","AbsorbedH", noBins, -histoRange/2, histoRange/2, noBins, -histoRange/2, histoRange/2);


  scintChain->GetEntry(entryNum);


  //	std::cout << "Num hits:\t" << numScintHits << "\n";
  if(numScintHits>3) {
    //Minimum requirement
    
    Int_t gotPlane[MAX_PLANES_PER_SIDE*2]={0};

    Double_t xPosTop[MAX_SCINT_HITS];
    Double_t yPosTop[MAX_SCINT_HITS];
    Double_t xErrTop[MAX_SCINT_HITS];
    Double_t yErrTop[MAX_SCINT_HITS];
    Double_t zPosTop[MAX_SCINT_HITS];
    Int_t rawNumTop=0;
    Int_t rawNumBot=0;

    Double_t xPosBot[MAX_SCINT_HITS];
    Double_t yPosBot[MAX_SCINT_HITS];
    Double_t xErrBot[MAX_SCINT_HITS];
    Double_t yErrBot[MAX_SCINT_HITS];
    Double_t zPosBot[MAX_SCINT_HITS];

    Int_t xPlaneTop[MAX_SCINT_HITS];
    Int_t xStripTop[MAX_SCINT_HITS];
    Int_t yPlaneTop[MAX_SCINT_HITS];
    Int_t yStripTop[MAX_SCINT_HITS];
    Double_t xPosTopReco[MAX_SCINT_HITS];
    Double_t yPosTopReco[MAX_SCINT_HITS];
    Double_t xErrTopReco[MAX_SCINT_HITS];
    Double_t yErrTopReco[MAX_SCINT_HITS];
    Double_t zxPosTopReco[MAX_SCINT_HITS];
    Double_t zyPosTopReco[MAX_SCINT_HITS];
    Int_t rawNumTopXReco=0;
    Int_t rawNumTopYReco=0;
    Int_t xPlaneBot[MAX_SCINT_HITS];
    Int_t yPlaneBot[MAX_SCINT_HITS];
    Int_t xStripBot[MAX_SCINT_HITS];
    Int_t yStripBot[MAX_SCINT_HITS];
    Double_t xPosBotReco[MAX_SCINT_HITS];
    Double_t yPosBotReco[MAX_SCINT_HITS];
    Double_t xErrBotReco[MAX_SCINT_HITS];
    Double_t yErrBotReco[MAX_SCINT_HITS];
    Double_t zxPosBotReco[MAX_SCINT_HITS];
    Double_t zyPosBotReco[MAX_SCINT_HITS];
    Int_t rawNumBotXReco=0;
    Int_t rawNumBotYReco=0;

    for(int hit=0;hit<numScintHits;hit++) {
      if(energyDep[hit]<MIN_ENERGY_THRESHOLD) continue;
      //      std::cout << plane[hit] << "\t" << strip[hit] << "\t" << truePos[hit][2] << "\t" << energyDep[hit] << "\n";
      Double_t recoXYZ[3];
      whichView[hit]=getStripCoords(plane[hit],strip[hit],recoXYZ);
 //      if(plane[hit]%2==0) {
	std::cout << plane[hit] << "\t" << strip[hit] << "\n";
// 	std::cout << "Reco: " << recoXYZ[0] << "\t" << recoXYZ[1] << "\t"
// 		  << recoXYZ[2] << "\n";
// 	std::cout << "True: " << truePos[hit][0] << "\t" << truePos[hit][1] << "\t"
// 		  << truePos[hit][2] << "\n";
// 	std::cout << "xDiff:\t" << truePos[hit][0]-recoXYZ[0] << "\n";
//       }
      if(TMath::IsNaN(recoXYZ[0])) {
	std::cout << "Nan at getStripCoords for plane " << plane[hit] << " strip " << strip[hit] << "\n";
      }

      gotPlane[plane[hit]]=1;
      
      if(plane[hit]<planesPerSide) {
	xPosTop[rawNumTop]=truePos[hit][0];
	xErrTop[rawNumTop]=1./energyDep[hit];
	yPosTop[rawNumTop]=truePos[hit][1];
	yErrTop[rawNumTop]=1./energyDep[hit];
	zPosTop[rawNumTop]=truePos[hit][2];
	rawNumTop++;
      }
      else {
	xPosBot[rawNumBot]=truePos[hit][0];
	xErrBot[rawNumBot]=1./energyDep[hit];
	yPosBot[rawNumBot]=truePos[hit][1];
	yErrBot[rawNumBot]=1./energyDep[hit];
	zPosBot[rawNumBot]=truePos[hit][2];
	rawNumBot++;
      }
      
      
      //Now the reco ones
      if(plane[hit]<planesPerSide) {
	if(whichView[hit]==PlaneView::kXView) {
	  xPlaneTop[rawNumTopXReco]=plane[hit];
	  xStripTop[rawNumTopXReco]=strip[hit];
	  xPosTopReco[rawNumTopXReco]=recoXYZ[0];
	  xErrTopReco[rawNumTopXReco]=1./energyDep[hit];
	  zxPosTopReco[rawNumTopXReco]=recoXYZ[2];
	  rawNumTopXReco++;
	}
	else {
	  yPlaneTop[rawNumTopYReco]=plane[hit];
	  xStripTop[rawNumTopXReco]=strip[hit];
	  yPosTopReco[rawNumTopYReco]=recoXYZ[1];
	  yErrTopReco[rawNumTopYReco]=1./energyDep[hit];
	  zyPosTopReco[rawNumTopYReco]=recoXYZ[2];
	  rawNumTopYReco++;
	}
      }
      else {
	if(whichView[hit]==PlaneView::kXView) {
	  xPlaneBot[rawNumBotXReco]=plane[hit];	  
	  xStripBot[rawNumTopXReco]=strip[hit];
	  xPosBotReco[rawNumBotXReco]=recoXYZ[0];
	  xErrBotReco[rawNumBotXReco]=1./energyDep[hit];
	  zxPosBotReco[rawNumBotXReco]=recoXYZ[2];
	  rawNumBotXReco++;
	}
	else {
	  yPlaneBot[rawNumBotYReco]=plane[hit];
	  xStripBot[rawNumTopXReco]=strip[hit];
	  yPosBotReco[rawNumBotYReco]=recoXYZ[1];
	  yErrBotReco[rawNumBotYReco]=1./energyDep[hit];
	  zyPosBotReco[rawNumBotYReco]=recoXYZ[2];
	  rawNumBotYReco++;
	}
      }
    }

    Double_t xTopAvgX[MAX_SCINT_HITS]={0};
    Double_t xTopAvgZ[MAX_SCINT_HITS]={0};
    Double_t xTopAvgW[MAX_SCINT_HITS]={0};
    Int_t countXTop=0;
    Double_t yTopAvgY[MAX_SCINT_HITS]={0};
    Double_t yTopAvgZ[MAX_SCINT_HITS]={0};
    Double_t yTopAvgW[MAX_SCINT_HITS]={0};
    Int_t countYTop=0;
    Double_t xBotAvgX[MAX_SCINT_HITS]={0};
    Double_t xBotAvgZ[MAX_SCINT_HITS]={0};
    Double_t xBotAvgW[MAX_SCINT_HITS]={0};
    Int_t countXBot=0;
    Double_t yBotAvgY[MAX_SCINT_HITS]={0};
    Double_t yBotAvgZ[MAX_SCINT_HITS]={0};
    Double_t yBotAvgW[MAX_SCINT_HITS]={0};
    Int_t countYBot=0;
    for(int testPlane=0;testPlane<planesPerSide;testPlane+=2) {
      Int_t countHits=0;
      Int_t stripsHit[MAX_SCINT_HITS]={0};
      Int_t index[MAX_SCINT_HITS]={0};
      Int_t indexIndex[MAX_SCINT_HITS]={0};
      for(int i=0;i<rawNumTopXReco;i++) {
	Int_t pl=xPlaneTop[i];
	Int_t st=xStripTop[i];
	if(pl==testPlane) {
	  stripsHit[countHits]=st;
	  index[countHits]=i;
	  countHits++;
	}
      }
      TMath::Sort(countHits,stripsHit,indexIndex);
      Int_t lastStrip=-1;
      for(int hit=0;hit<countHits;hit++) {
	Int_t i=index[indexIndex[hit]];
	Int_t st=xStripTop[i];
	Int_t pl=xPlaneTop[i];
	Double_t x=xPosTopReco[i];
	Double_t z=zxPosTopReco[i];
	Double_t w=xErrTopReco[i];
	if(st==lastStrip-1) {
	  Double_t newW=xTopAvgW[countXTop-1]+w;
	  Double_t newX=(xTopAvgX[countXTop-1]*xTopAvgW[countXTop-1])+(x*w);
	  Double_t newZ=(xTopAvgZ[countXTop-1]*xTopAvgW[countXTop-1])+(z*w);
	  newX/=newW;
	  newZ/=newW;
	  xTopAvgW[countXTop-1]=newW;
	  xTopAvgZ[countXTop-1]=newZ;
	  xTopAvgX[countXTop-1]=newX;	  
	}
	else {
	  xTopAvgX[countXTop]=x;
	  xTopAvgZ[countXTop]=z;
	  xTopAvgW[countXTop]=w;
	  countXTop++;
	}
	lastStrip=st;
	//	std::cout << pl << "\t" << st << "\t"  <<  x 
	//		<< "\t" << z << "\t" << w << "\n";
      }
    }

    //Next up same thing for yTop
    for(int testPlane=1;testPlane<planesPerSide;testPlane+=2) {
      Int_t countHits=0;
      Int_t stripsHit[MAX_SCINT_HITS]={0};
      Int_t index[MAX_SCINT_HITS]={0};
      Int_t indexIndex[MAX_SCINT_HITS]={0};
      for(int i=0;i<rawNumTopYReco;i++) {
	Int_t pl=yPlaneTop[i];
	Int_t st=yStripTop[i];
	if(pl==testPlane) {
	  stripsHit[countHits]=st;
	  index[countHits]=i;
	  countHits++;
	}
      }
      TMath::Sort(countHits,stripsHit,indexIndex);
      Int_t lastStrip=-1;
      for(int hit=0;hit<countHits;hit++) {
	Int_t i=index[indexIndex[hit]];
	Int_t st=yStripTop[i];
	Int_t pl=yPlaneTop[i];
	Double_t y=yPosTopReco[i];
	Double_t z=zyPosTopReco[i];
	Double_t w=yErrTopReco[i];
	if(st==lastStrip-1) {
	  Double_t newW=yTopAvgW[countYTop-1]+w;
	  Double_t newY=(yTopAvgY[countYTop-1]*yTopAvgW[countYTop-1])+(y*w);
	  Double_t newZ=(yTopAvgZ[countYTop-1]*yTopAvgW[countYTop-1])+(z*w);
	  newY/=newW;
	  newZ/=newW;
	  yTopAvgW[countYTop-1]=newW;
	  yTopAvgZ[countYTop-1]=newZ;
	  yTopAvgY[countYTop-1]=newY;	  
	}
	else {
	  yTopAvgY[countYTop]=y;
	  yTopAvgZ[countYTop]=z;
	  yTopAvgW[countYTop]=w;
	  countYTop++;
	}
	lastStrip=st;
	//	std::cout << pl << "\t" << st << "\t"  <<  y 
	//		<< "\t" << z << "\t" << w << "\n";
      }
    }

    //xBot
   for(int testPlane=planesPerSide;testPlane<2*planesPerSide;testPlane+=2) {
      Int_t countHits=0;
      Int_t stripsHit[MAX_SCINT_HITS]={0};
      Int_t index[MAX_SCINT_HITS]={0};
      Int_t indexIndex[MAX_SCINT_HITS]={0};
      for(int i=0;i<rawNumBotXReco;i++) {
	Int_t pl=xPlaneBot[i];
	Int_t st=xStripBot[i];
	if(pl==testPlane) {
	  stripsHit[countHits]=st;
	  index[countHits]=i;
	  countHits++;
	}
      }
      TMath::Sort(countHits,stripsHit,indexIndex);
      Int_t lastStrip=-1;
      for(int hit=0;hit<countHits;hit++) {
	Int_t i=index[indexIndex[hit]];
	Int_t st=xStripBot[i];
	Int_t pl=xPlaneBot[i];
	Double_t x=xPosBotReco[i];
	Double_t z=zxPosBotReco[i];
	Double_t w=xErrBotReco[i];
	if(st==lastStrip-1) {
	  Double_t newW=xBotAvgW[countXBot-1]+w;
	  Double_t newX=(xBotAvgX[countXBot-1]*xBotAvgW[countXBot-1])+(x*w);
	  Double_t newZ=(xBotAvgZ[countXBot-1]*xBotAvgW[countXBot-1])+(z*w);
	  newX/=newW;
	  newZ/=newW;
	  xBotAvgW[countXBot-1]=newW;
	  xBotAvgZ[countXBot-1]=newZ;
	  xBotAvgX[countXBot-1]=newX;	  
	}
	else {
	  xBotAvgX[countXBot]=x;
	  xBotAvgZ[countXBot]=z;
	  xBotAvgW[countXBot]=w;
	  countXBot++;
	}
	lastStrip=st;
	//	std::cout << pl << "\t" << st << "\t"  <<  x 
	//		<< "\t" << z << "\t" << w << "\n";
      }
    }

    //Next up same thing for yBot
    for(int testPlane=planesPerSide+1;testPlane<2*planesPerSide;testPlane+=2) {
      Int_t countHits=0;
      Int_t stripsHit[MAX_SCINT_HITS]={0};
      Int_t index[MAX_SCINT_HITS]={0};
      Int_t indexIndex[MAX_SCINT_HITS]={0};
      for(int i=0;i<rawNumBotYReco;i++) {
	Int_t pl=yPlaneBot[i];
	Int_t st=yStripBot[i];
	if(pl==testPlane) {
	  stripsHit[countHits]=st;
	  index[countHits]=i;
	  countHits++;
	}
      }
      TMath::Sort(countHits,stripsHit,indexIndex);
      Int_t lastStrip=-1;
      for(int hit=0;hit<countHits;hit++) {
	Int_t i=index[indexIndex[hit]];
	Int_t st=yStripBot[i];
	Int_t pl=yPlaneBot[i];
	Double_t y=yPosBotReco[i];
	Double_t z=zyPosBotReco[i];
	Double_t w=yErrBotReco[i];
	if(st==lastStrip-1) {
	  Double_t newW=yBotAvgW[countYBot-1]+w;
	  Double_t newY=(yBotAvgY[countYBot-1]*yBotAvgW[countYBot-1])+(y*w);
	  Double_t newZ=(yBotAvgZ[countYBot-1]*yBotAvgW[countYBot-1])+(z*w);
	  newY/=newW;
	  newZ/=newW;
	  yBotAvgW[countYBot-1]=newW;
	  yBotAvgZ[countYBot-1]=newZ;
	  yBotAvgY[countYBot-1]=newY;	  
	}
	else {
	  yBotAvgY[countYBot]=y;
	  yBotAvgZ[countYBot]=z;
	  yBotAvgW[countYBot]=w;
	  countYBot++;
	}
	lastStrip=st;
	//	std::cout << pl << "\t" << st << "\t"  <<  y 
	//		<< "\t" << z << "\t" << w << "\n";
      }
    }
 //    for(int i=0;i<countYBot;i++) {
//       std::cout << yBotAvgY[i] << "\t" << yBotAvgZ[i] << "\t"
// 		<< yBotAvgW[i] << "\n";
//     }
      



    Int_t topPlanes=0;
    Int_t botPlanes=0;
    for(int plane=0;plane<planesPerSide;plane++) {
      if(gotPlane[plane]) topPlanes++;
      if(gotPlane[plane+planesPerSide]) botPlanes++;
    }
    std::cout << "Planes:\t" << topPlanes << "\t" << botPlanes << "\n";
    if(topPlanes>=4) {

      //Now find the gradients using the true positions
      Double_t gradxz[2];
      Double_t gradyz[2];
      //      fitxyz(topPlanes,xTop,yTop,zTop,gradxz,gradyz);
      //      fitxyzw(rawNumTop,xPosTop,yPosTop,zPosTop,xErrTop,gradxz,gradyz);
      iterativelyFitXYZWithWeights(rawNumTop,xPosTop,yPosTop,zPosTop,xErrTop,gradxz,gradyz);
      xzGradTrue[0]=gradxz[0];
      yzGradTrue[0]=gradyz[0];
      xzCutTrue[0]=gradxz[1];
      yzCutTrue[0]=gradyz[1];
       
      //Now find the gradients using the reco positions
      Double_t recoGradXZ[2];
      Double_t recoGradYZ[2];
      iterativelyFitXZWithWeights(countXTop,xTopAvgX,xTopAvgZ,xTopAvgW,recoGradXZ);
      iterativelyFitXZWithWeights(countYTop,yTopAvgY,yTopAvgZ,yTopAvgW,recoGradYZ);
      //    findGradients(topPlanesXView,xTopReco,xzTopReco,recoGradXZ);
      //       findGradients(topPlanesYView,yTopReco,yzTopReco,recoGradYZ);
      xzGradReco[0]=recoGradXZ[0];
      yzGradReco[0]=recoGradYZ[0];
      xzCutReco[0]=recoGradXZ[1];
      yzCutReco[0]=recoGradYZ[1];
      


      TCanvas *canTop = new TCanvas("canTop","canTop",600,600);
      canTop->Divide(1,2);
      canTop->cd(1);
      TGraphErrors *grTopTruePosX = new TGraphErrors(rawNumTop,zPosTop,xPosTop,0,xErrTop);
      grTopTruePosX->SetMarkerStyle(22);
      grTopTruePosX->SetMarkerColor(50);
      grTopTruePosX->Draw("ap");
      TGraph *grTopRecoPosX = new TGraph(countXTop,xTopAvgZ,xTopAvgX);
      grTopRecoPosX->SetMarkerStyle(29);
      grTopRecoPosX->SetMarkerColor(30);
      grTopRecoPosX->Draw("p");
      canTop->cd(2);
      TGraphErrors *grTopTruePosY = new TGraphErrors(rawNumTop,zPosTop,yPosTop,0,yErrTop);
      grTopTruePosY->SetMarkerStyle(22);
      grTopTruePosY->SetMarkerColor(50);
      grTopTruePosY->Draw("ap");
      TGraph *grTopRecoPosY = new TGraph(countYTop,yTopAvgZ,yTopAvgY);
      grTopRecoPosY->SetMarkerStyle(29);
      grTopRecoPosY->SetMarkerColor(30);
      grTopRecoPosY->Draw("p");
      TF1 *topPol1X = new TF1("topPol1X","pol1",-8000,+8000);
      TF1 *topPol1Y = new TF1("topPol1Y","pol1",-8000,+8000);
      TF1 *topPol1XReco = new TF1("topPol1XReco","pol1",-8000,+8000);
      TF1 *topPol1YReco = new TF1("topPol1YReco","pol1",-8000,+8000);
      canTop->cd(1);
      topPol1X->SetParameters(xzCutTrue[0],xzGradTrue[0]);
      topPol1X->SetLineColor(50);
      topPol1X->Draw("same");
      topPol1XReco->SetParameters(xzCutReco[0],xzGradReco[0]);
      topPol1XReco->SetLineColor(8);
      topPol1XReco->Draw("same");
      canTop->cd(2);
      topPol1Y->SetParameters(yzCutTrue[0],yzGradTrue[0]);
      topPol1Y->SetLineColor(50);
      topPol1Y->Draw("same");
      topPol1YReco->SetParameters(yzCutReco[0],yzGradReco[0]);
      topPol1YReco->SetLineColor(8);
      topPol1YReco->Draw("same");
      


      if(botPlanes>3) {
	//Potential scatter
	//	fitxyz(botPlanes,xBot,yBot,zBot,gradxz,gradyz);
	iterativelyFitXYZWithWeights(rawNumBot,xPosBot,yPosBot,zPosBot,xErrBot,gradxz,gradyz);
	xzGradTrue[1]=gradxz[0];
	yzGradTrue[1]=gradyz[0];
	xzCutTrue[1]=gradxz[1];
	yzCutTrue[1]=gradyz[1];


	TCanvas *canBot = new TCanvas("canBot","canBot",600,600);
	canBot->Divide(1,2);
	canBot->cd(1);
	TGraphErrors *grBotTruePosX = new TGraphErrors(rawNumBot,zPosBot,xPosBot,0,xErrBot);
	grBotTruePosX->SetMarkerStyle(22);
	grBotTruePosX->SetMarkerColor(50);
	grBotTruePosX->Draw("ap");
	TGraph *grBotRecoPosX = new TGraph(rawNumBotXReco,zxPosBotReco,xPosBotReco);
	grBotRecoPosX->SetMarkerStyle(29);
	grBotRecoPosX->SetMarkerColor(30);
	grBotRecoPosX->Draw("p");
	canBot->cd(2);
	TGraphErrors *grBotTruePosY = new TGraphErrors(rawNumBot,zPosBot,yPosBot,0,yErrBot);
	grBotTruePosY->SetMarkerStyle(22);
	grBotTruePosY->SetMarkerColor(50);
	grBotTruePosY->Draw("ap");
	TGraph *grBotRecoPosY = new TGraph(rawNumBotYReco,zxPosBotReco,xPosBotReco);
	grBotRecoPosY->SetMarkerStyle(29);
	grBotRecoPosY->SetMarkerColor(30);
	grBotRecoPosY->Draw("p");
	TF1 *botPol1X = new TF1("botPol1X","pol1",-8000,+8000);
	TF1 *botPol1Y = new TF1("botPol1Y","pol1",-8000,+8000);
	canBot->cd(1);
	botPol1X->SetParameters(xzCutTrue[1],xzGradTrue[1]);
	botPol1X->Draw("same");
	canBot->cd(2);
	botPol1Y->SetParameters(yzCutTrue[1],yzGradTrue[1]);
	botPol1Y->Draw("same");





	{
	  float a1 = xzGradTrue[0];
	  float a2 = xzGradTrue[1];
	    
	  float modp = sqrt(1.0+a1*a1);
	  float modq = sqrt(1.0+a2*a2);
	    
	  float pdotq = 1.0/sqrt(modp*modq)+(a1*a2)/sqrt(modp*modq);
	    
	  thetaxzTrue = acos(pdotq/(modp*modq));
	    
	  a1 = yzGradTrue[0];
	  a2 = yzGradTrue[1];

	  modp = sqrt(1.0+a1*a1);
	  modq = sqrt(1.0+a2*a2);
	    
	  pdotq = 1.0/sqrt(modp*modq)+(a1*a2)/sqrt(modp*modq);
	    
	  thetayzTrue = acos(pdotq/(modp*modq));
	}

	//	findGradients(botPlanesXView,xBotReco,xzBotReco,recoGradXZ);
	//	findGradients(botPlanesYView,yBotReco,yzBotReco,recoGradYZ);
	xzGradReco[1]=recoGradXZ[0];
	yzGradReco[1]=recoGradYZ[0];
	xzCutReco[1]=recoGradXZ[1];
	yzCutReco[1]=recoGradYZ[1];
          
	{
	  float a1 = xzGradReco[0];
	  float a2 = xzGradReco[1];
	    
	  float modp = sqrt(1.0+a1*a1);
	  float modq = sqrt(1.0+a2*a2);
	    
	  float pdotq = 1.0/sqrt(modp*modq)+(a1*a2)/sqrt(modp*modq);
	    
	  thetaxzReco = acos(pdotq/(modp*modq));
	    
	  a1 = yzGradReco[0];
	  a2 = yzGradReco[1];

	  modp = sqrt(1.0+a1*a1);
	  modq = sqrt(1.0+a2*a2);
	    
	  pdotq = 1.0/sqrt(modp*modq)+(a1*a2)/sqrt(modp*modq);
	    
	  thetayzReco = acos(pdotq/(modp*modq));
	}
	//	  std::cout << "RecoXZ: " << xzGradReco[0] << "\t" << xzGradReco[1] << "\n";

	Double_t pcaTrue[3];
	PointOfCloseApproach(xzGradTrue,yzGradTrue,xzCutTrue,yzCutTrue,pcaTrue);
	//	CloseApproach(xzGradTrue,yzGradTrue,xzCutTrue,yzCutTrue,pcaTrue);
	xPosTrue = pcaTrue[0];
	yPosTrue = pcaTrue[1];
	zPosTrue = pcaTrue[2];


	TCanvas *canXZ = new TCanvas("canXZ","canXZ",600,600);
	TMultiGraph *mgX = new TMultiGraph();
	mgX->Add(grTopTruePosX,"p");
	mgX->Add(grBotTruePosX,"p");
	mgX->Draw("ap");
	topPol1X->SetLineColor(8);
	topPol1X->Draw("same");
	botPol1X->SetLineColor(9);
	botPol1X->Draw("same");
	TGraph *grPCAX = new TGraph(1,&zPosTrue,&xPosTrue);
	grPCAX->SetMarkerStyle(29);
	grPCAX->SetMarkerColor(1);
	grPCAX->Draw("p");


	TCanvas *canYZ = new TCanvas("canYZ","canYZ",600,600);
	TMultiGraph *mgY = new TMultiGraph();
	mgY->Add(grTopTruePosY,"p");
	mgY->Add(grBotTruePosY,"p");
	mgY->Draw("ap");
	topPol1Y->SetLineColor(8);
	topPol1Y->Draw("same");
	botPol1Y->SetLineColor(9);
	botPol1Y->Draw("same");
	TGraph *grPCAY = new TGraph(1,&zPosTrue,&yPosTrue);
	grPCAY->SetMarkerStyle(29);
	grPCAY->SetMarkerColor(1);
	grPCAY->Draw("p");

	TCanvas *canS = new TCanvas("canS","canS",600,600);
	TF1 *funcs = new TF1("funcs",getSSquared,-8000,8000,8);
	funcs->SetParameters(xzGradTrue[0],xzCutTrue[0],xzGradTrue[1],xzCutTrue[1],
			     yzGradTrue[0],yzCutTrue[0],yzGradTrue[1],yzCutTrue[1]);
	funcs->Draw("");

	TVector3 intialDirTrue(xzGradTrue[0],yzGradTrue[0],1);
	TVector3 finalDirTrue(xzGradTrue[1],yzGradTrue[1],1);
	thetaTrue=intialDirTrue.Angle(finalDirTrue);

	// 	  if(zPosTrue<-6000 && zPosTrue>-7000) {
	// 	    std::cout << "xzGradTrue: " << xzGradTrue[0] << "\t" << xzGradTrue[1] << std::endl;;
	// 	    std::cout << "xzCutTrue: " << xzCutTrue[0] << "\t" << xzCutTrue[1] << std::endl;;
	// 	    std::cout << "yzGradTrue: " << yzGradTrue[0] << "\t" << yzGradTrue[1] << std::endl;;
	// 	    std::cout << "yzCutTrue: " << yzCutTrue[0] << "\t" << yzCutTrue[1] << std::endl;;
	// 	    std::cout << "pca: " << xPosTrue << "\t" << yPosTrue << "\t" 
	// 		      << zPosTrue << std::endl;;

	// 	    Double_t ax=(xzGradTrue[0]-xzGradTrue[1])*(xzGradTrue[0]-xzGradTrue[1]);
	// 	    Double_t bx=2*(xzGradTrue[0]-xzGradTrue[1])*(xzCutTrue[0]-xzCutTrue[1]);
	// 	    std::cout << "ax = " << ax << " bx = " << bx << " zx = " << (-0.5*bx/ax)
	// 		      << "\n";

	// 	    Double_t ay=(yzGradTrue[0]-yzGradTrue[1])*(yzGradTrue[0]-yzGradTrue[1]);
	// 	    Double_t by=2*(yzGradTrue[0]-yzGradTrue[1])*(yzCutTrue[0]-yzCutTrue[1]);
	// 	    std::cout << "ay = " << ay << " by = " << by << " zy = " << (-0.5*by/ay)
	// 		      << "\n";
	// 	    std::cout << "Planes: " << topPlanes << "\t" << botPlanes << "\n";
	// 	    std::cout << "Theta: " << thetaTrue*TMath::RadToDeg() << "\n";


	// 	    for(int plane=0;plane<planesPerSide*2;plane++) {	  
	// 	      if(gotPlane[plane]) {
	// 		std::cout << "x[" << plane%8 << "]=" << xTrue[plane] << ";"
	// 			  << "\ty[" << plane%8 << "]=" << yTrue[plane] << ";"
	// 			  << "\tz[" << plane%8 << "]=" << zTrue[plane] << ";\n";
	// 	      }
	// 	    }

	// 	  }

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


	gotRecoPCA=1;
	if(TMath::IsNaN(xPosReco)) {
	  gotRecoPCA=0;
	  //Generally this happens when the gradients from the top and bottom are identical.
	  // 	    std::cout << "Oh no we have a Nan\n";
	  // 	    std::cout << xPosReco << "\t" << yPosReco << "\t" << zPosReco << "\n";
	  // 	    std::cout << "Gradients\t" << xzGradReco[0] << "\t" << xzGradReco[1] << "\n";
	  // 	    std::cout << "Intercepts\t" << xzCutReco[0] << "\t" << xzCutReco[1] << "\n";
	  // 	    std::cout << "True Gradients\t" << xzGradTrue[0] << "\t" << xzGradTrue[1] << "\n";
	  // 	    std::cout << "True Intercepts\t" << xzCutTrue[0] << "\t" << xzCutTrue[1] << "\n";
	}

	 
	//	  pcaTree->Fill();

	//	  PCAh->Fill(pcaTrue[0],pcaTrue[1],pcaTrue[2]);
	  
	
      }
      else {
	//Potential Absorbtion ... need to tweak them
	// 	  AbsorbedHisto->Fill(xTop[0],yTop[0]);

	AGradx=xzGradTrue[0];
	AGrady=yzGradTrue[0];
	ACutx=xzCutTrue[0];
	ACuty=yzCutTrue[0];
	// 	  Absorbed ->Fill();
      }

    }
  }

  

  //    newfile->cd();
  //  Absorbed->AutoSave();
  //  pcaTree->AutoSave();
  //  newfile->Close();

}


