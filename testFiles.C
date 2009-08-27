

void testFiles() {
   char inputFile[180];
   char *inputDir="/unix/anita1/creamtea/strips_650/fakecontainer_10cmtargetat_0p5_3_0p5/";
   char *rootName="fakecontainer_10cmtarget_";

  FileStat_t staty;
  
  ofstream Missing("missingFiles.txt");


  for(int fileNum=1;fileNum<=10000;fileNum++) {
    if(fileNum%100==0) cerr<< "*";
     sprintf(inputFile,"%s/%s%d.root",inputDir,rootName,fileNum);
     if(gSystem->GetPathInfo(inputFile,staty)) {
	std::cout << fileNum << "\n";
	Missing << inputFile << "\n";
	continue;
     }
     
     TFile fp(inputFile,"READ");
     if(fp.IsOpen()) {           
	TTree* scintTree = (TTree*) fp.Get("scintTree");
	if(scintTree) {
	   if(scintTree->GetEntries()>800) {
	      
	   }
	}
	else {
	   std::cout << fileNum << "\n";
	   Missing << inputFile << "\n";
	}
	fp.Close();	
     }
     else {
	std::cout << fileNum << "\n";
	Missing << inputFile << "\n";
     }
  }


}
