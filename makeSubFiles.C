

void makeSubFiles()
{
   char filename[180];
   char outname[180];
   ifstream Missing("missingFiles.txt");


   int tagCount=1;
   int fileCount=0;
   ofstream Out;
   while(Missing >> filename) {
      if(fileCount%25==0) {
	 if(Out.is_open()) {
	    Out.close();
	 }
	 sprintf(outname,"file%d.txt",tagCount);
	 Out.open(outname);
	 tagCount++;
      }
      Out << filename << endl;
      fileCount++;      
   }

}
