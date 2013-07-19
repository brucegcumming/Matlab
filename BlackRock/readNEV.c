#include <math.h>
#include "mex.h"
#include <stdio.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  char * filename, * wd;
 FILE * fid;
 int buflen;

 int bytesPerSample = 2;
 unsigned int numSamples;

 char fileType[8];
 unsigned char version[2];
 char fileFormatAdditional[2];
 unsigned int headerSize;
 unsigned int timeResTimeStamps;
 unsigned int timeResSamples;
 unsigned int packetSize;
 unsigned short timeOrigin[8];
 char application[32];
 unsigned char comment[256];
 int extendedHeaderNumber;
 int spikeCount,i;
 double* nevData;

 unsigned char junk[500];

 unsigned int time;
 unsigned char unit;
 short packetID;
 
 /* Get filename */
 buflen = (mxGetM(prhs[0]) * mxGetN(prhs[0]) * sizeof(mxChar))+1;
 filename = mxCalloc(buflen, sizeof(char));
 mxGetString(prhs[0], filename, buflen);

 if ((fid = fopen(filename,"rb")) == NULL)
 {
   wd = (char *)mxCalloc(100,1);
   getcwd(wd,100);
   mexPrintf("File %s/%s does not exist\n",wd,filename);
   plhs[0] = mxCreateNumericMatrix(0,0,mxDOUBLE_CLASS,mxREAL);
   return;
 }

 fseek(fid,0,SEEK_SET);

 /* Read the header */
 fread(fileType, 1, 8, fid);
 fread(version, 1, 2, fid);
 fread(fileFormatAdditional, 1, 2, fid);
 fread(&headerSize, 4, 1, fid);
 fread(&packetSize, 4, 1, fid);
 fread(&timeResTimeStamps, 4, 1, fid);
 fread(&timeResSamples, 4, 1, fid);
 fread(timeOrigin, 2, 8, fid);
 fread(application, 1, 32, fid);
 fread(comment, 1, 256, fid);
 fread(&extendedHeaderNumber,4,1,fid);
 numSamples = (packetSize-8)/bytesPerSample;

 fseek(fid,0,SEEK_END);
 spikeCount = (ftell(fid) - headerSize)/packetSize;
 
 mexPrintf("Reading %s, %li events...\n",filename,spikeCount);

 fseek(fid,headerSize,SEEK_SET);

 plhs[0] = mxCreateNumericMatrix(spikeCount,3,mxDOUBLE_CLASS,mxREAL);
 nevData = mxGetPr(plhs[0]);

 for (i = 0; i < spikeCount; i++)
 {
   fread(&time,4,1,fid);
   fread(&packetID,2,1,fid);
   fread(&unit,1,1,fid);
   
   if (packetID == 0)
   {
     fread(junk,1,1,fid);
     fread(&unit,1,1,fid);
     fread(junk,packetSize-9,1,fid);   
   }
   else
   {
     fread(junk,packetSize-7,1,fid);
   }
 
   *(nevData+i+spikeCount*2) = ((double)time)/timeResSamples;
   *(nevData+i) = packetID;
   *(nevData+i+spikeCount) = unit;     

   if (i % (spikeCount/10) == 0)
   {
     mexPrintf("  %li/%li\n",i,spikeCount);
     mexEvalString("drawnow;");
   }
 }

 fclose(fid);
}

