/*Compile the mex file:
 *
        mex -O cMultiFrameConnect_allcalcs3D.c
 *
 *      [coords]=cMultiFrameConnect_alcalcs3D(coords,sigma,maxdistxy,maxdistz,nframes)
 */

#include <stdlib.h>
#include <string.h>
#include "mex.h"
#include <math.h>

#define max(a,b) ( (a) >= (b) ? (a) : (b) )  
#define min(a,b) ( (a) < (b) ? (a) : (b) )  
#define FLOAT double
#define MLTYPENAME mxDOUBLE_CLASS

void mexFunction(int nlhs, mxArray *plhs[],	int	nrhs, const	mxArray	*prhs[])
{
    int N,tt;
    FLOAT sigma,maxdistxy,maxdistz;
    /* int outsize[3]; */
    long long unsigned int outsize[2];
    int nframes;
    FLOAT x1std,x2std,tmp;
    FLOAT y1std,y2std;
    FLOAT z1std,z2std;
    FLOAT xvar,xmean,yvar,ymean,zmean,zvar;
    FLOAT xii,xjj,yii,yjj,zii,zjj;
    FLOAT xvarii,xvarjj,yvarii,yvarjj,zvarii,zvarjj;
    FLOAT Nphii,Nphjj,bgii,bgjj;
    FLOAT runlengthii,runlengthjj;
    FLOAT locindexii,locindexjj;
    FLOAT xsumii,xsumjj,ysumii,ysumjj,zsumii,zsumjj;
    FLOAT xsqsumii,xsqsumjj,ysqsumii,ysqsumjj,zsqsumii,zsqsumjj;
    FLOAT crlbxii,crlbxjj,crlbx,crlbyii,crlbyjj,crlby,crlbzii,crlbzjj,crlbz;
    FLOAT *coords,*coordsin,*dataout;
    int cnt,kk,ii,jj;
    int numpars;
    const mwSize *datasize;
    
    /*input checks*/
    if (nrhs!=5)
        mexErrMsgTxt("Must input coords, sigma_connect, maxdistxy, maxdistz, nframes\n");
    
    datasize=mxGetDimensions(prhs[0]);
    
    /*input checks*/
    numpars = 20;
    if (datasize[1]!=numpars)
        mexErrMsgTxt("Size of coords must be Nx20.\n");
    
    N=datasize[0];
    
    coordsin=(FLOAT *) mxGetData(prhs[0]);
    sigma=(FLOAT)mxGetScalar(prhs[1]);
    maxdistxy=(FLOAT)mxGetScalar(prhs[2]);
    maxdistz=(FLOAT)mxGetScalar(prhs[3]);
    nframes=(FLOAT)mxGetScalar(prhs[4]);
    
    if (mxGetClassID(prhs[0])!=mxDOUBLE_CLASS)
        mexErrMsgTxt("data must be double\n");
    
    /*copy coordsin to coords*/
    coords=(FLOAT*)mxCalloc(N*numpars,sizeof(FLOAT));
    for (ii=0;ii<N*numpars;ii++)coords[ii]=coordsin[ii];
    
    cnt=0;
    for (ii=0;ii<N-1;ii++){
        tt=(int)coords[3*N+ii];
        for (jj=ii+1;jj<N;jj++){
           
            /*fast reject*/
            if ((int)coords[3*N+jj]==tt)continue;
            if ((int)coords[3*N+jj]>(tt+nframes))break;
            if (fabs(coords[0*N+jj]-coords[0*N+ii])>maxdistxy) continue;
            if (fabs(coords[1*N+jj]-coords[1*N+ii])>maxdistxy) continue;
            if (fabs(coords[2*N+jj]-coords[2*N+ii])>maxdistz) continue;
            if (sqrt( pow(coords[0*N+jj]-coords[0*N+ii],2)+ pow(coords[1*N+jj]-coords[1*N+ii],2) )> maxdistxy) continue;

            /*now compare to localization precision- separate dimensions*/
            x1std=sqrt(coords[4*N+ii]);
            x2std=sqrt(coords[4*N+jj]);
            /*tmp=sigma*(x1std+x2std)/2;*/
            tmp = sigma*fmax(x1std,x2std);
            if (fabs(coords[0*N+jj]-coords[0*N+ii])>tmp)continue;
            
            y1std=sqrt(coords[5*N+ii]);
            y2std=sqrt(coords[5*N+jj]);
            /*tmp=sigma*(y1std+y2std)/2;*/
            tmp = sigma*fmax(y1std,y2std);
            if (fabs(coords[1*N+jj]-coords[1*N+ii])>tmp)continue;
                    
            z1std=sqrt(coords[6*N+ii]);
            z2std=sqrt(coords[6*N+jj]);
            /*tmp=sigma*(z1std+z2std)/2;*/
            tmp = sigma*fmax(z1std,z2std);
            if (fabs(coords[2*N+jj]-coords[2*N+ii])>tmp)continue;
            
            /* It passes all tests, so combine this into the later coordinate
            (so it can be used later in the loop) */
            
            /*mexPrintf("found connection between %d and %d\n",ii,jj);*/
            
            /* calculation of the variables for the merged localization*/
            
            xii=coords[0*N+ii];
            xjj=coords[0*N+jj];
            yii=coords[1*N+ii];
            yjj=coords[1*N+jj];
            zii=coords[2*N+ii];
            zjj=coords[2*N+jj];
            xvarii=coords[4*N+ii];
            xvarjj=coords[4*N+jj];
            yvarii=coords[5*N+ii];
            yvarjj=coords[5*N+jj];
            zvarii=coords[6*N+ii];
            zvarjj=coords[6*N+jj];
            crlbxii=coords[7*N+ii];
            crlbxjj=coords[7*N+jj];
            crlbyii=coords[8*N+ii];
            crlbyjj=coords[8*N+jj];
            crlbzii=coords[9*N+ii];
            crlbzjj=coords[9*N+jj];
            Nphii=coords[10*N+ii];
            Nphjj=coords[10*N+jj];
            bgii=coords[11*N+ii];
            bgjj=coords[11*N+jj];
            runlengthii=coords[12*N+ii];
            runlengthjj=coords[12*N+jj];
            xsumii=coords[13*N+ii];
            xsumjj=coords[13*N+jj];
            ysumii=coords[14*N+ii];
            ysumjj=coords[14*N+jj];
            zsumii=coords[15*N+ii];
            zsumjj=coords[15*N+jj];
            xsqsumii=coords[16*N+ii];
            xsqsumjj=coords[16*N+jj];
            ysqsumii=coords[17*N+ii];
            ysqsumjj=coords[17*N+jj];
            zsqsumii=coords[18*N+ii];
            zsqsumjj=coords[18*N+jj];
            locindexii=coords[19*N+ii];
            locindexjj=coords[19*N+jj];
                                   
            /* compute the weighted mean position */
            xvar=1/(1/xvarii+ 1/xvarjj);
            yvar=1/(1/yvarii+ 1/yvarjj);
            zvar=1/(1/zvarii+ 1/zvarjj);
            xmean=xvar*(xii/xvarii+xjj/xvarjj);
            ymean=yvar*(yii/yvarii+yjj/yvarjj);
            zmean=zvar*(zii/zvarii+zjj/zvarjj);
            
            /* replace the variables for event jj with the ii & jj merged variables */
            coords[0*N+jj]=xmean;
            coords[1*N+jj]=ymean;
            coords[2*N+jj]=zmean;
            coords[4*N+jj]=xvar;
            coords[5*N+jj]=yvar;
            coords[6*N+jj]=zvar;
            coords[7*N+jj]=(runlengthii*crlbxii+runlengthjj*crlbxjj)/(runlengthii+runlengthjj);
            coords[8*N+jj]=(runlengthii*crlbyii+runlengthjj*crlbyjj)/(runlengthii+runlengthjj);
            coords[9*N+jj]=(runlengthii*crlbzii+runlengthjj*crlbzjj)/(runlengthii+runlengthjj);
            coords[10*N+jj]=Nphii+Nphjj;
            coords[11*N+jj]=bgii+bgjj;
            coords[12*N+jj]=runlengthii+runlengthjj;
            coords[13*N+jj]=xsumii+xsumjj;
            coords[14*N+jj]=ysumii+ysumjj;
            coords[15*N+jj]=zsumii+zsumjj;
            coords[16*N+jj]=xsqsumii+xsqsumjj;
            coords[17*N+jj]=ysqsumii+ysqsumjj;
            coords[18*N+jj]=zsqsumii+zsqsumjj;
            coords[19*N+jj]=locindexjj;
            
            /* set variables for event ii to dummy values */
            coords[0*N+ii]=0;
            coords[1*N+ii]=0;
            coords[2*N+ii]=0;
            coords[4*N+ii]=0;
            coords[5*N+ii]=0;
            coords[6*N+ii]=0;
            coords[7*N+ii]=0;
            coords[8*N+ii]=0;
            coords[9*N+ii]=0;
            coords[10*N+ii]=runlengthii;
            coords[11*N+ii]=0;
            coords[12*N+ii]=0;
            coords[13*N+ii]=0;
            coords[14*N+ii]=0;
            coords[15*N+ii]=0;
            coords[16*N+ii]=0;
            coords[17*N+ii]=0;
            coords[18*N+ii]=0;
            coords[19*N+ii]=locindexii;
            
            cnt++;
            break;
        }
    }
    
    outsize[0]=N-cnt;
    outsize[1]=numpars;
    plhs[0]=mxCreateNumericArray(2,outsize,mxDOUBLE_CLASS,mxREAL);
    dataout=(FLOAT *)mxGetData(plhs[0]);
    
    /*write into compressed coords array*/
    /*mexPrintf("Writing into output array: %d.\n",N-cnt);*/
    kk=0;
    for (ii=0;ii<N;ii++){
        if (coords[ii]){
            for (jj=0;jj<numpars;jj++)
                dataout[jj*(N-cnt)+kk]=coords[jj*N+ii];
                kk++;
            }
    }
    
     /*mexPrintf("kk: %d.\n",kk);*/
    
    mxFree(coords);
    return;
};

