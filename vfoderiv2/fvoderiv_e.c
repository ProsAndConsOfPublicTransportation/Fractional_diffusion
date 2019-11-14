
#define S_FUNCTION_NAME  fvoderiv_e
#define S_FUNCTION_LEVEL 2



#define MDL_START


#define N(I,J)  (*(mxGetPr(ssGetSFcnParam(S,0))+(mxGetM(ssGetSFcnParam(S,0)))*(J)+(I)))


#define Nu (mxGetM(ssGetSFcnParam(S,0)))
#define Nbuf ((int)(*(mxGetPr(ssGetSFcnParam(S,1)))))
#define Ts (*(mxGetPr(ssGetSFcnParam(S,2))))

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "simstruc.h"












static void mdlInitializeSizes(SimStruct *S)
{

ssSetNumSFcnParams(S, 3);
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
        return; /* Parameter mismatch will be reported by Simulink */
    }

   

    if (mxGetN(ssGetSFcnParam(S,0))!=1){
    printf("Matrix N must be an one column vector\n");
    return;
    }
    

    
    if (!ssSetNumInputPorts(S, 2)) return;
    ssSetInputPortWidth(S, 0, Nu);
    ssSetInputPortWidth(S, 1, Nu);
    ssSetInputPortDirectFeedThrough(S, 0, 1);
    ssSetInputPortDirectFeedThrough(S, 1, 1);
    if (!ssSetNumOutputPorts(S,1)) return;
    ssSetOutputPortWidth(S, 0, Nu);
    ssSetNumSampleTimes(S, 1);

    ssSetNumPWork(S,4); // X Delta wsp tempalf vectors 
    ssSetNumIWork(S,1); // px 

    
    
    
    
    /* Take care when specifying exception free code - see sfuntmpl_doc.c */
    ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE );
}



/* Function: mdlInitializeSampleTimes =========================================
 * Abstract:
 *    Specifiy that we inherit our sample time from the driving block.
 */
static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, Ts);
    ssSetOffsetTime(S, 0, 0.0);
}



static void mdlStart(SimStruct *S)
{
int_T i,j;
real_T **tempx;
real_T **tempalf;
real_T *Delta=calloc(Nu,sizeof(real_T));
real_T *wsp=calloc(Nu,sizeof(real_T));



  

    tempx=(real_T**)calloc(Nbuf,sizeof(real_T*));
    if (tempx==NULL){
    printf("Can't allocate memory");
    }
    for (i=0;i<Nbuf;i++){
    tempx[i]=(real_T*)calloc((Nu),sizeof(real_T));
    if (tempx[i]==NULL){
    printf("Can't allocate memory");
    }
    
    for (j=0;j<Nu;j++)
    tempx[i][j]=0; 
    }
    
    tempalf=(real_T**)calloc(Nbuf,sizeof(real_T*));
    if (tempalf==NULL){
    printf("Can't allocate memory");
    }
    for (i=0;i<Nbuf;i++){
    tempalf[i]=(real_T*)calloc((Nu),sizeof(real_T));
    if (tempalf[i]==NULL){
    printf("Can't allocate memory");
    }
    
    for (j=0;j<Nu;j++)
    tempalf[i][j]=0; 
    }
    
    ssSetPWorkValue(S,0,tempx);
    ssSetPWorkValue(S,1,Delta);
    ssSetPWorkValue(S,2,wsp);
    ssSetPWorkValue(S,3,tempalf);
    
    ssSetIWorkValue(S,0,0);

  
    
 
}



real_T xwy(SimStruct *S,real_T **X,int_T k,int_T px,int_T nr){
real_T xw;  
int_T p=px-k;


if (p>=Nbuf){
    xw=X[p-Nbuf][nr];
    }
else
if (p>=0){
    xw=X[p][nr];
    }
else {
    xw=X[Nbuf+p][nr];
    }

 return xw;
 }






/* Function: mdlOutputs =======================================================
 */
 
static void mdlOutputs(SimStruct *S, int_T tid)
{
 
    
    int_T               i,j,k;
    InputRealPtrsType   uPtrs = ssGetInputPortRealSignalPtrs(S,0);
    InputRealPtrsType   uPtrs2 = ssGetInputPortRealSignalPtrs(S,1);
    real_T              *y    = ssGetOutputPortRealSignal(S,0);
    int_T               width = ssGetOutputPortWidth(S,0);
    int_T               px=ssGetIWorkValue(S,0);
    real_T              **U,**Alf;
    real_T *Delta=(real_T*)ssGetPWorkValue(S,1);
    real_T *wsp=(real_T*)ssGetPWorkValue(S,2);
  
   U=(real_T**)ssGetPWorkValue(S,0);
    
    Alf=(real_T**)ssGetPWorkValue(S,3);
    
for (i=0;i<Nu;i++)
{
        Alf[px][i]=(*uPtrs2[i]);
   Delta[i]=*uPtrs[i]/pow(Ts,(*uPtrs2[i]));
   //Delta[i]=*uPtrs[i]/pow(Ts,xwy(S,Alf,1,px,i));
   
   for (j=1;j<Nbuf;j++)
   {wsp[i]=1;
    for (k=1;k<=j;k++){
            wsp[i]=wsp[i]*(-xwy(S,Alf,j,px,i)-k+1)/k;
    }
           Delta[i]-=(pow(Ts,xwy(S,Alf,j,px,i))/pow(Ts,(*uPtrs2[i])))*pow(-1,j)*wsp[i]*xwy(S,U,j,px,i);
   }
   
   U[px][i]=Delta[i];
   
}
    
    
    
    
for (i=0;i<Nu;i++)
{
    y[i]=Delta[i];

}




if (++px>=Nbuf){
    px=0;
   }
 
ssSetIWorkValue(S,0,px);

}


/* Function: mdlTerminate =====================================================
 * Abstract:
 *    No termination needed, but we are required to have this routine.
 */
static void mdlTerminate(SimStruct *S)
{
int_T i;
real_T **temp=(real_T**)ssGetPWorkValue(S,0);
real_T *Delta=(real_T*)ssGetPWorkValue(S,1);
real_T *wsp=(real_T*)ssGetPWorkValue(S,2);



for (i=0;i<Nbuf;i++)
free(temp[i]);
free(temp);
free(Delta);
free(wsp);

}



#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif

