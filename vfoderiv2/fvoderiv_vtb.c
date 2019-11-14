/*
C-MEX S-function code for fractional variable order derivative according to the second type of definition.
Authors: Dominik Sierociuk, Wiktor Malesza and Michal Macias
Institute of Control and Industrial Electronics,
Warsaw University of Technology, Koszykowa 75, Warsaw, Poland 
(e-mail: dsieroci,wmalesza,michal.macias@ee.pw.edu.pl)
 This file is a part of Fractional Variable Order Derivative Simulink Toolkit
 */
#define S_FUNCTION_NAME  fvoderiv_vtb
#define S_FUNCTION_LEVEL 2



#define MDL_START


#define N(I,J)  (*(mxGetPr(ssGetSFcnParam(S,0))+(mxGetM(ssGetSFcnParam(S,0)))*(J)+(I)))


#define Nu ((int)(*(mxGetPr(ssGetSFcnParam(S,0)))))
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
    

    
   if (!ssSetNumInputPorts(S, 4)) return;
    ssSetInputPortWidth(S, 0, Nu);
    ssSetInputPortWidth(S, 1, Nu);
    ssSetInputPortWidth(S, 2, Nu);
    ssSetInputPortWidth(S, 3, Nu);
    ssSetInputPortDirectFeedThrough(S, 0, 1);
    ssSetInputPortDirectFeedThrough(S, 1, 1);
    ssSetInputPortDirectFeedThrough(S, 2, 1);
    ssSetInputPortDirectFeedThrough(S, 3, 1);
    if (!ssSetNumOutputPorts(S,1)) return;
    ssSetOutputPortWidth(S, 0, Nu);
    ssSetNumSampleTimes(S, 1);

    ssSetNumPWork(S,5); // X Delta typ Alf param vectors
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
real_T **Alf;
real_T *Delta=calloc(Nu,sizeof(real_T));
real_T **typ,**param; //calloc(Nu,sizeof(real_T));

    typ=(real_T**)calloc(Nbuf,sizeof(real_T*));
    if (typ==NULL){
    printf("Can't allocate memory");
    }
    for (i=0;i<Nbuf;i++){
    typ[i]=(real_T*)calloc((Nu),sizeof(real_T));
    if (typ[i]==NULL){
    printf("Can't allocate memory");
    }
    for (j=0;j<Nu;j++)
    typ[i][j]=1; 
    }
    
    param=(real_T**)calloc(Nbuf,sizeof(real_T*));
    if (param==NULL){
    printf("Can't allocate memory");
    }
    for (i=0;i<Nbuf;i++){
    param[i]=(real_T*)calloc((Nu),sizeof(real_T));
    if (param[i]==NULL){
    printf("Can't allocate memory");
    }
    for (j=0;j<Nu;j++)
    param[i][j]=0; 
    }
    
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
    
    
    Alf=(real_T**)calloc(Nbuf,sizeof(real_T*));
    if (Alf==NULL){
    printf("Can't allocate memory");
    }
    for (i=0;i<Nbuf;i++){
    Alf[i]=(real_T*)calloc((Nu),sizeof(real_T));
    if (Alf[i]==NULL){
    printf("Can't allocate memory");
    }
    for (j=0;j<Nu;j++)
    Alf[i][j]=0; 
    }
    
    ssSetPWorkValue(S,0,tempx);
    ssSetPWorkValue(S,1,Delta);
    ssSetPWorkValue(S,2,typ);
    ssSetPWorkValue(S,3,Alf);
    ssSetPWorkValue(S,4,param);
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

void xwe(SimStruct *S,real_T **X,int_T k,int_T px,int_T nr,real_T xw){
int_T p=px-k;

if (p>=Nbuf){
    X[p-Nbuf][nr]=xw;
    }
else
if (p>=0){
    X[p][nr]=xw;
    }
else {
    X[Nbuf+p][nr]=xw;
    }
 }




/* Function: mdlOutputs =======================================================
 */
 
static void mdlOutputs(SimStruct *S, int_T tid)
{     
    int_T               i,j,k;
    InputRealPtrsType   uPtrs = ssGetInputPortRealSignalPtrs(S,0);
    InputRealPtrsType   uPtrs2 = ssGetInputPortRealSignalPtrs(S,1);
    InputRealPtrsType   uTyp = ssGetInputPortRealSignalPtrs(S,2);
    InputRealPtrsType   uParam = ssGetInputPortRealSignalPtrs(S,3);
    real_T              *y    = ssGetOutputPortRealSignal(S,0);
    int_T               width = ssGetOutputPortWidth(S,0);
    int_T               px=ssGetIWorkValue(S,0);
    real_T              **U,**Alf;
    real_T *Delta=(real_T*)ssGetPWorkValue(S,1);
    real_T **typ,wsp,**param;
  
   U=(real_T**)ssGetPWorkValue(S,0);
   Alf=(real_T**)ssGetPWorkValue(S,3);
   typ=(real_T**)ssGetPWorkValue(S,2);
   param=(real_T**)ssGetPWorkValue(S,4);
   
for (i=0;i<Nu;i++)
{
   U[px][i]=*uPtrs[i];
   Alf[px][i]=*uPtrs2[i];
   typ[px][i]=*uTyp[i];
   param[px][i]=*uParam[i];
   
    
   Delta[i]=param[px][i]*U[px][i]/pow(Ts,(Alf[px][i]));            
   for (j=1;j<Nbuf;j++){
     wsp=1;    
     if (xwy(S,typ,j,px,i)==2){    
     for (k=1;k<=j;k++){
        wsp=wsp*(xwy(S,Alf,j,px,i)-k+1)/k;
       }
     Delta[i]+=pow(-1,j)*wsp*xwy(S,param,j,px,i)*xwy(S,U,j,px,i)/pow(Ts,(xwy(S,Alf,j,px,i))); 
     }else{
       for (k=1;k<=j;k++){
        wsp=wsp*(Alf[px][i]-k+1)/k;
       }
       Delta[i]+=pow(-1,j)*wsp*param[px][i]*xwy(S,U,j,px,i)/pow(Ts,(Alf[px][i]));  
       //Delta[i]+=pow(-1,j)*wsp*xwy(S,param,j,px,i)*xwy(S,U,j,px,i)/pow(Ts,(Alf[px][i]));  
     }
     //xwe(S,wsp,j,px,i,xwy(S,wsp,j,px,i)*(xwy(S,Alf,j,px,i)-j+1)/(j));
               
   }
   
      
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
real_T **typ=(real_T**)ssGetPWorkValue(S,2);
real_T **Alf=(real_T**)ssGetPWorkValue(S,3);
real_T **param=(real_T**)ssGetPWorkValue(S,4);

for (i=0;i<Nbuf;i++)
{
    free(temp[i]);
    free(Alf[i]);
    free(typ[i]);
    free(param[i]);
}
    free(temp);
    free(Delta);
    free(typ);
    free(Alf);
    free(param);
}



#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif

