
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <stdlib.h>

#define PRECISION float

#define  max 999 
#define  maa   9 

/*Block size depends on maximum threads*/
#define BLOCK_SIZE_X  16
#define BLOCK_SIZE_Y  16

/* Parameters */
int imax , jmax ;
int imax2 , jmax2 ;
int ima  , jma  ;
int isd  , jsd  ;
int ied  , jed  ;
int numGPUs ;

long ncye , nwri ;

PRECISION reyn , rtau , csou , rcsu ;
PRECISION pini , uini , vini , uwui , runi ;

PRECISION ex[maa] , ey[maa] , we[maa] ;
PRECISION fn[max][max][maa] , fe[max][max][maa] , fp[max][max][maa] ;
PRECISION pn[max][max] , un[max][max] , vn[max][max] ;
PRECISION xg[max][max] , yg[max][max] ;


#define PTR(i, j) (imax2 * j + i)
#define PTRQ(i, j, k) ((j*imax2*maa)+(i*maa)+k) 


__constant__ PRECISION dex[maa];
__constant__ PRECISION dey[maa];
__constant__ PRECISION dwe[maa];

struct Conditions {
	int imax, jmax ;
	int imax2, jmax2 ;
	int isd, ied, jsd, jed ;
	PRECISION rtau, rcsu, runi, uwui ;
};

//ŠÖ”‚Ìƒvƒƒgƒ^ƒCƒvéŒ¾//
void iniset(void) ;
void inicon(void) ;
void solver(void) ;


//********************************‘æˆêˆ—‹æ‰æ**********************************//
//
//
//******************************************************************************//
/*__global__ void matrix1(PRECISION *dfe, PRECISION* dfp, PRECISION* dfn,
			            struct Conditions *d_cond) */ // original definition

__global__ void matrix1(PRECISION *dfe, PRECISION* dfp, PRECISION* dfn, PRECISION *dmarker,
						PRECISION *dun, PRECISION *dvn, PRECISION *dpn,
			            struct Conditions *d_cond)
{
	const int i = blockIdx.x * blockDim.x + threadIdx.x  ;
	const int j = blockIdx.y * blockDim.y + threadIdx.y  ;

	//ƒOƒ[ƒoƒ‹•Ï”‚ÍGPU“à‚ÉŠ±Â‚Å‚«‚È‚¢//
	const int imax2 = d_cond->imax2 ;
	const PRECISION rtau = d_cond->rtau  ;
	
	const PRECISION runi = d_cond->runi;

	const int imax = d_cond->imax ;
	const int jmax = d_cond->jmax ;

	if(i==0){
	}
	else if(i==imax){
	}
	else if(j==0){
	}
	else if(j==jmax){
	}
	else{
		for(int k=0;k<maa;k++){
			int ii = i - int(dex[k]) ;
			int jj = j - int(dey[k]) ;
			const int pointer = PTRQ(ii,jj,k) ;
			const int pointer2 = PTR(ii,jj);
//add external forcing here
			/*	euv = dex[k] * dun[PTR(i,j)] + dey[k] * dvn[PTR(i,j)] ;
		qau = 0.5 * (dun[PTR(i,j)] * dun[PTR(i,j)] + dvn[PTR(i,j)] * dvn[PTR(i,j)]) ;
		dfe[PTRQ(i, j, k)] = dwe[k] * (dpn[PTR(i, j)] + runi * (euv + 1.5 * euv * euv - qau)) ;
		dfn[PTRQ(i, j, k)] = dfp[PTRQ(i, j, k)] ; */


//dfp[PTRQ(i, j, k)] = dfn[pointer]+ rtau * (dfe[pointer] - dfn[pointer]) + dpn[pointer]*dmarker[pointer]*(dwe[k]*( -dex[k]*dun[pointer] -dey[k]*dvn[pointer] )) ;
dfp[PTRQ(i, j, k)] = dfn[pointer]+ rtau * (dfe[pointer] - dfn[pointer]) + runi*dmarker[pointer2]*(dwe[k]*( -dex[k]*dun[pointer2] -dey[k]*dvn[pointer2] )) ;

/* original line */
//dfp[PTRQ(i, j, k)] = dfn[pointer]+ rtau * (dfe[pointer] - dfn[pointer]);
		}
	}
	__syncthreads() ;
}


//********************************‘æ“ñˆ—‹æ‰æ**********************************//
//
//
//******************************************************************************//
__global__ void matrix2(PRECISION *dfp, 
						PRECISION *dpn, PRECISION *dun, PRECISION *dvn,
						PRECISION *ddps, PRECISION *dduv,
						struct Conditions *d_cond)
{
	const int i = blockIdx.x * blockDim.x + threadIdx.x  ;
	const int j = blockIdx.y * blockDim.y + threadIdx.y  ;

	const int imax2 = d_cond->imax2 ;
	const PRECISION rcsu  = d_cond->rcsu  ;

	const int imax = d_cond->imax ;
	const int jmax = d_cond->jmax ;

	if(i==0){
	}
	else if(i==imax){
	}
	else if(j==0){
	}
	else if(j==jmax){
	}
	else{
		PRECISION pss = 0.0 ;
		PRECISION uss = 0.0 ; 
		PRECISION vss = 0.0 ;

		for(int k=0;k<maa;k++){
			pss +=          dfp[PTRQ(i, j, k)] ;
			uss += dex[k] * dfp[PTRQ(i, j, k)] ;
			vss += dey[k] * dfp[PTRQ(i, j, k)] ;
		}
		pss =        pss ;
		uss = rcsu * uss ;
		vss = rcsu * vss ;
		
		ddps[PTR(i, j)] = (pss - dpn[PTR(i, j)]) * (pss - dpn[PTR(i, j)]) ;
		dduv[PTR(i, j)] = (uss - dun[PTR(i, j)]) * (uss - dun[PTR(i, j)]) 
						+ (vss - dvn[PTR(i, j)]) * (vss - dvn[PTR(i, j)]) ;
		dpn[PTR(i, j)] = pss ;
		dun[PTR(i, j)] = uss ;
		dvn[PTR(i, j)] = vss ;
	}
	__syncthreads() ;
}


//********************************‘æŽOˆ—‹æ‰æ**********************************//
//
//
//******************************************************************************//
__global__ void matrix3(PRECISION *dfe, PRECISION *dfp, PRECISION *dfn,
						PRECISION *dpn, PRECISION *dun, PRECISION *dvn,
						struct Conditions *d_cond)
{
	const int i = blockIdx.x * blockDim.x + threadIdx.x  ;
	const int j = blockIdx.y * blockDim.y + threadIdx.y  ;

	const int imax2 = d_cond->imax2 ;
	const PRECISION runi = d_cond->runi ;

	PRECISION euv, qau ;
	
	for(int k=0;k<maa;k++){
		euv = dex[k] * dun[PTR(i,j)] + dey[k] * dvn[PTR(i,j)] ;
		qau = 0.5 * (dun[PTR(i,j)] * dun[PTR(i,j)] + dvn[PTR(i,j)] * dvn[PTR(i,j)]) ;
		dfe[PTRQ(i, j, k)] = dwe[k] * (dpn[PTR(i, j)] + runi * (euv + 1.5 * euv * euv - qau)) ;
		dfn[PTRQ(i, j, k)] = dfp[PTRQ(i, j, k)] ;
	}
	__syncthreads() ;
}


//******************************‘æˆê‹«ŠEˆ—‹æ‰æ********************************//
//
//
//******************************************************************************//
__global__ void CUDAboundp(PRECISION *dpn, PRECISION *dun, PRECISION *dvn,
						   struct Conditions *d_cond)
{
	const int i = blockIdx.x * blockDim.x + threadIdx.x  ;
	const int j = blockIdx.y * blockDim.y + threadIdx.y  ;

	const int isd   = d_cond->isd ;
	const int ied   = d_cond->ied ;
	const int jsd   = d_cond->jsd ;
	const int jed   = d_cond->jed ;
	const int imax  = d_cond->imax ;
	const int jmax  = d_cond->jmax ;
	const int imax2 = d_cond->imax2 ;
	const PRECISION uwui = d_cond->uwui ;

	if(i == 0 && j == 0){//¶‰º//
		int i1 = i + 1 ; int j1 = j + 1 ;
		int i2 = i + 2 ; int j2 = j + 2 ;
		dpn[PTR(i, j)] = (9.0 * dpn[PTR(i1,j1)] - dpn[PTR(i2,j2)]) / 8.0 ;
		dun[PTR(i, j)] = 0.0 ;
		dvn[PTR(i, j)] = 0.0 ;
	}
	
	else if(i == imax && j == 0){//‰E‰º//
		int i1 = i - 1 ; int j1 = j + 1 ;
		int i2 = i - 2 ; int j2 = j + 2 ;
		dpn[PTR(i, j)] = (9.0 * dpn[PTR(i1,j1)] - dpn[PTR(i2,j2)]) / 8.0 ;
		dun[PTR(i, j)] = 0.0 ;
		dvn[PTR(i, j)] = 0.0 ;
	}
	
	else if(i == 0 && j == jmax){//¶ã//
		int i1 = i + 1 ; int j1 = j - 1 ;
		int i2 = i + 2 ; int j2 = j - 2 ;
		dpn[PTR(i, j)] = (9.0 * dpn[PTR(i1,j1)] - dpn[PTR(i2,j2)]) / 8.0 ;
		dun[PTR(i, j)] = 0.5 * uwui ;
		dvn[PTR(i, j)] = 0.0 ;
	}

	else if(i == imax && j == jmax){//‰Eã//
		int i1 = i - 1 ; int j1 = j - 1 ;
		int i2 = i - 2 ; int j2 = j - 2 ;
		dpn[PTR(i, j)] = (9.0 * dpn[PTR(i1,j1)] - dpn[PTR(i2,j2)]) / 8.0 ;
		dun[PTR(i, j)] = 0.5 * uwui ;
		dvn[PTR(i, j)] = 0.0 ;
	}
	
	else if(j == 0 && isd<=i && i<=ied){//‰º•Ó//
		int j1 = j + 1 ; int j2 = j + 2 ;
		dpn[PTR(i, j)] = (9.0 * dpn[PTR(i,j1)] - dpn[PTR(i,j2)]) / 8.0 ;
		dun[PTR(i, j)] = 0.0 ;
		dvn[PTR(i, j)] = 0.0 ;
	}

	else if(j == jmax && isd<=i && i<=ied){//ã•Ó//
		int j1 = j - 1 ; int j2 = j - 2 ;
		dpn[PTR(i, j)] = (9.0 * dpn[PTR(i,j1)] - dpn[PTR(i,j2)]) / 8.0 ;
		dun[PTR(i, j)] = uwui ;
		dvn[PTR(i, j)] = 0.0 ;
	}

	else if(i == 0 && jsd<=j && j<=jed){//¶•Ó//
		int i1 = i + 1 ; int i2 = i + 2 ;
		dpn[PTR(i, j)] = (9.0 * dpn[PTR(i1,j)] - dpn[PTR(i2,j)]) / 8.0 ;
		dun[PTR(i, j)] = 0.0 ;
		dvn[PTR(i, j)] = 0.0 ;
	}
	
	else if(i == imax && jsd<=j && j<=jed){//‰E•Ó//
		int i1 = i - 1 ; int i2 = i - 2 ;
		dpn[PTR(i, j)] = (9.0 * dpn[PTR(i1,j)] - dpn[PTR(i2,j)]) / 8.0 ;
		dun[PTR(i, j)] = 0.0 ;
		dvn[PTR(i, j)] = 0.0 ;
	}

	else{
		//’†S‹æ‰æ//
	}
	__syncthreads() ;
}


//******************************‘æ“ñ‹«ŠEˆ—‹æ‰æ********************************//
//
//
//******************************************************************************//
__global__ void CUDAboundf(PRECISION *dfe, PRECISION *dfp, PRECISION *dfn,
						   PRECISION *dpn, PRECISION *dun, PRECISION *dvn,
						   struct Conditions *d_cond)
{
	const int i = blockIdx.x * blockDim.x + threadIdx.x  ;
	const int j = blockIdx.y * blockDim.y + threadIdx.y  ;

	const int imax2 = d_cond->imax2 ;
	const int isd   = d_cond->isd ;
	const int ied   = d_cond->ied ;
	const int jsd   = d_cond->jsd ;
	const int jed   = d_cond->jed ;
	const int imax  = d_cond->imax ;
	const int jmax  = d_cond->jmax ;
	const PRECISION runi = d_cond->runi ;
	const PRECISION rcsu = d_cond->rcsu ;

	PRECISION pss , uss , vss ;
	PRECISION euv , qau ;

	if(i == 0 && j == 0){//¶‰º//
		pss = 0.0 ;
		uss = 0.0 ;
		vss = 0.0 ;
		for(int k=0;k<maa;k++){
			int ii = i + 1 ;
			int jj = j + 1 ;
			dfp[PTRQ(i,j,k)] = dfe[PTRQ(i,j,k)] + (dfn[PTRQ(ii,jj,k)] - dfe[PTRQ(ii,jj,k)]) ;
			dfp[PTRQ(i,j,k)] = 2.0 * dfp[PTRQ(i,j,k)] - dfp[PTRQ(ii,jj,k)] ;
			pss +=          dfp[PTRQ(i,j,k)] ;
			uss += dex[k] * dfp[PTRQ(i,j,k)] ;
			vss += dey[k] * dfp[PTRQ(i,j,k)] ;
		}
		dpn[PTR(i,j)] =        pss ;
		dun[PTR(i,j)] = rcsu * uss ;
		dvn[PTR(i,j)] = rcsu * vss ;
		for(int k=0;k<maa;k++){
			euv = dex[k] * dun[PTR(i,j)] + dey[k] * dvn[PTR(i,j)] ;
			qau = 0.5 * (dun[PTR(i,j)] * dun[PTR(i,j)] + dvn[PTR(i,j)] * dvn[PTR(i,j)]) ;
			dfe[PTRQ(i,j,k)] = dwe[k] * (dpn[PTR(i,j)] + runi * (euv + 1.5 * euv * euv - qau)) ;
			dfn[PTRQ(i,j,k)] = dfp[PTRQ(i,j,k)] ;
		}
	}
//
	else if(i == imax && j == 0){//‰E‰º//
		pss = 0.0 ;
		uss = 0.0 ;
		vss = 0.0 ;
		for(int k=0;k<maa;k++){
			int ii = i - 1 ;
			int jj = j + 1 ;
			dfp[PTRQ(i,j,k)] = dfe[PTRQ(i,j,k)] + (dfn[PTRQ(ii,jj,k)] - dfe[PTRQ(ii,jj,k)]) ;
			dfp[PTRQ(i,j,k)] = 2.0 * dfp[PTRQ(i,j,k)] - dfp[PTRQ(ii,jj,k)] ;
			pss +=          dfp[PTRQ(i,j,k)] ;
			uss += dex[k] * dfp[PTRQ(i,j,k)] ;
			vss += dey[k] * dfp[PTRQ(i,j,k)] ;
		}
		dpn[PTR(i,j)] =        pss ;
		dun[PTR(i,j)] = rcsu * uss ;
		dvn[PTR(i,j)] = rcsu * vss ;
		for(int k=0;k<maa;k++){
			euv = dex[k] * dun[PTR(i,j)] + dey[k] * dvn[PTR(i,j)] ;
			qau = 0.5 * (dun[PTR(i,j)] * dun[PTR(i,j)] + dvn[PTR(i,j)] * dvn[PTR(i,j)]) ;
			dfe[PTRQ(i,j,k)] = dwe[k] * (dpn[PTR(i,j)] + runi * (euv + 1.5 * euv * euv - qau)) ;
			dfn[PTRQ(i,j,k)] = dfp[PTRQ(i,j,k)] ;
		}
	}
//
	else if(i == 0 && j == jmax){//¶ã//
		pss = 0.0 ;
		uss = 0.0 ;
		vss = 0.0 ;
		for(int k=0;k<maa;k++){
			int ii = i + 1 ;
			int jj = j - 1 ;
			dfp[PTRQ(i,j,k)] = dfe[PTRQ(i,j,k)] + (dfn[PTRQ(ii,jj,k)] - dfe[PTRQ(ii,jj,k)]);
			dfp[PTRQ(i,j,k)] = 2.0 * dfp[PTRQ(i,j,k)] - dfp[PTRQ(ii,jj,k)] ;
			pss +=          dfp[PTRQ(i,j,k)] ;
			uss += dex[k] * dfp[PTRQ(i,j,k)] ;
			vss += dey[k] * dfp[PTRQ(i,j,k)] ;
		}
		dpn[PTR(i,j)] =        pss ;
		dun[PTR(i,j)] = rcsu * uss ;
		dvn[PTR(i,j)] = rcsu * vss ;
		for(int k=0;k<maa;k++){
			euv = dex[k] * dun[PTR(i,j)] + dey[k] * dvn[PTR(i,j)] ;
			qau = 0.5 * (dun[PTR(i,j)] * dun[PTR(i,j)] + dvn[PTR(i,j)] * dvn[PTR(i,j)]) ;
			dfe[PTRQ(i,j,k)] = dwe[k] * (dpn[PTR(i,j)] + runi * (euv + 1.5 * euv * euv - qau)) ;
			dfn[PTRQ(i,j,k)] = dfp[PTRQ(i,j,k)] ;
		}
	}
//
	else if(i == imax && j == jmax){//‰Eã//
		pss = 0.0 ;
		uss = 0.0 ;
		vss = 0.0 ;
		for(int k=0;k<maa;k++){
			int ii = i - 1 ;
			int jj = j - 1 ;
			dfp[PTRQ(i,j,k)] = dfe[PTRQ(i,j,k)] + (dfn[PTRQ(ii,jj,k)] - dfe[PTRQ(ii,jj,k)]) ;
			dfp[PTRQ(i,j,k)] = 2.0 * dfp[PTRQ(i,j,k)] - dfp[PTRQ(ii,jj,k)] ;
			pss +=          dfp[PTRQ(i,j,k)] ;
			uss += dex[k] * dfp[PTRQ(i,j,k)] ;
			vss += dey[k] * dfp[PTRQ(i,j,k)] ;
		}
		dpn[PTR(i,j)] =        pss ;
		dun[PTR(i,j)] = rcsu * uss ;
		dvn[PTR(i,j)] = rcsu * vss ;
		for(int k=0;k<maa;k++){
			euv = dex[k] * dun[PTR(i,j)] + dey[k] * dvn[PTR(i,j)] ;
			qau = 0.5 * (dun[PTR(i,j)] * dun[PTR(i,j)] + dvn[PTR(i,j)] * dvn[PTR(i,j)]) ;
			dfe[PTRQ(i,j,k)] = dwe[k] * (dpn[PTR(i,j)] + runi * (euv + 1.5 * euv * euv - qau)) ;
			dfn[PTRQ(i,j,k)] = dfp[PTRQ(i,j,k)] ;
		}
	}
//
	else if(j == 0 && isd<=i && i<=ied){//‰º•Ó//
		pss = 0.0 ;
		uss = 0.0 ;
		vss = 0.0 ;
		for(int k=0;k<maa;k++){
			int jj = j + 1 ;
			dfp[PTRQ(i,j,k)] = dfe[PTRQ(i,j,k)] + (dfn[PTRQ(i,jj,k)] - dfe[PTRQ(i,jj,k)]) ;
			dfp[PTRQ(i,j,k)] = 2.0 * dfp[PTRQ(i,j,k)] - dfp[PTRQ(i,jj,k)] ;
			pss +=          dfp[PTRQ(i,j,k)] ;
			uss += dex[k] * dfp[PTRQ(i,j,k)] ;
			vss += dey[k] * dfp[PTRQ(i,j,k)] ;
		}
		dpn[PTR(i,j)] =        pss ;
		dun[PTR(i,j)] = rcsu * uss ;
		dvn[PTR(i,j)] = rcsu * vss ;
		for(int k=0;k<maa;k++){
			euv = dex[k] * dun[PTR(i,j)] + dey[k] * dvn[PTR(i,j)] ;
			qau = 0.5 * (dun[PTR(i,j)] * dun[PTR(i,j)] + dvn[PTR(i,j)] * dvn[PTR(i,j)]) ;
			dfe[PTRQ(i,j,k)] = dwe[k] * (dpn[PTR(i,j)] + runi * (euv + 1.5 * euv * euv - qau)) ;
			dfn[PTRQ(i,j,k)] = dfp[PTRQ(i,j,k)] ;
		}
	}
//
	else if(j == jmax && isd<=i && i<=ied){//ã•Ó//
		pss = 0.0 ;
		uss = 0.0 ;
		vss = 0.0 ;
		for(int k=0;k<maa;k++){
			int jj = j - 1 ;
			dfp[PTRQ(i,j,k)] = dfe[PTRQ(i,j,k)] + (dfn[PTRQ(i,jj,k)] - dfe[PTRQ(i,jj,k)]) ;
			dfp[PTRQ(i,j,k)] = 2.0 * dfp[PTRQ(i,j,k)] - dfp[PTRQ(i,jj,k)] ;
			pss +=          dfp[PTRQ(i,j,k)] ;
			uss += dex[k] * dfp[PTRQ(i,j,k)] ;
			vss += dey[k] * dfp[PTRQ(i,j,k)] ;
		}
		dpn[PTR(i,j)] =        pss ;
		dun[PTR(i,j)] = rcsu * uss ;
		dvn[PTR(i,j)] = rcsu * vss ;
		for(int k=0;k<maa;k++){
			euv = dex[k] * dun[PTR(i,j)] + dey[k] * dvn[PTR(i,j)] ;
			qau = 0.5 * (dun[PTR(i,j)] * dun[PTR(i,j)] + dvn[PTR(i,j)] * dvn[PTR(i,j)]) ;
			dfe[PTRQ(i,j,k)] = dwe[k] * (dpn[PTR(i,j)] + runi * (euv + 1.5 * euv * euv - qau)) ;
			dfn[PTRQ(i,j,k)] = dfp[PTRQ(i,j,k)] ;
		}
	}
//
	else if(i == 0 && jsd<=j && j<=jed){//¶•Ó//
		pss = 0.0 ;
		uss = 0.0 ;
		vss = 0.0 ;
		for(int k=0;k<maa;k++){
			int ii = i + 1 ;
			dfp[PTRQ(i,j,k)] = dfe[PTRQ(i,j,k)] + (dfn[PTRQ(ii,j,k)] - dfe[PTRQ(ii,j,k)]) ;
			dfp[PTRQ(i,j,k)] = 2.0 * dfp[PTRQ(i,j,k)] - dfp[PTRQ(ii,j,k)] ;
			pss +=          dfp[PTRQ(i,j,k)] ;
			uss += dex[k] * dfp[PTRQ(i,j,k)] ;
			vss += dey[k] * dfp[PTRQ(i,j,k)] ;
		}
		dpn[PTR(i,j)] =        pss ;
		dun[PTR(i,j)] = rcsu * uss ;
		dvn[PTR(i,j)] = rcsu * vss ;
		for(int k=0;k<maa;k++){
			euv = dex[k] * dun[PTR(i,j)] + dey[k] * dvn[PTR(i,j)] ;
			qau = 0.5 * (dun[PTR(i,j)] * dun[PTR(i,j)] + dvn[PTR(i,j)] * dvn[PTR(i,j)]) ;
			dfe[PTRQ(i,j,k)] = dwe[k] * (dpn[PTR(i,j)] + runi * (euv + 1.5 * euv * euv - qau)) ;
			dfn[PTRQ(i,j,k)] = dfp[PTRQ(i,j,k)] ;
		}
	}
//
	else if(i == imax && jsd<=j && j<=jed){//‰E•Ó//
		pss = 0.0 ;
		uss = 0.0 ;
		vss = 0.0 ;
		for(int k=0;k<maa;k++){
			int ii = i - 1 ;
			dfp[PTRQ(i,j,k)] = dfe[PTRQ(i,j,k)] + (dfn[PTRQ(ii,j,k)] - dfe[PTRQ(ii,j,k)]) ;
			dfp[PTRQ(i,j,k)] = 2.0 * dfp[PTRQ(i,j,k)] - dfp[PTRQ(ii,j,k)] ;
			pss +=          dfp[PTRQ(i,j,k)] ;
			uss += dex[k] * dfp[PTRQ(i,j,k)] ;
			vss += dey[k] * dfp[PTRQ(i,j,k)] ;
		}
		dpn[PTR(i,j)] =        pss ;
		dun[PTR(i,j)] = rcsu * uss ;
		dvn[PTR(i,j)] = rcsu * vss ;
		for(int k=0;k<maa;k++){
			euv = dex[k] * dun[PTR(i,j)] + dey[k] * dvn[PTR(i,j)] ;
			qau = 0.5 * (dun[PTR(i,j)] * dun[PTR(i,j)] + dvn[PTR(i,j)] * dvn[PTR(i,j)]) ;
			dfe[PTRQ(i,j,k)] = dwe[k] * (dpn[PTR(i,j)] + runi * (euv + 1.5 * euv * euv - qau)) ;
			dfn[PTRQ(i,j,k)] = dfp[PTRQ(i,j,k)] ;
		}
	}
//
	else{
		//’†S‹æ‰æ//
	}
	__syncthreads() ;
}


//*****************************ƒƒCƒ“ƒvƒƒOƒ‰ƒ€*********************************//
//
//
//
//******************************************************************************//
int main()    
{

printf("This is my modification\n");

	cudaGetDeviceCount(&numGPUs) ;
	if(numGPUs == 0){
	printf("No GPU detected\n") ;
		return(0) ;
	}

	cudaDeviceProp dev;
	for (int i=0;i<numGPUs;i++){
		cudaGetDeviceProperties(&dev,i) ;
		printf("Device Number :     %d\n", i);
		printf("Using device :      %s\n",       dev.name);
		printf("totalGlobalMem      %d\n",        dev.totalGlobalMem);
		printf("sharedMemPerBlock   %d\n",     dev.sharedMemPerBlock);
		printf("regsPerBlock        %d\n",          dev.regsPerBlock);
		printf("warpSize            %d\n",              dev.warpSize);
		printf("memPitch            %d\n",              dev.memPitch);
		printf("maxThreadsPerBlock  %d\n",    dev.maxThreadsPerBlock);
		printf("maxThreadsDim[3]    %d,%d,%d\n",
				dev.maxThreadsDim[0], dev.maxThreadsDim[1], dev.maxThreadsDim[2]);
		printf("maxGridSize[3]      %d,%d,%d\n",  
				dev.maxGridSize[0], dev.maxGridSize[1], dev.maxGridSize[2]);
		printf("totalConstMem       %d\n",         dev.totalConstMem);
		printf("major.minor         %d.%d\n",        dev.major, dev.minor);
		printf("clockRate           %d\n",             dev.clockRate);
		printf("textureAlignment    %d\n",      dev.textureAlignment);
		printf("deviceOverlap       %d\n",         dev.deviceOverlap);
		printf("multiProcessorCount %d\n",   dev.multiProcessorCount);
	}


	iniset();

	if(imax2%BLOCK_SIZE_X != 0 || jmax2%BLOCK_SIZE_Y != 0){
		// printf("ŠiŽq”‚ÆBLOCK_SIZE‚Æ‚ÌŠÖŒW‚ªm9(^„D^)\n") ;
		printf("Mesh number and BLOCK_SIZE is m9\n") ;
		return(0) ;
	}

	inicon();
	solver();
//	printf("‚±‚±‚ÅI‚í‚è‚¾‚æ");
	printf("Program end\n");
	return(0);
}


//******************************ŒvŽZ—p‰ŠúÝ’è**********************************//
//
//
//******************************************************************************//
void iniset()
{                              
	int i , j ;
	PRECISION rnyi , taui ;
	PRECISION w0 , w1 , w2 ;
	int xsize, ysize ;
	PRECISION chle ; //‘ã•\’·‚³//
//
	xsize = 4 ;
	ysize = 4 ;
	ima = BLOCK_SIZE_X * xsize - 2 ;
	jma = BLOCK_SIZE_Y * ysize - 2 ;
	ncye = 50000 ;
	nwri = 1000 ;
	uwui = 0.1 ;
	reyn = 100 ; //ƒŒƒCƒmƒ‹ƒY”//
//
	isd   = 1 ;
	jsd   = 1 ;
	ied   = ima ;
	jed   = jma ;
	imax  = ima + 1 ;
	jmax  = jma + 1 ;
	imax2 = imax + 1 ; 
	jmax2 = jmax + 1 ;
//
	runi = 1.0 ;
    uini = 0.0 ;
	vini = 0.0 ;
	chle = jmax ;
	rnyi = chle * uwui / reyn ;
	taui = 0.5 * (6.0 * rnyi + 1.0) ;
	rtau = 1.0 / taui ;
	csou = 1.0 / sqrt(3.0) ;
	pini = runi * csou * csou ;
	rcsu = 1.0 / pini ;

printf("Mesh size, imax=%d, jmax=%d \n", imax, jmax);

//
	ex[0] =   0.0 ;
    ex[1] =   1.0 ;
    ex[2] =   0.0 ;
	ex[3] = - 1.0 ;
    ex[4] =   0.0 ;
    ex[5] =   1.0 ;
    ex[6] = - 1.0 ;
    ex[7] = - 1.0 ;
    ex[8] =   1.0 ;
//
	ey[0] =   0.0 ;
    ey[1] =   0.0 ;
    ey[2] =   1.0 ;
    ey[3] =   0.0 ;
    ey[4] = - 1.0 ;
    ey[5] =   1.0 ;
    ey[6] =   1.0 ;
    ey[7] = - 1.0 ;
    ey[8] = - 1.0 ;
//
    w0 = 4.0 /  9.0 ;
    w1 = 1.0 /  9.0 ;
    w2 = 1.0 / 36.0 ;
    we[0] = w0 ;
    we[1] = w1 ;
    we[2] = w1 ;
    we[3] = w1 ;
    we[4] = w1 ;
    we[5] = w2 ;
    we[6] = w2 ;
    we[7] = w2 ;
    we[8] = w2 ;
//
    for(j=0;j<=jmax;j++){
		for(i=0;i<=imax;i++){
			//xg[i][j] = ((i-0.5) / ima) ;
			xg[i][j] = (i / chle) ;
			yg[i][j] = (j / chle) ;
		}
	}
}             


//*******************************‰ŠúðŒÝ’è***********************************//
//
//
//******************************************************************************//
void inicon()
{
        int i , j , k ;
        PRECISION euv , qau ;
//
        for(j=0;j<=jmax;j++){
			for(i=0;i<=imax;i++){
				pn[i][j] = pini ;
				un[i][j] = uini ;
				vn[i][j] = vini ;
			}
		}
//
        j = jmax ;
        for(i=isd;i<=ied;i++){
			pn[i][j] = pini ;
			un[i][j] = uwui ;
			vn[i][j] = vini ;
        }
//
        for(k=0;k<maa;k++){
			for(j=0;j<=jmax;j++){
				for(i=0;i<=imax;i++){
					euv = ex[k] * un[i][j] + ey[k] * vn[i][j] ;
					qau = 0.5 * ( un[i][j] * un[i][j] + vn[i][j] * vn[i][j] ) ;
					fe[i][j][k] = we[k] * (pn[i][j] + runi * (euv + 1.5*euv*euv - qau)) ;
					fn[i][j][k] = fe[i][j][k] ;
					fp[i][j][k] = 0 ;
				}
			}
        }
}

//*****************************ŽåŒvŽZƒvƒƒOƒ‰ƒ€*********************************//
//
//
//******************************************************************************//
void solver()
{
	int i, j, k, nc ;
	struct Conditions cond;
	struct Conditions *d_cond;

	PRECISION dps, duv ;
	PRECISION rda, cda, rnp, dma ;
	FILE *fop ;
//
	PRECISION matrixsize1 = sizeof(PRECISION) * maa ;
	PRECISION matrixsize2 = sizeof(PRECISION) * imax2 * jmax2 ;
	PRECISION matrixsize3 = sizeof(PRECISION) * imax2 * jmax2 * maa ;
//	cda = 1.0e-06 ;
	cda = 1.0e-05 ;
	rnp = 1.0 / ((imax - 1.0)*(jmax - 1.0)) ;
	dma = 0.0 ;

	cond.imax = imax ;
	cond.jmax = jmax ;
	cond.imax2 = imax2 ;
	cond.jmax2 = jmax2 ;
	cond.isd = isd ;
	cond.ied = ied ;
	cond.jsd = jsd ;
	cond.jed = jed ;
	cond.rtau = rtau ;
	cond.rcsu = rcsu ;
	cond.runi = runi ;
	cond.uwui = uwui ;
//
	//ƒzƒXƒg‘¤‚Ì•Ï”Ý’è//
	PRECISION *hex, *hey, *hwe, *hfn, *hfe, *hfp, *hpn, *hun, *hvn, *hdps, *hduv ;

//Marker array
PRECISION *hmarker;

//
	//ƒzƒXƒg‘¤‚Ìƒƒ‚ƒŠŠm•Û//
	hex = (PRECISION*)malloc(matrixsize1);
	hey = (PRECISION*)malloc(matrixsize1);
	hwe = (PRECISION*)malloc(matrixsize1);

	hfn = (PRECISION*)malloc(matrixsize3);
	hfe = (PRECISION*)malloc(matrixsize3);
	hfp = (PRECISION*)malloc(matrixsize3);

	hpn = (PRECISION*)malloc(matrixsize2);
	hun = (PRECISION*)malloc(matrixsize2);
	hvn = (PRECISION*)malloc(matrixsize2);

//allocating array for marker at host
	hmarker = (PRECISION*)malloc(matrixsize2);

	hdps = (PRECISION*)malloc(matrixsize2);
	hduv = (PRECISION*)malloc(matrixsize2);
//
	if (hex == NULL) {
		printf("cannot allocate memory\n");
		return ;
	}
	if (hey == NULL) {
		printf("cannot allocate memory\n");
		return ;
	}
	if (hwe == NULL) {
		printf("cannot allocate memory\n");
		return ;
	}
	if (hfn == NULL) {
		printf("cannot allocate memory\n");
		return ;
	}
	if (hfe == NULL) {
		printf("cannot allocate memory\n");
		return ;
	}
	if (hfp == NULL) {
		printf("cannot allocate memory\n");
		return ;
	}
	if (hpn == NULL) {
		printf("cannot allocate memory\n");
		return ;
	}
	if (hun == NULL) {
		printf("cannot allocate memory\n");
		return ;
	}
	if (hvn == NULL) {
		printf("cannot allocate memory\n");
		return ;
	}
	if (hdps == NULL) {
		printf("cannot allocate memory\n");
		return ;
	}
	if (hduv == NULL) {
		printf("cannot allocate memory\n");
		return ;
	}

int ifront, iback, jdown, jtop;
ifront = 0.25*imax;
iback = 0.75*imax;
jdown = 0.25*jmax;
jtop = 0.75*jmax;

//Add initial value for hmarker
 for(i=0; i<=imax; i++){
	for(j=0; j<=jmax; j++){

if(i >= ifront && i<= iback  ){
	if(j >= jdown && j <= jtop){
	hmarker[PTR(i, j)] = 1.0;
		}
	}
else
	hmarker[PTR(i, j)] = 0.0;
	}
} 

/* FILE *fmarker;
fmarker = fopen("marker.txt","w");
for(i=0; i<= imax; i++){
	for(j=0; j<=jmax; j++){
	fprintf(fmarker, " %d %d %f\n", i, j, hmarker[PTR(i, j)]);
	}
}
fclose(fmarker); */

	//ƒzƒXƒg‚ÉŠe’l‚ð‘}“ü//
	for(j=0;j<=jmax;j++){
		for(i=0;i<=imax;i++){
			for(k=0;k<=8;k++){
				hfn[PTRQ(i, j, k)] = fn[i][j][k] ;
				hfe[PTRQ(i, j, k)] = fe[i][j][k] ;
				hfp[PTRQ(i, j, k)] = fp[i][j][k] ;
			}
		}
	}
	for(j=0;j<=jmax;j++){
		for(i=0;i<=imax;i++){
			hpn[PTR(i, j)] = pn[i][j] ;
			hun[PTR(i, j)] = un[i][j] ;
			hvn[PTR(i, j)] = vn[i][j] ;
			hdps[PTR(i, j)] = 0 ;
			hduv[PTR(i, j)] = 0 ;
		}
	}
	for(k=0;k<=8;k++){
		hex[k] = ex[k] ;
		hey[k] = ey[k] ;
		hwe[k] = we[k] ;
	}
//
	//ƒfƒoƒCƒX‘¤‚Ì•Ï”Ý’è//
	PRECISION *dfn, *dfe, *dfp, *dpn, *dun, *dvn, *ddps, *dduv ;
	PRECISION *dmarker;
//@
	//ƒfƒoƒCƒXƒƒ‚ƒŠŠm•Û‹y‚ÑƒRƒs[//
//	printf("(L¥ƒÖ¥M)\n") ; //cutilSafeCall‚Í‘‚¢‚¿‚áƒ_ƒ//
	printf("Yamazaki nuance, don't get it\n") ; //cutilSafeCall‚Í‘‚¢‚¿‚áƒ_ƒ//
	cudaMemcpyToSymbol(dwe, hwe, matrixsize1) ;
	cudaMemcpyToSymbol(dex, hex, matrixsize1) ;
	cudaMemcpyToSymbol(dey, hey, matrixsize1) ;
//
	cudaMalloc((void**)&dfn, matrixsize3) ;
	cudaMalloc((void**)&dfe, matrixsize3) ;
	cudaMalloc((void**)&dfp, matrixsize3) ;

	cudaMalloc((void**)&dpn, matrixsize2) ;
	cudaMalloc((void**)&dun, matrixsize2) ;
	cudaMalloc((void**)&dvn, matrixsize2) ;	
//allocating device memory for marker array
	cudaMalloc((void**)&dmarker, matrixsize2) ;	

	cudaMalloc((void**)&ddps, matrixsize2) ;
	cudaMalloc((void**)&dduv, matrixsize2) ;	
	cudaMalloc((void **)&d_cond, sizeof(struct Conditions)) ;
//
	cudaMemcpy(dfn, hfn, matrixsize3, cudaMemcpyHostToDevice) ;
	cudaMemcpy(dfe, hfe, matrixsize3, cudaMemcpyHostToDevice) ;
	cudaMemcpy(dfp, hfp, matrixsize3, cudaMemcpyHostToDevice) ;
	cudaMemcpy(dpn, hpn, matrixsize2, cudaMemcpyHostToDevice) ;
	cudaMemcpy(dun, hun, matrixsize2, cudaMemcpyHostToDevice) ;
	cudaMemcpy(dvn, hvn, matrixsize2, cudaMemcpyHostToDevice) ;
//copying memory from host to device for marker array
	cudaMemcpy(dmarker, hmarker, matrixsize2, cudaMemcpyHostToDevice) ;

	cudaMemcpy(ddps, hdps, matrixsize2, cudaMemcpyHostToDevice) ;
	cudaMemcpy(dduv, hduv, matrixsize2, cudaMemcpyHostToDevice) ;
	cudaMemcpy(d_cond, &cond, sizeof(struct Conditions), cudaMemcpyHostToDevice) ;
//
	//ƒuƒƒbƒNƒTƒCƒY‚ÆƒOƒŠƒbƒhƒTƒCƒY‚ÌÝ’è//
	dim3 threads(BLOCK_SIZE_X, BLOCK_SIZE_Y) ;
	dim3 block(imax2/BLOCK_SIZE_X, jmax2/BLOCK_SIZE_Y) ;
//
	//unsigned int timer = 0 ;
	//CUT_SAFE_CALL(cutCreateTimer(&timer)) ;
	//CUT_SAFE_CALL(cutStartTimer(timer)) ;
	PRECISION elapsed_time_ms = 0.0f ;
	cudaEvent_t start, stop;
	cudaEventCreate( &start );
	cudaEventCreate( &stop  );
 
	cudaEventRecord( start, 0 );
//
	for(nc=1;nc<=ncye;nc++){
//original thread
//		matrix1<<<block, threads>>>(dfe, dfp, dfn, d_cond) ;
		matrix1<<<block, threads>>>(dfe, dfp, dfn, dmarker, dun, dvn, dpn, d_cond) ;
		cudaThreadSynchronize() ;

		//‰¼‘z—¬‘©–@ŒÄ‚Ño‚µˆÊ’u//

		matrix2<<<block, threads>>>(dfp, dpn, dun, dvn, ddps, dduv, d_cond) ;
		cudaThreadSynchronize() ;

		cudaMemcpy(hduv, dduv, matrixsize2, cudaMemcpyDeviceToHost) ;
		cudaMemcpy(hdps, ddps, matrixsize2, cudaMemcpyDeviceToHost) ;

		dps = 0.0 ;
		duv = 0.0 ;

		for(j=0;j<=jmax;j++){
			for(i=0;i<=imax;i++){
				dps += hdps[PTR(i, j)] ;
				duv += hduv[PTR(i, j)] ;
			}
		}

		dps = sqrt(      rnp * dps) ;
		duv = sqrt(0.5 * rnp * duv) ;

		if(duv > dma){
			dma = duv ;
		}

		rda = duv / dma ;

		CUDAboundp<<<block, threads>>>(dpn, dun, dvn, d_cond) ;
		cudaThreadSynchronize() ;

		matrix3<<<block, threads>>>(dfe, dfp, dfn, dpn, dun, dvn, d_cond) ;
		cudaThreadSynchronize() ;

		CUDAboundf<<<block, threads>>>(dfe, dfp, dfn, dpn, dun, dvn, d_cond) ;
		cudaThreadSynchronize() ;

		if(nc%nwri == 0){
	
	printf("Resids = %d %e %e %e\n",nc,dps,duv,rda);	
		}

		if(rda < cda) break ; // original definition

	}
//
	printf("Residuals at timestep= %d, L2-pressure= %e, L2-uv= %e, L2-total= %e\n",nc,dps,duv,rda) ;
	//CUT_SAFE_CALL(cutStopTimer(timer)) ;
	//printf("Processing time : %f [msec]\n", cutGetTimerValue(timer)) ;
	//CUT_SAFE_CALL(cutDeleteTimer(timer)) ;

	cudaEventRecord( stop, 0 );
	cudaEventSynchronize( stop );
	cudaEventElapsedTime( &elapsed_time_ms, start, stop );
 
	printf("Processing time : %f [msec]\n", elapsed_time_ms) ;

	cudaEventDestroy( start );
	cudaEventDestroy( stop );

	cudaMemcpy(hpn, dpn, matrixsize2, cudaMemcpyDeviceToHost) ;
	cudaMemcpy(hun, dun, matrixsize2, cudaMemcpyDeviceToHost) ;
	cudaMemcpy(hvn, dvn, matrixsize2, cudaMemcpyDeviceToHost) ;
	cudaThreadSynchronize() ;
	
	printf("Writing result to file\n") ;
	fop = fopen("result.plt", "w") ;
	fprintf(fop, " variables=x,y,p,u,v\n" ); 
    fprintf(fop, "zone t = flowfield\n" ); 
    fprintf(fop, "i = %d , j = %d , f=point \n"  ,ima , jma) ;
	for(j=jsd;j<=jed;j++){
		for(i=isd;i<=ied;i++){
			fprintf (fop, "%e %e %e %e %e\n",xg[i][j], yg[i][j], hpn[PTR(i, j)],
				hun[PTR(i, j)], hvn[PTR(i, j)]) ;
		}
	}
/*
    fprintf(fop, " variables=x,y,p,u,v\n " ); 
    fprintf(fop, "zone t = marker\n" ); 
    fprintf(fop, "i = %d , j = %d , f=point \n"  ,ima , jma) ;	
    for(j=jsd;j<=jed;j++){
		for(i=isd;i<=ied;i++){
			fprintf (fop, "%e %e %e %e %e\n",xg[i][j], yg[i][j], hmarker[PTR(i, j)],
				hun[PTR(i, j)], hvn[PTR(i, j)]) ;
		}
	}
*/
	fclose(fop) ;


//
	free(hex) ;
	free(hey) ;
	free(hwe) ;
	free(hfn) ;
	free(hfe) ;
	free(hfp) ;
	free(hpn) ;
	free(hun) ;
	free(hvn) ;
	free(hdps) ;
	free(hduv) ;
	
	cudaFree(dfe) ;
	cudaFree(dfp) ;
	cudaFree(dfn) ;
	cudaFree(dpn) ;
	cudaFree(dun) ;
	cudaFree(dvn) ;
	cudaFree(ddps) ;
	cudaFree(dduv) ;
}
