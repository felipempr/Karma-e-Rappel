/********************************************************************************************************************************
*															        																
*		FORMA FINAL PARA SIMULAÇÃO EM DUAS DIMENSÕES DO CAMPO DE FASES COM ANISOTROPIA PARALELIZADO  		
*																
*																
*********************************************************************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda.h>
//#include <omp.h>
#include "entradas.h"
#include "funcC.cuh"

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

void inicializar(int I1,int I2,int J1,int J2,int K1,int K2, float ***Q, float B);

int main(void)
{
	int i;
	int j;
	int tempo;
	int numThreads = 512;
	int numBlocks = (TELX*TELY + numThreads - 1)/numThreads; //3052

	
	//char p;
		
	int telt=Ttot/dt;
	float dx=(float)compL/(float)TELX;
	float dy=(float)compL/(float)TELY;
	//const float sigma=E/sqrt(2*W);
	unsigned char* bmp; 
	bmp=(unsigned char*)malloc(3*TELX*TELY*sizeof(unsigned char));
	
	
	float ***P;
	P=(float***)malloc(2*sizeof(float**));
	for(int t=0;t<=1;t++){
		P[t]=(float**)malloc(TELX*sizeof(float*));
		for(int i=0;i<=TELX-1;i++){
			P[t][i]=(float*)malloc(TELY*sizeof(float));
		}
	}
	bool **FLAG;
	FLAG=(bool**)malloc(TELX*sizeof(bool*));
	for(int t=0;t<TELX;t++){
		FLAG[t]=(bool*)malloc(TELX*sizeof(bool));
	}
	
	float ***u;
	u=(float***)malloc(2*sizeof(float**));
	for(int t=0;t<=1;t++){
		u[t]=(float**)malloc(TELX*sizeof(float*));
		for(int i=0;i<TELX;i++){
			u[t][i]=(float*)malloc(TELY*sizeof(float));
						
		}
	}
	float ***X;
	X=(float***)malloc(2*sizeof(float**));
	for(int t=0;t<=1;t++){
		X[t]=(float**)malloc(TELX*sizeof(float*));
		for(int i=0;i<TELX;i++){
			X[t][i]=(float*)malloc(TELY*sizeof(float));
						
		}
	}
	
	float *Pc;	
	cudaMalloc((void **) &Pc, 2*TELX*TELY*sizeof(float));
	bool *FLAGc;	
	cudaMalloc((void **) &FLAGc, TELX*TELY*sizeof(bool));
	float *uc;	
	cudaMalloc((void **) &uc, 2*TELX*TELY*sizeof(float));
	float *Xc;	
	cudaMalloc((void **) &Xc, 2*TELX*TELY*sizeof(float));
	float *deltac;
	cudaMalloc((void **) &deltac, TELX*TELY*sizeof(float));
	

	//CÁLCULO DOS NÚMEROS DE FOURIER E BIOT
	float Fox=D*dt/(dx*dx);
	float Foy=D*dt/(dy*dy);
	
	float lambda=(1/a2)*((D*TAU0)/(E0*E0));//*(TAU0/(E0*E0));
	printf("lambda=%f\n",lambda);
printf("d0=%f\n",(a1*E0/lambda));
printf("E0/d0=%f\n",E0*(lambda/a1*E0));
	
	FILE *arq;
	FILE *arq2;
	FILE *arq3;
	
	
printf("\ntelt=%d\n",telt);

arq=fopen("Wd0_9_posy","w");//para binário "wb"
arq3=fopen("Wd0_9_posx","w");//para binário "wb"
arq2=fopen("Karma","wb");

//INSERINDO VALOR INICIAL NAS MATRIZES

inicializar(0,1,0,TELX,0,TELY,u,-Ui);
inicializar(0,1,0,TELX,0,TELY,P,LIQUIDO);
inicializar(0,2,0,TELX,0,TELY,X,-Ui);
/*
for (i=0;i<=R;i++)//(i=telx/2-r;i<=telx/2+r;i++)
{
	for (j=0;j<R;j++)//(j=tely/2-r;j<=tely/2+r;j++)
	{
		if((pow((i),2.0)+pow((j),2.0))<=(R*R))//((pow((i-telx/2.0),2)+pow((j-tely/2.0),2))<=(r*r))
		P[0][i][j]=SOLIDO;
	}
}
*/
readBMP(bmp);
cond_inicial(bmp,P,FLAG);

for(int t=0;t<2;t++)
{
	for(int i=0;i<TELX;i++)
	{
	cudaMemcpy(Pc+((t*TELX*TELY)+(i*TELX)), P[t][i], TELY*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(uc+((t*TELX*TELY)+(i*TELX)), u[t][i], TELY*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(Xc+((t*TELX*TELY)+(i*TELX)), X[t][i], TELY*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(FLAGc+(i*TELX), FLAG[i], TELY*sizeof(bool), cudaMemcpyHostToDevice);
	}
}


for(tempo=0;tempo<=telt;tempo++)//tempo<=telt
{
printf(".");
//CÁLCULO DA VARIÁVEL DE FASE		
P1<<<numBlocks,numThreads>>>(FLAGc,Pc,uc,dx,dy,lambda);
//cudaDeviceSynchronize();

//CALCULO DO CAMPO DE TEMPERATURAS
Temp<<<numBlocks, numThreads,4096>>>(uc, Xc, Pc, deltac, Fox, Foy);
cudaDeviceSynchronize();
/*cudaError_t erro = cudaGetLastError();        // Get error code

   if ( erro != cudaSuccess )
   {
      printf("CUDA Error: %s\n", cudaGetErrorString(erro));
      exit(-1);
   }
*/

//IMPRIMIR 		
	if (tempo%1000==0){
		printf("%d, ",tempo);
		for (int t=0;t<2;t++)
		{
			for (int i=0;i<TELX;i++)
			{
			cudaMemcpy(P[t][i], Pc+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			//cudaMemcpy(PS[t][i], PcS+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			//cudaMemcpy(PL[t][i], PcL+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			//cudaMemcpy(PI[t][i], PcI+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(u[t][i], uc+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(X[t][i], Xc+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(FLAG[i], FLAGc+(i*TELX), TELY*sizeof(bool), cudaMemcpyDeviceToHost);
			}
			//printf(" P0[5][5]=%f P1[30][30]=%f\n", P[t][5][5], P[t][30][30]);
			//printf(" X0[5][5]=%f X1[30][30]=%f\n", X[t][5][5], X[t][30][30]);
		}
	printf(" u0[5][0]=%f u0[30][0]=%f\n", u[0][5][0], u[0][30][0]);
	printf(" u0[0][5]=%f u0[0][30]=%f\n", u[0][0][5], u[0][0][30]);

	if (P[0][100][100]!=P[0][100][100]){
		break;
		}

/*	//IMPRESSÃO DAS POSIÇÕES E TEMPO DA PONTA DA DENDRITA NAS DIREÇÕES X E Y
		fprintf(arq,"phi:t=%.8f\n",(tempo*dt));
	for(int j=0;j<TELY-1;j++){
		if (P[1][0][j]*P[1][0][j+1]<0)
		fprintf(arq,"%d %f %d %f\n",j, P[1][0][j],j+1,P[1][0][j+1]);
		}
		fprintf(arq,"\n\n");

		fprintf(arq3,"phi:t=%.8f\n",(tempo*dt));
	for(int i=0;i<TELY-1;i++){
		if (P[1][i][0]*P[1][i+1][0]<0)
		fprintf(arq3,"%d %f %d %f\n",i, P[1][i][0],i+1,P[1][i+1][0]);
		}
		fprintf(arq3,"\n\n");
*/		
//IMPRESSÃO DO CAMPO DE FASES
	for (int i = 0; i < TELX; i++) {
		fwrite(P[0][i], sizeof(float), TELY, arq2);
		//fwrite(X[1][i], sizeof(float), TELY, arq2);
		}
}
//ATUALIZAR VARIÁVEIS PARA O PRÓXIMO CICLO

atualizaTudo<<<numBlocks, numThreads >>>(Pc,Pc);
atualizaTudo<<<numBlocks, numThreads >>>(uc,Xc);
}

printf("lambda=%f\n",lambda);
printf("d0=%f\n",(0.8839*E0/lambda));
printf("E0/d0=%f\n",(D/(0.8839*0.6267)));
//FECHAR E ARQUIVOS E LIBERAR MEMÓRIA

fclose(arq);
fclose(arq2);
fclose(arq3);

for(int t=0;t<2;t++){
	for(int i=0;i<TELX;i++){
		free(P[t][i]);
		//free(PS[t][i]);
		//free(PL[t][i]);
		//free(PI[t][i]);
		free(u[t][i]);
		free(X[t][i]);
		
	}
	free(P[t]);
	//free(PS[t]);
	//free(PL[t]);
	//free(PI[t]);
	free(u[t]);
	free(X[t]);
	
}
for(int i=0;i<TELX;i++){
free(FLAG[i]);
}
free(P);
//free(PS);
//free(PL);
//free(PI);
free(u);
free(X);
free(FLAG);
cudaFree(Pc);
//cudaFree(PcS);
//cudaFree(PcL);
//cudaFree(PcI);
cudaFree(uc);
cudaFree(Xc);
cudaFree(deltac);
cudaFree(FLAGc);

return 0;	
}

//FUNÇÃO QUE INICIALIZA AS MATRIZES
void inicializar(int I1,int I2,int J1,int J2,int K1,int K2, float ***Q, float B)
{
	for(int i=I1;i<=I2-1;i++)
	{	
		for(int j=J1;j<=J2-1;j++)
		{
			for(int k=K1;k<=K2-1;k++)
			{
				Q[i][j][k]=B;
			}
		}
	}
}

