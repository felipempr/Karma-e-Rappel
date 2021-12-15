/********************************************************************************************************************************
*															        																
*		SIMULAÇÃO EM DUAS DIMENSÕES DO CAMPO DE FASES PARALELIZADO REPRODUZINDO O MODELO DE JOHNSON DESCRITO NO ARTIGO
*																NÃO HÁ CAMPO DE TEMPERATURA
*																
*********************************************************************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda.h>
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
void cristal(float **cont);

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
	
	//VARIÁVEIS DAS FASES LÍQUIDA, SÓLIDA E IMÓVEL(?)
	float ***PL;
	PL=(float***)malloc(2*sizeof(float**));
	for(int t=0;t<=1;t++){
		PL[t]=(float**)malloc(TELX*sizeof(float*));
		for(int i=0;i<=TELX-1;i++){
			PL[t][i]=(float*)malloc(TELY*sizeof(float));
		}
	}
	float ***PS;
	PS=(float***)malloc(2*sizeof(float**));
	for(int t=0;t<=1;t++){
		PS[t]=(float**)malloc(TELX*sizeof(float*));
		for(int i=0;i<=TELX-1;i++){
			PS[t][i]=(float*)malloc(TELY*sizeof(float));
		}
	}
	float ***PI;
	PI=(float***)malloc(2*sizeof(float**));
	for(int t=0;t<=1;t++){
		PI[t]=(float**)malloc(TELX*sizeof(float*));
		for(int i=0;i<=TELX-1;i++){
			PI[t][i]=(float*)malloc(TELY*sizeof(float));
		}
	}
	float **contorno;
	contorno=(float**)malloc(TELX*sizeof(float*));
		for(int i=0;i<=TELX-1;i++){
			contorno[i]=(float*)malloc(TELY*sizeof(float));
	}
	//VARIÁVEIS REFERENTES À TEMPERATURA (SUPER-RESFRIAMENTO E TEMPERATURA PROVISÓRIA)
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
	
	//VARIÁVEIS DECLARADAS NA GPU
	float *PcL;	
	cudaMalloc((void **) &PcL, 2*TELX*TELY*sizeof(float));
	float *PcS;	
	cudaMalloc((void **) &PcS, 2*TELX*TELY*sizeof(float));
	float *PcI;	
	cudaMalloc((void **) &PcI, 2*TELX*TELY*sizeof(float));
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
	//FILE *arq3;
	
	
printf("\ntelt=%d\n",telt);

arq=fopen("JohnsonPS","w");
arq2=fopen("JohnsonPSbin","wb");
//arq=fopen("E_d08","w");//para binário "wb"
//arq3=fopen("E_d08i","w");//para binário "wb"
//arq2=fopen("DendritasParalE_d08","wb");

//INSERINDO VALOR INICIAL NAS MATRIZES

inicializar(0,1,0,TELX,0,TELY,u,-Ui);
inicializar(0,1,0,TELX,0,TELY,PL,1);
inicializar(0,1,0,TELX,0,TELY,PS,0);
inicializar(0,1,0,TELX,0,TELY,PI,0);
inicializar(0,2,0,TELX,0,TELY,X,-Ui);

for (i=0;i<TELX;i++)//(i=telx/2-r;i<=telx/2+r;i++)
{
	for (j=0;j<(3*TELY/4);j++)//(j=tely/2-r;j<=tely/2+r;j++)
	{
		//if((pow((i),2)+pow((j),2))<=(R*R))//((pow((i-telx/2.0),2)+pow((j-tely/2.0),2))<=(r*r))
		if(i<((j+((3.0*(TELY-1))/4.0))/(3.0*(TELY-1)/(TELX-1)))){
		//if(i<((j+((1.0*(TELY-1))/2.0))/(2.0*(TELY-1)/(TELX-1)))){
		PI[0][i][j]=0;
		PS[0][i][j]=1;
		PL[0][i][j]=0;
		}
		else if(i>((j-((9.0*(TELY-1))/4.0))/-((3.0*(TELY-1))/(TELX-1)))){
		//else if(i>((j-((3.0*(TELY-1))/2.0))/-((2.0*(TELY-1))/(TELX-1)))){
		PI[0][i][j]=0;
		PS[0][i][j]=0;
		PL[0][i][j]=1;
		}
		else{
		PI[0][i][j]=1;
		PS[0][i][j]=0;
		PL[0][i][j]=0;
		}

	}

	for(j=(3*TELY/4);j<TELY;j++)
		if(i<((TELX-1)/2)){
		PI[0][i][j]=0;
		PS[0][i][j]=1;
		PL[0][i][j]=0;
		}
		else{
		PI[0][i][j]=0;
		PS[0][i][j]=0;
		PL[0][i][j]=1;
		}
}
//printf(" PS[1][1]=%f PS[3][6]=%f PS[6][6]=%f\n", PS[0][1][1], PS[0][3][6], PS[0][6][6]);

//COPIANDO VALORES DA CPU (HOST) PARA AS VARIÁVEIS DA GPU (DEVICE)
for(int t=0;t<2;t++)
{
	for(int i=0;i<TELX;i++)
	{
	cudaMemcpy(PcS+((t*TELX*TELY)+(i*TELX)), PS[t][i], TELY*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(PcL+((t*TELX*TELY)+(i*TELX)), PL[t][i], TELY*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(PcI+((t*TELX*TELY)+(i*TELX)), PI[t][i], TELY*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(uc+((t*TELX*TELY)+(i*TELX)), u[t][i], TELY*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(Xc+((t*TELX*TELY)+(i*TELX)), X[t][i], TELY*sizeof(float), cudaMemcpyHostToDevice);
	}
}
/*for (int t=0;t<2;t++){
	for (int i=0;i<4;i++){
	cudaMemcpy(testec+((t*4*4)+(i*4)), teste[t][i], 4*sizeof(float), cudaMemcpyHostToDevice);}
}
for (int t=0;t<2;t++){
	for (int i=0;i<4;i++){
	cudaMemcpy(teste[t][i], testec+((t*4*4)+(i*4)), 4*sizeof(float), cudaMemcpyDeviceToHost);}
}*/

for(tempo=0;tempo<=telt;tempo++)//tempo<=telt
{
//printf("%d, ",tempo);
//CÁLCULO DA VARIÁVEL DE FASE		
P1<<<numBlocks,numThreads>>>(PcS,PcL,PcI,uc,dx,dy,lambda);
P1<<<numBlocks,numThreads>>>(PcL,PcS,PcI,uc,dx,dy,lambda);
P1<<<numBlocks,numThreads>>>(PcI,PcL,PcS,uc,dx,dy,lambda);
cudaDeviceSynchronize();
//CALCULO DO CAMPO DE TEMPERATURAS
//Temp<<<numBlocks, numThreads,4096>>>(uc, Xc, Pc, deltac, Fox, Foy);
//cudaDeviceSynchronize();
/*cudaError_t erro = cudaGetLastError();        // Get error code

   if ( erro != cudaSuccess )
   {
      printf("CUDA Error: %s\n", cudaGetErrorString(erro));
      exit(-1);
   }
*/

//IMPRIMIR 		
	if (tempo%1000==0){//%5000
		for (int t=0;t<2;t++)
		{
			for (int i=0;i<TELX;i++)
			{
			cudaMemcpy(PS[t][i], PcS+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(PL[t][i], PcL+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(PI[t][i], PcI+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(u[t][i], uc+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(X[t][i], Xc+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			}
		}
	//printf(" u0[5][0]=%f u0[30][0]=%f\n", u[0][5][0], u[0][30][0]);
	//printf(" u0[0][5]=%f u0[0][30]=%f\n", u[0][0][5], u[0][0][30]);
	/*fprintf(arq, "phi:t=%f\n",(tempo*dt));
	for(j=0;j<TELY;j++)
	{
		for(i=0;i<TELX;i++)
		{
			fprintf(arq,"%d %d %f\n", i,j,PS[0][i][j]);
		}
	}
	fprintf(arq,"\n\n");*/

	//IMPRESSÃO NO ARQUIVO BINÁRIO
	for (int j = 0; j < TELY; j++) {
		for(int i = 0; i < TELX; i++){
			//contorno[j][i]=PS[0][i][j]+PI[0][i][j]+PL[0][i][j];
			contorno[j][i]=PS[0][i][j]*PS[0][i][j]+PI[0][i][j]*PI[0][i][j]+PL[0][i][j]*PL[0][i][j];
		}
		fwrite(contorno[j], sizeof(float), TELX, arq2);//contorno[j]
		//fwrite(X[1][i], sizeof(float), TELY, arq2);
		}
//printf(" con[5][5]=%f con[70][63]=%f\n", contorno[5][5], contorno[70][63]);
}
//ATUALIZAR VARIÁVEIS PARA O PRÓXIMO CICLO

atualizaTudo<<<numBlocks, numThreads >>>(PcS,PcS);
atualizaTudo<<<numBlocks, numThreads >>>(PcL,PcL);
atualizaTudo<<<numBlocks, numThreads >>>(PcI,PcI);
atualizaTudo<<<numBlocks, numThreads >>>(uc,Xc);
}

printf("lambda=%f\n",lambda);
printf("d0=%f\n",(0.8839*E0/lambda));
printf("E0/d0=%f\n",(D/(0.8839*0.6267)));

cristal(contorno);
//FECHAR E ARQUIVOS E LIBERAR MEMÓRIA
fclose(arq);
fclose(arq2);
//fclose(arq3);
for(int t=0;t<2;t++){
	for(int i=0;i<TELX;i++){
		
		free(PS[t][i]);
		free(PL[t][i]);
		free(PI[t][i]);
		free(u[t][i]);
		free(X[t][i]);
	}
	free(PS[t]);
	free(PL[t]);
	free(PI[t]);
	free(u[t]);
	free(X[t]);
}
free(PS);
free(PL);
free(PI);
free(u);
free(X);

cudaFree(PcS);
cudaFree(PcL);
cudaFree(PcI);
cudaFree(uc);
cudaFree(Xc);
cudaFree(deltac);

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
//FUNÇÃO PARA OBTER VALORES DO CONTORNO
void cristal(float **cont)
{
	for(int i=TELX/2;i<TELX;i++)
	{
		for(int j=TELY-1;j>=0;j--)
		{
			if(0.49<cont[i][j]<0.51)
			printf("x=%d y=%d\n",i,j);
		}
	}
}
