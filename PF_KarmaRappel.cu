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

//void inicializar(int I1,int I2,int J1,int J2,int K1,int K2, float ***Q, float B);

int main(void)
{
	int i;
	int j;
	int tempo;
	int numThreads = 512;
	int numBlocks = (TELX*TELY + numThreads - 1)/numThreads; //3052
		
	int telt=Ttot/dt;
	float dx=(float)compL/(float)TELX;
	float dy=(float)compL/(float)TELY;

	unsigned char* bmp; 
	bmp=(unsigned char*)malloc(3*TELX*TELY*sizeof(unsigned char));
		
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

	float *rhs;
	rhs=(float*)malloc(sizeof(float));
	float *lamb;
	lamb=(float*)malloc(sizeof(float));
	/*float ***u;
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
	}*/

	float **contorno;
	contorno=(float**)malloc(TELX*sizeof(float*));
		for(int i=0;i<=TELX-1;i++){
			contorno[i]=(float*)malloc(TELY*sizeof(float));
	}
	
	float *PcL;	
	cudaMalloc((void **) &PcL, 2*TELX*TELY*sizeof(float));
	float *PcS;	
	cudaMalloc((void **) &PcS, 2*TELX*TELY*sizeof(float));
	float *PcI;	
	cudaMalloc((void **) &PcI, 2*TELX*TELY*sizeof(float));
	
	float *crhs;	
	cudaMalloc((void **) &crhs, sizeof(float));
	float *clamb;	
	cudaMalloc((void **) &clamb, sizeof(float));
	/*
	float *uc;	
	cudaMalloc((void **) &uc, 2*TELX*TELY*sizeof(float));
	float *Xc;	
	cudaMalloc((void **) &Xc, 2*TELX*TELY*sizeof(float));
	float *deltac;
	cudaMalloc((void **) &deltac, TELX*TELY*sizeof(float));*/
	

	//CÁLCULO E IMPRESSÃO DE VARIÁVEIS DE SUPORTE
	//float Fox=D*dt/(dx*dx);
	//float Foy=D*dt/(dy*dy);
	//float lambda=(1/a2)*((D*TAU0)/(E0*E0));//*(TAU0/(E0*E0));
	//printf("lambda=%f\n",lambda);
	//printf("d0=%f\n",(a1*E0/lambda));
	//printf("E0/d0=%f\n",E0*(lambda/a1*E0));
	printf("\ntelt=%d\n",telt);
		
	//ABERTURA DE ARQUIVOS EXTERNOS
	FILE *arq;
	FILE *arq2;
	arq=fopen("Nestler1","w");//para binário "wb"
	arq2=fopen("Nestler2","wb");//para binário "wb"
	
//INSERINDO VALOR INICIAL NAS MATRIZES
//inicializar(0,1,0,TELX,0,TELY,u,-Ui);
//inicializar(0,2,0,TELX,0,TELY,X,-Ui);
inicializar(0,1,0,TELX,0,TELY,PL,N_EXISTE);
inicializar(0,1,0,TELX,0,TELY,PS,N_EXISTE);
inicializar(0,1,0,TELX,0,TELY,PI,N_EXISTE);

readBMP(bmp);
cond_inicial(bmp,PI,PL,PS);
/*for (i=0;i<=R;i++)//(i=telx/2-r;i<=telx/2+r;i++)
{
	for (j=0;j<R;j++)//(j=tely/2-r;j<=tely/2+r;j++)
	{
		if((pow((i),2.0)+pow((j),2.0))<=(R*R))//((pow((i-telx/2.0),2)+pow((j-tely/2.0),2))<=(r*r))
		PS[0][i][j]=SOLIDO;
	}
}
*/
//COPIANDO VARIÁVEIS PARA A GPU
for(int t=0;t<2;t++)
{
	for(int i=0;i<TELX;i++)
	{
	cudaMemcpy(PcS+((t*TELX*TELY)+(i*TELX)), PS[t][i], TELY*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(PcL+((t*TELX*TELY)+(i*TELX)), PL[t][i], TELY*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(PcI+((t*TELX*TELY)+(i*TELX)), PI[t][i], TELY*sizeof(float), cudaMemcpyHostToDevice);
	//cudaMemcpy(uc+((t*TELX*TELY)+(i*TELX)), u[t][i], TELY*sizeof(float), cudaMemcpyHostToDevice);
	//cudaMemcpy(Xc+((t*TELX*TELY)+(i*TELX)), X[t][i], TELY*sizeof(float), cudaMemcpyHostToDevice);
	}
}
cudaMemcpy(crhs, rhs, sizeof(float), cudaMemcpyDeviceToHost);
cudaMemcpy(clamb, lamb, sizeof(float), cudaMemcpyDeviceToHost);
//LOOP PRINCIPAL
for(tempo=0;tempo<=telt;tempo++)
{
	printf("%d, ",tempo);

	//CÁLCULO DA VARIÁVEL DE FASE		
	P1<<<numBlocks,numThreads>>>(crhs,clamb,PcL,PcI,PcS,dx,dy,GAMMAAC,GAMMAAB,GAMMABC);
	P1<<<numBlocks,numThreads>>>(crhs,clamb,PcS,PcL,PcI,dx,dy,GAMMAAB,GAMMAAC,GAMMABC);
	P1<<<numBlocks,numThreads>>>(crhs,clamb,PcI,PcL,PcS,dx,dy,GAMMABC,GAMMAAC,GAMMAAB);
/*
	for (int t=0;t<2;t++)
		{
			for (int i=0;i<TELX;i++)
			{
			cudaMemcpy(PS[t][i], PcS+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(PL[t][i], PcL+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(PI[t][i], PcI+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			
			//cudaMemcpy(u[t][i], uc+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			//cudaMemcpy(X[t][i], Xc+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			}
		}
		cudaMemcpy(rhs, crhs, sizeof(float), cudaMemcpyDeviceToHost);
		cudaMemcpy(lamb, clamb, sizeof(float), cudaMemcpyDeviceToHost);
	printf("PLrhs=%f",*rhs);
	printf(" PLlamb=%f\n",*lamb);

	P1<<<numBlocks,numThreads>>>(crhs,clamb,PcS,PcL,PcI,dx,dy,GAMMAAB,GAMMAAC,GAMMABC,GAMMAABC);
	for (int t=0;t<2;t++)
		{
			for (int i=0;i<TELX;i++)
			{
			cudaMemcpy(PS[t][i], PcS+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(PL[t][i], PcL+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(PI[t][i], PcI+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			
			//cudaMemcpy(u[t][i], uc+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			//cudaMemcpy(X[t][i], Xc+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			}
		}
		cudaMemcpy(rhs, crhs, sizeof(float), cudaMemcpyDeviceToHost);
		cudaMemcpy(lamb, clamb, sizeof(float), cudaMemcpyDeviceToHost);
	printf("PSrhs=%f",*rhs);
	printf(" PSlamb=%f\n",*lamb);

	P1<<<numBlocks,numThreads>>>(crhs,clamb,PcI,PcL,PcS,dx,dy,GAMMABC,GAMMAAC,GAMMAAB,GAMMAABC);
	for (int t=0;t<2;t++)
		{
			for (int i=0;i<TELX;i++)
			{
			cudaMemcpy(PS[t][i], PcS+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(PL[t][i], PcL+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(PI[t][i], PcI+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			
			//cudaMemcpy(u[t][i], uc+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			//cudaMemcpy(X[t][i], Xc+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			}
		}
		cudaMemcpy(rhs, crhs, sizeof(float), cudaMemcpyDeviceToHost);
		cudaMemcpy(lamb, clamb, sizeof(float), cudaMemcpyDeviceToHost);
	printf("PIrhs=%f",*rhs);
	printf(" PIlamb=%f\n\n",*lamb);
	
	
	*/			
	cudaDeviceSynchronize();

	//CALCULO DO CAMPO DE TEMPERATURAS
	//Temp<<<numBlocks, numThreads,4096>>>(uc, Xc, Pc, deltac, Fox, Foy);
	//cudaDeviceSynchronize();
	
	//FUNÇÃO DE CAPTAÇÃO DE ERROS 
	/*cudaError_t erro = cudaGetLastError();        // Get error code
		if ( erro != cudaSuccess )
		{
			printf("CUDA Error: %s\n", cudaGetErrorString(erro));
			exit(-1);
		}
	*/

	//IMPRIMIR EM DETERMINADOS INTERVALOS DE TEMPO 	
	if (tempo%4000==0)
	{
		
		//COPIANDO VARIÁVEIS DA GPU
		for (int t=0;t<2;t++)
		{
			for (int i=0;i<TELX;i++)
			{
			cudaMemcpy(PI[t][i], PcI+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(PL[t][i], PcL+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(PS[t][i], PcS+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			
			//cudaMemcpy(u[t][i], uc+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			//cudaMemcpy(X[t][i], Xc+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			}
		}
		cudaMemcpy(rhs, crhs, sizeof(float), cudaMemcpyDeviceToHost);
		cudaMemcpy(lamb, clamb, sizeof(float), cudaMemcpyDeviceToHost);
		printf("rhs=%f",*rhs);
		printf(" lamb=%f",*lamb);
	fprintf(arq, "phi:t=%f\n",(tempo*dt));
	for(j=0;j<TELY;j++)
	{
		for(i=0;i<TELX;i++)
		{
			//contorno[j][i]=PS[0][i][j]+PI[0][i][j]+PL[0][i][j];
			contorno[j][i]=PL[0][i][j]*PL[0][i][j]+PS[0][i][j]*PS[0][i][j]+PI[0][i][j]*PI[0][i][j];//contorno[-j+TELY-1][i]
			fprintf(arq,"%d %d %f\n", i,j,contorno[j][i]); //ARQUIVO TEXTO
		}
		fwrite(contorno[j], sizeof(float), TELX, arq2);//ARQUIVO BINÁRIO
	}
	//IMPRESSÃO DA TEMPERATURA EM DETERMINADOS PONTOS PARA OBSERVAÇÃO DURANTE A EXECUÇÃO
		
		printf("\n PS[500][500]=%f PS[500][998]=%f PS[2][500]=%f\n",  PS[0][500][500], PS[0][500][998],PS[0][2][500]);
		printf(" PI[500][500]=%f PI[500][998]=%f PI[2][500]=%f\n", PI[0][500][500], PI[0][500][998],PI[0][2][500]);
		printf(" PL[500][500]=%f PL[500][998]=%f PL[2][500]=%f\n", PL[0][500][500], PL[0][500][998],PL[0][2][500]);
		//A leitura é feita de baixo pra cima. A coordenada x é ok, mas a y é invertida.
		printf(" contorno=%f\n", contorno[500][500]);
		//printf(" u0[0][5]=%f u0[0][30]=%f\n", u[0][0][5], u[0][0][30]);
		

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
	}

//ATUALIZAR VARIÁVEIS PARA O PRÓXIMO CICLO
atualizaTudo<<<numBlocks, numThreads >>>(PcS,PcS);
atualizaTudo<<<numBlocks, numThreads >>>(PcL,PcL);
atualizaTudo<<<numBlocks, numThreads >>>(PcI,PcI);
//atualizaTudo<<<numBlocks, numThreads >>>(uc,Xc);
}//FIM DO LOOP PRINCIPAL

//IMPRESSÃO DAS VARIÁVEIS DE SUPORTE DO FIM DA EXECUÇÃO
//printf("lambda=%f\n",lambda);
//printf("d0=%f\n",(0.8839*E0/lambda));
//printf("E0/d0=%f\n",(D/(0.8839*0.6267)));

//FECHAR E ARQUIVOS E LIBERAR MEMÓRIA
fclose(arq);
fclose(arq2);

for(int t=0;t<2;t++)
{
	for(int i=0;i<TELX;i++)
	{
		free(PS[t][i]);
		free(PL[t][i]);
		free(PI[t][i]);
		//free(u[t][i]);
		//free(X[t][i]);
	}
	free(PS[t]);
	free(PL[t]);
	free(PI[t]);
	//free(u[t]);
	//free(X[t]);
}
free(PS);
free(PL);
free(PI);
//free(u);
//free(X);

cudaFree(PcS);
cudaFree(PcL);
cudaFree(PcI);
//cudaFree(uc);
//cudaFree(Xc);
//cudaFree(deltac);

return 0;	
}

//FUNÇÃO QUE INICIALIZA AS MATRIZES
/*void inicializar(int I1,int I2,int J1,int J2,int K1,int K2, float ***Q, float B)
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
}*/

