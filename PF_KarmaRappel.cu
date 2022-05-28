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
	bool FLAG;
	int i;
	int j;
	int tempo;
	int numThreads = 512;
	int numBlocks = (TELX*TELY + numThreads - 1)/numThreads; //3052
		
	int telt=Ttot/dt;
	float dx=(float)compL/(float)TELX;
	float dy=(float)compL/(float)TELY;
	float ySim[TELX];

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

	/*float *rhs;
	rhs=(float*)malloc(sizeof(float));
	float *lamb;
	lamb=(float*)malloc(sizeof(float));
	*/
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
	
	/*
	float *crhs;	
	cudaMalloc((void **) &crhs, sizeof(float));
	float *clamb;	
	cudaMalloc((void **) &clamb, sizeof(float));
	*/
	float *uc;	
	cudaMalloc((void **) &uc, 2*TELX*TELY*sizeof(float));
	float *Xc;	
	cudaMalloc((void **) &Xc, 2*TELX*TELY*sizeof(float));
	float *deltac;
	cudaMalloc((void **) &deltac, TELX*TELY*sizeof(float));
	

	//CÁLCULO E IMPRESSÃO DE VARIÁVEIS DE SUPORTE
	float Fox=(constK/cv)*dt/(dx*dx);
	float Foy=(constK/cv)*dt/(dy*dy);
	//float lambda=(1/a2)*((D*TAU0)/(E0*E0));//*(TAU0/(E0*E0));
	//printf("lambda=%f\n",lambda);
	//printf("d0=%f\n",(a1*E0/lambda));
	//printf("E0/d0=%f\n",E0*(lambda/a1*E0));
	printf("\ntelt=%d\n",telt);
		
	//ABERTURA DE ARQUIVOS EXTERNOS
	FILE *arq;
	FILE *arqb;
	FILE *arqc;
	FILE *arq2;
	FILE *arq3;
	arq=fopen("VInt","w");//para binário "wb"
	arqb=fopen("TInt","w");//para binário "wb"
	arq2=fopen("NestlerT","wb");//para binário "wb"
	arq3=fopen("NestlerTemp","wb");
	
//INSERINDO VALOR INICIAL NAS MATRIZES
inicializar(0,2,0,TELX,0,TELY,u,TmS);
inicializar(0,2,0,TELX,0,TELY,X,TmS);
inicializar(0,1,0,TELX,0,TELY,PL,EXISTE);
inicializar(0,1,0,TELX,0,TELY,PS,N_EXISTE);
inicializar(0,1,0,TELX,0,TELY,PI,N_EXISTE);

//for(j=0;j<=10;j++){
for (i=0;i<TELX;i++)//(i=telx/2-r;i<=telx/2+r;i++)
{
	PS[0][i][0]=EXISTE;
	PL[0][i][0]=N_EXISTE;
	u[0][i][0]=T;
	u[1][i][0]=T;
	X[0][i][0]=T;
	X[1][i][0]=T;

}
//}

	/*for (j=0;j<(7*TELY/8);j++)//(j=tely/2-r;j<=tely/2+r;j++)
	{
		
		//if((pow((i),2)+pow((j),2))<=(R*R))//((pow((i-telx/2.0),2)+pow((j-tely/2.0),2))<=(r*r))
		//if(i<((j+((3.0*(TELY-1))/4.0))/(3.0*(TELY-1)/(TELX-1)))){
		if(i<((j+((7.0*(TELY-1))/8.0))/(7.0*(TELY-1)/(2.0*(TELX-1))))){
		PI[0][i][j]=0;
		PS[0][i][j]=1;
		PL[0][i][j]=0;
		}
		//else if(i>((j-((9.0*(TELY-1))/4.0))/-((3.0*(TELY-1))/(TELX-1)))){
		else if(i>((j-((21.0*(TELY-1))/8.0))/-((7.0*(TELY-1))/(2.0*(TELX-1))))){
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

	for(j=(7*TELY/8);j<TELY;j++)
		if(i<((TELX-1)/2)){
		PI[0][i][j]=0;
		PS[0][i][j]=1;
		PL[0][i][j]=0;
		}
		else{
		PI[0][i][j]=0;
		PS[0][i][j]=0;
		PL[0][i][j]=1;
		}*/

//readBMP(bmp);
//cond_inicial(bmp,PL,PI,PS);
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
	cudaMemcpy(uc+((t*TELX*TELY)+(i*TELX)), u[t][i], TELY*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(Xc+((t*TELX*TELY)+(i*TELX)), X[t][i], TELY*sizeof(float), cudaMemcpyHostToDevice);
	}
}
//cudaMemcpy(crhs, rhs, sizeof(float), cudaMemcpyDeviceToHost);
//cudaMemcpy(clamb, lamb, sizeof(float), cudaMemcpyDeviceToHost);
//LOOP PRINCIPAL
for(tempo=0;tempo<=telt;tempo++)
{
	printf("%d, ",tempo);

	//CÁLCULO DA VARIÁVEL DE FASE		
	//acrescentar antes crhs e clamb para monitorar esses valores
	P1<<<numBlocks,numThreads>>>(PcS,PcL,PcI,uc,dx,dy,GAMMAAB,GAMMAAC,GAMMABC,LS,TmS,LL,TmL,LI,TmI,0,2,1);
	P1<<<numBlocks,numThreads>>>(PcI,PcL,PcS,uc,dx,dy,GAMMABC,GAMMAAC,GAMMAAB,LI,TmI,LL,TmL,LS,TmS,1,2,0);
	P1<<<numBlocks,numThreads>>>(PcL,PcI,PcS,uc,dx,dy,GAMMAAC,GAMMAAB,GAMMABC,LL,TmL,LI,TmI,LS,TmS,2,1,0);
		cudaDeviceSynchronize();

	//CALCULO DO CAMPO DE TEMPERATURAS
	Temp<<<numBlocks, numThreads,4096>>>(uc, Xc, PcS,PcI, deltac, Fox, Foy,LS,LI);
	cudaDeviceSynchronize();
	
	//FUNÇÃO DE CAPTAÇÃO DE ERROS 
	/*cudaError_t erro = cudaGetLastError();        // Get error code
		if ( erro != cudaSuccess )
		{
			printf("CUDA Error: %s\n", cudaGetErrorString(erro));
			exit(-1);
		}
	*/

	//IMPRIMIR EM DETERMINADOS INTERVALOS DE TEMPO 	
	if (tempo%480==0)
	{
		
		FLAG=1;//flaf pra pegar o valor do contorno e calcular a velocidade da interface
		//COPIANDO VARIÁVEIS DA GPU
		for (int t=0;t<2;t++)
		{
			for (int i=0;i<TELX;i++)
			{
			cudaMemcpy(PI[t][i], PcI+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(PL[t][i], PcL+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(PS[t][i], PcS+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(u[t][i], uc+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(X[t][i], Xc+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			}
		}
		//fprintf(arq, "phi:t=%f\n",(tempo*dt));
	for(j=0;j<TELY;j++)
	{
		for(i=0;i<TELX;i++)
		{
			//contorno[j][i]=PS[0][i][j]+PI[0][i][j]+PL[0][i][j];
			contorno[j][i]=PL[0][i][j]*PL[0][i][j]+PS[0][i][j]*PS[0][i][j]+PI[0][i][j]*PI[0][i][j];//contorno[-j+TELY-1][i]
			if(PS[0][j][i]<0.4&&j==(TELX/2)&&FLAG==1){
			fprintf(arq,"%f %d \n", tempo*dt, i-1); //ARQUIVO TEXTO
			FLAG=0;}
		}
		if(tempo==38400){
			fprintf(arqb,"%d %f \n", j, u[0][TELX/2][j]); //ARQUIVO TEXTO
			}
		//fprintf(arq,"\n\n");
		fwrite(PS[0][j], sizeof(float), TELX, arq2);//ARQUIVO BINÁRIO
		fwrite(u[0][j], sizeof(float), TELX, arq3);//ARQUIVO BINÁRIO
	}
	
	//IMPRESSÃO DA TEMPERATURA EM DETERMINADOS PONTOS PARA OBSERVAÇÃO DURANTE A EXECUÇÃO
		printf("\n PS[10][10]=%f PS[190][10]=%f PS[100][100]=%f\n",  PS[0][10][10], PS[0][190][10],PS[0][100][100]);
		printf("\n PL[10][10]=%f PL[190][10]=%f PL[100][100]=%f\n",  PL[0][10][10], PL[0][190][10],PL[0][100][100]);
		printf("\n PI[10][10]=%f PI[190][10]=%f PI[100][100]=%f\n",  PI[0][10][10], PI[0][190][10],PI[0][100][100]);
		//A leitura é feita de baixo pra cima. A coordenada x é ok, mas a y é invertida.
		//printf(" contorno=%f\n", contorno[500][500]);
		printf(" u0[0][5]=%f u0[100][100]=%f\n", u[0][0][5], u[0][100][100]);
		

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
/*if (tempo==37632)//tem que ser um múltiplo do valor do loop acima
	{
	cristal(contorno,ySim);//(arqb, contorno, ySim);
	}*/
//ATUALIZAR VARIÁVEIS PARA O PRÓXIMO CICLO
atualizaTudo<<<numBlocks, numThreads >>>(PcS,PcS);
atualizaTudo<<<numBlocks, numThreads >>>(PcL,PcL);
atualizaTudo<<<numBlocks, numThreads >>>(PcI,PcI);
atualizaTudo<<<numBlocks, numThreads >>>(uc,Xc);


}//FIM DO LOOP PRINCIPAL
//curvas(arq,arqb,ySim);
//IMPRESSÃO DAS VARIÁVEIS DE SUPORTE DO FIM DA EXECUÇÃO
//printf("lambda=%f\n",lambda);
//printf("d0=%f\n",(0.8839*E0/lambda));
//printf("E0/d0=%f\n",(D/(0.8839*0.6267)));

//FECHAR E ARQUIVOS E LIBERAR MEMÓRIA
fclose(arq);
fclose(arqb);
fclose(arq2);
fclose(arq3);

for(int t=0;t<2;t++)
{
	for(int i=0;i<TELX;i++)
	{
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
