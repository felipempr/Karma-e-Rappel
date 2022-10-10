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

#ifdef NAN
/* NAN is supported */
#endif
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
	int numThreads = 512;//512
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
	FILE *arqS;
	FILE *arqL;
	FILE *arqI;
	FILE *arqT;

	FILE *curva1;
	FILE *curva2;
	FILE *curva3;
	FILE *curvaT1;
	FILE *curvaT2;
	FILE *curvaT3;

	arq=fopen("VInt","w");//para binário "wb"
	arqb=fopen("TInt","w");//para binário "wb"
	arqS=fopen("NestlerS","wb");//para binário "wb"
	arqL=fopen("NestlerL","wb");//para binário "wb"
	arqI=fopen("NestlerI","wb");//para binário "wb"
	arqT=fopen("NestlerTemp","wb");

	curva1=fopen("curva1","w");//para binário "wb"
	curva2=fopen("curva2","w");//para binário "wb"
	curva3=fopen("curva3","w");//para binário "wb"
	curvaT1=fopen("curvaT1","w");//para binário "wb"
	curvaT2=fopen("curvaT2","w");//para binário "wb"
	curvaT3=fopen("curvaT3","w");//para binário "wb"
	
//INSERINDO VALOR INICIAL NAS MATRIZES
inicializar(0,2,0,TELX,0,TELY,u,TmS);
inicializar(0,2,0,TELX,0,TELY,X,TmS);
inicializar(0,1,0,TELX,0,TELY,PL,N_EXISTE);
inicializar(0,1,0,TELX,0,TELY,PS,N_EXISTE);
inicializar(0,1,0,TELX,0,TELY,PI,0.0);


/*for (i=0;i<TELX;i++)//(i=telx/2-r;i<=telx/2+r;i++)
{
	PS[0][i][0]=0.5;
	PL[0][i][0]=0.5;
}*/


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

readBMP(bmp);
cond_inicial(bmp,PL,PI,PS);
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
printf("\n%f %f, ",dx,dy);
for(tempo=0;tempo<=telt;tempo++)
{
	printf("%d, ",tempo);
	
	//CÁLCULO DA VARIÁVEL DE FASE		
	//acrescentar antes crhs e clamb para monitorar esses valores
	P1<<<numBlocks,numThreads>>>(PcS,PcL,PcI,uc,dx,dy,EPSILONSL,EPSILONSI,EPSILONLI,WSL,WSI,WLI,MISL,MISI,MILI,LS,TmS,LL,TmL,LI,TmI,0,2,1);
	P1<<<numBlocks,numThreads>>>(PcL,PcI,PcS,uc,dx,dy,EPSILONLI,EPSILONSL,EPSILONSI,WLI,WSL,WSI,MILI,MISL,MISI,LL,TmL,LI,TmI,LS,TmS,2,1,0);
	P1<<<numBlocks,numThreads>>>(PcI,PcL,PcS,uc,dx,dy,EPSILONLI,EPSILONSI,EPSILONSL,WLI,WSI,WSL,MILI,MISI,MISL,LI,TmI,LL,TmL,LS,TmS,1,2,0);
	cudaDeviceSynchronize();
	//printf("Teste 1");
	//CALCULO DO CAMPO DE TEMPERATURAS
	Temp<<<numBlocks, numThreads,4096>>>(uc, Xc, PcS,PcL, deltac, Fox, Foy,LS,LL);
	cudaDeviceSynchronize();
	
	//FUNÇÃO DE CAPTAÇÃO DE ERROS 
	/*cudaError_t erro = cudaGetLastError();        // Get error code
		if ( erro != cudaSuccess )
		{
			printf("CUDA Error: %s\n", cudaGetErrorString(erro));
			exit(-1);
		}
	*/
	//printf("Teste 2");
	//IMPRIMIR EM DETERMINADOS INTERVALOS DE TEMPO 	
	if (tempo%10000==0)
	{
		
		FLAG=1;//flag pra pegar o valor do contorno e calcular a velocidade da interface
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
			//contorno[j][i]=PL[0][i][j]*PL[0][i][j]+PS[0][i][j]*PS[0][i][j]+PI[0][i][j]*PI[0][i][j];//contorno[-j+TELY-1][i]
			if(PS[0][j][i]<=0.6&&j==(TELX/2)&&FLAG==1){
			//printf("contorno=%f\n", contorno[i][j]);
			fprintf(arq,"%f %f \n", tempo*dt, i*dx); //ARQUIVO TEXTO
			//fprintf(curva1, "%f %f \n", tempo*dt, (2*1.8*sqrt((tempo*dt)*0.03)));
			//fprintf(curva2, "%f %f \n", tempo*dt, (2*1.8*sqrt((tempo*dt)*0.04)));
			//fprintf(curva3, "%f %f \n", tempo*dt, (2*1.8*sqrt((tempo*dt)*0.05)));
			FLAG=0;}
		}
		if(tempo==990000){
		fprintf(arqb,"%f %f \n", j*dx, X[0][TELX/2][j]); //ARQUIVO TEXTO
		//fprintf(curvaT1, "%f %f \n", j*dx, 1+1.0*erf((j*dx)/(2*sqrt(0.08*0.495)))/(erf(1.5)));
		//fprintf(curvaT2, "%f %f \n", j*dx, 1+1.0*erf((j*dx)/(2*sqrt(0.12*0.495)))/(erf(1.8)));
		//fprintf(curvaT3, "%f %f \n", j*dx, 1+1.0*erf((j*dx)/(2*sqrt(0.16*0.495)))/(erf(2.0)));
		}
		//fprintf(arq,"\n\n");
		fwrite(PS[0][j], sizeof(float), TELX, arqS);//ARQUIVO BINÁRIO
		fwrite(PL[0][j], sizeof(float), TELX, arqL);//ARQUIVO BINÁRIO
		fwrite(PI[0][j], sizeof(float), TELX, arqI);//ARQUIVO BINÁRIO
		fwrite(u[0][j], sizeof(float), TELX, arqT);//ARQUIVO BINÁRIO
	}
	
	//IMPRESSÃO DA TEMPERATURA EM DETERMINADOS PONTOS PARA OBSERVAÇÃO DURANTE A EXECUÇÃO
		printf("\n PS[50][10]=%f PS[50][50]=%f PS[50][20]=%f\n",  PS[0][10][5], PS[0][20][5],PS[0][30][5]);
		//printf("\n PL[10][10]=%f PL[190][10]=%f PL[100][100]=%f\n",  PL[0][10][10], PL[0][190][10],PL[0][100][100]);
		//printf("\n PI[10][10]=%f PI[190][10]=%f PI[100][100]=%f\n",  PI[0][10][10], PI[0][190][10],PI[0][100][100]);
		//A leitura é feita de baixo pra cima. A coordenada x é ok, mas a y é invertida.
		//printf(" contorno=%f\n", contorno[500][500]);
		printf(" u0[0][10]=%f u0[80][10]=%f\n", u[0][0][10], u[0][5][10]);
		
		

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
fclose(arqS);
fclose(arqL);
fclose(arqI);
fclose(arqT);
fclose(curva1);
fclose(curva2);
fclose(curva3);
fclose(curvaT1);
fclose(curvaT2);
fclose(curvaT3);

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
