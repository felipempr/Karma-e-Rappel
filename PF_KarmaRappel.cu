/********************************************************************************************************************************
*															        																
*		SIMULAÇÃO EM DUAS DIMENSÕES DO CAMPO DE FASES PARALELIZADO REPRODUZINDO O MODELO DE JOHNSON DESCRITO NO ARTIGO
*																NÃO HÁ CAMPO DE TEMPERATURA
*																
*********************************************************************************************************************************/


#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
void cristal(float **cont, float ySim[]);
void curvas(FILE *arquivo,FILE *arquivo2, float ySim[]);
void cond_inicial(unsigned char *data,float ***QA,float ***QB,float ***QC,float ***QG, float ***QD,float ***QE);
void readBMP(unsigned char *data);


int main(void)
{
printf("teste %d ",(11%9)/3);

	int i;
	int j;
	int tempo;
	int numThreads = 512;
	int numBlocks = (TELX*TELY + numThreads - 1)/numThreads; //3052
	int A[2]={0,0};
	//float delA;
	
	//char p;
		
	float ySim[TELX];
	int telt=Ttot/dt;
	float dx=(float)compL/(float)TELX;
	float dy=(float)compL/(float)TELY;
	//const float sigma=E/sqrt(2*W);

	unsigned char* bmp; 
	bmp=(unsigned char*)malloc(3*TELX*TELY*sizeof(unsigned char));
	
	//VARIÁVEIS DAS FASES LÍQUIDA, SÓLIDA E IMÓVEL(?)
	float ***PA;
	PA=(float***)malloc(2*sizeof(float**));
	for(int t=0;t<=1;t++){
		PA[t]=(float**)malloc(TELX*sizeof(float*));
		for(int i=0;i<=TELX-1;i++){
			PA[t][i]=(float*)malloc(TELY*sizeof(float));
		}
	}
	float ***PB;
	PB=(float***)malloc(2*sizeof(float**));
	for(int t=0;t<=1;t++){
		PB[t]=(float**)malloc(TELX*sizeof(float*));
		for(int i=0;i<=TELX-1;i++){
			PB[t][i]=(float*)malloc(TELY*sizeof(float));
		}
	}
	float ***PC;
	PC=(float***)malloc(2*sizeof(float**));
	for(int t=0;t<=1;t++){
		PC[t]=(float**)malloc(TELX*sizeof(float*));
		for(int i=0;i<=TELX-1;i++){
			PC[t][i]=(float*)malloc(TELY*sizeof(float));
		}
	}
	float ***PD;
	PD=(float***)malloc(2*sizeof(float**));
	for(int t=0;t<=1;t++){
		PD[t]=(float**)malloc(TELX*sizeof(float*));
		for(int i=0;i<=TELX-1;i++){
			PD[t][i]=(float*)malloc(TELY*sizeof(float));
		}
	}
	float ***PE;
	PE=(float***)malloc(2*sizeof(float**));
	for(int t=0;t<=1;t++){
		PE[t]=(float**)malloc(TELX*sizeof(float*));
		for(int i=0;i<=TELX-1;i++){
			PE[t][i]=(float*)malloc(TELY*sizeof(float));
		}
	}
	float ***PF;
	PF=(float***)malloc(2*sizeof(float**));
	for(int t=0;t<=1;t++){
		PF[t]=(float**)malloc(TELX*sizeof(float*));
		for(int i=0;i<=TELX-1;i++){
			PF[t][i]=(float*)malloc(TELY*sizeof(float));
		}
	}
	float ***PG;
	PG=(float***)malloc(2*sizeof(float**));
	for(int t=0;t<=1;t++){
		PG[t]=(float**)malloc(TELX*sizeof(float*));
		for(int i=0;i<=TELX-1;i++){
			PG[t][i]=(float*)malloc(TELY*sizeof(float));
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
	float *PcA;	
	cudaMalloc((void **) &PcA, 2*TELX*TELY*sizeof(float));
	float *PcB;	
	cudaMalloc((void **) &PcB, 2*TELX*TELY*sizeof(float));
	float *PcC;	
	cudaMalloc((void **) &PcC, 2*TELX*TELY*sizeof(float));
	float *PcD;	
	cudaMalloc((void **) &PcD, 2*TELX*TELY*sizeof(float));
	float *PcE;	
	cudaMalloc((void **) &PcE, 2*TELX*TELY*sizeof(float));
	float *PcF;	
	cudaMalloc((void **) &PcF, 2*TELX*TELY*sizeof(float));
	float *PcG;	
	cudaMalloc((void **) &PcG, 2*TELX*TELY*sizeof(float));
	float *uc;	
	cudaMalloc((void **) &uc, 2*TELX*TELY*sizeof(float));
	float *Xc;	
	cudaMalloc((void **) &Xc, 2*TELX*TELY*sizeof(float));
	float *deltac;
	cudaMalloc((void **) &deltac, TELX*TELY*sizeof(float));
	
	//CÁLCULO DOS NÚMEROS DE FOURIER E BIOT
	//float Fox=D*dt/(dx*dx);
	//float Foy=D*dt/(dy*dy);
	
	float lambda=(1/a2)*((D*TAU0)/(E0*E0));//*(TAU0/(E0*E0));
	printf("lambda=%f\n",lambda);
	printf("d0=%f\n",(a1*E0/lambda));
	printf("E0/d0=%f\n",E0*(lambda/a1*E0));
	
	FILE *arq;
	FILE *arqb;
	FILE *arq2;
	FILE *arq2b;
	
printf("\ntelt=%d\n",telt);

arq=fopen("JohnsonPSanal","w");
arqb=fopen("JohnsonPSsim","w");

arq2=fopen("JohnsonPSbin","wb");
arq2b=fopen("JohnsonPSDel","w");
//arq=fopen("E_d08","w");//para binário "wb"
//arq3=fopen("E_d08i","w");//para binário "wb"
//arq2=fopen("DendritasParalE_d08","wb");

//INSERINDO VALOR INICIAL NAS MATRIZES

inicializar(0,1,0,TELX,0,TELY,u,-Ui);
inicializar(0,1,0,TELX,0,TELY,PA,0);
inicializar(0,1,0,TELX,0,TELY,PB,0);
inicializar(0,1,0,TELX,0,TELY,PC,0);
inicializar(0,1,0,TELX,0,TELY,PD,0);
inicializar(0,1,0,TELX,0,TELY,PE,0);
inicializar(0,1,0,TELX,0,TELY,PF,0);
inicializar(0,1,0,TELX,0,TELY,PG,0);
inicializar(0,2,0,TELX,0,TELY,X,-Ui);

//TRECHO ANTIGO FEITO COM GEOMETRIA ANALITICA
/*for (i=0;i<TELX;i++)//(i=telx/2-r;i<=telx/2+r;i++)
{
	for (j=0;j<TELY;j++)//(j=tely/2-r;j<=tely/2+r;j++)
	{*/
		//TRECHO ABAIXO PARA UM QUADRADO
		/*if(i<TELX/2&&j<TELY/2)
		PA[0][i][j]=1;
		if(i<TELX/2&&j>TELY/2)
		PB[0][i][j]=1;
		if(i>TELX/2&&j<TELY/2)
		PC[0][i][j]=1;
		if(i>TELX/2&&j>TELY/2)
		PD[0][i][j]=1;
		if((pow((i-TELX/2),2)+pow((j-TELY/2),2))<=(R*R)){
		PG[0][i][j]=1;
		PA[0][i][j]=0;
		PB[0][i][j]=0;
		PC[0][i][j]=0;
		PD[0][i][j]=0;
		PE[0][i][j]=0;
		PF[0][i][j]=0;}*/

		//TRECHO ABAIXO PARA UM HEXÁGONO
		/*if(j>=(sqrt(3)/3.0)*i+TELY/2-(sqrt(3)/3.0)*(TELX/2)&&j<-(sqrt(3)/3.0)*i+TELY/2+(sqrt(3)/3.0)*(TELX/2))
		PA[0][i][j]=1;
		if(j<(sqrt(3)/3.0)*i+TELY/2-(sqrt(3)/3.0)*(TELX/2)&&i<TELX/2)
		PB[0][i][j]=1;
		if(j<-(sqrt(3)/3.0)*i+TELY/2+(sqrt(3)/3.0)*(TELX/2)&&i>=TELX/2)
		PC[0][i][j]=1;
		if(j>=-(sqrt(3)/3.0)*i+TELY/2+(sqrt(3)/3.0)*(TELX/2)&&j<(sqrt(3)/3.0)*i+TELY/2-(sqrt(3)/3.0)*(TELX/2))
		PA[0][i][j]=1;
		if(j>=(sqrt(3)/3.0)*i+TELY/2-(sqrt(3)/3.0)*(TELX/2)&&i>=TELX/2)
		PC[0][i][j]=1;
		if(j>=-(sqrt(3)/3.0)*i+TELY/2+(sqrt(3)/3.0)*(TELX/2)&&i<TELX/2)
		PB[0][i][j]=1;
		if((pow((i-TELX/2),2)+pow((j-TELY/2),2))<=(R*R)){//&&i>TELX/2
		PG[0][i][j]=1;
		PA[0][i][j]=0;
		PB[0][i][j]=0;
		PC[0][i][j]=0;
		PD[0][i][j]=0;
		PE[0][i][j]=0;
		PF[0][i][j]=0;}
	};
}*/
readBMP(bmp);
printf("bmp[0]=%d",(int)bmp[0]);
cond_inicial(bmp,PA,PB,PC,PG,PD,PE);		

//printf(" PS[1][1]=%f PS[3][6]=%f PS[6][6]=%f\n", PS[0][1][1], PS[0][3][6], PS[0][6][6]);

//COPIANDO VALORES DA CPU (HOST) PARA AS VARIÁVEIS DA GPU (DEVICE)
for(int t=0;t<2;t++)
{
	for(int i=0;i<TELX;i++)
	{
	cudaMemcpy(PcA+((t*TELX*TELY)+(i*TELX)), PA[t][i], TELY*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(PcB+((t*TELX*TELY)+(i*TELX)), PB[t][i], TELY*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(PcC+((t*TELX*TELY)+(i*TELX)), PC[t][i], TELY*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(PcD+((t*TELX*TELY)+(i*TELX)), PD[t][i], TELY*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(PcE+((t*TELX*TELY)+(i*TELX)), PE[t][i], TELY*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(PcF+((t*TELX*TELY)+(i*TELX)), PF[t][i], TELY*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(PcG+((t*TELX*TELY)+(i*TELX)), PG[t][i], TELY*sizeof(float), cudaMemcpyHostToDevice);
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
//CÁLCULO DA VARIÁVEL DE FASE		
P1<<<numBlocks,numThreads>>>(PcG,PcA,PcB,PcC,PcD,PcE,PcF,uc,dx,dy,lambda);
P1<<<numBlocks,numThreads>>>(PcA,PcB,PcC,PcD,PcE,PcF,PcG,uc,dx,dy,lambda);
P1<<<numBlocks,numThreads>>>(PcB,PcC,PcD,PcE,PcF,PcG,PcA,uc,dx,dy,lambda);
P1<<<numBlocks,numThreads>>>(PcC,PcD,PcE,PcF,PcG,PcA,PcB,uc,dx,dy,lambda);
P1<<<numBlocks,numThreads>>>(PcD,PcE,PcF,PcG,PcA,PcB,PcC,uc,dx,dy,lambda);
P1<<<numBlocks,numThreads>>>(PcE,PcF,PcG,PcA,PcB,PcC,PcD,uc,dx,dy,lambda);
P1<<<numBlocks,numThreads>>>(PcF,PcG,PcA,PcB,PcC,PcD,PcE,uc,dx,dy,lambda);
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
if (tempo%250==0)//%5000
	{printf("\t%d, ",tempo);
	A[1]=0;
		for (int t=0;t<2;t++)
		{
			for (int i=0;i<TELX;i++)
			{
			cudaMemcpy(PA[t][i], PcA+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(PB[t][i], PcB+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(PC[t][i], PcC+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(PD[t][i], PcD+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(PE[t][i], PcE+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(PF[t][i], PcF+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(PG[t][i], PcG+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(u[t][i], uc+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(X[t][i], Xc+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			}
		}
	//IMPRESSÃO NOS ARQUIVOS
	fprintf(arq, "phi:t=%f\n",(tempo*dt));
	for(j=0;j<TELY;j++)
	{
		for(i=0;i<TELX;i++)
		{
			//contorno[j][i]=PS[0][i][j]+PI[0][i][j]+PL[0][i][j];
			contorno[j][i]=PA[0][i][j]*PA[0][i][j]+PB[0][i][j]*PB[0][i][j]+PC[0][i][j]*PC[0][i][j]+PD[0][i][j]*PD[0][i][j]+PE[0][i][j]*PE[0][i][j]+PF[0][i][j]*PF[0][i][j]+PG[0][i][j]*PG[0][i][j];//contorno[-j+TELY-1][i]
			fprintf(arq,"%d %d %f\n", i,j,0); //ARQUIVO TEXTO
			
			if(PG[1][i][j]>0.9)
			A[1]=A[1]+1;
			
		}
		fwrite(contorno[j], sizeof(float), TELX, arq2);//ARQUIVO BINÁRIO
	}
	fprintf(arq,"\n\n");
	//delA=(A[1]-A[0])/(dt*1000);
	fprintf(arq2b,"%f %f\n",(tempo*dt),(A[1]*dx*dy));
	//fwrite(&A[1], sizeof(int),1, arq2b);
	A[0]=A[1];
}
//PERTO DO FINAL O AVANÇO JÁ É ESTÁVEL PONTO ONDE EXTRAIO O CONTORNO SIMULADO

//if (tempo==25000)//35000 paraLTJ 0.01
	//{
	//cristal(contorno, ySim);
	//}
//NO FINAL CALCULO A CURVA ANALÍTICA E COMPARO COM YSIM


//ATUALIZAR VARIÁVEIS PARA O PRÓXIMO CICLO
atualizaTudo<<<numBlocks, numThreads >>>(PcA,PcA);
atualizaTudo<<<numBlocks, numThreads >>>(PcB,PcB);
atualizaTudo<<<numBlocks, numThreads >>>(PcC,PcC);
atualizaTudo<<<numBlocks, numThreads >>>(PcD,PcD);
atualizaTudo<<<numBlocks, numThreads >>>(PcE,PcE);
atualizaTudo<<<numBlocks, numThreads >>>(PcF,PcF);
atualizaTudo<<<numBlocks, numThreads >>>(PcG,PcG);
atualizaTudo<<<numBlocks, numThreads >>>(uc,Xc);
}

//curvas(arq,arqb,ySim);

printf("lambda=%f\n",lambda);
printf("d0=%f\n",(0.8839*E0/lambda));
printf("E0/d0=%f\n",(D/(0.8839*0.6267)));

//FECHAR E ARQUIVOS E LIBERAR MEMÓRIA
fclose(arq);
fclose(arq2);
fclose(arq2b);
for(int t=0;t<2;t++){
	for(int i=0;i<TELX;i++){
		
		free(PA[t][i]);
		free(PB[t][i]);
		free(PC[t][i]);
		free(PD[t][i]);
		free(PE[t][i]);
		free(PF[t][i]);
		free(PG[t][i]);
		free(u[t][i]);
		free(X[t][i]);
	}
	free(PA[t]);
	free(PB[t]);
	free(PC[t]);
	free(PD[t]);
	free(PE[t]);
	free(PF[t]);
	free(PG[t]);
	free(u[t]);
	free(X[t]);
}
free(PA);
free(PB);
free(PC);
free(PD);
free(PE);
free(PF);
free(PG);
free(u);
free(X);
free(bmp);

cudaFree(PcA);
cudaFree(PcB);
cudaFree(PcC);
cudaFree(PcD);
cudaFree(PcE);
cudaFree(PcF);
cudaFree(PcG);
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
//FUNÇÃO PARA GERAR AS CURVAS ANALÍTICAS E FAZER A COMPARAÇÃO
void curvas(FILE *arquivo,FILE *arquivo2, float ySim[])
{
bool FLAG=true;
bool FLAG2=true;
int x;
int xoff;
int n2=0;
int nf;
int af;
float thetaf;
float y;
float eta;
float c1;
float c2;
float theta;
float sumquad=0;
float reg;
float yAna[TELX];
float dif[TELX];

//CALCULO A CURVA ANALÍTICA E IMPRES~SÃO DE AMBAS EM SEUS RESPECTIVOS ARQUIVOS
	for(int a=80;a<=200;a++)
	{
		for(theta=M_PI/30;theta<=M_PI/3;theta=theta+M_PI/300)
		{
		FLAG=true;
		sumquad=0;
		x=0;
		xoff=0;
		eta=a/(2*theta);
		c1=logf(sin(theta));
		c2=-eta*(M_PI/2-theta);
			while(xoff<TELX-2){
			x=x+1;
			//printf("x=%d xoff=%d\n",x,xoff);
			//for(int x=0;x<TELX;x++)
			y=eta*acos(exp((-x/eta)-c1))+c2;
			if(y>0&&FLAG==true)//Só pra estabelecer onde começa o x
				{
				n2=x;
				FLAG=false;
				}
				xoff=x-n2;
				if(xoff>=0){//condição para o vetor existir e para não dar os resultados no final da linha
				yAna[xoff]=y;
				if(yAna[xoff]>0&&ySim[xoff]>0){
				//fprintf(arquivo,"%f %f\n",x-n2,y);
				//printf("yAna[%d]=%f ySim[%d]=%d\n",d,yAna[d],d,ySim[d]);
				dif[xoff]=powf((yAna[xoff]-ySim[xoff]),2.0);
				sumquad=sumquad+dif[xoff];
				//printf("dif[%d]=%f\n",d,dif[d]);
				}
				}
				//printf("sumquad=%f ",sumquad);
			}
			if(FLAG2==true&&sumquad>0){//primeiro mvalor atribuido à reg
			reg=sumquad;
			FLAG2=false;
			}
			if(sumquad<reg&&sumquad>0){
			//printf("sum menor que reg ");
			reg=sumquad;
			nf=n2;
			af=a;
			thetaf=theta;
			printf("nf=%d n2=%d, a=%d af=%d, theta=%f reg=%f\n",nf,n2,a,af,theta*180/M_PI,reg);
			}
		printf(".");
		}
		printf("X");
	}
	//Depois que é decidido o theta e o a que fazem coincidir as duas curvas as imprimo
		eta=af/(2*thetaf);
		c1=logf(sin(thetaf));
		c2=-eta*(M_PI/2-thetaf);
		//printf("af=%d thetaf=%f eta=%f c1=%f c2=%f ln10=%f",af,thetaf,eta,c1,c2,log(10));
		for(int x=0;x<TELX;x++){
		//printf("n2=%d",nf);
		y=eta*acos(exp((-x/eta)-c1))+c2;
		//printf("[%d %f]",x-nf,y);
		fprintf(arquivo,"%d %f\n",x-nf,y);
		fprintf(arquivo2,"%d %f\n",x,ySim[x]); //-i+TELX -j+TELY/2
		}
}
//FUNÇÃO PARA OBTER VALORES DO CONTORNO SIMULADO
void cristal(float **cont, float ySim[]) //como não estava usando, retirei o parâmetro FILE
{
bool FLAG=true;

int j=0;
int n1=0;

	//OBTENHO A CURVA SIMULADA
	for(int i=0;i<TELX;i++)
	{
		j=TELX/2;
		do{
		j=j-1;
		}
		while(cont[i][j]>0.51&&j>0);
		if(j<((TELX/2)-1)&&FLAG==true){
		n1=i;
		FLAG=false;
		}
		if(j<((TELX/2)-1)){
		ySim[i-n1]=-j+TELY/2;
		//printf("n1=%d i=%d\n",n1,i);
		//printf("ySim[%d]=%d\n",i-n1,-j+TELY/2);
		//fprintf(arquivo,"%d %d\n",i-n1,-j+TELY/2); //-i+TELX -j+TELY/2
		}
	}
}
void readBMP(unsigned char* data)
{
    int v;
    FILE* f = fopen("C:\\Users\\Felipe Ribeiro\\Pictures\\teste.bmp", "rb");

    if(f == NULL)
        throw "Argument Exception";

    unsigned char info[54];
    fread(info, sizeof(unsigned char), 54, f); // read the 54-byte header

    // extract image height and width from header
    int size_of_header=*(int*)&info[14];
	int width = *(int*)&info[18];
    int height = *(int*)&info[22];
	int bits_per_pixel=*(int*)&info[28];
	int number_of_colors=*(int*)&info[46];
		
	printf("largura:%d altura:%d\ntamanho do cabeçalho:%d\nbits por pixel:%d\nnumero de cores:%d\n",width, height,size_of_header,bits_per_pixel,number_of_colors);
	
/*
    cout << endl;
    cout << "  Name: " << filename << endl;
    cout << " Width: " << width << endl;
    cout << "Height: " << height << endl;
	*/
    int size = 3 * TELX * TELY; //era width*height
    //unsigned char* data; 
	//data=(unsigned char*)malloc(size*sizeof(unsigned char));
	
    // read the rest of the data at once
    fread(data, sizeof(unsigned char), size, f); 
	for(v = 0; v < size; v=v+3)
    {
            // flip the order of every 3 bytes
            int tmp = data[v];
            data[v] = data[v+2];
            data[v+2] = tmp;
			
			//printf("R:%d G:%d B:%d\n",(int)data[v],data[v+1],data[v+2]);
	}
	printf("R:%d G:%d B:%d\n",(int)data[0],(int)data[1],(int)data[2]);
	printf("R:%d G:%d B:%d\n",(int)data[1497],(int)data[1498],(int)data[1499]);
	fclose(f);
}
void cond_inicial(unsigned char *data, float ***QA,float ***QB,float ***QC,float ***QG, float ***QD, float ***QE)
{
	for(int i=0;i<=3*TELX*TELY-3;i=i+3){
		//printf("%d ",i);	
		if((int)data[i]==0&&(int)data[i+1]==0&&(int)data[i+2]==0){//preto  BGR
		//printf("i=%d j=%d\n",(i/3)%TELX,i/(3*TELX));
		QA[0][(i/3)%TELX][i/(3*TELX)]=1;
		}
		else if((int)data[i]==255&&(int)data[i+1]==255&&(int)data[i+2]==255){//branco
		QB[0][(i/3)%TELX][i/(3*TELX)]=1;
		}
		else if((int)data[i]==63&&(int)data[i+1]==72&&(int)data[i+2]==204){//azul
		QC[0][(i/3)%TELX][i/(3*TELX)]=1;
		}
		else if((int)data[i]==255&&(int)data[i+1]==242&&(int)data[i+2]==0){//amarelo else if(*(int*)(data+i)==36&&*(int*)(data+(i+1))==28&&*(int*)(data+(i+2))==237)
		QD[0][(i/3)%TELX][i/(3*TELX)]=1;
		}
		else if((int)data[i]==34&&(int)data[i+1]==177&&(int)data[i+2]==76){//verde
		QE[0][(i/3)%TELX][i/(3*TELX)]=1;
		}
		else{
		QG[0][(i/3)%TELX][i/(3*TELX)]=1;
		}
	}
	printf("OKcond_inicial");
}