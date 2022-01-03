/********************************************************************************************************************************
*															        																
*		SIMULAÇÃO EM DUAS DIMENSÕES DO CAMPO DE FASES PARALELIZADO REPRODUZINDO O MODELO DE JOHNSON DESCRITO NO ARTIGO
*																NÃO HÁ CAMPO DE TEMPERATURA
*																
*********************************************************************************************************************************/


#define _USE_MATH_DEFINES
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
void cristal(FILE *arquivo, float **cont, float ySim[]);
void curvas(FILE *arquivo,FILE *arquivo2, float ySim[]);

int main(void)
{
	int i;
	int j;
	int tempo;
	int numThreads = 512;
	int numBlocks = (TELX*TELY + numThreads - 1)/numThreads; //3052
	
	//char p;
		
	float ySim[TELX];
	int telt=Ttot/dt;
	float dx=(float)compL/(float)TELX;
	float dy=(float)compL/(float)TELY;
	//const float sigma=E/sqrt(2*W);
	
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
	float ***P0;
	P0=(float***)malloc(2*sizeof(float**));
	for(int t=0;t<=1;t++){
		P0[t]=(float**)malloc(TELX*sizeof(float*));
		for(int i=0;i<=TELX-1;i++){
			P0[t][i]=(float*)malloc(TELY*sizeof(float));
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
	float *Pc0;	
	cudaMalloc((void **) &Pc0, 2*TELX*TELY*sizeof(float));
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
	//FILE *arq3;
	
	
printf("\ntelt=%d\n",telt);

arq=fopen("JohnsonPS","w");
arqb=fopen("JohnsonPSsim","w");
arq2=fopen("JohnsonPSbin","wb");
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
inicializar(0,1,0,TELX,0,TELY,P0,0);
inicializar(0,2,0,TELX,0,TELY,X,-Ui);

for (i=0;i<TELX;i++)//(i=telx/2-r;i<=telx/2+r;i++)
{
	for (j=0;j<TELY/2;j++)//(j=tely/2-r;j<=tely/2+r;j++)
	{
		if(j>(sqrt(3)/3)*i+TELY/2-(sqrt(3)/3)*(TELX/2))
		PA[0][i][j]=1;
		if(j<(sqrt(3)/3)*i+TELY/2-(sqrt(3)/3)*(TELX/2)&&i<TELX/2)
		PB[0][i][j]=1;
		if(j<(-sqrt(3)/3)*i+TELY/2-(sqrt(3)/3)*(TELX/2)&&i>TELX/2)
		PC[0][i][j]=1;
		if(j>(-sqrt(3)/3)*i+TELY/2-(sqrt(3)/3)*(TELX/2))
		PD[0][i][j]=1;
		if((pow((i),2)+pow((j),2))<=(R*R)){
		P0[0][i][j]=1;
		PA[0][i][j]=0;
		PB[0][i][j]=0;
		PC[0][i][j]=0;
		PD[0][i][j]=0;
		PE[0][i][j]=0;
		PF[0][i][j]=0;}
	}
		/*//if((pow((i),2)+pow((j),2))<=(R*R))//((pow((i-telx/2.0),2)+pow((j-tely/2.0),2))<=(r*r))
		if(i<((j+((3.0*(TELY-1))/4.0))/(3.0*(TELY-1)/(TELX-1)))){
		PI[0][i][j]=0;
		PS[0][i][j]=1;
		PL[0][i][j]=0;
		}
		else if(i>((j-((9.0*(TELY-1))/4.0))/-((3.0*(TELY-1))/(TELX-1)))){
		PI[0][i][j]=0;
		PS[0][i][j]=0;
		PL[0][i][j]=1;
		}
		else{
		PI[0][i][j]=1;
		PS[0][i][j]=0;
		PL[0][i][j]=0;
		}*/

	}

	for(j=TELY/2;j<TELY;j++)
	{
		if(j<(-sqrt(3)/3)*i+TELY/2-(sqrt(3)/3)*(TELX/2))
		PA[0][i][j]=1;
		if(j>(-sqrt(3)/3)*i+TELY/2-(sqrt(3)/3)*(TELX/2)&&i<TELX/2)
		PF[0][i][j]=1;
		if(j>(sqrt(3)/3)*i+TELY/2-(sqrt(3)/3)*(TELX/2)&&i>TELX/2)
		PE[0][i][j]=1;
		if(j<(sqrt(3)/3)*i+TELY/2-(sqrt(3)/3)*(TELX/2))
		PD[0][i][j]=1;
		if((pow((i),2)+pow((j),2))<=(R*R)){
		P0[0][i][j]=1;
		PA[0][i][j]=0;
		PB[0][i][j]=0;
		PC[0][i][j]=0;
		PD[0][i][j]=0;
		PE[0][i][j]=0;
		PF[0][i][j]=0;}
		/*if(i<((TELX-1)/2)){
		PI[0][i][j]=0;
		PS[0][i][j]=1;
		PL[0][i][j]=0;
		}
		else{
		PI[0][i][j]=0;
		PS[0][i][j]=0;
		PL[0][i][j]=1;
		}*/
	}
}
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
	cudaMemcpy(Pc0+((t*TELX*TELY)+(i*TELX)), P0[t][i], TELY*sizeof(float), cudaMemcpyHostToDevice);
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
printf("%d, ",tempo);
//CÁLCULO DA VARIÁVEL DE FASE		
P1<<<numBlocks,numThreads>>>(Pc0,PcA,PcB,PcC,PcD,PcE,PcF,uc,dx,dy,lambda);
P1<<<numBlocks,numThreads>>>(PcA,PcB,PcC,PcD,PcE,PcF,Pc0,uc,dx,dy,lambda);
P1<<<numBlocks,numThreads>>>(PcB,PcC,PcD,PcE,PcF,Pc0,PcA,uc,dx,dy,lambda);
P1<<<numBlocks,numThreads>>>(PcC,PcD,PcE,PcF,Pc0,PcA,PcB,uc,dx,dy,lambda);
P1<<<numBlocks,numThreads>>>(PcD,PcE,PcF,Pc0,PcA,PcB,PcC,uc,dx,dy,lambda);
P1<<<numBlocks,numThreads>>>(PcE,PcF,Pc0,PcA,PcB,PcC,PcD,uc,dx,dy,lambda);
P1<<<numBlocks,numThreads>>>(PcF,Pc0,PcA,PcB,PcC,PcD,PcE,uc,dx,dy,lambda);
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
if (tempo%5000==0)
	{//%5000
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
			cudaMemcpy(P0[t][i], Pc0+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
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
			contorno[-j+TELY-1][i]=PA[0][i][j]*PA[0][i][j]+PB[0][i][j]*PB[0][i][j]+PC[0][i][j]*PC[0][i][j]+PD[0][i][j]*PD[0][i][j]+PE[0][i][j]*PE[0][i][j]+PF[0][i][j]*PF[0][i][j]+P0[0][i][j]*P0[0][i][j];
		}
		fwrite(contorno[j], sizeof(float), TELX, arq2);//contorno[j]
	}
}
//PERTO DO FINAL O AVANÇO JÁ É ESTÁVEL PONTO ONDE EXTRAIO O CONTORNO SIMULADO

//if (tempo==25000)//35000 paraLTJ 0.01
	//{
	//cristal(arqb, contorno, ySim);
	//}
//NO FINAL CALCULO A CURVA ANALÍTICA E COMPARO COM YSIM


//ATUALIZAR VARIÁVEIS PARA O PRÓXIMO CICLO
atualizaTudo<<<numBlocks, numThreads >>>(PcA,PcA);
atualizaTudo<<<numBlocks, numThreads >>>(PcB,PcB);
atualizaTudo<<<numBlocks, numThreads >>>(PcC,PcC);
atualizaTudo<<<numBlocks, numThreads >>>(PcD,PcD);
atualizaTudo<<<numBlocks, numThreads >>>(PcE,PcE);
atualizaTudo<<<numBlocks, numThreads >>>(PcF,PcF);
atualizaTudo<<<numBlocks, numThreads >>>(Pc0,Pc0);
atualizaTudo<<<numBlocks, numThreads >>>(uc,Xc);
}

//curvas(arq,arqb,ySim);

printf("lambda=%f\n",lambda);
printf("d0=%f\n",(0.8839*E0/lambda));
printf("E0/d0=%f\n",(D/(0.8839*0.6267)));

//FECHAR E ARQUIVOS E LIBERAR MEMÓRIA
fclose(arq);
fclose(arq2);
//fclose(arq3);
for(int t=0;t<2;t++){
	for(int i=0;i<TELX;i++){
		
		free(PA[t][i]);
		free(PB[t][i]);
		free(PC[t][i]);
		free(PD[t][i]);
		free(PE[t][i]);
		free(PF[t][i]);
		free(P0[t][i]);
		free(u[t][i]);
		free(X[t][i]);
	}
	free(PA[t]);
	free(PB[t]);
	free(PC[t]);
	free(PD[t]);
	free(PE[t]);
	free(PF[t]);
	free(P0[t]);
	free(u[t]);
	free(X[t]);
}
free(PA);
free(PB);
free(PC);
free(PD);
free(PE);
free(PF);
free(P0);
free(u);
free(X);

cudaFree(PcA);
cudaFree(PcB);
cudaFree(PcC);
cudaFree(PcD);
cudaFree(PcE);
cudaFree(PcF);
cudaFree(Pc0);
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

//CALCULO A CURVA ANALÍTICA
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
void cristal(FILE *arquivo, float **cont, float ySim[])
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
