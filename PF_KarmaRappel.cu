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
	//int FLAGmobi;
	
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
printf("%d, ",tempo);
//CÁLCULO DA VARIÁVEL DE FASE		
P1<<<numBlocks,numThreads>>>(PcS,PcL,PcI,uc,dx,dy,lambda,1);
P1<<<numBlocks,numThreads>>>(PcL,PcS,PcI,uc,dx,dy,lambda,0);
P1<<<numBlocks,numThreads>>>(PcI,PcL,PcS,uc,dx,dy,lambda,1);
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
			cudaMemcpy(PS[t][i], PcS+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(PL[t][i], PcL+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(PI[t][i], PcI+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(u[t][i], uc+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			cudaMemcpy(X[t][i], Xc+((t*TELX*TELY)+(i*TELX)), TELY*sizeof(float), cudaMemcpyDeviceToHost);
			}
		}
	//printf(" u0[5][0]=%f u0[30][0]=%f\n", u[0][5][0], u[0][30][0]);
	//printf(" u0[0][5]=%f u0[0][30]=%f\n", u[0][0][5], u[0][0][30]);
	fprintf(arq, "phi:t=%f\n",(tempo*dt));
	for(j=0;j<TELY;j++)
	{
		for(i=0;i<TELX;i++)
		{
			fprintf(arq,"%d %d %f\n", i,j,PI[0][i][j]);
		}
	}
	fprintf(arq,"\n\n");

	//IMPRESSÃO NO ARQUIVO BINÁRIO
	for (int j = 0; j < TELY; j++) {
		for(int i = 0; i < TELX; i++){
			//contorno[j][i]=PS[0][i][j]+PI[0][i][j]+PL[0][i][j];
			contorno[-j+TELY-1][i]=PS[0][i][j]*PS[0][i][j]+PI[0][i][j]*PI[0][i][j]+PL[0][i][j]*PL[0][i][j];
		}
		fwrite(contorno[j], sizeof(float), TELX, arq2);//contorno[j]
	}
}
//PERTO DO FINAL O AVANÇO JÁ É ESTÁVEL PONTO ONDE EXTRAIO O CONTORNO SIMULADO

//if (tempo==17500)
	//{
	//cristal(arqb, contorno, ySim);
	//}
//NO FINAL CALCULO A CURVA ANALÍTICA E COMPARO COM YSIM


//ATUALIZAR VARIÁVEIS PARA O PRÓXIMO CICLO
atualizaTudo<<<numBlocks, numThreads >>>(PcS,PcS);
atualizaTudo<<<numBlocks, numThreads >>>(PcL,PcL);
atualizaTudo<<<numBlocks, numThreads >>>(PcI,PcI);
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
			while(xoff<TELX-1){
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
		if(j<=((TELX/2)-5)&&FLAG==true){
		n1=i;
		FLAG=false;
		}
		if(j<(TELX/2)-5){
		ySim[i-n1]=-j+TELY/2;
		//printf("n1=%d i=%d\n",n1,i);
		//printf("ySim[%d]=%d\n",i-n1,-j+TELY/2);
		//fprintf(arquivo,"%d %d\n",i-n1,-j+TELY/2); //-i+TELX -j+TELY/2
		}
		

	}
	

}
