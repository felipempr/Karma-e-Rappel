/********************************************************************************************************************************
*															        *	*																*
*		FORMA FINAL PARA SIMULAÇÃO EM DUAS DIMENSÕES DO CAMPO DE FASES COM ANISOTROPIA					*
*																*
*																*
*********************************************************************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "funcKR.h"



int main(void)
{
	int i;
	int j;
	int t;
	char p;
	int x1;
	int x0;
	
	float v;
	float constK;
	int h;
	//float alfa;
	float D;
	int telx;
	int tely;
	float Ttot;
	float dt;
	float tau0;
	float E0;//ao quadrado
	float La;
	//float TI;	
	//int TM;
	float compL;
	//float Ste;
	int r;
	float d0=0.554;

	
	
	//LEITURA DO INPUT	
	FILE *init;
	init=fopen("input3","r");
	while(p!='#')
	{fscanf(init,"%c",&p);}
	fscanf(init,"%f %d %f %f %d %d %d %f %f %f %f", &constK,&h,&E0,&compL,&r,&telx,&tely,&Ttot,&dt,&tau0,&D);
printf("\nK=%f\nh=%d\nE0=%f\nL=%f\nr=%d\ntelx=%d\ntely=%d\nTtot=%f\ndt=%.10f\ntau0=%f\nD=%f", constK,h,E0,compL,r,telx,tely,Ttot,dt,tau0,D);
	printf("\n");
	
	int telt=Ttot/dt;
	float dx=(float)compL/(float)telx;
	float dy=(float)compL/(float)tely;
	//const float sigma=E/sqrt(2*W);
	
	float ***P;
	P=(float***)malloc(2*sizeof(float**));
	for(int t=0;t<=1;t++){
		P[t]=(float**)malloc(telx*sizeof(float*));
		for(int i=0;i<=telx-1;i++){
			P[t][i]=(float*)malloc(tely*sizeof(float));
			//printf("%d %d\n",t,i);
			
		}

	}
	float ***u;
	u=(float***)malloc(2*sizeof(float**));
	for(int t=0;t<=1;t++){
		u[t]=(float**)malloc(telx*sizeof(float*));
		for(int i=0;i<=telx-1;i++){
			u[t][i]=(float*)malloc(tely*sizeof(float));
						
		}

	}
	float ***X;
	X=(float***)malloc(2*sizeof(float**));
	for(int t=0;t<=1;t++){
		X[t]=(float**)malloc(telx*sizeof(float*));
		for(int i=0;i<=telx-1;i++){
			X[t][i]=(float*)malloc(tely*sizeof(float));
						
		}
	}
	float *delta;
	delta=(float*)malloc((telx*tely)*sizeof(float));
	
	
	printf("\n ok1 \n");
	//float delta[telx*tely];
	printf("\n ok2 \n");
	float Fox=D*dt/(dx*dx);
	float Foy=D*dt/(dy*dy);
	float Bix=h*dx/constK;
	float Biy=h*dy/constK;
	
	float lambda=(D/0.6267)*(tau0/E0);
	printf("lambda=%f\n",lambda);
	
	FILE *arq;
	FILE *arq2;
	
	
printf("\ntelt=%d",telt);
//printf("\n dx=%.9f",dx);
//printf("\n Bi=%f\n",h*compL/constK);
//printf("\n Bix=%f\n",Bix);
//printf("\n Biy=%f\n",Biy);
//printf("\n Fox=%f\n",Fox);
//printf("\n Foy=%f\n",Foy);
//printf("\n Ttot=%d\n",Ttot);

arq=fopen("Dendritas2Dani","w");
arq2=fopen("TDendritas2Dani","w");

//INSERINDO VALOR INICIAL NAS MATRIZES

inicializar(0,1,0,telx,0,tely,u,-0.65);
inicializar(0,1,0,telx,0,tely,P,-1);
inicializar(0,2,0,telx,0,tely,X,-0.65);

for (i=0;i<=r;i++)//(i=telx/2-r;i<=telx/2+r;i++)
{
	for (j=0;j<=r;j++)//(j=tely/2-r;j<=tely/2+r;j++)
	{
		if((pow((i),2)+pow((j),2))<=(r*r))//((pow((i-telx/2.0),2)+pow((j-tely/2.0),2))<=(r*r))
		{
		P[0][i][j]=1;
		u[0][i][j]=0; //relacionado a Delta, super-resfriamento inicial
		}
	}
}

for(t=0;t<=telt;t++)
{

//CÁLCULO DA VARIÁVEL DE FASE		
printf("\nt=%d",t);

#pragma omp parallel shared(P,X) private(i,j)
{
#pragma omp for schedule (static)
	for(j=0;j<=tely-1;j++)
	{		
		for(i=0;i<=telx-1;i++)
		{
			P[1][i][j]=P1(0.05,P,u,i,j,telx,tely,dx,dy,dt,tau0,E0,lambda);
		}
	}
}			

//CALCULO DO CAMPO DE TEMPERATURAS
do{
#pragma omp parallel shared(X,P,u) private(i,j)
{
#pragma omp for schedule (static)
	for(j=0;j<=tely-1;j++)
	{		
		for(i=0;i<=telx-1;i++)
		{
			X[1][i][j]=X1(u,X,P,i,j,telx,tely,Fox,Foy);
		}
	}
}	
//CÁLCULO DO ERRO DA FORMULAÇÃO IMPLÍCITA DO CAMPO DE TEMPERATURAS	
		for(i=0;i<=telx-1;i++)
		{
			for(j=0;j<=tely-1;j++)
			{
				delta[j+i*telx]=X[1][i][j]-X[0][i][j];//vetor com os erros
				X[0][i][j]=X[1][i][j];	//nesse ponto X[1] se torna X[0]
			}
		}
		//printf(" %f",maior(telx*tely,delta));
}while(maior(telx*tely,delta)>0.00001);
//printf("P125 125=%f",P[0][125][125]);
//printf("T125 125=%f",u[0][125][125]);

//IMPRIMIR A CADA 100t			
if(t%1000==0){
	fprintf(arq,"phi:t=%.8f\n",(t*dt));
	fprintf(arq2,"u:t=%.3f\n",(t*dt));
	for(j=0;j<=tely-1;j++)
	{
		for(i=0;i<=telx-1;i++)
		{
		fprintf(arq,"%f %f %f\n",i*((compL/d0)/telx),j*((compL/d0)/tely), P[1][i][j]);
		fprintf(arq2,"%f %f %f\n",i*((compL/d0)/telx),j*((compL/d0)/tely), X[1][i][j]);
		}
		if (P[1][i][tely/2]==0.5)
		x1=P[1][i][tely/2]
		if (P[0][i][tely/2]==0.5)
		x0=P[0][i][tely/2]
	}
	v=x1-x0;
	x1=0;
	x0=0; //abrir outro arquivo para as velocidades
	fprintf(arq,"\n\n");
	fprintf(arq2,"\n\n");
	
	
}

//ATUALIZAR VARIÁVEIS PARA O PRÓXIMO CICLO
for(i=0;i<=telx-1;i++)
	{
		for(j=0;j<=tely-1;j++)
			{
			P[0][i][j]=P[1][i][j];	//atualizo p para o proximo ciclo
			u[0][i][j]=X[1][i][j];	//atualizo T para o proximo ciclo
			}
	}			

}

//FECHAR E ARQUIVOS E LIBERAR MEMÓRIA
fclose(arq);
fclose(arq2);

for(int t=0;t<=1;t++){
	for(int i=0;i<=telx-1;i++){
		free(P[t][i]);
		free(u[t][i]);
		free(X[t][i]);
	}
	free(P[t]);
	free(u[t]);
	free(X[t]);
}
free(P);
free(u);
free(X);
free(delta);

return 0;	
}

