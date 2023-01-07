//__global__ void testef(float *teste);

__device__ float maior(float *delta);
__device__ void deltaf(float *delta, float *X);
__device__ void atualiza(double *Q, double *K);
__global__ void atualizaTudo(double *Q, double *K);
__device__ void atualizafloat(float *Q, float *K);
__global__ void atualizaTudofloat(float *Q, float *K);
__global__ void P1(float *Pa, float *Pb, float *Pc, float *uf, double dx, double dy, double LAMBDAAB, double LAMBDAAC,double LAMBDABC, double SIGMAAB, double SIGMAAC,double SIGMABC, double miAB, double miAC,double miBC, double LAB, double LAC,double LBC, double TmAB, double TmAC,double TmBC, int FLAG1, int FLAG2, int FLAG3);
__device__ float rhs(float *Pa, float *Pb, float *Pc, float *uf, int idx, double dx, double dy, double LAMBDAAB, double LAMBDAAC,double LAMBDABC, double SIGMAAB,double SIGMAAC,double SIGMABC, double LAB, double LAC,double LBC, double TmAB, double TmAC,double TmBC, int FLAG1, int FLAG2, int FLAG3,double miAB, double miAC, double miBC);
__device__ float gradeng(float *Pa, float *Pb, float *Pc, int idx, double dx, double dy, double LAMBDAAB, double LAMBDAAC,double LAMBDABC, double SIGMAAB,double SIGMAAC,double SIGMABC, double miAB, double miAC,double miBC,int FLAG1);
//__device__ float DaDphi(float *Pa, double *Pb,double *Pc, int idx, float dx, float dy, float gammaConstAB, float gammaConstAC, int FLAG1, int FLAG2, int FLAG3);
//__device__ float divDaDNphi(float *Pa, double *Pb, double *Pc, int idx, float dx, float dy, float gammaConstAB, float gammaConstAC,float constmAB,float constmAC, int FLAG1, int FLAG2, int FLAG3);
//__device__ float MAB(float *Pa, float *Pb, int idx, double miAB);
//__device__ float MAC(float *Pa, float *Pc, int idx, double miAC);
//__device__ float MBC(float *Pb, float *Pc, int idx, double miBC);
//__device__ float DwDphi(float *Pa,float *Pb,float *Pc, int idx, float WAB,float WAC,float WBC,float MABf, float MACf, float MBCf);
//__device__ float lambda(float *Pa, float *Pb, float *Pc,float *uf, int idx, double dx, double dy, double epsilonAB, double epsilonAC, double epsilonBC, float WAB,float WAC,float WBC, float L1, float Tm1, float L2, float Tm2, float L3, float Tm3, int FLAG1, int FLAG2, int FLAG3,double miAB, double miAC, double miBC);
__device__ float freng(float *Pa, float *Pb, float *Pc, float *uf, int idx,double LAMBDAAB, double LAMBDAAC,double LAMBDABC, double SIGMAAB,double SIGMAAC,double SIGMABC, double LAB, double LAC,double LBC, double TmAB, double TmAC,double TmBC, double miAB, double miAC,double miBC,int FLAG);
__global__ void Temp(float *u, float *X, float *Pa, float *Pb, float *delta, float Fox, float Foy, double L1);
__device__ float X1(int idx,float *u, float *X, float *Pa,float *Pb, float Fox, float Foy, double L1);
//__device__ float tau(float *P, int idx, float dx, float dy);
//__device__ float E(float *P, int idx, float dx, float dy);
__device__ float Px(float *P, int idx, double dx);
__device__ float Py(float *P, int idx, double dy);
__device__ float Pxx(float *P, int idx, double dx);
__device__ float Pyy(float *P, int idx, double dy);
__device__ float Pxy(float *P, int idx, double dx, float dy);
__device__ float Pyx(float *P, int idx, double dx, float dy);
/*
__device__ float Epx(float *P, int idx, float dx, float dy);
__device__ float Epy(float *P, int idx, float dx, float dy);
__device__ float Ex(float *P, int idx, float dx, float dy);
__device__ float Ey(float *P, int idx, float dx, float dy);
__device__ float Epxpx(float *P, int idx, float dx, float dy);
__device__ float Epypy(float *P, int idx, float dx, float dy);
__device__ float Epxpy(float *P, int idx, float dx, float dy);
__device__ float Epypx(float *P, int idx, float dx, float dy);
__device__ float Expx(float *P,int idx,float dx, float dy);
__device__ float Eypy(float *P,int idx,float dx, float dy);
*/
void inicializar(int I1,int I2,int J1,int J2,int K1,int K2, double ***Q, double B);
void inicializarfloat(int I1,int I2,int J1,int J2,int K1,int K2, float ***Q, float B);
void readBMP(unsigned char* data);
void cond_inicial(unsigned char *data, float ***QA, float ***QB, float ***QC);//,float ***QG, float ***QD, float ***QE);
void curvas(FILE *arquivo,FILE *arquivo2, float ySim[]);
void cristal(float **cont, float ySim[]);//(FILE *arquivo, float **cont, float ySim[]);

//FUNÇÕES DE SUPORTE
void inicializar(int I1,int I2,int J1,int J2,int K1,int K2, double ***Q, double B)
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
void inicializarfloat(int I1,int I2,int J1,int J2,int K1,int K2, float ***Q, float B)
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

void readBMP(unsigned char* data)
{
    int v;
    FILE* f = fopen("C:\\Users\\Felipe Ribeiro\\Pictures\\teste1500.bmp", "rb");

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

void cond_inicial(unsigned char *data, float ***QA, float ***QB, float ***QC)//,float ***QG, float ***QD, float ***QE)
{
	for(int i=0;i<=3*TELX*TELY-3;i=i+3){
		//printf("%d ",i);	
		if((int)data[i]==255&&(int)data[i+1]==255&&(int)data[i+2]==255){//branco
		//printf("i=%d j=%d\n",(i/3)%TELX,i/(3*TELX));
		QA[0][(i/3)%TELX][i/(3*TELX)]=EXISTE;
		//u[0][(i/3)%TELX][i/(3*TELX)]=T;
		}
		else if((int)data[i]==0&&(int)data[i+1]==0&&(int)data[i+2]==0){//preto  BGR
		QB[0][(i/3)%TELX][i/(3*TELX)]=EXISTE;
		}
		/*else if((int)data[i]==63&&(int)data[i+1]==72&&(int)data[i+2]==204){//azul
		QC[0][(i/3)%TELX][i/(3*TELX)]=1;
		}
		else if((int)data[i]==255&&(int)data[i+1]==242&&(int)data[i+2]==0){//amarelo else if(*(int*)(data+i)==36&&*(int*)(data+(i+1))==28&&*(int*)(data+(i+2))==237)
		QD[0][(i/3)%TELX][i/(3*TELX)]=1;
		}
		else if((int)data[i]==34&&(int)data[i+1]==177&&(int)data[i+2]==76){//verde
		QE[0][(i/3)%TELX][i/(3*TELX)]=1;
		}*/
		else{
		QC[0][(i/3)%TELX][i/(3*TELX)]=EXISTE;
		}
	}
}
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
	for(int a=80;a<=300;a++)
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
			printf("\n nf=%d n2=%d, a=%d af=%d, theta=%f reg=%f\n",nf,n2,a,af,theta*180/M_PI,reg);
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
void cristal(float **cont, float ySim[])//(FILE *arquivo, float **cont, float ySim[])
{
bool FLAG=true;

int j;
int n1=0;

	//OBTENHO A CURVA SIMULADA
	for(int i=TELX-1;i>=0;i--)
	{
		j=(TELY-1)/2;
		do{
		j=j-1;
		}
		while(cont[i][j]>0.40&&j>0);
		if(j<=((TELX/2)-3)&&FLAG==true){
		n1=i;
		FLAG=false;
		}
		if(j<(TELX/2)-3)
		ySim[-i+n1]=-j+TELY/2;
		//printf("n1=%d i=%d\n",n1,i);
		//printf("ySim[%d]=%d\n",i-n1,-j+TELY/2);
		//fprintf(arquivo,"%d %d\n",i-n1,-j+TELY/2); //-i+TELX -j+TELY/2
	}
	
}

//FUNÇÕES PRINCIPAIS
extern __shared__ float cache[];
	
__device__ void deltaf(float *delta, float *X)
{
	int index = blockIdx.x*blockDim.x + threadIdx.x;
	int stride = blockDim.x*gridDim.x;
	for (int idx = index; idx < TELX*TELY; idx += stride) {
		delta[idx] = X[idx + TELX*TELY] - X[idx];
		}
}

__device__ float maior(float *delta)
{
int index = blockIdx.x*blockDim.x + threadIdx.x;
int stride = blockDim.x*gridDim.x;
int cacheIndex = threadIdx.x;
float m;
for (int idx = index; idx < TELX*TELY; idx += stride) {
	if (delta[idx] > m)
		m = delta[idx];
}
cache[cacheIndex] = m;
__syncthreads();

int ib = blockDim.x / 2;
while (ib != 0) {
	if (cacheIndex<ib&&cache[cacheIndex + ib]>cache[cacheIndex])
		cache[cacheIndex] = cache[cacheIndex + ib];
	__syncthreads();
	ib /= 2;
}
m = cache[0];
return m;
}
//printf("delta[%d]=%f erro=%f",idx,delta[idx],*m);}


//FUNÇÃO PARA ATUALIZAR OS VALORES DOS VETORES NA GPU
__global__ void atualizaTudo(double *Q, double *K)
{
			atualiza(Q, K);
}
__device__ void atualiza(double *Q, double *K)
{
//int idx=blockIdx.x * blockDim.x + threadIdx.x;
//if(idx<TELX*TELY){
	int index = blockIdx.x*blockDim.x + threadIdx.x;
	int stride = blockDim.x*gridDim.x;
	for (int idx = index; idx < TELX*TELY; idx += stride) {
	Q[idx]=K[idx+TELX*TELY];
	}
}
__global__ void atualizaTudofloat(float *Q, float *K)
{
			atualizafloat(Q, K);
}
__device__ void atualizafloat(float *Q, float *K)
{
//int idx=blockIdx.x * blockDim.x + threadIdx.x;
//if(idx<TELX*TELY){
	int index = blockIdx.x*blockDim.x + threadIdx.x;
	int stride = blockDim.x*gridDim.x;
	for (int idx = index; idx < TELX*TELY; idx += stride) {
	Q[idx]=K[idx+TELX*TELY];
	}
}

//FUNÇÃO QUE CALCULA NOVA VARIÁVEL DE FASE na GPU
__global__ void P1(float *Pa, float *Pb, float *Pc, float *uf, double dx, double dy, double LAMBDAAB, double LAMBDAAC,double LAMBDABC, double SIGMAAB, double SIGMAAC,double SIGMABC, double miAB, double miAC,double miBC, double LAB, double LAC,double LBC, double TmAB, double TmAC,double TmBC, int FLAG1, int FLAG2, int FLAG3)
{
	//float M;
	int index = blockIdx.x*blockDim.x+threadIdx.x;
	int stride = blockDim.x*gridDim.x;
	for (int idx = index; idx < TELX*TELY; idx += stride) {

	float rhsf=rhs(Pa,Pb,Pc,uf,idx,dx,dy,LAMBDAAB,LAMBDAAC,LAMBDABC,SIGMAAB,SIGMAAC,SIGMABC,LAB,LAC,LBC,TmAB,TmAC,TmBC,FLAG1,FLAG2,FLAG3,miAB,miAC,miBC);
	//float lambdaf=lambda(Pa,Pb,Pc,uf,idx,dx,dy,epsilonAB,epsilonAC,epsilonBC,WAB,WAC,WBC,L1,Tm1,L2,Tm2,L3,Tm3,FLAG1,FLAG2,FLAG3,miAB,miAC,miBC);
	
	//float Pxa=Px(Pa,idx,dx,FLAG1);
	//float Pya=Py(Pa,idx,dy);

	//float contorno=Pa[idx]*Pa[idx]+Pb[idx]*Pb[idx]+Pc[idx]*Pc[idx];
		
	//if(contorno>0.3067&&contorno<0.3667)
	//M=100.0*MISL;
	//else
	//M=miAB*(Pa[idx]/(1-Pa[idx]))*(Pb[idx]/(1-Pb[idx]))+miAC*(Pa[idx]/(1-Pa[idx]))*(Pc[idx]/(1-Pc[idx]))+miBC*(Pb[idx]/(1-Pb[idx]))*(Pc[idx]/(1-Pc[idx]));
	//M=(miAB*Pa[idx]*Pb[idx])/((1-Pa[idx])*(1-Pb[idx]))+(miAC*Pa[idx]*Pc[idx])/((1-Pa[idx])*(1-Pc[idx]))+(miBC*Pb[idx]*Pc[idx])/((1-Pb[idx])*(1-Pc[idx]));
	
	//float M=miAB*Pa[idx]*Pb[idx]+miAC*Pa[idx]*Pc[idx];
	//P1=P[0][i][j]-M*dt*((pow(P[0][i][j],3.0)-1.5*pow(P[0][i][j],2.0)+0.5*P[0][i][j])+(noise(P,i,j)+(0.9/M_PI)*atan(10*(T[0][i][j]-TM)))*(-pow(P[0][i][j],2.0)+P[0][i][j])-E0*(I+II+III));
	Pa[idx+TELX*TELY]=Pa[idx]+dt*(rhsf);//-(1.0/3.0)*lambdaf);
		
	//P[idx + TELX*TELY]=P[idx]+(dt/tauf)*((P[idx]-lambda*u[idx]*(1.0-powf(P[idx],2.0)))*(1.0-powf(P[idx],2.0))+E0*(I+II+III));
	}
}


__global__ void Temp(float *u, float *X, float *Pa,float *Pb, float *delta, float Fox, float Foy, double L1)
{
float m=0;
//printf(" x");
	do{
	int index = blockIdx.x*blockDim.x+threadIdx.x;
	int stride = blockDim.x*gridDim.x;
	for (int idx = index; idx < TELX*TELY; idx += stride) {	
		X[idx+TELX*TELY]=X1(idx,u, X, Pa,Pb, Fox, Foy, L1);
	}
			__threadfence();
			deltaf(delta,X);
			__threadfence();			
			atualizafloat(X, X);
			__threadfence();
			m=maior(delta);
			__threadfence();
	}while(m>err);
}


//FUNÇÃO QUE CALCULA NOVA TEMPERATURA NA GPU


__device__ float X1(int idx,float *u, float *X, float *Pa,float *Pb, float Fox, float Foy, double L1)
{
float X1;

if(idx>=0&&idx<TELX)//linha de baixo j==0
			{
				if(idx==0) //idx%TELX==0
				//X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx+1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx+1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx+TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[idx+TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(L1)*((Pa[idx+TELX*TELY]-Pa[idx])/(2.0*Fox+2.0*Foy+1));
				X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx+1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(T*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx+TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(T*(Foy/(2.0*Fox+2.0*Foy+1)))+(L1)*((Pa[idx+TELX*TELY]-Pa[idx])/(2.0*Fox+2.0*Foy+1));
				
				else if(idx==TELX-1)//(idx+1)%TELX==0
				//X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx-1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx-1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx+TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[idx+TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(L1)*((Pa[idx+TELX*TELY]-Pa[idx])/(2.0*Fox+2.0*Foy+1));
				X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx-1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(T*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx+TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(T*(Foy/(2.0*Fox+2.0*Foy+1)))+(L1)*((Pa[idx+TELX*TELY]-Pa[idx])/(2.0*Fox+2.0*Foy+1));
				
				else
				//X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx-1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx+1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx+TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[idx+TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(L1)*((Pa[idx+TELX*TELY]-Pa[idx])/(2.0*Fox+2.0*Foy+1));
				X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx-1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx+1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx+TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(T*(Foy/(2.0*Fox+2.0*Foy+1)))+(L1)*((Pa[idx+TELX*TELY]-Pa[idx])/(2.0*Fox+2.0*Foy+1));
			}	

else if(idx>=TELX*(TELY-1)&&idx<TELX*TELY)//linha de cima j==TELY-1
			{
				if(idx==TELX*(TELY-1))
				//X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx+1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx+1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx-TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[idx-TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(L1)*((Pa[idx+TELX*TELY]-Pa[idx])/(2.0*Fox+2.0*Foy+1));
				X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx+1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(T*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx-TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(T*(Foy/(2.0*Fox+2.0*Foy+1)))+(L1)*((Pa[idx+TELX*TELY]-Pa[idx])/(2.0*Fox+2.0*Foy+1));
				
				else if(idx==(TELX*TELY)-1)
				//X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx-1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx-1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx-TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[idx-TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(L1)*((Pa[idx+TELX*TELY]-Pa[idx])/(2.0*Fox+2.0*Foy+1));
				X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx-1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(T*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx-TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(T*(Foy/(2.0*Fox+2.0*Foy+1)))+(L1)*((Pa[idx+TELX*TELY]-Pa[idx])/(2.0*Fox+2.0*Foy+1));
				
				else
				//X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx-1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx+1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx-TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(T*(Foy/(2.0*Fox+2.0*Foy+1)))+(L1)*((Pa[idx+TELX*TELY]-Pa[idx])/(2.0*Fox+2.0*Foy+1));
				X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx-1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx+1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx-TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[idx-TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(L1)*((Pa[idx+TELX*TELY]-Pa[idx])/(2.0*Fox+2.0*Foy+1));
				
			}	
else//(1<=j<=tely-2)	// 1<=y<=tely-2  meio
			{
				if(idx%TELX==0)
				//X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx+1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx+1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx-TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[idx+TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(L1)*((Pa[idx+TELX*TELY]-Pa[idx])/(2.0*Fox+2.0*Foy+1));//0<x<L
				X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx+1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(T*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx-TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[idx+TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(L1)*((Pa[idx+TELX*TELY]-Pa[idx])/(2.0*Fox+2.0*Foy+1));//0<x<L
				
				else if((idx+1)%TELX==0)
				//X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx-1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx-1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx-TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[idx+TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(L1)*((Pa[idx+TELX*TELY]-Pa[idx])/(2.0*Fox+2.0*Foy+1));//x=L
				X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx-1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(T*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx-TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[idx+TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(L1)*((Pa[idx+TELX*TELY]-Pa[idx])/(2.0*Fox+2.0*Foy+1));//x=L
				
				else
				//X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx-1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx+1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx-TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[idx+TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(L1)*((Pa[idx+TELX*TELY]-Pa[idx])/(2.0*Fox+2.0*Foy+1));
				X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx-1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx+1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx-TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[idx+TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(L1)*((Pa[idx+TELX*TELY]-Pa[idx])/(2.0*Fox+2.0*Foy+1));
			}
			
//1D
/*
if(idx%TELX==0)
X1=(u[idx]/(2.0*Fox+1))+(T*(Fox/(2.0*Fox+1)))+(X[idx+1]*(Fox/(2.0*Fox+1)))+(L1)*((Pa[idx+TELX*TELY]-Pa[idx])/(2.0*Fox+1));//-(L2)*6.0*(Pa[idx]-powf(Pa[idx],2))*((Pb[idx+TELX*TELY]-Pb[idx])/(2.0*Fox+1));
else if((idx+1)%TELX==0)
X1=(u[idx]/(2.0*Fox+1))+(X[idx-1]*(Fox/(2.0*Fox+1)))+(X[idx-1]*(Fox/(2.0*Fox+1)))+(L1)*((Pa[idx+TELX*TELY]-Pa[idx])/(2.0*Fox+1));//-(L2)*6.0*(Pa[idx]-powf(Pa[idx],2))*((Pb[idx+TELX*TELY]-Pb[idx])/(2.0*Fox+1));
else
X1=(u[idx]/(2.0*Fox+1))+(X[idx-1]*(Fox/(2.0*Fox+1)))+(X[idx+1]*(Fox/(2.0*Fox+1)))+(L1)*((Pa[idx+TELX*TELY]-Pa[idx])/(2.0*Fox+1));//-(L2)*6.0*(Pa[idx]-powf(Pa[idx],2))*((Pb[idx+TELX*TELY]-Pb[idx])/(2.0*Fox+1));
*/
return X1;
}


//FUNCAO QUE CALCULA RHS NA GPU
__device__ float rhs(float *Pa, float *Pb, float *Pc, float *uf, int idx, double dx, double dy, double LAMBDAAB, double LAMBDAAC,double LAMBDABC, double SIGMAAB,double SIGMAAC,double SIGMABC, double LAB, double LAC,double LBC, double TmAB, double TmAC,double TmBC, int FLAG1, int FLAG2, int FLAG3,double miAB, double miAC, double miBC)
{
float rhs;

//if (FLAG1==2)
//miAB=-miAB;

//float DwDphif=DwDphi(Pa,Pb,Pc,idx,WAB,WAC,WBC,miAB,miAC,miBC);
float frengf=freng(Pa,Pb,Pc,uf,idx,LAMBDAAB,LAMBDAAC,LAMBDABC,SIGMAAB,SIGMAAC,SIGMABC,LAB,LAC,LBC,TmAB,TmAC,TmBC,miAB,miAC,miBC,FLAG1);
float gradengf=gradeng(Pa,Pb,Pc,idx,dx,dy,LAMBDAAB,LAMBDAAC,LAMBDABC,SIGMAAB,SIGMAAC,SIGMABC,miAB,miAC,miBC,FLAG1);

//rhs=gradengf-DwDphif-frengf;
rhs=gradengf+frengf;

//rhs=EPSILON*gradengf-(1.0/EPSILON)*DwDphif-(1.0/uf[idx])*frengf;
//rhs=EPSILON*(-divDaDNphif-DaDphif)-(1.0/EPSILON)*DwDphif-(1.0/uf[idx])*frengf;
return rhs;
}
/*__device__ float MAB(float *Pa, float *Pb, int idx, double miAB)
{
float M;
//M=miAB;
M=miAB*(Pa[idx]/(1-Pa[idx]))*(Pb[idx]/(1-Pb[idx]));
//M=(miAB*Pa[idx]*Pb[idx])/((1-Pa[idx])*(1-Pb[idx]))*(miAC*Pa[idx]*Pc[idx])/((1-Pa[idx])*(1-Pc[idx]));
return M;
}
__device__ float MAC(float *Pa, float *Pc, int idx, double miAC)
{
float M;
//M=miAC;
M=miAC*(Pa[idx]/(1-Pa[idx]))*(Pc[idx]/(1-Pc[idx]));
//M=(miAB*Pa[idx]*Pb[idx])/((1-Pa[idx])*(1-Pb[idx]))*(miAC*Pa[idx]*Pc[idx])/((1-Pa[idx])*(1-Pc[idx]));
return M;
}
__device__ float MBC(float *Pb, float *Pc, int idx, double miBC)
{
float M;
//M=miBC;
M=miBC*(Pb[idx]/(1-Pb[idx]))*(Pc[idx]/(1-Pc[idx]));
//M=(miAB*Pa[idx]*Pb[idx])/((1-Pa[idx])*(1-Pb[idx]))*(miAC*Pa[idx]*Pc[idx])/((1-Pa[idx])*(1-Pc[idx]));
return M;
}*/
//FUNÇÃO QUE CALCULA A ENERGIA DO GRADIENTE
__device__ float gradeng(float *Pa, float *Pb, float *Pc, int idx, double dx, double dy, double LAMBDAAB, double LAMBDAAC,double LAMBDABC, double SIGMAAB,double SIGMAAC,double SIGMABC, double miAB, double miAC,double miBC, int FLAG)
{
float gradeng;

float epsilonAB=LAMBDAAB*SIGMAAB;
float epsilonAC=LAMBDAAC*SIGMAAC;
float epsilonBC=LAMBDABC*SIGMABC;

//float Pxa=Px(Pa,idx,dx,FLAG1);
//float Pxb=Px(Pb,idx,dx,FLAG2);
//float Pxc=Px(Pc,idx,dx,FLAG3);
//float Pya=Py(Pa,idx,dy);
//float Pyb=Py(Pb,idx,dy);
//float Pyc=Py(Pc,idx,dy);

//float Pxxa=Pxx(Pa,idx,dx,FLAG1);
//float Pxxb=Pxx(Pb,idx,dx,FLAG2);
//float Pxxc=Pxx(Pc,idx,dx,FLAG3);

float N2A=Pxx(Pa,idx,dx)+Pyy(Pa,idx,dy);
float N2B=Pxx(Pb,idx,dx)+Pyy(Pb,idx,dy);
float N2C=Pxx(Pc,idx,dx)+Pyy(Pc,idx,dy);

if ((Pa[idx]*Pa[idx]+Pb[idx]*Pb[idx]+Pc[idx]*Pc[idx])<0.5)
miAB=miAC=miBC=MISL;
//gradeng=2.0*(epsilonAB)*(2.0*(Pa[idx]*(Pxb*Pxb+Pyb*Pyb)-Pb[idx]*(Pxa*Pxb+Pya*Pyb))+Pa[idx]*Pb[idx]*N2B-powf(Pb[idx],2)*N2A)+2.0*(epsilonAC)*(2.0*(Pa[idx]*(Pxc*Pxc+Pyc*Pyc)-Pc[idx]*(Pxa*Pxc+Pya*Pyc))+Pa[idx]*Pc[idx]*N2C-powf(Pc[idx],2));

//gradeng=MABf* epsilonAB*(N2B-N2A)+MACf* epsilonAC*(N2C-N2A)-MBCf* epsilonBC*(N2C+N2B);
gradeng= miAB* epsilonAB*(Pa[idx]*N2B-Pb[idx]*N2A)+ miAC* epsilonAC*(Pa[idx]*N2C-Pc[idx]*N2A)+Pa[idx]*Pb[idx]*Pc[idx];

//1D
//gradeng=2.0*(gammaConstAB/constmAB)*(2.0*(Pa[idx]*Pxb*Pxb-Pb[idx]*Pxa*Pxb)+Pa[idx]*Pb[idx]*Pxxb-powf(Pb[idx],2)*Pxxa)+2.0*(gammaConstAC/constmAC)*(2.0*(Pa[idx]*Pxc*Pxc-Pc[idx]*Pxa*Pxc)+Pa[idx]*Pc[idx]*Pxxc-powf(Pc[idx],2)*Pxxa);
//gradeng=epsilonAB*(Pxxb-Pxxa)+epsilonAC*(Pxxc-Pxxa)-epsilonBC*(Pxxc+Pxxb);

return -gradeng;

}

//FUNÇÃO QUE CALCULA DwDphi NA GPU
/*__device__ float DwDphi(float *Pa, float *Pb,float *Pc, int idx, float WAB,float WAC,float WBC,float MABf, float MACf, float MBCf)
{
float DwDphi;

//DwDphi=9.0*(2.0*gammaConstAB*constmAB*Pa[idx]*powf(Pb[idx],2)+2.0*gammaConstAC*constmAC*Pa[idx]*powf(Pc[idx],2))+72*(2.0*gammaConstAB*constmAB*Pa[idx]*powf(Pb[idx],2)*Pc[idx]+2.0*gammaConstAC*constmAC*Pa[idx]*powf(Pc[idx],2)*Pb[idx]);
//DwDphi=(9.0/EPSILON)*(2.0*WAB*Pa[idx]*powf(Pb[idx],2)+2.0*WAC*Pa[idx]*powf(Pc[idx],2))+(72.0/EPSILON)*(2.0*WAB*Pa[idx]*powf(Pb[idx],2)*Pc[idx]+2.0*WAC*Pa[idx]*powf(Pc[idx],2)*Pb[idx]);
//DwDphi=2.0*WAB*(Pa[idx]*powf(Pb[idx],2)-powf(Pa[idx],2)*Pb[idx])+2.0*WAC*(Pa[idx]*powf(Pc[idx],2)-powf(Pa[idx],2)*Pc[idx])-2.0*WBC*(Pb[idx]*powf(Pc[idx],2)+powf(Pb[idx],2)*Pc[idx]);

//DwDphi=MABf* WAB*(Pb[idx]-Pa[idx])+MACf* WAC*(Pc[idx]-Pa[idx])-MBCf* WBC*(Pc[idx]+Pb[idx]);
DwDphi=MABf*  WAB*Pa[idx]*Pb[idx]*(Pb[idx]-Pa[idx])+  MACf*   WAC*Pa[idx]*Pc[idx]*(Pc[idx]-Pa[idx])+MBCf*   WBC*Pb[idx]*Pc[idx]*(Pc[idx]-Pb[idx]);

return DwDphi;
}*/
//FUNÇÃO QUE CALCULA h`f() Por enquanto incluida dentro de rhs
__device__ float freng(float *Pa, float *Pb, float *Pc, float *uf, int idx,double LAMBDAAB, double LAMBDAAC,double LAMBDABC, double SIGMAAB,double SIGMAAC,double SIGMABC, double LAB, double LAC,double LBC, double TmAB, double TmAC,double TmBC, double miAB, double miAC,double miBC,int FLAG)
{
float freng;
float mAB;
float aAB;
float mAC;
float aAC;
float mBC;
float aBC;
float lab;
float lac;
float lbc;
if (FLAG==2||FLAG==1){
lab=-LAB;
lac=-LAC;
lbc=LBC;
}
else{
lab=LAB;
lac=LAC;
lbc=LBC;
}
aAB=LAMBDAAB/(72.0*SIGMAAB);
mAB=(6.0*aAB*lab*(TmAB-uf[idx]))/TmAB;

aAC=LAMBDAAC/(72.0*SIGMAAC);
mAC=(6.0*aAC*lac*(TmAC-uf[idx]))/TmAC;

//aBC=LAMBDABC/(72.0*SIGMABC);
//mBC=(6.0*aBC*lbc*(TmBC-uf[idx]))/TmBC;
if ((Pa[idx]*Pa[idx]+Pb[idx]*Pb[idx]+Pc[idx]*Pc[idx])<0.5)
miAB=miAC=miBC=MISL;

freng=-miAB*((Pa[idx]*Pb[idx])/(2.0*aAB))*(Pb[idx]-Pa[idx]-2.0*mAB*(TmAB-uf[idx]))-miAC*((Pa[idx]*Pc[idx])/(2.0*aAC))*(Pc[idx]-Pa[idx]-2.0*mAC*(TmAC-uf[idx]));//+miBC*((Pb[idx]*Pc[idx])/(2.0*aBC))*(Pc[idx]+Pb[idx]-2.0*mBC*(TmBC-uf[idx]));
//else
//freng=MISL*(L*((uf[idx]-Tm)/(Tm))*6.0*(Pa[idx]-powf(Pa[idx],2)));//+MACf*(L*((uf[idx]-Tm)/(Tm))*6.0*(Pa[idx]-powf(Pa[idx],2)));
//o correto seria um calor latente de formação do sólido a partir de cada uma das fases, mas como não se forma a aprtir de substrato, ta ok.

return freng;
}

//FUNÇÃO QUE CALCULA LAMBDA NA GPU
/*__device__ float lambda(float *Pa, float *Pb, float *Pc, float *uf, int idx, double dx, double dy, double epsilonAB, double epsilonAC, double epsilonBC, float WAB,float WAC,float WBC, float L1, float Tm1, float L2, float Tm2, float L3, float Tm3, int FLAG1, int FLAG2, int FLAG3,double miAB, double miAC, double miBC)
{
float lambda;

lambda=rhs(Pa,Pb,Pc,uf,idx,dx,dy,epsilonAB,epsilonAC,epsilonBC,WAB,WAC,WBC,L1,Tm1,FLAG1,FLAG2,FLAG3,miAB,miAC,miBC)+rhs(Pb,Pa,Pc,uf,idx,dx,dy,epsilonAB,epsilonBC,epsilonAC,WAB,WBC,WAC,L2,Tm2,FLAG2,FLAG1,FLAG3,miAB,miBC,miAC)+rhs(Pc,Pb,Pa,uf,idx,dx,dy,epsilonBC,epsilonAC,epsilonAB,WBC,WAC,WAB,L3,Tm3,FLAG3,FLAG2,FLAG1,miBC,miAC,miAB);
return lambda;
}*/
//FUNÇÕES QUE CALCULAM AS DERIVADAS DE P NA GPU
__device__ float Px(float *P, int idx, double dx)
{
float Px;

if(idx%TELX==0)
Px=(P[idx+1]-P[idx+1])/(2.0*dx);//Px=(P[idx+1]-boundE)/(2.0*dx);
else if((idx+1)%TELX==0)
Px=(P[idx-1]-P[idx-1])/(2.0*dx);
else
Px=(P[idx+1]-P[idx-1])/(2.0*dx);
return Px;
}

__device__ float Py(float *P, int idx, double dy)
{
float Py;
if(idx>=0&&idx<TELX)
Py=(P[idx+TELX]-P[idx+TELX])/(2.0*dy);
else if(idx>=TELX*(TELY-1)&&idx<TELX*TELY)
Py=(P[idx-TELX]-P[idx-TELX])/(2.0*dy);
else
Py=(P[idx+TELX]-P[idx-TELX])/(2.0*dy);
return Py;
}


__device__ float Pxx(float *P, int idx, double dx)
{
float Pxx;

if(idx%TELX==0)
Pxx=(P[idx+1]-2.0*P[idx]+P[idx+1])/(dx*dx);//Pxx=(P[idx+1]-2.0*P[idx]+boundE)/(dx*dx);Pxx=(P[idx+1]-2.0*P[idx]+P[idx+1])/(dx*dx);
else if((idx+1)%TELX==0)
Pxx=(P[idx-1]-2.0*P[idx]+P[idx-1])/(dx*dx);
else
Pxx=(P[idx+1]-2.0*P[idx]+P[idx-1])/(dx*dx);
return Pxx;
}
__device__ float Pyy(float *P, int idx, double dy)
{
float Pyy;

if(idx>=0&&idx<TELX)
Pyy=(P[idx+TELX]-2.0*P[idx]+P[idx+TELX])/(dy*dy);
else if(idx>=TELX*(TELY-1)&&idx<TELX*TELY)
Pyy=(P[idx-TELX]-2.0*P[idx]+P[idx-TELX])/(dy*dy);
else
Pyy=(P[idx+TELX]-2.0*P[idx]+P[idx-TELX])/(dy*dy);
return Pyy;
}

__device__ float Pxy(float *P, int idx, double dx, double dy)
{
float Pxy;
/*
float bound;
if(FLAG==0)
bound=1.0;
else
bound=0.0;
*/
if (idx>=0&&idx<TELX){
		if (idx==0)
			Pxy=(P[idx+TELX+1]-P[idx+TELX+1]-P[idx+TELX+1]+P[idx+TELX+1])/(4.0*dx*dy);
		else if (idx==TELX-1)
			Pxy=(P[idx-1+TELX]-P[idx-1+TELX]-P[idx-1+TELX]+P[idx-1+TELX])/(4.0*dx*dy);
		else
			Pxy=(P[idx+1+TELX]-P[idx-1+TELX]-P[idx+1+TELX]+P[idx-1+TELX])/(4.0*dx*dy);
	}
	else if (idx>=TELX*(TELY-1)&&idx<TELX*TELY){
		if (idx==TELX*(TELY-1))
			Pxy=(P[idx+1-TELX]-P[idx+1-TELX]-P[idx+1-TELX]+P[idx+1-TELX])/(4.0*dx*dy);
		else if (idx==(TELX*TELY)-1)
			Pxy=(P[idx-1-TELX]-P[idx-1-TELX]-P[idx-1-TELX]+P[idx-1-TELX])/(4.0*dx*dy);
		else
			Pxy=(P[idx+1-TELX]-P[idx-1-TELX]-P[idx+1-TELX]+P[idx-1-TELX])/(4.0*dx*dy);
	}
	else {
		if (idx%TELX==0)
			Pxy=(P[idx+1+TELX]-P[idx+1+TELX]-P[idx+1-TELX]+P[idx+1-TELX])/(4.0*dx*dy);
		else if ((idx+1)%TELX==0)
			Pxy=(P[idx-1+TELX]-P[idx-1+TELX]-P[idx-1-TELX]+P[idx-1-TELX])/(4.0*dx*dy);
		else
			Pxy=(P[idx+1+TELX]-P[idx-1+TELX]-P[idx+1-TELX]+P[idx-1-TELX])/(4.0*dx*dy);
	}
	return Pxy;
}


__device__ float Pyx(float *P, int idx, double dx, double dy)
{
	float Pxy;
	/*if(idx==0) //inferior esquerdo
	Pxy=(P[idx+1+TELX]-P[idx+1+TELX]-P[idx+1+TELX]+P[idx+1+TELX])/(4.0*dx*dy);	
	else if(idx==TELX-1)//inferior direito
	Pxy=(P[idx-1+TELX]-P[idx-1+TELX]-P[idx-1+TELX]+P[idx-1+TELX])/(4.0*dx*dy);	
	else if(idx==TELX*(TELY-1))//superior esquerdo
	Pxy=(P[idx+1-TELX]-P[idx+1-TELX]-P[idx+1-TELX]+P[idx+1-TELX])/(4.0*dx*dy);
	else if(idx==(TELX*TELY)-1)//superior direito
	Pxy=(P[idx-1-TELX]-P[idx-1-TELX]-P[idx-1-TELX]+P[idx-1-TELX])/(4.0*dx*dy);	
	else if(idx>0&&idx<(TELX-1))//contorno de baixo
	Pxy=(P[idx+1+TELX]-P[idx-1+TELX]-P[idx+1+TELX]+P[idx-1+TELX])/(4.0*dx*dy);
	else if(idx>TELX*(TELY-1)&&idx<((TELX*TELY)-1))//contorno de cima
	Pxy=(P[idx+1-TELX]-P[idx-1-TELX]-P[idx+1-TELX]+P[idx-1-TELX])/(4.0*dx*dy);	
	else if(idx%TELX==0)//contorno esquerdo
	Pxy=(P[idx+1+TELX]-P[idx+1+TELX]-P[idx+1-TELX]+P[idx+1-TELX])/(4.0*dx*dy);	
	else if((idx+1)%TELX==0)//contorno direito
	Pxy=(P[idx-1+TELX]-P[idx-1+TELX]-P[idx-1-TELX]+P[idx-1-TELX])/(4.0*dx*dy);	
	else//meio
	Pxy=(P[idx+1+TELX]-P[idx-1+TELX]-P[idx+1-TELX]+P[idx-1-TELX])/(4.0*dx*dy);
	*/
	/*
float bound;
if(FLAG==0)
bound=1.0;
else
bound=0.0;
*/
if (idx>=0&&idx<TELX){
		if (idx==0)
			Pxy=(P[idx+TELX+1]-P[idx+TELX+1]-P[idx+TELX+1]+P[idx+TELX+1])/(4.0*dx*dy);
		else if (idx==TELX-1)
			Pxy=(P[idx-1+TELX]-P[idx-1+TELX]-P[idx-1+TELX]+P[idx-1+TELX])/(4.0*dx*dy);
		else
			Pxy=(P[idx+1+TELX]-P[idx-1+TELX]-P[idx+1+TELX]+P[idx-1+TELX])/(4.0*dx*dy);
	}
	else if (idx>=TELX*(TELY-1)&&idx<TELX*TELY){
		if (idx==TELX*(TELY-1))
			Pxy=(P[idx+1-TELX]-P[idx+1-TELX]-P[idx+1-TELX]+P[idx+1-TELX])/(4.0*dx*dy);
		else if (idx==(TELX*TELY)-1)
			Pxy=(P[idx-1-TELX]-P[idx-1-TELX]-P[idx-1-TELX]+P[idx-1-TELX])/(4.0*dx*dy);
		else
			Pxy=(P[idx+1-TELX]-P[idx-1-TELX]-P[idx+1-TELX]+P[idx-1-TELX])/(4.0*dx*dy);
	}
	else {
		if (idx%TELX==0)
			Pxy=(P[idx+1+TELX]-P[idx+1+TELX]-P[idx+1-TELX]+P[idx+1-TELX])/(4.0*dx*dy);
		else if ((idx+1)%TELX==0)
			Pxy=(P[idx-1+TELX]-P[idx-1+TELX]-P[idx-1-TELX]+P[idx-1-TELX])/(4.0*dx*dy);
		else
			Pxy=(P[idx+1+TELX]-P[idx-1+TELX]-P[idx+1-TELX]+P[idx-1-TELX])/(4.0*dx*dy);
	}
	return Pxy;
}
/*
__device__ float Epx(float *P,int idx, float dx, float dy)
{
float Epx;
float Ef=E(P,idx,dx,dy);
float Pxf=Px(P,idx,dx);
float Pyf=Py(P,idx,dy);

if((Pxf<0.00000001&&Pxf>-0.00000001)&&(Pyf<0.00000001&&Pyf>-0.00000001))
	Epx=0;
	else
	Epx=16.0*S*((powf(Pxf,3)*powf(Pyf,2)-Pxf*powf(Pyf,4))/powf((powf(Pxf,2)+powf(Pyf,2)),3));
return Epx;
}



__device__ float Epy(float *P, int idx, float dx, float dy)//verificar se preciso declarar Pxf e Pyf
{
float Epy;
float Ef=E(P,idx,dx,dy);
float Pxf=Px(P,idx,dx);
float Pyf=Py(P,idx,dy);
if((Pxf<0.00000001&&Pxf>-0.00000001)&&(Pyf<0.00000001&&Pyf>-0.00000001))
	Epy=0;
	else
	Epy=16.0*S*((powf(Pyf,3)*powf(Pxf,2)-Pyf*powf(Pxf,4))/powf((powf(Pxf,2)+powf(Pyf,2)),3));
return Epy;
}

__device__ float Ex(float *P, int idx, float dx, float dy)
{
float Ex;
float Epxf=Epx(P,idx,dx,dy);
float Epyf=Epy(P,idx,dx,dy);
float Pxxf=Pxx(P,idx,dx);
float Pxyf=Pxy(P,idx,dx,dy);

Ex=Epxf*Pxxf+Epyf*Pxyf;
return Ex;
}


__device__ float Ey(float *P, int idx, float dx, float dy)
{
float Ey;
float Epxf=Epx(P,idx,dx,dy);
float Epyf=Epy(P,idx,dx,dy);
float Pyyf=Pyy(P,idx,dy);
float Pyxf=Pyx(P,idx,dx,dy);

Ey=Epyf*Pyyf+Epxf*Pyxf;
return Ey;
}

__device__ float Epxpx(float *P, int idx, float dx, float dy)
{
float Epxpx;
float Epxf=Epx(P,idx,dx,dy);
float Pxf=Px(P,idx,dx);
float Pyf=Py(P,idx,dy);
if((Pxf<0.00000001&&Pxf>-0.00000001)&&(Pyf<0.00000001&&Pyf>-0.00000001))
	Epxpx=0;
	else
	Epxpx=((16.0*S)/powf((powf(Pxf,2)+powf(Pyf,2)),3))*(3.0*powf(Pxf,2)*powf(Pyf,2)-powf(Pyf,4))-((6.0*Pxf)/(powf(Pxf,2)+powf(Pyf,2)))*Epxf;
return Epxpx;
}

__device__ float Epypy(float *P, int idx, float dx, float dy)
{
float Epypy;
float Epyf=Epy(P,idx,dx,dy);
float Pxf=Px(P,idx,dx);
float Pyf=Py(P,idx,dy);
if((Pxf<0.00000001&&Pxf>-0.00000001)&&(Pyf<0.00000001&&Pyf>-0.00000001))
	Epypy=0;
	else
	Epypy=((16.0*S)/powf((powf(Pxf,2)+powf(Pyf,2)),3))*(3.0*powf(Pyf,2)*powf(Pxf,2)-powf(Pxf,4))-((6.0*Pyf)/(powf(Pxf,2)+powf(Pyf,2)))*Epyf;
return Epypy;
}


__device__ float Epxpy(float *P, int idx, float dx, float dy)
{
float Epxpy;
float Epyf=Epy(P,idx,dx,dy);
float Pxf=Px(P,idx,dx);
float Pyf=Py(P,idx,dy);

if((Pxf<0.00000001&&Pxf>-0.00000001)&&(Pyf<0.00000001&&Pyf>-0.00000001))
	Epxpy=0;
	else
	Epxpy=((16.0*S)/powf((powf(Pxf,2)+powf(Pyf,2)),3))*(2.0*powf(Pyf,3)*Pxf-4.0*Pyf*powf(Pxf,3))-((6.0*Pxf)/(powf(Pxf,2)+powf(Pyf,2)))*Epyf;	

return Epxpy;
}

__device__ float Epypx(float *P,int idx, float dx, float dy)
{
float Epypx;
float Epxf=Epx(P,idx,dx,dy);
float Pxf=Px(P,idx,dx);
float Pyf=Py(P,idx,dy);

if((Pxf<0.00000001&&Pxf>-0.00000001)&&(Pyf<0.00000001&&Pyf>-0.00000001))
	Epypx=0;
	else
	Epypx=((16.0*S)/powf((powf(Pxf,2)+powf(Pyf,2)),3))*(2.0*powf(Pxf,3)*Pyf-4.0*Pxf*powf(Pyf,3))-((6.0*Pyf)/(powf(Pxf,2)+powf(Pyf,2)))*Epxf;	

return Epypx;
}

__device__ float Expx(float *P,int idx,float dx, float dy)
{
float Expx;
float Epxpxf=Epxpx(P,idx,dx,dy);
float Epypxf=Epypx(P,idx,dx,dy);
float Pxxf=Pxx(P,idx,dx);
float Pxyf=Pxy(P,idx,dx,dy);

Expx=Epxpxf*Pxxf+Epypxf*Pxyf;
return Expx;
}


__device__ float Eypy(float *P, int idx, float dx, float dy)
{
float Eypy;
float Epypyf=Epypy(P,idx,dx,dy);
float Epxpyf=Epxpy(P,idx,dx,dy);
float Pyyf=Pyy(P,idx,dy);
float Pyxf=Pyx(P,idx,dx,dy);

Eypy=Epypyf*Pyyf+Epxpyf*Pyxf;
return Eypy;
}
*/
