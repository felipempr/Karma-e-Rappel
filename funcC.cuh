//__global__ void testef(float *teste);

__device__ float maior(float *delta);
__device__ void deltaf(float *delta, float *X);
__device__ void atualiza(float *Q, float *K);
__device__ void atualizaDouble(float *Q, double *K);
__global__ void atualizaTudo(float *Q, float *K);
__global__ void atualizaTudoDouble(float *Q, double *K);
__global__ void P1(bool *FLAG,float *P, float *u, float dx, float dy);
__global__ void Temp(bool *FLAG, float *u, float *X, float *P,float *delta, float Fox, float Foy);
__device__ float X1(bool *FLAG, int idx,float *u, float *X, float *P, float Fox, float Foy);
//__device__ float tau(bool *FLAG,float *P, int idx, float dx, float dy);
__global__ void raio(float *k,bool *FLAG,float *P, float dx, float dy);
__device__ float E(bool *FLAG,float *P, int idx, float dx, float dy);
__device__ float Px(bool *FLAG,float *P, int idx, float dx);
__device__ float Py(bool *FLAG,float *P, int idx, float dy);
__device__ float Pxx(bool *FLAG,float *P, int idx, float dx);
__device__ float Pyy(bool *FLAG,float *P, int idx, float dy);
__device__ float Pxy(bool *FLAG,float *P, int idx, float dx, float dy);
__device__ float Pyx(bool *FLAG,float *P, int idx, float dx, float dy);
__device__ float Epx(bool *FLAG,float *P, int idx, float dx, float dy);
__device__ float Epy(bool *FLAG,float *P, int idx, float dx, float dy);
__device__ float Ex(bool *FLAG,float *P, int idx, float dx, float dy);
__device__ float Ey(bool *FLAG,float *P, int idx, float dx, float dy);
__device__ float Epxpx(bool *FLAG,float *P, int idx, float dx, float dy);
__device__ float Epypy(bool *FLAG,float *P, int idx, float dx, float dy);
__device__ float Epxpy(bool *FLAG,float *P, int idx, float dx, float dy);
__device__ float Epypx(bool *FLAG,float *P, int idx, float dx, float dy);
__device__ float Expx(bool *FLAG,float *P,int idx,float dx, float dy);
__device__ float Eypy(bool *FLAG,float *P,int idx,float dx, float dy);

void readBMP(unsigned char* data);
void cond_inicial(float ***QA, bool **FLAG);//,float ***QG, float ***QD, float ***QE);

void readBMP(unsigned char* data)
{
    int v;
    FILE* f = fopen("C:\\Users\\Felipe Ribeiro\\Pictures\\karma400.bmp", "rb");

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

void cond_inicial(float ***QA, bool **FLAG)//,float ***QG, float ***QD, float ***QE)

/*{
	for(int i=0;i<=3*TELX*TELY-3;i=i+3){
		//printf("%d ",i);	
		if((int)data[i]==255&&(int)data[i+1]==255&&(int)data[i+2]==255){//branco
		//printf("i=%d j=%d\n",(i/3)%TELX,i/(3*TELX));
		QA[0][(i/3)%TELX][i/(3*TELX)]=SOLIDO;
		FLAG[(i/3)%TELX][i/(3*TELX)]=true;//u[0][(i/3)%TELX][i/(3*TELX)]=T;
		}
		else if((int)data[i]==0&&(int)data[i+1]==0&&(int)data[i+2]==0){//preto  BGR
		QA[0][(i/3)%TELX][i/(3*TELX)]=LIQUIDO;
		FLAG[(i/3)%TELX][i/(3*TELX)]=true;
		}
		/*else if((int)data[i]==63&&(int)data[i+1]==72&&(int)data[i+2]==204){//azul
		QC[0][(i/3)%TELX][i/(3*TELX)]=1;
		}
		else if((int)data[i]==255&&(int)data[i+1]==242&&(int)data[i+2]==0){//amarelo else if(*(int*)(data+i)==36&&*(int*)(data+(i+1))==28&&*(int*)(data+(i+2))==237)
		QD[0][(i/3)%TELX][i/(3*TELX)]=1;
		}
		else if((int)data[i]==34&&(int)data[i+1]==177&&(int)data[i+2]==76){//verde
		QE[0][(i/3)%TELX][i/(3*TELX)]=1;
		}
		else{
		QA[0][(i/3)%TELX][i/(3*TELX)]=LIQUIDO;
		FLAG[(i/3)%TELX][i/(3*TELX)]=false;
		}
	}
}*/
{
for(int i=0;i<TELX;i++){
	for(int j=0;j<TELX;j++){
	if((pow((i-TELX/2),2)+pow((j-TELY/2),2))<=(R*R))
	QA[0][i][j]=SOLIDO;
	}
}
/*for(int i=0;i<TELX;i++){
	for(int j=0;j<TELX;j++){
	if((pow((i-0),2)+pow((j-0),2))<=(R*R))
	QA[0][i][j]=SOLIDO;
	}
}*/
}


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
__global__ void atualizaTudo(float *Q, float *K)
{
			atualiza(Q, K);
}

__global__ void atualizaTudoDouble(float *Q, double *K)
{
			atualizaDouble(Q, K);
}

__device__ void atualiza(float *Q, float *K)
{
//int idx=blockIdx.x * blockDim.x + threadIdx.x;
//if(idx<TELX*TELY){
	int index = blockIdx.x*blockDim.x + threadIdx.x;
	int stride = blockDim.x*gridDim.x;
	for (int idx = index; idx < TELX*TELY; idx += stride) {
	Q[idx]=K[idx+TELX*TELY];
}
	//printf(" atu.");
}
__device__ void atualizaDouble(float *Q, double *K)
{
//int idx=blockIdx.x * blockDim.x + threadIdx.x;
//if(idx<TELX*TELY){
	int index = blockIdx.x*blockDim.x + threadIdx.x;
	int stride = blockDim.x*gridDim.x;
	for (int idx = index; idx < TELX*TELY; idx += stride) {
	Q[idx]=K[idx+TELX*TELY];
}
	//printf(" atu.");
}

//FUNÇÃO QUE CALCULA NOVA VARIÁVEL DE FASE na GPU
__global__ void P1(bool *FLAG, float *P, float *u, float dx, float dy)
{
	//float I;
	//float II;
	//float III;

	//float Ef;
	//float Exf;
	//float Eyf;
	//float Epxf;
	//float Epyf;
	//float Expxf;
	//float Eypyf;
	//float Pxf;
	//float Pyf;
	float Pxxf;
	float Pyyf;
	//float Pxyf;
	//float tauf;

	int index = blockIdx.x*blockDim.x+threadIdx.x;
	int stride = blockDim.x*gridDim.x;
	for (int idx = index; idx < TELX*TELY; idx += stride) {
	/*
	if(FLAG[idx]==false){
	P[idx+TELX*TELY]=P[idx];
	}
	else{*/
	if(P[idx-1]==P[idx+1]&&P[idx+TELX]==P[idx-TELX])	
	P[idx+TELX*TELY]=P[idx];
	else{
	//int idx=blockIdx.x * blockDim.x + threadIdx.x;
	//if(idx<TELX*TELY){
	

		//Ef = E(FLAG, P, idx, dx, dy);
		//Exf = 0.0;//Ex(FLAG,P, idx, dx, dy);
		//Eyf = 0.0;//Ey(FLAG,P, idx, dx, dy);
		//Epxf = 0.0;//Epx(FLAG,P, idx, dx, dy);
		//Epyf = 0.0;//Epy(FLAG,P, idx, dx, dy);
		//Expxf = 0.0;//Expx(FLAG,P, idx, dx, dy);
		//Eypyf = 0.0;//Eypy(FLAG,P, idx, dx, dy);
		//Pxf = Px(FLAG,P, idx, dx);
		//Pyf = Py(FLAG,P, idx, dy);
		Pxxf = Pxx(FLAG,P, idx, dx);
		Pyyf = Pyy(FLAG,P, idx, dy);
		//Pxyf = Pxy(FLAG,P, idx, dx, dy);
		//tauf = TAU0;//tau(FLAG,P, idx, dx, dy);

		//I=powf(Ef, 2)*(Pxxf + Pyyf)+2.0*Ef*(Pxf*Exf+Pyf*Eyf);
		//II=2.0*Ef*Epxf*(Pxf*Pxxf+Pyf*Pxyf)+(powf(Pxf,2)+powf(Pyf,2))*(Epxf*Exf+Ef*Expxf);
		//III=2.0*Ef*Epyf*(Pyf*Pyyf+Pxf*Pxyf)+(powf(Pxf,2)+powf(Pyf,2))*(Epyf*Eyf+Ef*Eypyf);

		//P1=P[0][i][j]-M*dt*((pow(P[0][i][j],3.0)-1.5*pow(P[0][i][j],2.0)+0.5*P[0][i][j])+(noise(P,i,j)+(0.9/M_PI)*atan(10*(T[0][i][j]-TM)))*(-pow(P[0][i][j],2.0)+P[0][i][j])-E0*(I+II+III));

		//P[idx+TELX*TELY]=P[idx]+(dt/tauf)*((P[idx]-lambda*u[idx]*(1.0-powf(P[idx],2.0)))*(1.0-powf(P[idx],2.0))+(E0*E0)*(I+II+III));
		P[idx+TELX*TELY]=P[idx]+(dt*M0)*((E0*E0)*((Pxxf + Pyyf))-a1*E0*((u[idx]-Ueq))*(1.0-powf(P[idx],2.0))*(1.0-powf(P[idx],2.0))+P[idx]*(1.0-powf(P[idx],2.0)));
	//P[idx + TELX*TELY]=P[idx]+(dt/tauf)*((P[idx]-lambda*u[idx]*(1.0-powf(P[idx],2.0)))*(1.0-powf(P[idx],2.0))+E0*(I+II+III));

		//printf("P[5][5]=%f P[30][30]=%f\n", P[5 + 5 * TELX], P[30 + 30 * TELX]);
		//printf("P1[5][5]=%f P1[30][30]=%f\n", P[5 + 5 * TELX + TELX*TELY], P[30 + 30 * TELX + TELX*TELY]);
	//}
	}
	}
}
__global__ void raio(float *k, bool *FLAG,float *P, float dx, float dy)
{
	int index = blockIdx.x*blockDim.x+threadIdx.x;
	int stride = blockDim.x*gridDim.x;
	for (int idx = index; idx < TELX*TELY; idx += stride) {
if(FLAG[idx]==false)
k[idx]=0;
else{
float Pxf;
float Pyf;
float Pxxf;
float Pyyf;
float Pxyf;
//float theta;
Pxf = Px(FLAG,P, idx, dx);
Pyf = Py(FLAG,P, idx, dy);
Pxxf = Pxx(FLAG,P, idx, dx);
Pyyf = Pyy(FLAG,P, idx, dy);
Pxyf = Pxy(FLAG,P, idx, dx, dy);
float condicao=powf(abs(Pxf+Pyf),2);

//theta=atan2(Pyf,Pxf);
if((Pxf+Pyf)==0||condicao<powf(10.0,-10.0))
k[idx]=0;
else
//k[idx]=theta;
//k[idx]=((1.0/(2.0*fabs(Pxf+Pyf)))*((Pxxf+Pyyf)+cos(2.0*theta)*(Pyyf-Pxxf)-2.0*sin(2*theta)*Pxyf));
k[idx]=1.0/(-0.5*((Pxxf*powf(Pyf,2) + Pyyf*powf(Pxf,2) - 2.0*Pxf*Pyf*Pxyf)  /(powf((powf(Pxf,2)+powf(Pyf,2)),1.5)) ));
}
}
}






__global__ void Temp(bool *FLAG, float *u, float *X, float *P,float *delta, float Fox, float Foy)
{
float m=0;
//printf(" x");
	do{
	int index = blockIdx.x*blockDim.x+threadIdx.x;
	int stride = blockDim.x*gridDim.x;
	for (int idx = index; idx < TELX*TELY; idx += stride) {	

	/*if(FLAG[idx]==false){
	X[idx+TELX*TELY]=X[idx];
	}*/
	/*if(X[idx-1]==X[idx+1]&&X[idx+TELX]==X[idx-TELX])	//A NOITE TIRAR ESTE TRECHO
	X[idx+TELX*TELY]=X[idx];
	else{*/
		X[idx+TELX*TELY]=X1(FLAG,idx,u, X, P, Fox, Foy);
	}
	
			__threadfence();
			deltaf(delta,X);
			__threadfence();			
			atualiza(X, X);
			__threadfence();
			m=maior(delta);
			__threadfence();
	}while(m>err);
//printf(" m=%f", m);
}

//FUNÇÃO QUE CALCULA NOVA TEMPERATURA NA GPU

__device__ float X1(bool *FLAG, int idx,float *u, float *X, float *P, float Fox, float Foy)
{
float X1;
/*
	if(idx==0) //inferior esquerdo
		X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx+1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx+1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx+TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[idx+TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(0.5)*((P[idx+TELX*TELY]-P[idx])/(2.0*Fox+2.0*Foy+1));
	else if(idx==TELX-1)//inferior direito
		X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx-1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx-1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx+TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[idx+TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(0.5)*((P[idx+TELX*TELY]-P[idx])/(2.0*Fox+2.0*Foy+1));
	else if(idx==TELX*(TELY-1))//superior esquerdo
		X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx+1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx+1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx-TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[idx-TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(0.5)*((P[idx+TELX*TELY]-P[idx])/(2.0*Fox+2.0*Foy+1));
	else if(idx==(TELX*TELY)-1)//superior direito
		X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx-1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx-1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx-TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[idx-TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(0.5)*((P[idx+TELX*TELY]-P[idx])/(2.0*Fox+2.0*Foy+1));
	else if(idx>0&&idx<(TELX-1))//contorno de baixo
		X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx-1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx+1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx+TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[idx+TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(0.5)*((P[idx+TELX*TELY]-P[idx])/(2.0*Fox+2.0*Foy+1));
	else if(idx>TELX*(TELY-1)&&idx<((TELX*TELY)-1))//contorno de cima
		X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx-1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx+1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx-TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[idx-TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(0.5)*((P[idx+TELX*TELY]-P[idx])/(2.0*Fox+2.0*Foy+1)); 
	else if(idx%TELX==0)//contorno esquerdo
		X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx+1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx+1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx-TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[idx+TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(0.5)*((P[idx+TELX*TELY]-P[idx])/(2.0*Fox+2.0*Foy+1));
	else if((idx+1)%TELX==0)//contorno direito
		X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx-1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx-1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx-TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[idx+TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(0.5)*((P[idx+TELX*TELY]-P[idx])/(2.0*Fox+2.0*Foy+1));
	else//meio
		X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx-1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx+1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx-TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[idx+TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(0.5)*((P[idx+TELX*TELY]-P[idx])/(2.0*Fox+2.0*Foy+1));
	
return X1;*/

if(idx>=0&&idx<TELX)//linha de baixo j==0
			{
				if(idx==0) //idx%TELX==0
				X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx+1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx+1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx+TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[idx+TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(0.5)*((P[idx+TELX*TELY]-P[idx])/(dt+1));//(2.0*Fox+2.0*Foy+1));
				
				else if(idx==TELX-1)//(idx+1)%TELX==0
				X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx-1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx-1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx+TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[idx+TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(0.5)*((P[idx+TELX*TELY]-P[idx])/(dt+1));//(2.0*Fox+2.0*Foy+1));
				
				else
				X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx-1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx+1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx+TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[idx+TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(0.5)*((P[idx+TELX*TELY]-P[idx])/(dt+1));//(2.0*Fox+2.0*Foy+1));
			}	

else if(idx>=TELX*(TELY-1)&&idx<TELX*TELY)//linha de cima j==TELY-1
			{
				if(idx==TELX*(TELY-1))
				X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx+1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx+1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx-TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[idx-TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(0.5)*((P[idx+TELX*TELY]-P[idx])/(2*dt+1));//(2.0*Fox+2.0*Foy+1));
				
				else if(idx==(TELX*TELY)-1)
				X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx-1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx-1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx-TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[idx-TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(0.5)*((P[idx+TELX*TELY]-P[idx])/(2*dt+1));//(2.0*Fox+2.0*Foy+1));
				
				else
				X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx-1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx+1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx-TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[idx-TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(0.5)*((P[idx+TELX*TELY]-P[idx])/(2*dt+1));//(2.0*Fox+2.0*Foy+1));
				
			}	
else//(1<=j<=tely-2)	// 1<=y<=tely-2  meio
			{
				if(idx%TELX==0)
				X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx+1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx+1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx-TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[idx+TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(0.5)*((P[idx+TELX*TELY]-P[idx])/(2*dt+1));//(2.0*Fox+2.0*Foy+1));
				
				else if((idx+1)%TELX==0)
				X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx-1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx-1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx-TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[idx+TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(0.5)*((P[idx+TELX*TELY]-P[idx])/(2*dt+1));//(2.0*Fox+2.0*Foy+1));					//0<x<L
				
				else{
				if(FLAG[idx+1]==false)
				X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx-1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx-1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx-TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[idx+TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(0.5)*((P[idx+TELX*TELY]-P[idx])/(2*dt+1));//(2.0*Fox+2.0*Foy+1));
				else if(FLAG[idx-1]==false)
				X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx+1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx+1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx-TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[idx+TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(0.5)*((P[idx+TELX*TELY]-P[idx])/(2*dt+1));//(2.0*Fox+2.0*Foy+1));
				else if(FLAG[idx+TELX]==false)
				X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx-1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx+1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx-TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[idx-TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(0.5)*((P[idx+TELX*TELY]-P[idx])/(2*dt+1));//(2.0*Fox+2.0*Foy+1));
				else if(FLAG[idx-TELX]==false)
				X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx-1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx+1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx+TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[idx+TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(0.5)*((P[idx+TELX*TELY]-P[idx])/(2*dt+1));//(2.0*Fox+2.0*Foy+1));				 //x=L
				else
				X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx-1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx+1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx-TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[idx+TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(0.5)*((P[idx+TELX*TELY]-P[idx])/(2*dt+1));//(2.0*Fox+2.0*Foy+1));				 //x=L
				}
			}

return X1;
}

//FUNCAO QUE CALCULA TAU NA GPU
/*__device__ float tau(bool *FLAG, float *P, int idx, float dx, float dy)
{
float tau;
tau=TAU0*powf(E(FLAG,P,idx,dx,dy),2.0);
return tau;
}*/

//FUNÇÃO QUE CALCULA E NA GPU
__device__ float E(bool *FLAG, float *P, int idx, float dx, float dy)
{
float E;
float Pxf=Px(FLAG,P,idx,dx);
float Pyf=Py(FLAG,P,idx,dy);
	if((Pxf<0.00000001&&Pxf>-0.00000001)&&(Pyf<0.00000001&&Pyf>-0.00000001))
	E=(1.0-3.0*S);
	else
	E=(1.0-3.0*S)*(1.0+((4.0*S)/(1.0-3.0*S))*((powf(Pxf,4)+powf(Pyf,4))/powf((powf(Pxf,2)+powf(Pyf,2)),2)));	
return E;
}


//FUNÇÕES QUE CALCULAM AS DERIVADAS DE P NA GPU
__device__ float Px(bool *FLAG, float *P, int idx, float dx)
{
float Px;

if(idx%TELX==0)
Px=(P[idx+1]-P[idx+1])/(2.0*dx);
else if((idx+1)%TELX==0)
Px=(P[idx-1]-P[idx-1])/(2.0*dx);
else{
	if(FLAG[idx+1]==false)
	Px=(P[idx-1]-P[idx-1])/(2.0*dx);
	else if(FLAG[idx-1]==false)//POTENCIALMENTE PERIGOSO, COMO PRETENDO USAR SOMENTE UM RETANGULO NO MEIO NAO DEVE DAR PROBLEMA
	Px=(P[idx+1]-P[idx+1])/(2.0*dx);
	else
	Px=(P[idx+1]-P[idx-1])/(2.0*dx);
	}
return Px;
}

__device__ float Py(bool *FLAG, float *P, int idx, float dy)
{
float Py;

if(idx>=0&&idx<TELX)
Py=(P[idx+TELX]-P[idx+TELX])/(2.0*dy);
else if(idx>=TELX*(TELY-1)&&idx<TELX*TELY)
Py=(P[idx-TELX]-P[idx-TELX])/(2.0*dy);
else{
	if(FLAG[idx+TELX]==false)
	Py=(P[idx-TELX]-P[idx-TELX])/(2.0*dy);
	else if(FLAG[idx-TELX]==false)
	Py=(P[idx+TELX]-P[idx+TELX])/(2.0*dy);
	else
	Py=(P[idx+TELX]-P[idx-TELX])/(2.0*dy);
	}
return Py;
}


__device__ float Pxx(bool *FLAG, float *P, int idx, float dx)
{
float Pxx;
if(idx%TELX==0)
Pxx=(P[idx+1]-2.0*P[idx]+P[idx+1])/(dx*dx);
else if((idx+1)%TELX==0)
Pxx=(P[idx-1]-2.0*P[idx]+P[idx-1])/(dx*dx);
else{
	if(FLAG[idx+1]==false)
	Pxx=(P[idx-1]-2.0*P[idx]+P[idx-1])/(dx*dx);
	else if(FLAG[idx-1]==false)
	Pxx=(P[idx+1]-2.0*P[idx]+P[idx+1])/(dx*dx);
	else
	Pxx=(P[idx+1]-2.0*P[idx]+P[idx-1])/(dx*dx);
	}
return Pxx;
}
__device__ float Pyy(bool *FLAG, float *P, int idx, float dy)
{
float Pyy;
if(idx>=0&&idx<TELX)
Pyy=(P[idx+TELX]-2.0*P[idx]+P[idx+TELX])/(dy*dy);
else if(idx>=TELX*(TELY-1)&&idx<TELX*TELY)
Pyy=(P[idx-TELX]-2.0*P[idx]+P[idx-TELX])/(dy*dy);
else{
	if(FLAG[idx+TELX]==false)
	Pyy=(P[idx-TELX]-2.0*P[idx]+P[idx-TELX])/(dy*dy);
	else if(FLAG[idx-TELX]==false)
	Pyy=(P[idx+TELX]-2.0*P[idx]+P[idx+TELX])/(dy*dy);
	else
	Pyy=(P[idx+TELX]-2.0*P[idx]+P[idx-TELX])/(dy*dy);
	}
return Pyy;
}

__device__ float Pxy(bool *FLAG, float *P, int idx, float dx, float dy)
{//MANTER LONGE DAS BORDAS
float Pxy;
float DireitoCima;
float EsquerdoCima;
float DireitoBaixo;
float EsquerdoBaixo;



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
		else{
	if(FLAG[idx+1]==false){
	DireitoCima=P[idx-1+TELX];
	DireitoBaixo=P[idx-1-TELX];}
	else{
	DireitoCima=P[idx+1+TELX];
	DireitoBaixo=P[idx+1-TELX];}

	if(FLAG[idx-1]==false){
	EsquerdoCima=P[idx+1+TELX];
	EsquerdoBaixo=P[idx+1-TELX];}
	else{
	EsquerdoCima=P[idx-1+TELX];
	EsquerdoBaixo=P[idx-1-TELX];}

	if(FLAG[idx+TELX]==false){
	DireitoCima=P[idx+1-TELX];
	EsquerdoCima=P[idx-1-TELX];}
	else{
	DireitoCima=P[idx+1+TELX];
	EsquerdoCima=P[idx-1+TELX];}

	if(FLAG[idx-TELX]==false){
	DireitoBaixo=P[idx+1+TELX];
	EsquerdoBaixo=P[idx-1+TELX];}
	else{
	DireitoBaixo=P[idx+1-TELX];
	EsquerdoBaixo=P[idx-1-TELX];}
			//Pxy=(P[idx+1+TELX]-P[idx-1+TELX]-P[idx+1-TELX]+P[idx-1-TELX])/(4.0*dx*dy);
			Pxy=(DireitoCima-EsquerdoCima-DireitoBaixo+EsquerdoBaixo)/(4.0*dx*dy);}
	}
	return Pxy;
}


__device__ float Pyx(bool *FLAG, float *P, int idx, float dx, float dy)
{
	float Pxy;
			float DireitoCima;
	float EsquerdoCima;
	float DireitoBaixo;
	float EsquerdoBaixo;
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
	Pxy=(P[idx+1+TELX]-P[idx-1+TELX]-P[idx+1-TELX]+P[idx-1-TELX])/(4.0*dx*dy);*/
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
		else{
			//Pxy=(P[idx+1+TELX]-P[idx-1+TELX]-P[idx+1-TELX]+P[idx-1-TELX])/(4.0*dx*dy);
	
	if(FLAG[idx+1]==false){
	DireitoCima=P[idx-1+TELX];
	DireitoBaixo=P[idx-1-TELX];}
	else{
	DireitoCima=P[idx+1+TELX];
	DireitoBaixo=P[idx+1-TELX];}
	if(FLAG[idx-1]==false){
	EsquerdoCima=P[idx+1+TELX];
	EsquerdoBaixo=P[idx+1-TELX];}
	else{
	EsquerdoCima=P[idx-1+TELX];
	EsquerdoBaixo=P[idx-1-TELX];}

	if(FLAG[idx+TELX]==false){
	DireitoCima=P[idx+1-TELX];
	EsquerdoCima=P[idx-1-TELX];}
	else{
	DireitoCima=P[idx+1+TELX];
	EsquerdoCima=P[idx-1+TELX];}

	if(FLAG[idx-TELX]==false){
	DireitoBaixo=P[idx+1+TELX];
	EsquerdoBaixo=P[idx-1+TELX];}
	else{
	DireitoBaixo=P[idx+1-TELX];
	EsquerdoBaixo=P[idx-1-TELX];}
			Pxy=(DireitoCima-EsquerdoCima-DireitoBaixo+EsquerdoBaixo)/(4.0*dx*dy);}
	}
	return Pxy;
}

__device__ float Epx(bool *FLAG, float *P,int idx, float dx, float dy)
{
float Epx;
float Ef=E(FLAG,P,idx,dx,dy);
float Pxf=Px(FLAG,P,idx,dx);
float Pyf=Py(FLAG,P,idx,dy);

if((Pxf<0.00000001&&Pxf>-0.00000001)&&(Pyf<0.00000001&&Pyf>-0.00000001))
	Epx=0;
	else
	Epx=16.0*S*((powf(Pxf,3)*powf(Pyf,2)-Pxf*powf(Pyf,4))/powf((powf(Pxf,2)+powf(Pyf,2)),3));
return Epx;
}



__device__ float Epy(bool *FLAG, float *P, int idx, float dx, float dy)//verificar se preciso declarar Pxf e Pyf
{
float Epy;
float Ef=E(FLAG,P,idx,dx,dy);
float Pxf=Px(FLAG,P,idx,dx);
float Pyf=Py(FLAG,P,idx,dy);
if((Pxf<0.00000001&&Pxf>-0.00000001)&&(Pyf<0.00000001&&Pyf>-0.00000001))
	Epy=0;
	else
	Epy=16.0*S*((powf(Pyf,3)*powf(Pxf,2)-Pyf*powf(Pxf,4))/powf((powf(Pxf,2)+powf(Pyf,2)),3));
return Epy;
}

__device__ float Ex(bool *FLAG, float *P, int idx, float dx, float dy)
{
float Ex;
float Epxf=Epx(FLAG,P,idx,dx,dy);
float Epyf=Epy(FLAG,P,idx,dx,dy);
float Pxxf=Pxx(FLAG,P,idx,dx);
float Pxyf=Pxy(FLAG,P,idx,dx,dy);

Ex=Epxf*Pxxf+Epyf*Pxyf;
return Ex;
}


__device__ float Ey(bool *FLAG, float *P, int idx, float dx, float dy)
{
float Ey;
float Epxf=Epx(FLAG,P,idx,dx,dy);
float Epyf=Epy(FLAG,P,idx,dx,dy);
float Pyyf=Pyy(FLAG,P,idx,dy);
float Pyxf=Pyx(FLAG,P,idx,dx,dy);

Ey=Epyf*Pyyf+Epxf*Pyxf;
return Ey;
}

__device__ float Epxpx(bool *FLAG, float *P, int idx, float dx, float dy)
{
float Epxpx;
float Epxf=Epx(FLAG,P,idx,dx,dy);
float Pxf=Px(FLAG,P,idx,dx);
float Pyf=Py(FLAG,P,idx,dy);
if((Pxf<0.00000001&&Pxf>-0.00000001)&&(Pyf<0.00000001&&Pyf>-0.00000001))
	Epxpx=0;
	else
	Epxpx=((16.0*S)/powf((powf(Pxf,2)+powf(Pyf,2)),3))*(3.0*powf(Pxf,2)*powf(Pyf,2)-powf(Pyf,4))-((6.0*Pxf)/(powf(Pxf,2)+powf(Pyf,2)))*Epxf;
return Epxpx;
}

__device__ float Epypy(bool *FLAG, float *P, int idx, float dx, float dy)
{
float Epypy;
float Epyf=Epy(FLAG,P,idx,dx,dy);
float Pxf=Px(FLAG,P,idx,dx);
float Pyf=Py(FLAG,P,idx,dy);
if((Pxf<0.00000001&&Pxf>-0.00000001)&&(Pyf<0.00000001&&Pyf>-0.00000001))
	Epypy=0;
	else
	Epypy=((16.0*S)/powf((powf(Pxf,2)+powf(Pyf,2)),3))*(3.0*powf(Pyf,2)*powf(Pxf,2)-powf(Pxf,4))-((6.0*Pyf)/(powf(Pxf,2)+powf(Pyf,2)))*Epyf;
return Epypy;
}


__device__ float Epxpy(bool *FLAG, float *P, int idx, float dx, float dy)
{
float Epxpy;
float Epyf=Epy(FLAG,P,idx,dx,dy);
float Pxf=Px(FLAG,P,idx,dx);
float Pyf=Py(FLAG,P,idx,dy);

if((Pxf<0.00000001&&Pxf>-0.00000001)&&(Pyf<0.00000001&&Pyf>-0.00000001))
	Epxpy=0;
	else
	Epxpy=((16.0*S)/powf((powf(Pxf,2)+powf(Pyf,2)),3))*(2.0*powf(Pyf,3)*Pxf-4.0*Pyf*powf(Pxf,3))-((6.0*Pxf)/(powf(Pxf,2)+powf(Pyf,2)))*Epyf;	

return Epxpy;
}

__device__ float Epypx(bool *FLAG, float *P,int idx, float dx, float dy)
{
float Epypx;
float Epxf=Epx(FLAG,P,idx,dx,dy);
float Pxf=Px(FLAG,P,idx,dx);
float Pyf=Py(FLAG,P,idx,dy);

if((Pxf<0.00000001&&Pxf>-0.00000001)&&(Pyf<0.00000001&&Pyf>-0.00000001))
	Epypx=0;
	else
	Epypx=((16.0*S)/powf((powf(Pxf,2)+powf(Pyf,2)),3))*(2.0*powf(Pxf,3)*Pyf-4.0*Pxf*powf(Pyf,3))-((6.0*Pyf)/(powf(Pxf,2)+powf(Pyf,2)))*Epxf;	

return Epypx;
}

__device__ float Expx(bool *FLAG, float *P,int idx,float dx, float dy)
{
float Expx;
float Epxpxf=Epxpx(FLAG,P,idx,dx,dy);
float Epypxf=Epypx(FLAG,P,idx,dx,dy);
float Pxxf=Pxx(FLAG,P,idx,dx);
float Pxyf=Pxy(FLAG,P,idx,dx,dy);

Expx=Epxpxf*Pxxf+Epypxf*Pxyf;
return Expx;
}


__device__ float Eypy(bool *FLAG, float *P, int idx, float dx, float dy)
{
float Eypy;
float Epypyf=Epypy(FLAG,P,idx,dx,dy);
float Epxpyf=Epxpy(FLAG,P,idx,dx,dy);
float Pyyf=Pyy(FLAG,P,idx,dy);
float Pyxf=Pyx(FLAG,P,idx,dx,dy);

Eypy=Epypyf*Pyyf+Epxpyf*Pyxf;
return Eypy;
}

