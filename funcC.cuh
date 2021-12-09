//__global__ void testef(float *teste);

__device__ float maior(float *delta);
__device__ void deltaf(float *delta, float *X);
__device__ void atualiza(float *Q, float *K);
__global__ void atualizaTudo(float *Q, float *K);
__global__ void P1(float *P, float *u, float dx, float dy, float lambda);
__global__ void Temp(float *u, float *X, float *P,float *delta, float Fox, float Foy);
__device__ float X1(int idx,float *u, float *X, float *P, float Fox, float Foy);
__device__ float tau(float *P, int idx, float dx, float dy);
__device__ float E(float *P, int idx, float dx, float dy);
__device__ float Px(float *P, int idx, float dx);
__device__ float Py(float *P, int idx, float dy);
__device__ float Pxx(float *P, int idx, float dx);
__device__ float Pyy(float *P, int idx, float dy);
__device__ float Pxy(float *P, int idx, float dx, float dy);
__device__ float Pyx(float *P, int idx, float dx, float dy);
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

/*__global__ void testef(float *teste)
{
int index = blockIdx.x*blockDim.x + threadIdx.x;
	int stride = blockDim.x*gridDim.x;
	for (int idx = index; idx < 4*4; idx += stride) {
    teste[idx]=idx;

if(idx>=0&&idx<4)
		{
				if(idx==0) //idx%TELX==0
				teste[idx+4*4]=1;
				else if(idx==(4-1))//(idx+1)%TELX==0
				teste[idx+4*4]=2;
				else
				teste[idx+4*4]=3;
				
			}	

else if(idx>=4*(4-1)&&idx<4*4)//linha de cima j==TELY-1
			{
				if(idx==4*(4-1))
				teste[idx+4*4]=4;
				else if(idx==(4*4)-1)
				teste[idx+4*4]=5;
				else
				teste[idx+4*4]=6;
			}	
else//(1<=j<=tely-2)	// 1<=y<=tely-2  meio
			{
				if(idx%4==0)
				teste[idx+4*4]=7;
				else if((idx+1)%4==0)
				teste[idx+4*4]=8;
				else
				teste[idx+4*4]=9;
			}

}
}
*/
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

//FUNÇÃO QUE CALCULA NOVA VARIÁVEL DE FASE na GPU
__global__ void P1(float *P, float *u, float dx, float dy, float lambda)
{
	float I;
	float II;
	float III;

	float Ef;
	float Exf;
	float Eyf;
	float Epxf;
	float Epyf;
	float Expxf;
	float Eypyf;
	float Pxf;
	float Pyf;
	float Pxxf;
	float Pyyf;
	float Pxyf;
	float tauf;

	//int idx=blockIdx.x * blockDim.x + threadIdx.x;
	//if(idx<TELX*TELY){
	int index = blockIdx.x*blockDim.x+threadIdx.x;
	int stride = blockDim.x*gridDim.x;
	for (int idx = index; idx < TELX*TELY; idx += stride) {

		Ef = E(P, idx, dx, dy);
		Exf = Ex(P, idx, dx, dy);
		Eyf = Ey(P, idx, dx, dy);
		Epxf = Epx(P, idx, dx, dy);
		Epyf = Epy(P, idx, dx, dy);
		Expxf = Expx(P, idx, dx, dy);
		Eypyf = Eypy(P, idx, dx, dy);
		Pxf = Px(P, idx, dx);
		Pyf = Py(P, idx, dy);
		Pxxf = Pxx(P, idx, dx);
		Pyyf = Pyy(P, idx, dy);
		Pxyf = Pxy(P, idx, dx, dy);
		tauf = tau(P, idx, dx, dy);

		I=powf(Ef, 2)*(Pxxf + Pyyf)+2.0*Ef*(Pxf*Exf+Pyf*Eyf);
		II=2.0*Ef*Epxf*(Pxf*Pxxf+Pyf*Pxyf)+(powf(Pxf,2)+powf(Pyf,2))*(Epxf*Exf+Ef*Expxf);
		III=2.0*Ef*Epyf*(Pyf*Pyyf+Pxf*Pxyf)+(powf(Pxf,2)+powf(Pyf,2))*(Epyf*Eyf+Ef*Eypyf);

		//P1=P[0][i][j]-M*dt*((pow(P[0][i][j],3.0)-1.5*pow(P[0][i][j],2.0)+0.5*P[0][i][j])+(noise(P,i,j)+(0.9/M_PI)*atan(10*(T[0][i][j]-TM)))*(-pow(P[0][i][j],2.0)+P[0][i][j])-E0*(I+II+III));

		P[idx+TELX*TELY]=P[idx]+(dt/tauf)*((P[idx]-lambda*u[idx]*(1.0-powf(P[idx],2.0)))*(1.0-powf(P[idx],2.0))+(E0*E0)*(I+II+III));
	//P[idx + TELX*TELY]=P[idx]+(dt/tauf)*((P[idx]-lambda*u[idx]*(1.0-powf(P[idx],2.0)))*(1.0-powf(P[idx],2.0))+E0*(I+II+III));

		//printf("P[5][5]=%f P[30][30]=%f\n", P[5 + 5 * TELX], P[30 + 30 * TELX]);
		//printf("P1[5][5]=%f P1[30][30]=%f\n", P[5 + 5 * TELX + TELX*TELY], P[30 + 30 * TELX + TELX*TELY]);

	}
}
__global__ void Temp(float *u, float *X, float *P,float *delta, float Fox, float Foy)
{
float m=0;
//printf(" x");
	do{
	int index = blockIdx.x*blockDim.x+threadIdx.x;
	int stride = blockDim.x*gridDim.x;
	for (int idx = index; idx < TELX*TELY; idx += stride) {	
		X[idx+TELX*TELY]=X1(idx,u, X, P, Fox, Foy);
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

__device__ float X1(int idx,float *u, float *X, float *P, float Fox, float Foy)
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
				X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx+1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx+1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx+TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[idx+TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(0.5)*((P[idx+TELX*TELY]-P[idx])/(2.0*Fox+2.0*Foy+1));
				
				else if(idx==TELX-1)//(idx+1)%TELX==0
				X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx-1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx-1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx+TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[idx+TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(0.5)*((P[idx+TELX*TELY]-P[idx])/(2.0*Fox+2.0*Foy+1));
				
				else
				X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx-1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx+1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx+TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[idx+TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(0.5)*((P[idx+TELX*TELY]-P[idx])/(2.0*Fox+2.0*Foy+1));
			}	

else if(idx>=TELX*(TELY-1)&&idx<TELX*TELY)//linha de cima j==TELY-1
			{
				if(idx==TELX*(TELY-1))
				X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx+1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx+1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx-TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[idx-TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(0.5)*((P[idx+TELX*TELY]-P[idx])/(2.0*Fox+2.0*Foy+1));
				
				else if(idx==(TELX*TELY)-1)
				X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx-1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx-1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx-TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[idx-TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(0.5)*((P[idx+TELX*TELY]-P[idx])/(2.0*Fox+2.0*Foy+1));
				
				else
				X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx-1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx+1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx-TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[idx-TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(0.5)*((P[idx+TELX*TELY]-P[idx])/(2.0*Fox+2.0*Foy+1));
				
			}	
else//(1<=j<=tely-2)	// 1<=y<=tely-2  meio
			{
				if(idx%TELX==0)
				X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx+1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx+1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx-TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[idx+TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(0.5)*((P[idx+TELX*TELY]-P[idx])/(2.0*Fox+2.0*Foy+1));
				
				else if((idx+1)%TELX==0)
				X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx-1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx-1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx-TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[idx+TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(0.5)*((P[idx+TELX*TELY]-P[idx])/(2.0*Fox+2.0*Foy+1));					//0<x<L
				
				else
				X1=(u[idx]/(2.0*Fox+2.0*Foy+1))+(X[idx-1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx+1]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[idx-TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[idx+TELX]*(Foy/(2.0*Fox+2.0*Foy+1)))+(0.5)*((P[idx+TELX*TELY]-P[idx])/(2.0*Fox+2.0*Foy+1));				 //x=L
			}

return X1;
}

//FUNCAO QUE CALCULA TAU NA GPU
__device__ float tau(float *P, int idx, float dx, float dy)
{
float tau;
tau=TAU0*powf(E(P,idx,dx,dy),2.0);
return tau;
}

//FUNÇÃO QUE CALCULA E NA GPU
__device__ float E(float *P, int idx, float dx, float dy)
{
float E;
float Pxf=Px(P,idx,dx);
float Pyf=Py(P,idx,dy);
	if((Pxf<0.00000001&&Pxf>-0.00000001)&&(Pyf<0.00000001&&Pyf>-0.00000001))
	E=(1.0-3.0*S);
	else
	E=(1.0-3.0*S)*(1.0+((4.0*S)/(1.0-3.0*S))*((powf(Pxf,4)+powf(Pyf,4))/powf((powf(Pxf,2)+powf(Pyf,2)),2)));	
return E;
}


//FUNÇÕES QUE CALCULAM AS DERIVADAS DE P NA GPU
__device__ float Px(float *P, int idx, float dx)
{
float Px;
if(idx%TELX==0)
Px=(P[idx+1]-P[idx+1])/(2.0*dx);
else if((idx+1)%TELX==0)
Px=(P[idx-1]-P[idx-1])/(2.0*dx);
else
Px=(P[idx+1]-P[idx-1])/(2.0*dx);
return Px;
}

__device__ float Py(float *P, int idx, float dy)
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


__device__ float Pxx(float *P, int idx, float dx)
{
float Pxx;
if(idx%TELX==0)
Pxx=(P[idx+1]-2.0*P[idx]+P[idx+1])/(dx*dx);
else if((idx+1)%TELX==0)
Pxx=(P[idx-1]-2.0*P[idx]+P[idx-1])/(dx*dx);
else
Pxx=(P[idx+1]-2.0*P[idx]+P[idx-1])/(dx*dx);
return Pxx;
}
__device__ float Pyy(float *P, int idx, float dy)
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

__device__ float Pxy(float *P, int idx, float dx, float dy)
{
float Pxy;
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


__device__ float Pyx(float *P, int idx, float dx, float dy)
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
		else
			Pxy=(P[idx+1+TELX]-P[idx-1+TELX]-P[idx+1-TELX]+P[idx-1-TELX])/(4.0*dx*dy);
	}
	return Pxy;
}

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

