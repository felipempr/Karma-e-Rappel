void inicializar(int I1,int I2,int J1,int J2,int K1,int K2, float ***Q, float B);
float maior(int I,float *Q);
float P1(float s, float ***P, float ***u, int i, int j, int telx, int tely, float dx, float dy, float dt,float tau0, float E0, float lambda);
float X1(float ***u, float ***X, float ***P, int i, int j,int telx, int tely, float Fox, float Foy);
float Px(float ***P,int i,int j,float dx,int telx);
float Py(float ***P,int i,int j,float dy,int tely);
float Pxx(float ***P,int i,int j,float dx,int telx);
float Pyy(float ***P,int i,int j,float dy,int tely);
float Pxy(float ***P,int i,int j,float dx,float dy,int telx,int tely);
float Pyx(float ***P,int i,int j,float dx,float dy,int telx,int tely);
float Epx(float s,float ***P,int i,int j,float dx,float dy,int telx,int tely);
float Epy(float s,float ***P,int i,int j,float dx,float dy,int telx,int tely);
float Ex(float s,float ***P,int i,int j,float dx,float dy,int telx,int tely);
float Ey(float s,float ***P,int i,int j,float dx,float dy,int telx,int tely);
float Epxpx(float s,float ***P,int i,int j,float dx,float dy,int telx,int tely);
float Epypy(float s,float ***P,int i,int j,float dx,float dy,int telx,int tely);
float Epxpy(float s,float ***P,int i,int j,float dx,float dy,int telx,int tely);
float Epypx(float s,float ***P,int i,int j,float dx,float dy,int telx,int tely);
float Expx(float s,float ***P,int i,int j,float dx,float dy,int telx,int tely);
float Eypy(float s,float ***P,int i,int j,float dx,float dy,int telx,int tely);
float E(float s,float ***P,int i,int j,float dx,float dy,int telx,int tely);
float tau(float s,float ***P,int i,int j,float dx,float dy,int telx,int tely, float tau0);
float noise(float ***P,int i,int j);

//FUNÇÃO QUE inicializa AS MATRIZES
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
//FUNÇÃO PARA ESCOLHER O MAIOR VALOR EM UM VETOR
float maior(int I,float *Q)
{
float n=0;
for(int i=0;i<=I-1;i++)
{
	if(Q[i]>n)
	n=Q[i];
}
return n;
}
//FUNÇÃO QUE CALCULA TAU
float tau(float s,float ***P,int i,int j,float dx,float dy,int telx,int tely, float tau0)
{
float tau;
tau=tau0*pow(E(s,P,i,j,dx,dy,telx,tely),2.0);
return tau;
}

//FUNÇÃO QUE CALCULA NOVA VARIÁVEL DE FASE
float P1(float s, float ***P, float ***u, int i, int j,int telx, int tely, float dx, float dy,float dt,float tau0, float E0,float lambda)
{
float P1;
float I;
float II;
float III;

float Ef=E(s,P,i,j,dx,dy,telx,tely);
float Exf=Ex(s,P,i,j,dx,dy,telx,tely);
float Eyf=Ey(s,P,i,j,dx,dy,telx,tely);
float Epxf=Epx(s,P,i,j,dx,dy,telx,tely);
float Epyf=Epy(s,P,i,j,dx,dy,telx,tely);
float Expxf=Expx(s,P,i,j,dx,dy,telx,tely);
float Eypyf=Eypy(s,P,i,j,dx,dy,telx,tely);
float Pxf=Px(P,i,j,dx,telx);
float Pyf=Py(P,i,j,dy,tely);
float Pxxf=Pxx(P,i,j,dx,telx);
float Pyyf=Pyy(P,i,j,dy,tely);
float Pxyf=Pxy(P,i,j,dx,dy,telx,tely);
float tauf=tau(s,P,i,j,dx,dy,telx,tely,tau0);

I=pow(Ef,2)*(Pxxf+Pyyf)+2.0*Ef*(Pxf*Exf+Pyf*Eyf);
II= 2.0*Ef*Epxf*(Pxf*Pxxf+Pyf*Pxyf)+(pow(Pxf,2)+pow(Pyf,2))*(Epxf*Exf+Ef*Expxf);
III=2.0*Ef*Epyf*(Pyf*Pyyf+Pxf*Pxyf)+(pow(Pxf,2)+pow(Pyf,2))*(Epyf*Eyf+Ef*Eypyf);

	//P1=P[0][i][j]-M*dt*((pow(P[0][i][j],3.0)-1.5*pow(P[0][i][j],2.0)+0.5*P[0][i][j])+(noise(P,i,j)+(0.9/M_PI)*atan(10*(T[0][i][j]-TM)))*(-pow(P[0][i][j],2.0)+P[0][i][j])-E0*(I+II+III));
	P1=P[0][i][j]+(dt/tauf)*((P[0][i][j]-lambda*u[0][i][j]*(1.0-pow(P[0][i][j],2.0)))*(1.0-pow(P[0][i][j],2.0))+E0*(I+II+III));
		
return P1;
}
//FUNÇÃO QUE CALCULA E
float E(float s,float ***P,int i,int j,float dx,float dy,int telx,int tely)
{
float E;
float Pxf=Px(P,i,j,dx,telx);
float Pyf=Py(P,i,j,dy,tely);
	if((Pxf<0.00000001&&Pxf>-0.00000001)&&(Pyf<0.00000001&&Pyf>-0.00000001)){
	E=1.0-3.0*s;}
	else
	E=(1.0-3.0*s)*(1.0+((4.0*s)/(1.0-3.0*s))*((pow(Pxf,4)+pow(Pyf,4))/pow((pow(Pxf,2)+pow(Pyf,2)),2)));	
return E;
}
//FUNÇÃO QUE CALCULA NOVA TEMPERATURA
float X1(float ***u, float ***X, float ***P, int i, int j,int telx, int tely, float Fox, float Foy)
{
float X1;
	if(j==0)//linha de baixo
			{
				if(i==0)
				X1=(u[0][i][j]/(2.0*Fox+2.0*Foy+1))+(X[0][i+1][j]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[0][i+1][j]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[0][i][j+1]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[0][i][j+1]*(Foy/(2.0*Fox+2.0*Foy+1)))+(0.5)*((P[1][i][j]-P[0][i][j])/(2.0*Fox+2.0*Foy+1));
				//X1=(u[0][i][j]/(2.0*Fox+2.0*Foy+1))+(X[0][i+1][j]*(2.0*Fox/(2.0*Fox+2.0*Foy+1)))+(X[0][i][j+1]*(2.0*Foy/(2.0*Fox+2.0*Foy+1)))+(1/2)*((P[1][i][j]-P[0][i][j])/(2.0*Fox+2.0*Foy+1));//x=0,y=0
				//printf("==0 ");
				else if(i==telx-1)
				X1=(u[0][i][j]/(2.0*Fox+2.0*Foy+1))+(X[0][i-1][j]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[0][i-1][j]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[0][i][j+1]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[0][i][j+1]*(Foy/(2.0*Fox+2.0*Foy+1)))+(0.5)*((P[1][i][j]-P[0][i][j])/(2.0*Fox+2.0*Foy+1));
				//X1=(u[0][i][j]/(2.0*Fox+2.0*Foy+1))+(X[0][i-1][j]*(2.0*Fox/(2.0*Fox+2.0*Foy+1)))+(X[0][i][j+1]*(2.0*Foy/(2.0*Fox+2.0*Foy+1)))+(1/2)*((P[1][i][j]-P[0][i][j])/(2.0*Fox+2.0*Foy+1)); //y=0 parte de baixo
				else
				//printf("x=L y=%d\n",j);
				X1=(u[0][i][j]/(2.0*Fox+2.0*Foy+1))+(X[0][i-1][j]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[0][i+1][j]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[0][i][j+1]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[0][i][j+1]*(Foy/(2.0*Fox+2.0*Foy+1)))+(0.5)*((P[1][i][j]-P[0][i][j])/(2.0*Fox+2.0*Foy+1));
				//X1=(u[0][i][j]/(2.0*Fox+2.0*Foy+1))+(X[0][i-1][j]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[0][i+1][j]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[0][i][j+1]*(2.0*Foy/(2.0*Fox+2.0*Foy+1)))+(1/2)*((P[1][i][j]-P[0][i][j])/(2.0*Fox+2.0*Foy+1)); //y=0 parte de baixo
			}	
			else if(j==tely-1)//linha de cima
			{
				if(i==0)
				X1=(u[0][i][j]/(2.0*Fox+2.0*Foy+1))+(X[0][i+1][j]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[0][i+1][j]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[0][i][j-1]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[0][i][j-1]*(Foy/(2.0*Fox+2.0*Foy+1)))+(0.5)*((P[1][i][j]-P[0][i][j])/(2.0*Fox+2.0*Foy+1));
				//X1=u[0][i][j]/(2.0*Fox+2.0*Foy+1)+(X[0][i+1][j]*(2.0*Fox/(2.0*Fox+2.0*Foy+1)))+(X[0][i][j-1]*(2.0*Foy/(2.0*Fox+2.0*Foy+1)))+(1/2)*((P[1][i][j]-P[0][i][j])/(2.0*Fox+2.0*Foy+1));//x=0,y=tely-1
				//printf("==tely-1 ");
				else if(i==telx-1)
				X1=(u[0][i][j]/(2.0*Fox+2.0*Foy+1))+(X[0][i-1][j]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[0][i-1][j]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[0][i][j-1]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[0][i][j-1]*(Foy/(2.0*Fox+2.0*Foy+1)))+(0.5)*((P[1][i][j]-P[0][i][j])/(2.0*Fox+2.0*Foy+1));
				//X1=(u[0][i][j]/(2.0*Fox+2.0*Foy+1))+(X[0][i-1][j]*(2.0*Fox/(2.0*Fox+2.0*Foy+1)))+(X[0][i][j-1]*(2.0*Foy/(2.0*Fox+2.0*Foy+1)))+(1/2)*((P[1][i][j]-P[0][i][j])/(2.0*Fox+2.0*Foy+1));
				//y=tely-1 parte de cima
				else
				X1=(u[0][i][j]/(2.0*Fox+2.0*Foy+1))+(X[0][i-1][j]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[0][i+1][j]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[0][i][j-1]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[0][i][j-1]*(Foy/(2.0*Fox+2.0*Foy+1)))+(0.5)*((P[1][i][j]-P[0][i][j])/(2.0*Fox+2.0*Foy+1));
				//X1=u[0][i][j]/(2.0*Fox+2.0*Foy+1)+(X[0][i-1][j]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[0][i+1][j]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[0][i][j-1]*(2.0*Foy/(2.0*Fox+2.0*Foy+1)))+(1/2)*((P[1][i][j]-P[0][i][j])/(2.0*Fox+2.0*Foy+1));
				//x=telx-1, y=tely-1
				//printf("x=L y=%d\n",j);
			}	
			else//(1<=j<=tely-2)	// 1<=y<=tely-2  meio
			{
				if(i==0)
				X1=(u[0][i][j]/(2.0*Fox+2.0*Foy+1))+(X[0][i+1][j]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[0][i+1][j]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[0][i][j-1]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[0][i][j+1]*(Foy/(2.0*Fox+2.0*Foy+1)))+(0.5)*((P[1][i][j]-P[0][i][j])/(2.0*Fox+2.0*Foy+1));
				//X1=(u[0][i][j]/(2.0*Fox+2.0*Foy+1))+(X[0][i+1][j]*(2.0*Fox/(2.0*Fox+2.0*Foy+1)))+(X[0][i][j-1]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[0][i][j+1]*(Foy/(2.0*Fox+2.0*Foy+1)))+(1/2)*((P[1][i][j]-P[0][i][j])/(2.0*Fox+2.0*Foy+1));	//x=0
				//printf("outros 0 ");
				else if(i==telx-1)
				X1=(u[0][i][j]/(2.0*Fox+2.0*Foy+1))+(X[0][i-1][j]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[0][i-1][j]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[0][i][j-1]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[0][i][j+1]*(Foy/(2.0*Fox+2.0*Foy+1)))+(0.5)*((P[1][i][j]-P[0][i][j])/(2.0*Fox+2.0*Foy+1));					//0<x<L
				//X1=(u[0][i][j]/(2.0*Fox+2.0*Foy+1))+(X[0][i-1][j]*(2.0*Fox/(2.0*Fox+2.0*Foy+1)))+(X[0][i][j-1]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[0][i][j+1]*(Foy/(2.0*Fox+2.0*Foy+1)))+(1/2)*((P[1][i][j]-P[0][i][j])/(2.0*Fox+2.0*Foy+1));
				else
				X1=(u[0][i][j]/(2.0*Fox+2.0*Foy+1))+(X[0][i-1][j]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[0][i+1][j]*(Fox/(2.0*Fox+2.0*Foy+1)))+(X[0][i][j-1]*(Foy/(2.0*Fox+2.0*Foy+1)))+(X[0][i][j+1]*(Foy/(2.0*Fox+2.0*Foy+1)))+(0.5)*((P[1][i][j]-P[0][i][j])/(2.0*Fox+2.0*Foy+1));
				 //x=L
				//printf("x=19 y=%d\n",j);
			}
		
return X1;
}

//FUNÇÕES QUE CALCULAM AS DERIVADAS DE P
float Px(float***P, int i, int j,float dx,int telx)
{
float Px;
if(i==0)
//Px=0;
Px=(P[0][i+1][j]-P[0][i+1][j])/(2.0*dx);
else if (i==telx-1)
//Px=0;
Px=(P[0][i-1][j]-P[0][i-1][j])/(2.0*dx);
else
Px=(P[0][i+1][j]-P[0][i-1][j])/(2.0*dx);
return Px;
}

float Py(float***P, int i, int j,float dy,int tely)
{
float Py;
if(j==0)
//Py=0;
Py=(P[0][i][j+1]-P[0][i][j+1])/(2.0*dy);
else if(j==tely-1)
//Py=0;
Py=(P[0][i][j-1]-P[0][i][j-1])/(2.0*dy);
else
Py=(P[0][i][j+1]-P[0][i][j-1])/(2.0*dy);
return Py;
}

float Pxx(float ***P, int i, int j, float dx,int telx)
{
float Pxx;
if(i==0)
//Pxx=0;
Pxx=(P[0][i+1][j]-2.0*P[0][i][j]+P[0][i+1][j])/(dx*dx);
else if(i==telx-1)
//Pxx=0;
Pxx=(P[0][i-1][j]-2.0*P[0][i][j]+P[0][i-1][j])/(dx*dx);
else
Pxx=(P[0][i+1][j]-2.0*P[0][i][j]+P[0][i-1][j])/(dx*dx);
return Pxx;
}

float Pyy(float ***P, int i, int j, float dy, int tely)
{
float Pyy;
if(j==0)
//Pyy=0;
Pyy=(P[0][i][j+1]-2.0*P[0][i][j]+P[0][i][j+1])/(dy*dy);
else if(j==tely-1)
//Pyy=0;
Pyy=(P[0][i][j-1]-2.0*P[0][i][j]+P[0][i][j-1])/(dy*dy);
else
Pyy=(P[0][i][j+1]-2.0*P[0][i][j]+P[0][i][j-1])/(dy*dy);
return Pyy;
}

float Pxy(float ***P, int i, int j, float dx, float dy,int telx,int tely)
{
float Pxy;
/*if(i==0)
Pxy=(Py(P,i+1,j,dy,tely)-Py(P,i+1,j,dy,tely))/(2.0*dx);
else if(i==telx-1)
Pxy=(Py(P,i,j,dy,tely)-Py(P,i,j,dy,tely))/(2.0*dx);
else
Pxy=(Py(P,i+1,j,dy,tely)-Py(P,i-1,j,dy,tely))/(2.0*dx);*/

if(i==0)
	{
	if(j==0)
	//Pxy=0;
	Pxy=(P[0][i+1][j+1]-P[0][i+1][j+1]-P[0][i+1][j+1]+P[0][i+1][j+1])/(4.0*dx*dy);
	else if(j==tely-1)
	//Pxy=0;
	Pxy=(P[0][i+1][j-1]-P[0][i+1][j-1]-P[0][i+1][j-1]+P[0][i+1][j-1])/(4.0*dx*dy);
	else
	//Pxy=0;
	Pxy=(P[0][i+1][j+1]-P[0][i+1][j+1]-P[0][i+1][j-1]+P[0][i+1][j-1])/(4.0*dx*dy);
	}
else if(i==telx-1)
	{
	if(j==0)
	//Pxy=0;
	Pxy=(P[0][i-1][j+1]-P[0][i-1][j+1]-P[0][i-1][j+1]+P[0][i-1][j+1])/(4.0*dx*dy);
	else if(j==tely-1)
	//Pxy=0;
	Pxy=(P[0][i-1][j-1]-P[0][i-1][j-1]-P[0][i-1][j-1]+P[0][i-1][j-1])/(4.0*dx*dy);
	else
	//Pxy=0;
	Pxy=(P[0][i-1][j+1]-P[0][i-1][j+1]-P[0][i-1][j-1]+P[0][i-1][j-1])/(4.0*dx*dy);
	}
else
	{
	if(j==0)
	//Pxy=0;
	Pxy=(P[0][i+1][j+1]-P[0][i-1][j+1]-P[0][i+1][j+1]+P[0][i-1][j+1])/(4.0*dx*dy);
	else if(j==tely-1)
	//Pxy=0;
	Pxy=(P[0][i+1][j-1]-P[0][i-1][j-1]-P[0][i+1][j-1]+P[0][i-1][j-1])/(4.0*dx*dy);
	else
	Pxy=(P[0][i+1][j+1]-P[0][i-1][j+1]-P[0][i+1][j-1]+P[0][i-1][j-1])/(4.0*dx*dy);
	}
return Pxy;
}

float Pyx(float ***P, int i, int j, float dx, float dy,int telx,int tely)
{
float Pxy;
/*if(j==0)
Pyx=(Px(P,i,j+1,dx,telx)-Px(P,i,j+1,dx,telx))/(2.0*dy);
else if(j==tely-1)
Pyx=(Px(P,i,j-1,dx,telx)-Px(P,i,j-1,dx,telx))/(2.0*dy);
else
Pyx=(Px(P,i,j+1,dx,telx)-Px(P,i,j-1,dx,telx))/(2.0*dy);*/

if(i==0)
	{
	if(j==0)
	//Pxy=0;
	Pxy=(P[0][i+1][j+1]-P[0][i+1][j+1]-P[0][i+1][j+1]+P[0][i+1][j+1])/(4.0*dx*dy);
	else if(j==tely-1)
	//Pxy=0;
	Pxy=(P[0][i+1][j-1]-P[0][i+1][j-1]-P[0][i+1][j-1]+P[0][i+1][j-1])/(4.0*dx*dy);
	else
	//Pxy=0;
	Pxy=(P[0][i+1][j+1]-P[0][i+1][j+1]-P[0][i+1][j-1]+P[0][i+1][j-1])/(4.0*dx*dy);
	}
else if(i==telx-1)
	{
	if(j==0)
	//Pxy=0;
	Pxy=(P[0][i-1][j+1]-P[0][i-1][j+1]-P[0][i-1][j+1]+P[0][i-1][j+1])/(4.0*dx*dy);
	else if(j==tely-1)
	//Pxy=0;
	Pxy=(P[0][i-1][j-1]-P[0][i-1][j-1]-P[0][i-1][j-1]+P[0][i-1][j-1])/(4.0*dx*dy);
	else
	//Pxy=0;
	Pxy=(P[0][i-1][j+1]-P[0][i-1][j+1]-P[0][i-1][j-1]+P[0][i-1][j-1])/(4.0*dx*dy);
	}
else
	{
	if(j==0)
	//Pxy=0;
	Pxy=(P[0][i+1][j+1]-P[0][i-1][j+1]-P[0][i+1][j+1]+P[0][i-1][j+1])/(4.0*dx*dy);
	else if(j==tely-1)
	//Pxy=0;
	Pxy=(P[0][i+1][j-1]-P[0][i-1][j-1]-P[0][i+1][j-1]+P[0][i-1][j-1])/(4.0*dx*dy);
	else
	Pxy=(P[0][i+1][j+1]-P[0][i-1][j+1]-P[0][i+1][j-1]+P[0][i-1][j-1])/(4.0*dx*dy);
	}
return Pxy;
}
//FUNÇÃO QUE ADICIONA RUÍDO
float noise(float ***P, int i, int j)
{
float a;
if(P[0][i][j]>=-0.005&&P[0][i][j]<=0.005)
a=0.01;
else
a=0;
float r=((rand()%1000)/1000.0)-0.5;
return a*r;
}

//FUNÇÕES QUE CALCULAM AS DERIVADAS DE E
float Epx(float s,float ***P,int i,int j,float dx,float dy,int telx,int tely)
{
float Epx;
float Pxf=Px(P,i,j,dx,telx);
float Pyf=Py(P,i,j,dy,tely);

if((Pxf<0.00000001&&Pxf>-0.00000001)&&(Pyf<0.00000001&&Pyf>-0.00000001)){
	Epx=0;}
	else{
Epx=16.0*s*((pow(Pxf,3)*pow(Pyf,2)-Pxf*pow(Pyf,4))/pow((pow(Pxf,2)+pow(Pyf,2)),3));}
return Epx;
}

float Epy(float s,float ***P,int i,int j,float dx,float dy,int telx,int tely)
{
float Epy;
float Pxf=Px(P,i,j,dx,telx);
float Pyf=Py(P,i,j,dy,tely);

if((Pxf<0.00000001&&Pxf>-0.00000001)&&(Pyf<0.00000001&&Pyf>-0.00000001)){
	Epy=0;}
	else{
Epy=16.0*s*((pow(Pyf,3)*pow(Pxf,2)-Pyf*pow(Pxf,4))/pow((pow(Pxf,2)+pow(Pyf,2)),3));}
return Epy;
}

float Ex(float s,float ***P,int i,int j,float dx,float dy,int telx,int tely)
{
float Ex;
float Epxf=Epx(s,P,i,j,dx,dy,telx,tely);
float Epyf=Epy(s,P,i,j,dx,dy,telx,tely);
float Pxxf=Pxx(P,i,j,dx,telx);
float Pxyf=Pxy(P,i,j,dx,dy,telx,tely);

//Ex=Epx(s,P,i,j,dx,dy)*Pxx(P,i,j,dx)+Epy(s,P,i,j,dx,dy)*Pxy(P,i,j,dx,dy);
Ex=Epxf*Pxxf+Epyf*Pxyf;
return Ex;
}

float Ey(float s,float ***P,int i,int j,float dx,float dy,int telx,int tely)
{
float Ey;

float Epxf=Epx(s,P,i,j,dx,dy,telx,tely);
float Epyf=Epy(s,P,i,j,dx,dy,telx,tely);
float Pyyf=Pyy(P,i,j,dy,tely);
float Pyxf=Pyx(P,i,j,dx,dy,telx,tely);

//Ey=Epy(s,P,i,j,dx,dy)*Pyy(P,i,j,dy)+Epx(s,P,i,j,dx,dy)*Pxy(P,i,j,dx,dy);
Ey=Epyf*Pyyf+Epxf*Pyxf;
return Ey;
}

float Epxpx(float s,float ***P,int i,int j,float dx,float dy,int telx,int tely)
{
float Epxpx;
float Epxf=Epx(s,P,i,j,dx,dy,telx,tely);
float Pxf=Px(P,i,j,dx,telx);
float Pyf=Py(P,i,j,dy,tely);

if((Pxf<0.00000001&&Pxf>-0.00000001)&&(Pyf<0.00000001&&Pyf>-0.00000001)){
	Epxpx=0;}
	else{
	//Epxpx=16.0*s*(3.0*pow(Px(P,i,j,dx),2)*pow(Py(P,i,j,dy),2)-pow(Py(P,i,j,dy),4))/(pow((pow(Px(P,i,j,dx),2)+pow(Py(P,i,j,dy),2)),3))-(6.0*Px(P,i,j,dx)/(pow(Px(P,i,j,dx),2)+pow(Py(P,i,j,dy),2)))*Epx(s,P,i,j,dx,dy);
	Epxpx=((16.0*s)/pow((pow(Pxf,2)+pow(Pyf,2)),3))*(3.0*pow(Pxf,2)*pow(Pyf,2)-pow(Pyf,4))-((6.0*Pxf)/(pow(Pxf,2)+pow(Pyf,2)))*Epxf;
	}		
return Epxpx;
}

float Epypy(float s,float ***P,int i,int j,float dx,float dy,int telx,int tely)
{
float Epypy;
float Epyf=Epy(s,P,i,j,dx,dy,telx,tely);
float Pxf=Px(P,i,j,dx,telx);
float Pyf=Py(P,i,j,dy,tely);

if((Pxf<0.00000001&&Pxf>-0.00000001)&&(Pyf<0.00000001&&Pyf>-0.00000001)){
	Epypy=0;}
	else{
	Epypy=((16.0*s)/pow((pow(Pxf,2)+pow(Pyf,2)),3))*(3.0*pow(Pyf,2)*pow(Pxf,2)-pow(Pxf,4))-((6.0*Pyf)/(pow(Pxf,2)+pow(Pyf,2)))*Epyf;}
return Epypy;
}

float Epxpy(float s,float ***P,int i,int j,float dx,float dy,int telx,int tely)
{
float Epxpy;
float Epyf=Epy(s,P,i,j,dx,dy,telx,tely);
float Pxf=Px(P,i,j,dx,telx);
float Pyf=Py(P,i,j,dy,tely);
if((Pxf<0.00000001&&Pxf>-0.00000001)&&(Pyf<0.00000001&&Pyf>-0.00000001)){
	Epxpy=0;}
	else{
	Epxpy=((16.0*s)/pow((pow(Pxf,2)+pow(Pyf,2)),3))*(2.0*pow(Pyf,3)*Pxf-4.0*Pyf*pow(Pxf,3))-((6.0*Pxf)/(pow(Pxf,2)+pow(Pyf,2)))*Epyf;}
return Epxpy;
}

float Epypx(float s,float ***P,int i,int j,float dx,float dy,int telx,int tely)
{
float Epypx;
float Epxf=Epx(s,P,i,j,dx,dy,telx,tely);
float Pxf=Px(P,i,j,dx,telx);
float Pyf=Py(P,i,j,dy,tely);
if((Pxf<0.00000001&&Pxf>-0.00000001)&&(Pyf<0.00000001&&Pyf>-0.00000001)){
	Epypx=0;}
	else{
	Epypx=((16.0*s)/pow((pow(Pxf,2)+pow(Pyf,2)),3))*(2.0*pow(Pxf,3)*Pyf-4.0*Pxf*pow(Pyf,3))-((6.0*Pyf)/(pow(Pxf,2)+pow(Pyf,2)))*Epxf;}
//Epypx=((16.0*s)/ pow((pow(Px(P,i,j,dx),2)+pow(Py(P,i,j,dy),2)),3)) *(2.0*pow(Px(P,i,j,dx),3)*Py(P,i,j,dy)-4.0*Px(P,i,j,dx)*pow(Py(P,i,j,dy),3)) - (6.0*Py(P,i,j,dy)/(pow(Px(P,i,j,dx),2)+pow(Py(P,i,j,dy),2)))*Epx(s,P,i,j,dx,dy);
return Epypx;
}

float Expx(float s,float ***P,int i,int j,float dx,float dy,int telx,int tely)
{
float Expx;
float Epxpxf=Epxpx(s,P,i,j,dx,dy,telx,tely);
float Epypxf=Epypx(s,P,i,j,dx,dy,telx,tely);
float Pxyf=Pxy(P,i,j,dx,dy,telx,tely);
float Pxxf=Pxx(P,i,j,dx,telx);
//Expx=Epxpx(s,P,i,j,dx,dy,telx,tely)*Pxx(P,i,j,dx)+Epypx(s,P,i,j,dx,dy,telx,tely)*Pxy(P,i,j,dx,dy);
Expx=Epxpxf*Pxxf+Epypxf*Pxyf;
return Expx;
}

float Eypy(float s,float ***P,int i,int j,float dx,float dy,int telx,int tely)
{
float Eypy;
float Epypyf=Epypy(s,P,i,j,dx,dy,telx,tely);
float Epxpyf=Epxpy(s,P,i,j,dx,dy,telx,tely);
float Pyxf=Pyx(P,i,j,dx,dy,telx,tely);
float Pyyf=Pyy(P,i,j,dy,tely);
//Eypy=Epypy(s,P,i,j,dx,dy,telx,tely)*Pyy(P,i,j,dy)+Epxpy(s,P,i,j,dx,dy,telx,tely)*Pxy(P,i,j,dx,dy);
Eypy=Epypyf*Pyyf+Epxpyf*Pyxf;
return Eypy;
}

