/************************************************************************************************************************/
/*															*/
/*					ARQUIVO DE ENTRADAS DO CÓDIGO							*/
/*															*/
/************************************************************************************************************************/

//CONSTANTES
#define M_PI 3.14159265358979323846
#define KK 0.355 //78/220
#define aa 0.461 //1,49E-5/3,23E-5
//#define constK 81.0//0.01
//#define h 3000
//#define Cp 2920000.0
#define E0 2.0//1.0//W0 normalmente 1.0 COMO W0/d0<6 FIZ W0/d0=5
//#define TAU0 1.0
#define D 1.0//0.0000323//3.0
#define M0 0.2257
#define a1 0.8839
#define a2 0.6267
//PARÂMETROS DE GEOMETRIA
#define compL 1800	//comprimento do domínio x/d0=902  COMPRIMENTO NOVO =10X DIAMETRO DA PARTICULA E x*=x/d0 (MANTIVE ANTIGO)
#define R 21.2//16.84		//raio inicial
#define TELX 1500	//total de elementos na direção x 1250 deltax*=1,92, mantendo 800 eu calculo o tamanho do domínio (MANTIVE ANTIGO)
#define TELY 1500	//total de elementos na direção y 1250
//PARÂMETROS DE TEMPO
#define Ttot 150000.0//40	//tempo total da simulação 3750 depois 938
#define dt 0.2//0.02		//tamanho do intervalo entre os passos de tempo
//PARÂMETROS DE SIMULAÇÃO
#define SOLIDO 1.0	//valor da variável de fase no estado sólido
#define LIQUIDO -1.0	//valor da variável de fase no estado líquido
#define Ueq 0		//super-resfriamento de equilíbrio
#define Ui 0.05//0.45		//super-resfriamento inicial 0.65 normalmente
//#define L 1048000000.0//0.5
#define S 0.0		//grau de anisotropia
#define err 0.000001	//erro máximo da iteração da solução implicita 0.00001

#define l 60//183 //largura da particula
#define h 10.6//8.7 // altura acima da superficie da partícula
#define PHI0 0.0
#define TAN 1.732//0.84
#define TAN2 1.732//0.84
//#define SEN 0.866