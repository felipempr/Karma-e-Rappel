/************************************************************************************************************************/
/*															*/
/*					ARQUIVO DE ENTRADAS DO CÓDIGO							*/
/*															*/
/************************************************************************************************************************/

//CONSTANTES
#define M_PI 3.14159265358979323846
//#define constK 0.01
//#define h 3000
//#define E0 1.0//W0 normalmente 1.0
//#define TAU0 1.0
//#define D 3.0
//#define a1 0.8839
//#define a2 0.6267

#define GAMMAAB 0.25	
#define GAMMAAC 0.25
#define GAMMABC 0.25
#define GAMMAABC 3.0

#define TAU 1.0
#define EPSILON 3.0
#define M 1.0
//PARÂMETROS DE GEOMETRIA
#define compL 250	//comprimento do domínio x/d0=902
#define R 15		//raio inicial
#define TELX 1000	//total de elementos na direção x 1250
#define TELY 1000	//total de elementos na direção y 1250
//PARÂMETROS DE TEMPO
#define Ttot 10000	//tempo total da simulação 3750 depois 938
#define dt 0.05		//tamanho do intervalo entre os passos de tempo 0.02
//PARÂMETROS DE SIMULAÇÃO
#define EXISTE 1	//valor da variável de fase no estado sólido
#define N_EXISTE 0	//valor da variável de fase no estado líquido
#define Ueq 0		//super-resfriamento de equilíbrio
#define Ui 0.45		//super-resfriamento inicial 0.65 normalmente
#define S 0.05		//grau de anisotropia
#define err 0.000001	//erro máximo da iteração da solução implicita 0.00001

