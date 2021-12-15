/************************************************************************************************************************/
/*															*/
/*					ARQUIVO DE ENTRADAS DO CÓDIGO							*/
/*															*/
/************************************************************************************************************************/

//CONSTANTES
#define constK 0.01
#define h 3000
#define E0 1.0//W0 normalmente 1.0
#define TAU0 1.0
#define D 4.4312
#define a1 0.8839
#define a2 0.6267

#define Li -0.5 //considerando LGB=LTJ=1
#define LGB 1.0
#define LTJ 0.01
#define phiW 5000.0
#define phiMin 1.04

#define alfa 1.0
#define beta 1.0
#define gamma 1.5
#define kappa 1.0
//PARÂMETROS DE GEOMETRIA
#define compL 50	//comprimento do domínio era 500
#define R 15		//raio inicial
#define TELX 300	//total de elementos na direção x 1250
#define TELY 300	//total de elementos na direção y 1250
//PARÂMETROS DE TEMPO
#define Ttot 250	//tempo total da simulação 3750
#define dt 0.005		//tamanho do intervalo entre os passos de tempo 0.02
//PARÂMETROS DE SIMULAÇÃO
#define SOLIDO 1	//valor da variável de fase no estado sólido
#define LIQUIDO -1	//valor da variável de fase no estado líquido
#define Ueq 0		//super-resfriamento de equilíbrio
#define Ui 0.55		//super-resfriamento inicial 0.65 normalmente
#define S 0.05		//grau de anisotropia
#define err 0.0000001	//erro máximo da iteração da solução implicita 0.00001

