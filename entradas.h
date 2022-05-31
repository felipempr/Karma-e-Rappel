/************************************************************************************************************************/
/*															*/
/*					ARQUIVO DE ENTRADAS DO CÓDIGO							*/
/*															*/
/************************************************************************************************************************/

//CONSTANTES
#define M_PI 3.14159265358979323846
#define constK 1.0
#define cv 10.0
//#define h 3000
//#define E0 1.0//W0 normalmente 1.0
#define D 3.0
//#define a1 0.8839
//#define a2 0.6267

#define GAMMAAB 0.25 //era 0.25	
#define GAMMAAC 0.25
#define GAMMABC 0.25 
//#define GAMMAABC 3.0

#define LS 20.0 //preto
#define LI 0.0 //vermelho
#define LL 0.0 //liquido
#define TmS 110.0//300.0
#define TmI 110.0//300.0
#define TmL 110.0//300.0// não pode dividir por zero burrao

#define T 100.0 //100.0
//#define cv 0.1

#define TAU 1.0 //era 1.0
#define EPSILON 1.0 //era 3.0
#define M 1.0//era 1.0
//PARÂMETROS DE GEOMETRIA
#define compL 45	//comprimento do domínio, era 75 para TELX=300
#define R 15		//raio inicial
#define TELX 300	//total de elementos na direção x 1250
#define TELY 300	//total de elementos na direção y 1250
//PARÂMETROS DE TEMPO
#define Ttot 480	//tempo total da simulação 1920
#define dt 0.01		//tamanho do intervalo entre os passos de tempo 0.02
//PARÂMETROS DE SIMULAÇÃO
#define EXISTE 1	//valor da variável de fase no estado sólido
#define N_EXISTE 0	//valor da variável de fase no estado líquido
//#define Ueq 0		//super-resfriamento de equilíbrio
//#define Ui 0.45		//super-resfriamento inicial 0.65 normalmente
//#define S 0.05		//grau de anisotropia
#define err 0.00001	//erro máximo da iteração da solução implicita 0.00001
