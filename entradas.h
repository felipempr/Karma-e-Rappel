/************************************************************************************************************************/
/*															*/
/*					ARQUIVO DE ENTRADAS DO CÓDIGO							*/
/*															*/
/************************************************************************************************************************/

//CONSTANTES
#define M_PI 3.14159265358979323846
#define constK 200.0//1.0
#define cv 900.0//1.0
//#define h 3000
//#define E0 1.0//W0 normalmente 1.0
//#define D 3.0
//#define a1 0.8839
//#define a2 0.6267

//#define GAMMASL 0.001 //era 0.25 era 0.001
//#define GAMMASI 0.001
//#define GAMMALI 0.001
//#define GAMMAABC 0.001

#define MISL 0.0000000010//20.0 
#define MISI 0.0//4.0
#define MILI 0.0//0.0000000000010//20.0

//#define EPSILONSL 0.00000000000002//0.00000000000001
//#define EPSILONSI 0.00000000000002//0.00000000000002
//#define EPSILONLI 0.0000000000002//0.00000000000002

#define LAMBDASL 0.000000001
#define LAMBDASI 0.000000001
#define LAMBDALI 0.000000001

#define SIGMASL 0.0001
#define SIGMASI 0.000023
#define SIGMALI 0.0001

//#define WSL 45000.0//0.018
//#define WSI 45000.0//0.018
//#define WLI 450000.0//0.018

#define LSL 400000.0
#define LSI 0.0
#define LLI 0.0

#define TmSL 600.0
#define TmSI 1.0
#define TmLI 1.0

//#define LS 400000.0 //preto 1.0
//#define LI 0.0 //vermelho
//#define LL -400000.0//-400000.0 //liquido MUDEI
//#define TmS 600.0//600.0
//#define TmI 1.0//1.0
//#define TmL 600.0//1.0// não pode dividir por zero burrao

#define T 550.0 //500.0
//#define cv 0.1

//#define TAU 0.2 //era 1.0
//#define EPSILON 0.05 //era 3.0

//PARÂMETROS DE GEOMETRIA
#define compL 0.0000001//0.0000001
#define R 15		//raio inicial
#define TELX 1500	//total de elementos na direção x 1250
#define TELY 1500	//total de elementos na direção y 1250
//PARÂMETROS DE TEMPO
#define Ttot 10000	//tempo total da simulação 250000
#define dt 0.4 //0.001		//tamanho do intervalo entre os passos de tempo 0.02
//PARÂMETROS DE SIMULAÇÃO
#define EXISTE 1.0//0.99999	//valor da variável de fase no estado sólido
#define N_EXISTE 0.0//0.000005	//valor da variável de fase no estado líquido
//#define Ueq 0		//super-resfriamento de equilíbrio
//#define Ui 0.45		//super-resfriamento inicial 0.65 normalmente
//#define S 0.05		//grau de anisotropia
#define err 0.00001	//erro máximo da iteração da solução implicita 0.00001
