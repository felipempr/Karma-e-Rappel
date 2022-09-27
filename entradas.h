/************************************************************************************************************************/
/*															*/
/*					ARQUIVO DE ENTRADAS DO CÓDIGO							*/
/*															*/
/************************************************************************************************************************/

//CONSTANTES
#define M_PI 3.14159265358979323846
#define constK 1.0
#define cv 1.0
//#define h 3000
//#define E0 1.0//W0 normalmente 1.0
#define D 3.0
//#define a1 0.8839
//#define a2 0.6267

//#define GAMMASL 0.001 //era 0.25 era 0.001
//#define GAMMASI 0.001
//#define GAMMALI 0.001
//#define GAMMAABC 0.001

#define MISL 10.0 //era 0.1
#define MISI 10.0
#define MILI 10.0

#define EPSILONSL 0.0005
#define EPSILONSI 0.0005
#define EPSILONLI 0.0005

#define WSL 0.018
#define WSI 0.018
#define WLI 0.018

#define LS 5.0 //preto
#define LI 0.0 //vermelho
#define LL 0.0 //liquido MUDEI
#define TmS 2.0//300.0
#define TmI 1.0//300.0
#define TmL 2.0//300.0// não pode dividir por zero burrao

#define T 1.0 //100.0
//#define cv 0.1

#define TAU 0.2 //era 1.0
#define EPSILON 0.05 //era 3.0

//PARÂMETROS DE GEOMETRIA
#define compL 1.0	//comprimento do domínio, era 75 para TELX=300
#define R 15		//raio inicial
#define TELX 100	//total de elementos na direção x 1250
#define TELY 100	//total de elementos na direção y 1250
//PARÂMETROS DE TEMPO
#define Ttot 1	//tempo total da simulação 1920
#define dt 0.00001		//tamanho do intervalo entre os passos de tempo 0.02
//PARÂMETROS DE SIMULAÇÃO
#define EXISTE 0.99999	//valor da variável de fase no estado sólido
#define N_EXISTE 0.000005	//valor da variável de fase no estado líquido
//#define Ueq 0		//super-resfriamento de equilíbrio
//#define Ui 0.45		//super-resfriamento inicial 0.65 normalmente
//#define S 0.05		//grau de anisotropia
#define err 0.000001	//erro máximo da iteração da solução implicita 0.00001
