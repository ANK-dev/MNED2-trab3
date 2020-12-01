/******************************************************************************
 *        Métodos Numéricos para Equações Diferenciais II -- Trabalho 3       *
 *                          Ariel Nogueira Kovaljski                          *
 ******************************************************************************/

/*======================= Parâmetros a serem ajustados =======================*/

#define LX      (20.0)              /* Lx: comprimento do domínio (em m)      */
#define NX      (400)               /* nx: número de células                  */
#define DELTA_X (LX/NX)             /* Δx: largura de cada célula (em m)      */
#define U_BAR   (1.0)               /* ū: velocidade de escoamento (em m/s)   */
#define T_FINAL (14.0)               /* tempo final da simulação (em segundos) */
#define COURANT (0.8)               /* C: número de courant                   */

#define DELTA_T (COURANT*DELTA_X/U_BAR) /* Δt: passo de tempo (em segundos)   */

#define A (100.0)                   /*                                */
#define B (1.5)                     /*                                */
#define C (4.0)                     /* Parâmetros da condição inicial */
#define D (6.0)                     /*                                */
#define E (2.0)                     /*                                */

/*============================================================================*/

#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MAX3(a,b,c) (MAX(MAX((a),(b)), (c)))

/* Métodos utilizados para o cálculo de Q neste trabalho */
enum methods {UPWIND, SUPERBEE, VAN_ALBADA};

void listParameters();
void initializeArray(double arr[], int len, double a, double b, double c,
                     double d, double e);
double thetaPlusHalf(double arr[], int i);
double thetaMinusHalf(double arr[], int i);
void calculateQ(double old[], double new[], double (*psi)(double theta));
double upwind(double theta);
double superbee(double theta);
double vanAlbada(double theta);
void printAndSaveResults(double arr[], int len, int method);