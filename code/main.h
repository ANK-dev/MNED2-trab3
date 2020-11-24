/******************************************************************************
 *        Métodos Numéricos para Equações Diferenciais II -- Trabalho 2       *
 *                          Ariel Nogueira Kovaljski                          *
 ******************************************************************************/

/*======================= Parâmetros a serem ajustados =======================*/

#define LX      (10.0)              /* Lx: comprimento do domínio (em m)      */
#define NX      (200)               /* nx: número de células                  */
#define DELTA_X (LX/NX)             /* Δx: largura de cada célula (em m)      */
#define U_BAR   (2.0)               /* ū: velocidade de escoamento (em m/s)   */
#define T_FINAL (1.0)               /* tempo final da simulação (em segundos) */
#define COURANT (0.9)               /* C: número de courant                   */

#define DELTA_T (COURANT*DELTA_X/U_BAR) /* Δt: passo de tempo (em segundos)   */

#define A (100.0)                   /*                                */
#define B (1.5)                     /*                                */
#define C (4.0)                     /* Parâmetros da condição inicial */
#define D (6.0)                     /*                                */
#define E (2.0)                     /*                                */

/*============================================================================*/

/* Métodos utilizados para o cálculo de Q neste trabalho */
enum methods {FTBS, LF, LW, BW};

void listParameters();
void initializeArray(double arr[], int len, double a, double b, double c,
                     double d, double e);
void calculateQ_FTBS(double old_arr[], double new_arr[]);
void calculateQ_LF(double old_arr[], double new_arr[]);
void calculateQ_LW(double old_arr[], double new_arr[]);
void calculateQ_BW(double old_arr[], double new_arr[]);
void printAndSaveResults(double arr[], int len, int method);