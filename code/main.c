/******************************************************************************
 *        Métodos Numéricos para Equações Diferenciais II -- Trabalho 3       *
 *                          Ariel Nogueira Kovaljski                          *
 ******************************************************************************/

#include "main.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int main(void)
{
    double Q_new[NX];   /* Array de Q no tempo n+1 */
    double Q_old[NX];   /* Array de Q no tempo n */

    puts("MNED II - Trabalho 3\n"
         "====================\n"
         "por Ariel Nogueira Kovaljski\n");

    listParameters();

    /* Cálculo de Q via método Upwind */
    puts("Calculando Q via metodo Upwind...");
    initializeArray(Q_old, NX, A, B, C, D, E);
    initializeArray(Q_new, NX, A, B, C, D, E);
    calculateQ(Q_old, Q_new, upwind);
    printAndSaveResults(Q_new, NX, UPWIND);

    /* Cálculo de Q via método Superbee */
    puts("Calculando Q via metodo Superbee...");
    initializeArray(Q_old, NX, A, B, C, D, E);
    initializeArray(Q_new, NX, A, B, C, D, E);
    calculateQ(Q_old, Q_new, superbee);
    printAndSaveResults(Q_new, NX, SUPERBEE);

    /* Cálculo de Q via método Van Albada */
    puts("Calculando Q via metodo Van Albada...");
    initializeArray(Q_old, NX, A, B, C, D, E);
    initializeArray(Q_new, NX, A, B, C, D, E);
    calculateQ(Q_old, Q_new, vanAlbada);
    printAndSaveResults(Q_new, NX, VAN_ALBADA);

    return 0;
}

void listParameters()
{
    puts("Parametros\n----------");
    puts("Constantes da equacao:");
    printf("DELTA_T = %f, DELTA_X = %f, U_BAR = %3.2e\n\n",
           DELTA_T, DELTA_T, U_BAR);
    puts("Constantes da simulacao:");
    printf("LX = %f, NX = %d, T_FINAL = %f\n\n", LX, NX, T_FINAL);
}

/* Inicializa um array para um valor de entrada */
void initializeArray(double arr[], int len, double a, double b, double c,
                     double d, double e)
{
    /*************************************************************
     *
     *                    Condição de contorno
     *                c(x,0) = exp(-A(x-B)²) + s(x)
     * 
     *           onde s(x) = { E,    se  C <= x <= D
     *                       { 0,    c.c.
     */

    int i;
    double x, s;

    for (i = 0; i < len; ++i) {
        x = i * DELTA_X;
        s = (x >= c && x <= d ? e : 0);
        arr[i] = exp( -a * ((x - b)*(x - b)) ) + s;
    }
}

double thetaPlusHalf(double arr[], int i) 
{
    int a, b, c;

    /* Aplica condição de contorno / volume fantasma */
    a = (i - 1) < 0 ? 0 : (i - 1);
    b = i;
    c = i + 1;

    /* Workaround para divisão por zero */
    if (arr[c] - arr[b] == 0) {
        return (arr[b] - arr[a]) / 1e-10;
        /* return 0; */
    }

    return (arr[b] - arr[a]) / (arr[c] - arr[b]);
}

double thetaMinusHalf(double arr[], int i) 
{
    int a, b, c;
    
    /* Aplica condição de contorno / volume fantasma */
    a = (i - 2) < 0 ? 0 : (i - 2);
    b = (i - 1) < 0 ? 0 : (i - 1);
    c = i;

    /* Workaround para divisão por zero */
    if (arr[c] - arr[b] == 0) {
        return (arr[b] - arr[a]) / 1e-10;
        /* return 0; */
    }

    return (arr[b] - arr[a]) / (arr[c] - arr[b]);
}

/* calcula Q, recebendo como ponteiro a função `psi` */
void calculateQ( double old[], double new[], double (*psi)(double theta) )
{
    int i;
    double t = 0;

    do {
        /* 
         *                     Fronteira esquerda 
         *                |=@=|=@=|=@=|=@=|=@=|=@=|=@=|
         *                  ^ 
         */
        i = 0;
        new[i] = old[i] - COURANT/2 * (1-COURANT) * (
                     psi(thetaPlusHalf(old, i)) * (old[i+1] - old[i])
                 );

        /* 
         *                  Após a fronteira esquerda 
         *                |=@=|=@=|=@=|=@=|=@=|=@=|=@=|
         *                      ^ 
         */
        i = 1;
        new[i] = old[i] - COURANT * (
                        old[i] - old[i-1]
                    ) - COURANT/2 * (1-COURANT) * (
                          psi(thetaPlusHalf(old, i))  * (old[i+1] - old[i]  ) 
                        - psi(thetaMinusHalf(old, i)) * (old[i]   - old[i-1])
                    );

        /* 
         *                            Centro
         *                |=@=|=@=|=@=|=@=|=@=|=@=|=@=|
         *                          ^   ^   ^   ^
         */
        for (i = 2; i < NX - 1; ++i) {
            new[i] = old[i] - COURANT * (
                         old[i] - old[i-1]
                     ) - COURANT/2 * (1-COURANT) * (
                           psi(thetaPlusHalf(old, i))  * (old[i+1] - old[i]  ) 
                         - psi(thetaMinusHalf(old, i)) * (old[i]   - old[i-1])
                     );
        }

        /* 
         *                      Fronteira direita
         *                |=@=|=@=|=@=|=@=|=@=|=@=|=@=|
         *                                          ^ 
         */
        new[i] = old[i] - COURANT * (
                     old[i] - old[i-1]
                 ) - COURANT/2 * (1-COURANT) * (
                     - psi(thetaMinusHalf(old, i)) * (old[i]   - old[i-1])
                 );

       /* Atualiza array de valores antigos com os novos para o próximo
          passo de tempo */
        for (i = 0; i < NX; ++i) {
            old[i] = new[i];
        }      

    } while ( (t += DELTA_T) <= T_FINAL);
}

/* Método Upwind */
double upwind(double theta)
{
    return (theta = 0);
}

/* Método Superbee */
double superbee(double theta)
{
    return MAX3(0, MIN(1, 2*theta), MIN(2, theta));
}

/* Método Van Albada */
double vanAlbada(double theta)
{
    return ( (theta*theta) + theta ) / ( (theta*theta) + 1 );
}

/* Imprime na tela e salva os resultados num arquivo de saída */
void printAndSaveResults(double arr[], int len, int method)
{
    int i;
    char filename[50], file0[100], file1[100];
    FILE *results_file;     /* Ponteiro para o arquivo de resultados */

    /* Prepara nome do arquivo de saída */
    switch (method) {
        case UPWIND:
            sprintf(filename, "%s%d%s", "results", method, ".txt");
            break;
        case SUPERBEE:
            sprintf(filename, "%s%d%s", "results", method, ".txt");
            break;
        case VAN_ALBADA:
            sprintf(filename, "%s%d%s", "results", method, ".txt");
            break;
    }

    /* Imprime os resultados no console */
    printf("Q[%d] (tempo final: %.2fs) = [", NX, T_FINAL);
    for (i = 0; i < len - 1; ++i) {
        printf("%f, ", arr[i]);
    }
    printf("%f]\n\n", arr[i]);

    /* Error Handling -- Verifica se é possível criar/escrever o arquivo de
                         resultados */
    sprintf(file0, "%s%s", "./results/",    filename);
    sprintf(file1, "%s%s", "./../results/", filename);
    if (   ( results_file = fopen(file0, "w") ) == NULL
        && ( results_file = fopen(file1, "w") ) == NULL ) {
        fprintf(stderr, 
                "[ERR] Houve um erro ao escrever o arquivo \"%s\"! "
                "Os resultados nao foram salvos.\n", filename);
        exit(1);
    }

    /* Adiciona os resultados no arquivo "results.txt" */
    fprintf(results_file,
            "nx=%d\n"
            "Delta_t=%f\n"
            "Delta_x=%f\n"
            "t_final=%f\n"
            "u_bar=%f\n",
            NX, DELTA_T, DELTA_X, T_FINAL, U_BAR);
    fputs("********************\n", results_file);

    for (i = 0; i < len; ++i) {
        fprintf(results_file, "%f,%f\n", i * DELTA_X, arr[i]);
    }

    fclose(results_file);   /* Fecha o arquivo */

    printf("[INFO] Os resultados foram salvos no arquivo \"%s\" "
           "no diretorio \"results/\".\n\n", filename);
}