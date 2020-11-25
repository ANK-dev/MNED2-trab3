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
    initializeArray(Q_old, NX, A, B, C, D, E);
    initializeArray(Q_new, NX, A, B, C, D, E);
    calculateQ(Q_old, Q_new, upwind);
    printAndSaveResults(Q_new, NX, UPWIND);

    /* Cálculo de Q via método Superbee */
    initializeArray(Q_old, NX, A, B, C, D, E);
    initializeArray(Q_new, NX, A, B, C, D, E);
    calculateQ(Q_old, Q_new, superbee);
    printAndSaveResults(Q_new, NX, SUPERBEE);

    /* Cálculo de Q via método Van Albada */
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
    printf("NX = %d, T_FINAL = %f\n\n", NX, T_FINAL);
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

/* calcula Q, recebendo como ponteiro a função `psi`*/
void calculateQ( double old[], double new[], double (*psi)(double theta) )
{
    int i;
    double t = 0, theta_minus_half, theta_plus_half;

    do {
        /* Fronteira esquerda */
        theta_plus_half = 0;

        new[0] = old[0] - 1/2 * C * (1-C) * (
                     psi(theta_plus_half) * (old[2] - old[1])
                 );

        /* Centro */
        for (i = 1; i < NX; ++i) {
            theta_minus_half = (old[i-1] - old[i-2]) / (old[i]   - old[i-1]);
            theta_plus_half  = (old[i]   - old[i-1]) / (old[i+1] - old[i]  );

            new[i] = old[i] - C * (
                         old[i] - old[i-1]
                     ) - 1/2 * C * (1-C) * (
                           psi(theta_plus_half)  * (old[i+1] - old[i]  ) 
                         - psi(theta_minus_half) * (old[i]   - old[i-1])
                     );
        }

        /* Fronteira direita */
        theta_minus_half = (old[i-1] - old[i-2]) / (old[i] - old[i-1]);

        new[i] = old[i] - C * (
                     old[i] - old[i-1]
                 ) - 1/2 * C * (1-C) * (
                     - psi(theta_minus_half) * (old[i]   - old[i-1])
                 );
        
    } while ( (t += DELTA_T) < T_FINAL);
}


double upwind(double theta)
{
    return (theta = 0);
}

/* Método Lax-Friedrichs */
double superbee(double theta)
{
    return MAX3(0, MIN(1, 2*theta), MIN(2, theta));
}

/* Método Lax-Wendroff */
double vanAlbada(double theta)
{
    return ( (theta * theta) + theta ) / ( (theta * theta) + 1 );
}

/* Imprime na tela e salva os resultados num arquivo de saída */
void printAndSaveResults(double arr[], int len, int method)
{
    int i;
    char filename[50], file0[100], file1[100];
    FILE *results_file;     /* Ponteiro para o arquivo de resultados */

    /* Prepara nome do arquivo de saída */
    /* ATENÇÃO: verificar se não há buffer overflow! */
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
    printf("\n\nQ[%d] (tempo final: %.2fs) = [", NX, T_FINAL);
    for (i = 0; i < len - 1; ++i) {
        printf("%f, ", arr[i]);
    }
    printf("%f]\n\n", arr[i]);

    /* Error Handling -- Verifica se é possível criar/escrever o arquivo de
                         resultados */
    /* ATENÇÃO: verificar se não há buffer overflow! */
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