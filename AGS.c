/*
    Simple Genetic Algorith
    Artficial Inteligent, DICIS. 2017.
    April 2017

    Funcion F6

    Martinez Lona Veronica Montserrat

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define V 4189.829101

typedef unsigned int uint;

//Number of genes == number of variables of the equation
const uint NUMGENS = 10;
//Gen lenght (for simplicity and re-use all the gens has the same length )
const uint BITGEN = 9; // Para poder almacenar el numero 500
//Number of indiduals in the populations(has to be pair)
const uint NUMCHROMOSOME = 10;
//Limite inferior y superior del espacio de busqueda
const int LOWERLIM = -500;
const int UPPERLIM = 500;

typedef struct {
    char* gen; //Arreglo que contiene al gen
    //unsigned int bitgen; // Longitud del gen
    double fenotype;
}GEN;

typedef struct {
    //unsigned int numGens; //Numero de genes
    GEN* genes;
    double chromFit; //Fitness del gen
    double selProb; //Probabilidad de seleccion
} CHROMOSOME;

typedef struct {

    CHROMOSOME* chromosomes;
    //double bestFit;
    double bestChromosome;
} POPULATION;
//____________________________________ Init
//Inicializar gen
void Gen_Init(GEN *pGen);
//Inicializar chromosoma
void Chromosome_Init(CHROMOSOME *pChromosome);
//Inicalizar poblacion
POPULATION* Population_Init(void);
//____________________________________
//Decodificar gen
int Decode_Gen(GEN* pGen);
//Obtener el fenotipo
double Get_Fenotype(GEN* pGen, int binMax, int range);
//Evaluar individuo
void Evaluate_Chromosome(CHROMOSOME* pChromosome, int binMax, int range);
//seleccion

//Show poblacion
void Show_Population(POPULATION *pPopulation);
//Liberar poblacion
void Free_Population(POPULATION *pPopulation);

int main(void)
{
    srand(time(NULL));
    uint i, j, binMax, range;

    POPULATION *pPopulation;

    pPopulation = Population_Init();

    binMax = 2 << (BITGEN - 1);
    range = UPPERLIM - LOWERLIM;

    for(i = 0; i < NUMCHROMOSOME; i++)
    {
        Evaluate_Chromosome(&(pPopulation -> chromosomes[i]), binMax, range);
    }

    Show_Population(pPopulation);



    Free_Population(pPopulation);

    return 0;
}

void Gen_Init(GEN *pGen)
{
    uint i;

    pGen -> gen = (char*) malloc (BITGEN*sizeof(char));

    for(i = 0; i < BITGEN; i++)
    {
        pGen -> gen[i] = (char)round((double) rand() / RAND_MAX);
    }

    return;
}

void Chromosome_Init(CHROMOSOME *pChromosome)
{
    uint i;
    pChromosome -> chromFit = 0;
    pChromosome -> selProb = 0;

    pChromosome -> genes = (GEN*) malloc (NUMGENS * sizeof(GEN));

    for(i = 0; i < NUMGENS; i++)
    {
        Gen_Init(&(pChromosome -> genes[i]));
    }

    return;
}

POPULATION* Population_Init(void)
{
    uint i;
    POPULATION* pPopulation;

    pPopulation = (POPULATION*) malloc (sizeof(POPULATION));

    pPopulation -> bestChromosome = 0;
    pPopulation -> chromosomes = (CHROMOSOME*) malloc (NUMCHROMOSOME * sizeof(CHROMOSOME));

    for(i = 0; i < NUMCHROMOSOME; i++)
    {
        Chromosome_Init(&(pPopulation -> chromosomes[i]));
    }

    return pPopulation;
}

int Decode_Gen(GEN* pGen)
{
    uint i;
    uint base = 1;
    uint number = 0;

    for(i = BITGEN - 1; i >= 0; i --, base *= 2)
    {
        number += pGen -> gen[i] * base;
    }

    return number;
}

double Get_Fenotype(GEN* pGen, int binMax, int range)
{
    uint decodedGen;
    double fenotype;

    decodedGen = Decode_Gen(pGen);
    //printf("\n\tgen = %d\n", decodedGen);

    fenotype = (double)decodedGen / binMax * range + LOWERLIM;

    //printf("\n\t Fenotype = %f\n", fenotype);
    return fenotype;
}

void Evaluate_Chromosome(CHROMOSOME* pChromosome, int binMax, int range)
{
    uint i;
    double x;
    double F = 0;

    //Para cada gen
    for(i = 0; i < NUMGENS; i++)
    {
        x = Get_Fenotype(&(pChromosome -> genes[i]), binMax, range);
        if(x > 500 || x < -500)
        {
            pChromosome -> chromFit = 0;
            break;
        }

        else
        {

            F += -1 * x * sin(sqrt(fabs(x)));;
        }
    }

    pChromosome -> chromFit = F + (10 * V);

    return;
}




void Show_Population(POPULATION *pPopulation)
{
    uint i, j, k;

    //Para cada cromosoma
    for(i = 0; i < NUMCHROMOSOME; i++)
    {
        printf(" Individuo %3d", i + 1);
        //Para cada gen
        for(j = 0; j < NUMGENS; j++)
        {
            printf("\n\tGen %3d :  ", j + 1);
            //Para cada bit
            for(k = 0; k < BITGEN; k++)
            {
                printf("%d",
                    pPopulation -> chromosomes[i].genes[j].gen[k]);
            }

            //Get_Fenotype(&(pPopulation -> chromosomes[i].genes[j]), 512, 1000);
        }
        printf("\nchromosome's Fitness: %f", pPopulation -> chromosomes[i].chromFit);
        printf("\n-----\n");
    }

    return;
}

void Free_Population(POPULATION *pPopulation)
{
    uint i, j;

    for(i = 0; i < NUMCHROMOSOME; i++)
    {
        for(j = 0; j < NUMGENS; j++)
        {
            free(pPopulation -> chromosomes[i].genes[j].gen);
        }
    }

    for(i = 0; i < NUMCHROMOSOME; i++)
    {
        free(pPopulation -> chromosomes[i].genes);
    }

    free(pPopulation -> chromosomes);
    free(pPopulation);

    return;
}
