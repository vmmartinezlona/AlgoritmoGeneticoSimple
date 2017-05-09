/*
    Simple Genetic Algorith
    Artficial Inteligent, DICIS. 2017.
    April 2017

    Funcion F6

    Galvan Hernandez Armando
    Martinez Lona Veronica Montserrat

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// Variable perteneciente a la ecuacion
#define V 4189.829101

typedef unsigned int uint;

// Number of genes == number of variables of the equation
const uint NUMGENS = 10;
// Gen lenght (for simplicity and re-use all the gens has the same length )
const uint BITGEN = 9; // Para poder almacenar el numero 500
// Number of indiduals in the populations(has to be pair)
// 10, 20, 30, 50, 100, 200
const uint NUMCHROMOSOME = 20;
// Numero maximo de iteraciones
const uint NUMMAXIT = 20;
// Search space limits
const int LOWERLIM = -500;
const int UPPERLIM = 500;
//Probabilidad de mutacion 0.001 - 0.009
const double pMuta = 0.005;
//Probabilidad de cruza
const double pCross = 0.8;

typedef struct {
    char* gen;               // Arreglo que contiene al gen, bit por bit.
    // unsigned int bitgen;  // Longitud del gen.
    //double fenotype;         // Guarda el fenotipo, lo que es el valor numerico del gen completo
}GEN;

typedef struct {
    //unsigned int numGens; // Numero de genes
    GEN* genes;             // Conjunto de los genes que conforma el cromosoma
    double chromFit;        // Fitness del gen
    double bestFit;        // Mejor fitness historico
} CHROMOSOME;

//Estructura POPULATION, es para el manejo de la poblacion de cromosomas
typedef struct {
    CHROMOSOME* chromosomes; // Arreglo de cromosomas que conforman la ecuacion
    int bestChromosome;   // Index del cromosoma con mejor fitness
    // double bestGFit;        // Mejor fitness historico
} POPULATION;

// Estructura ARRAY, es auxiliar en la funcion Roulette_Metod.
// Guarda el valor del index (index) y el valor de la probabilidad de seleccion (value)
typedef struct {
    uint index;     // Index del cromosoma
    double value;   // Valor de probabilidad de seleccion
} ARRAY;

//____________________________________ Init
// Inicializar gen
void Gen_Init(GEN *pGen);
// Inicializar chromosoma
void Chromosome_Init(CHROMOSOME *pChromosome);
// Inicalizar poblacion
POPULATION* Population_Init(void);
//____________________________________ Evaluation
// Decodificar gen
int Decode_Gen(GEN* pGen);
// Obtener el fenotipo
double Get_Fenotype(GEN* pGen, int binMax, int range);
// Evaluar individuo
void Evaluate_Chromosome(CHROMOSOME* pChromosome, int binMax, int range);
//____________________________________ Selection
// Compute the totalFit
int Total_Fit(POPULATION *pPopulation);
// seleccion probability
ARRAY* Selection_Probability(POPULATION* pPopulation, int totalFit);
// For quickSort function
int cmpfunc (const void * a, const void * b);
// Roulette Method
int* Roulette_Metod(POPULATION *pPopulation);
//____________________________________ Evolucion
//Operador cruza
POPULATION* Cross_Population(POPULATION *pPopulation, int *selPair);
//Mutación
void Muta_Population(CHROMOSOME *pChromosome);
//
void UpdateBest(POPULATION *pPopulation);
//____________________________________ Finish
// Show poblacion
void Show_Population(POPULATION *pPopulation);
// Liberar poblacion
void Free_Population(POPULATION *pPopulation);

int main(void)
{
    srand(time(NULL));
    uint i, binMax, range, it;
    int bestChromosome;
    int *selPair;
    double error = 0.0001;
    POPULATION *pPopulation;

    // Inicalizamos poblacion
    pPopulation = Population_Init();

    binMax = 2 << (BITGEN - 1);
    range = UPPERLIM - LOWERLIM;

    for(it = 0; it < NUMMAXIT; it++)
    {
        // Evaluamos la poblacion
        for(i = 0; i < NUMCHROMOSOME; i++)
        {
            Evaluate_Chromosome(&(pPopulation -> chromosomes[i]), binMax, range);
        }
        UpdateBest(pPopulation);
        Show_Population(pPopulation);
        bestChromosome = pPopulation -> bestChromosome;
        printf("bestChromosome : %d , bestFitness : %f\n", bestChromosome, pPopulation -> chromosomes[bestChromosome].bestFit);
        printf("----------------------------Iteracion %d\n", it + 1);
        if(pPopulation -> chromosomes[bestChromosome].bestFit - 10 * V < error) break;
        // Seleccion
        selPair = Roulette_Metod(pPopulation);
        //Curza
        pPopulation = Cross_Population(pPopulation, selPair);
        // Mutacion
        for(i = 0; i < NUMCHROMOSOME; i++)
        {
            Muta_Population(&(pPopulation -> chromosomes[i]));
        }
    }

    Free_Population(pPopulation);

    return 0;
}


/*
    Inicializa cada bit del gen.
    Le asigna al gen un valor aleatorio (1 o 0) bit por bit.
    Inputs:
        Un gen vacio.
    Output:
        void
*/
void Gen_Init(GEN *pGen)
{
    uint i;

    pGen -> gen = (char*) malloc (BITGEN*sizeof(char));

    // Para cada bit del gen
    for(i = 0; i < BITGEN; i++)
    {
        // Alatoriamente asigna 1 o 0
        pGen -> gen[i] = (char)round((double) rand() / RAND_MAX);
    }

    return;
}

/*
    Inicializa el cromosoma
    Inicializa un cromosoma y su conjunto de genes .
    Tambien se le asigna un fitness igual a 0.
    Inputs:
        Un cromosoma
    Outputs:
        void
*/
void Chromosome_Init(CHROMOSOME *pChromosome)
{
    uint i;

    // El fitness se inicializa en 0
    pChromosome -> chromFit = 0;
    pChromosome -> bestFit = 0;

    //Se crea el conjunto de genes
    pChromosome -> genes = (GEN*) malloc (NUMGENS * sizeof(GEN));

    // Para cada gen
    for(i = 0; i < NUMGENS; i++)
    {
        //Se inicializa cada gen del cromosoma
        Gen_Init(&(pChromosome -> genes[i]));
    }

    return;
}

/*
    Inicializa la poblacion.
    Inicializa una poblacion y sus cromosamas.
    La variable de bestChromosome es igualada a 0.

    Inputs:
        Una poblacion vacia
    Outputs
        Una poblacion inicializada
*/
POPULATION* Population_Init(void)
{
    uint i;
    POPULATION* pPopulation;

    // Se reserva la memoria para la poblacion
    pPopulation = (POPULATION*) malloc (sizeof(POPULATION));

    // Por default, el mejor cromosoma es el 0
    pPopulation -> bestChromosome = 0;

    // Se reserva el espacio para el arreglo de cromosomas
    pPopulation -> chromosomes = (CHROMOSOME*) malloc (NUMCHROMOSOME * sizeof(CHROMOSOME));

    for(i = 0; i < NUMCHROMOSOME; i++)
    {
        // Se inicaliza cada cromosoma de la poblacion
        Chromosome_Init(&(pPopulation -> chromosomes[i]));
    }

    return pPopulation;
}

/*
    Decodificar el gen
    Calcula el valor numerico de cada gen, conviertiendo el conjunto binario a decimal

    Inputs:
        El gen a decodificar
    Outputs:
        El valor numerico calculado
*/
int Decode_Gen(GEN* pGen)
{
    int i;
    uint base = 1;
    uint number = 0;

    // Desde el bit menos significativo hasta el mas significativo
    // Por eso va "al reves" el conteo
    for(i = BITGEN - 1; i >= 0; i --, base *= 2)
    {
        number += pGen -> gen[i] * base;
    }

    return number;
}

/*
    Obtiene el fenotipo.
    Calcula el valor numerico de punto flotante que guarda cada gen
    Inputs:
        El gen del cual se obtendra el fenotipo.
        binMax: el numero maximo que se puede tener con ese numero de bits
        range: el rango donde se encuenra la solucion
    Outputs:
        El valor numerico calculado
*/
double Get_Fenotype(GEN* pGen, int binMax, int range)
{
    uint decodedGen;
    double fenotype;

    decodedGen = Decode_Gen(pGen);
    fenotype = (double)decodedGen / binMax * range + LOWERLIM;

    return fenotype;
}

/*
    Evalua el cromosoma en la ecuacion
    Fx = 10V + sum[-xi * sin(sqrt(abs(xi)))]
    Inputs:
        El cromosoma a evaluar
        binMax: el numero maximo que se puede tener con ese numero de bits
        range: el rango donde se encuenra la solucion
    Outputs:
        void
*/
void Evaluate_Chromosome(CHROMOSOME* pChromosome, int binMax, int range)
{
    uint i;
    double x;
    double F = 0;

    //Para cada gen
    for(i = 0; i < NUMGENS; i++)
    {
        x = Get_Fenotype(&(pChromosome -> genes[i]), binMax, range);
        // Si el fenotipo de un se sale de los limites asignados entonces el cromosoma no es una solucion
        if(x > UPPERLIM || x < LOWERLIM)
        {
            pChromosome -> chromFit = 0;
            return;
        }
        // Si el fenotipo tiene un valor correcto
        else
        {
            F += -1 * x * sin(sqrt(fabs(x)));;
        }
    }

    pChromosome -> chromFit = 10 * V + F;

    return;
}

/*
    Calcula el fitness total.
    Suma los fitness de cada uno de los cromosomas.
    Inputs:
        La poblacion
    Outputs:
        Valor del fitness total
*/
int Total_Fit(POPULATION *pPopulation)
{
    uint i;
    double totalFit = 0;

    //Por cada cromosoma
    for(i = 0; i < NUMCHROMOSOME; i++)
    {
        totalFit += pPopulation -> chromosomes[i].chromFit;
    }

    return totalFit;

}

/*
    Obtiene la probabilidad de seleccion de cada cromosoma
    Divide el fitness de cada cromosma entre el fitness total
    Inputs:
        La poblacion
        El fitness total
    Outputs:
        Un arreglo del tipo ARRAY, que contiene el valor de la probabilidad de seleccion y el indice del comosoma
*/
ARRAY* Selection_Probability(POPULATION* pPopulation, int totalFit)
{
    uint i;
    ARRAY *pSelProb;

    pSelProb  = (ARRAY*)  malloc (NUMCHROMOSOME * sizeof(ARRAY));

    for(i = 0; i < NUMCHROMOSOME; i++)
    {
        pSelProb[i].value = pPopulation -> chromosomes[i].chromFit / totalFit;
        pSelProb[i].index = i;
    }

    return pSelProb;
}

/*
    Funcion de comparacion.
    Regresa el elemento del tipo ARRAY cuyo .value sea mayor.
    Es para poder implementar la Funcion qsort de la libreria stdlib
    Inputs:
        Dos elemtos del arreglo
    Outputs:
        un valor entero
*/
int cmpfunc (const void * a, const void * b)
{
    const ARRAY *x;
    const ARRAY *y;
    x = a;
    y = b;

    if (x -> value < y -> value) return -1;
    else if (x -> value > y -> value) return 1;
    return 0;
}

/*
    Aplica el metodo de la ruleta
    Cada cromosoma tiene una probabilidad de seleccion,
    se calcula un numero random con el cual se selecciona un cromosoma y se van armando parejas
    Inputs:
        La poblacion
    Outputs:
        Un arreglo tipo entero con los pares de cromosomas que se cruzaran
*/
int* Roulette_Metod(POPULATION *pPopulation)
{
    uint i, j;
    double rNum, sumAux;
    ARRAY *pSelProb;
    int *pSelPair; //Guarda los indices de los cromosomas que seran pareja
    pSelPair = (int*) malloc (NUMCHROMOSOME * sizeof(int));

    // Se calcula la probabilida de seleccion de cada cromosoma de la poblacion
    pSelProb = Selection_Probability(pPopulation, Total_Fit(pPopulation));
    // Se ordenan de menor a mayor la probabilidad de seleccion de cada cromosoma
    qsort(pSelProb, NUMCHROMOSOME, sizeof(ARRAY), cmpfunc);

    for(i = 0; i < NUMCHROMOSOME; i++)
    {
        // Valor random que nos dice que parte de la ruleta se debe escoger
        rNum = (double) rand() / RAND_MAX;
        for(j = 0, sumAux = 0; j < NUMCHROMOSOME; j++)
        {
            //Va acumulando las probabilidades para avanzar de posicion en la ruleta
            sumAux += pSelProb[j].value;
            if(sumAux > rNum)
            {
                pSelPair[i] = pSelProb[j].index;
                break;
            }
        }
    }

    // Libera la memoria del arreglo de probabilidades
    free(pSelProb);

    return pSelPair;
}

/*
//----------------------------------------------------Cruza poblacion
Hace la cruza con 1 punto de cruza. (Intercambia bits entre dos cromosomas)
    Inputs:
        Poblacion
        Vector de parejas
    Output:
        Nueva generacion
*/
POPULATION* Cross_Population(POPULATION *pPopulation, int *selPair){
	int i, j, k, nCouples, nBits, bitCounter, crossPoint;
	double vCross; // Valor de cruza, determina si la pareja se cruzara o no
	POPULATION *newPopulation; // Guarda la nueva generacion

    // Numero de parejas, es la mitad del numero de cromosomas
	nCouples = NUMCHROMOSOME * 0.5;
    // Numero de bits totales del cromosoma
	nBits = NUMGENS * BITGEN;

    // Creamos e inicializamos la nueva poblacion
	newPopulation = Population_Init();

    // Para cada pareja
	for(i = 0; i < nCouples; i++){
		vCross = (double)rand()/RAND_MAX;
        // Si el valor de cruza es menor a la probabilidad de cruza hay cruza
		if(vCross <= pCross){
			bitCounter = 0;
            // El punto de cruza es aleatorio py diferente para cada par de parejas
			crossPoint = (int)(((double)rand() / RAND_MAX * (nBits - 1)) + 1);
            // Para cada gen
			for(j = 0; j < NUMGENS; j++){
                // Para cada bit
				for (k = 0; k < BITGEN; k++){
					if(bitCounter < crossPoint){ //Copia de forma normal hasta cierto punto
						newPopulation->chromosomes[2*i].genes[j].gen[k] = pPopulation->chromosomes[selPair[2*i]].genes[j].gen[k];
						newPopulation->chromosomes[(2*i)+1].genes[j].gen[k] = pPopulation->chromosomes[selPair[(2*i)+1]].genes[j].gen[k];
					}
					else{ // Hace la cruza (copia cruzada)
						newPopulation->chromosomes[2*i].genes[j].gen[k] = pPopulation->chromosomes[selPair[(2*i)+1]].genes[j].gen[k];
						newPopulation->chromosomes[(2*i)+1].genes[j].gen[k] = pPopulation->chromosomes[selPair[2*i]].genes[j].gen[k];
					}
                    // Para saber si llegamos al punto de cruza
					bitCounter++;
				}
			}
		}
		else{	//No hubo cruza, se pasan los padres igual.
            // Para cada gen
			for(j = 0; j < NUMGENS; j++){
                // Para cada bit
				for (k = 0; k < BITGEN; k++){
					newPopulation->chromosomes[2*i].genes[j].gen[k] = pPopulation->chromosomes[selPair[2*i]].genes[j].gen[k];
					newPopulation->chromosomes[(2*i)+1].genes[j].gen[k] = pPopulation->chromosomes[selPair[(2*i)+1]].genes[j].gen[k];
				}
			}
		}
	}

    free(selPair);
	Free_Population(pPopulation);//Libera poblacion anterior

	return newPopulation;
}

/*
//----------------------------------------------------Mutar poblacion
    Cambia el valor de un bit por su complemento.

    Inputs:
        Poblacion
    Output:
        void
*/
void Muta_Population(CHROMOSOME *pChromosome){
	int i, j;
	double vMuta; // Valor de mutacion, define si habra o no muta

    // Para cada gen
	for(i = 0; i < NUMGENS; i++){
        // Para cada bit
		for(j = 0; j < BITGEN; j++){
			vMuta = (double)rand()/RAND_MAX;
            // Si el valor de muta es menor a la probabilidad de muta hay mutacion
			if(vMuta <= pMuta){
                // Calculamos el complemento mediante modulos
				pChromosome -> genes[i].gen[j] = (pChromosome -> genes[i].gen[j] + 1) % 2;
			}
		}
	}

    return;
}

void UpdateBest(POPULATION *pPopulation)
{
  unsigned int i;
  //Peso del mejor cromosoma de la poblacion
  float best = pPopulation -> chromosomes[pPopulation -> bestChromosome].bestFit;

  //Para todos los cromosomas
  for(i = 0; i < NUMCHROMOSOME; i++)
  {
    //Si el fitness de actual es mayor que el del historico de ese chromosoma
    if(pPopulation -> chromosomes[i].chromFit > pPopulation -> chromosomes[i].bestFit)
    {
      pPopulation -> chromosomes[i].bestFit = pPopulation -> chromosomes[i].chromFit;
    }

    // Si el mejor peso actual de la particula es mejor que el de la mejor del ejambre
    if(pPopulation -> chromosomes[i].bestFit > best)
    {
      //Se actualiza el id y peso del mejor
      pPopulation -> bestChromosome = i;
      best = pPopulation -> chromosomes[i].chromFit;
    }
  }

  return;
}

/*
    Imprime la poblacion completa y otros datos de interes.
    Imprime el gen en binario y el fitness del cromosoma
    Inputs:
        La poblacion
    Outputs:
        void
*/
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
        }
        printf("\nchromosome's Fitness: %f", pPopulation -> chromosomes[i].chromFit);
        printf("\n-----\n");
    }

    return;
}

/*
    Libera el espacio de memoria que ocupa la poblacion
    Inputs:
        La poblacion
    Outputs:
        void
*/
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
