#include "cs-grn.h"
#include "ReadInputs.h"
#include "Cluster.h"
#include <thread>
using namespace std;

double abs(double x) {
    if (x < 0) return -x;
    return x;
}


int main (int argc, char * const argv[]) {
    
    if (argc < 7)
    {
        printf("usage: %s <numclusters> <numtrials> <couplingConstant> <orthfile> <numSpc> <spc1nw> <spc2nw> ... -s <startclustering> -t <maxtemp> -n <no arg, use if noise cluster needed>\n", argv[0]);
        exit(1);
    }
    
    clock_t begin = clock();
    /*====== read parameters =======*/
    int argbase = 1; // count for argv
    int numclusters = atoi(argv[argbase++]); // numclusters
    numclusters +=1; // add one more cluster: 'noise cluster 0'
    int numTrials = atoi(argv[argbase++]); // number of trials
    double couplingConstant = atof(argv[argbase++]); // couplingConstant
    char *orthfile = argv[argbase++]; // ortholog files
    int numSpc = atoi(argv[argbase++]); // number of species
    
    /*====== build network =========*/
    SpeciesNetwork **sns = new SpeciesNetwork*[numSpc]; // create a arr of size numSpc for network object
    for (int i=0; i<numSpc; i++) {
        char *spc = argv[argbase++]; // take a network file of a species
        sns[i] = ReadNetworkFile(spc); // initialize network for a species
    }
    
    /*====== read parameters: startclusfn '-s' and maxtemp '-t' =======*/
    char startclusfn[1024]; startclusfn[0] = 0;
    double maxtemp = -1;
    bool has_noise_cluster = false;
    while (argbase < argc)
    {
        // ground-truth clutering
        if (!strcmp(argv[argbase], "-s"))
        {
            argbase++;
            strcpy(startclusfn, argv[argbase++]);
            continue;
        }
        // starting tempreture
        if (!strcmp(argv[argbase], "-t"))
        {
            argbase++;
            maxtemp = atof(argv[argbase++]);
            continue;
        }
        // has noise cluster
        if (!strcmp(argv[argbase], "-n"))
        {
            argbase++;
            has_noise_cluster = true;
            continue;
        }
    }
    /*======= clustering ==========*/
    // create ortholog obj using (ortholog file, network arr, number of spe)
    Orthology *orth = ReadOrthologyFile(orthfile, sns, numSpc);
    // create arr of size numTrials for clustering obj
    Cluster **c = new Cluster *[numTrials];
    
    for (int i=0; i<numTrials; i++)
    {
        // initialize clustering parameters
        c[i] = new Cluster(sns, numSpc, orth, numclusters, couplingConstant, has_noise_cluster);
        
        // if ground-truth used, use ground-truth clustering label as start point
        if (startclusfn[0]!=0)
        {
            fprintf(stderr,"Using seed clustering file %s\n",startclusfn);
            c[i]->Preset(startclusfn);
        }
        
        // if starting tempreture used, set maxtemp
        if (maxtemp >= 0) c[i]->SetMaxTemp(maxtemp);
        
        // simulated annealing to find parameters with best cost
        c[i]->LearnGroundState();
    }
    
    /*==== go through all trails, and record the best cost ====*/
    double bestCost = 0.0;
    int bestTrial = -1;
    
    for (int i=0; i<numTrials; i++)
    {
        if (c[i]->BestCost < bestCost) // cost is negative
        {
            bestTrial = i;
            bestCost = c[i]->BestCost;
        }
    }
    
    /*====== print best clustering found =======*/
    printf("Best clustering, with cost %g is:\n", bestCost);
    c[bestTrial]->Print();
    
    /*====== print running time =========*/
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    fprintf(stderr,"Total time: %f s\n",elapsed_secs);
    
    for (int i=0; i<numTrials; i++)
    {
        delete c[i];
    }
    return 0;
}
