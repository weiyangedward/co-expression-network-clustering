#ifndef __CLUSTER__
#define __CLUSTER__

#include "cs-grn.h"
#define maxLogSize 1000

struct LogEntry
{
    int spc;
    int node;
    int oldstate;
};

class Cluster
{
    vector<double> spe_cost;
public:
    // constructor
    Cluster(SpeciesNetwork **specnws, int numS, Orthology *orthology, int numC, double cC);
    
    // clustering
    void LearnGroundState();
    
    SpeciesNetwork **sns;
    Orthology *orth;
    int numSpecies;
    int numClusters;
    double couplingConstant;
    double maxtemp;
    
    int ** snsClusterAssn;	// assignment of each node in each species to a cluster
    int ** bestClusterAssn;
    int * snsNumNodes; // number of nodes in each species network
    int * totalNumEdges;
    int ** degree;
    vector< vector< unordered_map<int, int> > > spe_cluster_nodes; // 3D arr to store nodes per cluster per spe
    
    double CurrentCost;
    double CurrentCostClus;
    
    double BestCost;
    
    struct LogEntry UndoLog[maxLogSize];
    int UndoLogSize;
    
    double OrthCost(int **ClusterAssn);
    double DeltaCost(int **ClusterAssn);
    double DeltaCostNew(int **ClusterAssn);
    double ScoreOfSpe(int spe);
    
    double Cost(int **ClusterAssn);
    int **Copy(int **ClusterAssn);
    void TotalNumEdges();
    void Delete(int **ClusterAssn);
    void CopyOver(int **ClusterAssn, int **oldClusterAssn);
    void Preset(char *filename);
    int Perturb(int **ClusterAssn, int numchanges);
    void UndoPerturb(int **ClusterAssnm, int new_state);
    void Print();
    void PrintBrief();
    void SetMaxTemp(double mt);
    
    int SA_counter;
    
};


#endif
