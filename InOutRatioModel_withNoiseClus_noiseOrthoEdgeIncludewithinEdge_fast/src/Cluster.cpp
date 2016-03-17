/*
 TODO:
 
 1) handle signal ctrl+c to output the best clustering resuts so far
 */

#include "Cluster.h"
#include "time.h"


/*==========
 constructor
 
 1) init 'snsNumNodes' to store total number of nodes
 2) init 'snsClusterAssn' to store current clustering with random assigned clusters
 3) init 'bestClusterAssn' to store best clustering
 4) init 'spe_cluster_nodes' to store nodes assigned to a cluster of a species
 
 */
Cluster::Cluster(SpeciesNetwork **specnws, int numS, Orthology *orthology, int numC, double cC, bool has_noiseCluster)
{
    //    srand((int)time(NULL)); // random seed for local time
    srand(1); // fix random seed
    
    sns = specnws;
    numSpecies = numS;
    orth = orthology;
    numClusters = numC;
    couplingConstant = cC;
    maxtemp = MAXTEMP;
    has_noise_cluster = has_noiseCluster;
    
    SA_counter = 0;
    UndoLogSize = 0;
    
    /* arr to store total number of nodes in a species */
    snsNumNodes = new int[numSpecies];
    for (int spe=0; spe<numSpecies; spe++)
    {
        snsNumNodes[spe] = sns[spe]->numNodes;
    }
    
    // init 3D arr to store nodes per cluster per species
    for (int spe = 0; spe<numSpecies; spe++)
    {
        vector< unordered_map<int, int> > tmp_cluster_node;
        for (int clus = 0; clus< numClusters; clus++)
        {
            unordered_map<int, int> tmp_node;
            tmp_cluster_node.push_back(tmp_node);
        }
        spe_cluster_nodes.push_back(tmp_cluster_node);
    }
    
    // 2D arr to store cluster label for each edge (i, j)
    snsClusterAssn = new int *[numSpecies];
    // 2D arr to store the best cluster label for each edge (i, j)
    bestClusterAssn = new int *[numSpecies];
    // initialize snsClusterAssn
    for (int spe=0; spe<numSpecies; spe++)
    {
        snsClusterAssn[spe] = new int[snsNumNodes[spe]];
        bestClusterAssn[spe] = new int[snsNumNodes[spe]];
        for (int node=0; node<snsNumNodes[spe]; node++)
        {
            // randomly assign clusters to each node
            int rand_index = 0;
            // noise cluster included, cluster0 will be assigned
            if (has_noise_cluster)
            {
                rand_index = rand() % numClusters;
            }
            // no noise cluster, then no cluster0, only cluster1 to clusterk will be assigned
            else
            {
               rand_index = (rand() % (numClusters-1))+1;
            }
            snsClusterAssn[spe][node] = rand_index; // [spe][node] = cluster
            spe_cluster_nodes[spe][rand_index][node] = 1; // [spe][cluster].emplace(node,1)
        }
    }
    
    // score the init cluster
    CurrentCostClus = Cost(snsClusterAssn);
    if (couplingConstant > 0)
    {
        CurrentCost = CurrentCostClus + couplingConstant * OrthCost(snsClusterAssn);
    }
    else
    {
        CurrentCost = CurrentCostClus;
    }
    
    // init BestCost to curr cost
    BestCost = CurrentCost;
}
/*==========
 Destructor
 
 1) de-allocate memory to prevent memory leak
 2) other memory de-allocation need to be done for structs in cs-grn and ReadInputs,
    work on this once the model is finalized
 */
Cluster::~Cluster()
{
    fprintf(stderr, "calling destructor ...\n");
    for (int spe=0; spe<numSpecies; spe++)
    {
        delete snsClusterAssn[spe];
        delete bestClusterAssn[spe];
    }
    delete snsNumNodes;
    delete snsClusterAssn;
    delete bestClusterAssn;
}

/*======================================
print final clustering result to output
 */
void Cluster::Print()
{
    for (int i=0; i<numSpecies; i++)
    {
        for (std::unordered_map<string, int>::iterator iter = sns[i]->nodeName2Id.begin(); iter != sns[i]->nodeName2Id.end(); ++iter)
        {
            string name = iter->first;
            int id = sns[i]->nodeName2Id[name];
            printf("Species %d\tGene %s\tCluster %d\n",i,name.c_str(),bestClusterAssn[i][id]);
        }
    }
}

void Cluster::PrintBrief()
{
    for (int i=0; i<numSpecies; i++) {
        printf("%d",snsClusterAssn[i][0]);
        for (int j=1; j<snsNumNodes[i]; j++) {
            printf("\t%d",snsClusterAssn[i][j]);
        }
        printf("\n");
    }
}

void Cluster::SetMaxTemp(double mt)
{
    maxtemp = mt;
}

void Cluster::Preset(char *filename)
{
    FILE *F = fopen(filename, "r");
    if (F==NULL) return;
    
    char line[linelen];
    fgets(line, linelen, F); // skip the single header line
    while (fgets(line, linelen, F)) {
        char spc[linelen];
        char id[linelen];
        char clusterid[linelen];
        char d1[linelen], d2[linelen], d3[linelen];
        // format of line is "Species <spc> Gene <gene> Cluster <cluster>"
        sscanf(line, "%s %s %s %s %s %s",(char *)&d1, (char *)&spc, (char *)&d2, (char *)&id, (char *)&d3, (char *)&clusterid);
        
        // convert spc id to numerics
        if (orth->spcName2Id.find(spc) == orth->spcName2Id.end()) {
            printf("Error: seed-clustering file %s has species name that wasnt seen in network files\n", filename);
            exit(1);
        }
        int spcid = orth->spcName2Id[spc];
        
        // convert gene id to numerics
        if (sns[spcid]->nodeName2Id.find(id) == sns[spcid]->nodeName2Id.end()) {
            continue; // this may happen because clusters file includes information about genes that are not there in
            // ... the network file (gene node with degree 0). It doesnt matter where (which cluster) we place such genes.
        }
        int geneid = sns[spcid]->nodeName2Id[id];
        
        // record information
        snsClusterAssn[spcid][geneid] = atoi(clusterid);
        fprintf(stderr,"Initial clustering has species %d gene %s assigned to cluster %d\n", spcid, id, snsClusterAssn[spcid][geneid]);
    }
    fclose(F);
    CurrentCostClus = Cost(snsClusterAssn);
    CurrentCost = CurrentCostClus + couplingConstant * OrthCost(snsClusterAssn);
}


/*============================================
 simulate annealing to find the best clustering
 
 1) for tempreture from max to min (max/1000),
    with step size = 0.9
        for 50*size_of_network
            compute delta_cost = new_cost - old_cost
            if delta_cost >= 0, accept new move with probability = exp(-delta_cost/temp)
            if delta_cost < 0, accept new move
    break for loop if reject new move for 500 times
 2) update best clustering when current clustering is better with a margine of 1
 */
void Cluster::LearnGroundState()
{
    fprintf(stderr,"Called LearnGroundState\n");
    double tempmax = maxtemp; // set start tempreture
    double tempmin = tempmax/1000; // set end tempreture
    
    // count total number of nodes in all species
    int totalNumNodes = 0;
    for (int i=0; i<numSpecies; i++)
        totalNumNodes += snsNumNodes[i];
    
    // counter for number of iterations that no better score is seen
    int nochangeiter = 0;
    
    /* simulated annealing (SA), each step = 0.9, 
     loop until tempreture is 1000 times smaller */
    for(double temp=tempmax; temp >= tempmin; temp *= 0.9)
    {
        // stop SA if > 500 iteration without change of cost
        if (nochangeiter > 500) break;
        
        // a "sweep" of the network does about 50 x as many changes as there are nodes overall
        for (int i=0; i < 20*totalNumNodes; i++)
        {
            // counter for SA iterations, used to determine when to print log
            SA_counter++;
            /* stop for loop if more than 500 iteration without change of cost, 
             here the counter 'nochangeiter' is the same as outer loop,
             which means once the inner loop break, the outer will also break */
            if (nochangeiter > 500) break;
            
            // record the old cost
            double OldCost = CurrentCost;
            double OldCostClus = CurrentCostClus;
            
            /* propose new assignment, by randomly perturbing the current assignment, 
             only one perturb occurs */
            int new_state = Perturb(snsClusterAssn, 1);
            
            // compute new cost of clusterng term
            double deltaCostCluster = DeltaCostNew(snsClusterAssn);
            double NewCostClus = OldCostClus + deltaCostCluster;
            
            double NewCost = 0.0;
            /* compute new cost adding ortho term
             only compute OrthCost when needed */
            if (couplingConstant > 0)
            {
                NewCost = NewCostClus + couplingConstant * OrthCost(snsClusterAssn);
            }
            else
            {
                NewCost = NewCostClus;
            }
            
            double deltaCost = NewCost - OldCost; // compute delta cost that includes orth term
            
//            if (SA_counter % 10000 == 0)
//                fprintf(stderr, "delta cost %g, old cost %g\n", deltaCost, OldCost);

            /* bad move, newCost >= oldCost */
            if (deltaCost >= 0)
            {
                /* reject bad move
                 when a random double (0~1) >= e^(-delta / t)
                 if new cost is more, i.e., worse,
                 delta is more positive, thus exp(-delta/T) is closer to 0,
                 thus reject probability is larger */
                if (double(rand())/RAND_MAX >= exp(-deltaCost/temp))
                {
                    // prev log
                    if (SA_counter % 10000 == 0)
                        fprintf(stderr, "REJECTED MOVE\t%g to %g at t %g, prob %g\n", OldCost, NewCost, temp, 1-exp(-deltaCost/temp));
                    
                    // undo pertubation
                    UndoPerturb(snsClusterAssn, new_state);
                    // increment counter no change of cost
                    nochangeiter++;
                }
                 /* accept bad move
                  with probability = exp(-deltaCost/temp) */
                else
                {
                    CurrentCost = NewCost;
                    CurrentCostClus = NewCostClus;
                    
                    if (SA_counter % 10000 == 0)
                        fprintf(stderr, "BAD MOVE\t%g to %g at t %g, prob %g\n", OldCost, NewCost, temp, exp(-deltaCost/temp)); // prev log
                    
                    // reset counter for no change of cost
                    nochangeiter = 0;
                    // reset counter if accept move
                    UndoLogSize = 0;
                }
            }
            /* good moves, accept it. 
             when deltaCost < 0, if newCost < oldCost */
            else
            {
                CurrentCost = NewCost;
                CurrentCostClus = NewCostClus;
                
                if (SA_counter % 10000 == 0)
                    fprintf(stderr, "GOOD MOVE\t%g to %g at t %g\n", OldCost, NewCost, temp); // prev log, print out good move
                
                /* any improvement less than this is not counted as an improvement, 
                 here set to 0, all improvements are counted */
                if (deltaCost < 0)
                    nochangeiter = 0;
                
                // reset log size
                UndoLogSize = 0;
                
                /* update BestCost
                 curr cost < BestCost - 1, here 1 is a margine of improvement */
                if (CurrentCost < BestCost -1)
                {
                    // update BestCost
                    BestCost = CurrentCost;
                    // update best clustering
                    for (int spe=0; spe<numSpecies; spe++)
                    {
                        for (int node=0; node<snsNumNodes[spe];node++)
                            bestClusterAssn[spe][node] = snsClusterAssn[spe][node];
                    }
                }
                
            }
            // print size of noise cluster
            if (SA_counter % 10000 == 0)
            {
                for (int spe = 0; spe<numSpecies; spe++)
                {
                    fprintf(stderr, "spe %d, noise size %lu\n", spe, spe_cluster_nodes[spe][0].size());
                }
                fprintf(stderr, "C = %g\n", CurrentCost); // prev log, print out curr cost
            }
        }
    }
}

/*============================================
 delta cost after perturb of a node in species
 
 1) if either the prev or new cluster of this node is noise cluster (0)
 then, compute the old and new score (delete this node from old cluster 
 and add to new cluster) of this species, define delta_score = new - old,
 and return -delta_score
 
 2) if both of prev and new cluster of this node are not noise cluster,
 then do the same thing as above to simplify computation steps. (it is possible
 to only score the new and old cluster, but that needs to take care of 
 adaptive nosie cost term)  
 */
double Cluster::DeltaCostNew(int **ClusterAssn)
{
    double delta_score = 0.0; // delta_cost = -delta_score
    
    int perturb_node = UndoLog[UndoLogSize-1].node;
    int spe = UndoLog[UndoLogSize-1].spc;
    int old_cluster = UndoLog[UndoLogSize-1].oldstate;
    int new_cluster = ClusterAssn[spe][perturb_node];
    
    double old_score = ScoreOfSpe(spe);
    
    spe_cluster_nodes[spe][old_cluster].erase(perturb_node);
    spe_cluster_nodes[spe][new_cluster][perturb_node] = 1;
    double new_score = ScoreOfSpe(spe);
    
    delta_score = new_score - old_score;
    return -delta_score;
}

/*===========================
 score clustering of a species
 1) used as a helper function for 'DeltaCostNew' to score
    new and old clusterings
 */
double Cluster::ScoreOfSpe(int spe)
{
    double spe_score = 0.0;
    
    unordered_map<int, int>::iterator node;
    list<int>::iterator adj_gene;
    unsigned long noise_cluster_size = spe_cluster_nodes[spe][0].size();
    // ignore noise cluster here
    for (int clus = 1; clus < numClusters; clus++)
    {
        int in_edges = 0;
        int out_edges = 0;
        unsigned long cluster_size = spe_cluster_nodes[spe][clus].size();
        for (node = spe_cluster_nodes[spe][clus].begin(); node!= spe_cluster_nodes[spe][clus].end(); node++)
        {
            for (adj_gene = sns[spe]->adjacencyList[(node->first)].begin(); adj_gene != sns[spe]->adjacencyList[node->first].end(); adj_gene++)
            {
                // adj_gene in the same cluster as cur gene
                if (spe_cluster_nodes[spe][clus].find(*adj_gene) !=spe_cluster_nodes[spe][clus].end())
                {
                    in_edges++;
                }
                else
                {
                    // adj_gene not in noise cluster
                    if (spe_cluster_nodes[spe][0].find(*adj_gene) == spe_cluster_nodes[spe][0].end())
                    {
                        out_edges++;
                    }
                }
            }
        }
        in_edges /=2;
        // cluster with size <= 1 will have possible_in_edge =0
        if (cluster_size >=2 && (in_edges > 0 || out_edges > 0))
        {
            unsigned long possible_in_edge = (cluster_size) * (cluster_size - 1) / 2;
            unsigned long possible_out_edge = ( snsNumNodes[spe] - (noise_cluster_size) - (cluster_size) ) * (cluster_size);
            
            double in_density = double(in_edges) /possible_in_edge;
            double out_density = double(out_edges) / possible_out_edge;
            
            double inOutRatio = (in_density / (in_density + out_density)) * double(cluster_size);
            
            spe_score += inOutRatio;
//            fprintf(stderr, "spe %d, cluster size %d, indensity %g, outdensity %g\n", spe, cluster_size, in_density, out_density);
        }
    }
    double noise_cost = (double)noise_cluster_size * (spe_score / snsNumNodes[spe]);
    spe_score +=noise_cost;
//    if (SA_counter % 10000 == 0)
//        fprintf(stderr, "spe %d, noise cluster size %d, spe_score %g\n", spe, noise_cluster_size, spe_score);
    return spe_score;
}


/*==============================
 score a clustering of a network
 1) loop through each node in a speices and assign them to corresponding clusters,
    count in-cluster edges, out-cluster edges, in-cluster nodes
 2) loop through each cluster of a species, compute size(cluster) * In-density/(In-density + Out-density)
 3) for each species, sum up InOutRatio score and compute adaptive noise score
 4) return -score as cost
 */
double Cluster::Cost(int **ClusterAssn)
{
    // init total score = 0
    double score = 0.0;

    for (int spe = 0; spe < numSpecies; spe++)
    {
        std::unordered_map<int, int> total_cluster_outEdges; // #out-cluster edges
        std::unordered_map<int, int> total_cluster_edge; // #in-cluster edges
        std::map<int, int> total_cluster_nodes; // #nodes in cluster, ordered
        int noise_clus_size = 0; // #nodes in noise cluster
        
        /* cost in normal clusters */
        // loop through node_i
        for (int i=0; i<snsNumNodes[spe]; i++)
        {
            // skip noise cluster
            if (ClusterAssn[spe][i] != 0)
            {
                total_cluster_edge[ClusterAssn[spe][i]] += 0; // initialize in-cluster total edge = 0
                // increment #nodes in clus
                total_cluster_nodes[ClusterAssn[spe][i]] ++;
                /* loop through node_j in a spe, start from j=0 
                 so that the out_edge can be double counted, which is needed */
                for (int j=0; j<snsNumNodes[spe]; j++)
                {
                    // skip noise cluster
                    if (ClusterAssn[spe][j] != 0)
                    {
                        // i and j in the same cluster
                        if (ClusterAssn[spe][i] == ClusterAssn[spe][j])
                        {
                            // if edge between i and j
                            if (sns[spe]->IsEdge(i,j))
                            {
                                // increment #in-cluster edges
                                total_cluster_edge[ClusterAssn[spe][i]] ++;
                            }
                        }
                        // i and j are in diff cluster
                        else
                        {
                            if (sns[spe]->IsEdge(i,j))
                            {
                                // increment #out-cluster edges
                                total_cluster_outEdges[ClusterAssn[spe][i]] ++;
                            }
                        }
                    }
                }
            }
            // if cur cluster is noise cluster
            else
            {
                // increment #nodes in noise cluster
                noise_clus_size ++;
            }
        }
        // init cost of spc
        double spe_score = 0.0;
        std::map<int, int>::iterator it;
        // loop through all clusters of a species, without noise cluster
        for (it=total_cluster_nodes.begin(); it != total_cluster_nodes.end(); it++)
        {
            int cluster = it->first, cluster_size =it->second;
            // #in-cluster edges > 0, do not need to worry about 'noise cluster' since it does not take into account above
            if (total_cluster_edge[cluster] > 0)
            {
                // #possible in-cluster edges
                int possible_in_edge = cluster_size * (cluster_size-1) / 2;
                // #possible out-cluster edges
                int possible_out_edge = ( snsNumNodes[spe] - noise_clus_size - (cluster_size) ) * (cluster_size);
                // in-density
                double total_cluster_inDensity = double(total_cluster_edge[cluster]/2) / double(possible_in_edge);
                // out-density
                double total_cluster_outDensity = double(total_cluster_outEdges[cluster]) / double(possible_out_edge);
                // cost of a spc = in-density * size / (in-density + out-density)
                spe_score += (total_cluster_inDensity / (total_cluster_inDensity + total_cluster_outDensity)) * double(cluster_size);
                
                fprintf(stderr, "cluster: %d, size: %d, in_edge: %d, out_edge: %d, possible_in_edge: %d\n", cluster, cluster_size, total_cluster_edge[cluster]/2, total_cluster_outEdges[cluster], possible_in_edge);
            }
        }
        // add spe_score to total score
        score += spe_score;
        //        cost += (double)noise_clus_size * 0.7;
        // cost of noise cluster = #nodes in noise cluster * (spc_cost / size) -> adaptive weight for noise cluster
        score += (double)noise_clus_size * (spe_score / snsNumNodes[spe]);
        
//        if (SA_counter % 1000 == 0)
            fprintf(stderr, "spe %d, noise cluster size %d, spe_score/#nodes %g\n", spe, noise_clus_size, (spe_score / snsNumNodes[spe]));
    }
    // return -score as cost
    return - score;
}

/*===============================
 OrthCost of all pairs of species
 1) for spe1 from 1 to total spe
        for spe2 from spe1+1 to total spe
            for each ortholog edge
                if both end nodes are in the same cluster, orth_score++
                if only one end node is in noise cluster, noise_score++
 2) return -score as cost
 */
double Cluster::OrthCost(int **ClusterAssn)
{
    double orthterms = 0.0; // orth cost without noise clusters
    double total_noise_nodes = 0.0;
    double orth_noise_term = 0.0; // orth cost of noise clusters
    double total_gene_among_spc = 0.0; // total number of genes among three spes
    double total_ortho_noise = 0.0;
    // loop through all spc
    for (int spc1=0; spc1 < numSpecies; spc1++)
    {
        // sum up total genes
        total_gene_among_spc += snsNumNodes[spc1];
        // sum up noise nodes
        total_noise_nodes += spe_cluster_nodes[spc1][0].size();
        // loop through spc + 1, avoid double-counting ortholog edges
        for (int spc2=spc1+1; spc2 < numSpecies; spc2++)
        {
            unordered_map<NodePair *, int, NodePairHasher, eqNodePair>::iterator iter; // hashmap iterator
            // loop through all ortholog gene pairs between two spc
            for ( iter = orth->orthtable[spc1][spc2].begin(); iter != orth->orthtable[spc1][spc2].end(); ++iter)
            {
                NodePair *np = (NodePair *)(iter->first);
                int n1 = np->id1; // node1
                int n2 = np->id2; // node2
                // count edges where both end nodes are in the same cluster, not counting nodes in noise cluster
                if (ClusterAssn[spc1][n1] == ClusterAssn[spc2][n2] && ClusterAssn[spc1][n1] != 0 && ClusterAssn[spc2][n2] != 0)
                {
                    orthterms += (orth->weighted_orth[spc1][spc2][n1] + orth->weighted_orth[spc1][spc2][n2]) / 2.0;
                }
                /* For noise nodes, only count ortholog edges that cross noise cluster. When both end of an ortho edge are in noise cluster, chances are these pair of nodes have co-expression edges, and so should be put in normal clusters.
                 This rewards:
                 1) all 3 ortho edge within a normal cluster
                 2) one end of ortho edge in noise cluster, and the other in normal cluster, max reward happens when end node of the other two species are in the same normal cluster: 2(cross noise) + 1(within normal) = 3
                 3) when 1 and 2 become a tie, then In/Out-density ratio will decide if a noise should be moved out of noise cluster
                 */
//                if (ClusterAssn[spc1][n1] != ClusterAssn[spc2][n2] && (ClusterAssn[spc1][n1] == 0 || ClusterAssn[spc2][n2] == 0))
                if ((ClusterAssn[spc1][n1] == 0 || ClusterAssn[spc2][n2] == 0))
                {
                    total_ortho_noise += (orth->weighted_orth[spc1][spc2][n1] + orth->weighted_orth[spc1][spc2][n2]) / 2.0;
                }
            }
        }
    }
    // adaptive weighting noise term by orthterm / #total genes
//    orth_noise_term = total_noise_nodes * orthterms / total_gene_among_spc;
//    orth_noise_term = total_ortho_noise * orthterms / total_gene_among_spc;
    orth_noise_term = total_ortho_noise;
    
    orthterms += orth_noise_term;
    
//    if (SA_counter % 10000 == 0)
//        fprintf(stderr, "total noise node %g, total ortho noise %g\n", total_noise_nodes, total_ortho_noise);
    //    fprintf(stderr, "cc = %g\tcccost = %g\n",couplingConstant, couplingConstant*orthterms); // prev log, print orth term
    return -orthterms;
}


int **Cluster::Copy(int **ClusterAssn)
{
    int **ca = new int *[numSpecies];
    for (int i=0; i<numSpecies; i++)
    {
        ca[i] = new int[snsNumNodes[i]];
        for (int j=0; j<snsNumNodes[i]; j++) ca[i][j] = ClusterAssn[i][j];
    }
    return ca;
}

void Cluster::Delete(int **ClusterAssn)
{
    for (int i=0; i<numSpecies; i++)
    {
        delete [] ClusterAssn[i];
    }
    delete [] ClusterAssn;
}

void Cluster::CopyOver(int **ClusterAssn, int **oldClusterAssn)
{
    for (int i=0; i<numSpecies; i++)
    {
        for (int j=0; j<snsNumNodes[i]; j++)
            ClusterAssn[i][j] = oldClusterAssn[i][j];
    }
}

/*================================================================
 UndoPerturb
 1) recover perturbed node to its old cluster
 2) for spe_cluster_nodes, delete perturbed node from new cluster,
    add perturbed node back to old cluster
 */
void Cluster::UndoPerturb(int **ClusterAssn, int new_state)
{
    for (int i = UndoLogSize; i>0; i--)
    {
        ClusterAssn[UndoLog[i-1].spc][UndoLog[i-1].node] = UndoLog[i-1].oldstate;
        spe_cluster_nodes[UndoLog[i-1].spc][new_state].erase(UndoLog[i-1].node);
        spe_cluster_nodes[UndoLog[i-1].spc][UndoLog[i-1].oldstate][UndoLog[i-1].node] = 1;
    }
    
    UndoLogSize = 0;
}

/*===========================
 Perturb
 1) randomly select a species
 2) randomly select a node
 3) randomly select a new cluster
 4) assign this to the new cluster
 */
int Cluster::Perturb(int **ClusterAssn, int numchanges = 1)
{
    int newstate = 0;
    for (int i=0; i<numchanges; i++)
    {
        // choose a random species
        int spc = rand() % numSpecies;
        // choose a random node
        int node = rand() % (snsNumNodes[spc]); // node is a uniq_ID instead of gene_ID
        int oldstate = ClusterAssn[spc][node];
        
        // choose a random cluster
        while (1)
        {
            // maxLogSize = 1000
            if (UndoLogSize > maxLogSize) {
                printf("Error: Undo Log Size Limit reached\n");
                exit(1);
            }
            
            // choose a randome new cluster
            int rand_index = 0;
            if (has_noise_cluster)
            {
                rand_index = rand() % numClusters; // if noise cluster included
            }
            else
            {
                rand_index = (rand() % (numClusters-1))+1; // no cluster0
            }
            newstate = rand_index;
            
            /* if new cluster is the same as old,
            jump to the end of loop, not increase UndoLogSize */
            if (newstate == oldstate) continue;
            
            // assign node to new cluster
            ClusterAssn[spc][node] = newstate;
            
            if (SA_counter % 100000 == 0)
                fprintf(stderr,"Per\tspe%d\tnode%s\told%d\tnew%d\n",spc,sns[spc]->nodeId2Name[node].c_str(),oldstate,newstate); // prev log
            
            // record the change for later undo
            UndoLog[UndoLogSize].spc = spc;
            UndoLog[UndoLogSize].node = node;
            UndoLog[UndoLogSize].oldstate = oldstate;
            UndoLogSize++;
            // stop the while loop
            break;
        }
    }
    return newstate;
}
