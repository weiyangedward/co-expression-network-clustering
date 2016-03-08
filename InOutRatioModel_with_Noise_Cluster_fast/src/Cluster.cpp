#include "Cluster.h"
#include "time.h"

/*===== constructor ======*/
Cluster::Cluster(SpeciesNetwork **specnws, int numS, Orthology *orthology, int numC, double cC)
{
    SA_counter = 0;
    sns = specnws;
    numSpecies = numS;
    orth = orthology;
    numClusters = numC;
    couplingConstant = cC;
    maxtemp = MAXTEMP;
    
    UndoLogSize = 0;
    snsNumNodes = new int[numSpecies];
    spe_cost.resize(numSpecies);
    for (int i=0; i<numSpecies; i++)
    {
        snsNumNodes[i] = sns[i]->numNodes;
        spe_cost[i] = 0.0;
    }
    
    // init 3D arr to store nodes per cluster per species
    for (int spe = 0; spe<numSpecies; spe++)
    {
        vector< unordered_set<int> > tmp_cluster_node;
        for (int clus = 0; clus< numClusters; clus++)
        {
            unordered_set<int> tmp_node;
            tmp_cluster_node.push_back(tmp_node);
        }
        spe_cluster_nodes.push_back(tmp_cluster_node);
    }
    
    // arr for total number of edges for each species
    totalNumEdges = new int[numSpecies];
    
    // 2D arr for degree of each node in each species
    degree = new int *[numSpecies];
    
    for (int i=0; i<numSpecies; i++)
    {
        degree[i] = new int[snsNumNodes[i]];
    }
    
    // count total number of edges and degree (do not need this )
//    TotalNumEdges();
    
    srand((int)time(NULL)); // random seed for local time
//         srand(1); // fix random seed
    
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
            snsClusterAssn[spe][node] = rand() % numClusters; // [spe][node] = cluster
            spe_cluster_nodes[spe][snsClusterAssn[spe][node]].emplace(node); // [spe][cluster].push_back(node)
        }
    }
    
    CurrentCostClus = Cost(snsClusterAssn); // compute the cost of clustering term given initial cluster assignment
    CurrentCost = CurrentCostClus + couplingConstant * OrthCost(snsClusterAssn); // compute cost = clustering term + ortho terms
    
    // initialize BestCost to curr cost
    BestCost = CurrentCost;
}

void Cluster::Print()
{
    for (int i=0; i<numSpecies; i++)
    {
        //        for (hash_map<const char*, int, hash<const char*>, eqstr>::iterator iter = sns[i]->nodeName2Id.begin(); iter != sns[i]->nodeName2Id.end(); ++iter) {
        for (std::unordered_map<string, int>::iterator iter = sns[i]->nodeName2Id.begin(); iter != sns[i]->nodeName2Id.end(); ++iter)
        {
            //            const char *name = (const char *)(iter->first);
            string name = iter->first;
            int id = sns[i]->nodeName2Id[name];
            // printf("Species %d\tGene %s\tCluster %d\n",i,name,snsClusterAssn[i][id]);
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

/*============
 TotalNumEdges (need to consider nodes in noise cluster!!!!!!!!!!)
 count total number of edges for each species
 */
void Cluster::TotalNumEdges()
{
    for (int spc = 0; spc < numSpecies; spc++)
    {
        // compute degree of each node
        // int *degree = new int[snsNumNodes[spc]];
        // initialize degree value
        for (int i=0; i<snsNumNodes[spc]; i++) degree[spc][i] = 0;
        
        for (int i=0; i<snsNumNodes[spc]; i++)
        {
            for (int j=i+1; j<snsNumNodes[spc]; j++)
            {
                if (sns[spc]->IsEdge(i,j)) { // if there is an edge between i and j
                    degree[spc][i]++;
                    degree[spc][j]++;
                }
            }
        }
        totalNumEdges[spc] = 0;
        for (int i=0; i<snsNumNodes[spc]; i++) {
            totalNumEdges[spc] += degree[spc][i];
        }
        // total number of edges, degree of both nodes added one each time
        totalNumEdges[spc] /= 2;
    }
}

/*===============
 LearnGroundState
 simulate annealing 
 */
void Cluster::LearnGroundState()
{
    
    fprintf(stderr,"Called LearnGroundState\n");
    double tempmax = maxtemp; // set start tempreture
    double tempmin = tempmax/1000; // set end tempreture
    
    /*====== count total node No. among species ========*/
    int totalNumNodes = 0;
    for (int i=0; i<numSpecies; i++)
        totalNumNodes += snsNumNodes[i];
    
    int nochangeiter = 0; // how many iterations have we seen no change in score
    
    /*======= simulated annealing (SA) =========*/
    for(double temp=tempmax; temp >= tempmin; temp *= 0.9)
    {
        if (nochangeiter > 500) break; // stop SA if > 500 iteration without change of cost
        
        // a "sweep" of the network does about 50 x as many changes as there are nodes overall
        for (int i=0; i < 50*totalNumNodes; i++)
        {
            SA_counter++;
            if (nochangeiter > 500) break; // stop for loop if more than 500 iteration without change of cost, here the counter 'nochangeiter' is the same as outter loop, which means once the inner loop break, the outter will also break
            
            double OldCost = CurrentCost; // set old cost = curr cost
            double OldCostClus = CurrentCostClus; // set old cost of clustering term to curr cost of clustering term
            
            // propose new assignment, by randomly perturbing the current assignment, Note that large number of perturb gives worse resutls
            int new_state = Perturb(snsClusterAssn, 1);
            
            // compute new cost of clusterng term
            double deltaCostCluster = DeltaCostNew(snsClusterAssn);
//            fprintf(stderr, "delta cost %g, old cost %g\n", deltaCostCluster, OldCost);
            double NewCostClus = OldCostClus + deltaCostCluster; // compute new clustering cost term
            double NewCost = NewCostClus + couplingConstant * OrthCost(snsClusterAssn); // compute new cost = clustering term + orth term (has noise cluster node been addressed? Yes it has)
            
            double deltaCost = NewCost - OldCost; // compute delta cost that includes orth term
            if (SA_counter % 10000 == 0)
                fprintf(stderr, "delta cost %g, old cost %g\n", deltaCost, OldCost);
            if (new_state == 0 and deltaCost > 0) // if node is moved to 'noise cluster' and has a worse cost
            {
                deltaCost -= 0.0; // make its cost better to 'help' a node to move into 'noise cluster', set to 0 by now
//                NewCost -= 0.0; // change the newCost so that there won't be a worse cost being accepted
            }
            
            // bad move
            if (deltaCost >= 0)  // if newCost >= oldCost
            {
                // reject bad move
                if (double(rand())/RAND_MAX >= exp(-deltaCost/temp)) // if a random double (0~1) >= e^(-delta / t), see explaination below
                {
                    /*==============================
                     if new cost is more, i.e., worse,
                     delta is more positive, thus exp(-delta/T) is closer to 0,
                     thus reject probability is larger
                     */
                    if (SA_counter % 10000 == 0)
                        fprintf(stderr, "REJECTED MOVE\t%g to %g at t %g, prob %g\n", OldCost, NewCost, temp, 1-exp(-deltaCost/temp)); // prev log
                    
                    UndoPerturb(snsClusterAssn, new_state); // undo pertubation
                    nochangeiter++; // increment counter no change of cost
                }
                else // accept bad move
                {
                    CurrentCost = NewCost;
                    CurrentCostClus = NewCostClus;
                    
                    if (SA_counter % 10000 == 0)
                        fprintf(stderr, "BAD MOVE\t%g to %g at t %g, prob %g\n", OldCost, NewCost, temp, exp(-deltaCost/temp)); // prev log
                    
                    nochangeiter = 0; // reset counter for no change of cost
                    UndoLogSize = 0; // reset counter if accept move
                }
            }
            else // deltaCost < 0, accept good moves, if newCost < oldCost
            {
                CurrentCost = NewCost;
                CurrentCostClus = NewCostClus;
                
                if (SA_counter % 10000 == 0)
                    fprintf(stderr, "GOOD MOVE\t%g to %g at t %g\n", OldCost, NewCost, temp); // prev log, print out good move
                
                if (deltaCost < 0) // any improvement less than this is not counted as an improvement
                    nochangeiter = 0; // reset counter if accept move
                
                UndoLogSize = 0; // reset log size
                
                // update BestCost
                if (CurrentCost < BestCost -1) // if curr cost < BestCost - 1, here 1 is a margine of improvement
                {
                    BestCost = CurrentCost; // update BestCost
                    for (int i=0; i<numSpecies; i++) // update best clustering
                    {
                        for (int j=0; j<snsNumNodes[i];j++)
                            bestClusterAssn[i][j] = snsClusterAssn[i][j];
                    }
                }
                
            }
            if (SA_counter % 10000 == 0)
                fprintf(stderr, "C = %g\n", CurrentCost); // prev log, print out curr cost
        }
    }
}

/* compute the delta cost after perturb of a node in species
 
 1) if either the prev or new cluster of this node is noise cluster (0)
 then, compute the old and new score (delete this node from old cluster 
 and add to new cluster) of this species, define delta_score = new - old,
 and return -delta_score
 2) if both of prev and new cluster of this node are not noise cluster,
 then do the same thing as above to simplify step so far. (it is possible
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
    spe_cluster_nodes[spe][new_cluster].emplace(perturb_node);
    double new_score = ScoreOfSpe(spe);
    
    delta_score = new_score - old_score;
    return -delta_score;
}


double Cluster::ScoreOfSpe(int spe)
{
    double spe_score = 0.0;
    
    unordered_set<int>::iterator node;
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
            for (adj_gene = sns[spe]->adjacencyList[*node].begin(); adj_gene != sns[spe]->adjacencyList[*node].end(); adj_gene++)
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
    if (SA_counter % 10000 == 0)
        fprintf(stderr, "spe %d, noise cluster size %d, spe_score %g\n", spe, noise_cluster_size, spe_score);
    return spe_score;
}


/*======================
 Cost of clustering term
 compute cost of clsutering term
 */
double Cluster::Cost(int **ClusterAssn)
{
    // init total cost = 0
    double cost = 0.0;
    // loop through all species
    for (int spc = 0; spc < numSpecies; spc++)
    {
        std::unordered_map<int, int> total_cluster_outEdges; // #out-cluster edges
        std::unordered_map<int, int> total_cluster_edge; // #in-cluster edges
        std::map<int, int> total_cluster_nodes; // #nodes in cluster
        int noise_clus_size = 0; // #nodes in noise cluster
        
        /* cost in normal clusters */
        // loop through each node in a spc
        for (int i=0; i<snsNumNodes[spc]; i++)
        {
            // skip noise cluster
            if (ClusterAssn[spc][i] != 0)
            {
                total_cluster_edge[ClusterAssn[spc][i]] += 0; // initialize in-cluster total edge = 0
                // increment #nodes in clus
                total_cluster_nodes[ClusterAssn[spc][i]] ++;
                // loop through each node in a spc, start from  j=0 so that the out_edge can be double counted, which is needed
                for (int j=0; j<snsNumNodes[spc]; j++)
                {
                    // skip noise cluster
                    if (ClusterAssn[spc][j] != 0)
                    {
                        // i and j in the same cluster
                        if (ClusterAssn[spc][i] == ClusterAssn[spc][j])
                        {
                            // if edge between i and j
                            if (sns[spc]->IsEdge(i,j))
                            {
                                // increment #in-cluster edges
                                total_cluster_edge[ClusterAssn[spc][i]] ++;
                            }
                        }
                        // i and j are in diff cluster
                        else
                        {
                            if (sns[spc]->IsEdge(i,j))
                            {
                                // increment #out-cluster edges
                                total_cluster_outEdges[ClusterAssn[spc][i]] ++;
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
        double spc_cost = 0.0;
        std::map<int, int>::iterator it;
        // loop through all clusters, without noise cluster
        for (it=total_cluster_nodes.begin(); it != total_cluster_nodes.end(); it++)
        {
            int cluster = it->first, cluster_size =it->second;
            // #in-cluster edges > 0, do not need to worry about 'noise cluster' since it does not take into account above
            if (total_cluster_edge[cluster] > 0)
            {
                // #possible in-cluster edges
                int possible_in_edge = cluster_size * (cluster_size-1) / 2;
                // #possible out-cluster edges
                int possible_out_edge = ( snsNumNodes[spc] - noise_clus_size - (cluster_size) ) * (cluster_size);
                // in-density
                double total_cluster_inDensity = double(total_cluster_edge[cluster]/2) / double(possible_in_edge);
                // out-density
                double total_cluster_outDensity = double(total_cluster_outEdges[cluster]) / double(possible_out_edge);
                // cost of a spc = in-density * size / (in-density + out-density)
                spc_cost += (total_cluster_inDensity / (total_cluster_inDensity + total_cluster_outDensity)) * double(cluster_size);
                
                fprintf(stderr, "cluster: %d, size: %d, in_edge: %d, out_edge: %d, possible_in_edge: %d\n", cluster, cluster_size, total_cluster_edge[cluster]/2, total_cluster_outEdges[cluster], possible_in_edge);
            }
        }
        // add spc_cost to total cost
        cost += spc_cost;
        spe_cost[spc] = -spc_cost;
        //        cost += (double)noise_clus_size * 0.7;
        // cost of noise cluster = #nodes in noise cluster * (spc_cost / size) -> adaptive weight for noise cluster
        cost += (double)noise_clus_size * (spc_cost / snsNumNodes[spc]);
        
//        if (SA_counter % 1000 == 0)
            fprintf(stderr, "spc %d, noise cluster size %d, cost/#nodes %g\n", spc, noise_clus_size, (spc_cost / snsNumNodes[spc]));
    }
    // return -cost since here cost is a score
    return - cost;
}

/*=======
 OrthCost (just for now set cc = 0)
 compute cost of orth term
 */
double Cluster::OrthCost(int **ClusterAssn)
{
    double orthterms = 0.0; // orth cost without noise clusters
    double total_noise_nodes = 0.0;
    double orth_noise_term = 0.0; // orth cost of noise clusters
    double total_gene_among_spc = 0.0; // total number of genes among three spes
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
                if (ClusterAssn[spc1][n1] == ClusterAssn[spc2][n2] and ClusterAssn[spc1][n1] != 0 and  ClusterAssn[spc2][n2] != 0)
                {
                    orthterms += (orth->weighted_orth[spc1][spc2][n1] + orth->weighted_orth[spc1][spc2][n2]) / 2.0;
                }
            }
        }
    }
    // adaptive weighting noise term by orthterm / #total genes
    orth_noise_term = total_noise_nodes * orthterms / total_gene_among_spc;
    //    orth_noise_term = total_noise_nodes;
    // ortholog term = orthlog term in network + #noise nodes
    orthterms += orth_noise_term;
    //    fprintf(stderr, "cc = %g\tcccost = %g\n",couplingConstant, couplingConstant*orthterms); // prev log, print orth term
    return -orthterms;
}


int **Cluster::Copy(int **ClusterAssn) {
    
    int **ca = new int *[numSpecies];
    for (int i=0; i<numSpecies; i++) {
        ca[i] = new int[snsNumNodes[i]];
        for (int j=0; j<snsNumNodes[i]; j++) ca[i][j] = ClusterAssn[i][j];
    }
    return ca;
}

void Cluster::Delete(int **ClusterAssn) {
    
    for (int i=0; i<numSpecies; i++) {
        delete [] ClusterAssn[i];
    }
    delete [] ClusterAssn;
}

void Cluster::CopyOver(int **ClusterAssn, int **oldClusterAssn) {
    
    for (int i=0; i<numSpecies; i++) {
        for (int j=0; j<snsNumNodes[i]; j++) ClusterAssn[i][j] = oldClusterAssn[i][j];
    }
}

/*==========
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
        spe_cluster_nodes[UndoLog[i-1].spc][UndoLog[i-1].oldstate].emplace(UndoLog[i-1].node);
    }
    
    UndoLogSize = 0;
}

/*======
 Perturb
 randomly assign a cluster to a random node
 */
int Cluster::Perturb(int **ClusterAssn, int numchanges = 1)
{
    int newstate = 0;
    for (int i=0; i<numchanges; i++)
    {
        // choose a random species
        int spc = rand() % numSpecies;
        // chose a random node
        int node = rand() % (snsNumNodes[spc]); // node is a uniq_ID instead of gene_ID
        int oldstate = ClusterAssn[spc][node];
        
        while (1)
        {
            newstate = rand() % numClusters; // choose a randome new cluster
            /*================*/
//            newstate = (rand() % (numClusters-1))+1;
            /*================*/
            if (newstate == oldstate) continue;
            if (SA_counter % 10000 == 0)
                fprintf(stderr,"Per\tspe%d\tnode%s\told%d\tnew%d\n",spc,sns[spc]->nodeId2Name[node].c_str(),oldstate,newstate); // prev log
            // assign node to new cluster
            ClusterAssn[spc][node] = newstate;
            
            if (UndoLogSize > maxLogSize) { // maxLogSize = 1000
                printf("Error: Undo Log Size Limit reached\n");
                exit(1);
            }
            
            // record the change for later undo
            UndoLog[UndoLogSize].spc = spc;
            UndoLog[UndoLogSize].node = node;
            UndoLog[UndoLogSize].oldstate = oldstate;
//            fprintf(stderr, "Pertub - species: %d, node: %d, old cluster: %d, new cluster: %d\n", spc, node, oldstate, newstate);
            UndoLogSize++;
            break; // stop the while loop
        }
    }
    return newstate;
}
