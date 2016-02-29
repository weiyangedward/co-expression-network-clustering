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
    num_node_noise_in_cluster = new int [numS]; // arr of numS to count number of nodes in 'noise cluster'
    snsNumNodes = new int[numSpecies];
    
    for (int i=0; i<numSpecies; i++)
    {
        snsNumNodes[i] = sns[i]->numNodes;
    }
    
    // arr for total number of edges for each species
    totalNumEdges = new int[numSpecies];
    
    // 2D arr for degree of each node in each species
    degree = new int *[numSpecies];
    
    for (int i=0; i<numSpecies; i++)
    {
        degree[i] = new int[snsNumNodes[i]];
        num_node_noise_in_cluster[i] = 0; // initialize number of nodes in noise cluster to be 0
    }
    
    // count total number of edges and degree (do not need this )
    //    TotalNumEdges();
    
    srand((int)time(NULL)); // random seed for local time
    // srand(1); // fix random seed
    
    // 2D arr to store cluster label for each edge (i, j)
    snsClusterAssn = new int *[numSpecies];
    
    // 2D arr to store the best cluster label for each edge (i, j)
    bestClusterAssn = new int *[numSpecies];
    
    // initialize snsClusterAssn
    for (int i=0; i<numSpecies; i++)
    {
        snsClusterAssn[i] = new int[snsNumNodes[i]];
        bestClusterAssn[i] = new int[snsNumNodes[i]];
        for (int j=0; j<snsNumNodes[i]; j++)
        {
            snsClusterAssn[i][j] = rand() % numClusters;
            /*================*/
            //            snsClusterAssn[i][j] = (rand() % (numClusters-1))+1;
            /*================*/
            if (snsClusterAssn[i][j] == 0)
                num_node_noise_in_cluster[i] ++;
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
    // m.unlock();
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
            pair<int, int> new_old_state = Perturb(snsClusterAssn, 1);
            
            // score new assignment
//            double deltaCostCluster = DeltaCost(snsClusterAssn); // compute delta cost of clusterng term between new cost and curr cost
//            double NewCostClus = OldCostClus + deltaCostCluster; // compute new clustering cost term
//            double NewCost = NewCostClus + couplingConstant * OrthCost(snsClusterAssn); // compute new cost = clustering term + orth term (has noise cluster node been addressed? Yes it has)
            
            double NewCostClus = Cost(snsClusterAssn);
            double NewCost = NewCostClus + couplingConstant * OrthCost(snsClusterAssn);
            
            double deltaCost = NewCost - OldCost; // compute delta cost that includes orth term
//            fprintf(stderr, "deltaCost in learnGroundState: %g\n", deltaCost);
            if (new_old_state.first == 0 and deltaCost > 0) // if node is moved to 'noise cluster' and has a worse cost
            {
                deltaCost -= 0; // make its cost better to 'help' a node to move into 'noise cluster', set to 0 by now
                //                NewCost -= 0.5; // change the newCost so that there won't be a worse cost being accepted
            }
            
            if (new_old_state.second == 0 and deltaCost > 0) // if node is moved to 'noise cluster' and has a worse cost
            {
                deltaCost -= 0; // make its cost better to 'help' a node to move into 'noise cluster', set to 0 by now
                //                NewCost -= 0.5; // change the newCost so that there won't be a worse cost being accepted
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
                    if (SA_counter % 1000 == 0)
                        fprintf(stderr, "REJECTED MOVE\t%g to %g at t %g, prob %g\n", OldCost, NewCost, temp, 1-exp(-deltaCost/temp)); // prev log
                    UndoPerturb(snsClusterAssn, new_old_state.first); // undo pertubation
                    nochangeiter++; // increment counter no change of cost
                }
                else // accept bad move
                {
                    CurrentCost = NewCost;
                    CurrentCostClus = NewCostClus;
                    if (SA_counter % 1000 == 0)
                        fprintf(stderr, "BAD MOVE\t%g to %g at t %g, prob %g\n", OldCost, NewCost, temp, exp(-deltaCost/temp)); // prev log
                    nochangeiter = 0; // reset counter for no change of cost
                    UndoLogSize = 0; // reset counter if accept move
                }
            }
            else // deltaCost < 0, accept good moves, if newCost < oldCost
            {
                CurrentCost = NewCost;
                CurrentCostClus = NewCostClus;
                
                if (SA_counter % 1000 == 0)
                    fprintf(stderr, "GOOD MOVE\t%g to %g at t %g\n", OldCost, NewCost, temp); // prev log, print out good move
                
                if (deltaCost < 0) // any improvement less than this is not counted as an improvement
                    nochangeiter = 0; // reset counter if a good move is accepted
                
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
            if (SA_counter % 1000 == 0)
                fprintf(stderr, "C = %g\n", CurrentCost); // prev log, print out curr cost
        }
    }
}
/*============================
 DeltaCost of clustering term
 compute delta cost of clustering term
 */
double Cluster::DeltaCost(int **ClusterAssn)
{
    
    double delta_cost = 0.0;
    
    int perturb_node = UndoLog[UndoLogSize-1].node;
    int spc = UndoLog[UndoLogSize-1].spc;
    int oldstate = UndoLog[UndoLogSize-1].oldstate;
    int curstate = ClusterAssn[spc][perturb_node];
    
    /*====== get nodes in the old and new clusters ========*/
    std::vector<int> new_cluster_node;
    std::vector<int> old_cluster_node;
    
    for (int i=0; i<snsNumNodes[spc]; i++) // here needs to go visit all nodes in a species, can improve by only looking at the old and new clusters
    {
        if (ClusterAssn[spc][i] == curstate)
            new_cluster_node.push_back(i);
        
        if (ClusterAssn[spc][i] == oldstate)
            old_cluster_node.push_back(i);
    }
    
    old_cluster_node.push_back(perturb_node); // add the perturbed node back to the old cluster
    
    int old_cluster_size = (int)old_cluster_node.size(); // old cluster size is +=1 here since it includes perturb_node
    int new_cluster_size = (int)new_cluster_node.size();
//    if (SA_counter % 10000 == 0)
//    {
        fprintf(stderr, "old_cluster: %d, old_cluster_size: %d, new_cluster: %d, new_cluster_size: %d\n", oldstate, old_cluster_size, curstate, new_cluster_size);
        fprintf(stderr, "spe: %d, total_num_nodes: %d, num_node_noise_cluster: %d\n", spc, snsNumNodes[spc], num_node_noise_in_cluster[spc]);
//    }
    
    
    int in_edge_before_perturb = 0, in_edge_after_perturb = 0;
    int out_edge_before_perturb = 0, out_edge_after_perturb = 0;
    
    
    if (curstate == 0) // new cluster is noise cluster, only compute delta cost in old cluster
    {
//        if (SA_counter % 10000 == 0)
            fprintf(stderr, "new cluster is noise cluster ======================\n");
        //delta cost = decrease of cost in the old cluster (only need to compute delta cost for the old cluster)
        in_edge_before_perturb = 0, in_edge_after_perturb = 0;
        out_edge_before_perturb = 0, out_edge_after_perturb = 0;
        
        for (int i = 0; i < old_cluster_size; i++)
        {
            for (list<int>::iterator it = sns[spc]->adjacencyList[old_cluster_node[i]].begin(); it!=sns[spc]->adjacencyList[old_cluster_node[i]].end(); it++)
            {
                int adj_node = *it;
                if (ClusterAssn[spc][adj_node] == oldstate or // nodes in old cluster
                    adj_node == perturb_node) // adj_node had perturb_node though it is now in new cluster
                    in_edge_before_perturb ++;
                
                if (ClusterAssn[spc][adj_node] == oldstate and // nodes in old cluster
                    old_cluster_node[i] != perturb_node) // exclude perturb_node
                    in_edge_after_perturb ++;
                
                if (ClusterAssn[spc][adj_node] != oldstate and // adj_node is not in old cluster
                    adj_node != perturb_node and // perturb_node is now in new state, so excluding it
                    ClusterAssn[spc][adj_node] != 0) // adj_node is not in noise cluster
                    out_edge_before_perturb++;
                
                if (old_cluster_node[i] != perturb_node and // inside cluster does not have perturb_node
                    ClusterAssn[spc][adj_node] != oldstate and // adj_node is not in old cluster
                    ClusterAssn[spc][adj_node] != 0) // adj_node is not in noise cluster (here includes perturb_node)
                    out_edge_after_perturb++;
            }
        }
        in_edge_before_perturb /= 2;
        in_edge_after_perturb /= 2;
        /*====== before perturb in old cluster ==========*/
        if (in_edge_before_perturb > 0)
        {
            int possible_in_edge_before_perturb = (old_cluster_size) * (old_cluster_size - 1) / 2;
            int possible_out_edge_before_perturb = ( snsNumNodes[spc] - num_node_noise_in_cluster[spc] - (old_cluster_size) ) * (old_cluster_size);
            
            double in_density_before_perturb = double(in_edge_before_perturb) /possible_in_edge_before_perturb;
            
            double out_density_before_perturb = double(out_edge_before_perturb) / possible_out_edge_before_perturb;
            
            double inOutRatio_before_perturb = (in_density_before_perturb / (in_density_before_perturb + out_density_before_perturb)) * double(old_cluster_size);
            
            delta_cost -= inOutRatio_before_perturb;
            
//            if (SA_counter % 10000 == 0)
                fprintf(stderr, "in_edge_before_perturb: %d, possible_in_edge_before_perturb: %d, out_edge_before_perturb: %d, possible_out_edge_before_perturb: %d, inOutRatio_before_perturb: %g, delta_cost: %g\n", in_edge_before_perturb, possible_in_edge_before_perturb, out_edge_before_perturb, possible_out_edge_before_perturb, inOutRatio_before_perturb, delta_cost);
        }
        /*====== after perturb in old cluster ==========*/
        if (in_edge_after_perturb > 0)
        {
            int possible_in_edge_after_perturb = (old_cluster_size - 1) * (old_cluster_size - 2) / 2;
            int possible_out_edge_after_perturb = ( snsNumNodes[spc] - num_node_noise_in_cluster[spc] - (old_cluster_size-1) ) * (old_cluster_size-1);
            
            double in_density_after_perturb = double(in_edge_after_perturb) / possible_in_edge_after_perturb;
            
            double out_density_after_perturb = double(out_edge_after_perturb) / possible_out_edge_after_perturb;
            
            double inOutRatio_after_perturb = (in_density_after_perturb / (in_density_after_perturb + out_density_after_perturb)) * double(old_cluster_size-1);
            
            delta_cost += inOutRatio_after_perturb;
            
//            if (SA_counter % 10000 == 0)
                fprintf(stderr, "in_edge_after_perturb: %d, possible_in_edge_after_perturb: %d, out_edge_after_perturb: %d, possible_out_edge_after_perturb: %d, inOutRatio_after_perturb: %g, delta_cost: %g\n", in_edge_after_perturb, possible_in_edge_after_perturb, out_edge_after_perturb, possible_out_edge_after_perturb, inOutRatio_after_perturb, delta_cost);
        }
    }
    else if (oldstate == 0) // the old cluster is noise cluster, only need to compute delta cost for the new cluster
    {
//        if (SA_counter % 10000 == 0)
            fprintf(stderr, "old cluster is noise cluster ======================\n");
        //delta cost = increase of cost in the new cluster (only need to compute delta cost for the new cluster)
        in_edge_before_perturb = 0, in_edge_after_perturb = 0;
        out_edge_before_perturb = 0, out_edge_after_perturb = 0;
        
        for (int i = 0; i < new_cluster_size; i++) // go through each node in new cluster
        {
            for (list<int>::iterator it = sns[spc]->adjacencyList[new_cluster_node[i]].begin(); it!=sns[spc]->adjacencyList[new_cluster_node[i]].end(); it++)
            {
                int adj_node = *it;
                if (ClusterAssn[spc][adj_node] == curstate)
                    in_edge_after_perturb ++;
                
                if (ClusterAssn[spc][adj_node] == curstate and
                    new_cluster_node[i] != perturb_node and // perturb_node was not in the new cluster
                    adj_node != perturb_node)
                    in_edge_before_perturb ++;
                
                if (ClusterAssn[spc][adj_node] != curstate and
                    ClusterAssn[spc][adj_node] != 0)
                    out_edge_after_perturb ++;
                
                if (ClusterAssn[spc][adj_node] != curstate and
                    new_cluster_node[i] != perturb_node and // perturb_node was not in the new cluster
                    ClusterAssn[spc][adj_node] != 0 and
                    adj_node != perturb_node) // adj_node is not perturb_node since it is now in new cluster instead of noise cluster, but before it was in noise cluster
                    out_edge_before_perturb ++;
            }
        }
        in_edge_before_perturb /= 2;
        in_edge_after_perturb /= 2;
        /*====== before perturb in new cluster ==========*/
        if (in_edge_before_perturb > 0) // new_cluster is not noise cluster
        {
            int possible_in_edge_before_perturb = (new_cluster_size - 1) * (new_cluster_size - 2) / 2;
            int possible_out_edge_before_perturb = ( snsNumNodes[spc] - num_node_noise_in_cluster[spc] - (new_cluster_size-1) ) * (new_cluster_size-1);
            
            double in_density_before_perturb = double(in_edge_before_perturb) / possible_in_edge_before_perturb;
            
            double out_density_before_perturb = double(out_edge_before_perturb) / possible_out_edge_before_perturb;
            
            double inOutRatio_before_perturb = ( in_density_before_perturb / (in_density_before_perturb + out_density_before_perturb) ) * double(new_cluster_size - 1);
            
            delta_cost -= inOutRatio_before_perturb;
            
//            if (SA_counter % 10000 == 0)
                fprintf(stderr, "in_edge_before_perturb: %d, possible_in_edge_before_perturb: %d, out_edge_before_perturb: %d, possible_out_edge_before_perturb: %d, inOutRatio_before_perturb: %g, delta_cost: %g\n", in_edge_before_perturb, possible_in_edge_before_perturb, out_edge_before_perturb, possible_out_edge_before_perturb, inOutRatio_before_perturb, delta_cost);
        }
        /*====== after perturb in new cluster ==========*/
        if (in_edge_after_perturb > 0)
        {
            int possible_in_edge_after_perturb = (new_cluster_size) * (new_cluster_size - 1) / 2;
            int possible_out_edge_after_perturb = ( snsNumNodes[spc] - num_node_noise_in_cluster[spc] - (new_cluster_size) ) * (new_cluster_size);
            
            double in_density_after_perturb = double(in_edge_after_perturb) / possible_in_edge_after_perturb;
            
            double out_density_after_perturb = double(out_edge_after_perturb) / possible_out_edge_after_perturb;
            
            double inOutRatio_after_perturb = (in_density_after_perturb / (in_density_after_perturb + out_density_after_perturb)) * double(new_cluster_size);
            
            delta_cost += inOutRatio_after_perturb;
            
//            if (SA_counter % 10000 == 0)
                fprintf(stderr, "in_edge_after_perturb: %d, possible_in_edge_after_perturb: %d, out_edge_after_perturb: %d, possible_out_edge_after_perturb: %d, inOutRatio_after_perturb: %g, delta_cost: %g\n", in_edge_after_perturb, possible_in_edge_after_perturb, out_edge_after_perturb, possible_out_edge_after_perturb, inOutRatio_after_perturb, delta_cost);
        }
        
    }
    else // when neither old or new cluster is noise cluster, both delta cost need to be computed
    {
//        if (SA_counter % 10000 == 0)
            fprintf(stderr, "neither old or new cluster is noise cluster ======================\n");
        /*====================== delta cost in new cluster ==========================*/
        //        printf("delta cost in new cluster======================\n");
        in_edge_before_perturb = 0, in_edge_after_perturb = 0;
        out_edge_before_perturb = 0, out_edge_after_perturb = 0;
        
        for (int i = 0; i < new_cluster_size; i++) // go through each node in new cluster
        {
            for (list<int>::iterator it = sns[spc]->adjacencyList[new_cluster_node[i]].begin(); it!=sns[spc]->adjacencyList[new_cluster_node[i]].end(); it++)
            {
                int adj_node = *it;
                if (ClusterAssn[spc][adj_node] == curstate)
                    in_edge_after_perturb ++;
                
                if (ClusterAssn[spc][adj_node] == curstate and
                    new_cluster_node[i] != perturb_node and
                    adj_node != perturb_node)
                    in_edge_before_perturb ++;
                
                if (ClusterAssn[spc][adj_node] != curstate and
                    ClusterAssn[spc][adj_node] != 0)
                    out_edge_after_perturb ++;
                
                if (ClusterAssn[spc][adj_node] != curstate and
                    ClusterAssn[spc][adj_node] != 0 and
                    new_cluster_node[i] != perturb_node)
                    out_edge_before_perturb ++;
                
                if (adj_node == perturb_node and
                    new_cluster_node[i] != perturb_node) // special case when adj_node is perturb_node and in new cluster now
                    out_edge_before_perturb ++;
            }
        }
        in_edge_before_perturb /= 2;
        in_edge_after_perturb /= 2;
        /*====== before perturb in new cluster ==========*/
        if (in_edge_before_perturb > 0) // new_cluster is not noise cluster
        {
            int possible_in_edge_before_perturb = (new_cluster_size - 1) * (new_cluster_size - 2) / 2;
            int possible_out_edge_before_perturb = ( snsNumNodes[spc] - num_node_noise_in_cluster[spc] - (new_cluster_size-1) ) * (new_cluster_size-1);
            
            double in_density_before_perturb = double(in_edge_before_perturb) / possible_in_edge_before_perturb;
            
            double out_density_before_perturb = double(out_edge_before_perturb) / possible_out_edge_before_perturb;
            
            double inOutRatio_before_perturb = ( in_density_before_perturb / (in_density_before_perturb + out_density_before_perturb) ) * double(new_cluster_size - 1);
            
            delta_cost -= inOutRatio_before_perturb;
            
//            if (SA_counter % 10000 == 0)
                fprintf(stderr, "in_edge_before_perturb: %d, possible_in_edge_before_perturb: %d, out_edge_before_perturb: %d, possible_out_edge_before_perturb: %d, inOutRatio_before_perturb: %g, delta_cost: %g\n", in_edge_before_perturb, possible_in_edge_before_perturb, out_edge_before_perturb, possible_out_edge_before_perturb, inOutRatio_before_perturb, delta_cost);
        }
        /*====== after perturb in new cluster ==========*/
        if (in_edge_after_perturb > 0)
        {
            int possible_in_edge_after_perturb = (new_cluster_size) * (new_cluster_size - 1) / 2;
            int possible_out_edge_after_perturb = ( snsNumNodes[spc] - num_node_noise_in_cluster[spc] - (new_cluster_size) ) * (new_cluster_size);
            
            double in_density_after_perturb = double(in_edge_after_perturb) / possible_in_edge_after_perturb;
            
            double out_density_after_perturb = double(out_edge_after_perturb) / possible_out_edge_after_perturb;
            
            double inOutRatio_after_perturb = (in_density_after_perturb / (in_density_after_perturb + out_density_after_perturb)) * double(new_cluster_size);
            
            delta_cost += inOutRatio_after_perturb;
            
//            if (SA_counter % 10000 == 0)
                fprintf(stderr, "in_edge_after_perturb: %d, possible_in_edge_after_perturb: %d, out_edge_after_perturb: %d, possible_out_edge_after_perturb: %d, inOutRatio_after_perturb: %g, delta_cost: %g\n", in_edge_after_perturb, possible_in_edge_after_perturb, out_edge_after_perturb, possible_out_edge_after_perturb, inOutRatio_after_perturb, delta_cost);
        }
        
        
        /*============================= delta cost in old cluster =========================*/
        //        printf("delta cost in old cluster======================\n");
        in_edge_before_perturb = 0, in_edge_after_perturb = 0;
        out_edge_before_perturb = 0, out_edge_after_perturb = 0;
        
        for (int i = 0; i < old_cluster_size; i++)
        {
            for (list<int>::iterator it = sns[spc]->adjacencyList[ old_cluster_node[i] ].begin(); it!=sns[spc]->adjacencyList[ old_cluster_node[i] ].end(); it++)
            {
                int adj_node = *it;
                if (ClusterAssn[spc][adj_node] == oldstate or
                    adj_node == perturb_node) // perturb_node is no longer in oldstate
                    in_edge_before_perturb ++;
                
                if (ClusterAssn[spc][adj_node] == oldstate and
                    old_cluster_node[i] != perturb_node)
                    in_edge_after_perturb ++;
                
                if (ClusterAssn[spc][adj_node] != oldstate and
                    adj_node != perturb_node and
                    ClusterAssn[spc][adj_node] != 0)
                    out_edge_before_perturb++;
                
                if (old_cluster_node[i] != perturb_node and
                    ClusterAssn[spc][adj_node] != oldstate and
                    ClusterAssn[spc][adj_node] != 0) // adj_node is not in noise cluster, but now perturb_node is in a new cluster instead of noise cluster
                    out_edge_after_perturb++;
            }
        }
        in_edge_before_perturb /= 2;
        in_edge_after_perturb /= 2;
        /*====== before perturb in old cluster ==========*/
        if (in_edge_before_perturb > 0)
        {
            int possible_in_edge_before_perturb = (old_cluster_size) * (old_cluster_size - 1) / 2;
            int possible_out_edge_before_perturb = ( snsNumNodes[spc] - num_node_noise_in_cluster[spc] - (old_cluster_size) ) * (old_cluster_size);
            
            double in_density_before_perturb = double(in_edge_before_perturb) /possible_in_edge_before_perturb;
            
            double out_density_before_perturb = double(out_edge_before_perturb) / possible_out_edge_before_perturb;
            
            double inOutRatio_before_perturb = (in_density_before_perturb / (in_density_before_perturb + out_density_before_perturb)) * double(old_cluster_size);
            
            delta_cost -= inOutRatio_before_perturb;
            
//            if (SA_counter % 10000 == 0)
                fprintf(stderr, "in_edge_before_perturb: %d, possible_in_edge_before_perturb: %d, out_edge_before_perturb: %d, possible_out_edge_before_perturb: %d, inOutRatio_before_perturb: %g, delta_cost: %g\n", in_edge_before_perturb, possible_in_edge_before_perturb, out_edge_before_perturb, possible_out_edge_before_perturb, inOutRatio_before_perturb, delta_cost);
        }
        /*====== after perturb in old cluster ==========*/
        if (in_edge_after_perturb > 0)
        {
            int possible_in_edge_after_perturb = (old_cluster_size - 1) * (old_cluster_size - 2) / 2;
            int possible_out_edge_after_perturb = ( snsNumNodes[spc] - num_node_noise_in_cluster[spc] - (old_cluster_size-1) ) * (old_cluster_size-1);
            
            double in_density_after_perturb = double(in_edge_after_perturb) / possible_in_edge_after_perturb;
            
            double out_density_after_perturb = double(out_edge_after_perturb) / possible_out_edge_after_perturb;
            
            double inOutRatio_after_perturb = (in_density_after_perturb / (in_density_after_perturb + out_density_after_perturb)) * double(old_cluster_size-1);
            
            delta_cost += inOutRatio_after_perturb;
            
//            if (SA_counter % 10000 == 0)
                fprintf(stderr, "in_edge_after_perturb: %d, possible_in_edge_after_perturb: %d, out_edge_after_perturb: %d, possible_out_edge_after_perturb: %d, inOutRatio_after_perturb: %g, delta_cost: %g\n", in_edge_after_perturb, possible_in_edge_after_perturb, out_edge_after_perturb, possible_out_edge_after_perturb, inOutRatio_after_perturb, delta_cost);
        }
    }
    return -delta_cost;
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
                double clus_score = 0.0;
                // #possible in-cluster edges
                int possible_in_edge = cluster_size * (cluster_size-1) / 2;
                // #possible out-cluster edges
//                int possible_out_edge = ( snsNumNodes[spc] - num_node_noise_in_cluster[spc] - (cluster_size) ) * (cluster_size);
                int possible_out_edge = ( snsNumNodes[spc] - noise_clus_size - (cluster_size) ) * (cluster_size);
                // in-density
                double total_cluster_inDensity = double(total_cluster_edge[cluster]/2) / double(possible_in_edge);
                // out-density
                double total_cluster_outDensity = double(total_cluster_outEdges[cluster]) / double(possible_out_edge);
                // cost of a spc = in-density * size / (in-density + out-density)
                clus_score = (total_cluster_inDensity / (total_cluster_inDensity + total_cluster_outDensity)) * double(cluster_size);
                spc_cost += clus_score;
                
//                fprintf(stderr, "cluster: %d, size: %d, in_edge: %d, out_edge: %d, possible_in_edge: %d, possible_out_edge: %d, cluster score: %g\n", cluster, cluster_size, total_cluster_edge[cluster]/2, total_cluster_outEdges[cluster], possible_in_edge, possible_out_edge, clus_score);
            }
        }
        // add spc_cost to total cost
        cost += spc_cost;
//        cost += (double)noise_clus_size * 0.7;
        
        // adaptive cost of noise cluster = #nodes in noise cluster * (spc_cost / size) -> adaptive weight for noise cluster
        cost += (double)noise_clus_size * (spc_cost / snsNumNodes[spc]);
        
        if (SA_counter % 1000 == 0)
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
    double orthterms = 0.0;
    double total_noise_nodes = 0.0;
    double orth_noise_term = 0.0;
    double total_gene_among_spc = 0.0;
    // loop through all spc
    for (int spc1=0; spc1 < numSpecies; spc1++)
    {
        // sum up total genes
        total_gene_among_spc += snsNumNodes[spc1];
        // sum up noise nodes
        total_noise_nodes += num_node_noise_in_cluster[spc1];
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
 assign perturbed node with its old state
 */
void Cluster::UndoPerturb(int **ClusterAssn, int new_state)
{
    for (int i = UndoLogSize; i>0; i--)
    {
        ClusterAssn[UndoLog[i-1].spc][UndoLog[i-1].node] = UndoLog[i-1].oldstate;
        if (new_state == 0) // if new cluster is noise cluster, to undo it, decrease count of nodes in noise cluster by 1
            num_node_noise_in_cluster[UndoLog[i-1].spc] --;
        if (UndoLog[i-1].oldstate == 0) // if old cluster is noise cluster, to undo it, increase count of nodes in noise cluster by 1
            num_node_noise_in_cluster[UndoLog[i-1].spc] ++;
    }
    
    UndoLogSize = 0;
}

/*======
 Perturb
 randomly assign a cluster to a random node
 */
pair<int, int> Cluster::Perturb(int **ClusterAssn, int numchanges = 1)
{
    int newstate = 0;
    int oldstate;
    for (int i=0; i<numchanges; i++)
    {
        // choose a random species
        int spc = rand() % numSpecies;
        // chose a random node
        int node = rand() % (snsNumNodes[spc]); // node is a uniq_ID instead of gene_ID
        oldstate = ClusterAssn[spc][node];
        
        while (1)
        {
            newstate = rand() % numClusters; // choose a randome new cluster
            /*================*/
            //            newstate = (rand() % (numClusters-1))+1;
            /*================*/
            if (newstate == oldstate) continue;
//            fprintf(stderr,"Per\tspe%d\tnode%s\told%d\tnew%d\n",spc,sns[spc]->nodeId2Name[node].c_str(),oldstate,newstate); // prev log
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
        if (oldstate == 0) // if old cluster is noise cluster, then decrease count of nodes in noise cluster by 1
            num_node_noise_in_cluster[UndoLog[UndoLogSize-1].spc] --;
        if (newstate == 0) // if new cluster is noise cluster, then increase count of nodes in noise cluster by 1
            num_node_noise_in_cluster[UndoLog[UndoLogSize-1].spc] ++;
    }
    return make_pair(newstate, oldstate);
}
