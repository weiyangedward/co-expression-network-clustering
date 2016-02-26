#include "cs-grn.h"
#include "ReadInputs.h"
//#include <cmath>

/*===========
 IsEdge(i, j)
 return true if edge weight between i and j > edgeWeightThreshold
 */
bool SpeciesNetwork::IsEdge(int i, int j)
{
    if (i >= numNodes || j >= numNodes)
        return false;
    
    if (std::abs(network[i][j]) > edgeWeightThreshold)
        return true;
    
    return false;
}

/*==================
 GetEdgeWeight(i, j)
 return edge weight between i and j
 */
double SpeciesNetwork::GetEdgeWeight(int i, int j)
{
    if (i >= numNodes || j >= numNodes)
        return false;
    
    return network[i][j];
}

/*==========
 hashGeneID
 function to convert gene_ID to uniq_ID
 */
void hashGeneID(SpeciesNetwork *sn, string gene_ID, int & uniqueNodes, char *fn)
{
    if (sn->nodeName2Id.find(gene_ID) == sn->nodeName2Id.end()) // if gene1 is not seen
    {
        fprintf(stderr,"Assigning id %d to node %s in species %s\n",uniqueNodes, gene_ID.c_str(), fn); // print out which gene is assigned what unique ID
        sn->nodeId2Name[uniqueNodes] = gene_ID; // create a new entry for (uniq_ID, gene_ID) in hashmap
        sn->nodeName2Id[gene_ID] = uniqueNodes; // create a new entry for (gene_ID, gene_ID) in hashmap
        uniqueNodes++; // increase uniq_ID by 1
        
        if (uniqueNodes > sn->numNodes) // error: if uniq_ID > number of nodes
        {
            printf("Error: In file %s, read more unique node ids (%d) than the %d expected\n", fn, uniqueNodes, sn->numNodes);
            exit(1);
        }
    }
}

/*============================
 ReadNetworkFile(network_file)
 read network file from each spe
 */
SpeciesNetwork *ReadNetworkFile(char *fn) // fn is a file name
{
    FILE *F = fopen(fn, "r"); // open file
    if (F==NULL) return NULL;
    
    SpeciesNetwork *sn = new SpeciesNetwork; // create a new network object 'sn' for a species
    
    char line[linelen]; // buffer for each line
    
    fgets(line, linelen, F); // first line: species name
    sscanf(line, "%s", (char *)&sn->spcname);
    
    
    fgets(line, linelen, F); // second line: number of nodes
    sscanf(line, "%d", &sn->numNodes);
    
    /*===== initialize network using a adj matrix filled by 0s ======*/
    sn->network = new double *[sn->numNodes];
    for (int i=0; i<sn->numNodes; i++)
    {
        sn->network[i] = new double[sn->numNodes];
        for (int j=0; j<sn->numNodes; j++)
            sn->network[i][j] = 0;
    }
    
    /*==== initialize network using adjacency list (arr of list) =====*/
    for (int i=0; i<sn->numNodes;i++)
    {
        list<int> adj_genes;
        sn->adjacencyList.push_back(adj_genes);
    }
    
    int uniqueNodes = 0; // a unique ID assigned each gene so that the all gene IDs will be mapped to 1,2,3,...,n
    /*====== read the rest of network file ======*/
    while (fgets(line, linelen, F))
    {
        char *id1 = new char[linelen]; // buff for gene1
        char *id2 = new char[linelen]; // buff for gene2
        double wt; // edge weight
        // wt = 1
        sscanf(line, "%s\t%s\t%lf",id1, id2, &wt);
        string str_id1(id1); // copy id1 to a str
        string str_id2(id2); // copy id2 to a str
        
        /*======== convert gene_ID to uniq_ID =========*/
        hashGeneID(sn, str_id1, uniqueNodes, fn); // convert gene1
        hashGeneID(sn, str_id2, uniqueNodes, fn); // convert gene2
        
        /*====== build network using adj matrix and adj list =======*/
        sn->network[sn->nodeName2Id[str_id1]][sn->nodeName2Id[str_id2]] = wt; // store edge weight between gene1 and gene2 in network adj matrix
        sn->adjacencyList[sn->nodeName2Id[str_id1]].push_back(sn->nodeName2Id[str_id2]); // store gene2 in the adj list of gene1, only do this for gene1 assuming that network file has both edges of: gene1->gene2 and gene2->gene1
    }
    fclose(F); // close file
    
    return sn; // return network 'sn' for a species
}


/*================
 ReadOrthologyFile
 read orth file
 */
Orthology *ReadOrthologyFile(char *orthfn, SpeciesNetwork **sns, int numSpc)
{
    // create space for the table
    Orthology *o = new Orthology;
    o->numSpecies = numSpc;
    o->orthtable = new unordered_map<NodePair *, int, NodePairHasher, eqNodePair> *[numSpc];
    // o->weighted_orth = new hash_map<int , double, NodePairHasher, eqNodePair> *[numSpc]; // hash table to store 1/d(#ortho genes from spe1 to spe2)
    o->weighted_orth = new unordered_map<int, double> *[numSpc];
    for (int i=0; i<numSpc; i++)
    {
        o->orthtable[i] = new unordered_map<NodePair *, int, NodePairHasher, eqNodePair> [numSpc];
        // o->weighted_orth[i] = new hash_map<int , double, NodePairHasher, eqNodePair> [numSpc];
        o->weighted_orth[i] = new unordered_map<int, double> [numSpc];
    }
    
    // create numerical ids for species names
    o->spcId2Name = new char *[numSpc];
    for (int i=0; i<numSpc; i++) {
        o->spcId2Name[i] = new char[linelen];
        strcpy(o->spcId2Name[i], sns[i]->spcname);
        o->spcName2Id[(char *)sns[i]->spcname] = i;
    }
    
    // read the file and populate table
    FILE *F = fopen(orthfn, "r");
    if (F==NULL) return NULL;
    
    char line[linelen];
    /* ortholog files have double-counted edges: spc1.gene1->spc2.gene2 and spc2.gene2->spc1.gene1 */
    while (fgets(line, linelen, F))
    {
        char spc1[linelen];
        char spc2[linelen];
        char id1[linelen];
        char id2[linelen];
        sscanf(line, "%s %s %s %s",(char *)&spc1, (char *)&spc2, (char *)&id1, (char *)&id2);
        string str_spc1(spc1);
        string str_spc2(spc2);
        string str_id1(id1);
        string str_id2(id2);
        // convert spc ids to numerics
        if (o->spcName2Id.find(str_spc1) == o->spcName2Id.end())
        {
            printf("Error: orthology file %s has species name that wasnt seen in network files\n", orthfn);
            exit(1);
        }
        int spc1id = o->spcName2Id[str_spc1];
        if (o->spcName2Id.find(str_spc2) == o->spcName2Id.end())
        {
            printf("Error: orthology file %s has species name that wasnt seen in network files\n", orthfn);
            exit(1);
        }
        int spc2id = o->spcName2Id[str_spc2];
        
        // convert gene ids to numerics
        if (sns[spc1id]->nodeName2Id.find(str_id1) == sns[spc1id]->nodeName2Id.end()) { continue; }
        if (sns[spc2id]->nodeName2Id.find(str_id2) == sns[spc2id]->nodeName2Id.end()) { continue; }
        int gene1id = sns[spc1id]->nodeName2Id[str_id1];
        int gene2id = sns[spc2id]->nodeName2Id[str_id2];
        
        NodePair *np = new NodePair;
        np->id1 = gene1id;
        np->id2 = gene2id;
        o->orthtable[spc1id][spc2id][np] = 1;
        // only count edges from one direction: spc1->spc2
        o->weighted_orth[spc1id][spc2id][gene1id] += 1.0;
        // fprintf(stderr, "gene id: %d,  orth effect: %f\n", gene1id,  o->weighted_orth[spc1id][spc2id][gene1id]);
    }
    fclose(F);
    
    /* update ortholog edges from counts to ratio */
    // loop through all spc
    for (int spc1=0; spc1 < numSpc; spc1++)
    {
        // loop through spc+1, avoid spc2->spc1
        for (int spc2=spc1+1; spc2 < numSpc; spc2++)
        {
            unordered_map<int, double>::iterator iter;
            for ( iter = o->weighted_orth[spc1][spc2].begin(); iter != o->weighted_orth[spc1][spc2].end(); ++iter)
            {
                // second is ortholog edges counts
                double orth_effect = 1.0 / iter->second;
                o->weighted_orth[spc1][spc2][iter->first] = orth_effect;
                // fprintf(stderr, "orth effect %f\n", o->weighted_orth[spc1][spc2][iter->first]);
            }
        }
    }
    
    return o;
}


