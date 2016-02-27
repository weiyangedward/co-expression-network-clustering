#!/usr/bin/perl

## This script evaluates clustering results with 'noise' cluster presented

use strict;

die "Usage: perl $0 ClusteringOutput RealCluster RealNetwork\n" unless @ARGV==3;

my $output = shift (@ARGV); ## one file with given clustering (output from program)
my $real = shift(@ARGV); ## one file with real clustering 
my $realnw = shift(@ARGV); ## one file with real network

##===== read clustering Output =====##
my %outclus;
open (OUT, "<$output");
my $hdr = <OUT>;
while (<OUT>) 
{
    my ($d1, $species, $d2, $gene, $d3, $cluster) = split(/\s+/);
    # warn "$species\t$gene\t$cluster\n";
    $outclus{$species}{$gene} = $cluster; ## outclus{0}{200} = 29
}
close(OUT);

##===== read Real Cluster =====##
my %realclus;
open (REAL, "<$real");
$hdr = <REAL>;
while (<REAL>) 
{
    my ($d1, $species, $d2, $gene, $d3, $cluster) = split(/\s+/);
    # warn "$species\t$gene\t$cluster\n";
    $realclus{$species}{$gene} = $cluster; ## realclus{0}{559} = 17
}
close(REAL);

##===== Random cluster =====##
my %randclus;
foreach my $spc (keys %outclus) 
{
    my %rgene;
    my @genes = (keys %{$outclus{$spc}});
    fisher_yates_shuffle(\@genes);
    my @ogenes = (keys %{$outclus{$spc}});
    for (my $i=0; $i<=$#ogenes; $i++) 
    {
		$rgene{$ogenes[$i]} = $genes[$i];
    }
    foreach my $g (keys %{$outclus{$spc}}) 
    {
		$randclus{$spc}{$rgene{$g}} = $outclus{$spc}{$g};
    }
}

##===== read Real Network =====##
my %network_total_edge_2x;
my %realnetwork;
my %nodeDegree; ## degree of node in a species
my %total_spc_node_degree;
open (REAL, "<$realnw");
$hdr = <REAL>;
while (<REAL>) 
{
    my ($species, $g1, $g2, $wt) = split(/\s+/);
    $realnetwork{$species}{$g1}{$g2} = $wt; ## realnetwork{0}{559}{365} = 1, this network only counts within species edges, NOT orthologous edges
    $nodeDegree{$species}{$g1} ++;
    $network_total_edge_2x{$species} ++;
}
close(REAL);

for my $spc (keys %nodeDegree)
{
	for my $g (keys %{$nodeDegree{$spc}})
    {
		$total_spc_node_degree{$spc} += $nodeDegree{$spc}{$g};
	}
}


print "Examining if real co-clustered pairs are predicted as such\n";
ComputeSensitivity(\%realclus, \%outclus);
print "Examining if predicted co-clustered pairs are really so\n";
ComputeSensitivity(\%outclus, \%realclus);
print "Examining if real co-clustered pairs are predicted as such by chance\n";
ComputeSensitivity(\%realclus, \%randclus);

print "Descriptive statistics of real clusters\n";
DescriptiveStats(\%realclus);
print "Descriptive statistics of predicted clusters\n";
DescriptiveStats(\%outclus);
##============================ 
#evaluate clustering using NMI
sub ComputeSensitivity {
	## clus1ptr = real; clus2ptr = pred
    my ($clus1ptr, $clus2ptr) = @_;

    my ($trueIn, $falseIn, $trueCross, $falseCross) = (0,0,0,0);
    my $NMI = 0;
    ## for each species
    foreach my $spc (keys %{$clus1ptr}) 
    {
    	my $total_node = 0;
    	my %clusterIntercept_real_to_found = ();
    	my %clusterIntercept_found = ();
    	my %clusterIntercept_real = ();
		## within cluster gene pairs
		my @genes = keys %{$$clus1ptr{$spc}};
		foreach my $g1 (@genes) 
        {
            ## cluste1
            my $real_cluster = $$clus1ptr{$spc}{$g1};
            ## cluster2
			my $found_cluster = $$clus2ptr{$spc}{$g1};
			$clusterIntercept_real_to_found{$real_cluster}{$found_cluster} ++;
			$clusterIntercept_found{$found_cluster} ++;
			$clusterIntercept_real{$real_cluster} ++;
			$total_node ++;
		    foreach my $g2 (@genes) 
            {
				if ($g1 eq $g2) { next; }
				if ($g1 > $g2) { next; }
				## in real cluster, g1 and g2 are in the same cluster
				if ($$clus1ptr{$spc}{$g1} eq $$clus1ptr{$spc}{$g2}) 
                {
					## in pred cluster, g1 and g2 are in the same cluster
				    if ($$clus2ptr{$spc}{$g1} eq $$clus2ptr{$spc}{$g2}) 
                    { 
				    	$trueIn ++ ; 
				    }
				    else 
                    { 
				    	$falseIn ++; 
				    }
				}
				else 
                {	## in real cluster, g1 and g2 are in diff cluster
					## in pred cluster, g1 and g2 are in the same cluster
				    if ($$clus2ptr{$spc}{$g1} eq $$clus2ptr{$spc}{$g2}) 
                    { 
                        $falseCross ++ ; 
                    }
				    else 
                    { 
                        $trueCross ++; 
                    }
				}
		    }
		}
		my $infomation = 0;
		for my $real (sort keys %clusterIntercept_real_to_found)
        {
			for my $found (sort keys %{$clusterIntercept_real_to_found{$real}})
            {
				my $N_ij = $clusterIntercept_real_to_found{$real}{$found};
				if ($N_ij > 0)
                {
					$infomation += $N_ij * log(($N_ij * $total_node) / ($clusterIntercept_real{$real} * $clusterIntercept_found{$found}));
				}
			}
		}
		$infomation = -2 * $infomation;

		my $sum_real = 0;
		for my $real (sort keys %clusterIntercept_real)
        {
			$sum_real += $clusterIntercept_real{$real} * log($clusterIntercept_real{$real} / $total_node);
		}

		my $sum_found = 0;
		for my $found (sort keys %clusterIntercept_found)
        {
			$sum_found += $clusterIntercept_found{$found} * log($clusterIntercept_found{$found} / $total_node);
		}
		$NMI += $infomation / ($sum_real + $sum_found);
    }
    $NMI = $NMI / 3;
    print "TrueIn\tFalseIn\tTrueCross\tFalseCross\tNMI\n";
    my $score = (($trueIn/($trueIn + $falseIn)) + ($trueCross/($trueCross+$falseCross)))/2;
    my $rindex = (($trueIn+$trueCross)/($trueIn+$trueCross+$falseIn+$falseCross));
    print "$trueIn\t$falseIn\t$trueCross\t$falseCross\t$NMI\n";
    print "Score = $score\tRand Index = $rindex\n";
}

##====================================== 
#print in/out density and other measures 
# noise cluster is not counted in computing these measures 
sub DescriptiveStats 
{
    
    ## format: clus1ptr{spc}{gene} = clus_id
    my ($clus1ptr) = @_;
    my %cs;
    my $clusterNum = 0;
    ## loop through all spc
    foreach my $spc (keys %{$clus1ptr}) 
    {
		## get all clus_id
		my @clusters = values %{$$clus1ptr{$spc}};
		## loop through all clus_id
		foreach my $c (@clusters) 
		{ 
			## format: cs{c} = 1, where c = 1,2,3
			$cs{$c} = 1;
		}
		## total #clust
		$clusterNum = keys %cs; ## cluster number
    }
    
    my $sumindensity = 0;
    my $countindensity = 0;
    my $sumoutdensity = 0;
    my $countoutdensity = 0;
    my $suminoutratio = 0;
    my $sum_ave_modularity = 0;
    my %total_cluster_edge = ();
	my %total_cluster_modularity = ();
	my %total_cluster_modularity_weighted = ();
	my %total_cluster_modularity_bySize = ();
	my %total_cluster_degree = ();
	my %cluster_size = (); ## cluster size
	my $sum_modularity_cluster_spc = 0;
	my $sum_modularity_cluster_spc_weighted = 0;
	my $sum_modularity_cluster_spc_bySize = 0;

	my %indensity = ();
	my %outdensity = ();
	my %outdegree = ();
	my %incount = ();
	my %outcount = ();
	my %weighted_in_out_ratio = ();
	my $sum_weighted_in_out_ratio = 0;

	# my $sum_product_verDen_verIntro = 0;

	## for each Cluster
    foreach my $c (keys %cs) 
    {
		## for each Species
		foreach my $spc (keys %{$clus1ptr}) 
        {

		    my @genes = keys %{$$clus1ptr{$spc}};
		    $cluster_size{$spc}{$c} = 0;

		    $indensity{$spc}{$c} = 0;
		    $outdensity{$spc}{$c} = 0;
		    $incount{$spc}{$c} = 0;
		    $outcount{$spc}{$c} = 0;
		   
		    # my $verIntrovert = 0;
		    # my $product_verDen_verIntro = 0;
		    # my $ave_modularity = 0;

		    foreach my $g1 (@genes) {
		    	## g1 in c
				if ($$clus1ptr{$spc}{$g1} eq $c) { 
					## sum genes in a cluster
					$cluster_size{$spc}{$c}++;
					## sum total degree in cluster
			    	$total_cluster_degree{$spc}{$c} += $nodeDegree{$spc}{$g1};

			    	my $inDegree = 0;
			    	foreach my $g2 (@genes) 
                    {
			    		## i /= j
						if ($g1 == $g2) { next; } 
						
						## g2 in c
						if ($$clus1ptr{$spc}{$g2} eq $c) 
                        { 

				    		## if there is an edge i,j
				    		if (defined($realnetwork{$spc}{$g1}{$g2})) 
                            {
				    			## increase total edges in cluster  
				    			$indensity{$spc}{$c}++;
				    			$inDegree ++;
				    			## sum total edges in cluster
				    			$total_cluster_edge{$spc}{$c} ++;
				    			## sum total modularity in cluster
				  				$total_cluster_modularity{$spc}{$c} += (1 / $network_total_edge_2x{$spc}) - (($nodeDegree{$spc}{$g1} * $nodeDegree{$spc}{$g2}) / ($network_total_edge_2x{$spc} ** 2));
				    		}
				    		else{
				    			## add penalty to nodes without edges
				    			$total_cluster_modularity{$spc}{$c} += - (($nodeDegree{$spc}{$g1} * $nodeDegree{$spc}{$g2}) / ($network_total_edge_2x{$spc} ** 2));
				    		}
				    		## sum total potential edges in-cluster
				    		$incount{$spc}{$c}++; 
						}
						else { ## g2 not in c
				    		if (defined($realnetwork{$spc}{$g1}{$g2})) {
				    			## sum total out-edges
				    			$outdensity{$spc}{$c}++;
				    			$outdegree{$spc}{$c} ++; 
				    		}
				    		## sum total potential edges out-cluster
				    		$outcount{$spc}{$c}++; 
						}
			    	}
	
				}
		    }
		    ## weight modularity by "#total edges in cluster" / "#total degree in cluster"
		    if ($total_cluster_degree{$spc}{$c} >= 1) 
            {
		    	$total_cluster_modularity_bySize{$spc}{$c} = $total_cluster_modularity{$spc}{$c} / $cluster_size{$spc}{$c};
		    	$total_cluster_modularity_weighted{$spc}{$c} = $total_cluster_modularity{$spc}{$c} * ($total_cluster_edge{$spc}{$c} / $total_cluster_degree{$spc}{$c});
		    	
		    }
		    else 
            {
		    	$total_cluster_modularity_bySize{$spc}{$c} = $total_cluster_modularity{$spc}{$c} / 1;
		    	$total_cluster_modularity_weighted{$spc}{$c} = $total_cluster_modularity{$spc}{$c} * ($total_cluster_edge{$spc}{$c} / 1);
		    }
		    ## sum modularity over clusters and spc
		    $sum_modularity_cluster_spc += $total_cluster_modularity{$spc}{$c};
		    $sum_modularity_cluster_spc_bySize += $total_cluster_modularity_bySize{$spc}{$c};
		    $sum_modularity_cluster_spc_weighted += $total_cluster_modularity_weighted{$spc}{$c};

		    ## total nodes in a cluster >= 2
		    if ($incount{$spc}{$c} > 0) 
            {
		    	$indensity{$spc}{$c} = $indensity{$spc}{$c} / $incount{$spc}{$c};
		    	$outdegree{$spc}{$c} = $outdegree{$spc}{$c} / $incount{$spc}{$c};
		    } 

		    if ($outcount{$spc}{$c} > 0) 
            { 
		    	$outdensity{$spc}{$c} = $outdensity{$spc}{$c} / $outcount{$spc}{$c}; 
		    } 

		    if ($incount{$spc}{$c} > 0)
            {
		    	$weighted_in_out_ratio{$spc}{$c} = ($indensity{$spc}{$c} / ($indensity{$spc}{$c} + $outdensity{$spc}{$c})) * ($cluster_size{$spc}{$c});
		    }
		    print "Cluster $c\tSpecies $spc\t Size $cluster_size{$spc}{$c} \t In-density $indensity{$spc}{$c}\tOut-density $outdensity{$spc}{$c}\t Cluster modularity $total_cluster_modularity{$spc}{$c} \t Cluster modularity weighted $total_cluster_modularity_weighted{$spc}{$c} \t Cluster modularity by Size $total_cluster_modularity_bySize{$spc}{$c} \t Total in-cluster edge $total_cluster_edge{$spc}{$c}\t Total in-cluster degree $total_cluster_degree{$spc}{$c} \t Weighted in/out density ratio $weighted_in_out_ratio{$spc}{$c}\n";
            
            ## skip noise cluster when computing the average
            if ($c != 0)
            {
		        $sumindensity += $cluster_size{$spc}{$c} * $indensity{$spc}{$c};
            
		        $countindensity += $cluster_size{$spc}{$c};
		        $sumoutdensity += $cluster_size{$spc}{$c} * $outdensity{$spc}{$c};
		        $countoutdensity += $cluster_size{$spc}{$c};
		        $sum_weighted_in_out_ratio += $weighted_in_out_ratio{$spc}{$c};
            }
		}
    }
    $sumindensity /= $countindensity;
    $sumoutdensity /= $countoutdensity;
    my $out_in_ratio = $sumoutdensity/$sumindensity;
    # $sum_product_verDen_verIntro /= ($clusterNum * 3);
    print "All clusters on average have in-density of $sumindensity and out-density of $sumoutdensity and Cluster num = $clusterNum\n";
}


sub fisher_yates_shuffle {
    my $array = shift;
    my $i;
    for ($i = @$array; --$i; ) {
	my $j = int rand ($i+1);
	next if $i == $j;
	@$array[$i,$j] = @$array[$j,$i];
    }
}
