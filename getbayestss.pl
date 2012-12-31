#!/usr/bin/perl
#given a set of transcription start site estimates (or really any set of measures of central tendency)
#given prior prob distribution (e.g., empirically estimated)
#gives the posterior prob distribution of the estimate taking into account the set of measures and prior

use strict;
use warnings;

#declare variables
my (@genelist, @logprior) = ();
my $genelist = '';


open(FILE, "<co2_tss_preds_with_strands_test.txt" ) or die "Can't open co2_tsss.txt";
    while( <FILE> ) {
        @genelist = <FILE>;
    }
close (FILE);

open(FILE, "<bayes_gld_tss_log_prior.txt" ) or die "Can't open prior.txt";
    while( <FILE> ) {
        @logprior = <FILE>;
		chomp (@logprior);
    }
close (FILE);

foreach my $line (@genelist){
	my (@TSS, @TSS_w, @loglik, @logpost, @post, @nonzero_TSS, @nonzero_TSS_w) = ();
	my ($index, $TSS_estimate, $position, $weight, $nonzero_TSS, $estimate, $call, $nonzero_TSS_w, $maxvalue) = 0;
	#estimated from sample variance...
	my $sigma = 20;
	#upstream-most 5' utr position we're considering, note "really" 400 but perl indexes from 0
	my $upend = 399;
	#number of samples, note really 3 for now but perl indexes from 0...
	my $samples = 2;
	#list here if looking to visualize the distributions of a particular gene's TSS probs
	my $gene_of_interest = "GBAA0789";
	#TSS = empirical transcription start site, w = weight (sharpness), $TLSS = ncbi translation start site
	my ($genename, $TSS1, $TSS2, $TSS3, $TSS_w1, $TSS_w2, $TSS_w3, $TLSS, $strand) = split (' ', $line);
	@TSS = ($TSS1, $TSS2, $TSS3);
	@TSS_w = ($TSS_w1, $TSS_w2, $TSS_w3);
	
	#remove nonzero TSS's and weights, and put them into reference positions
	for ($weight = 0; $weight <= $samples; $weight++) {
		if ($TSS_w[$weight] > 0){
			my $outcome_TSS = abs ($TLSS - $TSS[$weight]);
			push (@nonzero_TSS, $outcome_TSS);
			push (@nonzero_TSS_w, $TSS_w[$weight]);
		}
	}
	
	#get length of nonzero_TSS array. if less than zero, no elements to consider, so skip the rest.
	$nonzero_TSS = @nonzero_TSS;
	if ($nonzero_TSS > 0){
		
		my ($k, $sample_var, $mu, $i, $sum, $S_real) = (0,0,0,0,0,0);
		
		#get the mean of the TSS positions
		for ($i = 0; $i < $nonzero_TSS; $i++){
			$sum += $nonzero_TSS[$i];
			}
		$mu = ($sum / $nonzero_TSS);	
		
		print "initially S is $sample_var\n";
	
		#get sigma estimate based on shrinkage estimator of the pop var, this is empirical bayes
		for ($k = 0; $k < $nonzero_TSS; $k++){
			$sample_var += (($nonzero_TSS[$k] - $mu)**2);
			#print "$nonzero_TSS[$k]\n";
			#print "$sample_var\n";
			#print "$k\n";
			#print "$nonzero_TSS\n";
			}		
		$S_real = sqrt($sample_var * (1 / ($nonzero_TSS - 1)));
		my $sigma_est = (((1 / $nonzero_TSS) * $sigma) + ((($nonzero_TSS - 1) / $nonzero_TSS) * $S_real)); 
		#print "$mu\t$sample_var\t$S_real\t$sigma_est\n";
	
		#for each possible TSS estimate, determine likelihood (w/ weighted avg of weights), add prior, get posterior
		for ($estimate = 0; $estimate <= $upend; $estimate++) {
			my ($loglik, $loglik_unweighted, $sum_weights, $logpost) = 0;
			for ($call = 0; $call <= $nonzero_TSS; $call++){
				$loglik_unweighted -= (($nonzero_TSS_w[$call] * (($nonzero_TSS[$call] - $estimate)**2)) / (2 * ($sigma_est**2)));
				#$sum_weights += $nonzero_TSS_w[$call];
			}
			#$loglik = ($loglik_unweighted / $sum_weights); 
			$loglik = $loglik_unweighted; #get rid of if go back to weights, obvi
			$logpost = $loglik_unweighted + $logprior[$estimate];
			push (@loglik, $loglik);
			push (@logpost, $logpost);
		}	
		
		#find max of log posterior and note the index, although of course it *could* be unchanged from 0, too
		$maxvalue = $logpost[0];
		for ($estimate = 0; $estimate <= $upend; $estimate++) {
			if ($logpost[$estimate] > $maxvalue){
				$maxvalue = $logpost[$estimate];
				$index = $estimate;
			}	
		}
	
		#depending on strand, add or subtract index from TLSS to find estimate that maximizes posterior density
		if ($strand == 1){
			$TSS_estimate = $TLSS - $index;
		}	
		else {
			$TSS_estimate = $TLSS + $index;
		}
	
		#if looking for a particular gene to visualize the distributions of, prep that here
		if ($genename =~ /$gene_of_interest/){
			for ($position = 0; $position <= $upend; $position++){
				my ($real_position, $real_logprior, $real_loglik, $real_logpost) = 0;
				$real_position = $position + 1;
				$real_logprior = $logprior[$position];
				$real_loglik = $loglik[$position];
				$real_logpost = $logpost[$position];
				open (OUT, ">>bayes_dstn_gene_of_interest.txt");
				print (OUT "$real_position\t$real_logprior\t$real_loglik\t$real_logpost\n");
				close (OUT);
			}
		}
	}
	
	#print stuff to output file	
	open (OUT, ">>bayes_tss_estimates.txt");
	print (OUT "$genename\t$TSS_estimate\t$maxvalue\n");
	close (OUT);
	
}