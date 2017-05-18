# synthetic-promoter-design
Various approaches for designing complex synthetic promoters using PWM from the CISBP database

During my post-doctoral training, I came accross and complex problem: we need to have a reporter construct that would be
recognized and bound by  a set of 78 transcription factors in vivo. This "super promoter" would then be a perfect tool
for a high-throughput reporter assay.

So we set to develop the optimal approach for this task.

The promoter would need to satisfy a few criteria:
- Have sequences likely to be bound by the 78 TFs 
- satisfy any known positional requirements for specific TFs
- Be as short as possible: TF motifs far from the TSS are less likely to drive efficient transcription.
- Have the motifs in the right orientation: every motif has a reverse complement, but not every TF is bidirectional!
- Not be digested by the restriction enzymes used for cloning the synthetic promoter in desired vector

Using PWM files from the CISBP database, one can derive the most optimal sequence for a given promoter. 
However, each of these files are derived from different experimental conditions, and thus have variations in
size, probabilities, and orientations. 

Using known CORE motifs from the litterature , I was able to make scripts to generate a PWM in a correct orientation
for each TF (pwmorient.py), and then determine which is the most representative PWM among all for a given TF (pwmunif.py)

Using these unified PWMs, I could generate promoters based on 3 approaches:
1. Most probable sequence from each PWM (highest LOG-likelihood) in a random order. The script can generate many versions (topsequence
2. Similar to 1, but goes through iterations to find possible overlaps between motifs to make a shorter promoter. (topsequenceoverlap.py)
3. Use of a probabilistic approach to derive a promoter with "good enough" sequences, with a reduction in length.  (montecarlominiseed /montecarloseeded.py)

I also generated a script (promoterscan.py) to quantify how a given generated promoter performs. 

We then also set out to generate simpler promoters containing only a singlet TFBS 
