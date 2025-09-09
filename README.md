# op

Calculates the Billera-Holmes-Vogtmann geodesic distance 
(using the Owens-Provan algorithm) distance between trees:

It assumes you want tree distances between some reference tree and at 
least one other tree.

If you want geodesic distances (the default) between the first tree in 
trees.tre and all the others,

./op --treefile trees.tre

If you want geodesic distances (the default) between the first tree in 
trees.tre and all the others, and you want to see all the details,

./op --treefile trees.tre --quiet no

If you want to calculate the Frechet mean tree, variance, and 95% HPD interval, specify a file name prefix (e.g. "meantree") and the parameters that
determine the stopping criterion for determining the mean:
* k is the maximum number of iterations
* n  is the number of recent iterations to compare
* epsilon is the maximum amount by which every pairwise comparison of recent iterations must differ in order to stop. 
The following example uses values for k, n, and epsilon similar to those used by Brown and Owen (2020):

./op --treefile trees.tre --frechet meantree --frechet-k 1000000 --frechet-n 10 --frehet-epsilon 0.00001

The resulting output file will, in this case, be named "meantree.txt" and will contain comments (lines beginning with #) specifying the mean tree newick description. tree length, variance, and the 95% HPD interval lower and upper bounds. The file can also be executed in R to plot the kernel density of tree distances with 95% HPD shaded.

Other options can be seen by asking for help:

./op --help

Literature Cited

LJ Billera, SP Holmes, and K Vogtmann. 2001. Geometry of the space of phylogenetic trees.
Advances in Applied Mathematics 27:733-767.
[DOI:10.1006/aama.2001.0759](https://doi.org/10.1006/aama.2001.0759)

DG Brown and M Owen. 2020. Mean and variance of phylogenetic trees. Systematic Biology 69:139-154. [DOI:10.1093/sysbio/syz041](https://doi.org/10.1093/sysbio/syz041)

M Owens and JS Provan. 2011. A fast algorithm for computing geodesic distances
in tree space. IEEE/ACM Transactions on Computational Biology and Bioinformatics
8:2-12. [DOI:10.1109/TCBB.2010.3](https://doi.org/10.1109/TCBB.2010.3)

