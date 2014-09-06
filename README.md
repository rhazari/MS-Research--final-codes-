Master's Thesis Research codes
==============================

Analyzes the load-based performances of 3 different Connected Dominating Set algorithm. 

1. ECDS (Essential Connected Dominating set) algorithm is an approxiamte MCDS (Minimum Connected Dominating Set) algorithm.

2. 2-CDS (2 connectivity CDS) algorithm attempts to provide atleast 2 node disjoint paths between any two source-destination terminals.

3. LoB-CDS (Load based CDS) algorithm attempts to improve on the ECDS algorithm, by identifying terminals having load-factor values (bottleneck terminals) and attempts to add terminals from its neigbourhood to reduce its traffic load.

All the 3 CDS algorithms needs to be complemented with appropriate routing metrics to move traffic away from the bottleneck terminals.
