Empirical likelihood based inference for functional means with application to wearable device data 
by Hsin-wen Chang and Ian W. McKeague


OVERVIEW
The following code is named after different (sub)sections to reproduce the results in those (sub)sections of the paper.


REQUIREMENTS
The attached code uses `fdEL' package from http://github.com/news11/fdEL, which includes documented R code for implementing the proposed EL methods, along with the NHANES data used in the paper. To install this package from Github, R package `devtool' is required. The attached code installs and loads these packages.


RECOMMENDATION
Although it takes little time to implement the methods in the paper, repeating those multiple times for simulation study purpose can be time consuming. We recommend using parallel computing through a server or cluster.


SECTION 3.1_1
This section constructs simultaneous confidence bands based on simulated data in the first example of Section 3.1 of the paper. The parameters for simulations can be adjusted to obtain the numbers in the first four columns of Table 1 of the paper.


SECTION 3.1_2
This section constructs simultaneous confidence bands based on simulated data in the second example of Section 3.1 of the paper. The parameters for simulations can be adjusted to obtain the numbers in the last four columns of Table 1 of the paper.


SECTION 3.2
This section conducts functional ANOVA tests based on simulated data in Section 3.2 of the paper. The parameters for simulations can be adjusted to obtain the numbers in Table 2 of the paper.


SECTION 4
This section analyzes real data from the NHANES study. The results in Section 4 of the paper can be reproduced using this file.


Hsin-wen Chang, Ph.D. 
Assistant Research Fellow 
Institute of Statistical Science 
Academia Sinica 
Phone: +886-2-27875715 
Fax: +886-2-27886833 
Email: hwchang@stat.sinica.edu.tw

