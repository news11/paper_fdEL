Empirical likelihood based inference for functional means with application to wearable device data 
by Hsin-wen Chang and Ian W. McKeague


OVERVIEW
The documented R code for implementing our method and the NHANES occupation time data used in this paper are posted on Github (http://github.com/news11/fdEL) as an R package. This package is used in the following code named after different (sub)sections to reproduce the results in those (sub)sections of the paper.


REQUIREMENTS
The attached code uses `fdEL' package (http://github.com/news11/fdEL) hosted on Github. To install this package from Github, R package `devtool' is required. Other than these, R package `refund' and `fda' are used. The attached code installs and loads these packages.


RECOMMENDATION
***use similar write-up: Although it takes little time to achieve confidence bands and perform hypothesis testings suggested in the paper, reapeating those multiple times for simulation study purpose can be time consuming. Since the methods depend heavily on vector/matrix operations, we recommend using better BLAS libraries, like openblas or Intel(R) MKL to reduce computation time.


SECTION 4.1
This section performs hypothesis testings on simulated data. The simulation size, sample size, smoothness of curves should be adjusted as desired. 
***use similar write-up: Note that "Bs", the parametric bootstrap based band -- which is not a novel contribution of the paper but is included as a comparison -- takes the most computation time and therefore is turned off by default. It can be easily tuned on by including "Bs" for the types of bands to achieve.


SECTION 4.2
This section compares shapes of confidence bands and local coverage rates. The code was modified to be more concise than the original code used in the paper, by utilizing functions written later on.


SECTION 4.3
This section compares effects of different smoothing technique and those of more complicated data structure. This is the same original code used is the paper, and therefore, is lengthier than the code in Section 4.2 above.


SECTION 5
This section includes a real data example. Since the raw activity data from the NHANESs website is huge, the R package contains only the cleaned data and the code for the cleaning (in data-raw/DATASET.R).


Hsin-wen Chang, Ph.D. 
Assistant Research Fellow 
Institute of Statistical Science 
Academia Sinica 
Phone: +886-2-27875715 
Fax: +886-2-27886833 
Email: hwchang@stat.sinica.edu.tw

