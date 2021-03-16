# Computing gene knockout strategies (MCS) for different degrees of growth-coupled production

Supplemetary code for:

Systematizing different notions of growth-coupled product synthesis and a single framework for computing corresponding strain designs

Philipp Schneider, Radhakrishnan Mahadevan, Steffen Klamt

2021/03/16

corresponding author: klamt@mpi-magdeburg.mpg.de


[Download supplementary code](https://github.com/klamt-Lab/MCS_growth-coupling/releases)

[Download CNA](https://www2.mpi-magdeburg.mpg.de/projects/cna/download.html)

Added Features (CellNetAnalyzer):
---------------

1.  Definition of optimality constraints for describing desired and undesired behavior
    (in addition to static linear constraints of the form **T**·**r** ≤ **t** and **D**·**r** ≤ **d**).
    Such constraints can be used to define MCS setups for potentially and weakly growth-coupled product synthesis.
    For potentially growth-coupled production, one defines desired flux states such that growth rate optimality is demanded
    together minima for the attainable growth rate and product synthesis rate. For weakly growth-coupled product synthesis, 
    one uses a target system in which growth rate optimality is selected together with the static constraint r_P = 0. 
    MCS will seek to eliminate where maximum growth and zero production coincide.

Software Requirements:
----------------------

1.  MATLAB 2016b® or later

2.  IBM ILOG® CPLEX® 12.7, 12.8, 12.9 or 12.10 (Make sure to use compatible CPLEX® and MATLAB® versions. Version 12.10 is recommended. 
    CPLEX® 20.1 DOES NOT WORK PROVIDE A MATLAB API AND CANNOT BE USED HERE)

3.  [CellNetAnalyzer](https://www2.mpi-magdeburg.mpg.de/projects/cna/download.html)2021.1 or later 

4.  Set up *CellNetAnalyzer* to access the CPLEX-Matlab-API (as described by *CellNetAnalyzer* manual)
    
Getting Started:
----------------------
1. Download this project to your computer (see [release page](https://github.com/klamt-Lab/MCS_growth-coupling/releases)) and extract all files.
2. Download the *E. coli* model iML1515 from the [Bigg-database](http://bigg.ucsd.edu/static/models/iML1515.mat). Place it in the project folder.
3. Start MATLAB.
4. Add the main directory of your installation of the *CellNetAnalyzer* toolbox to your MATLAB path.
5. Add the project folder and the subdirectory 'functions' to your MATLAB path.

Script Files:
-------------

1. **MCS_1_coupling_degrees.m**

   Computes **exhaustively all minimal cut sets up to three gene knockouts for the potentially, weakly and strongly growth-coupled and substrate-uptake coupled 
   production of ethanol** with *E. coli*. A subnetwork/core network (ca. 600 reactions) of the iML1515 is used to shorten the computation runtime.
   Finally, the relationship between the MCS sets of the different coupling types is plotted. The user can set/unset the flag to de-/activate the 
   minimum ATP maintenance demand.

2. **MCS_2_smallest.m** 

   Computes **the smallest minimal cut set** for the potentially, weakly and strongly growth-coupled and substrate-uptake coupled synthesis of 12 different products
   with *E. coli*. The genome-scale model iML1515 is used. For heterologous products, the pathways are added automatically.
   As the runtime strongly depends on the random seed used in the MCS computation, multiple computations are run with a time
   limit. To avoid memory problems, the computation is run in a seperate MATLAB instance. By default, the script runs with the following
   settings: (1) ethanol production (2) weak growth-coupled production (3) 6 iterations (4) 4 hours per computation (5) run each computation
   in a new process each (6) maximum of 60 knockouts (7) return after finding 1 solution. The user can change these settings to compute
   MCS for different products or coupling strengths.

3. **MCS_3_any.m** 

   Computes **a random minimal cut set** for the potentially, weakly and strongly growth-coupled and substrate-uptake coupled synthesis of 12 different products
   with *E. coli*. The genome-scale model iML1515 is used. For heterologous products, the pathways are added automatically.
   As the runtime strongly depends on the random seed used in the MCS computation, multiple computations are run with a time
   limit. To avoid memory problems, the computation is run in a seperate MATLAB instance. By default, the script runs with the following
   settings: (1) ethanol production (2) weak growth-coupled production (3) 12 iterations (4) 2 hours per computation (5) run each computation
   in a new process each (6) maximum of 60 knockouts (7) return after finding 1 solution. The user can change these settings to compute
   MCS for different products or coupling strengths.

4. **MCS_4_ACP.m** 

   Computes **exhaustively all minimal cut sets up to three gene knockouts for the strongly growth-coupled, ATP-coupled and substrate-uptake coupled 
   production of ethanol** with *E. coli*. A subnetwork/core network (ca. 600 reactions) of the iML1515 is used to shorten the computation runtime.
   The minimum ATP maintenace demand from the original *iML1515* model is omitted. Finally, the relationships between the MCS sets of the different 
   coupling types are plotted.

Other functions required for the scripts above:
-----------------------------------------------

5. functions/**load_pathway.m** 
   Contains the pathways for all heterologous products. This function returns the species and reactions that need to be added to the iML1515.

6. functions/**block_non_standard_products.m**
   Shut metabolite exchanges for all atypical products of *E. coli*.

7. functions/**check_mass_balance.m**

8. functions/**compare_mcs_sets.m**
   Between two sets of MCS, identify matching MCS, MCS subsets and supersets.

9. functions/**compute_MCS_ext**
   MCS computation function that is run by a new MATLAB instance.

10. functions/**plot_mcs_relationships.m**

Model files:
-------------

11. **core.mat** - Indicates species and reactions within the iML1515 that are part of *E. coli*'s core metabolism

12. **iML1515geneNames.mat** - contains a map of gene names and b-numbers. This helps to generate a user friendly output.

13. **iML1515.mat** - **Required for computation, but not provided in this repository. Please download from the **
                        [Bigg-database](http://bigg.ucsd.edu/static/models/iML1515.mat)


Relevant new (API) functions included in the most recent release (2021.1) of the *CellNetAnalyzer* toolbox :
------------------------------------------------------------------------------------------------------------

* **CNAgeneMCSEnumerator3**

   Function wrapper for CNAMCSEnumerator3 that allows the computation of gene-MCS with all features introduced in
   CNAgeneMCSEnumerator2 and additionally the option to define optimality constraints.

* **CNAMCSEnumerator3**

   Features all functions of CNAMCSEnumerator2 and additionally the option to define optimality constraints to describe
   desired or undesired flux vectors.

Remarks:
--------

-   As their runtimes are usually very long, scripts MCS_2_smallest.m and MCS_3_any.m compute only one specific setup (by default: ethanol, wGCP). 
    The computation for other products or coupling strengths must be defined by the user. The according lines are marked in the comments
    of the scripts.
    
