This directory includes a demonstration of the fitness function 
and optimization method described in

    Weaver, CM and Wearne, SL.  The role of action potential shape
    and parameter constraints in optimization of compartment models.
    Neurocomputing 69: 1053-1057, (2006).



I.  To run this demo:

% nrngui main_demoMRF.hoc

A.  Different Objective Functions

This demo evaluates the difference betweeen model and target data by two different functions, specified under the Multiple Run Fitter Generators:

(1) RegionFitness, i.e. mean squared error between voltage traces
(2) APShpFRCVFitness:  a linear combination of the errors in time-aligned AP shape, mean FR, and irregularity measured by the CV.  (See Weaver and Wearne reference above).  The weights of each component can be changed by displaying the APShpFRCVFitness generator, then choosing Select, followed by the appropriate drop-down menu choice.

Click on one of the buttons to the left.  This will initialize the model, run it, and evaluate the error between that model and the target data shown in the MRF windows.

Summary data are written to matlab files (.m extension) that print the model and target voltage traces, plus a summary of the AP Shape / FR / CV error.  The name of this output file can be changed by bringing up the appropriate MRF Generator and choosing Select -> Set output info.

B.  Simulated annealing for parameter optimization

   If you go to the MRF window and choose Parameters -> Select Optimizer, there are two new options:  

    - Simulated Annealing
    - Constrained Simulated Annealing with Recentering

    The first, implemented by Andrew Davison, performs the Simplex-based simulated annealing optimization algorithm found in Numerical Recipes in C.

    The second, implemented by Christina Weaver, performs a version of this algorithm that incorporates constraints on the parameters to be optimized.  See Weaver and Wearne (2006) for details, and for the original reference (Cardoso et al. 1996) which introduced this algorithm.  See the enclosed PDF help file for details of this NEURON implementation.


   In addition to a new set of Optimizer menu items specific to these optimization techniques, one change has been made to the general MRF Optimize menu.  In the past the output from a parameter search was always output to the same filename, so that one would often run the risk of overwriting previous searches.  The name of the output file can be changed by clicking on "Change output filename" and entering a new string, e.g. "filename".  Data will be written to "filename".tmp during the search, and "filename".fit once the search has finished.




II.  Files in the model/ directory:

biophys_params.hoc
error_testshapes.hoc
instrumentation.hoc
main_demoMRF.hoc
modelAII.hoc
morph.hoc
numerics.hoc
testshapesGUI.hoc

cad.mod
cahi.mod
fn.mod
ka.mod
kca.mod
nap.mod

MRFdemos.ses
MRFdemos.ses.ft1
MRFdemos.ses.fd1
Vrun.ses

sample .m output files:

eWin105const_dflt.m
eWin105const_burst.m

III. Files in the optmz/ directory

[TOM & MICHAEL -- These are the files for Simulated Annealing / APShpFRCVFitness that need to be integrated appropriately into NEURON.]

mulfit.hoc:  goes just above mulfit/ directory
feature_weaver.mod:  to incorporate new Vector functions

mulfit/:
     mulfit1.hoc
     fitparm.hoc
     eparmlst.hoc
    simanneal_seq_weaver_Feb07.hoc  (Unconstrained SA)
     simanneal_cardoso.hoc          (Constrained SA + Recenter)
     eonerun.hoc
     e_apwinfrcv.hoc




IV.  Disclaimer

I think this is everything a user would need to reproduce my functions.  I removed a lot of comments & debugging stuff from the e_apwinfrcv.hoc file - I hope I didn't break something in the meantime!

Other files in this directory:

weaver06_Optimization.pdf:  The Weaver & Wearne reference cited above.
simanneal_unconst.html:  HTML help file written by Andrew Davison (Unconstrained SA)
simannrctr_help.pdf:  Constrained Simulated Annealing help file.