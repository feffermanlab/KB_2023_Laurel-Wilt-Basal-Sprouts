# KB_2023_Vector-Borne-Tree-Disease
This is the READ ME file for all code related to 
Mathematical Model of Basal Sprout Production in Vector-Borne Tree Disease 
By Kelly Buch and Nina Fefferman 
Published in Forests in February 2023.

Programing Language: Matlab 2022a
There are two pertinent files of code: chapter1_finalcode.m and stabilityanalusis_condensed. They are described below.
chapter1_finalcode.m : This file produces the figures and results presented in the main text of the paper. The file begins by defining our parameter sweep, running the short-term ODE over all parameter combos, and then analyzing changes in trends through the parameter sweep as well as make a short- term plot. It then considers the long term effects of rho and sigma.
Supporting files: Laurel_Model1_Equations.m (defines the differential equations as needed by ode45) and params_refigured_simplified (defines the base parameter values)
Necessary toolboxes: Statistics and Machine Learning (v12.3)
stabilityanalysis_condensed: This provides support for all results within the appendix, specifically related to stability of the equilibrium. It begins by defining the equations and initializing variables. In section 2, it considers the stability of the disease free equilibrium algebraically. In section 3, the endemic equilibrium is considered algebraically. Please note that you must run section 1, then section 2, then section 1, and then section 3. That is, the set up must be reran before section 3. Section 4 defines a parameter sweep on which to consider equilibrium. Section 5 calculates the equilibrium value according to the trajectories created by ode45 (5a), calculates the theoretical equilibrium according to the algebraic expressions in section 3 (5b), corrects a common issue manually (5c),  and computes an executives summary of stability analysis by comparing the theoretical and numerical equilibrium and by determining the sign of eigenvalues. If the Stability report determines instability somewhere according to the ode45 (inf) calculations, there is probably more correcting to do in section 5c.
Supporting Files: None
Required Toolboxes: Parallel Computing Toolbox (v7.6)  and Symbolic Math Toolbox (v9.1) 
