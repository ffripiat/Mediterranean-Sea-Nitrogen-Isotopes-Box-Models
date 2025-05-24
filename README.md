Prognostic multi-box (4) model to simulate nitrogen isotope dynamics in the Mediterranean Sea

The model is run on MATLAB.

Reference: 


Wald, T., F. Fripiat, A.D. Foreman, Y. Ryu, D. Marconi, T. Tanhua, G. Sisma-Ventura, D.M. Sigman, G.H. Haug, and A. Martínez-García (2025). Origins of the nitrate 15N depletion in the Mediterranean Sea. Accepted in Global Biogeochemical Cycles.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FILE DESCRIPTION for the four-box model:

med_main_template_fourbox_steadystate.m : master file for the sensitivity experiments for the four-box model, assuming that we reach a steady state (1,000 year)
med_main_template_fourbox_transient.m : master file for the sensitivity experiments for the four-box model, assuming that transient simulation (70 years)
med_fourbox_mixingWE_ode.m: ODE solver for both the "med_main_template_fourbox_steadystate.m" and "med_main_template_fourbox_transient.m"

HOW TO for the four-box model:

to run "sensitivity_generation.m" to generate the "sensitivity" matrix by randomly varying parameters over a range well beyond literature estimates (see table 1 in the reference manuscript)
to run either "med_main_template_fourbox_steadystate.m" or "med_main_template_fourbox_transient.m"
