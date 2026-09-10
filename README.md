# preparedness-network-design
This project contains the input data, the source code, and the output results for the paper "Emergency Preparedness Network Design: Joint Chance-Constrained Large-Scale Optimization and a Case Study of Typhoon Maria".

# For the input data:

* Fi.txt records the fixed location and operation cost for each relief facility.

* Hi.txt records the pre-positioning cost of each relief package in each relief facility, including the acquisition and holding costs of each relief package.

* Cij.txt records the allocation cost of establishing service linkage between each relief faciity and each affected area.

* Dij.txt records the distance between each relief facility and each affected area.

* Qi.txt records the storage capacity for each relief facility.

* MeanDj.txt records the forecasted demand of each affected area.

# For the source code and output results:

* case study experiment.cpp is the code for case study, with output in folder results/13&14(6hours)

* large-scale experiment.cpp is the code for large-scale instances, with output in folder results/Large scale

* The sensitivity analysis on facility coverage is conducted with source code case study experiment.cpp, with output in folder results/13&14(109876)

* lognormal.cpp is the code for the sensitivity analysis on skewed distributed demand, with output in folder results/lognormal

* extended model - K.cpp and extended model - T_k.cpp is the code for responsiveness differentiation, with output in folder results/k = (1,2,3) and T_K = 1 and folder results/K= 1 and T_k = (1,2,3)
