# co-expression-probability

These MatLab scripts calculate the probability of a given neuron co-expressing two neurotransmitters. The script calculates co-expression probabilities for a 
population of 360 neurons that individually express GABA, Allatotropin (Mas-AT), Myoinhibitory peptide (MIP), Tachykinin (TKK), FMRFamide.  
Code can be altered to include the total # of neurons in a given neural population and the avg # (and stdev) of neurons that express each transmitter. 
The random_probability_script assumes no expression dependencies, and thus the predicitons of co-expression are based purely on chance co-expression based 
on independent expression probability. 

Similarly, the script_Predict_FMRF_ATR script allows you to explicity set the co-expression probability of two given transmitters to its observed value. 
For example, you may know that 30% of FMRF expressing neurons also express the transmitter ATR based on physical data from immunocytochemistry. This relationship is then set as an explicit rule, 
leaving the remaining co-expression relationships to emerge based on independent probablility of expression. This allows you to determine if there are key co-expression relationships
in your population of neurons that may be predictive of other relationships in the population in an unbiased manner. 

averaging_all is a simple script used to calculate the average # of neurons that co-express a given transmitter pair from the output of the above scripts. 
These scripts are fully described in Lizbinski et al. 2017 in BioRXIV (https://doi.org/10.1101/167403)
Please cite this paper if you use these scripts.

This project is licensed under the terms of the MIT license.
