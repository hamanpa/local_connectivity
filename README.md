This repository contains scripts to analyze the data on local connectivity (lateral distance from soma < 500 um>) in cat primary visual cortex from

Armen Stepanyants, Judith A. Hirsch, Luis M. Martinez, Zoltán F. Kisvárday, Alex S. Ferecskó, Dmitri B. Chklovskii, Local Potential Connectivity in Cat Primary Visual Cortex, Cerebral Cortex, Volume 18, Issue 1, January 2008, Pages 13–28, https://doi.org/10.1093/cercor/bhm027


### Data
Data were privately shared. It should be added into subfolder `data`
```'./data/Cat_Np.mat'```

Raw data are a dictionary with following items:


| key | value | shape | comments |
| --- | ---   | ---   | ---      |
| Np_int_ee | numpy array | (150, 150, 51) | Expected number of potential synapses Exc -> Exc |
| Np_int_ei | numpy array | (150, 150, 51) | Expected number of potential synapses Exc -> Inh |
| Np_int_ie | numpy array | (150, 150, 51) | Expected number of potential synapses Inh -> Exc |
| Np_int_ii | numpy array | (150, 150, 51) | Expected number of potential synapses Inh -> Inh |
| | | | |
| P_int_ee | numpy array | (150, 150, 51) | Probability of potential connectivity Exc -> Exc |
| P_int_ei | numpy array | (150, 150, 51) | Probability of potential connectivity Exc -> Inh |
| P_int_ie | numpy array | (150, 150, 51) | Probability of potential connectivity Inh -> Exc |
| P_int_ii | numpy array | (150, 150, 51) | Probability of potential connectivity Inh -> Inh |
| | | | |
| roi | numpy array | (1, 5) | |
| roe | numpy array | (1, 5) | |
| htz | numpy array | (1, 5) | borders of layers [um] |
| z_int | numpy array | (1, 150) | depth, 5 -- 1495 um |
| Rxy_int | numpy array | (1, 51) | lateral position, 0 -- 500 um |

For `Np_int_XY` and `P_int_XZ` axis are indexing 
    somatic depth of presynaptic cell (`z_pre`), 
    somatic depth of postsynaptic cell (`z_post`),
    relative lateral (horizontal) position (`rho`).

#### Comments
data are sufficient/insufficient for

Missing data
- `Np_int_XY` has the same missing data as `P_int_XY`
- for a given combination of `z_pre` and `z_post`there are either data available for all values of `rho` or no values at all.
- thus one can check the missing data with a fixed value of `rho`
- Figs 8, 10, 11, 12 show there are data missing

Projection percentage in the data
- 




### Describe file and how to use it 
Layer 1 is omitted from the analysis, since no data for Layer 1 were provided


### Describe fitting process (how it works)

It can fit following functions:


Jak resim normalizaci???
