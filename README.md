# Unpacking sources of transmission in HIV prevention trials
###
Welcome to the GitHub repository for [*Unpacking sources of transmission in HIV prevention trials with deep-sequence pathogen data – BCPP/ Ya Tsie study*](https://www.medrxiv.org/content/10.1101/2024.08.30.24312845v1).

## Repo organization 

This repo is organized as follows:
- [bcpp_spillovers](bcpp_spillovers):
    - [data](bcpp_spillovers/data/): input data
    - [scripts](bcpp_spillovers/scripts/): scripts to estimate the risk of transmission between communities and the preventable transmissions with a nationwide intervention
    - [results](bcpp_spillovers/results/): preprocessed outputs
 
 
## Dependencies

| Software         | Version  |
|:-----------------|:-------- |
| R                | 4.1.2    |
| Stata (optional) | 13.1     |

Note that Stata is required only if you would like to run the script: estimate_transmission_risk.do

## Installation

Clone the spillover-infections repo and change directory into it

```
git clone git@github.com:magosil86/spillover-infections.git
cd spillover-infections
```

You're all set! Now for some examples

## Examples

####  1. To estimate the risk of transmission between communities we run the script:

*R version*
```
estimate_transmission_risk.r
```

Notes

* This script will fit models to estimate the risk of transmission between communities and save the models in the workspace: `estimate_transmission_risk.rda`
* Pre-run models are available at: `bcpp_spillovers/results/estimate_transmission_risk.rda`
* Tip! To view the pre-run models in R:

```   
      # Load workspace
      load(bcpp_spillovers/results/estimate_transmission_risk.rda)

      # Load libraries
      library(dplyr)         # needed for sorting and merging data.frames
      library(tidyr)         # needed for combining or splitting fields
      library(MASS)          # for predicting fitted values and fitting negative binomial models

      # View post-baseline model
      nb_glm_offset_ln_sample_oc_same_comm_linear_post_baseline

      # View baseline model
      nb_glm_offset_ln_sample_oc_same_comm_linear_baseline_2
      
```

*Stata version*
```
estimate_transmission_risk.do
```

#### 2. To estimate preventable transmissions with a nationwide intervention we run the script:

```
estimate_preventable_trms_with_nationwide_intervention.r
```
Notes
* This script produces a .csv file named: `preventable_trms_nationwide_intervention.csv`
* The .csv file contains three columns:
    * mean estimate of preventable transmissions
    * 95% lower confidence bound
    * 95% upper confidence bound
* Pre-processed output is available at: `bcpp_spillovers/results/preventable_trms_nationwide_intervention.csv`


And that's it, you are all set!


## Getting help
- To ask questions, learn about updates, or just interact with other users, please use the [Github Discussions](https://github.com/magosil86/spillover-infections/discussions) forum.
- To suggest new features or report bugs, please file an [issue](https://github.com/magosil86/spillover-infections/issues).

## Code of conduct
Contributions are welcome. Please observe the [Contributor Code of Conduct](https://github.com/magosil86/spillover-infections/blob/master/CONDUCT.md) when participating in this project.

## Citation
Lerato E. Magosi, Eric Tchetgen Tchetgen, Vlad Novitsky, Molly Pretorius Holme, Janet Moore, 
Pam Bachanas, Refeletswe Lebelonyane, Christophe Fraser, Sikhulile Moyo, Kathleen E. Hurwitz, 
Tendani Gaolathe, Ravi Goyal, Joseph Makhema, Shahin Lockman, Max Essex, Victor De Gruttola, & Marc Lipsitch 
(2024) [Unpacking sources of transmission in HIV prevention trials with deep-sequence pathogen data – BCPP/ Ya Tsie study](https://www.medrxiv.org/content/10.1101/2024.08.30.24312845v1)


## Authors
Lerato E. Magosi, Victor De Gruttola and Marc Lipsitch.

## Maintainer
Lerato E. Magosi magosil86@gmail.com or lmagosi@hsph.harvard.edu

## License

See the [LICENSE](https://github.com/magosil86/spillover-infections/blob/master/LICENSE) file.