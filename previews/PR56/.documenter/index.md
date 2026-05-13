


# TransmissionLinelist.jl {#TransmissionLinelist.jl}

A Julia + Turing model fitted to the Epuyén 2018–19 Andes hantavirus outbreak ([Martínez et al. 2020, NEJM](https://doi.org/10.1056/NEJMoa2009040)).

Four quantities are estimated jointly from the line list: the incubation period, the transmission timing of each secondary infection relative to its source&#39;s symptom onset, a weekly time-varying reproduction number, and offspring dispersion. Exposure and onset dates are interval-censored. Each case is given a continuous latent infection time and a continuous latent onset time, each sampled within its recorded window. Generation interval and serial interval are derived from the fitted distributions in post-processing.

## Pages {#Pages}
- [Model](model.md) — priors, data augmentation, GI / SI derivation.
  
- [Limitations](limitations.md) — known caveats around exposure encoding, late R(t) bins, dispersion identifiability, and right-truncation.
  
- [Analysis walkthrough](analysis.md) — runs the full analysis end to end.
  
- [API Reference](api.md) — exported functions.
  

## Citing {#Citing}
> 
> Martínez VP, Di Paola N, Alonso DO, et al. _&quot;Super-spreaders&quot; and person-to-person transmission of Andes virus in Argentina._ N Engl J Med 2020;383:2230–41. [doi:10.1056/NEJMoa2009040](https://doi.org/10.1056/NEJMoa2009040)
> 


The reporting follows:
> 
> Charniga K, et al. _Best practices for estimating and reporting epidemiological delay distributions of infectious diseases._ 2024. [arXiv:2405.08841](https://arxiv.org/abs/2405.08841)
> 

