# Can Species Distribution Model predictions correlate with local abundances of non-native species?
The project aimed to assess whether relative local abundances for non-native species correlated with habitat suitability projections from Species Distribution Models [SDMs].
Global occurrence data was used for two non-native species [_Rattus rattus_ (Linnaeus, 1758) and _Neogale vison_ (Schreber, 1777)], alongside environmental predictors, 
to produce habitat suitability projections. Four models were produced, two for each species. For both species, there was a Climate Only model that used on predictors 
from the WorldClim set (Fick and Hijmans, 2017) and a Footprint model that also included a human footprint predictor (Venter et al, 2016; Venter et al., 2018).
The results of these models was tested for correlation with the relative local abundance of each non-native species.

## Repository
The Code folder contains the annotated scripts used for data processing, model building and model evaluation.
The Data folder contains the raw global occurrence data from GBIF (GBIF, 2026) and the local abundance data from published literature, for each species.
In the Evaluations folder, is the model output and evaluations for each of the four models. 
Finally, the Outputs folder contains the results of the correlation tests for each model. As well as the habitat suitability predictions used for the tests.

## FAIR Data Principles
**Findable**: The data and associated scripts are locatable through this GitHub repository, as referenced in the dissertation.

**Accessible**: All datasets are in either of the commonly used .csv or .xlsx formats. Scripts are formatted for use in R, a widely available software. 
Whilst the model evaluations and outputs in the Evaluations folder are in an .rds format that can be accessed through the R software.

**Interoperable**: All naming conventions remain consistent throughout the project and datasets. 
The accompanying R scripts are annotated and describe the processing and analysation of the data.

**Reuse**: All analytical processes are annotated within the scripts. Alongside the raw data, this enables the repetition of any stage and 
for easy incorporation into future research.

## Contact
For any queries about the data in this repository, or its associated results, please contact Malcolm Stewart at s2193313@ed.ac.uk

#### References
Espinoza, F., 2026 Personal communication: Email to Malcolm Stewart, 24 July.

Fick, S.E. and Hijmans, R.J., 2017. WorldClim 2: new 1km spatial resolution climate surfaces for global land areas. International Journal of Climatology 37 (12): 4302-4315.

GBIF.org. 2026. GBIF Home Page. Available online at: [https://www.gbif.org/](https://www.gbif.org/). Accessed: 13/07/2026.

Venter, O., Sanderson, E.W., Magrach, A., Allan, J.R., Beher, J., Jones, K.R., Possingham, H.P., Laurance, W.F., Wood, P., Fekete, B.M. and Levy, M.A., 2016. Global terrestrial Human Footprint maps for 1993 and 2009. Scientific data, 3(1), p.160067.

Venter, O., Sanderson, E.W., Magrach, A., Allan, J.R., Beher, J., Jones, K.R., Possingham, H.P., Laurance, W.F., Wood, P., Fekete, B.M., Levy, M.A., and Watson, J.E., 2018. Last of the Wild Project, Version 3 (LWP-3): 2009 Human Footprint, 2018 Release. Palisades, New York: NASA Socioeconomic Data and Applications Center (SEDAC). Available online at: [https://doi.org/10.7927/H46T0JQ4]. Accessed: 24/06/2026.
