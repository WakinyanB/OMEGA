# Readme

Data & Codes for:<br>
**The interplay between migration and selection on the dynamics of pathogen variants**

<!--Wakinyan Benhamou<sup>1,2</sup>, Rémi Choquet<sup>1</sup> and Sylvain Gandon<sup>1</sup>

*<sup>1</sup> CEFE, Univ Montpellier, CNRS, EPHE, IRD, Montpellier, France*<br>
*<sup>2</sup> High Meadows Environmental Institute, Princeton University, Princeton, NJ, USA*<br>-->

## Data (sources and short descriptions)

All data used in the scripts are in the 'Data' folder. 

### VOC_202012_01_Technical_Briefing_5_Data_England.ods

[DOWNLOAD](https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/957631/Variant_of_Concern_VOC_202012_01_Technical_Briefing_5_Data_England.ods) (ods file)

&#11169;&emsp; Fequencies of S-Gene Target Failures (SGTF) in the 9 regions of England; underlying data from *Public Health England* [*Technical Briefing 5: Investigation of novel SARS-CoV-2 variant - Variant of Concern 202012/01*](https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/959426/Variant_of_Concern_VOC_202012_01_Technical_Briefing_5.pdf).

<sub>All PHE (now  replaced by UKHSA) technical briefings dealing with the investigation of SARS-CoV-2 variants may be found on *GOV.UK*:</sub><br>
<sub>Technical briefings 1 to 23: https://www.gov.uk/government/publications/investigation-of-novel-sars-cov-2-variant-variant-of-concern-20201201</sub><br>
<sub>From technical briefing 24: https://www.gov.uk/government/publications/investigation-of-sars-cov-2-variants-technical-briefings</sub>

### Delta_data.rds

&#11169;&emsp; Data for the sweep of the Delta variant. The data were shared by Erik Volz who used them in [Volz, E. Fitness, growth and transmissibility of SARS-CoV-2 genetic variants. Nat Rev Genet 24, 724–734 (2023)](https://doi.org/10.1038/s41576-023-00610-z) (see Fig. 1-a).

### changes-visitors-covid

&#11169;&emsp; Google mobility reports. This dataset was downloaded from the website [*Our World in Data*](https://ourworldindata.org/covid-mobility-trends) on April 8, 2025.<br>

<sub>*Hannah Ritchie (2020) - “Google Mobility Trends: How has the pandemic changed the movement of people around the world?” Published online at OurWorldinData.org. Retrieved from: ['https://ourworldindata.org/covid-mobility-trends'](https://ourworldindata.org/covid-mobility-trends) [Online Resource]*</sub>

## R codes

All the scripts are in the root directory of this repositery.

Scripts:
- **Functions.R** (script with custom functions, including the ODE models)
- **Theoretical_figures.Rmd** (plot all the theoretical figures)
- **Analyses_VOC_Alpha_Delta.R** (analyses and plots using the data for the Alpha and Delta variants in England)

## Outputs

- **long_term_diff.RDS** (first panel, N<sup>B</sup>/N<sup>A</sup>=1)
- **long_term_diff2.RDS** (second panel, N<sup>B</sup>/N<sup>A</sup>=10)

This two files are generated in **Theoretical_figures.Rmd** when we compute the long-term differentiation. The code takes some time to run so results are saved here.
