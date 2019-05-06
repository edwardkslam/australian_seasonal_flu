# Code for manuscript

Reproduces the analyses and figures in the Main Text and Supplementary figures

## Data
### ``raw_data.csv``
This file contains the full dataset of 18,250 laboratory confirmed and antigenically characterised influenza virus infections from Adelaide, Brisbane, Melbourne, Perth and Sydney from 01/01/2000 to 31/12/15.

The column names are as follows:
- ``"city"`` = City
- ``"type"`` = Virus type
- ``"subtype"`` = Virus subtype/lineage
- ``"specimen_date"`` = Collection Date of submitted specimen
- ``"reference_strain"`` = Antigenic variant of submitted specimen, based on comparison against current vaccine reference strain by HI assay
- ``"agree_with_phylogenetic_analyses"`` = Y/N
- ``"assumed_antigenic_variant"`` = Assumed antigenic variant of submitted specimen, given potential errors in HI testing.

As mentioned in the Materials and Methods, since the Southern Hemisphere’s influenza vaccine composition is updated every September (typically after the end of the influenza season in Australia), viruses characterised by HI during the preceding season may have been misidentified.  To account for this potential source of bias, we compared the antigenic characterization data against phylogenetic data. In particular, we identified two potential instances of misidentification of A/H3 viruses in 2004 and 2005.  In these seasons, there is evidence to suggest that the newly emerged antigenic variant had already replaced its predecessor so all A/H3 cases are assumed to be attributable to the new variant.

### ``epi_table.csv``
Here ``assumed_antigenic_variant`` from ``raw_data.csv`` is used to account for potential cases of misidentification.  ``reference_strain`` in ``epi_table.csv`` is thus derived from ``assumed_antigenic_variant`` in ``raw_data.csv``.

For each of antigenic variants recorded in a given season and city, we used the algorithm presented in the Materials and Methods to ascertain whether or not above-baseline levels of epidemic activity could be detected.

``epi_table.csv`` contains the summary statistics for each of the epidemics identified, over the study period from 01/01/2000 to 31/12/15.

- ``city`` = City 
- ``strain_year`` = Concatenated string formed from reference_strain and year 
- ``year`` = Year 
- ``subtype`` = Virus Subtype 
- ``reference_strain`` = ``assumed_antigenic_variant`` of submitted specimen 
- ``epi_alarm`` = Y (above-baseline levels of epidemic activity detected) / N 
- ``start`` = Fortnight in which above-baseline levels of activity is first detected \
- ``epi_counts`` = Total Number of specimens collected over the period of epidemic-levels of activity (as defined by algorithm) 
- ``incidence_per_mil`` = ``epi_counts`` divided by the city-specific Annual Estimated Resident Population 
- ``epi_fractional_counts`` = ``epi_counts`` divided by total ``epi_counts`` across epidemics of all subtypes in that ``city`` and ``year`` 
- ``new_ag_marker`` = 1 (first epidemic of an antigenic variant) / 0 (subsequent epidemics) / NA (antigenic variant emerged before the start of study period) 
- ``first_detection_of_new_ag`` = 1 (first detection of an antigenic variant irrespective of whether or not there was above-baseline levels of activity) / 0 (subsequent detections) / NA (antigenic variant emerged before the start of study period) 
- ``first_n_biggest`` = Y (this antigenic variant caused the epidemic with earliest ``start`` and largest ``incidence_per_mil`` in a particular City and Year) / N 
- ``relative_to_first_n_biggest`` = ``epi_counts`` divided by the largest ``epi_counts`` value for that particular ``year`` and ``city`` 
- ``delay`` = difference between ``start`` and the earliest ``start`` value for that particular ``year`` and ``city`` 
- ``prior_everything_scaled`` = sum of all specimens of other subtypes with ``specimen_dates`` prior to ``start`` 
- ``mean_epi_size`` = mean of "incidence_per_mil" across all epidemics of all subtypes, within that particular "city"

### ``epi_table_no_corrections.csv``
Here, we make no assumptions about the existence of potential mis-identifications: ``reference_strain`` in ``epi_table_no_corrections.csv`` is derived from ``reference_strain`` in ``raw_data.csv``. Colnames are the same as ``epi_table.csv``.

### ``mean_fortnightly_climate_30years.csv``
For each of the five cities, daily mean temperature (°C) and relative humidity (%) values from 1985 to 2015 were retrieved from TuTiempo (https://en.tutiempo.net/).  We then calculated the mean absolute humidity ($g/m^3$) values for each two-week period over the 31 years.  The "historic average" mean climatic values for each of the 26 two-week periods of the year were also calculated.

- ``city`` = City 
- ``year`` = Year 
- ``fortnights_since_start_of_year`` = The fortnight  of the year (ranges from 1 to 26)
- ``mean_AH`` = mean absolute humidity (g/m3) for that particular ``fortnights_since_start_of_year``, ``year`` and ``city`` 
- ``mean_RH`` = mean relative humidity (%) for that particular ``fortnights_since_start_of_year``, ``year`` and ``city`` 
- ``mean_temp`` = mean temperature humidity (oC) for that particular ``fortnights_since_start_of_year``, ``year`` and ``city``
- ``mean_AH_for_that_fortnight_of_year`` = The historic average of ``mean_AH`` for that particular ``fortnights_since_start_of_year`` and ``city`` 
- ``mean_RH_for_that_fortnight_of_year`` = The historic average of ``mean_RH`` for that particular ``fortnights_since_start_of_year`` and ``city`` 
- ``mean_temp_for_that_fortnight_of_year`` = = The historic average of ``mean_temp`` for that particular ``fortnights_since_start_of_year`` and ``city``
- ``d.AH`` = Anomalous absolute humidity: defined as ``mean_AH`` - ``mean_AH_for_that_fortnight_of_year`` 
- ``d.RH`` = Anomalous relative humidity: defined as ``mean_RH`` - ``mean_RH_for_that_fortnight_of_year`` 
- ``d.temp`` = Anomalous temperature: defined as ``mean_temp`` - ``mean_temp_for_that_fortnight_of_year``

### ``robustness vs geoghegan timings/Geoghegan_2018_estimated_A_onsets.csv``
The onset timings of Influenza A epidemic activity across the five cities from 2007 - 2015 were derived from Figure 3b of Geoghegan et al. (2018).

- ``city`` = City 
- ``year`` = Year 
- ``start`` = Fortnight in which above-baseline levels of activity is first detected

### ``robustness vs geoghegan timings/epi_table_with_geoghegan_estimates.csv``
Additional columns were appended to "epi_table.csv": \
type = Virus Type \
earliest_A = Y (this antigenic variant is of Type A and has the epidemic with the earliest "start" amongst Influenza A antigenic variants active within that particular "year" and "city") / N \
largest_A = Y (this antigenic variant is of Type A and has the epidemic with the largest "incidence_per_mil" amongst Influenza A antigenic variants active within that particular "year" and "city") / N \
poor_timeseries = Y (it is difficult to discern between periods of below and above baseline activity, based upon inspection by eye) / N

earliest_geog_start = if "earliest_A" == Y, "earliest_geog_start" takes the value of "start" from "Geoghegan_2018_estimated_A_onsets.csv".  Otherwise, "earliest_geog_start" takes the original, unchanged value of "start" from "epi_table.csv". \
firstNbiggest_earliest_geog = Y (this antigenic variant caused the epidemic with the earliest "earliest_geog_start" and largest "incidence_per_mil" in a particular City and Year) / N

largest_geog_start = if "largest_A" == Y, "largest_geog_start" takes the value of "start" from "Geoghegan_2018_estimated_A_onsets.csv".  Otherwise, "largest_geog_start" takes the original, unchanged value of "start" from "epi_table.csv". \
firstNbiggest_largest_geog = Y (this antigenic variant caused the epidemic with the earliest "largest_geog_start" and largest "incidence_per_mil" in a particular City and Year) / N

poor_timeseries_geog_start = if "poor_timeseries" == Y, "poor_timeseries_geog_start" takes the value of "start" from "Geoghegan_2018_estimated_A_onsets.csv".  Otherwise, "poor_timeseries_geog_start" takes the original, unchanged value of "start" from "epi_table.csv". \
firstNbiggest_poor_geog = Y (this antigenic variant caused the epidemic with the earliest "poor_timeseries_geog_start" and largest "poor_timeseries_geog_start" in a particular City and Year) / N

## Analysis scripts
### main text analysis scripts/Climatic factors and epidemic onset.R
This set of analyses considers whether or not the onset of epidemics are preceded with fluctuations in climatic factors that are:
1) anomalous when compared with against "typical wintertime" fluctuations using the bootstrap sampling method presented by Shaman (2010). The output is a data frame: final_results (Table S1).
2) anomalous when compared against historical averages that typical for that time of the year. The output is two figures: T_plot and AH-plot (Figure 2).

### main text analysis scripts/Effect of antigenic change.R
Here we assess the effect of antigenic change on:
1) the size of epidemics (Figure 3)
2) the onset timing of epidemics (Figure S8)
3) the temporal synchrony of epidemics across all five cities (Figure S9)

### main text analysis scripts/Effect of prior immunity.R
Here we assess the effect of accumulated antigenic-variant specific immunity on:
1) epidemic incidence (Figure 4)
2) probability of successful epidemic initiation (Figure S16 and Table S6)

### main text analysis scripts/Effect of prior immunity.R
Here we assess the effect of heterosubtypic competition on the size of subsequent epidemic:
1)  the relationship between prior epidemic activity by other subtypes within the same season and city and the relative size of an epidemic (Figure 5).
2)  the relationship between delay in onset timing and relative size of an epidemic (Figure 5).

### main text analysis without phylogenetic informed corrections/Effect of antigenic change no corrections.R
Here we assess (whilst making no corrctions for potential mis-identification during antigenic characterisation) the effect of antigenic change on:
1) the size of epidemics (Figure S17)
2) the onset timing of epidemics (Figure S18)
3) the temporal synchrony of epidemics across all five cities (Figure S19)

### main text analysis without phylogenetic informed corrections/Effect of prior immunity no corrections.R
Here we assess (whilst making no corrctions for potential mis-identification during antigenic characterisation) the effect of accumulated antigenic-variant specific immunity on:
1) epidemic incidence (Figure S20)
2) probability of successful epidemic initiation (Figure S21 and Table S9)



### robustness vs geoghegan timings/Robustness of climatic analyses using Geoghegan et al (2018) timings.R
Here we assess the robustness towards potential inaccuracies in our epidemic onset timing estimates, of our analyses on climatic fluctuations preceding epidemic onset. We considered whether or not these preceding climatic fluctuations were: \
i)  anomalous when compared with against "typical wintertime" fluctuations using the bootstrap sampling method presented by Shaman (2010). \
ii) anomalous when compared against historical averages that typical for that time of the year.

We incorporate alternative estimates for the onset timing of Influenza A epidemic activity from Geoghegan et al. (2018), based on the following series of assumptions: \
1) For each of the seasons between 2007 and 2015, we assumed that our timing estimate for the DOMINANT influenza A subtype was incorrect and replaced it with estimates from Geoghegan et al. (2018). \
    i) Table S2; ii) Figure S4
2) For each of the seasons between 2007 and 2015, we assumed that our timing estimate for the EARLIEST influenza A subtype was incorrect and replaced it with estimates from Geoghegan et al. (2018). \
    i) Table S3; ii) Figure S5
3) For seasons between 2007 and 2015 in which the number of cases for the dominant influenza A subtype were small or it was difficult to discern the period of epidemic from background activity, we assumed that our timing estimate was incorrect and replaced it with estimates from Geoghegan et al. (2018). \
    i) Table S4; ii) Figure S6

We also considered whether or not, more generally, Influenza A epidemic activity was preceded by anomalous fluctuations. \
4) Utilising only the timing estimates by Geoghegan et al. (2018), we assessed if more generally, the onset of influenza A epidemic activity in the seasons from 2007 to 2015 was preceded by periods of anomalous climatic conditions. \
    i) Table S5; ii) Figure S7

### robustness vs geoghegan timings/Robustness of ag change analyses using Geoghegan et al (2018) timings.R
Here we assess the robustness towards potential inaccuracies in our epidemic onset timing estimates, of our analyses on the effect of antigenic change on: \
i) the onset timing of epidemics \
ii) the the temporal synchrony of epidemics across all five cities

We incorporate alternative estimates for the onset timing of Influenza A epidemic activity from Geoghegan et al. (2018), based on the following series of assumptions:
1) For each of the seasons between 2007 and 2015, we assumed that our timing estimate for the DOMINANT influenza A subtype was incorrect and replaced it with estimates from Geoghegan et al. (2018). \
    i) Figure S10 ; ii) Figure S11
2) For each of the seasons between 2007 and 2015, we assumed that our timing estimate for the EARLIEST influenza A subtype was incorrect and replaced it with estimates from Geoghegan et al. (2018). \
    i) Figure S12 ; ii) Figure S13
3) For seasons between 2007 and 2015 in which the number of cases for the dominant influenza A subtype were small or it was difficult to discern the period of epidemic from background activity, we assumed that our timing estimate was incorrect and replaced it with estimates from Geoghegan et al. (2018).
    i) Figure S14 ; ii) Figure S15
