
# Description of output files

## LE models

For more details see notebook `life_expectancy`

## GDP only models or individual models (no stacking)

#### 1900

- `shift_1900_gdponly.csv`: shifts using data since 1900 with GDP as the only varying coefficient. See notebook for more details.
- `shift_1900_gdponly_year.csv`: shifts using data since 1900 with GPD and year (standardized variable) as the only covariates. See notebooks for more details.
- `shift_1900_gdponly_autocorrelation.csv`: shifts using data since 1900 with GPD as only varying coefficient and autocorrelation adjustment (one lag of LE variable). See notebooks for more details.

- `shift_1900_full_autocorrelation.csv`: shifts using data since 1900 using all covariates available + adjustment autocorrelation See notebook for more details.

- `shift_1900_stacking.csv`: shifts using data since 1900 with stacking. GDP is the only varying coefficient, all the rest are constant. Stacking consists of 4 models:
  - wy ~ 1 + igdp_log + (igdp_log|ctry_year)
  - wy ~ 1 + igdp_log  + iurban_log +  (igdp_log|ctry_year)
  - wy ~ 1 + igdp_log  + ilit_log +  (igdp_log|ctry_year)
  - wy ~ 1 + igdp_log + ilit_log + iurban_log + (igdp_log|ctry_year)

- `shift_1900_stacking_year.csv`: shifts using data since 1900 with stacking adjusting for year (standardized variable). GDP is the only varying coefficient, all the rest are constant. Stacking consists of 4 models:
  - wy ~ 1 + igdp_log + zyear + (igdp_log|ctry_year)
  - wy ~ 1 + igdp_log  + iurban_log +  zyear + (igdp_log|ctry_year)
  - wy ~ 1 + igdp_log  + ilit_log +  zyear + (igdp_log|ctry_year)
  - wy ~ 1 + igdp_log + ilit_log + iurban_log + zyear + (igdp_log|ctry_year)

#### 1950

- `shift_1950_gdponly.csv`: shifts using data since 1950 with GDP as the only varying coefficient. See notebook for more details.
- `shift_1950_gdponly_year.csv`: shifts using data since 1950 with GPD and year (standardized variable) as the only covariates. See notebooks for more details.
- `shift_1950_gdponly_autocorrelation.csv`: shifts using data since 1950 with GPD as only varying coefficient and autocorrelation adjustment (one lag of LE variable). See notebooks for more details.

- `shift_1950_full_autocorrelation.csv`: shifts using data since 1950 using all covariates available + adjustment outcorrelation. See notebook for more details.See notebook for more details.

- `shift_1950_stacking.csv`: shifts using data since 1950 with stacking. GDP is the only varying coefficient, all the rest are constant. Stacking consists of 4 models:
  - wy ~ 1 + igdp_log + (igdp_log|ctry_year)
  - wy ~ 1 + igdp_log + ilit_log + iurban_log + (igdp_log|ctry_year)
  - wy ~ 1 + igdp_log + iwater_log + isewage_log + ielec_log + (igdp_log|ctry_year)
  - wy ~ 1 + igdp_log + ilit_log + iurban_log + iwater_log + isewage_log + ielec_log + (igdp_log|ctry_year)


- `shift_1950_stacking_year.csv`: shifts using data since 1950 with stacking adjusting for year (standardized variable). GDP is the only varying coefficient, all the rest are constant. Stacking consists of 4 models:
  - wy ~ 1 + igdp_log + zyear + (igdp_log|ctry_year)
  - wy ~ 1 + igdp_log + ilit_log + iurban_log + zyear + (igdp_log|ctry_year)
  - wy ~ 1 + igdp_log + iwater_log + isewage_log + ielec_log + zyear  +(igdp_log|ctry_year)
  - wy ~ 1 + igdp_log + ilit_log + iurban_log + iwater_log + isewage_log + ielec_log + zyear + (igdp_log|ctry_year)


## Cause of death models

For more details see notebook `cause_of_death`.

For each cause we get:

- `shift_cause_X_full.csv`:
  - cause ~ 1 + igdp_log + ilit_log + iurban_log + iwater_log + isewage_log + ielec_log + (igdp_log|ctry_year)
- `shift_cause_X_full_leinitial_year.csv`:
  - cause ~ 1 + igdp_log + ilit_log + iurban_log + iwater_log + isewage_log + ielec_log + LE_initial + zyear + (igdp_log|ctry_year)



