## Empirical Moments for "Self-Fulfilling Debt Crises Revisited"

This sub-directory contains data and code for the calibration of the model.

Several data series are from Haver Analytics, and therefore cannot be posted publicly.  Those wishing to replicate the estimation can either download the data directly from Haver using the codes provided below, or can contact the authors.

### The Process for GDP

The stata do files, state_space_estimation_mexico.do and state_space_estimation_italy.do, contain the code to estimate the process for GDP for Mexico and Italy, respectively.  The model is:
$$\Delta y_t= g_t+z_t-z_{t-1} \\
g_t=\mu_g+\rho_g g_{t-1}+\sigma_g \epsilon_t \\
z_t = \sigma_z \nu_t$$
where $\epsilon,\nu$ are standard Normal and $y$ is log real GDP.

The estimation for Mexico uses constant USD GDP from 1980Q through 2020Q3 (series code S273NGCD@EMERGE) and is estimated with state_space_estimation_mexico.do. 

For Mexico, the results are:

| Parameter      | Estimate     | Standard Error     |
| :------------- | :----------: | :-----------: |
|  $\rho_g$ | .3308571  | (.2577268)     |
| $\sigma_g$   | .0138298 | (.0039443)  |
| $\mu_g$ | .0033107 | (.0017096) |
| $\sigma_z$ |  .0113592 | (.002094) |

Estimation of same model (using state_space_estimation_italy.do) on Italian real GDP 1960Q1-2020Q3 from the OECD contained in the spreadsheet italy_gdp.xlsx:

| Parameter      | Estimate     | Standard Error     |
| :------------- | :----------: | :-----------: |
|  $\rho_g$ |.579109   | (.1677794)     |
| $\sigma_g$   |  .0082939  | (.0022173)   |
| $\mu_g$ | .0022693  | (.0010752) |
| $\sigma_z$ | .0096982  | (.0010449) |


#### Debt and Revenue


The calibration involves targeting the ratio of debt payments to tax revenue.  Using Haver Analytics codes N273FGR@EMERGE[A] for nominal revenue and S273NGDP@@EMERGE for nominal GDP, we compute that tax revenue is 0.176 of GDP in 1995Q1.  In Table 1 of Cole-Kehoe JIE 1996, debt due in 1995Q1 is 3,015 Cetes, 9,873.94 Tesobonos, or total 12,888.94, all in millions USD.  Haver reports a (current) USD quarterly GDP (Haver code H273NGPD@EMERGE) for Mexico of 84.247 billion.  Restricting attention to Tesobonos, debt due to GDP is  0.117.  Dividing this by the tax-to-gdp ratio of 0.176, we obtain that debt due in 1995Q1 is 66% of tax revenue.  The face value of debt relative to annualized revenue is obtained by dividing the face value of debt (N273FDG@EMERGE) by 4* quarterly revenue (N273FGR@EMERGE[A]).  This ratio is 1.78 for 1995Q1, which is the 97.5th percentile of the 1990Q1-2020Q3 sample.  Hence, we target a debt due to tax revenue ratio of 0.66 at the 97.5th percentile of the simulated ergodic distribution.  

For Italy, tax data are from the Bocola-Dovis AER data set as well as italy_tax.xls, which was downloaded from OECDstat on 1/16/21.  For 2011, total tax revenue was 686,500 million euros.  Social security contributions in that year was 211,637.  Tax revenue net of SSC was 474,863 million. GDP in 2011Q4 at an annualized rate was 1,635,142.  T/Y was therefore 42.0% and (T-SSC)/Y= 29.0%. Face value of total debt was 1,800,287, of which Bocola-Dovis report 1,588,796.5 was in securities.  D/T was 236.1% for securities, where T is an annual number.  D/(T-SSC) was 334.6%. Bocola-Dovis report debt due in the year following 2011Q4 of 373,000 million euros. Debt due to GDP is therefore  22.8% of GDP and 78.5% of tax revenue net of SSC.  2011Q4 represented a debt burden relative to tax revenue that is the 93rd percentile of the Bocola-Dovis dataset. We therefore target debt payments of 78.5% of resources at the 93rd percentile of the model’s ergodic distribution of debt. 

### Spreads

The standard deviation of the Mexican EMBI sovereign spread, for the available sample period of 1994Q1-2020Q3, is 2.3%.  This is computed using the EMBI+ spread is from Haver Analytics (code P273S@EMBI[A]).  

For Italy, the spread is Italian government bond yields, 5-year  (average), Haver code N136RG5@G10[A], minus the interest rate on German federal government securities of 5-6 year maturity (N134RG5@G10[A]).  We target the standard deviation of the spread for the period 1999Q1-2020Q4.
