import localprojections as lp
import pandas as pd
from scripts.utilities import renameColumnOfDataframe
from scripts.variable_preprocessing import makeLogColumns, makeDiffColumns

df_trimestral_crudo = pd.read_excel("https://github.com/mbarrena/tesis/blob/main/data/Data%20trimestral%201950%20a%202023%20con%20DUMMIES%20outliers%20(por%20trimestre).xlsx?raw=true")
df_Arg = df_trimestral_crudo[['año', 'trimestre', 'ipc_ajust', 'E', 'Ebc', 'pbird', 'impp_usa', 'Psoja_USA', 'Pmaíz_USA', 'Ptrigo_USA', 'TOTfmi', "D1", "D2"]].copy()
df_Arg = renameColumnOfDataframe(df_Arg, "ipc_ajust", "ipc")
df_Arg59 = df_Arg.copy().iloc[36:].reset_index(drop=True)[['año', 'trimestre', 'ipc', 'E', 'Ebc', 'pbird', 'impp_usa', 'Psoja_USA', 'Pmaíz_USA', 'Ptrigo_USA', 'TOTfmi', "D1", "D2"]].copy()
df_ERPT_Arg59 = df_Arg59.copy()
df_ERPT_Arg59 = makeLogColumns(['ipc','E', 'Ebc', 'pbird', 'impp_usa', 'Psoja_USA', 'Pmaíz_USA', 'Ptrigo_USA', 'TOTfmi'],df_ERPT_Arg59)
df_ERPT_Arg59 = makeDiffColumns(['ipc','E', 'Ebc', 'pbird', 'impp_usa', 'Psoja_USA', 'Pmaíz_USA', 'Ptrigo_USA', 'TOTfmi'],df_ERPT_Arg59)


df = df_ERPT_Arg59
endog=["impp_usa", "E", "ipc", "pbird"]
exog=[]
response = endog+exog 
response = [x for x in response if x not in ["impp_usa","Psoja_USA","Pmaiz_USA","Ptrigo_USA"]] # estimate the responses of all variables to shocks from all variables
irf_horizon = 10 # estimate IRFs up to 10 periods ahead
opt_lags = 5 # include 5 lags in the local projections model
opt_cov = 'robust' # HAC standard errors
opt_ci = 0.95 # 95% confidence intervals

if len(exog) == 0:
    irf = lp.TimeSeriesLP(data=df, 
                    Y=endog, 
                    response=response, 
                    horizon=irf_horizon, 
                    lags=opt_lags, 
                    newey_lags=None, 
                    ci_width=0.95
                    )
else:
    irf = lp.TimeSeriesLPX(data=df, 
                    Y=endog, 
                    X=exog, 
                    response=response, 
                    horizon=irf_horizon, 
                    lags=opt_lags, 
                    newey_lags=None, 
                    ci_width=0.95
                    )

irfplot = lp.IRFPlot(irf=irf, # take output from the estimated model
                     response=['ipc'], # plot only response of invest ...
                     shock=["E"], # ... to shocks from all variables
                     n_columns=2, # max 2 columns in the figure
                     n_rows=2, # max 2 rows in the figure
                     maintitle=f'LP: IRFs E -> ipc', # self-defined title of the IRF plot
                     show_fig=True, # display figure (from plotly)
                     save_pic=False # don't save any figures on local drive
                     )

