import localprojections as lp
import pandas as pd
import sys 

sh_path = '/home/mako/Documents/personal/tesis_mariana/code/scripts'
sys.path.insert(0, sh_path)

from utilities import renameColumnOfDataframe
from variable_preprocessing import makeLogColumns, makeDiffColumns
from regression import run_irf, LPResults
from visualization import plotIrfWithSignifLP

df_trimestral_crudo = pd.read_excel("https://github.com/mbarrena/tesis/blob/main/data/Data%20trimestral%201950%20a%202023%20con%20DUMMIES%20outliers%20(por%20trimestre).xlsx?raw=true")
df_Arg = df_trimestral_crudo[['año', 'trimestre', 'ipc_ajust', 'E', 'Ebc', 'pbird', 'impp_usa', 'Psoja_USA', 'Pmaíz_USA', 'Ptrigo_USA', 'TOTfmi', "D1", "D2"]].copy()
df_Arg = renameColumnOfDataframe(df_Arg, "ipc_ajust", "ipc")
df_Arg59 = df_Arg.copy().iloc[36:].reset_index(drop=True)[['año', 'trimestre', 'ipc', 'E', 'Ebc', 'pbird', 'impp_usa', 'Psoja_USA', 'Pmaíz_USA', 'Ptrigo_USA', 'TOTfmi', "D1", "D2"]].copy()
df_ERPT_Arg59 = df_Arg59.copy()
df_ERPT_Arg59 = makeLogColumns(['ipc','E', 'Ebc', 'pbird', 'impp_usa', 'Psoja_USA', 'Pmaíz_USA', 'Ptrigo_USA', 'TOTfmi'],df_ERPT_Arg59)
df_ERPT_Arg59 = makeDiffColumns(['ipc','E', 'Ebc', 'pbird', 'impp_usa', 'Psoja_USA', 'Pmaíz_USA', 'Ptrigo_USA', 'TOTfmi'],df_ERPT_Arg59)


signifs = [0.05,0.32]
df = df_ERPT_Arg59
endog=["impp_usa", "E", "ipc", "pbird"]
exog=[]
response = endog+exog 
response = [x for x in response if x not in ["impp_usa","Psoja_USA","Pmaiz_USA","Ptrigo_USA"]] # estimate the responses of all variables to shocks from all variables
irf_horizon = 10 # estimate IRFs up to 10 periods ahead
opt_lags = 5 # include 5 lags in the local projections model
opt_cov = 'robust' # HAC standard errors

lp_results = LPResults([])
lp_results.names = endog+exog
lp_results.signifs = signifs

for signif in signifs:
    opt_ci = 1-signif
    print(f"Signif. {opt_ci}")

    if len(exog) == 0:
        irf = lp.TimeSeriesLP(data=df, 
                        Y=endog, 
                        response=response, 
                        horizon=irf_horizon, 
                        lags=opt_lags, 
                        newey_lags=None, 
                        ci_width=opt_ci
                        )
    else:
        irf = lp.TimeSeriesLPX(data=df, 
                        Y=endog, 
                        X=exog, 
                        response=response, 
                        horizon=irf_horizon, 
                        lags=opt_lags, 
                        newey_lags=None, 
                        ci_width=opt_ci
                        )
    lp_results.append(irf)

run_irf(lp_results, lp_results.signifs)

