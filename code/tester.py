import sys
sh_path = './code/scripts'
sys.path.insert(0, sh_path)

from tests import *
from regression import *
#%matplotlib inline
# Para elegir que modelo usar para regresión, descomentar la línea con el modelo
# que se quiera usar y comentar las otras dos. Despues correr la celda
REGRESS_MODEL = "VAR" #o "VECM" o "OLS"|

df_trimestral_crudo = pd.read_excel("https://github.com/mbarrena/tesis/blob/main/data/Data%20trimestral%201950%20a%202023%20con%20DUMMIES%20outliers%20(por%20trimestre).xlsx?raw=true")
df_Arg = df_trimestral_crudo.iloc[36:].reset_index(drop=True)[['año', 'trimestre', 'ipc_ajust', 'E', 'Ebc', 'pbird', 'impp_usa', 'Psoja_USA', 'Pmaíz_USA', 'Ptrigo_USA', 'TOTfmi']].copy()
df_Arg = df_trimestral_crudo.iloc[36:].reset_index(drop=True)[['año', 'trimestre', 'ipc_ajust', 'E', 'Ebc', 'pbird', 'impp_usa', 'Psoja_USA', 'Pmaíz_USA', 'Ptrigo_USA', 'TOTfmi']].copy()
df_Arg = renameColumnOfDataframe(df_Arg, "ipc_ajust", "ipc")
df_ERPT_Arg = df_Arg.copy()
df_ERPT_Arg = makeLogColumns(['ipc','E', 'Ebc', 'pbird', 'impp_usa', 'Psoja_USA', 'Pmaíz_USA', 'Ptrigo_USA', 'TOTfmi'],df_ERPT_Arg)
df_ERPT_Arg = makeDiffColumns(['ipc','E', 'Ebc', 'pbird', 'impp_usa', 'Psoja_USA', 'Pmaíz_USA', 'Ptrigo_USA', 'TOTfmi'],df_ERPT_Arg)

#ChowTest(df_Arg59,19,21,y_vars=["ipc"],variables=removeInplace(df_Arg59.columns[2:-4],"tot"))
#BaiPerron(df_Arg59, variables=['ipc'])

import regression
regression.REGRESS_MODEL = 'LocalProjections'
regress(endog=["ipc","E","pbird"], data=df_ERPT_Arg, exog=["impp_usa"], maxlags=2, run_other_tests_on=False, lp_threshold=121)