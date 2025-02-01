import sys
sh_path = './code/scripts'
sys.path.insert(0, sh_path)

from tests import *
from regression import *
#%matplotlib inline
# Para elegir que modelo usar para regresión, descomentar la línea con el modelo
# que se quiera usar y comentar las otras dos. Despues correr la celda
REGRESS_MODEL = "VAR" #o "VECM" o "OLS"|

df_mensual = pd.read_excel("https://github.com/mbarrena/tesis/blob/main/data/Mensuales%202004%20a%202023.xlsx?raw=true")
df_Arg_mens = df_mensual.copy()
df_Arg = df_Arg_mens.iloc[:240].reset_index(drop=True)[['Año','Mes','ipc_ajust', 'Ebc', 'Efmi','emae', 'impp_usa', 'Psoja_USA', 'Pmaíz_USA', 'Ptrigo_USA', 'tot_04', 'TOTfmi']].copy()
df_Arg = renameColumnOfDataframe(df_Arg, "ipc_ajust", "ipc")
df_Arg = renameColumnOfDataframe(df_Arg, "Efmi", "E")
df_ERPT_Arg = df_Arg #ARGENTINA, mensual 2004 a 2023. Modelo cálculo ERPT de corto plazo con E nominal, PBI real desestacionalizado, precios externos (P de importados de USA) e IPC.
df_ERPT_Arg = makeLogColumns(['ipc', 'Ebc', 'E', 'emae', 'impp_usa', 'Psoja_USA', 'Ptrigo_USA', 'Pmaíz_USA', 'tot_04', 'TOTfmi'],df_ERPT_Arg)
df_ERPT_Arg = makeDiffColumns(['ipc', 'Ebc', 'E', 'emae', 'impp_usa', 'Psoja_USA', 'Ptrigo_USA', 'Pmaíz_USA', 'tot_04', 'TOTfmi'],df_ERPT_Arg)

#ChowTest(df_Arg59,19,21,y_vars=["ipc"],variables=removeInplace(df_Arg59.columns[2:-4],"tot"))
#BaiPerron(df_Arg59, variables=['ipc'])

import regression
regression.REGRESS_MODEL = 'LocalProjections'
regress(["ipc", "Ebc", "impp_usa", "emae"], df_ERPT_Arg, maxlags=2)