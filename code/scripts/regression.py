from statsmodels.tsa.api import VAR
from statsmodels.tsa.vector_ar.vecm import VECM
import pandas as pd
import statsmodels.api as sm
from statsmodels.tsa.stattools import adfuller

from variable_preprocessing import *
from utilities import *
from visualization import *

# Para elegir que modelo usar para regresión, descomentar la línea con el modelo
# que se quiera usar y comentar las otras dos. Despues correr la celda
#REGRESS_MODEL = "OLS"
REGRESS_MODEL = "VAR"
#REGRESS_MODEL = "VECM"


def regress(endog, data, exog=[], maxlags=3, rModel=None, estacionalidad=True):
  """
  Funcion que realiza la regresion OLS
  - `endog`: la lista de variables endogenas
  - `data`: es la base de datos
  - `exog`: la lista de variables exogenas
  - `maxlags`: máximo número de lags para la selección de orden al hacer modelo.fit()
  - `rModel`: parámetro opcional, permite cambiar el modelo de regresión (se usa
  el valor de REGRESS_MODEL elegido al principio de la notebook por defecto)
  """
  if rModel is None:
    rModel = REGRESS_MODEL
  datos=data[exog+endog].copy().reset_index(drop=True)
  # el metodo dropna() me permite eliminar las filas que tienen algun valor missing
  datos=datos.dropna().astype(float)
  if rModel == "OLS":
    if len(exog) > 1 or len(endog) == 0:
      raise Exception(f'''OLS requiere a lo sumo una variable exógena (dependiente)
      y que el resto sean endógenas (independientes).
      Se están pasando {len(endog)} endógenas y {len(exog)} exógenas''')
    Y_endog=datos[endog].astype(float)
    X_exog=datos[exog] # MA; 25 de febrero 2024, aca podes agregar los nombres de las dmmies que quieras. Asi:  datos[exog] + ["D1","D2", "D3"]
    X_exog = X_exog.astype(float)
  elif len(endog) >=2:
    Y_endog=datos[endog].astype(float)
    if len(exog)==0:
      X_exog=None
    else:
      X_exog= datos[exog]
      X_exog = X_exog.astype(float)
  else:
    raise Exception(f"VAR y VEC requieren al menos dos variables endógenas, pero sólo se está pasando {len(endog)}.")

  print(f"Realizando regresión con modelo {rModel}")
  if rModel == "OLS":
    modelo=sm.OLS(Y_endog,X_exog)
    modelo.select_order(8)
    resultados=modelo.fit(maxlag=8, ic='aic')
    print(str(resultados.summary()))

  elif rModel == "VAR":
    #27 de febrero- copiar desde aca
    if estacionalidad:
      dates, freq = generatePeriodIndex(data) #primera original
      dummies = generatePeriodDummies(data, index=dates)
      data_exog = pd.concat([X_exog, dummies], axis=1, ignore_index=True)
      #MA NUEVO 18/06
      if X_exog is not None:
        data_exog.columns = np.concatenate((X_exog.columns.values, dummies.columns.values))
      else:
        data_exog.columns = dummies.columns.values
      #MA NUEVO 18/06 ES HASTA ACA
    else:
      data_exog = X_exog
      dates=None
      freq=None
    print(data_exog)
    modelo=VAR(Y_endog,data_exog,dates=dates,freq=freq) #ultima original
    # 27 de ferero - HASTA LA ANTERIOR A ACA

    lor = modelo.select_order(6)

    print(f"Selected orders are: {lor.selected_orders}")
    display(lor.summary()) #Imprime tabla VAR Order Selection

    resultados=modelo.fit(maxlags)  # Acá cambiar el fit para los rezagos

    print(str(resultados.summary()))
    #print(str(resultados.summary()).split("=\n\nResults")[0])
    printmd(bold("Pvalues:"))
    display(resultados.pvalues[[resultados.names[0]]].T)

    # Saco esto porque para normalizar tengo que calcular los irf a mano
    #irf = resultados.irf(periods=10)

    printmd(bold("Impulso-respuesta:"))
    var_respuesta = ["ipc"] if "ipc" in resultados.names else []
    var_respuesta += ["E" if "E" in resultados.names else "Ebc"]
    #var_respuesta += ["E" if "E" in resultados.names else "Efmi"]
    #var_respuesta = resultados.names #OJO! MUY LENTO! Grafica impulso-respuesta de todas las variables contra todas.
    #var_respuesta = ['ipc', 'E', 'Ebc', 'emae', 'immp_usa', 'Psoja_USA']:

    for var in resultados.names:
      for var_res in var_respuesta: #Grafica los impulso-respuesta de todas las variables contra ipc y contra E
        if var != var_res:
          sigma = resultados.sigma_u.loc[var,var]
          print(f"Sigma {var}: {sigma}")
          plotIrfWithSignif(signifs=[0.05, 0.32], impulse=var, response=var_res, var_results=resultados, periods=10, orth=True, cumulative=False, figsize=(3,2))
          plotIrfWithSignif(signifs=[0.05, 0.32], impulse=var, response=var_res, var_results=resultados, periods=10, orth=True, cumulative=True, figsize=(3,2))
        
    vars_contra_emae_o_pbi = ["Psoja_USA","Pmaíz_USA","Ptrigo_USA","tot_04","TOTfmi","impp_usa","E","Ebc", "ipc"]
    emae_o_pbi = "emae" if "emae" in resultados.names else "pbird"
    for var_contra in vars_contra_emae_o_pbi: #Para todas las variables de la lista anterior
      if var_contra in resultados.names: #Si son variables del dataset sobre el que estamos trabajando
        sigma = resultados.sigma_u.loc[var_contra,var_contra]
        print(f"Sigma {var_contra}: {sigma}")
        plotIrfWithSignif(signifs=[0.05, 0.32], impulse=var_contra, response=emae_o_pbi, var_results=resultados, periods=10, orth=True, cumulative=False, figsize=(3,2))
        plotIrfWithSignif(signifs=[0.05, 0.32], impulse=var_contra, response=emae_o_pbi, var_results=resultados, periods=10, orth=True, cumulative=True, figsize=(3,2))
        
    printmd(bold("Resultados descomposición de varianza (FEVD)"))
    try:
      fevd = resultados.fevd()
      fevd.plot(figsize=(5,10))
      fevd.summary()
    except Exception as inst:
      print("[Fallo el test!] motivo:")
      print(inst)

    printmd(bold("Test normality:"))
    tryrun(lambda: print(resultados.test_normality().summary()))
    #print(resultados.test_causality('ipc', ["E", "emae", "impp_usa"], kind='f'))
    printmd(bold("Test causality:"))
    tryrun(lambda: print(resultados.test_causality("ipc", removeInplace(resultados.names, "ipc"), kind='f').summary()))

    printmd(bold("Test causality:"))
    tryrun(lambda: print(resultados.test_causality("ipc", "emae", kind='f').summary()))

    printmd(bold("Test causality:"))
    tryrun(lambda: print(resultados.test_causality("emae", "ipc", kind='f').summary()))


    printmd(bold("Test whiteness:"))
    tryrun(lambda: print(resultados.test_whiteness(nlags=6, signif=0.05, adjusted=True).summary()))
    print("")
    printmd(bold("Raíces:"))
    tryrun(lambda: print(resultados.roots))
    print("")

  elif rModel == "VECM":
    modelo=VECM(Y,X)

    lor = modelo.select_order(6)
    print(f"Selected orders are: {lor.selected_orders}")
    display(lor.summary()) #Imprime tabla VAR Order Selection

    resultados=modelo.fit(maxlags)
    printRes = f"Det. terms outside the coint {str(resultados.summary()).split('Det. terms outside the coint')[1]}"
    printRes += f"\n                Loading coefficients (alpha){str(resultados.summary()).split('                Loading coefficients (alpha)')[1]}"
    printRes += f"\n          Cointegration relations {str(resultados.summary()).split('          Cointegration relations')[1]}"
    print(printRes)
    irf = resultados.irf(periods=10)

    for var in resultados.names:
        for var_res in ['ipc', 'E']:
            irf.plot(orth=True, impulse=var, response=var_res, figsize=(3,2))
            irf.plot_cum_effects(orth=True, impulse=var, response=var_res, figsize=(3,2))
    #fevd = resultados.fevd()
    #fevd.summary()
  return modelo, resultados


def correrTestsADF(df, maxlag=None, regression="c", variables=None):
  """
  Realiza test ADF y muestra los resultados en una tabla legible.
  Parametros:
  `df`: el dataframe sobre el que aplicar ADF
  `maxlag` (opcional): máximo lag incluido en el test, valor por defecto el mismo de la funciòn adfuller de statsmodels
  `regression` (opcional): orden de constante y tendencia a incluir en ADF. Valor por defecto mismo de la funcion adfuller de statsmodels ("c")
    - "c" : solo constante (default).
    - "ct" : constante y tendencia.
    - "ctt" : constante y tendencia lineal y cuadrática.
    - "n" : sin constante ni tendencia.
  `variables` (opcional): las variables sobre las que aplicar ADF. Valor por defecto: df.columns[2:].
  """
  if variables is None:
    variables = df.columns[2:]

  df_adf = pd.DataFrame({})
  df_temp = df.dropna() #.dropna() funcion saca las filas que tienen algun valor NaN, que en este caso es solo la primera

  if (df.isna().sum().sum() > 0): #cuento la cantidad de nans que aparecen en el dataframe
    only_na = df[~df.index.isin(df_temp.index)]
    print("WARNING: hay NaNs en el dataframe! Se excluyen las siguientes filas de los tests:")
    display(only_na)
    print("===")

  for variable in variables:
    adf = adfuller(df_temp[variable], maxlag=maxlag, regression=regression)
    df_adf.loc["result",variable] = adf[0]
    df_adf.loc["1%",variable] = adf[4]["1%"]
    df_adf.loc["5%",variable] = adf[4]["5%"]
    df_adf.loc["10%",variable] = adf[4]["10%"]
    df_adf.loc["rechaza 1%", variable] = adf[0] < adf[4]["1%"]
    df_adf.loc["rechaza 5%", variable] = adf[0] < adf[4]["5%"]
    df_adf.loc["rechaza 10%", variable] = adf[0] < adf[4]["10%"]
    df_adf.loc["pvalue",variable] = adf[1]
    df_adf.loc["lag",variable] = adf[2]
  return df_adf

def ADFmultipleLags(df, lags=[0,1,2,None], regression="c", variables=None):
  """
  Corre ADFtable muchas veces, con distintos valores de maxlag.
  `df`: el dataframe sobre el que aplicar ADF
  `lags`: los valores de maxlag a correr. Por defecto son [0,1,2,None]. None es el valor por defecto de adfuller de statsmodels, que es uno bueno calculado con cierta ecuacion.
  `regression` (opcional): orden de constante y tendencia a incluir en ADF. Valor por defecto mismo de la funcion adfuller de statsmodels ("c"). Opciones:
    - "c" : solo constante (default).
    - "ct" : constante y tendencia.
    - "ctt" : constante y tendencia lineal y cuadrática.
    - "n" : sin constante ni tendencia.
  `variables` (opcional): las variables sobre las que aplicar ADF. Valor por defecto: df.columns[2:].
  """
  for lag in lags:
    printmd(bold(f"Maxlag {lag}, regression {regression}"))
    display(correrTestsADF(df, maxlag=lag, regression=regression, variables=variables))