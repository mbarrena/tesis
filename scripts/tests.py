from statsmodels.tools.sm_exceptions import InterpolationWarning
from utilities import *
import warnings
import numpy as np
import pandas as pd
from statsmodels.tsa.stattools import range_unit_root_test, zivot_andrews, kpss
from arch.unitroot import PhillipsPerron

# Funciones para hacer que el resultado de los tests RUR, Zivot-Andrews y KPSS sea más legible.
def formatInterpolationWarning(fun):
  """Catch p-value InterpolationWarning that may arise during root testing and reformats p-value accordingly."""
  warnings.filterwarnings("error", category=InterpolationWarning)
  try:
    y = fun()
    warnings.resetwarnings()
    return y
  except InterpolationWarning as err:
    warnings.filterwarnings("ignore", category=InterpolationWarning)
    y = list(fun())
    warnings.resetwarnings()
    if "actual p-value is larger" in str(err):
      y[1] = f">={y[1]}"
    elif "actual p-value is greater" in str(err):
      y[1] = f">={y[1]}"
    elif "actual p-value is smaller" in str(err):
      y[1] = f"<={y[1]}"
    else:
      warnings.warn(str(err))
    return y

def _runAndFormatKPSSOutput(datos, var):
  """Reformats KPSS test output to be more legible

  Args:
      var: nombre de la variable evaluada, para imprimirlo con los resultados.
      fun: la funcion a ejecutar para evaluar la variable.
  """
  y = formatInterpolationWarning(lambda: kpss(datos[var]))
  res = {
      "kpss_stat": y[0],
      "p_value": y[1],
      "lags" : y[2],
      "crit" : y[3]
  }
  return pd.DataFrame.from_dict(res, orient='index')
  #return f"{bold(var)}: {mono(str(res))}"

def _runAndFormatZivotAndrewsOutput(datos, var):
  """Reformats Zivot-Andrews test output to be more legible

  Args:
      var: nombre de la variable evaluada, para imprimirlo con los resultados.
      fun: la funcion a ejecutar para evaluar la variable.
  """
  y = formatInterpolationWarning(lambda: zivot_andrews(datos[var]))
  res = {
      "za_stat": y[0],
      "p_value": y[1],
      "cvdict" : y[2],
      "baselag" : y[3],
      "bpidx" : y[4]
  }
  return pd.DataFrame.from_dict(res, orient='index')
  #return f"{bold(var)}: {mono(str(res))}"

def _runAndFormatRUROutput(datos, var):
  """Runs RUR test on data and reformats the output to be more legible

  Args:
      var: nombre de la variable evaluada, para imprimirlo con los resultados.
      fun: la funcion a ejecutar para evaluar la variable.
  """
  y = formatInterpolationWarning(lambda: range_unit_root_test(datos[var]))
  res = {
      "rur_stat": y[0],
      "p_value": y[1],
      "crit" : y[2]
  }
  return pd.DataFrame.from_dict(res, orient='index')
  #return f"{bold(var)}: {mono(str(res))}"

def _formatPPOutput(var, pp):
  """Reformats Phillips-Perron test output to be more legible

  Args:
      var: nombre de la variable evaluada, para imprimirlo con los resultados.
      pp: resultado de correr test de Phillips-Perron sobre variable.
  """
  cv_string = ""
  cv = pp.critical_values.keys()
  cv_numeric = np.array([float(x.split("%")[0]) for x in cv])
  cv_numeric = np.sort(cv_numeric)
  for val in cv_numeric:
      p = str(int(val)) + "%"
      cv_string += f"{pp.critical_values[p]:0.2f}"
      cv_string += " (" + p + ")"
      if val != cv_numeric[-1]:
          cv_string += ", "
  res = {
    "Test Statistic": f"{pp.stat:0.3f}",
    "P-value": f"{pp.pvalue:0.3f}",
    "Lags": f"{pp.lags:d}",
    "Critical Values": cv_string
  }
  return pd.DataFrame.from_dict(res, orient='index')
  #return f"{bold(var)}: {mono(str(res))}"

def _runTestAndFormat(fun, header, datos, variables):
  variables = datos.columns[2:] if variables is None else variables
  printmd(bold(header))
  dfs = []
  for var in variables:
    res = fun(datos, var)
    dfs.append(res)
  df = pd.concat(dfs, axis= 1, ignore_index=True)
  df.columns = variables
  display(df)

def RUR(datos,variables=None):
  return _runTestAndFormat(_runAndFormatRUROutput, "Range Unit Root Test (RUR):", datos, variables)
  #variables = datos.columns[2:] if variables is None else variables
  #printmd(bold("Range Unit Root Test (RUR):"))
  #for var in variables:
  #  printmd(runAndFormatRUROutput(var, lambda: range_unit_root_test(datos[var])))

def ZivotAndrews(datos, variables=None):
  return _runTestAndFormat(_runAndFormatZivotAndrewsOutput, "Test Zivot-Andrews:", datos, variables)
  #variables = datos.columns[2:] if variables is None else variables
  #printmd(bold("Test Zivot-Andrews:"))
  #for var in variables:
  #  printmd(runAndFormatZivotAndrewsOutput(var, lambda: zivot_andrews(datos[var])))

def KPSS(datos, variables=None):
  return _runTestAndFormat(_runAndFormatKPSSOutput, "Test KPSS:", datos, variables)
  #variables = datos.columns[2:] if variables is None else variables
  #printmd(bold("Test KPSS:"))
  #for var in variables:
  #  printmd(runAndFormatKPSSOutput(var, lambda: kpss(datos[var])))

def PP(datos, variables=None):
  variables = datos.columns[2:] if variables is None else variables
  printmd(bold("Phillips-Perron Test (PP):"))
  dfs = []
  for var in variables:
    pp = PhillipsPerron(datos[var])
    dfs.append(_formatPPOutput(var, pp))
    #printmd(formatPPOutput(var, pp))
  df = pd.concat(dfs, axis = 1, ignore_index=True)
  df.columns = variables
  print(f"Hipótesis nula: {pp.null_hypothesis}")
  print(f"Hipótesis alternativa: {pp.alternative_hypothesis}")
  return df

# MA: el test de chow necesita que uno elija cual cree que es el punto donde ocurre el quiebre estructural,
# y le pase a la función el índice del último dato antes del quiebre (last_index_in_model_1)
# y del primero después del quiebre (first_index_in_model_2).
# Como no sabia que poner aca, lo dejo comentado y descomentalo y ponele los valores que te parezcan.
# Referencias: https://pypi.org/project/chowtest/
# https://www.statology.org/chow-test-in-python/
# Tampoco estoy segura si hice bien en hacer que y sea cada una de las variables y x sean todas las demas.

def ChowTest(datos, last_index, first_index, significance, variables=None):
  """The chow test is computed and assessed against the significance argument. The chow_test value and p_value are returned
    from the function.
    
  **Nota:**
  El test de chow necesita que uno elija cual cree que es el punto donde ocurre el quiebre estructural,
  y le pase a la función el índice del último dato antes del quiebre (last_index_in_model_1)
  y del primero después del quiebre (first_index_in_model_2).

  Como no sabia que poner aca, lo dejo comentado y descomentalo y ponele los valores que te parezcan.
    
  **Referencias:** 
  https://pypi.org/project/chowtest/

  https://www.statology.org/chow-test-in-python/
    
  Tampoco estoy segura si hice bien en hacer que y sea cada una de las variables y x sean todas las demas.

  Args:
      datos: los datos. Se realiza test de chow sobre cada una de las variables especificadas de a una (y), contra todas las demas (x).
      last_index (int): the final index value to be included before the data split.
      first_index (int): the first index value to be included after the data split.
      significance_level (float): the significance level against which the p-value is assessed.

  Returns:
      chow_value (float): the chow test output value.
      p_value (float): the associated p-value for the chow test.
  """
  from chow_test import chow_test
  variables = datos.columns[2:] if variables is None else variables
  printmd(f"{bold('Test de Chow:')} idx before break {last_index}, idx after break {first_index}")
  for var in variables:
    printmd(f"{bold(var)}:")
    ct = chow_test(y_series=datos[var], X_series=datos.drop(var, axis=1),
          last_index=last_index,
          first_index=first_index,
          significance=significance)
  print(f"chow_test value: {ct[0]}, p_value: {ct[1]}")
