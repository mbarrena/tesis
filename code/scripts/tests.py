from statsmodels.tools.sm_exceptions import InterpolationWarning
from statsmodels.tsa.stattools import range_unit_root_test, zivot_andrews, kpss
from sklearn.metrics import mean_squared_error
from scripts.utilities import *
import warnings
import numpy as np
import pandas as pd
from arch.unitroot import PhillipsPerron
import chowtest
import ruptures as rpt
import matplotlib.pyplot as plt


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

def ChowTest(datos, last_index, first_index, y_vars=None, variables=None):
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
      last_index (int): el último índice (número de fila) a incluir PREVIO al quiebre estructural.
      first_index (int): el último índice (número de fila) a incluir POSTERIOR al quiebre estructural.
      y_vars (list[str]): una lista de columnas que representa las variables (y) contra las cuales evaluar. Por defecto se hace contra todas las variables. 
      variables (list[str]): subset of variables of the dataframe to use.

  Returns:
      chow_value (float): the chow test output value.
      p_value (float): the associated p-value for the chow test.
  """
  variables = datos.columns[2:] if variables is None else variables
  datos = datos[variables]
  y_vars = variables if y_vars is None else y_vars
  printmd(f"{bold('Test de Chow:')} idx before break {last_index}, idx after break {first_index}")
  for var in y_vars:
    printmd(f"{bold(var)}:")
    ct = chowtest.ChowTest(datos.drop(var, axis=1), datos[var],
          last_index_in_model_1=last_index,
          first_index_in_model_2=first_index)
  print(f"ChowTest value: {ct[0]:.4f}, p_value: {ct[1]:4f}")

def BaiPerron(datos, variables=None, cost_fun="l2", penalty=3, n_breaks=None, **kwargs):
  """
  Implementa el test de Bai-Perron para detectar múltiples cambios estructurales en series temporales.

  Args:
      datos (pd.DataFrame): El conjunto de datos.
      variables (list[str], opcional): Lista de variables a analizar. Si es None, se analizan todas las columnas excepto las dos primeras.
      cost_fun (str): La función de costo a usar. Por defecto es 'l2'. Opciones: 'l2', 'l1', 'rbf', etc.
      penalty (int): penalidad para busqueda de nùmero de quiebres (solo vale si no se eligion n_breaks). 3 por defecto.
      n_breaks (int): Número de rupturas, si es conocido. Por defecto None (no es conocido), en cuyo caso se usa Pelt para esncontrar quiebres. Si es conocido, se usa segmentacion binaria. None por defecto.

  Returns:
      pd.DataFrame: DataFrame con los resultados (índices de las rupturas y el número de rupturas detectadas).
  """
  variables = datos.columns[2:] if variables is None else variables
  resultados = []

  print(f"Test de Bai-Perron (Algoritmo de detecciòn: {'Pelt' if n_breaks is None else 'Prog Din.'}, Nro. quiebres: {'no especificado' if n_breaks is None else n_breaks}, Función de costo: {cost_fun} ):")

  for var in variables:
    print(f"Variable: {var}")
    # Extraer la serie
    serie = datos[var].dropna().values
    # Configurar el modelo de búsqueda de rupturas
    if n_breaks is None:
      algo = rpt.Pelt(model=cost_fun, **kwargs).fit(serie)
      bkps = algo.predict(pen=penalty)  # Puedes ajustar el penalizador
    else:
      algo = rpt.Dynp(model=cost_fun, **kwargs).fit(serie)
      bkps = algo.predict(n_bkps=n_breaks)

    errores = []
    for i in range(len(bkps) - 1):
      inicio, fin = bkps[i], bkps[i + 1]
      segmento = serie[inicio:fin]
      prediccion = np.mean(segmento)  # Usar la media como modelo
      error = mean_squared_error(segmento, [prediccion] * len(segmento))
      errores.append(error)

    error_total = np.mean(errores)
    print(f"Puntos de ruptura detectados: {bkps[:-1]}")
    #print(f"Errores por segmento: {errores}")
    print(f"Error total (MSE): {error_total:.4f}")

    # Guardar resultados
    resultados.append({
      "Variable": var,
      "Breakpoints": bkps,
      "Number of Breaks": len(bkps) - 1,
      "Error Total (MSE)": error_total
    })

    # Graficar la serie y los puntos de ruptura
    plt.figure(figsize=(10, 4))
    plt.plot(serie, label=f"Serie: {var}", color="blue")
    for bkp in bkps[:-1]:  # Evitar incluir el índice final
      plt.axvline(bkp, color="red", linestyle="--", label=f"Ruptura: {bkp}")
      plt.text(bkp, serie[bkp], f"{bkp}", color="red", fontsize=10, verticalalignment='bottom')

    plt.title(f"Rupturas detectadas en {var}")
    plt.xlabel("Índice")
    plt.ylabel("Valor")
    plt.legend()
    plt.grid()
    plt.show()

  # Crear DataFrame de resultados
  resultados_df = pd.DataFrame(resultados)
  return resultados_df
