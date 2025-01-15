from statsmodels.tsa.deterministic import Seasonality
import pandas as pd 
import numpy as np

def encontrarTrimestre(trim: str):
  res = 1
  if("2" in trim):
    res = 2
  if("3" in trim):
    res = 3
  if("4" in trim):
    res = 4
  return res

def findDateCols(data: pd.DataFrame):
  cols = [x.lower() for x in data.columns]
  try:
    año_k = data.columns[cols.index("año")]
  except ValueError:
    año_k = None
  try:
    mes_k = data.columns[cols.index("mes")]
  except ValueError:
    mes_k = None
  try:
    trimestre_k = data.columns[cols.index("trimestre")]
  except ValueError:
    trimestre_k = None

  print(año_k, mes_k, trimestre_k)
  return año_k, mes_k, trimestre_k

# return  -> (pd.PeriodIndex, str)
def generatePeriodIndex(data: pd.DataFrame):
  año_k, mes_k, trimestre_k = findDateCols(data)
  if trimestre_k is not None:
    trim_inicio = encontrarTrimestre(data[trimestre_k].iloc[0])
    trim_fin = encontrarTrimestre(data[trimestre_k].iloc[-1])
    anios_unicos = data[año_k].unique()
    idx_anio = 0

    trim = []
    anios = []
    if len(anios_unicos) == 1:
      trim = [i for i in range(trim_inicio, trim_fin+1)]
      anios = [int(anios_unicos[0]) for _ in range(len(trim))]
    else:
      #primer anio
      trim = [i for i in range(trim_inicio, 5)]    #Trimestres son de 1 a 4, así que pongo range con limite 5, ya que limite superior es exclusive
      anios = [int(anios_unicos[0]) for _ in range(len(trim))]
      idx_anio += 1
      while(idx_anio < len(anios_unicos)-1):
        #mientras no estemos en el ultimo año
        trim += [1,2,3,4]
        anios += [int(anios_unicos[idx_anio]) for _ in range(4)]
        idx_anio += 1
      #ultimo anio
      trim += [i for i in range(1, trim_fin+1)]
      anios += [int(anios_unicos[idx_anio]) for _ in range(len(range(1, trim_fin+1)))]

    return pd.PeriodIndex(year=anios, quarter=trim, freq="Q"), "Q"
  else:
    if año_k is not None and mes_k is not None:
      return pd.PeriodIndex(year=data[año_k], month=data[mes_k], freq="M"), "M"
    else:
      raise NotImplementedError("No se encontro año, mes ni trimestre en el dataset. Asegurese de que las columnas correpsondientes esten bien nombradas.")

def generatePeriodDummies(data: pd.DataFrame, index=None):
  if index is None:
    index, _ = generatePeriodIndex(data)
  dummies = Seasonality.from_index(index).in_sample(index).reset_index(drop=True)
  return dummies

def makeLogColumns(lista: list[any], data: pd.DataFrame):
  df = data.copy()
  for c in lista:
    df[c] = np.log(df[c].astype(float))
  return df

def makeDiffColumns(lista: list[any], data: pd.DataFrame, factor: int=1, drop: bool=True):
  df = data.copy()
  for c in lista:
    df[c] = df[c].diff(factor)
  if drop == True:
    df = df.iloc[factor: , :]
  return df

def makeMovingAverageColumns(columns: list[str], dataframe: pd.DataFrame, windowSize: int=12):
  df_res = pd.DataFrame({})

  periods = []
  for i in range(dataframe.shape[0] - windowSize + 1):
      period = str(dataframe.iloc[i]['trimestre'][0]) + "Tr" + str(dataframe.iloc[i]["año"])
      period += " a "
      period += str(dataframe.iloc[i+windowSize-1]['trimestre'][0]) + "Tr" + str(dataframe.iloc[i+windowSize-1]["año"])
      periods.append(period)
  df_res["trimestres"] = periods

  for var in columns:
    means = []
    for i in range(dataframe.shape[0] - windowSize + 1):
      mean = dataframe.iloc[i:i+windowSize][var].mean()
      means.append(mean)
    df_res[var] = means

  return df_res

