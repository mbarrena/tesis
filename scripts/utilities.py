from IPython.display import Markdown, display
import pandas as pd

def tryrun(func):
  try:
    func()
  except Exception as inst:
    print("[Fallo el test!] motivo:")
    print(inst)

def printmd(string):
    display(Markdown(string))

def bold(str):
  return f"**{str}**"

def mono(str):
  return f"`{str}`"

def renameColumnOfDataframe(dataframe: pd.DataFrame, oldColumnName: str, newColumnName: str):
  """Funcion para renombrar una columna de un dataframe, "oldColumnName". 
  
  Le pone el nombre "newColumnName". 
  No modifica el dataframe que se pasa a la funcion, solamente devuelve un nuevo dataframe igual pero con la columna renombrada.
  
  **Modo de uso:**
  ```python
  df = renameColumnOfDataframe(df, "mi columna mal nombrada", "mi buen nombre de columna")
  ```
  **Nota:** si la llamas asi:
  ```python
  renameColumnOfDataframe(df, "mi columna mal nombrada", "mi buen nombre de columna")
  ```
  (sin el "df =" adelante) el cambio no se va a guardar en el dataframe en el que estas trabajando.
  """
  renamedDataframe = dataframe.copy()
  renamedDataframe[newColumnName] = renamedDataframe[oldColumnName]
  return renamedDataframe.drop([oldColumnName], axis=1)

def removeInplace(list: list[any], val: any):
  """Devuelve una lista igual a list pero sin ninguna de las apariciones de val.

  **NO MODIFICA LA LISTA DE ENTRADA** (nombre mal puesto). 
  """
  return [x for x in list if x != val]
