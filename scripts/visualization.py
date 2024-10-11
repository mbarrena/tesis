from utilities import printmd, bold
import matplotlib.pyplot as plt
import numpy as np
from statsmodels.tsa.vector_ar.irf import IRAnalysis
import statsmodels.tsa.vector_ar.plotting as plotting

#NOTA: svar no está implementado

def getNormalizedIRFs(var_results, impulse, periods, orth, cumulative):
  # Basado en https://stackoverflow.com/questions/62269621/statsmodels-vector-ar-and-iranalysis
  sigma = var_results.sigma_u.loc[impulse,impulse]
  if orth:
    J = var_results.orth_ma_rep(periods)
  else:
    J = var_results.ma_rep(periods)
  J = np.array(J)/sigma
  if cumulative:
    J = J.cumsum(axis=0)
  return J

def getNormalizedIRFerrband(var_results, impulse, periods, signif, orth, cum=False, stderr='asym'):
  sigma = var_results.sigma_u.loc[impulse,impulse]
  if stderr == "mc":
    G,H = var_results.irf_errband_mc(orth=orth, repl=100, steps=periods, signif=signif, seed=None, burn=100, cum=cum)
    G = G/sigma
    H = H/sigma
    return (G,H)
  elif stderr =="asym":
    ir_analysis = var_results.irf(periods)
    cov = ir_analysis.cov(orth=orth)
    cov = cov/sigma
    return cov

def plotIrfWithSignif(signifs, impulse, response, var_results, periods, orth, cumulative=False, **kwargs):
  stderr_type = 'mc'
  # Basado en la función plot de https://www.statsmodels.org/dev/_modules/statsmodels/tsa/vector_ar/irf.html#IRAnalysis
  irfs = getNormalizedIRFs(var_results, impulse, periods, orth, cumulative)
  title = f'Impulse Responses'
  if cumulative:
    title = 'Cumulative' + title
  if orth:
    title += '(orthogonalized)'
  for signif in signifs:
    title_s = title + f'({(1-signif)*100:.0f}%)'
    stderr = getNormalizedIRFerrband(var_results, impulse, periods, signif, orth, cum=cumulative, stderr=stderr_type)
    fig = plotting.irf_grid_plot(irfs, stderr, impulse, response,
                                     var_results.names, title=title_s, 
                                     signif=signif,
                                     stderr_type=stderr_type,
                                     **kwargs)
  if cumulative:
    vals = plt.gca().lines[0].get_ydata()
    print(f"Cum effect {impulse} (impulso) - {response} (respuesta): {vals}")
  return fig

  # Obtengo la matriz de covarianza (?) del VAR, y de ahí saco el desvío correspondiente a la variable impulso
  # Normalizo los IRFs con dicho desvío
  irf.irfs = irf.irfs/sigma
    
  for signif in signifs:
    fig = irf.plot(signif=signif, impulse=impulse, **kwargs)
    for plot in fig.get_axes(): plot.title.set_text(f'{plot.title.get_text()} ({(1-signif)*100:.0f}%)')

def plotCumulativeIrfWithSignif(irf, signifs, impulse, response, orth, var_results, **kwargs):
  # Obtengo la matriz de covarianza (?) del VAR, y de ahí saco el desvío correspondiente a la variable impulso
  sigma = var_results.sigma_u.loc[impulse,impulse]
  # Normalizo los IRFs cumulativos con dicho desvío
  irf.irfs = irf.irfs/sigma
  
  for signif in signifs:
      fig = irf.plot_cum_effects(signif=signif, impulse=impulse, response=response, orth=orth, **kwargs)
      for plot in fig.get_axes():
          plot.title.set_text(f'{plot.title.get_text()} ({(1-signif)*100:.0f}%)')
      vals = plt.gca().lines[0].get_ydata()
  
  print(f"Cum effect {impulse} (impulso) - {response} (respuesta): {vals}")
