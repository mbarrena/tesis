#from utilities import printmd, bold
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
#from statsmodels.tsa.vector_ar.irf import IRAnalysis
import statsmodels.tsa.vector_ar.plotting as plotting
import localprojections as lp
import plotly.subplots as sp
import plotly.graph_objects as go

#NOTA: svar no está implementado

def getNormalizedIRFs(var_results, impulse, periods, orth, cumulative):
    # Basado en https://stackoverflow.com/questions/62269621/statsmodels-vector-ar-and-iranalysis
    sigma = var_results.sigma_u.loc[impulse,impulse]
    sigma = math.sqrt(sigma)
    if orth:
        J = var_results.orth_ma_rep(periods)
    else:
        J = var_results.ma_rep(periods)

    if J=="LocalProjections":
        impulse_coeff = var_results.params[impulse]
        J = np.zeros(periods)
        J[:] = impulse_coeff / sigma
    else:
        J = np.array(J)/sigma
        
    if cumulative:
        J = J.cumsum(axis=0)
    return J

def getNormalizedIRFerrband(var_results, impulse, periods, signif, orth, cum=False, stderr='asym'):
    sigma = var_results.sigma_u.loc[impulse,impulse]
    sigma = math.sqrt(sigma)
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

def plotIrfWithSignif(var_results, signifs: list, impulse: str, response: str, periods: int, orth: bool, cumulative: bool=False, **kwargs):
    """Abstract interface for specific plotIrf functions for VAR or LocalProjections"""
    raise NotImplementedErrorpl

def plotIrfWithSignifVAR(var_results, signifs: list, impulse: str, response: str, periods: int, orth: bool, cumulative: bool=False, **kwargs):
    stderr_type = 'mc'
    # Basado en la función plot de https://www.statsmodels.org/dev/_modules/statsmodels/tsa/vector_ar/irf.html#IRAnalysis
    irfs = getNormalizedIRFs(var_results, impulse, periods, orth, cumulative)
    title = f'Impulse responses'
    if cumulative:
        title = 'Cumulative ' + title.lower()
    if orth:
        title += ' (orthogonalized)'
    
    for signif in signifs:
        title_s = title + f' ({(1-signif)*100:.0f}%)'
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

def plotIrfWithSignifLP(var_results, signifs: list, impulse: str, response: str, periods: int, orth: bool, cumulative: bool=False, **kwargs):
    # Afaik los resultados de LP ya vienen normalizados (ya son porcentajes) así que normalizar no haría falta.
    title = f'{impulse} -> {response} ({{signif:.0f}}%)'
    title += '\nCumulative impulse responses' if cumulative else '\nImpulse responses'
    title += ' (orthogonalized)' if orth else ""
    
    for signif, irf in zip(signifs, var_results):
        title_s = title.format(signif=(1-signif)*100)
        df = irf[(irf["Shock"] == impulse) & (irf["Response"] == response)]
        if cumulative:
            df.loc[:,"Mean"] = df.loc[:,"Mean"].cumsum(axis=0)
            df.loc[:,"LB"] = df.loc[:,"LB"].cumsum(axis=0)
            df.loc[:,"UB"] = df.loc[:,"UB"].cumsum(axis=0)
        
        print(f"{'Cumulative ' if cumulative else ''}IRF values for significance {(1-signif)*100}")
        print(irf)
        
        lp.IRFPlot(irf=df, # take output from the estimated model
                    response=[response], # plot only response of invest ...
                    shock=[impulse], # ... to shocks from all variables
                    n_columns=1, # max 2 columns in the figure
                    n_rows=1, # max 2 rows in the figure
                    maintitle=title_s, # self-defined title of the IRF plot
                    show_fig=True, # display figure (from plotly)
                    save_pic=False, # don't save any figures on local drive
                    font_size=16
                    )

def plotIrfWithSignifLPthr(var_results, signifs: list, impulse: str, response: str, periods: int, orth: bool, cumulative: bool=False, **kwargs):
    # Afaik los resultados de LP ya vienen normalizados (ya son porcentajes) así que normalizar no haría falta.
    title = f'{impulse} -> {response} ({{signif:.0f}}%)'
    title += '\nCumulative impulse responses' if cumulative else '\nImpulse responses'
    title += ' (orthogonalized)' if orth else ""

    for signif, irf in zip(signifs, var_results):
        title_s = title.format(signif=(1-signif)*100)
        df_on = irf[0][(irf[0]["Shock"] == impulse) & (irf[0]["Response"] == response)]
        df_off = irf[1][(irf[1]["Shock"] == impulse) & (irf[1]["Response"] == response)] 
        
        if cumulative:
            df_off.loc[:,"Mean"] = df_off.loc[:,"Mean"].cumsum(axis=0)
            df_off.loc[:,"LB"] = df_off.loc[:,"LB"].cumsum(axis=0)
            df_off.loc[:,"UB"] = df_off.loc[:,"UB"].cumsum(axis=0)

        print(f"{'Cumulative ' if cumulative else ''}IRF values for significance {(1-signif)*100}")
        print(pd.concat([df_off, df_on], keys=["Pre-threshold", "Post-threshold"]))
        
        fig1 = lp.IRFPlot(irf=df_off, # take output from the estimated model
                    response=[response], # plot only response of invest ...
                    shock=[impulse], # ... to shocks from all variables
                    n_columns=1, # max 2 columns in the figure
                    n_rows=1, # max 2 rows in the figure
                    maintitle=title_s + ", pre-threshold", # self-defined title of the IRF plot
                    show_fig=False, # display figure (from plotly)
                    save_pic=False, # don't save any figures on local drive
                    font_size=16
                    )
        
        if cumulative:
            df_on.loc[:,"Mean"] = df_on.loc[:,"Mean"].cumsum(axis=0)
            df_on.loc[:,"LB"] = df_on.loc[:,"LB"].cumsum(axis=0)
            df_on.loc[:,"UB"] = df_on.loc[:,"UB"].cumsum(axis=0)
        
        fig2 = lp.IRFPlot(irf=df_on, # take output from the estimated model
                    response=[response], # plot only response of invest ...
                    shock=[impulse], # ... to shocks from all variables
                    n_columns=1, # max 2 columns in the figure
                    n_rows=1, # max 2 rows in the figure
                    maintitle=title_s + ", post-threshold", # self-defined title of the IRF plot
                    show_fig=False, # display figure (from plotly)
                    save_pic=False, # don't save any figures on local drive
                    font_size=16
                    )
        
        fig = sp.make_subplots(rows=1, cols=2, subplot_titles=["Pre-threshold", "Post-threshold"], 
            horizontal_spacing=0.05  # Reduce the spacing between subplots
            )

        # Add traces from fig1 to the left subplot (col=1)
        for trace in fig1.data:
            trace.line.color = 'black'  # Set line color to black for the left plot
            fig.add_trace(trace, row=1, col=1)

        # Add traces from fig2 to the right subplot (col=2)
        for trace in fig2.data:
            trace.line.color = 'crimson'  # Set line color to crimson for the right plot
            fig.add_trace(trace, row=1, col=2)

        # Update the layout
        fig.update_layout(
            title=title_s,
            plot_bgcolor="white",
            hovermode="x unified",
            showlegend=False,
            font=dict(color="black", size=16),
        )
        fig.show()