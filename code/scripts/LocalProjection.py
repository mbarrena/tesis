from visualization import plotIrfWithSignifLP
from scipy.stats import chi2
from statsmodels.stats.diagnostic import acorr_ljungbox
import pandas as pd

class LPResults(list):
        names: list[str] = None
        signifs: list[float] = None
        lags: int = None
        sigma_u: None = None #LocalProjection results are already normalised
        roots = "NOT YET IMPLEMENTED"
        
        def plotIrfWithSignif(signifs: list, impulse: str, response: str, var_results, periods: int, orth: bool, cumulative: bool=False, **kwargs):
            return plotIrfWithSignifLP(signifs=signifs, impulse=impulse, response=response, var_results=var_results, periods=periods, orth=orth, cumulative=cumulative, **kwargs)
        
        def __repr__(self):
            return f"LPResults([{', '.join([f'irf_df (signif {s})' for s in self.signifs])}])"
        
        def test_whiteness(self, nlags: int, signif: float, adjusted: bool):
            """
            Perform the Ljung-Box test for whiteness of residuals.
                residuals: Residuals from the model.
                lags: Number of lags to consider in the test.
                alpha: Significance level for the test.

            :return: DataFrame summarizing the Ljung-Box test results.
            """
            class whitenessResults:
                summary = lambda *args, **kwargs: "NOT YET IMPLEMENTED"
                plot = lambda *args, **kwargs: "NOT YET IMPLEMENTED"
            return whitenessResults()
            results = []
            alpha = signif
            
            irf = self[self.signifs.index(alpha)]

            residuals = irf['Mean'] - irf['Mean'].shift(1)  # Example residual calculation
            lb_test = acorr_ljungbox(residuals, lags=6, return_df=True)
            #lb_test['Critical Value'] = chi2.ppf(1 - alpha, df=lb_test.index)  # df = lags
            lb_test['Significant'] = lb_test['lb_pvalue'] < alpha

            results = []
            for lag in lb_test.index:
                results.append({
                    'Significance': alpha,
                    'Lag': lag,
                    'Test Statistic': lb_test.loc[lag, 'lb_stat'],
                #    'Critical Value': lb_test.loc[lag, 'Critical Value'],
                    'p-value': lb_test.loc[lag, 'lb_pvalue'],
                    'Is significant': lb_test.loc[lag, 'Significant']
                })
                
            return pd.DataFrame(results)
        
        def fevd(*args, **kwargs):
            class fevdResults:
                summary = lambda *args, **kwargs: "NOT YET IMPLEMENTED"
                plot = lambda *args, **kwargs: "NOT YET IMPLEMENTED"
            return fevdResults()
        
        def test_normality(*args, **kwargs):
            class normalityResults:
                summary = lambda *args, **kwargs: "NOT YET IMPLEMENTED"
            return normalityResults()

        def test_causality(*args, **kwargs):
            class causalityResults:
                summary = lambda *args, **kwargs: "NOT YET IMPLEMENTED"
            return causalityResults()
        
        @property
        def pvalues(self):
            """
            Return the p-values from the local projection results. 
            This assumes `irf` is a DataFrame with coefficients and p-values.
            """
            return "NOT SURE HOW TO IMPLEMENT..."