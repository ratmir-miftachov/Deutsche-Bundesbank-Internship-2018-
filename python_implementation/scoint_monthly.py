import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
from statsmodels.tsa.seasonal import seasonal_decompose
from statsmodels.tsa.stattools import adfuller
import warnings
warnings.filterwarnings('ignore')

class SeasonalCointegration:
    """
    Seasonal Cointegration Test - Monthly Frequency
    Extension of EGHL test from quarterly to monthly frequency
    
    Deutsche Bundesbank Internship Project 2018
    """
    
    def __init__(self):
        # Initialize with pre-computed distributions if available
        # In practice, these would be loaded from saved files
        self.distr = None  # HEGY distributions
        self.distr2 = None  # EGHL distributions
    
    def diff(self, x, lag=1):
        """Difference function equivalent to R's diff()"""
        if isinstance(x, pd.Series):
            result = x.diff(lag)
            return result
        elif isinstance(x, np.ndarray):
            result = np.concatenate([np.full(lag, np.nan), np.diff(x, lag)])
            return result
        else:
            raise ValueError("Input must be pandas Series or numpy array")
    
    def lag(self, x, k):
        """Lag function equivalent to R's lag()"""
        if isinstance(x, pd.Series):
            return x.shift(k)
        elif isinstance(x, np.ndarray):
            result = np.concatenate([np.full(k, np.nan), x[:-k]])
            return result
        else:
            raise ValueError("Input must be pandas Series or numpy array")
    
    def ftest(self, model_restricted, model_unrestricted):
        """F-test between nested models"""
        rss_r = model_restricted.ssr
        rss_u = model_unrestricted.ssr
        df_r = model_restricted.df_resid
        df_u = model_unrestricted.df_resid
        
        f_stat = ((rss_r - rss_u) / (df_r - df_u)) / (rss_u / df_u)
        return f_stat
    
    def pval(self, test_value, distribution, two_sided=True, side="left"):
        """Calculate p-values from empirical distribution"""
        if isinstance(test_value, (list, np.ndarray)):
            test_value = np.array(test_value)
        else:
            test_value = np.array([test_value])
        
        distribution = np.array(distribution)
        
        if two_sided:
            p_values = []
            for tv in test_value:
                p = np.mean(np.abs(distribution) > np.abs(tv))
                p_values.append(p)
        else:
            p_values = []
            for tv in test_value:
                if side == "left":
                    p = np.mean(distribution < tv)
                else:  # side == "right"
                    p = np.mean(distribution > tv)
                p_values.append(p)
        
        return p_values[0] if len(p_values) == 1 else p_values
    
    def interpolation(self, t_value, distributions, side, column_index):
        """Interpolate p-values based on sample size"""
        t = len(self.current_series)
        
        # Sample size intervals for interpolation
        intervals = [(0, 60), (60, 120), (120, 240), (240, 1200), (1200, 100000)]
        
        # Find appropriate interval
        interval_idx = 0
        for i, (lower, upper) in enumerate(intervals):
            if lower < t <= upper:
                interval_idx = i
                break
        
        if 1 < interval_idx < 4:
            # Interpolate between adjacent distributions
            lower_dist = distributions[interval_idx - 1][:, column_index]
            upper_dist = distributions[interval_idx][:, column_index]
            
            pval_lower = self.pval(t_value, lower_dist, two_sided=False, side=side)
            pval_upper = self.pval(t_value, upper_dist, two_sided=False, side=side)
            
            lower_bound, upper_bound = intervals[interval_idx]
            weight = (t - lower_bound) / (upper_bound - lower_bound)
            pval_interpolated = pval_lower + (pval_upper - pval_lower) * weight
        else:
            # Use boundary distributions
            if interval_idx <= 1:
                dist = distributions[0][:, column_index]
            else:
                dist = distributions[3][:, column_index]
            pval_interpolated = self.pval(t_value, dist, two_sided=False, side=side)
        
        return pval_interpolated
    
    def hegy_filter(self, series):
        """HEGY linear filters for monthly data (12 frequencies)"""
        s = pd.Series(series) if not isinstance(series, pd.Series) else series
        
        filters = []
        
        # Filter 1: y1 (frequency 0)
        y1 = s
        for i in range(1, 12):
            y1 = y1 + self.lag(s, i)
        filters.append(y1)
        
        # Filter 2: y2 (frequency 6/12)
        y2 = s
        for i in range(1, 12):
            y2 = y2 + ((-1)**i) * self.lag(s, i)
        y2 = -y2
        filters.append(y2)
        
        # Filter 3: y3 (frequency 3/12, 9/12)
        y3 = -sum([self.lag(s, i) for i in [1, 3, 5, 7, 9, 11]])
        filters.append(y3)
        
        # Filter 4: y4 (frequency 3/12, 9/12)
        y4 = s - self.lag(s, 2) + self.lag(s, 4) - self.lag(s, 6) + self.lag(s, 8) - self.lag(s, 10)
        y4 = -y4
        filters.append(y4)
        
        # Filters 5-12: Complex frequency filters
        # Filter 5
        y5 = s + self.lag(s, 1) - 2*self.lag(s, 2) + self.lag(s, 3) + self.lag(s, 4) - 2*self.lag(s, 5) + \
             self.lag(s, 6) + self.lag(s, 7) - 2*self.lag(s, 8) + self.lag(s, 9) + self.lag(s, 10) - 2*self.lag(s, 11)
        y5 = -0.5 * y5
        filters.append(y5)
        
        # Filter 6
        y6 = np.sqrt(3)/2 * (s - self.lag(s, 1) + self.lag(s, 3) - self.lag(s, 4) + 
                             self.lag(s, 6) - self.lag(s, 7) + self.lag(s, 9) - self.lag(s, 10))
        filters.append(y6)
        
        # Filter 7
        y7 = s - self.lag(s, 1) - 2*self.lag(s, 2) - self.lag(s, 3) + self.lag(s, 4) + 2*self.lag(s, 5) + \
             self.lag(s, 6) - self.lag(s, 7) - 2*self.lag(s, 8) - self.lag(s, 9) + self.lag(s, 10) + 2*self.lag(s, 11)
        y7 = 0.5 * y7
        filters.append(y7)
        
        # Filter 8
        y8 = -np.sqrt(3)/2 * (s + self.lag(s, 1) - self.lag(s, 3) - self.lag(s, 4) + 
                              self.lag(s, 6) + self.lag(s, 7) - self.lag(s, 9) - self.lag(s, 10))
        filters.append(y8)
        
        # Filter 9
        y9 = np.sqrt(3)*s - self.lag(s, 1) + self.lag(s, 3) - np.sqrt(3)*self.lag(s, 4) + 2*self.lag(s, 5) - \
             np.sqrt(3)*self.lag(s, 6) + self.lag(s, 7) - self.lag(s, 9) + np.sqrt(3)*self.lag(s, 10) - 2*self.lag(s, 11)
        y9 = -0.5 * y9
        filters.append(y9)
        
        # Filter 10
        y10 = s - np.sqrt(3)*self.lag(s, 1) + 2*self.lag(s, 2) - np.sqrt(3)*self.lag(s, 3) + self.lag(s, 4) - \
              self.lag(s, 6) + np.sqrt(3)*self.lag(s, 7) - 2*self.lag(s, 8) + np.sqrt(3)*self.lag(s, 9) - self.lag(s, 10)
        y10 = 0.5 * y10
        filters.append(y10)
        
        # Filter 11
        y11 = np.sqrt(3)*s + self.lag(s, 1) - self.lag(s, 3) - np.sqrt(3)*self.lag(s, 4) - 2*self.lag(s, 5) - \
              np.sqrt(3)*self.lag(s, 6) - self.lag(s, 7) + self.lag(s, 9) + np.sqrt(3)*self.lag(s, 10) + 2*self.lag(s, 11)
        y11 = 0.5 * y11
        filters.append(y11)
        
        # Filter 12
        y12 = s + np.sqrt(3)*self.lag(s, 1) + 2*self.lag(s, 2) + np.sqrt(3)*self.lag(s, 3) + self.lag(s, 4) - \
              self.lag(s, 6) - np.sqrt(3)*self.lag(s, 7) - 2*self.lag(s, 8) - np.sqrt(3)*self.lag(s, 9) - self.lag(s, 10)
        y12 = -0.5 * y12
        filters.append(y12)
        
        # Filter 13: First difference at lag 12
        y13 = s - self.lag(s, 12)
        filters.append(y13)
        
        return filters
    
    def eghl_filter(self, series):
        """EGHL linear filters for monthly cointegration"""
        s = pd.Series(series) if not isinstance(series, pd.Series) else series
        
        filters = []
        
        # Filter 1: Sum filter (frequency 0)
        y1 = s + sum([self.lag(s, i) for i in range(1, 12)])
        filters.append(y1)
        
        # Filter 2: Alternating sum (frequency 6/12)
        y2 = s
        for i in range(1, 12):
            y2 = y2 + ((-1)**i) * self.lag(s, i)
        y2 = -y2
        filters.append(y2)
        
        # Filter 3: Even lags (frequency 3/12, 9/12)
        y3 = s - self.lag(s, 2) + self.lag(s, 4) - self.lag(s, 6) + self.lag(s, 8) - self.lag(s, 10)
        y3 = -y3
        filters.append(y3)
        
        # Filter 4: Complex pattern 1
        y4 = s - self.lag(s, 1) + self.lag(s, 3) - self.lag(s, 4) + self.lag(s, 6) - \
             self.lag(s, 7) + self.lag(s, 9) - self.lag(s, 10)
        y4 = -y4
        filters.append(y4)
        
        # Filter 5: Complex pattern 2
        y5 = s + self.lag(s, 1) - self.lag(s, 3) - self.lag(s, 4) + self.lag(s, 6) + \
             self.lag(s, 7) - self.lag(s, 9) - self.lag(s, 10)
        y5 = -y5
        filters.append(y5)
        
        # Filter 6: Complex pattern with sqrt(3)
        y6 = s + 2*self.lag(s, 2) + self.lag(s, 4) - self.lag(s, 6) - 2*self.lag(s, 8) - \
             self.lag(s, 10) - np.sqrt(3)*self.lag(s, 1) - np.sqrt(3)*self.lag(s, 3) + \
             np.sqrt(3)*self.lag(s, 7) + np.sqrt(3)*self.lag(s, 9)
        y6 = -y6
        filters.append(y6)
        
        # Filter 7: Complex pattern with sqrt(3) variation
        y7 = s + 2*self.lag(s, 2) + self.lag(s, 4) - self.lag(s, 6) - 2*self.lag(s, 8) - \
             self.lag(s, 10) + np.sqrt(3)*self.lag(s, 1) + np.sqrt(3)*self.lag(s, 3) - \
             np.sqrt(3)*self.lag(s, 7) - np.sqrt(3)*self.lag(s, 9)
        y7 = -y7
        filters.append(y7)
        
        # Filter 8: Annual difference
        y8 = s - self.lag(s, 12)
        filters.append(y8)
        
        return filters
    
    def lag_selection_procedure(self, dependent_var, regressors, max_lags=24, significance_level=0.05):
        """Automatic lag selection procedure"""
        # Create lag matrix
        lags_data = pd.DataFrame()
        for i in range(1, max_lags + 1):
            lags_data[f'Lag{i}'] = self.lag(dependent_var, i)
        
        # Combine with regressors
        if isinstance(regressors, pd.DataFrame):
            full_data = pd.concat([regressors, lags_data], axis=1)
        else:
            full_data = lags_data.copy()
        
        # Remove rows with NaN values
        full_data = full_data.dropna()
        dependent_clean = dependent_var.dropna()
        
        # Align lengths
        min_len = min(len(dependent_clean), len(full_data))
        dependent_clean = dependent_clean.iloc[-min_len:]
        full_data = full_data.iloc[-min_len:]
        
        # Backward elimination
        while len(full_data.columns) > 0:
            # Fit regression
            X = sm.add_constant(full_data)
            try:
                model = sm.OLS(dependent_clean, X).fit()
                
                # Get p-values for lag terms (exclude constant and initial regressors)
                if isinstance(regressors, pd.DataFrame):
                    lag_start_idx = len(regressors.columns) + 1  # +1 for constant
                else:
                    lag_start_idx = 1  # Just constant
                
                if len(model.pvalues) > lag_start_idx:
                    lag_pvalues = model.pvalues.iloc[lag_start_idx:]
                    
                    # Find highest p-value
                    max_pval_idx = lag_pvalues.idxmax()
                    max_pval = lag_pvalues.max()
                    
                    # Remove if not significant
                    if max_pval > significance_level and len(full_data.columns) > 1:
                        full_data = full_data.drop(columns=[max_pval_idx])
                    elif max_pval > significance_level and len(full_data.columns) == 1:
                        # Remove last lag if not significant
                        full_data = pd.DataFrame()
                        break
                    else:
                        break
                else:
                    break
            except:
                break
        
        # Extract selected lag numbers
        if len(full_data.columns) > 0:
            selected_lags = [int(col.replace('Lag', '')) for col in full_data.columns if col.startswith('Lag')]
        else:
            selected_lags = []
        
        return selected_lags, full_data
    
    def hegy_test(self, series):
        """HEGY test for seasonal unit roots"""
        self.current_series = series
        filters = self.hegy_filter(series)
        
        # Dependent variable: annual difference
        series_diff = filters[12]  # y13 = (1-L^12)
        
        # Independent variables: lagged filtered series
        regressors = pd.DataFrame()
        for i in range(12):
            regressors[f'Reg{i+1}'] = self.lag(filters[i], 1)
        
        # Lag selection
        selected_lags, lag_data = self.lag_selection_procedure(series_diff, regressors)
        
        # Final regression
        if len(lag_data) > 0:
            X = sm.add_constant(pd.concat([regressors, lag_data], axis=1).dropna())
        else:
            X = sm.add_constant(regressors.dropna())
        
        y = series_diff.dropna()
        
        # Align data
        min_len = min(len(y), len(X))
        y = y.iloc[-min_len:]
        X = X.iloc[-min_len:]
        
        try:
            model = sm.OLS(y, X).fit()
            
            # Test statistics
            t1 = model.tvalues.iloc[1]  # First regressor
            t2 = model.tvalues.iloc[2]  # Second regressor
            
            # F-tests for pairs of coefficients
            def joint_test(model, indices):
                """Joint significance test for multiple coefficients"""
                coef_subset = model.params.iloc[indices]
                cov_subset = model.cov_params().iloc[indices, indices]
                
                if len(coef_subset) == 2:
                    chi2_stat = coef_subset.T @ np.linalg.inv(cov_subset) @ coef_subset
                    f_stat = chi2_stat / 2 * (model.df_resid) / (model.df_resid)
                    return f_stat
                return None
            
            f34 = joint_test(model, [3, 4]) if len(model.params) > 4 else 0
            f56 = joint_test(model, [5, 6]) if len(model.params) > 6 else 0
            f78 = joint_test(model, [7, 8]) if len(model.params) > 8 else 0
            f910 = joint_test(model, [9, 10]) if len(model.params) > 10 else 0
            f1112 = joint_test(model, [11, 12]) if len(model.params) > 12 else 0
            
            test_stats = {
                't1': t1, 't2': t2, 'f34': f34 or 0, 'f56': f56 or 0, 
                'f78': f78 or 0, 'f910': f910 or 0, 'f1112': f1112 or 0
            }
            
            return test_stats, selected_lags
            
        except Exception as e:
            print(f"Error in HEGY test: {e}")
            return {
                't1': 0, 't2': 0, 'f34': 0, 'f56': 0, 
                'f78': 0, 'f910': 0, 'f1112': 0
            }, []
    
    def eghl_cointegration_test(self, x, y):
        """EGHL cointegration test between two series"""
        x_series = pd.Series(x) if not isinstance(x, pd.Series) else x
        y_series = pd.Series(y) if not isinstance(y, pd.Series) else y
        
        filx = self.eghl_filter(x_series)
        fily = self.eghl_filter(y_series)
        
        test_results = {}
        selected_lags_all = []
        
        # Test 1: Frequency 0
        try:
            reg1 = sm.OLS(fily[0].dropna(), sm.add_constant(filx[0].dropna())).fit()
            res1 = pd.Series(reg1.resid, index=fily[0].dropna().index)
            
            # Lag selection for residuals
            res1_diff = self.diff(res1)
            lag_res1 = self.lag(res1, 1)
            
            selected_lags1, lag_data1 = self.lag_selection_procedure(res1_diff, pd.DataFrame())
            
            # Final regression
            if len(lag_data1) > 0:
                X1 = sm.add_constant(pd.concat([lag_data1, lag_res1], axis=1).dropna())
            else:
                X1 = sm.add_constant(lag_res1.dropna())
            
            y1 = res1_diff.dropna()
            min_len = min(len(y1), len(X1))
            y1, X1 = y1.iloc[-min_len:], X1.iloc[-min_len:]
            
            model1 = sm.OLS(y1, X1).fit()
            t1 = model1.tvalues.iloc[-1]  # Coefficient on lagged residual
            
            test_results['t1'] = t1
            selected_lags_all.append(selected_lags1)
            
        except Exception as e:
            test_results['t1'] = 0
            selected_lags_all.append([])
        
        # Test 2: Frequency 6/12
        try:
            reg2 = sm.OLS(fily[1].dropna(), sm.add_constant(filx[1].dropna())).fit()
            res2 = pd.Series(reg2.resid, index=fily[1].dropna().index)
            
            v = self.eghl_filter(res2)
            regres2 = -v[7] / v[1]  # -(1-L^12) / alternating sum
            
            selected_lags2, lag_data2 = self.lag_selection_procedure(regres2, pd.DataFrame())
            
            if len(lag_data2) > 0:
                X2 = sm.add_constant(pd.concat([lag_data2, self.lag(-res2, 1)], axis=1).dropna())
            else:
                X2 = sm.add_constant(self.lag(-res2, 1).dropna())
            
            y2 = regres2.dropna()
            min_len = min(len(y2), len(X2))
            y2, X2 = y2.iloc[-min_len:], X2.iloc[-min_len:]
            
            model2 = sm.OLS(y2, X2).fit()
            t2 = model2.tvalues.iloc[-1]
            
            test_results['t2'] = t2
            selected_lags_all.append(selected_lags2)
            
        except Exception as e:
            test_results['t2'] = 0
            selected_lags_all.append([])
        
        # Tests 3-7: Pair-wise frequency tests
        frequency_pairs = [3, 4, 5, 6, 7]  # Corresponding to filters 2-6
        f_stats = []
        
        for k in frequency_pairs:
            try:
                # Cointegrating regression with lag
                X_reg = pd.concat([filx[k-1], self.lag(filx[k-1], 1)], axis=1).dropna()
                X_reg.columns = ['x', 'x_lag']
                X_reg = sm.add_constant(X_reg)
                
                y_reg = fily[k-1].dropna()
                min_len = min(len(y_reg), len(X_reg))
                y_reg, X_reg = y_reg.iloc[-min_len:], X_reg.iloc[-min_len:]
                
                reg3 = sm.OLS(y_reg, X_reg).fit()
                res3 = pd.Series(reg3.resid, index=y_reg.index)
                
                w = self.eghl_filter(res3)
                regres = -w[7] / w[k-1]
                
                # Lag selection
                selected_lags3, lag_data3 = self.lag_selection_procedure(regres, pd.DataFrame())
                
                # Build regression
                base_regressors = pd.concat([self.lag(-res3, 2), self.lag(-res3, 1)], axis=1)
                base_regressors.columns = ['lag2', 'lag1']
                
                if len(lag_data3) > 0:
                    X3 = sm.add_constant(pd.concat([lag_data3, base_regressors], axis=1).dropna())
                else:
                    X3 = sm.add_constant(base_regressors.dropna())
                
                y3 = regres.dropna()
                min_len = min(len(y3), len(X3))
                y3, X3 = y3.iloc[-min_len:], X3.iloc[-min_len:]
                
                model3 = sm.OLS(y3, X3).fit()
                f_stat = model3.fvalue
                
                f_stats.append(f_stat)
                selected_lags_all.append(selected_lags3)
                
            except Exception as e:
                f_stats.append(0)
                selected_lags_all.append([])
        
        # Add F-statistics to results
        test_results.update({
            'f34': f_stats[0] if len(f_stats) > 0 else 0,
            'f56': f_stats[1] if len(f_stats) > 1 else 0,
            'f78': f_stats[2] if len(f_stats) > 2 else 0,
            'f910': f_stats[3] if len(f_stats) > 3 else 0,
            'f1112': f_stats[4] if len(f_stats) > 4 else 0
        })
        
        return test_results, selected_lags_all
    
    def simulate_distribution(self, x, y, lags_list, num_simulations=1000):
        """Simulate null distribution for cointegration test"""
        results = []
        
        for i in range(num_simulations):
            # Generate integrated series
            np.random.seed(i)
            sim_x = np.cumsum(np.random.normal(0, 1, len(x)))
            sim_y = np.cumsum(np.random.normal(0, 1, len(y)))
            
            # Apply seasonal integration
            for j in range(11):
                sim_x = np.cumsum(sim_x)
                sim_y = np.cumsum(sim_y)
            
            # Normalize
            sim_x = sim_x / (10 ** int(np.log10(np.abs(np.mean(sim_x)) + 1e-10)))
            sim_y = sim_y / (10 ** int(np.log10(np.abs(np.mean(sim_y)) + 1e-10)))
            sim_x = (sim_x - np.min(sim_x)) + 100
            sim_y = (sim_y - np.min(sim_y)) + 100
            
            # Run EGHL test
            try:
                test_stats, _ = self.eghl_cointegration_test(sim_x, sim_y)
                results.append([
                    test_stats['t1'], test_stats['t2'], test_stats['f34'],
                    test_stats['f56'], test_stats['f78'], test_stats['f910'], test_stats['f1112']
                ])
            except:
                results.append([0, 0, 0, 0, 0, 0, 0])
            
            if (i + 1) % 100 == 0:
                print(f"Simulation progress: {i + 1}/{num_simulations}")
        
        return np.array(results)
    
    def scoint(self, x, y, option="default", num_simulations=1000):
        """
        Main seasonal cointegration test function
        
        Parameters:
        -----------
        x, y : array-like
            Time series to test for cointegration
        option : str
            "default" to use pre-computed distributions, "manual" to simulate
        num_simulations : int
            Number of simulations for manual option
        
        Returns:
        --------
        dict : Dictionary containing p-values for X, Y integration and XY cointegration
        """
        x_series = pd.Series(x) if not isinstance(x, pd.Series) else x
        y_series = pd.Series(y) if not isinstance(y, pd.Series) else y
        
        # Test individual series for seasonal unit roots
        print("Testing X for seasonal unit roots...")
        teststat_x, lags_x = self.hegy_test(x_series)
        
        print("Testing Y for seasonal unit roots...")
        teststat_y, lags_y = self.hegy_test(y_series)
        
        print("Testing XY for cointegration...")
        teststat_xy, lags_xy = self.eghl_cointegration_test(x_series, y_series)
        
        if option == "manual":
            print("Simulating distributions...")
            # Simulate distributions
            distr_xy = self.simulate_distribution(x_series, y_series, lags_xy, num_simulations)
            distr_x = self.simulate_distribution(x_series, x_series, lags_x, num_simulations)
            distr_y = self.simulate_distribution(y_series, y_series, lags_y, num_simulations)
            
            # Calculate p-values using simulated distributions
            pval_x = [
                self.pval(teststat_x['t1'], distr_x[:, 0], two_sided=False, side="left"),
                self.pval(teststat_x['t2'], distr_x[:, 1], two_sided=False, side="left"),
                self.pval(teststat_x['f34'], distr_x[:, 2], two_sided=False, side="right"),
                self.pval(teststat_x['f56'], distr_x[:, 3], two_sided=False, side="right"),
                self.pval(teststat_x['f78'], distr_x[:, 4], two_sided=False, side="right"),
                self.pval(teststat_x['f910'], distr_x[:, 5], two_sided=False, side="right"),
                self.pval(teststat_x['f1112'], distr_x[:, 6], two_sided=False, side="right")
            ]
            
            pval_y = [
                self.pval(teststat_y['t1'], distr_y[:, 0], two_sided=False, side="left"),
                self.pval(teststat_y['t2'], distr_y[:, 1], two_sided=False, side="left"),
                self.pval(teststat_y['f34'], distr_y[:, 2], two_sided=False, side="right"),
                self.pval(teststat_y['f56'], distr_y[:, 3], two_sided=False, side="right"),
                self.pval(teststat_y['f78'], distr_y[:, 4], two_sided=False, side="right"),
                self.pval(teststat_y['f910'], distr_y[:, 5], two_sided=False, side="right"),
                self.pval(teststat_y['f1112'], distr_y[:, 6], two_sided=False, side="right")
            ]
            
            pval_xy = [
                self.pval(teststat_xy['t1'], distr_xy[:, 0], two_sided=False, side="left"),
                self.pval(teststat_xy['t2'], distr_xy[:, 1], two_sided=False, side="left"),
                self.pval(teststat_xy['f34'], distr_xy[:, 2], two_sided=False, side="right"),
                self.pval(teststat_xy['f56'], distr_xy[:, 3], two_sided=False, side="right"),
                self.pval(teststat_xy['f78'], distr_xy[:, 4], two_sided=False, side="right"),
                self.pval(teststat_xy['f910'], distr_xy[:, 5], two_sided=False, side="right"),
                self.pval(teststat_xy['f1112'], distr_xy[:, 6], two_sided=False, side="right")
            ]
            
        else:  # option == "default"
            # Use pre-computed distributions (interpolation)
            if self.distr is None or self.distr2 is None:
                print("Warning: Pre-computed distributions not available. Using approximate p-values.")
                # Use standard normal approximation (this is a simplification)
                pval_x = [min(stats.norm.cdf(teststat_x['t1']), 0.999) for _ in range(7)]
                pval_y = [min(stats.norm.cdf(teststat_y['t1']), 0.999) for _ in range(7)]
                pval_xy = [min(stats.norm.cdf(teststat_xy['t1']), 0.999) for _ in range(7)]
            else:
                # Use interpolation with pre-computed distributions
                pval_x = [
                    self.interpolation(teststat_x['t1'], self.distr, "left", 0),
                    self.interpolation(teststat_x['t2'], self.distr, "left", 1),
                    self.interpolation(teststat_x['f34'], self.distr, "right", 2),
                    self.interpolation(teststat_x['f56'], self.distr, "right", 3),
                    self.interpolation(teststat_x['f78'], self.distr, "right", 4),
                    self.interpolation(teststat_x['f910'], self.distr, "right", 5),
                    self.interpolation(teststat_x['f1112'], self.distr, "right", 6)
                ]
                
                pval_y = [
                    self.interpolation(teststat_y['t1'], self.distr, "left", 0),
                    self.interpolation(teststat_y['t2'], self.distr, "left", 1),
                    self.interpolation(teststat_y['f34'], self.distr, "right", 2),
                    self.interpolation(teststat_y['f56'], self.distr, "right", 3),
                    self.interpolation(teststat_y['f78'], self.distr, "right", 4),
                    self.interpolation(teststat_y['f910'], self.distr, "right", 5),
                    self.interpolation(teststat_y['f1112'], self.distr, "right", 6)
                ]
                
                pval_xy = [
                    self.interpolation(teststat_xy['t1'], self.distr2, "left", 0),
                    self.interpolation(teststat_xy['t2'], self.distr2, "left", 1),
                    self.interpolation(teststat_xy['f34'], self.distr2, "right", 2),
                    self.interpolation(teststat_xy['f56'], self.distr2, "right", 3),
                    self.interpolation(teststat_xy['f78'], self.distr2, "right", 4),
                    self.interpolation(teststat_xy['f910'], self.distr2, "right", 5),
                    self.interpolation(teststat_xy['f1112'], self.distr2, "right", 6)
                ]
        
        # Round p-values
        pval_x = [round(p, 3) for p in pval_x]
        pval_y = [round(p, 3) for p in pval_y]
        pval_xy = [round(p, 3) for p in pval_xy]
        
        results = {
            'X': pval_x,
            'Y': pval_y,
            'XY': pval_xy,
            'test_statistics': {
                'X': teststat_x,
                'Y': teststat_y,
                'XY': teststat_xy
            }
        }
        
        return CointegrationResults(results)

class CointegrationResults:
    """Class to handle and display cointegration test results"""
    
    def __init__(self, results_dict):
        self.results = results_dict
        self.frequencies = ["0", "6/12", "3/12 and 9/12", "4/12 and 8/12", 
                           "2/12 and 10/12", "5/12 and 7/12", "1/12 and 11/12"]
    
    def print_results(self):
        """Print integration and cointegration results"""
        print("\n" + "="*60)
        print("SEASONAL COINTEGRATION TEST RESULTS")
        print("="*60)
        
        # Integration results for X
        print("\nSeries X - Integration Results:")
        print("-" * 30)
        for i, freq in enumerate(self.frequencies):
            if self.results['X'][i] > 0.05:
                print(f"X is integrated at frequency {freq}")
        
        # Integration results for Y
        print("\nSeries Y - Integration Results:")
        print("-" * 30)
        for i, freq in enumerate(self.frequencies):
            if self.results['Y'][i] > 0.05:
                print(f"Y is integrated at frequency {freq}")
        
        # Cointegration results
        print("\nCointegration Results:")
        print("-" * 30)
        for i, freq in enumerate(self.frequencies):
            x_integrated = self.results['X'][i] > 0.05
            y_integrated = self.results['Y'][i] > 0.05
            cointegrated = self.results['XY'][i] <= 0.05
            
            if cointegrated and x_integrated and y_integrated:
                print(f"X and Y are cointegrated at frequency {freq}")
            elif cointegrated and (not x_integrated or not y_integrated):
                print(f"(X and Y are cointegrated at frequency {freq})")
    
    def summary(self):
        """Detailed summary table"""
        print("\n" + "="*80)
        print("DETAILED SUMMARY")
        print("="*80)
        
        # X series summary
        print("\nX Series:")
        print(f"{'Frequency':<15} {'p-value':<10} {'Integration':<15}")
        print("-" * 40)
        for i, freq in enumerate(self.frequencies):
            integration = "present" if self.results['X'][i] > 0.05 else "not present"
            print(f"{freq:<15} {self.results['X'][i]:<10} {integration:<15}")
        
        # Y series summary
        print("\nY Series:")
        print(f"{'Frequency':<15} {'p-value':<10} {'Integration':<15}")
        print("-" * 40)
        for i, freq in enumerate(self.frequencies):
            integration = "present" if self.results['Y'][i] > 0.05 else "not present"
            print(f"{freq:<15} {self.results['Y'][i]:<10} {integration:<15}")
        
        # Cointegration summary
        print("\nCointegration (XY):")
        print(f"{'Frequency':<15} {'p-value':<10} {'Cointegration':<15}")
        print("-" * 40)
        for i, freq in enumerate(self.frequencies):
            x_integrated = self.results['X'][i] > 0.05
            y_integrated = self.results['Y'][i] > 0.05
            cointegrated = self.results['XY'][i] <= 0.05
            
            if cointegrated and x_integrated and y_integrated:
                coint_status = "present"
            elif cointegrated:
                coint_status = "(present)"
            else:
                coint_status = "not present"
            
            print(f"{freq:<15} {self.results['XY'][i]:<10} {coint_status:<15}")
    
    def __repr__(self):
        return f"CointegrationResults(X={len(self.results['X'])} frequencies, Y={len(self.results['Y'])} frequencies)"

# Example usage and testing
if __name__ == "__main__":
    # Create example data (similar to the R script)
    np.random.seed(42)
    
    # Generate integrated time series with seasonal patterns
    n = 200
    x = np.cumsum(np.random.normal(0, 1, n))
    y = 0.8 * x + np.cumsum(np.random.normal(0, 0.5, n))
    
    # Add seasonal patterns
    seasonal_x = 5 * np.sin(2 * np.pi * np.arange(n) / 12)
    seasonal_y = 3 * np.cos(2 * np.pi * np.arange(n) / 12)
    
    x = x + seasonal_x
    y = y + seasonal_y
    
    # Create cointegration test object
    coint_test = SeasonalCointegration()
    
    # Run the test
    print("Running seasonal cointegration test...")
    print("This may take a few minutes...")
    
    results = coint_test.scoint(x, y, option="default")
    
    # Display results
    results.print_results()
    results.summary() 