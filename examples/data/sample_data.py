"""
Sample Data Generation for Seasonal Cointegration Testing
Deutsche Bundesbank Internship Project 2018

This module provides functions to generate various types of sample data
for testing the seasonal cointegration procedures.
"""

import numpy as np
import pandas as pd
from typing import Tuple, Optional


def generate_tourism_data(n: int = 567, seed: Optional[int] = None) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generate sample tourism expenditure data similar to Austria/Spain data
    
    Parameters:
    -----------
    n : int
        Number of monthly observations (default: 567, ~47 years)
    seed : int, optional
        Random seed for reproducibility
        
    Returns:
    --------
    Tuple[np.ndarray, np.ndarray]
        Two time series representing tourism expenditures
    """
    if seed is not None:
        np.random.seed(seed)
    
    # Base integrated processes (unit roots)
    x_base = np.cumsum(np.random.normal(0, 1, n))
    
    # Cointegrating relationship with error correction
    error_correction = np.cumsum(np.random.normal(0, 0.3, n))
    y_base = 0.8 * x_base + error_correction
    
    # Monthly seasonal patterns (tourism peaks in summer)
    months = np.arange(n) % 12
    
    # Austria tourism: strong winter/summer peaks
    seasonal_x = 15 * np.sin(2 * np.pi * months / 12) + \
                 8 * np.cos(2 * np.pi * months / 12 + np.pi/3) + \
                 5 * np.sin(4 * np.pi * months / 12)  # Semi-annual component
    
    # Spain tourism: strong summer peaks, some winter
    seasonal_y = 12 * np.sin(2 * np.pi * months / 12 + np.pi/4) + \
                 6 * np.cos(2 * np.pi * months / 12) + \
                 3 * np.sin(4 * np.pi * months / 12 + np.pi/2)
    
    # Long-term trends (growing tourism industry)
    trend_x = 0.08 * np.arange(n) + 0.001 * np.arange(n)**1.1
    trend_y = 0.06 * np.arange(n) + 0.0008 * np.arange(n)**1.1
    
    # Combine components with base levels
    x = x_base + seasonal_x + trend_x + 150
    y = y_base + seasonal_y + trend_y + 120
    
    return x, y


def generate_production_data(n: int = 336, seed: Optional[int] = None) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generate sample production and sales data
    
    Parameters:
    -----------
    n : int
        Number of monthly observations (default: 336, ~28 years)
    seed : int, optional
        Random seed for reproducibility
        
    Returns:
    --------
    Tuple[np.ndarray, np.ndarray]
        Production and sales time series
    """
    if seed is not None:
        np.random.seed(seed)
    
    # Base integrated processes
    production_base = np.cumsum(np.random.normal(0, 0.8, n))
    
    # Sales follow production with some independent variation
    sales_innovation = np.cumsum(np.random.normal(0, 0.4, n))
    sales_base = 0.9 * production_base + sales_innovation
    
    # Industrial seasonal patterns
    months = np.arange(n) % 12
    
    # Production: lower in summer holidays, higher in spring/fall
    seasonal_prod = -8 * np.sin(2 * np.pi * months / 12 + np.pi/2) + \
                    4 * np.cos(2 * np.pi * months / 12) + \
                    3 * np.sin(4 * np.pi * months / 12 + np.pi/3)
    
    # Sales: peak before holidays, lower during holidays
    seasonal_sales = -6 * np.sin(2 * np.pi * months / 12 + np.pi/3) + \
                     5 * np.cos(2 * np.pi * months / 12 + np.pi/6) + \
                     2 * np.sin(4 * np.pi * months / 12)
    
    # Productivity trends
    trend_prod = 0.03 * np.arange(n)
    trend_sales = 0.025 * np.arange(n)
    
    # Combine components
    production = production_base + seasonal_prod + trend_prod + 100
    sales = sales_base + seasonal_sales + trend_sales + 95
    
    return production, sales


def generate_energy_data(n: int = 240, seed: Optional[int] = None) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generate sample energy consumption data (heating vs cooling)
    
    Parameters:
    -----------
    n : int
        Number of monthly observations (default: 240, ~20 years)
    seed : int, optional
        Random seed for reproducibility
        
    Returns:
    --------
    Tuple[np.ndarray, np.ndarray]
        Heating and cooling energy consumption time series
    """
    if seed is not None:
        np.random.seed(seed)
    
    # Base processes
    heating_base = np.cumsum(np.random.normal(0, 0.5, n))
    cooling_base = np.cumsum(np.random.normal(0, 0.6, n))
    
    # Strong seasonal patterns (opposite phases)
    months = np.arange(n) % 12
    
    # Heating: peak in winter (months 0, 1, 11)
    seasonal_heating = 20 * np.cos(2 * np.pi * months / 12) + \
                      8 * np.cos(4 * np.pi * months / 12) + \
                      5 * np.cos(6 * np.pi * months / 12)
    
    # Cooling: peak in summer (months 6, 7, 8)
    seasonal_cooling = -18 * np.cos(2 * np.pi * months / 12) + \
                      10 * np.sin(2 * np.pi * months / 12) + \
                      6 * np.cos(4 * np.pi * months / 12 + np.pi)
    
    # Mild trends (efficiency improvements)
    trend_heating = -0.01 * np.arange(n)  # Declining due to efficiency
    trend_cooling = 0.02 * np.arange(n)   # Increasing due to climate change
    
    # Combine components
    heating = heating_base + seasonal_heating + trend_heating + 80
    cooling = cooling_base + seasonal_cooling + trend_cooling + 40
    
    return heating, cooling


def generate_agricultural_data(n: int = 300, seed: Optional[int] = None) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generate sample agricultural price data (wheat vs corn)
    
    Parameters:
    -----------
    n : int
        Number of monthly observations (default: 300, ~25 years)  
    seed : int, optional
        Random seed for reproducibility
        
    Returns:
    --------
    Tuple[np.ndarray, np.ndarray]
        Wheat and corn price time series
    """
    if seed is not None:
        np.random.seed(seed)
    
    # Commodity prices often cointegrated but with different seasonal patterns
    wheat_base = np.cumsum(np.random.normal(0, 0.12, n))
    
    # Corn prices cointegrated with wheat (substitutes)
    cointegration_error = np.cumsum(np.random.normal(0, 0.08, n))
    corn_base = 0.85 * wheat_base + cointegration_error
    
    # Agricultural seasonal patterns
    months = np.arange(n) % 12
    
    # Wheat: harvest in summer (lower prices), planting costs in spring
    seasonal_wheat = -0.15 * np.sin(2 * np.pi * months / 12 + np.pi/3) + \
                     0.08 * np.cos(2 * np.pi * months / 12) + \
                     0.05 * np.sin(4 * np.pi * months / 12)
    
    # Corn: different harvest timing and storage patterns
    seasonal_corn = -0.12 * np.sin(2 * np.pi * months / 12 + np.pi/2) + \
                    0.06 * np.cos(2 * np.pi * months / 12 + np.pi/4) + \
                    0.04 * np.sin(4 * np.pi * months / 12 + np.pi/3)
    
    # Price inflation trends
    trend_wheat = 0.001 * np.arange(n)
    trend_corn = 0.0008 * np.arange(n)
    
    # Combine components (log prices)
    wheat_log = wheat_base + seasonal_wheat + trend_wheat + 2.5
    corn_log = corn_base + seasonal_corn + trend_corn + 2.3
    
    # Convert to price levels
    wheat_price = np.exp(wheat_log)
    corn_price = np.exp(corn_log)
    
    return wheat_price, corn_price


def create_sample_datasets() -> dict:
    """
    Create a collection of sample datasets for testing
    
    Returns:
    --------
    dict
        Dictionary containing different sample datasets
    """
    datasets = {}
    
    # Tourism data (like original R script)
    tourism_x, tourism_y = generate_tourism_data(567, seed=42)
    datasets['tourism'] = {
        'x': tourism_x,
        'y': tourism_y,
        'names': ['Austria_Tourism', 'Spain_Tourism'],
        'description': 'Tourism expenditure data (Austria vs Spain)'
    }
    
    # Production data  
    prod_x, prod_y = generate_production_data(336, seed=123)
    datasets['production'] = {
        'x': prod_x,
        'y': prod_y,
        'names': ['Production', 'Sales'],
        'description': 'Industrial production and sales data'
    }
    
    # Energy data
    energy_x, energy_y = generate_energy_data(240, seed=456)
    datasets['energy'] = {
        'x': energy_x,
        'y': energy_y,
        'names': ['Heating_Consumption', 'Cooling_Consumption'],
        'description': 'Seasonal energy consumption data'
    }
    
    # Agricultural data
    agri_x, agri_y = generate_agricultural_data(300, seed=789)
    datasets['agriculture'] = {
        'x': agri_x,
        'y': agri_y,
        'names': ['Wheat_Price', 'Corn_Price'],
        'description': 'Agricultural commodity prices'
    }
    
    return datasets


def save_datasets_to_csv(output_dir: str = "sample_datasets"):
    """
    Generate and save sample datasets to CSV files
    
    Parameters:
    -----------
    output_dir : str
        Directory to save the CSV files
    """
    import os
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Generate datasets
    datasets = create_sample_datasets()
    
    # Save each dataset
    for name, data in datasets.items():
        df = pd.DataFrame({
            data['names'][0]: data['x'],
            data['names'][1]: data['y']
        })
        
        filename = os.path.join(output_dir, f"{name}_data.csv")
        df.to_csv(filename, index=False)
        print(f"Saved {name} dataset to {filename}")
        print(f"  Description: {data['description']}")
        print(f"  Observations: {len(data['x'])}")
        print()


if __name__ == "__main__":
    # Example usage
    print("Sample Data Generation for Seasonal Cointegration Testing")
    print("=" * 60)
    
    # Generate tourism data (like original R script)
    x, y = generate_tourism_data(seed=42)
    print(f"Generated tourism data: {len(x)} observations")
    print(f"X (Austria) range: [{x.min():.2f}, {x.max():.2f}]")
    print(f"Y (Spain) range: [{y.min():.2f}, {y.max():.2f}]")
    
    print("\nSaving all sample datasets to CSV files...")
    save_datasets_to_csv()
    
    print("All datasets generated and saved successfully!") 