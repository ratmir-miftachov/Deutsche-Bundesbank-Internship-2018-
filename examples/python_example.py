#!/usr/bin/env python3
"""
Example usage of the Seasonal Cointegration Test (Python version)
Deutsche Bundesbank Internship Project 2018

This script demonstrates how to use the scoint_monthly.py module 
to test for seasonal cointegration between two time series.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os

# Add the python_implementation directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'python_implementation'))

from scoint_monthly import SeasonalCointegration

def generate_sample_data():
    """Generate sample tourism/economic data similar to the original R script"""
    np.random.seed(42)
    n = 567  # Similar to the Austrian/Spanish tourism data length
    
    # Generate integrated time series with seasonal patterns
    # Base integrated series
    x_base = np.cumsum(np.random.normal(0, 1, n))
    y_base = 0.7 * x_base + np.cumsum(np.random.normal(0, 0.5, n))
    
    # Add monthly seasonal patterns (tourism seasonality)
    months = np.arange(n) % 12
    seasonal_x = 10 * np.sin(2 * np.pi * months / 12) + 5 * np.cos(4 * np.pi * months / 12)
    seasonal_y = 8 * np.sin(2 * np.pi * months / 12 + np.pi/4) + 3 * np.cos(4 * np.pi * months / 12)
    
    # Add trend
    trend_x = 0.05 * np.arange(n)
    trend_y = 0.03 * np.arange(n)
    
    # Combine components
    x = x_base + seasonal_x + trend_x + 100  # Add base level similar to R script
    y = y_base + seasonal_y + trend_y + 100
    
    return x, y

def plot_data(x, y):
    """Plot the time series data"""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))
    
    ax1.plot(x, label='Series X (e.g., Austria Tourism)', color='blue')
    ax1.set_title('Time Series X')
    ax1.set_ylabel('Values')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    ax2.plot(y, label='Series Y (e.g., Spain Tourism)', color='red')
    ax2.set_title('Time Series Y')
    ax2.set_xlabel('Time (months)')
    ax2.set_ylabel('Values')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('sample_data.png', dpi=150, bbox_inches='tight')
    plt.show()

def main():
    """Main example execution"""
    print("="*60)
    print("SEASONAL COINTEGRATION TEST - PYTHON VERSION")
    print("Deutsche Bundesbank Internship Project 2018")
    print("="*60)
    
    # Generate sample data
    print("\n1. Generating sample tourism data...")
    x, y = generate_sample_data()
    print(f"   Generated {len(x)} monthly observations")
    
    # Plot the data
    print("\n2. Plotting time series data...")
    try:
        plot_data(x, y)
        print("   Plot saved as 'sample_data.png'")
    except:
        print("   Warning: Could not create plot (matplotlib not available)")
    
    # Initialize the cointegration test
    print("\n3. Initializing seasonal cointegration test...")
    coint_test = SeasonalCointegration()
    
    # Run the test with default option (faster)
    print("\n4. Running seasonal cointegration test (default mode)...")
    print("   This tests for integration and cointegration at 7 different frequencies:")
    print("   - Frequency 0: Annual/long-run relationship")
    print("   - Frequency 6/12: Semi-annual patterns")  
    print("   - Frequencies 1/12-5/12: Monthly seasonal patterns")
    
    results = coint_test.scoint(x, y, option="default")
    
    # Display results
    print("\n5. RESULTS:")
    results.print_results()
    
    print("\n6. DETAILED SUMMARY:")
    results.summary()
    
    # Interpretation
    print("\n" + "="*60)
    print("INTERPRETATION GUIDE:")
    print("="*60)
    print("Integration Test:")
    print("- If p-value > 0.05: Series has a unit root at that frequency")
    print("- If p-value ≤ 0.05: Series is stationary at that frequency")
    print()
    print("Cointegration Test:")
    print("- 'present': Both series integrated at frequency, cointegration detected")
    print("- '(present)': Cointegration detected but integration pattern unclear")  
    print("- 'not present': No cointegration at that frequency")
    print()
    print("Economic Interpretation:")
    print("- Frequency 0: Long-run equilibrium relationship")
    print("- Frequency 6/12: Semi-annual seasonal cointegration")
    print("- Other frequencies: Monthly seasonal relationships")
    
    # Optional: Run with manual simulation (slower but more accurate)
    run_manual = input("\nRun manual simulation for more accurate p-values? (y/n): ")
    if run_manual.lower() == 'y':
        print("\n7. Running manual simulation (this may take several minutes)...")
        results_manual = coint_test.scoint(x, y, option="manual", num_simulations=500)
        print("\nManual simulation results:")
        results_manual.summary()

if __name__ == "__main__":
    main() 