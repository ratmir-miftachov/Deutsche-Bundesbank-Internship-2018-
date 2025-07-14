# Deutsche Bundesbank Internship Project 2018

## Extension of EGHL Seasonal Cointegration Test to Monthly Frequency

**A 12-week research project extending the Engle-Granger-Hylleberg-Lee seasonal cointegration test from quarterly to monthly frequency for enhanced analysis of seasonal economic relationships.**

<table>
<tr>
<td><img src="season1.png" alt="Seasonal Analysis 1" width="400"/></td>
<td><img src="season2.png" alt="Seasonal Analysis 2" width="400"/></td>
</tr>
</table>

---

## Abstract

This project presents an extension of the Engle-Granger-Hylleberg-Lee (EGHL) seasonal cointegration methodology from quarterly to monthly frequency. The standard EGHL approach, originally designed for quarterly data, tests for cointegration relationships at seasonal frequencies by decomposing the seasonal difference operator (1-L⁴) into its constituent roots. This work adapts the methodology to handle monthly data by extending the seasonal difference operator to (1-L¹²) and implementing the corresponding 12 seasonal frequencies.

The monthly extension requires significant theoretical and computational adaptations. The factorization of (1-L¹²) yields multiple complex seasonal frequencies that must be tested individually for unit roots and cointegration relationships. This research develops the complete statistical framework including: (1) HEGY-type tests for seasonal unit roots at all 12 monthly frequencies, (2) EGHL cointegration tests adapted for monthly seasonal patterns, (3) automatic lag selection procedures for optimal model specification, and (4) Monte Carlo simulation methods for generating critical values.

The methodology is implemented in both R and Python, providing researchers with robust tools for analyzing seasonal economic relationships at monthly frequency. Applications include tourism expenditure relationships, production-sales dynamics, energy consumption patterns, and agricultural commodity price linkages. The enhanced temporal resolution enables detection of more nuanced seasonal cointegration patterns that are masked in quarterly analysis, offering valuable insights for economic policy and forecasting.

---

## Repository Structure

```
Deutsche-Bundesbank-Internship-2018-/
├── README.md                          # Project overview (this file)
├── .gitignore                         # Git ignore patterns
│
├── docs/                              # Documentation
│   ├── Documentation.pdf             # Theoretical framework and methodology
│   └── Presentation.pdf              # Project presentation and results
│
├── r_implementation/                  # Original R Implementation
│   └── scoint_monthly.R              # Complete R package (667 lines)
│
├── python_implementation/             # Python Translation
│   ├── scoint_monthly.py             # Main Python implementation (785 lines)
│   └── requirements.txt              # Python dependencies
│
└── examples/                          # Usage Examples
    ├── python_example.py             # Complete Python usage example
    └── data/
        └── sample_data.py            # Sample data generation utilities
```

---

## Quick Start

**Python:**
```bash
cd python_implementation
pip install -r requirements.txt
cd ../examples && python python_example.py
```

**R:**
```bash
R -e "install.packages('smooth')" 
# Then: source("r_implementation/scoint_monthly.R")
```

---
---

## Methodology Overview

### Theoretical Framework

The implementation extends the **Engle-Granger-Hylleberg-Lee (EGHL)** seasonal cointegration test from quarterly (4 frequencies) to **monthly (12 frequencies)** data:

**Monthly Seasonal Difference Operator:**
```
(1 - L¹²) = (1 - L)(1 + L)(1 + L²)(1 + L + L²)(1 - L + L²)(1 + √3L + L²)(1 - √3L + L²)
```

**Tested Frequencies:**
1. **Frequency 0**: Long-run/annual relationship
2. **Frequency π (6/12)**: Semi-annual patterns
3. **Frequencies π/2, 3π/2**: Quarterly patterns  
4. **Frequencies 2π/3, 4π/3**: Tri-annual patterns
5. **Frequencies π/3, 5π/3**: Bi-monthly patterns
6. **Frequencies 5π/6, 7π/6**: Five-month patterns
7. **Frequencies π/6, 11π/6**: Monthly patterns

### Statistical Tests

#### HEGY Test (Individual Series)
- Tests for **seasonal unit roots** at each of the 12 frequencies
- Uses **13 linear filters** for monthly decomposition
- Automatic lag selection based on statistical significance

#### EGHL Test (Cointegration)
- Tests for **seasonal cointegration** between two series
- Uses **8 specialized filters** for cointegration analysis
- Monte Carlo simulation for critical value generation

### Key Features

- **Comprehensive Frequency Coverage**: Tests all 12 monthly seasonal frequencies
- **Automatic Lag Selection**: Statistical procedures for optimal model specification
- **Robust Implementation**: Handles near-singular matrices and numerical stability
- **Dual Language Support**: Complete implementations in both R and Python
- **Sample Data Generation**: Built-in utilities for testing and validation

---

## Applications

The monthly seasonal cointegration test is particularly valuable for:

### 🏖️ **Tourism Economics**
- Seasonal tourism demand relationships between countries/regions
- Monthly visitor flow patterns and expenditure analysis
- Tourism policy impact assessment at monthly frequency

### 🏭 **Industrial Analysis**
- Production-sales relationships with monthly seasonal components
- Supply chain dynamics and inventory management
- Monthly capacity utilization patterns

### ⚡ **Energy Economics**
- Monthly energy consumption patterns (heating vs cooling)
- Seasonal electricity demand relationships
- Renewable energy integration analysis

### 🌾 **Agricultural Economics**
- Commodity price relationships with monthly harvest cycles
- Seasonal arbitrage opportunities in agricultural markets
- Food security and price stability analysis

---

## Example Usage

### Python
```python
from python_implementation.scoint_monthly import SeasonalCointegration

# Initialize the test
coint_test = SeasonalCointegration()

# Load your monthly time series data
# x, y = your_monthly_data()

# Run seasonal cointegration test
results = coint_test.scoint(x, y, option="default")

# Display results
results.print_results()
results.summary()
```

### R
```r
source("r_implementation/scoint_monthly.R")

# Load your monthly time series data
# x <- ts(your_data_x, frequency=12)
# y <- ts(your_data_y, frequency=12)

# Run seasonal cointegration test
results <- scoint(x, y, option="default")
print(results)
```

---


