# Leveraged ETFs Project

This project models the performance of **Leveraged Exchange-Traded Funds (LETFs)** using a **complex stochastic model**. By simulating market conditions and incorporating leverage adjustments, the project aims to analyze the risks and long-term performance of LETFs, focusing on their compounding effect under different market conditions.

---

## Dataset Overview

The dataset used in this project consists of historical price data for leveraged ETFs:
- **Data Files**:
  - `dailyPrices.csv`: Contains daily price data for selected ETFs.
  - `weeklyPrices.csv`: Contains aggregated weekly prices for trend analysis.
- **Data Processing**: The raw data has been cleaned and preprocessed to ensure compatibility for the stochastic model.

---

## Methodology

This study utilizes a **stochastic differential equation (SDE)** framework to model the behavior of leveraged ETFs under varying market conditions. The approach leverages **Geometric Brownian Motion (GBM)** to capture the price dynamics of the underlying index and incorporates leverage factors to simulate the effects on ETF performance.

### 1. Stochastic Modeling

We use a **Geometric Brownian Motion (GBM)** model to simulate the price dynamics of the underlying index (e.g., DAX). The model is represented as follows:

$$
dS_t = \mu S_t \, dt + \sigma S_t \, dW_t
$$

Where:
- $\( S_t \)$ is the price of the underlying index.
- $\( \mu \)$ is the drift term (expected return).
- $\( \sigma \)$ is the volatility.
- $\( W_t \)$ is a Wiener process, representing random market fluctuations.

#### Leveraged ETFs Simulation:
For LETFs, the model is adjusted to account for leverage multipliers, where the return dynamics are given by:

$$
dL_t = L_t (\mu dt + \sigma dW_t)
$$

Where:
- $\( L_t \)$ represents the leveraged ETF price at time \( t \).
- $\( \mu \)$ and $\( \sigma \)$ are the drift and volatility of the underlying index.
- Leverage (e.g., 2x or 3x) is incorporated to scale the price dynamics of the ETF.

### 2. Numerical Solutions
We solved these SDEs using **Monte Carlo simulations** and MATLAB’s numerical solvers, running simulations over thousands of paths to estimate the ETF's performance under various market conditions. This approach allows us to capture the randomness in the price evolution and analyze different leverage factors.

### 3. Data Processing and Validation
- The raw data from `dailyPrices.csv` and `weeklyPrices.csv` was processed and transformed to match the required input format for the stochastic model.
- Simulated results were validated against real-world historical performance to check the accuracy of the model.

---

## Performance Evaluation

The model's performance was evaluated using the following key metrics:

### 1. Return Deviation
- We calculated the return deviation as the difference between the actual ETF return and the expected return based on the modeled leverage factor.
- The Monte Carlo simulation was run for 5,000 iterations over multiple time horizons (weekly, monthly, yearly) to assess how leverage and time impact return deviations.

### 2. Risk Assessment
- **Volatility**: The standard deviation of return deviations was analyzed to assess risk.
- **Drawdown Analysis**: We calculated the maximum drawdown across different market conditions (bull vs bear markets).
  
### 3. Simulation Results

The simulation results demonstrated that higher leverage amplifies the returns during favorable market conditions, but also increases the downside risk significantly during volatile or unfavorable markets.

---

## Key Visualizations

- **Leverage Impact on Return Distributions**:
  ![Leverage Impact](figures/leverage_impact.png)

- **Simulated Paths of Leveraged ETFs**:
  ![Simulated ETF Paths](figures/simulated_paths.png)

- **Risk Analysis (Standard Deviation of Return Deviations)**:
  ![Risk Analysis](figures/risk_analysis.png)

---
