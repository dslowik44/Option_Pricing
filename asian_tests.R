library(RQuantLib)

AmericanOption(type="put", underlying=50, strike=50, dividendYield=0.0, 
riskFreeRate=0.1, maturity=5/12, volatility=0.4, engine="CrankNicolson")

EuropeanOption("put", 40.2, 43.0, 0.2, 0.07, 3/12, 0.32)

AsianOption("arithmetic", "put", 49.0, 50.0, 0.0, .07, 0.5, 0.3, length=0.5, fixings=1000)

AsianOption("geometric", "put", 49.0, 50.0, 0.0, .07, 0.5, 0.3, length=0.5, fixings=10000)

AsianOption("geometric", "put", 49.0, 50.0, 0.0, .07, 0.5, 0.3, length=0.05, fixings=10)

AsianOption("geometric", "put", 49.0, 50.0, 0.0, .07, 0.5, 0.3)
