library(RQuantLib)

AsianOption("arithmetic", "put", 49.0, 50.0, 0.0, .07, 0.5, 0.3)

AsianOption("arithmetic", "put", 49.0, 50.0, 0.0, .07, 0.5, 0.3, length=1.0)


AsianOption("geometric", "put", 49.0, 50.0, 0.0, .07, 0.5, 0.3)
