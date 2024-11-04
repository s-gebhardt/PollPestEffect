# Pollination model with pesticide effects

## Description
In the first part of our model, the pollinator flights and the subsequent exposure to pesticide hazard are simulated with a modelling approach presented by Lonsdorf et al. (2023). 
Using the published R code for their “Spatially Explicit Model of Landscape Exposure (SEMLE)”, we first calculate the exposure of pollinators at a given site through simulating their foraging behavior. Through a kernel smoothing calculation, the hazard within the foraging range is summarized based on the probability of a pollinator visiting a location, which exhibits a distance-decay function. While the original model in Lonsdorf et al. (2023) simulates landscape pesticide exposure as described, we treat the calculated numbers as the incurred mortality hazard. As such, the resulting values are substracted from the initial pollinator capacity (100%) in the respective habitat raster cell. The result is the remaining pollinator capacity, which we refer to as pollinator health in our results to reflect how much damages the pesticide exposure has caused. In our model, we further calculate pollination service as the visitation of pollinators to the agricultural fields. For this, our model utilizes again a kernel-smoothing approach to estimate the probability of pollinators arriving at an agricultural field location based on the reduced capacity or pollinator health in the habitats. 

## Citing Information
This model is based on the pesticide exposure modelling approach presented in 
Lonsdorf, E. V., Nicholson, C. C., Rundlöf, M., & Williams, N. M. (2023). A spatially explicit model of landscape pesticide exposure to bees: development, exploration, and evaluation. Science of the Total Environment doi: 10.1016/j.scitotenv.2023.168146.
