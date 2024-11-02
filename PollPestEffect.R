# License of original code:

# MIT License
# Copyright (c) 2023 Charlie
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#   
#   The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.



# code adapted from https://github.com/CCNicholson/SEMLE/blob/main/SEMLE.R
# publication for original code: Lonsdorf, E. V., Nicholson, C. C., Rundlöf, M., & Williams, N. M. (2023). A spatially explicit model of landscape pesticide exposure to bees: development, exploration, and evaluation. Science of the Total Environment doi: 10.1016/j.scitotenv.2023.168146.



#### DEPENDENCIES ####
require(smoothie)
require(raster)

#####  Create test landscape ######
set.seed(2024)

m_m <-matrix( 
  c(rep(1,125000), 
    rep( c(rep (1,425),rep(2,150),rep(1,425)), 150),
    rep( c(rep (1,275),rep(2,150),rep (1,150),rep(2,150),rep(1,275)), 150),
    rep(1,150000),
    rep( c(rep (1,125),rep(2,150),rep (1,150),rep(2,150),rep (1,150),rep(2,150),rep(1,125)), 150),
    rep( c(rep (1,275),rep(2,150),rep (1,150),rep(2,150),rep(1,275)), 150),
    rep(1,125000)
  ),  1000, 1000)

LU.m <- raster(m_m, xmn=0, xmx=1000, ymn=0, ymx=1000)

#### PollPestEffect - Pollination model with pesticide effects ####

### Define land use scenario and masks for habitat and agricultural fields
  #agriculture = 1  #nature = 2  
  field_mask<-reclassify(LU.m,matrix(c(1,2,1,NA),nrow=2) ) 
  habitat_mask<-reclassify(LU.m, matrix(c(1,2,NA,1),nrow=2) ) 
  
### Define pesticide hazard on agricultural land uses and pollination capacity in habitat
  load.df <- data.frame(LU_Code = c( 1,2), mortality = c( 0.25, 0)) # mortality hazard 25%
  load.ras <- reclassify(LU.m, load.df)
  capa.df <- data.frame(LU_Code = c(1,2), capacity = c( 0, 1)) # initial habitat suitability 100%
  capa.ras <- reclassify(LU.m, capa.df)

### Make moving window object based on foraging range
  gamma_param = 150    # set example foraging range 
  
  {maxforage.dist = 2*gamma_param
  c.size = res(LU.m)[1]
  radius <- round(maxforage.dist/c.size) # matrix boundary for moving window
  
  # create matrices for weight arguments in focal() and setting moving window distance matrix
  weight.m = dist.m = matrix(data=1, nrow= radius * 2 + 1,   ncol= radius * 2 + 1) 
  focal.c = median(1:nrow(dist.m)) # focal.c is defined as row number of focal cell in the distance matrix
  
  # calculating distance of all cells from the focal cell in the nested matrix
  for(i in 1:nrow(dist.m)) {
    for (j in 1:ncol(dist.m)) {
      dist.m[i,j] = sqrt(((i-0.5)*c.size - (focal.c - 0.5)*c.size)^2 +
                           ((j-0.5)*c.size - (focal.c - 0.5)*c.size)^2 )
    }
  }
  
  # calculate effective distance matrix 
  effdist.m = effdist.m = exp(-dist.m / gamma_param)  # possible change of exponential decay formula  
  effdist.m[which(dist.m > (2*gamma_param) | dist.m > maxforage.dist)] <- 0 # set effective distance to 0 at distances larger than 2x foraging range
  }
  
### Perform FFT kernel smoothing
  FFT_matrix <- smoothie::kernel2dsmooth(x=raster::as.matrix(load.ras), K=effdist.m,   setup=T)
  load_dw <- smoothie::kernel2dsmooth(x=raster::as.matrix(load.ras),   W=FFT_matrix)
  pestiLoad <- raster::raster(load_dw, template=load.ras)  
  
  #convert all land uses into value 1 to clip output raster and calculate moving window normalization values
  mask_land <- raster::crop(raster::reclassify(LU.m, cbind(-1,  max(raster::values(LU.m), na.rm=T), 1)), extent(LU.m)) 
  #calculate moving window weight sums across landscape (will be the same value for most of the raster except along edges)
  mask_dw <- smoothie::kernel2dsmooth(x=as.matrix(mask_land),  W=FFT_matrix )
  window_sum <- raster::raster(mask_dw, template=mask_land)  
  
  #divide the summarized mortality values by the moving window sum & clip result to boundary of land use raster
  mortality <- (pestiLoad/window_sum)* mask_land  #plot(mortality, col = heat.colors(n = 100, rev = T))

### Reduce pollination capacity in habitat patches according to incurred mortality hazard
  mortality_habitat<-mask(mortality, habitat_mask) # Clip mortality to habitat patches
  polli_habi<- capa.ras - mortality_habitat  # subtract mortality from initial pollination capacity
  polli_habi<-reclassify(polli_habi,cbind(-Inf,0,0)) # in case the capacity is reduced to below 0, set it back to 0

### Calculate pollination from habitats with reduced pollination capacity
  # use FFT matrix, mask and window sum from before
  polli_dw <- smoothie::kernel2dsmooth(x=raster::as.matrix(polli_habi),  W=FFT_matrix)
  polliLoad <- raster::raster(polli_dw, template=polli_habi)
  polliLoad <- (polliLoad/window_sum)* mask_land
  polli_field<-mask(polliLoad, field_mask)   # clip calculated pollination to agricultural fields

