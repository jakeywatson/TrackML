# ML Final Project: Particle Tracking
# Jake Watson,20/06/2020
# UPC

# Libraries
require(data.table)
require(bit64)
require(dbscan)
require(doParallel)
require(dplyr)
require(rBayesianOptimization)
require(plotly)
require(DescTools)
require(viridis)

path='../data/train_1/'

###############################################################################################################
# ANALYSE THE HELIX PARAMETERS FOR THE TRUTH HITS FOR ONE EVENT
###############################################################################################################
# Approach taken from CPMP - Kaggle LB #12
helix_params =  data.frame()
for (ii in 10:99){
  print(ii)
  truth = data.table(read.csv(paste0(path, "event0000010", ii, "-truth.csv"), header=TRUE, sep=","))
  hits = data.table(read.csv(paste0(path, "event0000010", ii, "-hits.csv"), header=TRUE, sep=","))
  cells = data.table(read.csv(paste0(path, "event0000010", ii, "-cells.csv"), header=TRUE, sep=","))
  particles = data.table(read.csv(paste0(path, "event0000010", ii, "-particles.csv"), header=TRUE, sep=","))
  
  data = hits
  data = merge(data, truth, by="hit_id")
  data = merge(data, particles, by="particle_id")
  data[, rv:=sqrt(vx*vx + vy*vy)]
  data = data[rv <= 1 & vz <= 50 & vz >=-50,]
  data = data[weight > 0]
  data[, event_id:=1000+ii]
  data[, pt:=sqrt(px*px+py*py)]
  data[, alpha:=exp(-8.115 - log(pt))]

  if (ii == 0) {
    helix_params = as.data.frame(data[, .(event_id, particle_id, alpha, vz)])
  } else {
    helix_params = rbind(helix_params, as.data.frame(data[, .(event_id, particle_id, alpha, vz)]))
  }
}

df = helix_params %>% group_by(event_id, particle_id) %>%   filter(row_number()==1)
df <- na.omit(df)

write.csv(df, "../data/helix_params.csv", row.names=FALSE)
df <- read.csv("../data/helix_params.csv")


###############################################################################################################
# DATA EXPLORATION AND VISUALISATION
###############################################################################################################


#Read single event
hits = data.table(read.csv(paste0(path, "event000001000-hits.csv"), header=TRUE, sep=","))
hits[, r:=sqrt(x*x+y*y+z*z)]
hits[, rt:=sqrt(x*x+y*y)]
hits[,x1:=x/r]
hits[,y1:=y/r]
hits[,z1:=z/rt]
summary(hits)
hits

cells = data.table(read.csv(paste0(path, "event000001000-cells.csv"), header=TRUE, sep=","))
summary(cells)

particles = data.table(read.csv(paste0(path, "event000001000-particles.csv"), header=TRUE, sep=","))
summary(particles)

truth=data.table(read.csv(paste0(path,'event000001000-truth.csv'),header = TRUE, sep=","))
summary(truth)

detectors=data.table(read.csv("../data/detectors.csv",header = TRUE, sep=","))
summary(detectors)

PlotFdist(hits$x, "Distribution of X co-ord position")
PlotFdist(hits$y, "Distribution of Y co-ord position")
PlotFdist(hits$z, "Distribution of Z co-ord position")

PlotFdist(cells$value, "Distribution of charge deposition energy")

PlotFdist(particles$nhits)

PlotFdist(truth$tx, "Distribution of X-position of interactions")
PlotFdist(truth$ty, "Distribution of X-position of interactions")
PlotFdist(truth$tz, "Distribution of X-position of interactions")

PlotFdist(truth$tpx, "Distribution of X-momentum")
PlotFdist(truth$tpy, "Distribution of Y-momentum")
PlotFdist(truth$tpz, "Distribution of Z-momentum")

PlotFdist(truth$weight, "Distribution of hit weights")

hist(hits$r,breaks=500)
hist(hits$rt,breaks=500)
plot(hits[,.(x,x1)],pch='.',col=rgb(abs(hits$z)/3100,0,0,.01))
plot(hits[,.(z,z1)],pch='.',col=rgb(abs(hits$z)/3100,0,0,.01))
plot(hits[,.(x,r)],pch='.',col=rgb(abs(hits$z)/3100,0,0,.01))

plot(hits$x, hits$y)
plot(hits$z, hits$y)

particles %>% 
  plot_ly(x = ~vx, y = ~vy, color = ~vz,
          frame = ~nhits)

particles %>%
  arrange(desc(nhits)) %>%
  mutate(nhits = as.factor(nhits)) %>%
  plot_ly(x = ~vx, y = ~vy, z = ~vz, color = ~nhits, colors = viridis(19)) %>%
  add_markers(opacity = 0.8) %>%
  layout(title = "Particles by vertex and number of hits",
         annotations=list(yref='paper',xref="paper",y=1.05,x=1.1, text="Number of hits",showarrow=F), 
         scene = list(xaxis = list(title = 'vx'),
                      yaxis = list(title = 'vy'),
                      zaxis = list(title = 'vz')))


particles[, p:=sqrt(px*px +py*py +pz*pz)]
particles %>% plot_ly(x = ~px, y = ~py, z=~pz, color=~nhits, type = 'scatter', mode = 'markers')

particles %>%
  arrange(desc(nhits)) %>%
  mutate(nhits = as.factor(nhits)) %>%
  plot_ly(x = ~px, y = ~py, z = ~pz, color = ~nhits, colors = viridis(19)) %>%
  add_markers(opacity = 0.8) %>%
  layout(title = "Particles by momentum and number of hits",
         annotations=list(yref='paper',xref="paper",y=1.05,x=1.1, text="Number of hits",showarrow=F), 
         scene = list(xaxis = list(title = 'px'),
                      yaxis = list(title = 'py'),
                      zaxis = list(title = 'pz')))

# Visualise true hits in both representations

hits[, r:=sqrt(x*x+y*y+z*z)]
hits[, rt:=sqrt(x*x+y*y)]
hits[, a0:=atan2(y, x)]
hits[, ydivr:=y/r]
hits[, xdivr:=x/r]
hits[, dtheta:=(rt+0.000005*rt^2)/1000*(100/2)/180*pi]
hits[, theta_ratio:= 1 - (abs(z + 200) / 6000)**2.4 + 0.005]
hits[, a1:=a0+dtheta*theta_ratio]
hits[, zdivrt:=z/rt]
hits[, zdivr:=z/r]
hits[,sina1:=sin(a1)]
hits[,cosa1:=cos(a1)]

data = hits
data = merge(data, truth, by="hit_id")
data = merge(data, particles, by="particle_id")
first_particles = unique(data$particle_id)[1:8]
colours <- row_number(first_particles)
combined <- cbind(first_particles, colours)
colnames(combined) <- c("particle_id", "colour")

data = merge(data, combined, by="particle_id")
tracks = data[which(particle_id %in% first_particles), .(colour, x, y, z, zdivr, sina1, cosa1)]
scaled <- data.table(scale(tracks[, .(zdivr, sina1, cosa1)]))
tracks$cosa1 <- scaled$cosa1
tracks$sina1 <- scaled$sina1
tracks$zdivr <- scaled$zdivr
scaled <- data.table(scale(tracks[, .(x, y, z)]))
tracks$x <- scaled$x
tracks$y <- scaled$y
tracks$z <- scaled$z

zpt <- tracks %>% plot_ly(x= ~zdivr*0.25, y= ~cosa1*3, z= ~sina1*3, color= ~colour)
xyz <- tracks %>% plot_ly(x= ~x, y= ~y, z= ~z, color= ~colour)


###############################################################################################################
# TRANSFORM AND CLUSTER
###############################################################################################################

# PARALLEL LOOP
nrCores <- detectCores()
c1 <- makeCluster(nrCores)
registerDoParallel(c1)
clusterEvalQ(c1, c(library(data.table), library(dbscan)))

scores <- foreach (event=0:9, .combine=rbind) %dopar% {
  n_event = 1001 + event
  path = paste0("C:/Users/jaked/OneDrive/Documents/UPC/Semester 2/ML/project/data/train_1/event00000", n_event)
  truth = data.table(read.csv(paste0(path, "-truth.csv"), header=TRUE, sep=","))
  
  clustered <- cluster_transform3(path, 
                                  w1=3.0, w2=3.0, w3=1.0, w4=0.2521, w5=0.0211, w6=0.0474,
                                  n_iters = 190, 
                                  ep = 0.0088)
  
  labelled <- clustered[,.(n_event, hit_id, track_id)]
  s = score(labelled, truth)
  return(s)
}
scores
# SEQUENTIAL LOOP
scores = data.frame(0:0)
for (event in 0:0) {
  n_event = 1000
  path = paste0("C:/Users/jaked/OneDrive/Documents/UPC/Semester 2/ML/project/data/train_1/event00000", n_event)
  truth = data.table(read.csv(paste0(path, "-truth.csv"), header=TRUE, sep=","))
  print("Clustering")
  clustered <- cluster_transform3(path, 
                w1=3.0, w2=3.0, w3=1.0, w4=0.2521, w5=0.0211, w6=0.0474,
                n_iters = 190, 
                ep = 0.0088)
  
  labelled <- clustered[,.(n_event, hit_id, track_id)]
  s = score(labelled, truth)
  print(paste0("Event score: ", s))
  scores[event+1, ] = data.frame(s)
}
scores

###############################################################################################################
# OPTIMISATION OF PARAMETERS
###############################################################################################################


single4opt <- function(w1, w2, w3, w4, w5, w6, n_iters, ep) {
  n_event = 1000
  path = paste0("C:/Users/jaked/OneDrive/Documents/UPC/Semester 2/ML/project/data/train_1/event00000", n_event)
  truth = data.table(read.csv(paste0(path, "-truth.csv"), header=TRUE, sep=","))
  
  clustered <- cluster_transform3(path, w1, w2, w3, w4, w5, w6, n_iters, ep)
  labelled <- clustered[,.(n_event, hit_id, track_id)]
  s = score(labelled, truth)
  list(Score=s, Pred=0)
}

OPT <- BayesianOptimization(single4opt,
  bounds = list(w1 = c(1.0, 3.0), # sina1 - 3.0
                w2 = c(1.0, 3.0), # cosa1 - 3.0
                w3 = c(0.5, 1.0), # zdivrt - 1.0
                w4 = c(0.1, 0.5), # zdivr - 0.2521
                w5 = c(0.01, 0.05), # xdivr - 0.0211
                w6 = c(0.01, 0.05), # ydivr - 0.0474
                n_iters = c(150L, 200L), # iters - 190
                ep = c(0.001, 0.01)), # dbscan epsilon - 0.0088 
  init_points = 3, 
  n_iter = 20,
  acq = "ucb",
  kappa = 2.576,
  eps = 0.0,
  verbose = TRUE)


cluster_transform <- function(path, n){
  truth = data.table(read.csv(paste0(path, "event00000", n, "-truth.csv"), header=TRUE, sep=","))
  
  hits = data.table(read.csv(paste0(path, "event00000", n, "-hits.csv"), header=TRUE, sep=","))
  hits[, r:=sqrt(x*x+y*y+z*z)]
  hits[, rt:=sqrt(x*x+y*y)]
  hits[,x1:=x/r]
  hits[,y1:=y/r]
  hits[,z1:=z/rt]
  scaled = scale((hits[,.(x1, y1, z1)])) # Scale the transformed coords: subtract mean and divide by s.d.
  clusters = dbscan(scaled, eps = 0.007, minPts = 1)
  
  labelled = data.table(hit_id=hits$hit_id, track_id=clusters$cluster)
  s = score(labelled, truth)
  
  return(s)
}


# Transforms data into helix parameters, then performs a clustering on the transformed data.
# Helix params are unknown, so check multiple values, and choose best clustering.
# Transforms:
# angle a1: atan2(y, x), tangential angle in x-y, which is robustified by separating into sin and cos a0 (avoids singularities).
# z1: z/r, speed in z against total speed.
# z2: z/rt, slope of the helix (speed in z against speed in x-y)
# For each clustering, params z0 and R for the helix are unknown, so test many values.
# Best clusters are chosen by length. Longer is better.
cluster_transform2 <- function(path, n){
  
  # Read truth data
  truth = data.table(read.csv(paste0(path, "event00000", n, "-truth.csv"), header=TRUE, sep=","))
  
  # Read hits data
  # Transform initial hit coords to helix parameters, for use by dbscan
  # Helix params are invariant for points on helix, so dbscan will cluster them.
  hits = data.table(read.csv(paste0(path, "event00000", n, "-hits.csv"), header=TRUE, sep=","))
  hits[, r:=sqrt(x*x+y*y+z*z)]
  hits[, rt:=sqrt(x*x+y*y)]
  hits[, a0:=atan2(y, x)]
  hits[, x1:=x/r]
  hits[, y1:=y/r]
  hits[, z1:=z/rt]
  hits[, z2:=z/r]
  
  # We don't know which helix a point belongs to
  # So we try multiple helix radii and initial z-pos
  dz0    <- -0.00070
  stepdz <-  0.00001
  stepeps<-  0.000005
  mm     <-  1
  
  # Loop for n_iters
  for (ii in 0:100) {
    
    #Increment params
    mm <- mm*(-1)
    dz <- mm*(dz0 + ii*stepdz)
    hits[,a1:=a0+dz*z*sign(z)]
    
    # Use robust params
    hits[,sina1:=sin(a1)]
    hits[,cosa1:=cos(a1)]
    
    # Scale data
    scaled=scale(hits[,.(sina1,cosa1,z1,z2)])
    cx <- c(1, 1, 0.4, 0.4)
    for (jj in 1:ncol(scaled)) scaled[,jj] <- scaled[,jj]*cx[jj]
    
    # Do clustering
    clusts=dbscan(scaled,eps=0.0035+ii*stepeps,minPts = 1)
    
    # If initial run, just label each hit to a track
    # and count track populations
    if (ii==0) {
      hits[,s1:=clusts$cluster]
      hits[,N1:=.N, by=s1]
      
    }else{
      # If not initial run, label each hit to a new_track
      # count each new_track population
      hits[,s2:=clusts$cluster]
      hits[,N2:=.N, by=s2]
      
      # Find unique large number
      maxs1 <- max(hits$s1)
      
      # If the old_track population is less than the new_track population,
      # and the new_track is less than 20, then relabel hit to new_track using unique id
      hits[,s1:=ifelse(N2>N1 & N2<20,s2+maxs1,s1)]
      hits[,s1:=as.integer(as.factor(s1))]
      
      # Count new population
      hits[,N1:=.N, by=s1]
    }
  }
  # Merge hits with labels
  labelled = data.table(event_id=1001+event, hit_id=hits$hit_id, track_id=hits$s1)
  s = score(labelled, truth)
  
  return(s)
}


###############################################################################################################
# MAIN CLUSTERING ALGORITHM
###############################################################################################################
# Transforms XYZ hit data into invariant helix parameters, then performs a clustering on the transformed data.
# Multiple clusterings are performed for multiple helix parameters, favouring the parameters that produce the 
# largest clusters.
#
# Transforms:
# angle a1: unrolling angle for helix, robustified by separating into sin and cos a0 (avoids singularities).
# z1: z/r, speed in z against total speed.
# z2: z/rt, slope of the helix (speed in z against speed in x-y)
# zdivr: z/r, speed in z against total speed.
# zdivrt: z/rt, slope of the helix (speed in z against speed in x-y)
# xdivr: x/r, x-position of the hit divided by the z-position
# ydivr: y/r, y-position of the hit divided by the z-position

# Inputs:
# path: full path to input event file
# w1 - w6: weights for the features, 
# n_iters: number of iterations to try clustering
# ep: epsilon parameter for dbscan (radius of clustering)

# Outputs:
# labelled hits: adds columns to hits with cluster id and size of the cluster.

cluster_transform3 <- function(path, w1, w2, w3, w4, w5, w6, n_iters, ep){

  # Read hits data
  # Transform initial hit coords to helix parameters, for use by dbscan
  # Helix params are invariant for points on helix, so dbscan will cluster them.
  hits = data.table(read.csv(paste0(path, "-hits.csv"), header=TRUE, sep=","))
  hits[, r:=sqrt(x*x+y*y+z*z)]
  hits[, rt:=sqrt(x*x+y*y)]
  hits[, a0:=atan2(y, x)]
  hits[, ydivr:=y/r]
  hits[, xdivr:=x/r]

  mm     <-  1
  
  # Loop for n_iters
  for (ii in 0:n_iters) {

    # Helix goes either direction
    mm <- mm*(-1)

    # Calculate transformed coords
    hits[, dtheta:=mm*(rt+0.000005*rt^2)/1000*(ii/2)/180*pi]
    hits[, theta_ratio:= 1 - (abs(z + 200) / 6000)**2.4 + 0.005]
    hits[, a1:=a0+dtheta*theta_ratio]
    hits[, zdivrt:=z/rt]
    hits[, zdivr:=z/r]
    
    # Use robust params
    hits[,sina1:=sin(a1)]
    hits[,cosa1:=cos(a1)]
    
    # Scale data
    scaled=scale(hits[,.(sina1,cosa1,zdivrt,zdivr, xdivr, ydivr)])
    weights <- c(w1, w2, w3, w4, w5, w6)
    for (jj in 1:ncol(scaled)) scaled[,jj] <- scaled[,jj]*weights[jj]
    
    # Do clustering
    clusts=dbscan(scaled,eps=ep,minPts = 1)
    
    # If initial run, just label each hit to a track
    # and count track populations
    if (ii==0) {
      hits[,track_id:=clusts$cluster]
      hits[,N1:=.N, by=track_id]
    }else{
      # If not initial run, label each hit to a new_track
      # count each new_track population
      hits[,s2:=clusts$cluster]
      hits[,N2:=.N, by=s2]
      
      # Find unique large number
      maxs1 <- max(hits$track_id)
      
      # If the old_track population is less than the new_track population,
      # and the new_track is less than 20, then relabel hit to new_track using unique id
      hits[,track_id:=ifelse(N2>N1 & N2<20,s2+maxs1,track_id)]
      hits[,track_id:=as.integer(as.factor(track_id))]
      
      # Count new population
      hits[,N1:=.N, by=track_id]
    }
  }
  return(hits)
}



###############################################################################################################
# CLUSTER POSTPROCESSING
###############################################################################################################
# Attempts to remove poor clusters. Loops over each cluster, removes those of size above or below thresholds.
# Uses a least-squares regression to determine how well the cluster fits to a cylindrical shape.
# 
# Inputs: 
# hits: clustered hits data
# min: minimum size of clusters
# max: maximum size of clusters
#
# Outputs:
# hits: clustered hits data, with spurious clusters removed

find_norms <- function(hits, min, max) {
  hits[, track_id:=ifelse(N1 < min | N1 > max, 0, track_id)]
  hits[, nm:=0]
  
  clusters <- unique(hits$track_id)

  counter = 0
  for (cluster_id in clusters) {
    print(paste0("Analysing cluster ", counter, " out of ", length(clusters)))
    cluster <- hits[which(track_id==cluster_id),]
    size <- nrow(cluster)
    counter <- counter + 1
    
    if (cluster_id==0) next

    xyz <- as.matrix(cluster[, .(x, y, z)])
    Z <- matrix(nrow=size, ncol=10)
    x <- xyz[,1]; y <- xyz[, 2]; z <- xyz[, 3];

    Z[, 1] <- x*x
    Z[, 2] <- 2*x*y
    Z[, 3] <- 2*x*z
    Z[, 4] <- 2*x    
    Z[, 5] <- y*y
    Z[, 6] <- 2*y*z
    Z[, 7] <- 2*y
    Z[, 8] <- z*z
    Z[, 9] <- 2*z
    Z[, 10] <- 1
    
    sv <- svd(Z, nu=nrow(Z), nv=ncol(Z)) 
    min_index = which.min(sv$d)
    T <- as.matrix(sv$v[min_index,])
    fit <- norm(Z %*% T, type="2")**2
    hits[, nm:=ifelse(track_id==cluster_id, fit, nm)]
  }
}

# Find unique cluster size 
analyse_clusters <- function(hits) {
  
  clusters <- hits %>% group_by(track_id)
  return(group_size(clusters))
  
}


###############################################################################################################
# SCORING FUNCTION
###############################################################################################################
# Gives the accuracy of the clustering algorithm.
# It is the intersection between the reconstructed tracks and
# the ground truth tracks. For each event, this is normalized to one.
# The score is then the average of all events computed.

# Only those clusters with more than 50% of its hits correctly matched to a single particle
# have score above zero. Additionally, the cluster must contain more than 50% of the matching 
# particle's hits. This is known as the double majority.
# The score of the track that matches this rule is then the sum of the weights of the points
# the intersection between the track and the particle.
# This per-track score is then summed over the event, normalized to one.

# Inputs: 
# labelled: clustered hits
# truth: ground truth hits
#
# Outputs:
# score: value from 0 to 1 indicating the accuracy of the clustering according to the above metric.

score <- function(labelled, truth) {
  tmp = merge(labelled, truth[,.(hit_id, particle_id, weight)])
  tmp[, Np:=.N, by=particle_id]
  tmp[, Nt:=.N, by=track_id]
  tmp[, Ntp:=.N, by=list(track_id, particle_id)]
  tmp[, r1:=Ntp/Nt]
  tmp[, r2:=Ntp/Np]
  sum(tmp[r1>0.5 & r2 > 0.5, weight])
}
