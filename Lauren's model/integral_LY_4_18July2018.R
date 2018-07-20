require(deSolve)
library(ggplot2)
library(dplyr)

##Blue rockfish
##derivative of how population changes in time based on von bertalanffy growth and mortality.  i.e., this reflects the size distribution at equilibrium under constant recruitment. Look at stand distribution WFC222 Lecture 4.
eqn=function(l,y,parms) {
  N1=y[1]
  k=parms[1] 
  M=parms[2] 
  Linf=parms[3] 
  
  dN=numeric(1); # Vector to hold the integrands
  dN[1]=-N1*((M-k)/(k*(Linf-l)))
  return(list(dN))
}

##If need initial abundance N1 to match the mockdata
# F0.05 <- read.csv("/R_stuff/R_Research/MPAs/Blue_rockfish_sizes_F0.05.csv")
# N1 <- sum(F0.05$Year1)
# y=c(N1=20); 
##

##gopher rockfish parms
#y=c(N1=20);
#parms=c(k=0.2256, M=0.2, Linf=34.1);
#l=seq(1,parms[3]-0.01,length=10000); #sequence of lengths, previously 10000; another option: length=round(parms[3]

##blue rockfish parms
y=c(N1=20);
parms=c(k=0.172, M=0.14, Linf=38.15);
l=seq(1,parms[3]-0.01,length=10000); #sequence of lengths, previously 10000; another option: length=round(parms[3]

#Integrate the derivative to get abundance N as a function of length l
out=ode(y,l,eqn,parms) #out[,2] are the abundances at sizes 1-34
t0=-0.5
max_age=35 #was 35; try to lower this to 10, 15, 20 to reduce # of cohorts
bin_w <- 1

#Plot the integral of abundance N as a function of length l
plot(out, col=1, xlab = "Length", ylab = "Abundance (N)", main = "Blue rockfish", type = "p")

colnames(out) <- c("length", "N")
abundance_at_size <- data.frame(out)

#Convert ages to mean size-at-age based on von bertalanffy eq.
#t0 = -0.5, so x0 = 3.6375
Linf <- parms[3]
k <- parms[1]

size_at_age <- vector()
for (j in 0:max_age){
size_at_age[j+1] <- Linf*(1-exp(-k*(j-t0)))
}

mean_size <- vector()
lower_size <- vector()
upper_size <- vector()
sd <- vector()

for (s in 1:length(size_at_age)){
  mean_size[s] <- size_at_age[s]
  sd[s] <- (mean_size[s])*.1
  lower_size[s] <- mean_size[s]-(3*sd[s])
  upper_size[s] <- mean_size[s]+(3*sd[s])
  #upper_size[s] <- ifelse(upper_size[s] > Linf, 34, upper_size[s])
}

# Take the VB size-at-age, and add variability around each mean size-at-age equal to 3 standard deviations
Blue_sizes <- data.frame(mean_size = size_at_age, sd = sd, lower_size = lower_size, upper_size = upper_size)

Blue_sizes <- Blue_sizes %>% mutate_all(.,funs(round(.,digits = 2))) #mutate_each

#abundance_at_size <- data.frame(out)
# Round VB size_at_age to 2 decimal places to make it easier to match with sizes from abundance as a function of length from the ODE
size_at_age <- round(size_at_age, digits = 2)
lower_size <- round(lower_size, digits = 2)
upper_size <-round(upper_size, digits = 2)

#Determine how many individuals are in each size class associated with the ODE.
hold <- abundance_at_size %>% mutate(., length_round = round(length, digits = 2)) %>% group_by(length_round) %>% summarise(max = max(N))
abundance_at_size <- data.frame("length"=hold$length_round, "N"=hold$max)

# Match each abundance_at_size$length from ODE to the Blue_sizes$mean_size.  Note that because we did not look at ages on a finer scale than seq(1:max_age, by = 1) for the VB based sizes-at-age, it is hard to match these exactly with the ODE derived abundance_at_size$length's
N_mean_size <- subset(abundance_at_size, length %in% c(Blue_sizes$mean_size))
colnames(N_mean_size) <- c("mean_size", "N_mean_size")
Blue_sizes <- left_join(Blue_sizes, N_mean_size, by = "mean_size")
# Deal with the annoying NA's for Blue_sizes$N_mean_size that arise from no individuals in those size classes (likely because they exceed Linf)
Blue_sizes$N_mean_size <- ifelse(is.na(Blue_sizes$N_mean_size), 0, Blue_sizes$N_mean_size)

df_len <- nrow(Blue_sizes)
for(i in 1:df_len){
Blue_sizes$size_obs[i] <- list(rnorm(Blue_sizes$N_mean_size[i], mean = Blue_sizes$mean_size[i], sd = Blue_sizes$sd[i]))
}

Blue_sizeobs_F0 <- unlist(Blue_sizes$size_obs)

hist(Blue_sizeobs_F0, breaks = 100, xlab = "Size (cm)", main = "Blue F=0", ylim = c(0,25))

########### Higher values of F for Blue

#F=0.05--------------------------------------------------------
#Solve ODE up to Lfish:
#Lfish <- 25 # Blue Lfish
Lfish <- 21.02 #Blue rockfish Lfish
len_noF = round(10000 *(Lfish/Linf))
l=seq(1,Lfish,length=len_noF) #Lfish
parms=c(k=0.172, M=0.14, Linf=38.15);
outF0=ode(y,l,eqn,parms)
#outF0

#Then add the M+F into it after Lfish
##add M+ F to eqn
eqn1=function(l,y,parms) {
  N1=y[1]
k=parms[1]
M=parms[2]
Linf=parms[3]
Fi=parms[4]
 dN=numeric(1); # Vector to hold the derivatives
  dN[1]=-N1*(((M+Fi)-k)/(k*(Linf-l)))
  return(list(dN))
}

y=c(N1=20);
len_F <- round(10000 *((Linf-Lfish)/Linf))
l=seq(Lfish, round(parms[3]),length=len_F); #round(Linf-Lfish)
parms=c(k=0.172, M=0.14, Linf=38.15,Fi=0.05);
outF0.05=ode(y,l,eqn1,parms)
#outF0.05
out.s005=rbind(outF0,outF0.05)
#out.s005
#plot(out.s005[,2]~out.s005[,1],lwd=2,col=1,type='l',xlab="size class",ylab="N",ylim=c(20,25), main="Blue rockfish")

colnames(out.s005) <- c("length", "N")
abundance_at_size <- data.frame(out.s005)

#Convert ages to mean size-at-age based on von bertalanffy eq.
#t0 = -0.5, so x0 = 3.6375
Linf <- parms[3]
k <- parms[1]

size_at_age <- vector()
for (j in 0:max_age){
size_at_age[j+1] <- Linf*(1-exp(-k*(j-t0)))
}

mean_size <- vector()
lower_size <- vector()
upper_size <- vector()
sd <- vector()
#
for (s in 1:length(size_at_age)){
mean_size[s] <- size_at_age[s]
sd[s] <- (mean_size[s])*.1
lower_size[s] <- mean_size[s]-(3*sd[s])
upper_size[s] <- mean_size[s]+(3*sd[s])
upper_size[s] <- ifelse(upper_size[s] > Linf, 34, upper_size[s])
}

Blue_sizes <- data.frame(mean_size = size_at_age, sd = sd, lower_size = lower_size, upper_size = upper_size)
#
Blue_sizes <- Blue_sizes %>% mutate_all(.,funs(round(.,digits = 2)))
  
hold <- abundance_at_size %>% mutate(., length_round = round(length, digits = 2)) %>% group_by(length_round) %>% summarise(max = max(N))
abundance_at_size <- data.frame("length"=hold$length_round, "N"=hold$max)

N_mean_size <- subset(abundance_at_size, length %in% c(Blue_sizes$mean_size))
colnames(N_mean_size) <- c("mean_size", "N_mean_size")
Blue_sizes <- left_join(Blue_sizes, N_mean_size, by = "mean_size")
# Deal with the annoying NA's for Blue_sizes$N_mean_size that arise from no individuals in those size classes (likely because they exceed Linf)
Blue_sizes$N_mean_size <- ifelse(is.na(Blue_sizes$N_mean_size), 0, Blue_sizes$N_mean_size)

df_len <- nrow(Blue_sizes)
for(i in 1:df_len){
Blue_sizes$size_obs[i] <- list(rnorm(Blue_sizes$N_mean_size[i], mean = Blue_sizes$mean_size[i], sd = Blue_sizes$sd[i]))
}

Blue_sizeobs_F0.05 <- unlist(Blue_sizes$size_obs)

hist(Blue_sizeobs_F0.05, breaks = 100, xlab = "Size (cm)", main = "Blue F=0.05", ylim =c(0,15))

#F=0.1 ----------------------------------------------------------
y=c(N1=20);
len_F <- round(10000 *((Linf-Lfish)/Linf))
l=seq(Lfish, round(parms[3]),length=len_F); #round(Linf-Lfish)
parms=c(k=0.172, M=0.14, Linf=38.15, Fi=0.1);
outF0.1=ode(y,l,eqn1,parms)
#outF0.1
out.s01=rbind(outF0,outF0.1)
#out.s01
plot(out.s01[,2]~out.s01[,1],lwd=2,col=1,type='l',xlab="size class",ylab="N",ylim=c(20,25), main="Blue rockfish")

colnames(out.s01) <- c("length", "N")
abundance_at_size <- data.frame(out.s01)

#Convert ages to mean size-at-age based on von bertalanffy eq.
#t0 = -0.5, so x0 = 3.6375
Linf <- parms[3]
k <- parms[1]

size_at_age <- vector()
for (j in 0:max_age){
size_at_age[j+1] <- Linf*(1-exp(-k*(j-t0)))
}

mean_size <- vector()
lower_size <- vector()
upper_size <- vector()
sd <- vector()

for (s in 1:length(size_at_age)){
  mean_size[s] <- size_at_age[s]
  sd[s] <- (mean_size[s])*.1
  lower_size[s] <- mean_size[s]-(3*sd[s])
  upper_size[s] <- mean_size[s]+(3*sd[s])
  #upper_size[s] <- ifelse(upper_size[s] > Linf, 34, upper_size[s])
}
  
Blue_sizes <- data.frame(mean_size = size_at_age, sd = sd, lower_size = lower_size, upper_size = upper_size)

Blue_sizes <- Blue_sizes %>% mutate_all(.,funs(round(.,digits = 2)))

# round abundance_at_size for mean, lower_size, upper_size
#abundance_at_size <- data.frame(out)

size_at_age <- round(size_at_age, digits = 2)
lower_size <- round(lower_size, digits = 2)
upper_size <-round(upper_size, digits = 2)
  
hold <- abundance_at_size %>% mutate(., length_round = round(length, digits = 2)) %>% group_by(length_round) %>% summarise(max = max(N))
abundance_at_size <- data.frame("length"=hold$length_round, "N"=hold$max)

N_mean_size <- subset(abundance_at_size, length %in% c(Blue_sizes$mean_size))
colnames(N_mean_size) <- c("mean_size", "N_mean_size")
Blue_sizes <- left_join(Blue_sizes, N_mean_size, by = "mean_size")
# Deal with the annoying NA's for Blue_sizes$N_mean_size that arise from no individuals in those size classes (likely because they exceed Linf)
Blue_sizes$N_mean_size <- ifelse(is.na(Blue_sizes$N_mean_size), 0, Blue_sizes$N_mean_size)
df_len <- nrow(Blue_sizes)
for(i in 1:df_len){
Blue_sizes$size_obs[i] <- list(rnorm(Blue_sizes$N_mean_size[i], mean = Blue_sizes$mean_size[i], sd = Blue_sizes$sd[i]))
}

Blue_sizeobs_F0.1 <- unlist(Blue_sizes$size_obs)

hist(Blue_sizeobs_F0.1, breaks = 100, xlab = "Size (cm)", main = "Blue F=0.1", ylim =c(0,15))

#F=0.2--------------------------------------------------------------------------------------
y=c(N1=20);
len_F <- round(10000 *((Linf-Lfish)/Linf))
l=seq(Lfish, round(parms[3]),length=len_F); #round(Linf-Lfish)
parms=c(k=0.172, M=0.14, Linf=38.15, Fi=0.2);
outF0.2=ode(y,l,eqn1,parms)
#outF0.2
out.s02=rbind(outF0,outF0.2)
#out.s02
# plot(out.s[,2]~out.s[,1],lwd=2,col=1,type='l',xlab="size class",ylab="N",ylim=c(20,25), main="Blue rockfish")

colnames(out.s02) <- c("length", "N")
abundance_at_size <- data.frame(out.s02)

#Convert ages to mean size-at-age based on von bertalanffy eq.
#t0 = -0.5, so x0 = 3.6375
Linf <- parms[3]
k <- parms[1]

size_at_age <- vector()
for (j in 0:max_age){
size_at_age[j+1] <- Linf*(1-exp(-k*(j-t0)))
}

mean_size <- vector()
lower_size <- vector()
upper_size <- vector()
sd <- vector()

for (s in 1:length(size_at_age)){
  mean_size[s] <- size_at_age[s]
  sd[s] <- (mean_size[s])*.1
  lower_size[s] <- mean_size[s]-(3*sd[s])
  upper_size[s] <- mean_size[s]+(3*sd[s])
  #upper_size[s] <- ifelse(upper_size[s] > Linf, 34, upper_size[s])
}
  
Blue_sizes <- data.frame(mean_size = size_at_age, sd = sd, lower_size = lower_size, upper_size = upper_size)

Blue_sizes <- Blue_sizes %>% mutate_all(.,funs(round(.,digits = 2)))

# round abundance_at_size for mean, lower_size, upper_size
#abundance_at_size <- data.frame(out)

size_at_age <- round(size_at_age, digits = 2)
lower_size <- round(lower_size, digits = 2)
upper_size <-round(upper_size, digits = 2)
  
hold <- abundance_at_size %>% mutate(., length_round = round(length, digits = 2)) %>% group_by(length_round) %>% summarise(max = max(N))
abundance_at_size <- data.frame("length"=hold$length_round, "N"=hold$max)

N_mean_size <- subset(abundance_at_size, length %in% c(Blue_sizes$mean_size))
colnames(N_mean_size) <- c("mean_size", "N_mean_size")
Blue_sizes <- left_join(Blue_sizes, N_mean_size, by = "mean_size")
# Deal with the annoying NA's for Blue_sizes$N_mean_size that arise from no individuals in those size classes (likely because they exceed Linf)
Blue_sizes$N_mean_size <- ifelse(is.na(Blue_sizes$N_mean_size), 0, Blue_sizes$N_mean_size)

df_len <- nrow(Blue_sizes)
for(i in 1:df_len){
Blue_sizes$size_obs[i] <- list(rnorm(Blue_sizes$N_mean_size[i], mean = Blue_sizes$mean_size[i], sd = Blue_sizes$sd[i]))
}

Blue_sizeobs_F0.2 <- unlist(Blue_sizes$size_obs)

hist(Blue_sizeobs_F0.2, breaks = 100, xlab = "Size (cm)", main = "Blue F=0.2", ylim =c(0,15))

#F=0.25----------------------------------------------------------------------
y=c(N1=20);
len_F <- round(10000 *((Linf-Lfish)/Linf))
l=seq(Lfish, round(parms[3]),length=len_F); #round(Linf-Lfish)
parms=c(k=0.172, M=0.14, Linf=38.15, Fi=0.25);
outF0.25=ode(y,l,eqn1,parms)
#outF0.25
out.s025=rbind(outF0,outF0.25)
#out.s025
# plot(out.s[,2]~out.s[,1],lwd=2,col=1,type='l',xlab="size class",ylab="N",ylim=c(20,25), main="Blue rockfish")

colnames(out.s025) <- c("length", "N")
abundance_at_size <- data.frame(out.s025)

#Convert ages to mean size-at-age based on von bertalanffy eq.
#t0 = -0.5, so x0 = 3.6375
Linf <- parms[3]
k <- parms[1]

size_at_age <- vector()
for (j in 0:max_age){
size_at_age[j+1] <- Linf*(1-exp(-k*(j-t0)))
}

mean_size <- vector()
lower_size <- vector()
upper_size <- vector()
sd <- vector()

for (s in 1:length(size_at_age)){
  mean_size[s] <- size_at_age[s]
  sd[s] <- (mean_size[s])*.1
  lower_size[s] <- mean_size[s]-(3*sd[s])
  upper_size[s] <- mean_size[s]+(3*sd[s])
  #upper_size[s] <- ifelse(upper_size[s] > Linf, 34, upper_size[s])
}
  
Blue_sizes <- data.frame(mean_size = size_at_age, sd = sd, lower_size = lower_size, upper_size = upper_size)

Blue_sizes <- Blue_sizes %>% mutate_all(.,funs(round(.,digits = 2)))

# round abundance_at_size for mean, lower_size, upper_size
#abundance_at_size <- data.frame(out)

size_at_age <- round(size_at_age, digits = 2)
lower_size <- round(lower_size, digits = 2)
upper_size <-round(upper_size, digits = 2)
  
hold <- abundance_at_size %>% mutate(., length_round = round(length, digits = 2)) %>% group_by(length_round) %>% summarise(max = max(N))
abundance_at_size <- data.frame("length"=hold$length_round, "N"=hold$max)

N_mean_size <- subset(abundance_at_size, length %in% c(Blue_sizes$mean_size))
colnames(N_mean_size) <- c("mean_size", "N_mean_size")
Blue_sizes <- left_join(Blue_sizes, N_mean_size, by = "mean_size")

# Deal with the annoying NA's for Blue_sizes$N_mean_size that arise from no individuals in those size classes (likely because they exceed Linf)
Blue_sizes$N_mean_size <- ifelse(is.na(Blue_sizes$N_mean_size), 0, Blue_sizes$N_mean_size)

df_len <- nrow(Blue_sizes)
for(i in 1:df_len){
Blue_sizes$size_obs[i] <- list(rnorm(Blue_sizes$N_mean_size[i], mean = Blue_sizes$mean_size[i], sd = Blue_sizes$sd[i]))
}

Blue_sizeobs_F0.25 <- unlist(Blue_sizes$size_obs)

hist(Blue_sizeobs_F0.25, breaks = 100, xlab = "Size (cm)", main = "Blue F=0.25", ylim =c(0,15))

#Call all histograms of F values for a species at once
#lines(density(FBsizeF0.0), col="black", lwd=2) #for continuous line over histogram
#par(mfrow = c(3,2))

#To match with PISCO and mockdata, change to bin_w = 3
bins_F0 <- seq(0,(ceiling(max(Blue_sizeobs_F0)/bin_w))*bin_w, by = bin_w)
bins_F0.05 <- seq(0,(ceiling(max(Blue_sizeobs_F0.05)/bin_w))*bin_w, by = bin_w)
bins_F0.1 <- seq(0,(ceiling(max(Blue_sizeobs_F0.1)/bin_w))*bin_w, by = bin_w)
bins_F0.2 <- seq(0,(ceiling(max(Blue_sizeobs_F0.2)/bin_w))*bin_w, by = bin_w)
bins_F0.25 <- seq(0,(ceiling(max(Blue_sizeobs_F0.25)/bin_w))*bin_w, by = bin_w)

hist(Blue_sizeobs_F0, breaks = c(bins_F0), xlab = "Size (cm)", main = "Blue F=0", ylim =c(0,120)) #Also will need to add to hist function:  prob = TRUE
hist(Blue_sizeobs_F0.05, breaks = c(bins_F0.05), xlab = "Size (cm)", main = "Blue F=0.05", ylim =c(0,120), add = TRUE, col = 2)
hist(Blue_sizeobs_F0.1, breaks = c(bins_F0.1), xlab = "Size (cm)", main = "Blue F=0.1", ylim =c(0,120), add = TRUE, col = "blue")
hist(Blue_sizeobs_F0.2, breaks = c(bins_F0.2), xlab = "Size (cm)", main = "Blue F=0.2", ylim =c(0,120), add = TRUE, col = "green")
hist(Blue_sizeobs_F0.25, breaks = c(bins_F0.25), xlab = "Size (cm)", main = "Blue F=0.25", ylim =c(0,120), add = TRUE, col = "black")
legend("topleft", c("F=0", "F=0.05", "F=0.1","F=0.2","F=0.25"), fill=c("white", "red", "blue", "green", "black"), bty = "n")


