library(rethinking)

##As per https://www.science.org/doi/10.1126/science.adg8116?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed


d_all <- read.csv("hMPXV1_root_to_tip.data.csv",header=T)


year_lower <- 2012
y_upper <- 80
interval_prob <- 0.90
background_colour <- alpha("#AFDBF5", 0.1)
foreground_colour <- "#3B3B53"
lineageA_colour <- "#FFA000"
lineageB_colour <- "#4997D0"
density_colour <- "#002FA7" 



show_simulated <- FALSE
show_density <- TRUE

#most recent
year_max <- max(d_all$year)

#most recent
year_min <- min(d_all$year)

d_all$age = year_max - d_all$year
non_apobec_min <- min(d_all$non_apobec_snps)
d_all$all_mut = d_all$all_snps - non_apobec_min

d_all$mut = d_all$apobec_snps

mut_min <- min(d_all$mut)
mut_max <- max(d_all$mut)

#mean year of sampling
xbar <- mean(d_all$year)

#linear model of root to tip vs time
rtt <- quap(
  alist(
    mut ~ dnorm(mu, sigma),
    mu <- a + b * (year - xbar),
    a ~ dnorm(mut_min, 100),
    b ~ dlnorm(0, 1),
    sigma ~ dnorm(0, 20)
  ),
  data = d_all)

#print parameters
print(precis(rtt, prob = interval_prob))

#make a grid of times
years.seq <- seq(from=year_lower, to=year_max, by=1/52)
mu <- link(rtt, data=data.frame(year=years.seq))

# get means and intervals for the grid
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI, prob = 0.9)

#plot the data
par(fin=c(7, 7), col=foreground_colour, family = "Helvetica") # plot=c(0, 1, 0, 1), 

par(tck = -0.01)
#par(mgp = c(3, 0.5, 0))  

# Plot background
plot(mut ~ year, data=d_all, pch = 20, col='black', cex = 1, 
     #axes = FALSE, 
     xlab = "", ylab = "APOBEC3 mutations", las=1,
     xlim = c(year_lower, year_max), ylim = c(0, y_upper), 
     cex.lab = 1.5) 



rect(par("usr")[1], par("usr")[3],
     par("usr")[2], par("usr")[4],
     xlim = c(year_lower, year_max), ylim = c(0, y_upper),
     col = background_colour)

#points(mut ~ year, data=d, pch = 16, col='black', cex = 0.8)
points(mut ~ year, data=d_all, pch = 16, col='white', cex = 0.6)
points(mut ~ year, data=d_all, pch = 16, col=alpha(lineageA_colour, 0.5), cex = 0.6)
#text(d$year, d$mut-1, labels=d$name)

axis(side = 1, lwd = 1, las=1)
axis(side = 2, lwd = 1, las=1)

line_df<-as.data.frame(cbind(years.seq, mu.mean))

#plot the mean line
clip(x1=year_lower, x2=year_max, y1=0, y2=par("usr")[4])
lines(years.seq, mu.mean, col=lineageA_colour, lwd=2, xlab = "", ylab = "", 
      xlim = c(year_lower, year_max), ylim = c(0, y_upper))

#plot the intervals
shade(mu.PI, years.seq, col=alpha(lineageA_colour, 0.2), ylab = "", 
      xlim = c(year_lower, year_max), ylim = c(0, y_upper))


lower<-mu.PI[1,]
upper<-mu.PI[2,]


shade_df<-as.data.frame(cbind(years.seq, lower, upper))



if (show_simulated) {
  #simulate data on the grid times
  sim.muts <- sim(rtt, data=list(year=years.seq), N=1E10)
  muts.PI <- apply(sim.muts, 2, PI, prob = interval_prob)
  
  #plot the simulation intervals
  shade(muts.PI, years.seq, col=alpha(lineageA_colour, 0.2))
}


#extract a whole lot of sample lines
post <- extract.samples(rtt, N=1E200)


x_df<-as.data.frame(x_at_0)



#get the distribution of intersects with x-axis
x_at_0 <- ((0 - post$a) / post$b) + xbar
x_at_0_density <- density(x_at_0, adjust = 0.5)
mean_x_at_0 <- mean(x_at_0)

#try different numbers of pre-APOBEC3 era mutations
x_at_1 <- ((1 - post$a) / post$b) + xbar
x_at_1_density <- density(x_at_1, adjust = 0.5)
mean_x_at_1 <- mean(x_at_1)
x_at_2 <- ((2 - post$a) / post$b) + xbar
x_at_2_density <- density(x_at_2, adjust = 0.5)
mean_x_at_2 <- mean(x_at_2)
x_at_3 <- ((3 - post$a) / post$b) + xbar
x_at_3_density <- density(x_at_3, adjust = 0.5)
mean_x_at_3 <- mean(x_at_3)

#lines(x_at_0_density$x, x_at_0_density$y * 10)
if (show_density) {
  par(fg = alpha("#4997D0", 0.1))
  polygon(x_at_0_density$x, x_at_0_density$y * 10, col = alpha("#4997D0", 0.3), border = NULL)
  par(fg=foreground_colour)
}

par(mar=c(1,1,1,1))


  
#print the intervals
origin <- HPDI(x_at_0, prob = interval_prob)
cat("\n\nmean_x at y = 0: ", mean(x_at_0), " [", origin[1], ", ", origin[2], "]\n")

origin1 <- HPDI(x_at_1, prob = interval_prob)
cat("mean_x at y = 1: ", mean(x_at_1), " [", origin1[1], ", ", origin1[2], "]\n")

origin2 <- HPDI(x_at_2, prob = interval_prob)
cat("mean_x at y = 2: ", mean(x_at_2), " [", origin2[1], ", ", origin2[2], "]\n")
origin3 <- HPDI(x_at_3, prob = interval_prob)
cat("mean_x at y = 3: ", mean(x_at_3), " [", origin3[1], ", ", origin3[2], "]\n\n")

  
dev.off()
  

  
  