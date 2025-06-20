#
# Tandem of C M/M/1 queues


rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(simmer)

Empirical_AAoI = function(td) {
  
  # Input
  #   > td : simmer object
  # 
  # Output
  #   > Empirical Average Age of Information
  
  td.arr = get_mon_arrivals(td) %>% arrange(start_time)
  
  Yn = td.arr %>% 
    filter(finished) %>% 
    mutate(Yn=c(end_time[1],diff(end_time))) %>% 
    pull(Yn) %>%
    lead(default=0)
  Tn = td.arr %>% 
    filter(finished) %>% 
    mutate(Tn=end_time-start_time) %>% 
    pull(Tn)
  Hn = Yn*(Tn + Yn/2)
  maxT = max(td.arr$end_time)
  
  return(sum(Hn)/maxT)
}

AAoI0 = function(la,mu) {
  
  # Input
  #   > la : arrival rate
  #   > mu : service rate
  # 
  # Output
  #   > Theoretical Average Age of Information of a M/M/1 queue with zero-aged updates
  
  t1 = 1/la + 1/mu
  t2= (la**2)/(mu**2)/(mu-la)
  return(t1+t2)
}

CorrTerm = function(la,ga,mu) {
  
  # Input
  #   > la : arrival rate
  #   > mu : service rate (queue of interest, second one in the tandem)
  #   > ga : service rates (first queue in the tandem)
  # 
  # Output
  #   > Exact Correction Term for a 2-queued tandem (equation in Lemma 2 multiplied by lambda)
  
  t1 = 1/ga
  t2 = (la**2)/(ga**2)/(ga-la)
  t3 = (la**2)/ga/mu/(ga+mu-la)
  return(t1+t2+t3)
}

Bounds_CT = function(la,gammas) {
  
  # Input
  #   > la     : arrival rate
  #   > gammas : vector of services rate
  # 
  # Output
  #   > Lower and upper bounds for the correction term in a tandem of M/M/1 queues, by Corollary 1.
  
  C = length(gammas) - 1
  t1 = round(sum(1/(gammas[1:C]-la)),3)
  t2 = round(sqrt(sum(1/(gammas[1:C]-la)**2)),3)
  return(c(t1-t2,t1+t2))
}

AAoI_Bounds_Tandem_MM1 = function(la=1,Ns=2,rhos=seq(.1,.9,length.out=Ns),seed=1234,T=10000,niter=100) {
  
  # Input
  #   > la   : arrival rate
  #   > Ns   : number of servers (if the load rates are not provided)
  #   > rhos : vector of load rates (if the number of servers is not provided)
  #   > seed : random seed
  #   > T    : time horizon in queueing simulation
  #   > iter : iterations per model
  #
  # Output
  # List with the following elements:
  #   > service_rates : vector of service rates (of length Ns)
  #   > arrival_rate  : arrival rate
  #   > load_rates    : vector of load rates (of length Ns)
  #   > bounds        : AAoI bounds of Equation (16)
  #   > AAoI_emp      : empirical AAoI (vector of niter elements)
  
  mu = la/rhos
  Ns = length(rhos)
  mm1 = function(., id) {
    seize(., paste0("mm1_q", id), 1) %>% 
      timeout(function() rexp(1,mu[id])) %>%
      release(paste0("mm1_q", id), 1)
  }
  
  Results = list()
  Results[['service_rates']] = mu
  Results[['arrival_rate']]  = la
  Results[['load_rates']]    = rhos
  Results[['bounds']]        = AAoI0(la,mu[Ns]) + Bounds_CT(la,mu)
  Results[['AAoI_emp']]      = rep(0,niter)
  
  # set up the point of attachment for the generator
  traj = trajectory()
  for (n in 1:Ns) traj = traj %>% mm1(n)
  
  set.seed(seed)
  
  for (i in 1:niter) {
    # prepare and run the simulation environment with the resources and generators required:
    tandem = simmer()
    for (n in 1:Ns) tandem %>% 
      add_resource(paste0("mm1_q", n))
    
    tandem %>%
      add_generator('update ', traj, function() {rexp(1,la)}, mon=2) %>%
      run(T)
    
    df.arr    = get_mon_arrivals(tandem) %>% arrange(start_time)
    df.arr.pr = get_mon_arrivals(tandem, per_resource=TRUE) %>% arrange(start_time)
    
    Results[['AAoI_emp']][i] = Empirical_AAoI(tandem)
  }
  return(Results)
}

la = 1
set.seed(2143)

# 3 servers (3!=6 different orders)
o1 = c(1,2,3)
o2 = c(1,3,2)
o3 = c(2,1,3)
o4 = c(2,3,1)
o5 = c(3,1,2)
o6 = c(3,2,1)

rhos3.1 = seq(.1,.9,length.out=3)
rhos3.2 = rhos3.1[o2]
rhos3.3 = rhos3.1[o3]
rhos3.4 = rhos3.1[o4]
rhos3.5 = rhos3.1[o5]
rhos3.6 = rhos3.1[o6]

AAoI.results.3servers.1 = AAoI_Bounds_Tandem_MM1(rhos=rhos3.1) 
AAoI.results.3servers.2 = AAoI_Bounds_Tandem_MM1(rhos=rhos3.2)
AAoI.results.3servers.3 = AAoI_Bounds_Tandem_MM1(rhos=rhos3.3)
AAoI.results.3servers.4 = AAoI_Bounds_Tandem_MM1(rhos=rhos3.4)
AAoI.results.3servers.5 = AAoI_Bounds_Tandem_MM1(rhos=rhos3.5)
AAoI.results.3servers.6 = AAoI_Bounds_Tandem_MM1(rhos=rhos3.6)

#
# Table II
#
Summary.df.3servers = tibble(s1=0, s2=0, s3=0, 
                             av.AAoI_emp=0, sd.AAoI_emp=0, AAoI_LB=0, AAoI_UB=0)
for (i in 1:6) {
  name = paste0('AAoI.results.3servers.',i)
  df = get(name)
  Summary.df.3servers[i,1:3] = t(df$load_rates)
  Summary.df.3servers[i,4]   = mean(df$AAoI_emp)
  Summary.df.3servers[i,5]   = sd(df$AAoI_emp)
  Summary.df.3servers[i,6]   = df$bounds[1]
  Summary.df.3servers[i,7]   = df$bounds[2]
}

  #      s1    s2    s3 av.AAoI_emp sd.AAoI_emp AAoI_LB AAoI_UB
  #   <dbl> <dbl> <dbl>       <dbl>       <dbl>   <dbl>   <dbl>
  # 1   0.1   0.5   0.9       10.1         1.53    9.29    11.3
  # 2   0.1   0.9   0.5        9.86        1.40    1.86    19.9
  # 3   0.5   0.1   0.9       10.1         1.87    9.29    11.3
  # 4   0.5   0.9   0.1       10.3         1.67    2.05    20.2
  # 5   0.9   0.1   0.5       10.1         1.73    1.86    19.9
  # 6   0.9   0.5   0.1       10.2         1.76    2.05    20.2


# 6 servers (6!=720 different orders)
rhos6.1 = seq(.1,.9,length.out=6)
rhos6.2 = rev(rhos6.1)
rhos6.3 = sample(rhos6.1,6,replace=F)
rhos6.4 = sample(rhos6.1,6,replace=F)
rhos6.5 = sample(rhos6.1,6,replace=F)
rhos6.6 = sample(rhos6.1,6,replace=F)

AAoI.results.6servers.1 = AAoI_Bounds_Tandem_MM1(rhos=rhos6.1) 
AAoI.results.6servers.2 = AAoI_Bounds_Tandem_MM1(rhos=rhos6.2)
AAoI.results.6servers.3 = AAoI_Bounds_Tandem_MM1(rhos=rhos6.3)
AAoI.results.6servers.4 = AAoI_Bounds_Tandem_MM1(rhos=rhos6.4)
AAoI.results.6servers.5 = AAoI_Bounds_Tandem_MM1(rhos=rhos6.5)
AAoI.results.6servers.6 = AAoI_Bounds_Tandem_MM1(rhos=rhos6.6)

#
# Table III
#
Summary.df.6servers = tibble(s1=0, s2=0, s3=0, s4=0, s5=0, s6=0, 
                             av.AAoI_emp=0, sd.AAoI_emp=0, AAoI_LB=0, AAoI_UB=0)
for (i in 1:6) {
  name = paste0('AAoI.results.6servers.',i)
  df = get(name)
  Summary.df.6servers[i,1:6] = t(df$load_rates)
  Summary.df.6servers[i,7]   = mean(df$AAoI_emp)
  Summary.df.6servers[i,8]   = sd(df$AAoI_emp)
  Summary.df.6servers[i,9]   = df$bounds[1]
  Summary.df.6servers[i,10]  = df$bounds[2]
}

  #     s1    s2    s3    s4    s5    s6 av.AAoI_emp sd.AAoI_emp AAoI_LB AAoI_UB
  #  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>       <dbl>       <dbl>   <dbl>   <dbl>
  # 1  0.1   0.26  0.42  0.58  0.74  0.9         14.4        1.79   11.3     17.9
  # 2  0.9   0.74  0.58  0.42  0.26  0.1         14.4        1.95    5.83    25.0
  # 3  0.26  0.1   0.9   0.58  0.74  0.42        14.3        2.01    5.69    24.8
  # 4  0.1   0.74  0.26  0.58  0.9   0.42        14.5        1.72    5.69    24.8
  # 5  0.74  0.58  0.9   0.1   0.42  0.26        14.5        1.50    5.78    24.9
  # 6  0.1   0.9   0.42  0.26  0.74  0.58        14.4        1.52    5.60    24.6
 
# 10 servers (10!=3628800 different orders)
rhos10.1 = seq(.1,.9,length.out=10)
rhos10.2 = rev(rhos10.1)
rhos10.3 = sample(rhos10.1,10,replace=F)
rhos10.4 = sample(rhos10.1,10,replace=F)
rhos10.5 = sample(rhos10.1,10,replace=F)
rhos10.6 = sample(rhos10.1,10,replace=F)

AAoI.results.10servers.1 = AAoI_Bounds_Tandem_MM1(rhos=rhos10.1) 
AAoI.results.10servers.2 = AAoI_Bounds_Tandem_MM1(rhos=rhos10.2)
AAoI.results.10servers.3 = AAoI_Bounds_Tandem_MM1(rhos=rhos10.3)
AAoI.results.10servers.4 = AAoI_Bounds_Tandem_MM1(rhos=rhos10.4)
AAoI.results.10servers.5 = AAoI_Bounds_Tandem_MM1(rhos=rhos10.5)
AAoI.results.10servers.6 = AAoI_Bounds_Tandem_MM1(rhos=rhos10.6)

#
# Table IV
#
Summary.df.10servers = tibble(s1=0, s2=0, s3=0, s4=0, s5=0, s6=0, s7=0, s8=0, s9=0, s10=0, 
                              av.AAoI_emp=0, sd.AAoI_emp=0, AAoI_LB=0, AAoI_UB=0)
for (i in 1:6) {
  name = paste0('AAoI.results.10servers.',i)
  df = get(name)
  Summary.df.10servers[i,1:10] = t(df$load_rates)
  Summary.df.10servers[i,11]   = mean(df$AAoI_emp)
  Summary.df.10servers[i,12]   = sd(df$AAoI_emp)
  Summary.df.10servers[i,13]   = df$bounds[1]
  Summary.df.10servers[i,14]   = df$bounds[2]
}

  #    s1    s2    s3    s4    s5    s6    s7    s8    s9   s10 av.AAoI_emp sd.AAoI_emp AAoI_LB AAoI_UB
  # <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>       <dbl>       <dbl>   <dbl>   <dbl>
  # 1 0.1   0.189 0.278 0.367 0.456 0.544 0.633 0.722 0.811 0.9          20.9        1.88    15.6    26.7
  # 2 0.9   0.811 0.722 0.633 0.544 0.456 0.367 0.278 0.189 0.1          21.0        2.06    11.4    32.5
  # 3 0.1   0.544 0.367 0.722 0.811 0.278 0.9   0.189 0.456 0.633        21.1        2.12    11.1    32.0
  # 4 0.367 0.456 0.1   0.9   0.633 0.278 0.189 0.722 0.544 0.811        20.5        1.84    11.6    31.0
  # 5 0.811 0.544 0.367 0.456 0.722 0.278 0.9   0.189 0.633 0.1          20.8        1.84    11.4    32.5
  # 6 0.9   0.278 0.367 0.1   0.544 0.189 0.633 0.722 0.456 0.811        20.8        1.91    11.6    31.0


#
# Assess the validity of the AAoI0 and CorrTerm functions in a tandem of two queues.
# (not in article)
# 

rhos2.1 = c(0.1,0.9)
rhos2.2 = c(0.9,0.1)
rhos2.3 = c(0.4,0.6)
rhos2.4 = c(0.6,0.4)
rhos2.5 = c(0.25,0.3)
rhos2.6 = c(0.3,0.25)

AAoI.results.2servers.1 = AAoI_Bounds_Tandem_MM1(rhos=rhos2.1) 
AAoI.results.2servers.2 = AAoI_Bounds_Tandem_MM1(rhos=rhos2.2)
AAoI.results.2servers.3 = AAoI_Bounds_Tandem_MM1(rhos=rhos2.3)
AAoI.results.2servers.4 = AAoI_Bounds_Tandem_MM1(rhos=rhos2.4)
AAoI.results.2servers.5 = AAoI_Bounds_Tandem_MM1(rhos=rhos2.5)
AAoI.results.2servers.6 = AAoI_Bounds_Tandem_MM1(rhos=rhos2.6)

Summary.df.2servers = tibble(s1=0, s2=0, th.AAoI=0,
                             av.AAoI_emp=0, sd.AAoI_emp=0, AAoI_LB=0, AAoI_UB=0)
for (i in 1:6) {
  name = paste0('AAoI.results.2servers.',i)
  df = get(name)
  mu = df$service_rates[2]
  ga = df$service_rates[1]
  la = df$arrival_rate
  Summary.df.2servers[i,1:2] = t(df$load_rates)
  Summary.df.2servers[i,3]   = AAoI0(la,mu) + CorrTerm(la,ga,mu)
  Summary.df.2servers[i,4]   = mean(df$AAoI_emp)
  Summary.df.2servers[i,5]   = sd(df$AAoI_emp)
  Summary.df.2servers[i,6]   = df$bounds[1]
  Summary.df.2servers[i,7]   = df$bounds[2]
}

  #     s1    s2 th.AAoI av.AAoI_emp sd.AAoI_emp AAoI_LB AAoI_UB
  #  <dbl> <dbl>   <dbl>       <dbl>       <dbl>   <dbl>   <dbl>
  # 1  0.1   0.9     9.30        9.33      1.80      9.19    9.41
  # 2  0.9   0.1     9.30        8.87      1.25      1.10   19.1 
  # 3  0.4   0.6     2.72        2.73      0.0476    2.14    3.47
  # 4  0.6   0.4     2.72        2.72      0.0515    1.51    4.51
  # 5  0.25  0.3     1.62        1.62      0.0163    1.34    2.00
  # 6  0.3   0.25    1.62        1.62      0.0158    1.27    2.13

  # Empirical AAoI = Theoretical AAoI

# 
# Homogeneous servers (not in the article)
# 

# 3 servers

Ns = 3

AAoI.results.3eq.1 = AAoI_Bounds_Tandem_MM1(rhos=rep(0.1,Ns))
AAoI.results.3eq.2 = AAoI_Bounds_Tandem_MM1(rhos=rep(0.2,Ns))
AAoI.results.3eq.3 = AAoI_Bounds_Tandem_MM1(rhos=rep(0.3,Ns))
AAoI.results.3eq.4 = AAoI_Bounds_Tandem_MM1(rhos=rep(0.4,Ns))
AAoI.results.3eq.5 = AAoI_Bounds_Tandem_MM1(rhos=rep(0.5,Ns))
AAoI.results.3eq.6 = AAoI_Bounds_Tandem_MM1(rhos=rep(0.6,Ns))
AAoI.results.3eq.7 = AAoI_Bounds_Tandem_MM1(rhos=rep(0.7,Ns))
AAoI.results.3eq.8 = AAoI_Bounds_Tandem_MM1(rhos=rep(0.8,Ns))
AAoI.results.3eq.9 = AAoI_Bounds_Tandem_MM1(rhos=rep(0.9,Ns))

Summary.df.3eq = tibble(mu=0, rho=0, av.AAoI_emp=0, sd.AAoI_emp=0, AAoI_LB=0, AAoI_UB=0)
for (i in 1:9) {
  name = paste0('AAoI.results.3eq.',i)
  df = get(name)
  Summary.df.3eq[i,1] = df$service_rates[1]
  Summary.df.3eq[i,2] = df$load_rates[1] 
  Summary.df.3eq[i,3] = mean(df$AAoI_emp)
  Summary.df.3eq[i,4] = sd(df$AAoI_emp)
  Summary.df.3eq[i,5] = df$bounds[1]
  Summary.df.3eq[i,6] = df$bounds[2]
}

  #      mu   rho av.AAoI_emp sd.AAoI_emp AAoI_LB AAoI_UB
  #   <dbl> <dbl>       <dbl>       <dbl>   <dbl>   <dbl>
  # 1 10      0.1        1.31      0.0149    1.17    1.48
  # 2  5      0.2        1.64      0.0161    1.36    2.06
  # 3  3.33   0.3        2.06      0.0178    1.59    2.80
  # 4  2.5    0.4        2.62      0.0247    1.90    3.78
  # 5  2      0.5        3.46      0.0418    2.34    5.16
  # 6  1.67   0.6        4.80      0.102     3.02    7.26
  # 7  1.43   0.7        7.13      0.250     4.21   10.8 
  # 8  1.25   0.8       12.0       0.671     6.70   18.0 
  # 9  1.11   0.9       26.5       3.20     14.5    39.9 

# 6 servers

Ns = 6

AAoI.results.6eq.1 = AAoI_Bounds_Tandem_MM1(rhos=rep(0.1,Ns))
AAoI.results.6eq.2 = AAoI_Bounds_Tandem_MM1(rhos=rep(0.2,Ns))
AAoI.results.6eq.3 = AAoI_Bounds_Tandem_MM1(rhos=rep(0.3,Ns))
AAoI.results.6eq.4 = AAoI_Bounds_Tandem_MM1(rhos=rep(0.4,Ns))
AAoI.results.6eq.5 = AAoI_Bounds_Tandem_MM1(rhos=rep(0.5,Ns))
AAoI.results.6eq.6 = AAoI_Bounds_Tandem_MM1(rhos=rep(0.6,Ns))
AAoI.results.6eq.7 = AAoI_Bounds_Tandem_MM1(rhos=rep(0.7,Ns))
AAoI.results.6eq.8 = AAoI_Bounds_Tandem_MM1(rhos=rep(0.8,Ns))
AAoI.results.6eq.9 = AAoI_Bounds_Tandem_MM1(rhos=rep(0.9,Ns))

Summary.df.6eq = tibble(mu=0, rho=0, av.AAoI_emp=0, sd.AAoI_emp=0, AAoI_LB=0, AAoI_UB=0)
for (i in 1:9) {
  name = paste0('AAoI.results.6eq.',i)
  df = get(name)
  Summary.df.6eq[i,1] = df$service_rates[1]
  Summary.df.6eq[i,2] = df$load_rates[1] 
  Summary.df.6eq[i,3] = mean(df$AAoI_emp)
  Summary.df.6eq[i,4] = sd(df$AAoI_emp)
  Summary.df.6eq[i,5] = df$bounds[1]
  Summary.df.6eq[i,6] = df$bounds[2]
}


# 10 servers

Ns = 10

AAoI.results.10eq.1 = AAoI_Bounds_Tandem_MM1(rhos=rep(0.1,Ns))
AAoI.results.10eq.2 = AAoI_Bounds_Tandem_MM1(rhos=rep(0.2,Ns))
AAoI.results.10eq.3 = AAoI_Bounds_Tandem_MM1(rhos=rep(0.3,Ns))
AAoI.results.10eq.4 = AAoI_Bounds_Tandem_MM1(rhos=rep(0.4,Ns))
AAoI.results.10eq.5 = AAoI_Bounds_Tandem_MM1(rhos=rep(0.5,Ns))
AAoI.results.10eq.6 = AAoI_Bounds_Tandem_MM1(rhos=rep(0.6,Ns))
AAoI.results.10eq.7 = AAoI_Bounds_Tandem_MM1(rhos=rep(0.7,Ns))
AAoI.results.10eq.8 = AAoI_Bounds_Tandem_MM1(rhos=rep(0.8,Ns))
AAoI.results.10eq.9 = AAoI_Bounds_Tandem_MM1(rhos=rep(0.9,Ns))

Summary.df.10eq = tibble(mu=0, rho=0, av.AAoI_emp=0, sd.AAoI_emp=0, AAoI_LB=0, AAoI_UB=0)
for (i in 1:9) {
  name = paste0('AAoI.results.10eq.',i)
  df = get(name)
  Summary.df.10eq[i,1] = df$service_rates[1]
  Summary.df.10eq[i,2] = df$load_rates[1] 
  Summary.df.10eq[i,3] = mean(df$AAoI_emp)
  Summary.df.10eq[i,4] = sd(df$AAoI_emp)
  Summary.df.10eq[i,5] = df$bounds[1]
  Summary.df.10eq[i,6] = df$bounds[2]
}

cbind(cbind(Summary.df.3eq,Summary.df.6eq[,3:6]),Summary.df.10eq[,3:6])

#          mu rho av.AAoI_emp sd.AAoI_emp   AAoI_LB   AAoI_UB av.AAoI_emp sd.AAoI_emp   AAoI_LB   AAoI_UB av.AAoI_emp sd.AAoI_emp   AAoI_LB    AAoI_UB
# 1 10.000000 0.1    1.307717  0.01493273  1.166111  1.480111    1.612039  0.01666842  1.409111  1.905111    2.027470  0.01577019  1.768111   2.434111
# 2  5.000000 0.2    1.642588  0.01611794  1.356000  2.064000    2.305686  0.01632451  1.901000  3.019000    3.212688  0.01862752  2.710000   4.210000
# 3  3.333333 0.3    2.059067  0.01781798  1.589571  2.801571    3.192422  0.02100933  2.523571  4.439571    4.764777  0.02803454  3.909571   6.481571
# 4  2.500000 0.4    2.624014  0.02465593  1.896667  3.782667    4.422197  0.03863015  3.348667  6.330667    6.931932  0.05999278  5.506667   9.506667
# 5  2.000000 0.5    3.463250  0.04183661  2.336000  5.164000    6.235496  0.09183593  4.514000  8.986000   10.122243  0.14128011  7.750000  13.750000
# 6  1.666667 0.6    4.799025  0.10171661  3.019000  7.261000    9.076379  0.21433249  6.286000 12.994000   15.027200  0.32028360 11.140000  20.140000
# 7  1.428571 0.7    7.129532  0.25043249  4.210333 10.810333   13.976980  0.41852256  9.293333 19.727333   23.315515  0.61124549 16.843333  30.843333
# 8  1.250000 0.8   11.954959  0.67115062  6.703000 18.017000   23.618754  1.13137303 15.416000 33.304000   39.917613  1.97327151 28.360000  52.360000
# 9  1.111111 0.9   26.524081  3.20165018 14.462000 39.918000   51.798659  5.16111793 34.065000 74.315000   87.140203  8.73295198 63.190000 117.190000