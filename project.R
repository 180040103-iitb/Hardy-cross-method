#Declaring Variables

P0       <- 75000        #pump delivering pressure
dia      <- 0.05         #diameter of the original pipe
rho      <- 997          #density of the water
gval     <- 9.81         #gravitational acc of earth
pi       <- 22/7         #constant
kin_vis  <- 0.801*10^-6  #kinematic viscosity of water
epi      <- 0.26*10^-3   #roughness of the pipe
k        <-20            #no of pipes
l1       <- 0.6          #length of section from pump to pipe 1
l0       <-0.2           #length of section between consecutive pipes

#assume friction factor(ss) of pipe from pump to 1st pipe
ss = 0.02

#By mood digram equation 
ss <- 1/((-2*log10((epi/(3.7*dia))+((2.51*kin_vis)/sqrt(((P0/(rho*gval))*2*gval*(dia)^3*ss)/(dia+ss*l1)))))^2)


#V1 is the velocity of water in original pipe from pump to 1st pipe
v1 <- sqrt(((P0/(rho*gval))*2*gval*dia)/(dia+ss*l1))

#Q0 is the discharge in section from pump to 1st pipe
Q0 = (v1*pi*dia^2)/4

#for loop vector
n <- 1:(k-1)

#Declaring vector of friction factor
ff <- numeric(k-1)

#Assuming the initial friction factor values for the section of pipes 
ff[1:(k-1)] <-0.02

#Evaluating friction factor in each section
for (i in n) {
  
  ff[i] = 1/((-2*log10((epi/(3.7*dia))+((1.97*dia*kin_vis)/(Q0*(1-i/k)*sqrt(ff[i])))))^2)
  
}

#Evaluating discharge in each section
Q <- numeric(k-1)
for (i in n) {
  Q[i]=Q0*(1-i/k)
}

#Evaluating head loss in each section
hf <- numeric(k-1)
 
for (i in n) {
  hf[i]=(ff[i]*8*l0*(Q[i])^2)/(gval*pi^2*dia^5)
}

#Evaluating head loss in section from pump to pipe 1
hf0=(8*ss*l1*Q0^2)/(gval*pi^2*dia^5)


hff=sum(hf[1:(k-1)])
#Evaluating sum of head loss whole pipe
final_hf= hff+hf0

#calculating respective velocites 
vel <- numeric(k-1)

vel[k-1]=sqrt((P0/(rho*gval)-final_hf)*2*gval)

diam <- numeric(k-1)

diam[k-1] = sqrt((4*Q[k-1])/(vel[k-1]*pi))


r<-(k-2):1
for ( i in r) {
  diam[i]=((l0/dia^5)*(ff[i]*(k-1-i)^2)+1/(diam[i+1])^4)^-0.25
}

diam0=((l1/dia^5)*(ss*(k-1)^2)+1/(diam[1])^4)^-0.25
print(diam0)
for (i in n) {
  print(diam[i])
}

plot(diam)





