#Declaring Variables

P0       <- 750000        #pump delivering pressure
rho      <- 997          #density of the water
gval     <- 9.81         #gravitational acc of earth
pi       <- 22/7         #constant
kin_vis  <- 0.801*10^-6  #kinematic viscosity of water
epi      <- 0.26*10^-3   #roughness of the pipe
k        <- 25           #no of pipes


Q_pipe <-numeric(k)
vel_pipe <-numeric(k)
dia_pipe <-numeric(k)

Q_sec  <-numeric(k)
vel_sec <-numeric(k)
dia_sec <-numeric(k)
area_sec <-numeric(k)
hf <- numeric(k)

ff <-numeric(k)
l <-numeric(k)
len_pump <-numeric(k)


l[1:k]<-0.2
l[1]<-0.2

n<-(k-1):1
m<-2:k
o<-1:k

for (i in o) {
  len_pump[i] = sum(l[1:i])
}

dia_sec[1:k] <-0.05
for (i in o) {
  area_sec[i] <- (pi*dia_sec[i]^2)/4
}

ff[1:k]=0.02
ff[1]  =  1/((-2*log10((epi/(3.7*dia_sec[1]))+((2.51*kin_vis*area_sec[1])/(dia_sec[1]*sqrt(ff[1])*(sqrt((P0*2*dia_sec[1]*area_sec[1]^2)/(rho*(dia_sec[1]+ff[1]*l[1]))))))))^2)

Q_sec[1] = (sqrt((P0*2*dia_sec[1]*area_sec[1]^2)/(rho*(dia_sec[1]+ff[1]*l[1]))))

Q_pipe[1:k] <- Q_sec[1]/k
#Q_pipe[10] = 1.5*Q_sec[1]/k
#Q_pipe[7] = 0.5*Q_sec[1]/k
#Q_pipe[15] = 0
#Q_pipe[17] = 4*Q_sec[1]/k


for (i in m) {
  Q_sec[i] <- Q_sec[1]-sum(Q_pipe[1:i-1])
}

for (i in o) {
  vel_sec[i] = Q_sec[i]/area_sec[i]
}

for (i in m) {
  ff[i]  =  1/((-2*log10((epi/(3.7*dia_sec[i]))+((2.51*kin_vis*area_sec[i])/(dia_sec[i]*sqrt(ff[i])*Q_sec[i]))))^2)
}

for (i in o) {
  hf[i] = (ff[i]*l[i]*vel_sec[i]^2)/(2*gval*dia_sec[i])
}

total_hf = sum(hf[1:k])

#dia_pipe[k] = ((8*Q_pipe[k]^2*rho)/(pi^2*(P0-gval*rho*total_hf)))^0.25
#print(dia_pipe[k])

vel_pipe[k]= sqrt(2*P0/rho-2*gval*total_hf)

for (i in n) {
 vel_pipe[i] = sqrt(((ff[i]*l[i]*vel_sec[i]^2)/dia_sec[i])+vel_pipe[i+1]^2)
}

for (i in o) {
  dia_pipe[i] = sqrt((4*Q_pipe[i])/(vel_pipe[i]*pi))
  print(dia_pipe[i])
}

plot(len_pump,dia_pipe,type = "o",main = "Variation of pipe diameters with length for same flow",xlab = "Length from pump in meter",ylab = "Diameter of pipes in meter",col="red",cex.lab="1.2",cex.main="1.4",cex.axis="1.2",ps="2",pch="0", lwd = "2",xaxs="r",yaxs="r")

