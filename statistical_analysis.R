#Statistical analysis of dip directions

library("circular")
#Calculation of azimuth
calculate_azimuth<- function(x)
{
  
    angle <- atan2(x[2], x[1])
    angle <- angle * (180/pi)
    if (angle < 0)
    {
      rezultat <- (angle + 360)
    }
    else
    {
      rezultat <- angle
    }
    return (rezultat)
    
}

#Vector normalization
norm_vec <- function(x)
{
  len <- sqrt(sum(x^2))
  y <- x/len
  return (y)
}

#Combinatorial results
surface <- read.table(".txt", header=TRUE, sep = ";", dec=".")
nrow(surface)

surface<-dplyr::filter(surface, DOC<1)#deleting collinear configurations
surface<-dplyr::filter(surface, Z_N<1)#deleting horizontal samples
surface<-dplyr::filter(surface, Z_N>0)#deleting vertical samples
N=nrow(surface)    #number of valid samples
N

###3D analysis
mean_dir3D=atan2(mean(surface$Y_N),mean(surface$X_N))*180/pi+360
mean_dir3D



###2D analysis
dip_direction=c(rep(0,N))
sumx=0
sumy=0


for(i in 1:nrow(surface)){
  my_vector=c(surface$X_N[i],surface$Y_N[i])
  new_vector=norm_vec(my_vector)
  sumx=sumx+new_vector[1]
  sumy=sumy+new_vector[2]
  dip_direction[i]=calculate_azimuth(new_vector)
}
sumx
sumy
mean_dir_atan2=atan2(sumy,sumx)*180/pi+360
mean_dir_atan2

C=sum(cos(dip_direction*pi/180)) #dip directions in angles are converted into radians
S=sum(sin(dip_direction*pi/180))
C
S

#Mean direction
mean_dir=atan(S/C)*180/pi+2*180  #C>0 S<0 (checking Fisher's conditions in Eq. 2.9)
mean_dir
mean_dir_v2=atan2(S,C)*180/pi
mean_dir_v2+360

#Median dip direction
x=circular::circular(dip_direction, units="degrees")
x[87:121] #example data
median_circular<-circular::median.circular(x)
median_circular

#Resultant length
resultant_length_squared=C^2+S^2
resultant_length=sqrt(resultant_length_squared)
resultant_length

#Mean resultant length
mean_resultant_length=sqrt((C/N)^2+(S/N)^2)
mean_resultant_length

#Verification of calculation of resultant length
mean_resultant_length*N
resultant_length

#Sample circular variance (V=1-mean_resultant_length)
circular_variance=1-mean_resultant_length
circular_variance

#Sample circular standard deviation
circular_stdev=sqrt(-2*log(mean_resultant_length))
circular_stdev
stdev<- circular::sd.circular(x = dip_direction*pi/180 )
stdev


#Sample circular dispersion
m2=(1/N)*sum(cos(2*(dip_direction*pi/180-mean_dir*pi/180))) #p-th (p=2) trigonometric moment about the mean direction
m2

circ_disp=(1-m2^2)/(2*mean_resultant_length^2)
circ_disp

#Circular standard error
circ_ster=sqrt(circ_disp/N)
circ_ster

#Confidence intervals using non-parametric methods (Fisher (1993) "Statistical analysis of circular data", p. 76)
left=mean_dir-asin(1.9604*circ_ster)*180/pi #left bound of the confidence interval for the mean direction
left
right=mean_dir+asin(1.9604*circ_ster)*180/pi #right bound of the confidence interval for the mean direction
right

