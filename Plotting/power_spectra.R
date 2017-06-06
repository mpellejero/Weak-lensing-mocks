library(ggplot2)

dir <- "../Output_Rsp/"
dir <- "../Output_3deg_Rsp/"
dir <- "../Output_3deg_1Rvir/"
dir <- "../Output_3deg_Rsp-Alias/"
#dir <- "../Output_3deg_1Rvir-Alias/"
#dir <- "../Output_3deg_1Rvir_large/"
#dir <- "../Output_1Rvir/"
#dir <- "../Output_Rand_1Rvir/"

dir <- "../Output_particles_6arcs/"
dir <- "../Output_particles_6arcs_alias/"
dir <- "../Output_particles_random/"
dir <- "../Output_particles_maps/"

factor = 3*pi/2
#factor = 1.75
dir2 <- "../Output_3deg_1Rvir-Alias/"
dir2 <- "../Output_3deg_1Rvir/"
dir2 <- "../Output_particles_6arcs_alias/"

L_limit = 1.0e5

plt <- ggplot() +
  scale_x_log10(limit=c(100,L_limit)) + scale_y_log10() + ggtitle(paste(dir,factor)) +
  xlab("l") + ylab( expression( l^2~P(l) ) ) +
  theme(axis.title.x=element_text(face="italic"))

names <- c("testmap0_","testmap1_","testmap2_","testmap3_",
           "testmap4_","testmap5_","testmap6_","testmap7_"
           ,"testmap8_","testmap9_")

ps_data3 <- read.csv(paste0(dir2,"testmap0_2.297000PS.csv"))
ps_data3 <- subset(ps_data3, PS > 0 & l > 100)
ps_data3$llPS <- ps_data3$PS*ps_data3$l**2*factor# /pi
plt <- plt + geom_line(data=ps_data3,aes(x=l,y=llPS,color=dir2))

ps_ave1 <- read.csv(paste0(dir,names[1],"2.297000PS.csv") )

#ps_ave1 <- ps_data3
#ps_ave1$PS <- 0*ps_ave1$PS
for(n in names[2:length(names)]){
  print(n)
  ps_data3 <- read.csv(paste0(dir,n,"2.297000PS.csv") )
  ps_ave1$PS <- ps_ave1$PS + ps_data3$PS   
}
ps_ave1$PS <- ps_ave1$PS/length(names)
ps_ave1$llPS <- ps_ave1$PS*ps_ave1$l**2*factor

plt <- plt + geom_line(data=subset(ps_ave1, PS > 0 & l > 100)
                       ,aes(x=l,y=llPS,color='halos'))

#ps_ave1$PS = ps_ave1$PS - ps_ave1$PS[nrow(ps_ave1)]
#ps_ave1$llPS <- ps_ave1$PS*ps_ave1$l**2*factor
#plt <- plt + geom_line(data=subset(ps_ave1,PS != 0),aes(x=l,y=llPS,color='shot subtracted'))

ps_ave <- ps_ave1
ps_ave$PS <- 0*ps_ave$PS
for(n in names){
  print(n)
  ps_data3 <- read.csv(paste0(dir,n,"1.075000PS.csv") )
  ps_ave$PS <- ps_ave$PS + ps_data3$PS   
}
ps_ave$PS <- ps_ave$PS/length(names)
ps_ave$llPS <- ps_ave$PS*ps_ave$l**2*factor

plt <- plt + geom_line(data=subset(ps_ave, PS > 0 & l > 100)
                       ,aes(x=l,y=llPS,color='halos'))

#ps_ave$PS = ps_ave$PS - ps_ave$PS[nrow(ps_ave)]
#ps_ave$llPS <- ps_ave$PS*ps_ave$l**2*factor
#plt <- plt + geom_line(data=subset(ps_ave,PS != 0),aes(x=l,y=llPS,color='shot subtracted'))

ps_ave$PS <- 0*ps_ave$PS
for(n in names){
  print(n)
  ps_data3 <- read.csv(paste0(dir,n,"0.489200PS.csv") )
  ps_ave$PS <- ps_ave$PS + ps_data3$PS   
}
ps_ave$PS <- ps_ave$PS/length(names)
ps_ave$llPS <- ps_ave$PS*ps_ave$l**2*factor

plt <- plt + geom_line(data=subset(ps_ave,PS > 0 & l > 100),aes(x=l,y=llPS,color='halos')) +
  geom_point(data=subset(ps_ave,PS != 0),aes(x=l,y=llPS))


#ps_ave$PS = ps_ave$PS - ps_ave$PS[nrow(ps_ave)]
#ps_ave$llPS <- ps_ave$PS*ps_ave$l**2*factor
#plt <- plt + geom_line(data=subset(ps_ave,PS != 0),aes(x=l,y=llPS,color='shot subtracted'))
  

##########  Read in CAMB power spectra

MDps <- read.table("../MDpowerSpectra/BPl_nonlinear_CAMB_zs0.5_Omega_M0.31_Omega_DE0.69_scalarAmp2.157e-09.dat"
                   , sep="",header = FALSE)
colnames(MDps) <- c("l","PS")

MDps$llPS <- MDps$PS*MDps$l**2/4/pi/pi
MDps <- subset(MDps,l>100)

#plt <- plt + geom_line(data=MDps
#                       ,aes(x=l,y=llPS,color='CAMB'))

MDps <- read.table("../MDpowerSpectra/BPl_nonlinear_CAMB_zs1.0_Omega_M0.31_Omega_DE0.69_scalarAmp2.157e-09.dat",
                   header=FALSE, sep="")
colnames(MDps) <- c("l","PS")

MDps$llPS <- MDps$PS*MDps$l**2/4/pi/pi
MDps <- subset(MDps,l>100)

#plt <- plt + geom_line(data=MDps
#                       ,aes(x=l,y=llPS,color='CAMB'))

MDpower <- read.table("../w4vipers/1/lensing2.0.10/mapPowerSpectrum.dat",
                header=FALSE,sep="")
colnames(MDpower) <- c("l","PS")

for (i in 2:99) {
  p <- read.table(paste0("../w4vipers/",i,"/lensing2.0.10/mapPowerSpectrum.dat"),
                  header = FALSE,sep="")
  
  MDpower$PS <- MDpower$PS + p$V2
}
MDpower$PS <- MDpower$PS/99.
MDpower$llPS <- MDpower$PS*MDpower$l**2
MDpower <- subset(MDpower,l>100)

plt <- plt + geom_line(data=MDpower
                  ,aes(x=l,y=llPS,color='particles'))

MDpower1 <- MDpower
MDpower1$PS <- 0*MDpower$PS
for (i in 1:99) {
  p <- read.table(paste0("../w4vipers/",i,"/lensing2.0.1/mapPowerSpectrum.dat"),
                  header = FALSE,sep="")
  p <- subset(p,V1>100)
  
  MDpower1$PS <- MDpower1$PS + p$V2
}
MDpower1$PS <- MDpower1$PS/99.
MDpower1$llPS <- MDpower1$PS*MDpower1$l**2
MDpower1 <- subset(MDpower1,l>100)

plt <- plt + geom_line(data=MDpower1
                       ,aes(x=l,y=llPS,color='particles'))

MDpower$PS <- 0*MDpower$PS
for (i in 1:99) {
  p <- read.table(paste0("../w4vipers/",i,"/lensing2.0.17/mapPowerSpectrum.dat"),
                  header = FALSE,sep="")
  p <- subset(p,V1>100)
  
  MDpower$PS <- MDpower$PS + p$V2
}
MDpower$PS <- MDpower$PS/99.
MDpower$llPS <- MDpower$PS*MDpower$l**2
MDpower <- subset(MDpower,l>100)

plt <- plt + geom_line(data=MDpower
                       ,aes(x=l,y=llPS,color='particles'))

plt

#y <- approx(MDpower1$l,MDpower1$llPS,ps_ave1$l)
#plot(log10(ps_ave1$l),ps_ave1$llPS/y$y)
