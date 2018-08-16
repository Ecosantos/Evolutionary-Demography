#################################################################
####	   	SCRIPT TO SHOW GRAPHICALLY THE IMPORTANCE TO 	      
####		    CORRECT VARIANCE ACCORDING WITH MAXIMUM 	
####				          VARIANCE POSSIBLE				
####		     (only for vital rates range [0,1])		 	
#################################################################


# As discussed in my qualification, as better detailed by Morris and Doak 2004
# and followed by Morris et al. 2011, vital rates range from 0 to 1 need to 
# by corrected they respective maximum variance

### References
#Morris, W. F., and D. F. Doak. 2004. Buffering of life histories against environmental stochasticity: accounting for a spurious correlation between the variabilities of vital rates and their contributions to fitness. The American naturalist 163:579–590.
#Morris, W. F., J. Altmann, D. K. Brockman, M. Cords, L. M. Fedigan, A. E. Pusey, T. S. Stoinski, A. M. Bronikowski, S. C. Alberts, and K. B. Strier. 2011. Low Demographic Variability in Wild Primate Populations: Fitness Impacts of Variation, Covariation, and Serial Correlation in Vital Rates. The American Naturalist 177:E14–E28.

# Maximum variance correction was also discussed in this post: 
https://math.stackexchange.com/questions/769923/maximum-possible-variance

#In this presentation I compare difference between corrected and uncorrected variance 
# Also compare Deterministic elasticity x Stochastic elasticity
#	this last comparing are checkable by difference between "Y axis" 
#		while difference between variance is checkable by differences in "X axis"

# Code
rm(list=ls())

library(popbio)
data(monkeyflower)

data(monkeyflower)

mim<- subset(monkeyflower, species == "cardinalis" &
site == "Carlon" & year != "pooled", select = c(4:19))

## convert data frame to list of matrices using split
mim1<-split(mim, 2000:2002)
mim2<-lapply(mim1, matrix, nrow=4, byrow=TRUE)
vr1<- pfister.plot(mim2)
length(vr1[,1])	#one value was deleted (a24) because mean and variance = 0 

A1<-mim2[1]
A2<-mim2[2]
A3<-mim2[3]

elas<-
popbio::stoch.sens(list(
mean(c(eval(A1),eval(A1),eval(A1))),
mean(c(eval(A2),eval(A2),eval(A2))),
mean(c(eval(A3),eval(A3),eval(A3)))))$elasticities



splitA(var2(mim2))$T*(mean(mim2)*(1-mean(mim2)))	#variance corrected

#then, we return to algamatade form

varmxcorrected<-splitA(var2(mim2))$T*(mean(mim2)*(1-mean(mim2))) +
 splitA(var2(mim2))$F

CV<-sqrt(varmxcorrected)/mean(mim2)*100


plot(vr1$cv, vr1$elas, xlab="CV", ylab="Elasticity",xlim=c(0.01,130),log="y",,type="n")
y2<- expression(S[11],G[21],G[31],G[41],
F[12],S[22],G[32],G[42],
                 F[13],R[23],S[33],G[43],
                 F[14],R[34],S[44])

text(CV, elas, y2,col="blue")
text(vr1$cv, vr1$elas, y2,col="red")

a<- cor.test(vr1$cv, vr1$elas, method="spearman")
acorrected<- cor.test(CV, elas, method="spearman")
a
text(100, .0015, substitute(rho == x, list(x=round(a$estimate,2))), col="red")
text(100, .0010, substitute(rho("corrected") == x, list(x=round(acorrected$estimate,2))), col="blue")

arrows(vr1$cv,elasticity(mean(mim2)),CV,elas,length = .1,code=2)


