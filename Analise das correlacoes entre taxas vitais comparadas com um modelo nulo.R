###############################################################################
###                 ANÁLISE DO NÍVEL DE CORRELAÇÃO ENTRE AS TAXAS VITAIS	# 
###				DOS DADOS UTILIZADOS POR PFISTER 1998  		      #
###############################################################################
###                                                                           #
###      As análises que seguem são a tentativa de avaliar o quanto os dados  #
###    obtidos por Catherine Pfister em seu artigo clássico publicado na PNAS #
###  "Patterns of variance in stage-structured populations: evolutionary      #
### predictions and ecological implications" são auto-correlacionados entre 	#
###              si, uma vez que isso é de grande importância para		#
###		     sabermos como analisar de forma adequada esses ados          #         
###                                                                           #
###                                                                           #
###                                                                           #
###############################################################################
###   		Para facilitar o processo, eu utilizei algumas 			#
###	funções previamente criadas criadas por mim, com apenas algumas 		#
###	modificações para adequação à proposta desse script, que é a avaliação  #
###	do grau de correlação entre as taxas vitais analisadas			#		                                 
###                                                                           #
###############################################################################
     
library(popbio)
library(reshape)
library(corrplot)
setwd(choose.dir())
dir()


filelist = list.files(,pattern = ".*.txt")
filelist 


mylist <- lapply(filelist, read.table) 
mylist

summary(mylist)
as.numeric(unlist(mylist))

a<-lapply(mylist, as.numeric)

num<-lapply(list(mylist[[1]],mylist[[2]]),as.matrix)
num
summary(num)

nummean <-mean(num)
nummean

#================================================================================
###############         ==== CRIANDO AS FUNÇÕES UTILZIADAS====         	#########
#================================================================================


PLOTDOC<- function (A) 
{
    n <- length(A)
    if (class(A) != "list") {
        stop("A LIST of annual matrices is required")
    }
    if (length(unique(lapply(A, dim))) > 1) {
        stop("Matrices have different dimensions")
    }
    if (n < 2) {
        stop("A list of TWO or more annual matrices is required input")
    }
    col <- names(A)
    x <- dim(A[[1]])[1]
    row <- paste("a", 1:x, rep(1:x, each = x), sep = "")
    vr <- data.frame(matrix(unlist(A), ncol = n, dimnames = list(row, 
        col)))
    vr$mean <- apply(vr, 1, mean)
    vr$var <- apply(vr, 1, var)
    vr$cv <- vr$var^0.5/vr$mean * 100
    meanA <- matrix(vr$mean, nrow = x)
    eigA <- eigen.analysis(meanA)
    vr$sens <- as.vector(eigA$sensitivities)
    vr$elas <- as.vector(eigA$elasticities)

###################################################################################
#Comando utilizado para adicionar o tipo de taxa vital referente ao elemento da matrix
###################################################################################

	a=rep(rep(1:x,1),x)
	vr$VR<- VitalRate<-as.factor(as.vector((ifelse(a == 1, "Fertility", "Survival"))))

###################################################################################
							#Aqui é necessário retirar os valores 
    							#referentes à medias zero pois elas 
  vr1 <- subset(vr, mean > 0 & var > 0)	#Se tornarão valores NaN e interrogações 
							#no resultado das análises de correlação 
							#e no gráfico gerado pela função corrplot
######################################################################################
 					#Aqui se encontra a principal modificação desse script
					#Na versão original a funçao pfister.plot retorna apenas os 
	   vr1[, 1:n ]		# valores já computados de média, elasticidade e sensibilidade
					# Aqui eu consegui fazer retornar em forma de lista as taxas
					# vitais referente à todas as observações na matrix	
######################################################################################
}

#================================================================================
###  ==== CRIANDO MATRIZES ALEATÓRIAS COM DISTRIBUIÇÃO BETA A PARTIR    ====  ###
#	==		DA MÉDIA E DA VARIÂNCIA DE CADA CLASSE ETÁRIA		==	  #
#================================================================================


BNM<-BetaNullModel<-function(A){
vr<-PLOTDOC(A)
	vr$mean <- apply(vr, 1, mean)
	vr$var <- apply(vr, 1, var)
	vr$cv <- vr$var^0.5/vr$mean * 100

vr1 <- subset(vr, mean > 0 & var > 0)

vr1$beta<-replicate(length(n),(vr1$xusbus<-as.vector(matrix(numeric(1*length(vr1$mean))))))
	for (i in 1:(length(vr1$mean)*length(n)))
{
	vr1$beta[i]<-betaval(vr1$mean,sqrt(vr1$var))
} 

return(as.data.frame(vr1$beta,row.names=row.names(vr1)))

}




#================================================================================
##      ==== CRIANDO A FUNÇÃO PARA DETECTAR VALORES CRÍTICOS ====            ####
#================================================================================
# Função retirada na integra do link: XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#================================================================================

d<-cbind(rnorm(1:100),rnorm(1:10000))

critical.r <- function( n, alpha = .05 ) {
  df <- n - 2
  critical.t <- qt(alpha/2, df, lower.tail = F)
  critical.r <- sqrt( (critical.t^2) / ( (critical.t^2) + df ) )
  return(critical.r)
}

critical.r (length(d))



#================================================================================
# Função Para comparar a frequência de valores significativos
#================================================================================

ComparedPairs <-function(A,n){
#Sendo A uma matriz de correlação qualquer no presente estudo representada por cor(t(PLOTDOC(n)))
#e n o número de matrizes existentes que irá representar o tempo
	A[upper.tri(A)] <- Inf
	b<-critical.r(length(n))
return(list(
#Tamanho total da matriz de correlação (elementos em "cor(cor(t(PLOTDOC(n))"
	"Correlation matrix length" =	length(A),
#Número de elemetos na porção triangular superior (sem a diagonal =1)
	"Upper length" = 	sum(A == Inf),
#Valor crítico
	"Critical value of significance" = b,
#Elementos SIGNITIVCATIVOS (+ & -)
	"Signitive values"= sum(A >= b & A <1) + sum(A <= (b*-1)),
#Elementos ACIMA acima do valor crítico de significancia
	"POSITIVE significative values" = sum(A >= b & A <1),
#Elementos ACIMA acima do valor crítico de significancia
	"NEGATIVE significative values" = sum(A <= (b*-1)),
#Elementos não significativos
	"NOT Sinigficative (p > 0.05)" = sum(A > (b*-1) & A < b)))
}




#==================================================================================
#USAGE:   ComparedPairs (cor(t(PLOTDOC(n))),n)
#==================================================================================


#================================================================================
###############                 ==== ANALISE POR TRABALHO ====         	#########
#================================================================================
#### PRESTAR ATENÇÂO QUE TODOS VIRARÃO n PARA OTIMIZAR AS ANÁLISES

#				n

####

##############################################################################
###BENGTSSON1993
##############################################################################

n<-bengtsson1993<-lapply(list(mylist[[3]],mylist[[4]],mylist[[5]],mylist[[6]],mylist[[7]],mylist[[8]]),as.matrix)
bengtsson1993

#Gráfico de correlações comparativo
par(mfrow=c(2,2))

#Matriz de correlações
corrplot(cor(t(PLOTDOC(n))),type = "upper")
corrplot(cor(t(BNM(n))),type="upper")

#Histogramas de distribuição de valores de coeficiente de correlação (rho)
#Dados originais
hist((cor(t(PLOTDOC(n))))
[!cor(PLOTDOC(n))==1]
,main="Correlation with alpha(0.05) significance level",xlab="rho")

abline(v=(critical.r(length(n))),lty=2,col="red")
abline(v=(critical.r(length(n))*-1),lty=2,col="red")

#Modelo nulo
hist((cor(t(PLOTDOC(n))))
[!cor(PLOTDOC(n))==1]
,main="Null Model Matrix correlation",xlab="rho")


#Comparando a proporção de valores significativos positivos, signif. negativos e não significativos
ComparedPairs (cor(t(PLOTDOC(n))),n)
ComparedPairs (cor(t(BNM(n))),n)

#Comparando se a proporção de valores significativos é estatisticamente diferente do esperado ao caso
#Segundo as correlações geradas a partir de uma matriz criada aleatóriamente


chisq.test(matrix(cbind(ComparedPairs (cor(t(PLOTDOC(n))),n)$"Signitive values",
ComparedPairs (cor(t(BNM(n))),n)$"Signitive values",
ComparedPairs (cor(t(PLOTDOC(n))),n)$"NOT Sinigficative (p > 0.05)",
ComparedPairs (cor(t(BNM(n))),n)$"NOT Sinigficative (p > 0.05)"),2,2))

a





##############################################################################
###BIERZYCHUDEK1982
##############################################################################

####BROOKTONDALE & FALL CREEK

n<-Bierzychudek1982Brook<-lapply(list(mylist[[9]],mylist[[10]],mylist[[11]],mylist[[12]]),as.matrix)
Bierzychudek1982Brook

par(mfrow=c(2,1))

cor(t(PLOTDOC(n)))

corrplot(cor(t(PLOTDOC(n))),type = "upper")

hist((cor(t(PLOTDOC(n))))
[!cor(PLOTDOC(n))==1]
,main="Correlation with alpha(0.05) significance level",xlab="rho")

abline(v=(critical.r(length(n))),lty=2,col="red")
abline(v=(critical.r(length(n))*-1),lty=2,col="red")




##############################################################################
###CHARRON E GAGNON 1991
##############################################################################

###TOTAL    Todas as populações amostradas

n<-CG1991<-lapply(list(mylist[[13]],mylist[[14]],mylist[[15]],mylist[[16]],mylist[[17]],mylist[[18]],mylist[[19]],mylist[[20]]),as.matrix)
CG1991

par(mfrow=c(2,1))

cor(t(PLOTDOC(n)))

corrplot(cor(t(PLOTDOC(n))),type = "upper")

hist((cor(t(PLOTDOC(n))))
[!cor(PLOTDOC(n))==1]
,main="Correlation with alpha(0.05) significance level",xlab="rho")

abline(v=(critical.r(length(n))),lty=2,col="red")
abline(v=(critical.r(length(n))*-1),lty=2,col="red")

ComparedPairs (cor(t(PLOTDOC(n))),n)


##############################################################################
###DOAK ET AL 1994
##############################################################################

n<-doak1994<-lapply(list(mylist[[21]],mylist[[22]],mylist[[23]],mylist[[24]]),as.matrix)
doak1994

par(mfrow=c(2,1))

cor(t(PLOTDOC(n)))

corrplot(cor(t(PLOTDOC(n))),type = "upper")

hist((cor(t(PLOTDOC(n))))
[!cor(PLOTDOC(n))==1]
,main="Correlation with alpha(0.05) significance level",xlab="rho")

abline(v=(critical.r(length(n))),lty=2,col="red")
abline(v=(critical.r(length(n))*-1),lty=2,col="red")

ComparedPairs (cor(t(PLOTDOC(n))),n)

##############################################################################
###FIELDER 1987
##############################################################################

#COLOCHORTUS ALBUS
n<-FilderCA<-lapply(list(mylist[[27]],mylist[[28]]),as.matrix)
FilderCA

par(mfrow=c(2,1))

cor(t(PLOTDOC(n)))

corrplot(cor(t(PLOTDOC(n))),type = "upper")

hist((cor(t(PLOTDOC(n))))
[!cor(PLOTDOC(n))==1]
,main="Correlation with alpha(0.05) significance level",xlab="rho")

abline(v=(critical.r(length(n))),lty=2,col="red")
abline(v=(critical.r(length(n))*-1),lty=2,col="red")

#COLOCHORTUS OBIPONENSIS
n<-FilderCO<-lapply(list(mylist[[29]],mylist[[30]]),as.matrix)
FilderCO

par(mfrow=c(2,1))

cor(t(PLOTDOC(n)))

corrplot(cor(t(PLOTDOC(n))),type = "upper")

hist((cor(t(PLOTDOC(n))))
[!cor(PLOTDOC(n))==1]
,main="Correlation with alpha(0.05) significance level",xlab="rho")

abline(v=(critical.r(length(n))),lty=2,col="red")
abline(v=(critical.r(length(n))*-1),lty=2,col="red")

ComparedPairs (cor(t(PLOTDOC(n))),n)

#COLOCHORTUS PULCHELLUS 
n<-fielderCP<-lapply(list(mylist[[31]],mylist[[32]]),as.matrix)
fielderCP

par(mfrow=c(2,1))

cor(t(PLOTDOC(n)))

corrplot(cor(t(PLOTDOC(n))),type = "upper")

hist((cor(t(PLOTDOC(n))))
[!cor(PLOTDOC(n))==1]
,main="Correlation with alpha(0.05) significance level",xlab="rho")

abline(v=(critical.r(length(n))),lty=2,col="red")
abline(v=(critical.r(length(n))*-1),lty=2,col="red")

ComparedPairs (cor(t(PLOTDOC(n))),n)


#COLOCHORTUS TIBURONENSIS
n<-fielderCT<-lapply(list(mylist[[33]],mylist[[34]]),as.matrix)
fielderCT

par(mfrow=c(2,1))

cor(t(PLOTDOC(n)))

corrplot(cor(t(PLOTDOC(n))),type = "upper")

hist((cor(t(PLOTDOC(n))))
[!cor(PLOTDOC(n))==1]
,main="Correlation with alpha(0.05) significance level",xlab="rho")

abline(v=(critical.r(length(n))),lty=2,col="red")
abline(v=(critical.r(length(n))*-1),lty=2,col="red")

ComparedPairs (cor(t(PLOTDOC(n))),n)


##############################################################################
###HORVITZ E SCHEMSKE 1995
##############################################################################


###PLOT1

n<-HS1995plot1<-lapply(list(mylist[[36]],mylist[[40]],mylist[[44]],mylist[[48]]),as.matrix)
HS1995plot1

par(mfrow=c(2,1))

cor(t(PLOTDOC(n)))

corrplot(cor(t(PLOTDOC(n))),type = "upper")

hist((cor(t(PLOTDOC(n))))
[!cor(PLOTDOC(n))==1]
,main="Correlation with alpha(0.05) significance level",xlab="rho")

abline(v=(critical.r(length(n))),lty=2,col="red")
abline(v=(critical.r(length(n))*-1),lty=2,col="red")

###PLOT2

n<-HS1995plot2<-lapply(list(mylist[[37]],mylist[[41]],mylist[[45]],mylist[[49]]),as.matrix)
HS1995plot2

par(mfrow=c(2,1))

cor(t(PLOTDOC(n)))

corrplot(cor(t(PLOTDOC(n))),type = "upper")

hist((cor(t(PLOTDOC(n))))
[!cor(PLOTDOC(n))==1]
,main="Correlation with alpha(0.05) significance level",xlab="rho")

abline(v=(critical.r(length(n))),lty=2,col="red")
abline(v=(critical.r(length(n))*-1),lty=2,col="red")

###PLOT3

n<-HS1995plot3<-lapply(list(mylist[[38]],mylist[[42]],mylist[[46]],mylist[[50]]),as.matrix)
HS1995plot3

par(mfrow=c(2,1))

cor(t(PLOTDOC(n)))

corrplot(cor(t(PLOTDOC(n))),type = "upper")

hist((cor(t(PLOTDOC(n))))
[!cor(PLOTDOC(n))==1]
,main="Correlation with alpha(0.05) significance level",xlab="rho")

abline(v=(critical.r(length(n))),lty=2,col="red")
abline(v=(critical.r(length(n))*-1),lty=2,col="red")

###PLOT4

n<-HS1995plot4<-lapply(list(mylist[[39]],mylist[[43]],mylist[[47]],mylist[[51]]),as.matrix)
HS1995plot4

par(mfrow=c(2,1))

cor(t(PLOTDOC(n)))

corrplot(cor(t(PLOTDOC(n))),type = "upper")

hist((cor(t(PLOTDOC(n))))
[!cor(PLOTDOC(n))==1]
,main="Correlation with alpha(0.05) significance level",xlab="rho")

abline(v=(critical.r(length(n))),lty=2,col="red")
abline(v=(critical.r(length(n))*-1),lty=2,col="red")


###TODOS OS PLOTS!

n<-HS1995plot4<-lapply(list(mylist[[36]],
mylist[[37]],
mylist[[38]],
mylist[[39]],
mylist[[40]],
mylist[[41]],
mylist[[42]],
mylist[[43]],
mylist[[44]],
mylist[[45]],
mylist[[46]],
mylist[[47]],
mylist[[48]],
mylist[[49]],
mylist[[50]],
mylist[[51]]),as.matrix)

par(mfrow=c(2,1))

cor(t(PLOTDOC(n)))

corrplot(cor(t(PLOTDOC(n))),type = "upper")

hist((cor(t(PLOTDOC(n))))
[!cor(PLOTDOC(n))==1]
,main="Correlation with alpha(0.05) significance level",xlab="rho")

abline(v=(critical.r(length(n))),lty=2,col="red")
abline(v=(critical.r(length(n))*-1),lty=2,col="red")

ComparedPairs (cor(t(PLOTDOC(n))),n)



##############################################################################
###HUENKE E MARKS 1987
##############################################################################

###DIVIDE SWAMP

n<-HM1987DS<-lapply(list(mylist[[52]],mylist[[55]],mylist[[56]]),as.matrix)
HM1987DS

par(mfrow=c(2,1))

cor(t(PLOTDOC(n)))

corrplot(cor(t(PLOTDOC(n))),type = "upper")

hist((cor(t(PLOTDOC(n))))
[!cor(PLOTDOC(n))==1]
,main="Correlation with alpha(0.05) significance level",xlab="rho")

abline(v=(critical.r(length(n))),lty=2,col="red")
abline(v=(critical.r(length(n))*-1),lty=2,col="red")

ComparedPairs (cor(t(PLOTDOC(n))),n)

###OLD FIELDER

n<-HM1987OF<-lapply(list(mylist[[53]],mylist[[55]],mylist[[57]]),as.matrix)
HM1987OF

par(mfrow=c(2,1))

cor(t(PLOTDOC(n)))

corrplot(cor(t(PLOTDOC(n))),type = "upper")

hist((cor(t(PLOTDOC(n))))
[!cor(PLOTDOC(n))==1]
,main="Correlation with alpha(0.05) significance level",xlab="rho")

abline(v=(critical.r(length(n))),lty=2,col="red")
abline(v=(critical.r(length(n))*-1),lty=2,col="red")

ComparedPairs (cor(t(PLOTDOC(n))),n)


###TODOS!

n<-HM1987TODOS<-lapply(list(
mylist[[52]],
mylist[[53]],
mylist[[54]],
mylist[[55]],
mylist[[56]],
mylist[[57]]),as.matrix)

HM1987TODOS

par(mfrow=c(2,1))

cor(t(PLOTDOC(n)))

corrplot(cor(t(PLOTDOC(n))),type = "upper")

hist((cor(t(PLOTDOC(n))))
[!cor(PLOTDOC(n))==1]
,main="Correlation with alpha(0.05) significance level",xlab="rho")

abline(v=(critical.r(length(n))),lty=2,col="red")
abline(v=(critical.r(length(n))*-1),lty=2,col="red")

ComparedPairs (cor(t(PLOTDOC(n))),n)

##############################################################################
###HUGHES 1984
##############################################################################

n<-hughes1984<-lapply(list(mylist[[58]],mylist[[59]]),as.matrix)
hughes1984

par(mfrow=c(2,1))

cor(t(PLOTDOC(n)))

corrplot(cor(t(PLOTDOC(n))),type = "upper")

hist((cor(t(PLOTDOC(n))))
[!cor(PLOTDOC(n))==1]
,main="Correlation with alpha(0.05) significance level",xlab="rho")

abline(v=(critical.r(length(n))),lty=2,col="red")
abline(v=(critical.r(length(n))*-1),lty=2,col="red")

ComparedPairs (cor(t(PLOTDOC(n))),n)


##############################################################################
###MENGES 1990
##############################################################################

n<-menges1990<-lapply(list(mylist[[60]],mylist[[61]],mylist[[62]]),as.matrix)
menges1990

par(mfrow=c(2,1))

cor(t(PLOTDOC(n)))

corrplot(cor(t(PLOTDOC(n))),type = "upper")

hist((cor(t(PLOTDOC(n))))
[!cor(PLOTDOC(n))==1]
,main="Correlation with alpha(0.05) significance level",xlab="rho")

abline(v=(critical.r(length(n))),lty=2,col="red")
abline(v=(critical.r(length(n))*-1),lty=2,col="red")


ComparedPairs (cor(t(PLOTDOC(n))),n)


##############################################################################
####MOLONEY1988
##############################################################################

### TODAS AS LOCALIDADES!

n<-Moloney1988<-lapply(
list(mylist[[63]],
mylist[[64]],
mylist[[65]],
mylist[[66]],
mylist[[67]],
mylist[[68]],
mylist[[69]],
mylist[[70]],
mylist[[71]],
mylist[[72]]),as.matrix)
Moloney1988

par(mfrow=c(2,1))

cor(t(PLOTDOC(n)))

corrplot(cor(t(PLOTDOC(n))),type = "upper")

hist((cor(t(PLOTDOC(n))))
[!cor(PLOTDOC(n))==1]
,main="Correlation with alpha(0.05) significance level",xlab="rho")

abline(v=(critical.r(length(n))),lty=2,col="red")
abline(v=(critical.r(length(n))*-1),lty=2,col="red")


ComparedPairs (cor(t(PLOTDOC(n))),n)


##############################################################################
#### MENGES 1990
##############################################################################

### TODAS AS LOCALIDADES!

n<-Menges1990<-lapply(list(mylist[[60]],mylist[[61]],mylist[[62]]),as.matrix)
Menges1990
n

Menges1990

par(mfrow=c(2,1))

cor(t(PLOTDOC(n)))

corrplot(cor(t(PLOTDOC(n))),type = "upper")

hist((cor(t(PLOTDOC(n))))
[!cor(PLOTDOC(n))==1]
,main="Correlation with alpha(0.05) significance level",xlab="rho")

abline(v=(critical.r(length(n))),lty=2,col="red")
abline(v=(critical.r(length(n))*-1),lty=2,col="red")


ComparedPairs (cor(t(PLOTDOC(n))),n)

##############################################################################
#### Nault e Gagnon 1993
##############################################################################


n<-NG1993<-lapply(
list(mylist[[73]],
mylist[[74]],
mylist[[75]],
mylist[[76]]),as.matrix)
NG1993
n

par(mfrow=c(2,1))

cor(t(PLOTDOC(n)))

corrplot(cor(t(PLOTDOC(n))),type = "upper")

hist((cor(t(PLOTDOC(n))))
[!cor(PLOTDOC(n))==1]
,main="Correlation with alpha(0.05) significance level",xlab="rho")

abline(v=(critical.r(length(n))),lty=2,col="red")
abline(v=(critical.r(length(n))*-1),lty=2,col="red")


ComparedPairs (cor(t(PLOTDOC(n))),n)


##############################################################################
#### Nault e Gagnon 1993
##############################################################################


n<-NG1993<-lapply(list(mylist[[73]],mylist[[74]],mylist[[75]],mylist[[76]]),as.matrix)
NG1993
n

par(mfrow=c(2,1))

cor(t(PLOTDOC(n)))

corrplot(cor(t(PLOTDOC(n))),type = "upper")

hist((cor(t(PLOTDOC(n))))
[!cor(PLOTDOC(n))==1]
,main="Correlation with alpha(0.05) significance level",xlab="rho")

abline(v=(critical.r(length(n))),lty=2,col="red")
abline(v=(critical.r(length(n))*-1),lty=2,col="red")


ComparedPairs (cor(t(PLOTDOC(n))),n)

##############################################################################
#### McFadden 1967
##############################################################################


n<-McFadden1967<-lapply(list(mylist[[77]],
mylist[[78]],
mylist[[79]],
mylist[[80]],
mylist[[81]],
mylist[[82]],
mylist[[83]],
mylist[[84]],
mylist[[85]],
mylist[[86]],
mylist[[87]]),as.matrix)


McFadden1967

n

par(mfrow=c(2,1))

cor(t(PLOTDOC(n)))

corrplot(cor(t(PLOTDOC(n))),type = "upper")

hist((cor(t(PLOTDOC(n))))
[!cor(PLOTDOC(n))==1]
,main="Correlation with alpha(0.05) significance level",xlab="rho")

abline(v=(critical.r(length(n))),lty=2,col="red")
abline(v=(critical.r(length(n))*-1),lty=2,col="red")


ComparedPairs (cor(t(PLOTDOC(n))),n)
