# Experimental analysis of Kendall correction to sampling variation

# In this example, I run Kendall functions of popbio package to estimate Kendall's mean and Kendall's unbiased variance.
# The experimental approach of this variations were compare differences between absolute values of individuals with percentage or individuals analysed

# Data analysis was performed with Hudosnia montana plant available in Box 8.5 of Morris e Doak 2002. Quantitative Conservation Biology pag. 275 

# Other example used were Desert Tortoise and are  available in popbio package with "tor"


#Let's go

library(popbio)

#Exemplo do popbio
## desert tortoise input from Box 8.2 - compare results to Table 8.3
tor<-data.frame(rate=rep(c("g4","g5","g6"),each=3),
year=rep(1:3,3),      ## representing 70s, early 80s, late 80s
start=c(17,15,7,22,19,4,32,31,10),
grow=c(8,1,0,5,5,0,2,1,0))
## use fewer grades for faster loop
tor.est<-Kendall(tor, grades=200)
tor.est



# EXEMPLO Hudsonia montana
# Os objetos foram transformados com a função "dput" para serem exportado para cá


#OBJETO 1 - Apenas copiar e colar, relaxa!
bla<- structure(list(rate = structure(c(1L, 1L, 1L, 1L, NA, NA, NA, 
NA, NA), .Label = c("g4", "g5", "g6"), class = "factor"), year = c(1, 
2, 3, 4, NA, NA, NA, NA, NA), start = c(25, 19, 10, 13, NA, NA, 
NA, NA, NA), grow = c(4, 7, 2, 3, 0, 0, 0, 0, 0)), .Names = c("rate", 
"year", "start", "grow"), row.names = c(NA, 9L), class = "data.frame")

bla

# O objeto "bla" segue a lógica proposta para as correções por amostragem proposta por Morris e Doak 2002. 
# No entanto, nosso objetivo aqui é comparar se a diferença entre as correções de Kendall com valores absolutos e os valores em porcentagem
# São realmente muito diferentes.

#Para isso vamos trabalhar o objeto bla2 como sendo a porcentagem (vezes 100) 
bla2<-structure(list(rate = structure(c(1L, 1L, 1L, 1L, NA, NA, NA, 
NA, NA), .Label = c("g4", "g5", "g6"), class = "factor"), year = c(1, 
2, 3, 4, NA, NA, NA, NA, NA), start = c(100, 100, 100, 100, NA, 
NA, NA, NA, NA), grow = c(16, 37, 20, 23, 0, 0, 0, 0, 0)), .Names = c("rate", 
"year", "start", "grow"), row.names = c(NA, 9L), class = "data.frame")

#Agora no objeto "bla3" vamos trabalhar com a porcentagem de 0-1,uma vez que sugere-se que a análise desenvolvida para variáveis que variam
# de 0-1 (não sei onde li isso, conferir).

bla3<-structure(list(rate = structure(c(1L, 1L, 1L, 1L, NA, NA, NA, 
NA, NA), .Label = c("g4", "g5", "g6"), class = "factor"), year = c(1, 
2, 3, 4, NA, NA, NA, NA, NA), start = c(1, 1, 1, 1, NA, NA, NA, 
NA, NA), grow = c(0.16, 0.37, 0.2, 0.23, 0, 0, 0, 0, 0)), .Names = c("rate", 
"year", "start", "grow"), row.names = c(NA, 9L), class = "data.frame")



set.seed(1)

Kendall(bla[1:4,],500)
Kendall(bla2[1:4,],500)
Kendall(bla3[1:4,],500)


#Notar as diferenças entre eles, bla 2 parece diferir bastante, porém bla e bla 3 são muito semelhantes, diferindo apenas na 3ª casa decimal, geralmente.
# Quando comparamos os valores obtidos por Morris e Doak nos box 8.3 e 8.5 para Desert tortoise e Hudsonia montana, respectivamente, não teriamos motivo para 
# Acreditar que bla e bla3 são significativamente diferente entre si, mostrando que pode ser adequado utilizarmos as correções de Kendall para as proporções das "vital rates"

