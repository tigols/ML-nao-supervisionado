# Atividade 1 ME921 #

rm(list = ls(all.names = TRUE))
par(mar = c(2.5,2.5,2.5,2.5) + 1.5)

#Bibliotecas
library(readr)
library(tidyverse)
library(corrplot)

# Funcoes #
panel.hist <- function(x, ...)
{
  usr <- par("usr")
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "blue", ...)
}

#Web scraping feito do seguinte site usando a extensao "Web Scraper"
url = "https://www.transfermarkt.com.br/campeonato-brasileiro-serie-a/torschuetzenliste/wettbewerb/BRA1/saison_id/2022/altersklasse/alle/detailpos//plus/1"

dados = read_csv("transfermarkt.csv")

# Tratamento dos dados #

dados = dados[-c(1:25),-c(1,2,3)] #primeiras linhas estao repetidas

dados$`minutos em campo` = gsub("[\\.']","",dados$`minutos em campo`)#retirar pontuacoes
dados$`minutos para gol` = gsub("'","",dados$`minutos para gol`)#retirar pontuacao de miutagem
dados$`gols por jogo` = gsub(",",".",dados$`gols por jogo`)#trocar , por .
dados$idade = str_sub(dados$idade, 1,2)#pegar apenas a idade que o jogador tinha na epoca e nao idade atual

dados = dados %>% 
  mutate_at(vars(idade,jogos,assistencias,penalti,`minutos em campo`,`minutos para gol`,`gols por jogo`,gols), as.numeric)

dados = dados %>% mutate_at(vars(posicao), as.factor)
# definindo tipo das variaveis


# Analise Descritiva #
pairs(dados[,c(4,5,10)],
       diag.panel = panel.hist,cex = 1.5)


# Cluster Hierarquico #

#Complete Linkage

analise_centroavante = dados[dados$posicao == "Centroavante",c(10,4)]
# analise_centroavante = as.data.frame(apply(analise_centroavante, 2, scale))


dist(analise_centroavante)
plot.new()
plot(analise_centroavante, pch = 20)

cluster_CA = hclust(dist(analise_centroavante),method = "complete")

# Definindo numero de clusters (altura do corte) (Elbow)

n_CA = nrow(analise_centroavante)
h_teste = seq(0, 40, by = 1)
k_CA = numeric(length(h_teste))
y_CA = with(dados[dados$posicao == "Centroavante",], cbind(jogos, gols))
total_var = numeric(length(h_teste))

for( i in seq_along(h_teste)){
  groups = factor(cutree(cluster_CA, h = h_teste[i]))
  k_CA[i] = length(levels(groups))
  B = summary(manova(y_CA ~ groups), tol = 0)$SS$groups
  W = summary(manova(y_CA ~ groups), tol = 0)$SS$Residuals
  total_var[i] = det(B)/det(B+W)
}

plot(total_var)
abline(v = 13, col = 'red')

h_teste[13]

plot(cluster_CA)
abline(h = h_teste[13], col = "Red", lty = 2)

grupos = cutree(cluster_CA, h = h_teste[13])

complete_CA = cbind(dados[dados$posicao == "Centroavante",],grupos)
complete_CA[,c(10,4)] = analise_centroavante
complete_CA = complete_CA %>% rename(cluster = grupos)
complete_CA = complete_CA %>% group_by(cluster) %>%
  mutate( media_assist = as.factor(round(mean(assistencias),2)) ) 
complete_CA = complete_CA %>% group_by(cluster) %>%
                    mutate( centro_x = mean(gols), centro_y = mean(jogos))
complete_CA = complete_CA %>% arrange(cluster, desc(`gols por jogo`))

plot(analise_centroavante, col = grupos, pch = 20)

for(i in unique(grupos)){
  objeto = analise_centroavante[grupos == i,]
  xts = chull(objeto)
  xts = c(xts,xts[1])
  lines(objeto[xts,], col = i)
}

for(i in unique(grupos)){
  text(complete_CA$centro_x[grupos == i] - 0.1, complete_CA$centro_y[grupos == i], labels = complete_CA$media_assist[grupos == i])
}

# Zagueiros
analise_zag = dados[dados$posicao == "Zagueiro",c(10,4)]

dist(analise_zag)
plot.new()
plot(analise_zag, pch = 20)

cluster_ZAG = hclust(dist(analise_zag),method = "complete")

# Definindo numero de clusters (altura do corte) (Elbow)

n_ZAG = nrow(analise_zag)
h_teste = seq(0, 35, by = 0.5)
k_ZAG = numeric(length(h_teste))
y_ZAG = with(dados[dados$posicao == "Zagueiro",], cbind(jogos, gols))
total_var = numeric(length(h_teste))

for( i in seq_along(h_teste)){
  groups = factor(cutree(cluster_ZAG, h = h_teste[i]))
  k_ZAG[i] = length(levels(groups))
  B = summary(manova(y_ZAG ~ groups), tol = 0)$SS$groups
  W = summary(manova(y_ZAG ~ groups), tol = 0)$SS$Residuals
  total_var[i] = det(B)/det(B+W)
}

plot(total_var)
abline(v = 14, col = 'red')
h_teste[14]

plot(cluster_ZAG)
abline(h = h_teste[14], col = "Red", lty = 2)
grupos = cutree(cluster_ZAG, h = h_teste[14])

complete_ZAG = cbind(dados[dados$posicao == "Zagueiro",],grupos)
complete_ZAG[,c(10,4)] = analise_zag

complete_ZAG = complete_ZAG %>% rename(cluster = grupos)
complete_ZAG = complete_ZAG %>% group_by(cluster) %>%
  mutate( media_assist = as.factor(round(mean(assistencias),2)) ) 
complete_ZAG = complete_ZAG %>% group_by(cluster) %>%
  mutate( centro_x = mean(gols), centro_y = mean(jogos))
complete_ZAG = complete_ZAG %>% arrange(cluster, desc(`gols por jogo`))

plot(analise_zag, col = grupos, xlim= c(0,5), pch = 20)

for(i in unique(grupos)){
  objeto = analise_zag[grupos == i,]
  xts = chull(objeto)
  xts = c(xts,xts[1])
  lines(objeto[xts,], col = i)
}

for(i in unique(grupos)){
  text(complete_ZAG$centro_x[grupos == i] - 0.25, complete_ZAG$centro_y[grupos == i], labels = complete_ZAG$media_assist[grupos == i])
}

#Meias
analise_mei = dados[dados$posicao == "Meia Ofensivo",c(10,4)]

dist(analise_mei)
plot.new()
plot(analise_mei, pch = 20)

cluster_MEI = hclust(dist(analise_mei),method = "complete")

# Definindo numero de clusters (altura do corte) (Elbow)

n_MEI = nrow(analise_mei)
h_teste = seq(0.05, 30, by = 0.5)
k_MEI = numeric(length(h_teste))
y_MEI = with(dados[dados$posicao == "Meia Ofensivo",], cbind(jogos, gols))
total_var = numeric(length(h_teste))

for( i in seq_along(h_teste)){
  groups = factor(cutree(cluster_MEI, h = h_teste[i]))
  k_MEI[i] = length(levels(groups))
  B = summary(manova(y_MEI ~ groups), tol = 0)$SS$groups
  W = summary(manova(y_MEI ~ groups), tol = 0)$SS$Residuals
  total_var[i] = det(B)/det(B+W)
}

plot(total_var)
abline(v = 15, col = 'red')

h_teste[15]

plot(cluster_MEI)
abline(h = h_teste[15], col = "Red", lty = 2)

grupos = cutree(cluster_MEI, h = h_teste[15])

complete_MEI = cbind(dados[dados$posicao == "Meia Ofensivo",],grupos)
complete_MEI[,c(10,4)] = analise_mei
complete_MEI = complete_MEI %>% rename(cluster = grupos)
complete_MEI = complete_MEI %>% group_by(cluster) %>%
  mutate( media_assist = as.factor(round(mean(assistencias),2)) ) 
complete_MEI = complete_MEI %>% group_by(cluster) %>%
  mutate( centro_x = mean(gols), centro_y = mean(jogos))
complete_MEI = complete_MEI %>% arrange(cluster, desc(`gols por jogo`))

plot(analise_mei, col = grupos, pch = 20)

for(i in unique(grupos)){
  objeto = analise_mei[grupos == i,]
  xts = chull(objeto)
  xts = c(xts,xts[1])
  lines(objeto[xts,], col = i)
}

for(i in unique(grupos)){
  text(complete_MEI$centro_x[grupos == i] - 0.05, complete_MEI$centro_y[grupos == i], labels = complete_MEI$media_assist[grupos == i])
}

###################################################################################################################

analise = dados[,c(10,4)]
# analise_mei = as.data.frame(apply(analise_mei,2 ,scale))

dist(analise)
plot.new()
plot(analise, pch = 20)

cluster = hclust(dist(analise),method = "complete")
plot(cluster)
abline(h = 0.5, col = "Red", lty = 2)

# Definindo numero de clusters (altura do corte) (Elbow)

n = nrow(analise)
h_teste = seq(0.05, 41, by = 0.5)
k = numeric(length(h_teste))
y = with(dados, cbind(jogos, gols))
total_var = numeric(length(h_teste))

for( i in seq_along(h_teste)){
  groups = factor(cutree(cluster, h = h_teste[i]))
  k[i] = length(levels(groups))
  B = summary(manova(y ~ groups), tol = 0)$SS$groups
  W = summary(manova(y ~ groups), tol = 0)$SS$Residuals
  total_var[i] = det(B)/det(B+W)
}

plot(total_var)
axis(1, at = seq(0, length(h_teste), by = 1), labels = seq(0, length(h_teste), by = 1))
abline(v = 16, col = 'red')



grupos = cutree(cluster, h = h_teste[14])

complete = cbind(dados,grupos)
complete[,c(10,4)] = analise
complete = complete %>% rename(cluster = grupos)
complete = complete %>% group_by(cluster) %>%
  mutate( media_assist = as.factor(round(mean(assistencias),2)) ) 
complete = complete %>% group_by(cluster) %>%
  mutate( centro_x = mean(gols), centro_y = mean(jogos))
complete = complete %>% arrange(cluster, desc(`gols por jogo`))

plot(analise, col = grupos, pch = 20)

for(i in unique(grupos)){
  objeto = analise[grupos == i,]
  xts = chull(objeto)
  xts = c(xts,xts[1])
  lines(objeto[xts,], col = i)
}

for(i in unique(grupos)){
  text(complete$centro_x[grupos == i] - 0.05, complete$centro_y[grupos == i], labels = complete$media_assist[grupos == i])
}

