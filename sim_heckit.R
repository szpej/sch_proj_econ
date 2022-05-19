library(tidymodels) # crossing z tidyr, więc tidyverse też wystarczy
library(sampleSelection)
library(foreach)
library(doParallel)
library(stringr)



# Generowanie pre-populacji -- przed ustaleniem docelowej macierzy korelacji
# Y, X są normalne, Z jest gamma
set.seed(123)

y_1 <- rnorm(1000, 20, 5)
x_1 <- matrix(rnorm(2000), ncol = 2)
z_1 <- rgamma(1000, 0.7)


d_1 <- cbind(y_1,z_1, x_1)
d_1 <- scale(d_1)
WK <- var(d_1)
W <- solve(chol(WK))
T_1 <- d_1 %*% W

### Stworzenie 10000 możliwości do losowania
range_cor <- seq(-1,1, by = 0.01)[seq(-1,1, by = 0.01) != 0]


list_cor <- vector(mode = "list", length = 10)
for (i in 1:10) {
  list_cor[[i]] <- range_cor[(i*20-19):(20*i)]
}


name_cor <- levels(cut(range_cor, 10))
comb_cor <- crossing(yz = name_cor, yx =name_cor,  zx1 = name_cor, zx2 = name_cor)


### Pętla tworzaca populacje o korelacji z danego przedziału
# Zakładam, że korelacja pomiędzy X jest na stałym poziomie 0.21


losowanie <- function(x) {
  if (x == "(-1,-0.8]") 
    sample(list_cor[[1]],1)
  else if (x == "(-0.8,-0.6]")
    sample( list_cor[[2]],1)
  else if (x == "(-0.6,-0.4]")
    sample(list_cor[[3]],1)
  else if (x == "(-0.4,-0.2]")
    sample(list_cor[[4]],1)
  else if (x == "(-0.2,0]")
    sample(list_cor[[5]],1)
  else if (x == "(0,0.2]")
    sample(list_cor[[6]],1)
  else if (x == "(0.2,0.4]" )
    sample(list_cor[[7]],1)
  else if (x == "(0.4,0.6]")
    sample(list_cor[[8]],1)
  else if (x == "(0.6,0.8]")
    sample(list_cor[[9]],1)
  else 
    sample(list_cor[[10]],1)
}

# generowanie populacji - nie więcej niż 20 minut przy 2.9
population_list <- vector(mode = "list", 10000)
proportion <- tibble(id = 1:10000,
                     corYZ = character(10000),
                     corYX = character(10000),
                     corZX1 = character(10000),
                     corZX2 = character(10000)
)
index = 0

for (i in 1:10000) {
  set.seed(i)
  yz <- losowanie(comb_cor[i,1])
  yx1 <- losowanie(comb_cor[i,2])
  yx2 <- losowanie(comb_cor[i,2])
  zx1<- losowanie(comb_cor[i,3])
  zx2<- losowanie(comb_cor[i,4])
  
  C_2 <- matrix(
    c(1, yz, yx1, yx2,
      yz, 1, zx1, zx2,
      yx1, zx1, 1, 0.21,
      yx2, zx2, 0.21, 1),
    ncol = 4) 
  
  # nie zawsze można uzyskać dodatnio określoną macierz
  # więc sprawdzana są wartości własne
  # problem jest taki, że warunek nie może być > 0
  # bo w jednym wypadku przy wartości własnej = 1.332268e-15
  # nie jest możliwa dekompozycja
  # ze względu na to, że któryś minor nie jest większy od zera
  k = 1
  index = index + 1
  while (sum(eigen(C_2)$values > 1.332268e-15) != 4) { 
    set.seed(k + 10000)
    
    yz <- losowanie(comb_cor[i,1])
    yx1 <- losowanie(comb_cor[i,2])
    yx2 <- losowanie(comb_cor[i,2])
    zx1<- losowanie(comb_cor[i,3])
    zx2<- losowanie(comb_cor[i,4])
    xx <- 0.21
    
    C_2 <- matrix(
      c(1, yz, yx1, yx2,
        yz, 1, zx1, zx2,
        yx1, zx1, 1, xx,
        yx2, zx2, xx, 1),
      ncol = 4)  
    
    k = k +1
    if (k >= 50) 
      break
  }
  
  if (k >= 50) {
    index = index - 1
    next 
  }
  
  proportion[index,2] <- comb_cor[i,1]
  proportion[index,3] <- comb_cor[i,2]
  proportion[index,4] <- comb_cor[i,3]
  proportion[index,5] <- comb_cor[i,4]
  
  
  E <- chol(C_2)
  D_2 <- T_1 %*% E
  colnames(D_2) <- c("variableY", "variableZ", "variableX1", "variableX2")
  # prob - prawdopodobieństwo przynależności do big data
  # prob_2 - prawdopodobieństwo selekcji
  D_2 <- as_tibble(D_2) %>% mutate(prob = 1-1/(1+exp(1.5*-variableZ)),
                                   prob_2 = 1-1/(1+exp(-variableZ)))
  
  population_list[[index]] <- D_2
  
} 


# population_list[[4608]] - liczba populacji
table(proportion[,2])
table(proportion[,3])
table(proportion[,4])
table(proportion[,5])


################# obliczenia równoległe


### Do analizy wyników
# największą swartość bezwzględną sumy korelacji ma populacja 4604
# korelacja pomiędzy x1 a x2 jest elemnetem neutralnym - taka sama dla każderj populacji
dist_corr <- foreach(i=1:4608, .combine = rbind) %do% {
  tibble(id = i,
         sum_corr = abs(sum(c(cor(population_list[[i]][,1:4]))[c(2:4,7:8)]))
  )
}


dist_corr <- foreach(i=1:4608, .combine = rbind) %do% {
  y <- cor(population_list[[i]][,1:4])[c(2:4,7:8)]
  tibble(
    id = i,
    yz = y[1],
    yx1 = y[2],
    yx2 = y[3],
    zx1 = y[4],
    zx2 = y[5],
    dist = 0
  )
}


for(i in 1:4608) {
  dist_corr$dist[i] <- 
    dist(rbind(dist_corr[i,2:6], dist_corr[4604,2:6]))
}


numCores <- detectCores()
registerDoParallel(numCores)



wynik <- foreach (i=1:4608,.combine = rbind, .errorhandling = "remove") %dopar% {
  
  test <- population_list[[i]]
  lm_coef <- unname(coef(lm(variableY~variableX1, data = test)))
  real_coef <- lm_coef[2]
  realy <- mean(test$variableY)
  
  bias_tibble <- 
    tibble(naive_x1bias = numeric(100),
           naive_ybias = numeric(100),
           heckit_x1bias = numeric(100),
           heckit_ybias = numeric(100)
    )
  
  for (j in 1:100) {
    set.seed(j)
    test$nielosowa <- rbinom(n = 1000, prob = test$prob, size = 1)
    test$selection <- rbinom(n = 1000, prob = test$prob_2, size = 1)
    x <- selection(selection = selection ~ variableX2,
                   outcome = variableY ~variableX1,
                   method = "ml",
                   subset = which(test$nielosowa == 1),
                   data = test,
                   maxMethod = "NR"
    )
    
    if(!str_detect(x$message, "succes")) {
      x <- selection(selection = selection ~ variableX2,
                     outcome = variableY ~variableX1,
                     method = "ml",
                     subset = which(test$nielosowa == 1),
                     data = test,
                     maxMethod = "BHHH"
      )
    } else if(!str_detect(x$message, "succes")) {
      x <- selection(selection = selection ~ variableX2,
                     outcome = variableY ~variableX1,
                     method = "ml",
                     subset = which(test$nielosowa == 1),
                     data = test,
                     maxMethod = "SANN"
      )
    }
    
    naive_lm <- lm(variableY~variableX1, data = test, subset = which(test$nielosowa == 1))
    
    bias_tibble[j,1] <- real_coef- unname(coef(naive_lm))[2]
    bias_tibble[j,2] <- realy - mean(naive_lm$fitted.values)
    bias_tibble[j,3] <- real_coef - unname(x$estimate[4])
    bias_tibble[j,4] <- realy - mean(test$variableX1 * unname(x$estimate[4]) + unname(x$estimate[3]))
    
  }
  
  tibble(id = proportion[[i,1]],
         corYZ = proportion[[i,2]],
         corYX = proportion[[i,3]],
         corZX1 = proportion[[i,4]],
         corZX2 = proportion[[i,5]],
         distance = dist_corr$dist[i],
         naive_x1bias = mean(bias_tibble$naive_x1bias),
         naive_ybias = mean(bias_tibble$naive_ybias),
         naive_x1var = var(bias_tibble$naive_x1bias),
         naive_yvar = var(bias_tibble$naive_ybias),
         naive_x1mse = naive_x1bias^2 + naive_x1var,
         naive_ymse = naive_ybias^2 + naive_yvar,
         heckit_x1bias = mean(bias_tibble$heckit_x1bias),
         heckit_ybias = mean(bias_tibble$heckit_ybias),
         heckit_x1var = var(bias_tibble$heckit_x1bias),
         heckit_yvar = var(bias_tibble$heckit_ybias),
         heckit_x1mse = heckit_x1bias^2 + heckit_x1var,
         heckit_ymse = heckit_ybias^2 + heckit_yvar
    )
  

  
}
stopImplicitCluster()
write.csv(wynik, "wynik.csv")


cut_wynik <- wynik
cut_wynik <- cut_wynik %>% filter(abs(heckit_x1bias) <= 10)
cut_wynik <- cut_wynik %>% filter(abs(heckit_ybias) <= 100)
cut_wynik <- cut_wynik %>% filter(abs(heckit_x1var) <= 20)
cut_wynik <- cut_wynik %>% filter(abs(heckit_yvar) <= 10)

cut_wynik %>% group_by(corYZ) %>% 
  summarise(naive_x1bias = mean(naive_x1bias),
            heckit_x1bias = mean(heckit_x1bias),
            naive_ybias = mean(naive_ybias),
            heckit_ybias = mean(heckit_ybias),
            naive_yvar = mean(naive_yvar),
            heckit_yvar = mean(heckit_yvar),
            naive_x1var = mean(naive_x1var),
            heckit_x1var = mean(heckit_x1var),
            n = n()
  )

cut_wynik %>% group_by(corYX) %>% 
  summarise(naive_x1bias = mean(naive_x1bias),
            heckit_x1bias = mean(heckit_x1bias),
            naive_ybias = mean(naive_ybias),
            heckit_ybias = mean(heckit_ybias),
            naive_yvar = mean(naive_yvar),
            heckit_yvar = mean(heckit_yvar),
            naive_xvar = mean(naive_x1var),
            heckit_xvar = mean(heckit_x1var),
            n = n()
  )

cut_wynik %>% group_by(corZX1) %>% 
  summarise(naive_x1bias = mean(naive_x1bias),
            heckit_x1bias = mean(heckit_x1bias),
            naive_ybias = mean(naive_ybias),
            heckit_ybias = mean(heckit_ybias),
            naive_yvar = mean(naive_yvar),
            heckit_yvar = mean(heckit_yvar),
            naive_xvar = mean(naive_x1var),
            heckit_xvar = mean(heckit_x1var),
            n = n()
  )

cut_wynik %>% group_by(corZX2) %>% 
  summarise(naive_x1bias = mean(naive_x1bias),
            heckit_x1bias = mean(heckit_x1bias),
            naive_ybias = mean(naive_ybias),
            heckit_ybias = mean(heckit_ybias),
            naive_yvar = mean(naive_yvar),
            heckit_yvar = mean(heckit_yvar),
            naive_xvar = mean(naive_x1var),
            heckit_xvar = mean(heckit_x1var),
            n = n()
  )


cut_wynik %>% 
  filter(abs(heckit_x1bias) <= 1) %>%
  ggplot(aes(x = distance, y = heckit_x1bias, color = corZX2)) + geom_point(size = 1.5) +
  scale_color_brewer(palette="Paired") + theme_bw()
cut_wynik %>% 
  filter(abs(heckit_ybias) <= 10) %>%
  ggplot(aes(x = distance, y = heckit_ybias, color = corZX2)) + geom_point(size = 1.5) +
  scale_color_brewer(palette="Paired") + theme_bw()

cut_wynik[,c(1:6,8,14)] %>%
  filter(abs(heckit_ybias) <= 10) %>%
  pivot_longer(cols = 7:8, names_to = "model") %>%
  ggplot(aes(x = distance, y = value, color = model), data = .) + 
  geom_point(alpha = 0.2)+ theme_bw() +
  scale_color_brewer(palette="Dark2")


cut_wynik[,c(1:6,7,13)]  %>% 
  filter(abs(heckit_x1bias) <= 1) %>%
  pivot_longer(cols = 7:8, names_to = "model") %>%
  ggplot(aes(x = corZX2, y = value, fill = model)) + 
  geom_boxplot(orientation = "x") +
  theme_bw() +
  scale_fill_brewer(palette="Dark2") +
  ylab("x1bias")
cut_wynik[,c(1:6,8,14)]  %>% 
  filter(abs(heckit_ybias) <= 10) %>%
  pivot_longer(cols = 7:8, names_to = "model") %>%
  ggplot(aes(x = corZX2, y = value, fill = model)) + 
  geom_boxplot(orientation = "x") +
  theme_bw() +
  scale_fill_brewer(palette="Dark2") +
  ylab("ybias")
  











