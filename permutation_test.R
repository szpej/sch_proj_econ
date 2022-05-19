library(tidymodels)
library(gtools)

#################################    Symulacja    ##############################
my_data <- 
  tibble(
    value = c(27,33,37,52,53,57,69,70,71,77,
              6,9,14,16,29,43,45,47,50,55),
    grupa = rep(c("pierwsza", "druga"), each = 10),
    ranga = c(5,7,8,13,14,16,17,18,19,20,
              1,2,3,4,6,9,10,11,12,15)
  )

benchmark <-
  mean(my_data$ranga[which(my_data$grupa == "pierwsza")])-
  mean(my_data$ranga[which(my_data$grupa == "druga")])

M <- 1000000
wynik <- numeric(1000000)
for (i in 1:M){
  x <- sample(my_data$ranga, 20)
  wynik[i] <- 
    abs(mean(x[1:10]) - mean(x[11:20]))
}
sum(wynik >= benchmark)/1000000

#################################    exact p value    ##########################
my_data <- 
  tibble(
    value = c(27,33,37,52,53,57,69,70,71,77,
              6,9,14,16,29,43,45,47,50,55),
    grupa = rep(c("pierwsza", "druga"), each = 10),
    ranga = c(5,7,8,13,14,16,17,18,19,20,
              1,2,3,4,6,9,10,11,12,15)
  )

mianownik <- choose(20,10)
benchmark <-
  mean(my_data$ranga[which(my_data$grupa == "pierwsza")])-
  mean(my_data$ranga[which(my_data$grupa == "druga")])

list_comb <- combn(my_data$ranga,10)
list_comb <- t(list_comb)
wynik <- numeric(mianownik)
for (i in 1:mianownik) {
  grupa_1 <- list_comb[i,]
  grupa_2 <- my_data$ranga[!(my_data$ranga %in% grupa_1)]
  wynik[i] <-  abs(mean(grupa_1)- mean(grupa_2))
}

mean(wynik >= benchmark)

#################################    exact p value    ##########################
mianownik <- choose(20,10)

benchmark <- mean(my_data$ranga[which(my_data$grupa == "pierwsza")])-
mean(my_data$ranga[which(my_data$grupa == "druga")])

list_comb <- permutations(20,20,my_data$ranga,repeats.allowed = FALSE)
list_comb <- t(list_comb)

wynik <- numeric(184756)

mean(list_comb[1,1:5])
mean(list_comb[1,6:10])
list_comb[1,]

for (i in 1:184756) {
  wynik[i] <- 
    abs(mean(list_comb[i,1:5]) - mean(list_comb[i,6:10]))
}

(length(wynik[which(wynik >= benchmark)])-1)/mianownik

wynik %>% view()


for (i in 1:184756) {
  wynik[i] <- 
    t.test(list_comb[i,1:5], list_comb[i,6:10], 
           alternative = "two.sided", var.equal = TRUE)$statistic
}
pt(-2.21,6, lower = TRUE)*2
pt(2.21,6,lower = FALSE) *2
?t.test

t_true <- t.test(
  my_data$ranga[which(my_data$grupa == "pierwsza")],
  my_data$ranga[which(my_data$grupa == "druga")],
  alternative = "two.sided",
  var.equal = TRUE)$statistic

t_false <- t.test(
  my_data$ranga[which(my_data$grupa == "pierwsza")],
  my_data$ranga[which(my_data$grupa == "druga")],
  alternative = "two.sided",
  var.equal = FALSE)$statistic
