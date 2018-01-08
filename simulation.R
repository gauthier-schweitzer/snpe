
# Simulations -------------------------------------------------------------
# On veut constituer un echantillon de taille N
N = 10000

# On commence par simuler les erreurs
library(mvtnorm)
sigma_mat = matrix(c(1, 0.5, 0.5, 1), nrow=2, ncol=2, byrow = TRUE)
UV = rmvnorm(n = N, mean = c(0,0), sigma = sigma_mat)

# On simule X
X = matrix(data = 1, nrow = N, ncol = 1)

# On definit ensuite les constantes de nos variables latentes
beta = 3
gamma = 0.2
pi = matrix (c(1,0.25), nrow = 1, ncol = 2)

# On assigne les indiv a un groupe de traitement (proportion p) et de controle (1-p)
p = 0.5
D = rbinom(n = N, size = 1, prob = p)

# On introduit les variables latentes et l'observee
Y_hat = D*beta + X*pi[1,1] + UV[,1]
Z_hat = D*gamma + X*pi[1,2] + UV[,2]
Y = 0 + (Z_hat>=0)*Y_hat

# On verifie en comparant les taux de chomage dans les deux groupes
# Chez les non traites
sum((Y == 0)&(D == 0))/length(D ==0)
# Chez les traites
sum((Y == 0)&(D == 1))/length(D ==1)


# Estimation --------------------------------------------------------------
# On estime d'abord la proportion de trimming
p = (sum((Y > 0)&(D == 1))/length(D ==1) - sum((Y > 0)&(D == 0))/length(D ==0))/(sum((Y > 0)&(D == 1))/length(D == 1))

# On estime les quantiles associes
y_p = quantile(x = Y[(Y > 0)&(D == 1)], probs = c(p))
y_1p = quantile(x = Y[(Y > 0)&(D == 1)], probs = c(1-p))

# On calcule la lower bound
lb = mean(Y[(Y>0)&(D==1)&(Y<=y_1p)]) - mean(Y[(Y>0)&(D==0)])
ub = mean(Y[(Y>0)&(D==1)&(Y>=y_p)]) - mean(Y[(Y>0)&(D==0)])


# Comparaison -------------------------------------------------------------
# On compare avec une simple difference de moyenne
mean(Y[(Y>0)&(D==1)]) - mean(Y[(Y>0)&(D==0)])
