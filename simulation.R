library(mvtnorm)
#library(Hmisc)

# Definition des parametres -----------------------------------------------
sigma_mat = matrix(c(0.6, 0.25, 0.25, 1), nrow=2, ncol=2, byrow = TRUE) # Matrice de var-cov des erreurs
beta = 0.2 # Effet traitement sur salaire
gamma = 0.2 # Effet traitement sur employabilite 
pi = matrix (c(9.5,0.25), nrow = 1, ncol = 2) # Constantes pour les 2 eq
p = 0.5  # Proportion de l'echantillon traitee

# Simulations -------------------------------------------------------------
# On veut constituer un echantillon de taille N
N = 100000
set.seed(123)

# On commence par simuler les erreurs
UV = rmvnorm(n = N, mean = c(0,0), sigma = sigma_mat)

# On simule X constant
X = matrix(data = 1, nrow = N, ncol = 1)

# On assigne les indiv a un groupe de traitement (proportion p) et de controle (1-p)
D = rbinom(n = N, size = 1, prob = p)

# On introduit les variables latentes et l'observee
Y_hat = D*beta + X*pi[1,1] + UV[,1]
Z_hat = D*gamma + X*pi[1,2] + UV[,2]
Y = 0 + (Z_hat>=0)*Y_hat

# On introduit S tel que S=1 si salaire observe
S = vector(length = N)
S = 0 + (Y>0)*1

# On verifie en comparant les taux de chomage et salaires dans les deux groupes
# Chez les non traites
sum((Y == 0)&(D == 0))/length(D ==0)
summary(exp(Y[(S==1)&(D==1)]))
# Chez les traites
sum((Y == 0)&(D == 1))/length(D ==1)
summary(exp(Y[(S==1)&(D==0)]))


# Estimation --------------------------------------------------------------
# On estime d'abord la proportion de trimming
p = (sum((S==1)&(D == 1))/length(D ==1) - sum((S==1)&(D == 0))/length(D ==0))/(sum((S==1)&(D == 1))/length(D == 1))

# On estime les quantiles associes
y_p = quantile(x = Y[(S==1)&(D == 1)], probs = c(p))
y_1p = quantile(x = Y[(S==1)&(D == 1)], probs = c(1-p))
y_p = as.numeric(y_p)
y_1p = as.numeric(y_1p)

# On calcule la lower bound
lb = mean(Y[(S==1)&(D==1)&(Y<=y_1p)]) - mean(Y[(S==1)&(D==0)])
ub = mean(Y[(S==1)&(D==1)&(Y>=y_p)]) - mean(Y[(S==1)&(D==0)])


# Comparaison -------------------------------------------------------------
# On compare avec une simple difference de moyenne
mean(Y[(S==1)&(D==1)]) - mean(Y[(S==1)&(D==0)])


# Calcul de la variance ---------------------------------------------------
# On calcule d'abord a part les mu et les variances incluses dans les formules de variance des bornes
mu_lb = mean(Y[(S==1)&(D==1)&(Y<=y_1p)])
mu_ub = mean(Y[(S==1)&(D==1)&(Y>=y_p)])
v1_temp = sum((Y[(S==1)&(D==1)&(Y<=y_1p)]-mu_lb)^2)/(length(Y[(S==1)&(D==1)&(Y<=y_1p)])-1)
v2_temp = sum((Y[(S==1)&(D==1)&(Y>=y_p)]-mu_lb)^2)/(length(Y[(S==1)&(D==1)&(Y>=y_p)])-1)

# Puis on calcule les variances desirees
v_lb = 1/(mean(S*D)*(1-p))*(v1_temp+(y_1p-mu_lb)^2*p)+(y_1p-mu_lb)^2*((1-mean(S[D==0])-p*(1-mean(D)))/(mean(D)*mean(S[D==0])*(1-mean(D))))
v_ub = 1/(mean(S*D)*(1-p))*(v2_temp+(y_p-mu_ub)^2*p)+(y_p-mu_ub)^2*((1-mean(S[D==0])-p*(1-mean(D)))/(mean(D)*mean(S[D==0])*(1-mean(D))))

# Et enfin les IC
ic_lb = lb -1.96*(sqrt(v_lb/N))
ic_ub = ub +1.96*(sqrt(v_ub/N))
