library(mvtnorm)
#library(Hmisc)

# Definition des parametres -----------------------------------------------
sigma_mat = matrix(c(0.25, 0.5, 0.5, 1), nrow=2, ncol=2, byrow = TRUE) # Matrice de var-cov des erreurs
beta = 0.2 # Effet traitement sur salaire
gamma = 0.3 # Effet traitement sur employabilite 
pi = matrix (c(9.8,0.35,0.3,1), nrow = 2, ncol = 2) # Constantes pour les 2 eq
p = 0.5  # Proportion de l'echantillon traitee

# Simulations -------------------------------------------------------------
# On veut constituer un echantillon de taille N
N = 100000
set.seed(123)

# On commence par simuler les erreurs
UV = rmvnorm(n = N, mean = c(0,0), sigma = sigma_mat)

# On simule X: constante et propensity score
X = matrix(data = 1, nrow = N, ncol = 2)
# Pour le score, on veut une assignation uniforme dans 5 groupes
#runifdisc<-function(n, min, max) sample(min:max, n, replace=T)
#X[,2] = runifdisc(n = N, min = -2, max = 2)
X[,2] = rnorm(n = N, mean = 0, sd = 1)

# On assigne les indiv a un groupe de traitement (proportion p) et de controle (1-p)
D = rbinom(n = N, size = 1, prob = p)

# On introduit les variables latentes et l'observee
Y_hat = D*beta + X%*%pi[,1] + UV[,1]
Z_hat = D*gamma + X%*%pi[,2] + UV[,2]
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

# Estimation sans covariates ---------------------------------------------------------
# On estime d'abord la proportion de trimming
p = (sum((S==1)&(D == 1))/length(D ==1) - sum((S==1)&(D == 0))/length(D ==0))/(sum((S==1)&(D == 1))/length(D == 1))

# On estime les quantiles associes
y_p = as.numeric(quantile(x = Y[(S==1)&(D == 1)], probs = c(p)))
y_1p = as.numeric(quantile(x = Y[(S==1)&(D == 1)], probs = c(1-p)))

# On calcule la lower bound
lb = mean(Y[(S==1)&(D==1)&(Y<=y_1p)]) - mean(Y[(S==1)&(D==0)])
ub = mean(Y[(S==1)&(D==1)&(Y>=y_p)]) - mean(Y[(S==1)&(D==0)])

# Calcul de la variance sans covariates -----------------------------------------------
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

# On supprime les variables devenues inutiles
rm(y_p, y_1p, v1_temp, v2_temp, mu_lb, mu_ub)

# Utilisation des covariates pour affiner les bornes ----------------------
# On definit 5 groupes, on commence par regresser le salaire sur X sans constante (deja dans X)
# pour ceux dont on observe le salaire
data = data.frame(Y[S==1],X[S==1,])
colnames(data) = c("Y","X1", "X2")
fit <- lm(Y~X1 + X2 -1, data = data)
summary(fit) 

# On applique ensuite les coefficients pour predire un salaire a tout l'echantillon
Y_predicted = X %*% fit$coefficients

# On calcule les quintiles empiriques pour classer dans 5 groupes
quantiles = as.numeric(quantile(x = Y_predicted, probs = c(1/5,2/5,3/5,4/5)))
group = vector(mode = "numeric", length = N)
group = (Y_predicted < quantiles[1])*1
for (i in 1:3) {
  group = group + (i+1)*((Y_predicted>=quantiles[i])&(Y_predicted<quantiles[i+1]))
}
group = group + (Y_predicted > quantiles[4])*5

# On definit 5 vecteurs pour les variables d'interet
p_cov = vector(mode = "numeric", length = 5)
y_p_cov = vector(mode = "numeric", length = 5)
y_1p_cov = vector(mode = "numeric", length = 5)
lb_cov = vector(mode = "numeric", length = 5)
ub_cov = vector(mode = "numeric", length = 5)

for (i in 1:5){
  # On estime d'abord la proportion de trimming
  p_cov[i] = (sum((S==1)&(D==1)&(group==i))/length((D ==1)&(group==i)) - sum((S==1)&(D==0)&(group==i))/length((D ==0)&(group==i)))/(sum((S==1)&(D==1)&(group==i))/length((D==1)&(group==i)))
  
  # On estime les quantiles associes
  y_p_cov[i] = as.numeric(quantile(x = Y[(S==1)&(D==1)&(group==i)], probs = c(p_cov[i])))
  y_1p_cov[i] = as.numeric(quantile(x = Y[(S==1)&(D==1)&(group==i)], probs = c(1-p_cov[i])))
  
  # On calcule la lower bound
  lb_cov[i] = mean(Y[(S==1)&(D==1)&(group==i)&(Y<=y_1p_cov[i])]) - mean(Y[(S==1)&(D==0)&(group==i)])
  ub_cov[i] = mean(Y[(S==1)&(D==1)&(group==i)&(Y>=y_p_cov[i])]) - mean(Y[(S==1)&(D==0)&(group==i)])
}

# Pour calculer la moyenne des ub et lw, on pondere par 1- trimming proportion
weight = (1-p_cov)/sum((1-p_cov))
ub_final = as.numeric(ub_cov %*% weight)
lb_final = as.numeric(lb_cov %*% weight)

# On supprime les variables inutiles
rm (quantiles, i, Y_predicted)


# Calcul des variances dans chaque groupe ------------------------------------

# On definit 2 vecteurs pour les variables d'interet
v_lb_cov = vector(mode = "numeric", length = 5)
v_ub_cov = vector(mode = "numeric", length = 5)

for (i in 1:5){
  y_1p = y_1p_cov[i]
  y_p = y_p_cov[i]
  p = p_cov[i]
  # On calcule d'abord a part les mu et les variances incluses dans les formules de variance des bornes
  mu_lb = mean(Y[(S==1)&(D==1)&(Y<=y_1p)&(group==i)])
  mu_ub = mean(Y[(S==1)&(D==1)&(Y>=y_p)&(group==i)])
  v1_temp = sum((Y[(S==1)&(D==1)&(Y<=y_1p)&(group==i)]-mu_lb)^2)/(length(Y[(S==1)&(D==1)&(Y<=y_1p)&(group==i)])-1)
  v2_temp = sum((Y[(S==1)&(D==1)&(Y>=y_p)&(group==i)]-mu_lb)^2)/(length(Y[(S==1)&(D==1)&(Y>=y_p)&(group==i)])-1)
  
  # Puis on calcule les variances desirees
  v_lb_cov[i] = 1/(mean(S[group==i]*D[group==i])*(1-p))*(v1_temp+(y_1p-mu_lb)^2*p)+(y_1p-mu_lb)^2*((1-mean(S[(D==0)&(group==1)])-p*(1-mean(D[group==i])))/(mean(D[group==i])*mean(S[(D==0)&(group==i)])*(1-mean(D[group==i]))))
  v_ub_cov[i] = 1/(mean(S[group==i]*D[group==i])*(1-p))*(v2_temp+(y_p-mu_ub)^2*p)+(y_p-mu_ub)^2*((1-mean(S[(D==0)&(group==1)])-p*(1-mean(D[group==i])))/(mean(D[group==i])*mean(S[(D==0)&(group==i)])*(1-mean(D[group==i]))))
}


# Calcul de la variance totale et IC --------------------------------------------
# Variances
v_lb_final = as.numeric(v_lb_cov %*% weight + ((lb_cov - lb_final)^2)%*% weight)
v_ub_final = as.numeric(v_ub_cov %*% weight + ((ub_cov - ub_final)^2)%*% weight)

# Et enfin les IC
ic_lb_final = lb_final -1.96*(sqrt(v_lb_final/N))
ic_ub_final = ub_final +1.96*(sqrt(v_ub_final/N))

# On supprime les variables inutiles
rm (weight,v1_temp, v2_temp, mu_lb, mu_ub)


# Presentation des resultats ----------------------------------------------
# Vraie valeur
beta

# Simple difference de moyenne
mean(Y[(S==1)&(D==1)]) - mean(Y[(S==1)&(D==0)])

# Methode sans covariate
lb # lower bound
ub # upper bound
ic_lb # lower bound de l'ic de la lower bound
ic_ub # upper bound de l'ic de la upper bound

# Methode avec covariate
lb_final # lower bound
ub_final # upper bound
ic_lb_final # lower bound de l'ic de la lower bound
ic_ub_final # upper bound de l'ic de la upper bound
