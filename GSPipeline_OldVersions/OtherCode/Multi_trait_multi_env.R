#download package
library(BMTME)

#download data
data("MaizeToy")

#rename the rows of the phenotypic dataset
phenoMaizeToy<-(phenoMaizeToy[order(phenoMaizeToy$Env,
                                    phenoMaizeToy$Line),])
rownames(phenoMaizeToy)=1:nrow(phenoMaizeToy)
head(phenoMaizeToy)

#Then the design matrices for the line effects, the environment and the genotype·environment interaction are generated
LG <- cholesky(genoMaizeToy)
ZG <- model.matrix(~0 + as.factor(phenoMaizeToy$Line))
Z.G <- ZG %*% LG
Z.E <- model.matrix(~0 + as.factor(phenoMaizeToy$Env))
ZEG <- model.matrix(~0 + as.factor(phenoMaizeToy$Line):as.factor(phenoMaizeToy$Env))
G2 <- kronecker(diag(length(unique(phenoMaizeToy$Env))),
                data.matrix(genoMaizeToy))
LG2 <- cholesky(G2)
Z.EG <- ZEG %*% LG2
Y <- as.matrix(phenoMaizeToy[, -c(1, 2)])

#fit the model
fm <- BMTME(Y = Y, X = Z.E, Z1 = Z.G, Z2 = Z.EG,
            nIter =15000, burnIn =10000, thin = 2,bs = 50)

names(fm)

#To extract thematrix of covariances between traits
COV_TraitGenetic <- fm$varTrait
COV_TraitGenetic
COR_TraitGenetic <- cov2cor(COV_TraitGenetic)
COR_TraitGenetic
COV_ResGenetic <- fm$vare
COV_ResGenetic

#The observed and predicted values of each trait can be plotted using the following plot(), function 
plot(fm, trait='Yield')

# fivefold CV strategy and its implementation with the BMTME function:
pheno <- data.frame(GID = phenoMaizeToy[, 1], Env =
                      phenoMaizeToy[, 2], Response = phenoMaizeToy[, 3])
CrossV <- CV.KFold(pheno, DataSetID = 'GID', K = 5,
                   set_seed = 123)
pm <- BMTME(Y = Y, X = Z.E, Z1 = Z.G, Z2 = Z.EG,
            nIter = 10000, burnIn = 2000, thin = 2,bs = 50, testingSet = CrossV)
summary(pm)

boxplot(pm, select='Pearson', las = 2)