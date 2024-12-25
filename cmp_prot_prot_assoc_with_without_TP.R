prot1 = 'CA4'
prot2 <- 'FGF5'

prots1 <- sample(colnames(d_wide)[4:ncol(d_wide)], 100)
prots2 <- sample(colnames(d_wide)[4:ncol(d_wide)], 100)

res <- data.frame(matrix(nrow = length(prots1) * length(prots2), ncol = 7))
cnt <- 1
for (prot1 in prots1) {
  for (prot2 in prots2){
    
    d_subs <- inner_join(d_wide[,c("ID", "TP", prot1, prot2)], covariates, by = c('ID'))
    
    colnames(d_subs)[1:4] <- c("ID", "TP", "prot1", "prot2")
    d_subs$TP <- as.numeric(d_subs$TP)
    d_subs <- na.omit(d_subs)
    
    model1 <- lmer(prot1 ~ prot2 + (1|ID), data = d_subs)
    model10 <- lmer(prot1 ~  (1|ID), data = d_subs)
    est1 <- summary(model1)$coefficients["prot2", "Estimate"]
    an1 <- suppressMessages(anova(model1, model10))
    pval1 <- an1$`Pr(>Chisq)`[2]
       
    model2 <- lmer(prot1 ~ prot2 + TP + (1|ID), data = d_subs)
    model20 <- lmer(prot1 ~ TP +(1|ID), data = d_subs)

    est2 <- summary(model2)$coefficients["prot2", "Estimate"]
    an2 <- suppressMessages(anova(model2, model20))
    pval2 <- an2$`Pr(>Chisq)`[2]
    an22 <- suppressMessages(anova(model2, model1))
    pval22 <- an22$`Pr(>Chisq)`[2]
    res[cnt, ] <- c(prot1, prot2, pval1, est1, pval2, est2, pval22)
    cnt <- cnt + 1

  }
}
colnames(res) <- c("prot1", "prot2", "pval1", "est1", "pval2", "est2", 'pval_tp')
res <- na.omit(res) %>%
  mutate(across(-c(prot1, prot2), as.numeric)) 
res$logp1 <- -log10(res$pval1)
res$logp2 <- -log10(res$pval2)

plot(res$est1, res$est2, pch = 16, xlab = 'no visit in predictors', ylab = 'with visit in predictors')
plot(res$logp1, res$logp2, pch = 16, xlab = 'no visit in predictors', ylab = 'with visit in predictors')
