
d_wide <- read.delim("/Users/Dasha/work/Sardinia/W4H/olink/data/olink_clean_CVD+INF.txt", as.is = T, check.names = F, sep = "\t")
ID <- gsub("_.*", "", d_wide$SampleID)
TP <- gsub(".*_", "", d_wide$SampleID)
d_wide <- cbind(ID, TP, d_wide)

prot1 <- 'FGF5'
prot2 <- 'CA4'
prot3 <- 'MMP10'

prot1 <- 'DBNL'
prot2 <- 'CRKL'


d_subs <- d_wide[,c("ID", "TP", prot1, prot2, prot3)]
colnames(d_subs) <- c("ID", "TP", "prot1", "prot2", 'prot3')
d_subs$TPcat <- as.factor(d_subs$TP)
d_subs <- na.omit(d_subs)


ggplot(d_subs, aes(y = prot2, x = TP, group = ID)) + 
  geom_line(aes(color = ID), alpha = 0.5) + 
  geom_point() + 
  theme_bw() +
  theme(legend.position="none") 

d_subs2 <- d_subs[d_subs$ID %in% sample(unique(d_subs$ID), 20),]
ggplot(d_subs2, aes(y = prot2, x = TP, group = ID)) + 
  geom_line(aes(color = ID)) + 
  geom_point() + 
  theme_bw() +
  theme(legend.position="none") 


# 1. Do protein levels change with time

# lmm
m <- lmer(prot3 ~ poly(TP,3) + (1|ID), data = d_subs)
m0 <- lmer(prot3 ~ (1|ID), data = d_subs)
anova(m, m0)

ggplot(d_subs, aes(x = TP, y = prot3, color = ID)) + geom_point() + geom_line()+ theme_bw() +
  theme(legend.position="none") 
ggplot(d_subs, aes(x = TP, y = prot3, color = ID)) + geom_point() + geom_line(aes(y = predict(m))) + theme_bw() +
  theme(legend.position="none") 

# gls
gls_fit <- gls(prot3 ~ poly(TP, 3), data = d_subs, correlation = corCompSymm(form = ~ TP | ID))
anova(gls_fit)

ggplot(d_subs, aes(x = TP, y = prot3)) + geom_jitter(width = 0.2, alpha = 0.4) + geom_line(aes(y = predict(gls_fit))) + theme_bw()


# 2. Is prot1 assoicated with prot2?
m <- lmer(prot1 ~  prot2 + TP + (1|ID), data = d_subs)
m0 <- lmer(prot1 ~  TP + (1|ID), data = d_subs)

anova(m,m0)


ggplot(d_subs) + 
  geom_jitter(aes(x = TP, y = prot1), color = 'blue', alpha = 0.2, width = 0.2) +
  geom_line(aes(x = TP, y = prot1), color = 'blue', stat="smooth",method = "lm", formula =y ~ poly(x, 3, raw=TRUE), se = F) +
  geom_jitter(aes(x = TP, y = prot2), color = 'red', alpha = 0.2, width = 0.2) +
  geom_line(aes(x = TP, y = prot2), color = 'red', stat="smooth",method = "lm", formula =y ~ poly(x, 3, raw=TRUE), se = F) +
  theme_bw()

plot_together(d_wide, d_wide, prot1, prot2)

scatter_col_tp(d_wide, d_wide , prot1, prot2)

m_inter <- lmer(prot1 ~  prot2 + TP + prot2*TP + (1|ID), data = d_subs)
anova(m_inter, m)

m2 <- lmer(prot1 ~  prot2 + TP + (TP|ID), data = d_subs)
m02 <- lmer(prot1 ~  TP + (TP|ID), data = d_subs)
anova(m2,m02)

m3 <- lmer(prot1 ~  prot2 + TP + (1|ID) + (1|TP), data = d_subs)
m03 <- lmer(prot1 ~   TP + (1|ID) + (1|TP), data = d_subs)
anova(m3,m03)




#
#
#

gls_fit1 <- gls(prot1 ~ poly(TP, 1), data = d_subs, correlation = corCompSymm(form = ~ TP | ID))
gls_fit2 <- gls(prot1 ~ poly(TP, 2), data = d_subs, correlation = corCompSymm(form = ~ TP | ID))
gls_fit3 <- gls(prot1 ~ poly(TP, 3), data = d_subs, correlation = corCompSymm(form = ~ TP | ID))

gls_fit1$coefficients
gls_fit2$coefficients
gls_fit3$coefficients


prot = 'CHRDL2'
get_simplified_coefs <- function(d_wide, prot, scale = T){

  d_subs <- d_wide[,c("ID", "TP", prot)]
  colnames(d_subs) <- c("ID", "TP", "prot")
  d_subs$TP <- as.numeric(d_subs$TP)
  d_subs <- na.omit(d_subs)

  if (scale) d_subs$prot <- scale(d_subs$prot)

  gls_fit <- gls(prot ~ poly(TP, 3), data = d_subs, correlation = corSymm(form = ~ TP | ID))
  coefs <- gls_fit$coefficients
  
  coefs <- coefs[2,3,4]
  coefs[coefs > 2] <- 3
  coefs[coefs < -2] <- -3
  
  coefs[coefs > 1 & coefs < 2] <- 2
  coefs[coefs < -1 & coefs > -2] <- -2
  
  coefs[coefs < 1 & coefs > 0 ] <- 1
  coefs[coefs > -1 & coefs < 0 ] <- -1
  
  return(coefs)
} 
