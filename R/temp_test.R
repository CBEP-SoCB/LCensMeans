
library(ggplot2)

df = data.frame(sim = sort(rlnorm(25,2,3)),
                cens=c(rep.int(TRUE,5), rep.int(FALSE,20)))
df$sim[1:4] <- df$sim[5]
est <- simulate_mean_censored(lmu = 2, lsigma = 3, cutoff = df$sim)


ggplot(df, aes(x = 1:25)) +
  geom_line(aes(y = sim, color=cens)) +
  geom_point(aes(y = est)) +
  scale_y_log10() +
  scale_color_discrete(name = 'Censored') +
  theme_minimal() +
  xlab('Rank Order') +
  ylab('Raw Data (Line) and Conditional Means (points)')

vals <- sub_conditional_means(cc = df$sim, flg =df$cens)
ggplot(df, aes(x = 1:25)) +
  geom_line(aes(y = sim, color=cens)) +
  geom_point(aes(y = vals)) +
  scale_y_log10() +
  scale_color_discrete(name = 'Censored') +
  theme_minimal() +
  xlab('Rank Order') +
  ylab('Raw Data (Line) and Data with Substitutions (points)')


