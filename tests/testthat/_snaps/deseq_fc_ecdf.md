# Expect printed output from running function

    Code
      SeedMatchR::deseq_fc_ecdf(res, list(Background = background.list, mer7m8 = mer7m8.list))$
        stats
    Output
      
      	Exact two-sample Kolmogorov-Smirnov test
      
      data:  target_fc and background_fc
      D^+ = 0.51035, p-value = 0.0001038
      alternative hypothesis: the CDF of x lies above that of y
      

# Expect a change in line color

    Code
      SeedMatchR::deseq_fc_ecdf(res, list(Background = background.list, mer7m8 = mer7m8.list),
      palette = c("red", "green"))$stats
    Output
      
      	Exact two-sample Kolmogorov-Smirnov test
      
      data:  target_fc and background_fc
      D^+ = 0.51035, p-value = 0.0001038
      alternative hypothesis: the CDF of x lies above that of y
      

# Expect a change in plot ordering of lines

    Code
      SeedMatchR::deseq_fc_ecdf(res, list(Background = background.list, mer7m8 = mer7m8.list),
      palette = c("red", "green"), factor_order = c("mer7m8", "Background"))$stats
    Output
      
      	Exact two-sample Kolmogorov-Smirnov test
      
      data:  target_fc and background_fc
      D^+ = 0.51035, p-value = 0.0001038
      alternative hypothesis: the CDF of x lies above that of y
      

# Test parameters stats_test = Wilcoxen, alternative = less.

    Code
      SeedMatchR::ecdf_stat_test(res, background.list, mer7m8.list, stats_test = "Wilcoxen",
        alternative = "less")
    Output
      
      	Wilcoxon rank sum test with continuity correction
      
      data:  target_fc and background_fc
      W = 2150, p-value = 0.9996
      alternative hypothesis: true location shift is less than 0
      

# Check Wass statistic calculation.

    Code
      SeedMatchR::ecdf_stat_test(res, background.list, mer7m8.list, stats_test = "Wass")$
        statistic
    Output
      [1] 0

