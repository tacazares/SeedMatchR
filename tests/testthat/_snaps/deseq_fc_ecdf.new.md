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
      

# Expect a change in plot stats due to target name and background name changes.

    Code
      SeedMatchR::deseq_fc_ecdf(res, list(Background = background.list, mer7m8 = mer7m8.list),
      palette = c("red", "green"), factor_order = c("mer7m8", "Background"),
      null_name = 2, target_name = 1)$stats
    Output
      
      	Exact two-sample Kolmogorov-Smirnov test
      
      data:  target_fc and background_fc
      D^+ = 1.1015e-15, p-value = 1
      alternative hypothesis: the CDF of x lies above that of y
      

# Expect a change in plot stats due to alternative hypothesis change.

    Code
      SeedMatchR::deseq_fc_ecdf(res, list(Background = background.list, mer7m8 = mer7m8.list),
      palette = c("red", "green"), factor_order = c("mer7m8", "Background"),
      null_name = 2, target_name = 1, alternative = "less")$stats
    Output
      
      	Exact two-sample Kolmogorov-Smirnov test
      
      data:  target_fc and background_fc
      D^- = 0.51035, p-value = 0.0001038
      alternative hypothesis: the CDF of x lies below that of y
      

# Expect a change in plot stats due different stats test.

    Code
      SeedMatchR::deseq_fc_ecdf(res, list(Background = background.list, mer7m8 = mer7m8.list),
      palette = c("red", "green"), factor_order = c("mer7m8", "Background"),
      null_name = 2, target_name = 1, stats_test = "Wilcoxen", alternative = "less")$
        stats
    Output
      
      	Wilcoxon rank sum test with continuity correction
      
      data:  target_fc and background_fc
      W = 2150, p-value = 0.9996
      alternative hypothesis: true location shift is less than 0
      

---

    Code
      SeedMatchR::deseq_fc_ecdf(res, list(Background = background.list, mer7m8 = mer7m8.list),
      palette = c("red", "green"), factor_order = c("mer7m8", "Background"),
      null_name = 2, target_name = 1, stats_test = "Wilcoxen", alternative = "greater")$
        stats
    Output
      
      	Wilcoxon rank sum test with continuity correction
      
      data:  target_fc and background_fc
      W = 2150, p-value = 0.00039
      alternative hypothesis: true location shift is greater than 0
      

# Test parameters stats_test = Wilcoxen, alternative = less.

    Code
      SeedMatchR::ecdf_stat_test(res, background.list, mer7m8.list, stats_test = "Wilcoxen",
        alternative = "less")
    Output
      
      	Wilcoxon rank sum test with continuity correction
      
      data:  target_fc and background_fc
      W = 2150, p-value = 0.9996
      alternative hypothesis: true location shift is less than 0
      

# Test parameters stats_test = KS, alternative = less.

    Code
      SeedMatchR::ecdf_stat_test(res, background.list, mer7m8.list, stats_test = "KS",
        alternative = "less")
    Output
      
      	Exact two-sample Kolmogorov-Smirnov test
      
      data:  target_fc and background_fc
      D^- = 0.51035, p-value = 0.0001038
      alternative hypothesis: the CDF of x lies below that of y
      

# Test parameters stats_test = KS, alternative = greater

    Code
      SeedMatchR::ecdf_stat_test(res, background.list, mer7m8.list, stats_test = "KS",
        alternative = "greater")
    Output
      
      	Exact two-sample Kolmogorov-Smirnov test
      
      data:  target_fc and background_fc
      D^+ = 1.1015e-15, p-value = 1
      alternative hypothesis: the CDF of x lies above that of y
      

# Test parameters stats_test = KS, alternative = two.sided

    Code
      SeedMatchR::ecdf_stat_test(res, background.list, mer7m8.list, stats_test = "KS",
        alternative = "less")
    Output
      
      	Exact two-sample Kolmogorov-Smirnov test
      
      data:  target_fc and background_fc
      D^- = 0.51035, p-value = 0.0001038
      alternative hypothesis: the CDF of x lies below that of y
      

# Expect different results based on list order.

    Code
      SeedMatchR::ecdf_stat_test(res, mer7m8.list, background.list, stats_test = "KS",
        alternative = "less")
    Output
      
      	Exact two-sample Kolmogorov-Smirnov test
      
      data:  target_fc and background_fc
      D^- = 1.1015e-15, p-value = 1
      alternative hypothesis: the CDF of x lies below that of y
      

# Check Wass statistic calculation.

    Code
      SeedMatchR::ecdf_stat_test(res, background.list, mer7m8.list, stats_test = "Wass")$
        statistic
    Output
      [1] 0

# Check DTS statistic calculation.

    Code
      SeedMatchR::ecdf_stat_test(res, background.list, mer7m8.list, stats_test = "DTS")$
        statistic
    Condition
      Warning in `scored$statistic <- scored[[1]]`:
      Coercing LHS to a list
    Output
      [1] 74.91157

