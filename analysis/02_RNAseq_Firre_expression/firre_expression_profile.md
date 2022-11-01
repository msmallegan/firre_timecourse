Firre’s expression profile
================
Michael Smallegan

October 31, 2022

## Firre expression at time zero

``` r
# Firre's expression in WT
firre_tpm <- tpm %>%
  filter(gene_name == "Firre") %>%
  pivot_longer(3:ncol(.), names_to = "sample_id", values_to = "tpm") %>%
  left_join(rnaseq_samples, by = "sample_id") %>%
  filter(cell_type == "ESC", timepoint == 0)

firre_mean <- firre_tpm %>%
  group_by(firre_induced, firre_ko) %>%
  summarize(sd = sd(tpm),
            tpm = mean(tpm))

ggplot(firre_tpm, aes(x = firre_induced, y = tpm)) +
  geom_beeswarm(cex = 6, alpha = 0.8) +
  stat_summary(aes(y = tpm, group = firre_induced), fun=mean, colour="red", geom="point") +
  stat_summary(fun='mean', geom='text', 
               aes(label=signif(..y..,3),x=firre_induced)) +
  facet_grid(~firre_ko) +
  ggtitle("Firre expression at 0 min")
```

![](firre_expression_profile_files/figure-gfm/firre-at-zero-1.png)<!-- -->

``` r
ggsave("figures/firre_expression_at_zero.pdf", height = 2, width = 2)

(firre_wt_level <- firre_mean$tpm[firre_mean$firre_induced == "control" & firre_mean$firre_ko == "WT"])
```

    ## [1] 25.13169

## ESC KO (rescue)

### ESC KO long timecourse

``` r
firre_tpm <- tpm %>%
  filter(gene_name == "Firre") %>%
  pivot_longer(3:ncol(.), names_to = "sample_id", values_to = "tpm") %>%
  left_join(rnaseq_samples, by = "sample_id") %>%
  filter(cell_type == "ESC", firre_ko == "KO", timecourse_length == "long")

firre_means <- firre_tpm %>%
  group_by(timepoint, firre_induced) %>%
  summarize(tpm = mean(tpm)) %>%
  mutate(timepoint = as.numeric(as.character(timepoint)))

## FIG 1B
ggplot(firre_tpm, 
       aes(x = as.numeric(as.character(timepoint)), y = tpm, color = firre_induced)) +
  geom_line(data = firre_means, aes(x = timepoint, y = tpm, color = firre_induced,
                                    group = firre_induced)) +
  geom_point() + 
  theme(legend.position = "none") +
  xlab("T (min)") +
  ggtitle("FirreRescue Firre expr. long",
          subtitle = "Fig. 1B")
```

![](firre_expression_profile_files/figure-gfm/esc_ko_long_firre_profile-1.png)<!-- -->

``` r
ggsave("figures/firre_esc_ko_profile_long.pdf", height = 3, width = 3, useDingbats = FALSE)

firre_means <- firre_means %>%
  mutate(rel_wt = tpm / firre_wt_level)
ggplot(firre_means, 
       aes(x = timepoint, y = rel_wt, color = firre_induced)) +
  geom_line(data = firre_means, aes(x = timepoint, y = rel_wt, color = firre_induced,
                                    group = firre_induced)) +
  geom_point() + 
  theme(legend.position = "none") +
  xlab("T (min)") +
  ggtitle("FirreRescue Firre expr. rel. WT Firre")
```

![](firre_expression_profile_files/figure-gfm/esc_ko_long_firre_profile-2.png)<!-- -->

``` r
ggsave("figures/firre_esc_ko_profile_long_relwt.pdf", height = 2, width = 2, useDingbats = FALSE)
```

### ESC KO short timecourse

``` r
firre_tpm <- tpm %>%
  filter(gene_name == "Firre") %>%
  pivot_longer(3:ncol(.), names_to = "sample_id", values_to = "tpm") %>%
  left_join(rnaseq_samples, by = "sample_id") %>%
  filter(cell_type == "ESC", firre_ko == "KO", timecourse_length == "short")

firre_means <- firre_tpm %>%
  group_by(timepoint, firre_induced) %>%
  summarize(tpm = mean(tpm)) %>%
  mutate(timepoint = as.numeric(as.character(timepoint)))

## Fig. 2A
ggplot(firre_tpm, 
       aes(x = as.numeric(as.character(timepoint)), y = tpm, color = firre_induced)) +
  geom_line(data = firre_means, aes(x = timepoint, y = tpm, color = firre_induced,
                                    group = firre_induced)) +
  geom_point() + 
  theme(legend.position = "none") +
  xlab("T (min)") +
  ggtitle("FirreRescue short Firre expr.",
          subtitle = "Fig. 2A")
```

![](firre_expression_profile_files/figure-gfm/esc_ko_short_firre_profile-1.png)<!-- -->

``` r
ggsave("figures/firre_esc_ko_profile_short.pdf", height = 3, width = 3, useDingbats = FALSE)

firre_means <- firre_means %>%
  mutate(rel_wt = tpm / firre_wt_level)

## Fig. S4A
ggplot(firre_means, 
       aes(x = timepoint, y = rel_wt, color = firre_induced)) +
  geom_line(data = firre_means, aes(x = timepoint, y = rel_wt, color = firre_induced,
                                    group = firre_induced)) +
  geom_point() + 
  theme(legend.position = "none") +
  geom_text(data = firre_means %>% filter(rel_wt == max(rel_wt)), aes(y = rel_wt + 2, label = round(rel_wt,1))) +
  scale_x_continuous(breaks = seq(0,330, by = 30), labels = c("0", "", "1", "", "2", "", "3", "", "4", "", "5", "")) +
  xlab("Time (h)") +
  ggtitle("FirreRescue Firre rel. WT",
          subtitle = "Fig. S4A")
```

![](firre_expression_profile_files/figure-gfm/esc_ko_short_firre_profile-2.png)<!-- -->

``` r
ggsave("figures/firre_esc_ko_profile_short_relwt.pdf", height = 2, width = 2, useDingbats = FALSE)
```

### ESC KO Combined

``` r
firre_tpm <- tpm %>%
  filter(gene_name == "Firre") %>%
  pivot_longer(3:ncol(.), names_to = "sample_id", values_to = "tpm") %>%
  left_join(rnaseq_samples, by = "sample_id") %>%
  filter(cell_type == "ESC", firre_ko == "KO")

firre_means <- firre_tpm %>%
  group_by(timepoint, firre_induced) %>%
  summarize(tpm = mean(tpm)) %>%
  mutate(timepoint = as.numeric(as.character(timepoint)))

ggplot(firre_tpm, 
       aes(x = as.numeric(as.character(timepoint)), y = tpm, color = firre_induced)) +
  geom_line(data = firre_means, aes(x = timepoint, y = tpm, color = firre_induced,
                                    group = firre_induced)) +
  geom_point() + 
  theme(legend.position = "none") +
  xlab("T (min)") +
  ggtitle("FirreRescue full Firre expr.")
```

![](firre_expression_profile_files/figure-gfm/esc_ko_combined_firre_profile-1.png)<!-- -->

``` r
ggsave("figures/firre_esc_ko_profile.pdf", height = 3, width = 3, useDingbats = FALSE)
```

## ESC WT (overexpression)

### ESC WT long timecourse

``` r
firre_tpm <- tpm %>%
  filter(gene_name == "Firre") %>%
  pivot_longer(3:ncol(.), names_to = "sample_id", values_to = "tpm") %>%
  left_join(rnaseq_samples, by = "sample_id") %>%
  filter(cell_type == "ESC", firre_ko == "WT", timecourse_length == "long")

firre_means <- firre_tpm %>%
  group_by(timepoint, firre_induced) %>%
  summarize(tpm = mean(tpm)) %>%
  mutate(timepoint = as.numeric(as.character(timepoint)))

## Fig. S2C
ggplot(firre_tpm, 
       aes(x = as.numeric(as.character(timepoint)), y = tpm, color = firre_induced)) +
  geom_line(data = firre_means, aes(x = timepoint, y = tpm, color = firre_induced,
                                    group = firre_induced)) +
  geom_point() + 
  theme(legend.position = "none") +
  xlab("T (min)") +
  ggtitle("FirreOE Firre expr.",
          subtitle = "Fig. S2C")
```

![](firre_expression_profile_files/figure-gfm/esc_wt_long_firre_profile-1.png)<!-- -->

``` r
ggsave("figures/firre_esc_wt_profile_long.pdf", height = 3, width = 3, useDingbats = FALSE)

firre_means <- firre_means %>%
  mutate(rel_wt = tpm / firre_wt_level)

ggplot(firre_means, 
       aes(x = timepoint, y = rel_wt, color = firre_induced)) +
  geom_line(data = firre_means, aes(x = timepoint, y = rel_wt, color = firre_induced,
                                    group = firre_induced)) +
  geom_point() + 
  theme(legend.position = "none") +
  geom_text(data = firre_means %>% filter(rel_wt == max(rel_wt)), aes(y = rel_wt + 0.5, label = round(rel_wt,1))) +
  scale_x_continuous(breaks = c(0, 720, 1440, 2880, 5670), labels = c("0", "12", "24", "48", "96")) +
  xlab("Time (h)") +
  ylim(0,6) +
  ggtitle("FirreOE Firre rel. WT")
```

![](firre_expression_profile_files/figure-gfm/esc_wt_long_firre_profile-2.png)<!-- -->

``` r
ggsave("figures/firre_esc_wt_profile_long_relwt.pdf", height = 2, width = 2, useDingbats = FALSE)
```

### ESC WT short timecourse

``` r
firre_tpm <- tpm %>%
  filter(gene_name == "Firre") %>%
  pivot_longer(3:ncol(.), names_to = "sample_id", values_to = "tpm") %>%
  left_join(rnaseq_samples, by = "sample_id") %>%
  filter(cell_type == "ESC", firre_ko == "WT", timecourse_length == "short")

firre_means <- firre_tpm %>%
  group_by(timepoint, firre_induced) %>%
  summarize(tpm = mean(tpm)) %>%
  mutate(timepoint = as.numeric(as.character(timepoint)))


ggplot(firre_tpm, 
       aes(x = as.numeric(as.character(timepoint)), y = tpm, color = firre_induced)) +
  geom_line(data = firre_means, aes(x = timepoint, y = tpm, color = firre_induced,
                                    group = firre_induced)) +
  geom_point() + 
  theme(legend.position = "none") +
  xlab("T (min)") +
  ggtitle("FirreOE short Firre expr.")
```

![](firre_expression_profile_files/figure-gfm/esc_wt_short_firre_profile-1.png)<!-- -->

``` r
ggsave("figures/firre_esc_wt_profile_short.pdf", height = 3, width = 3, useDingbats = FALSE)

firre_means <- firre_means %>%
  mutate(rel_wt = tpm / firre_wt_level)

## Fig. S4B
ggplot(firre_means, 
       aes(x = timepoint, y = rel_wt, color = firre_induced)) +
  geom_line(data = firre_means, aes(x = timepoint, y = rel_wt, color = firre_induced,
                                    group = firre_induced)) +
  geom_point() + 
  theme(legend.position = "none") +
  geom_text(data = firre_means %>% filter(rel_wt == max(rel_wt)), aes(y = rel_wt + 2, label = round(rel_wt,1))) +
  scale_x_continuous(breaks = seq(0,330, by = 30), labels = c("0", "", "1", "", "2", "", "3", "", "4", "", "5", "")) +
  xlab("Time (h)") +
  ylim(0,50) +
  ggtitle("FirreOE short Firre rel. WT",
          subtitle = "Fig. S4B")
```

![](firre_expression_profile_files/figure-gfm/esc_wt_short_firre_profile-2.png)<!-- -->

``` r
ggsave("figures/firre_esc_wt_profile_short_relwt.pdf", height = 2, width = 2, useDingbats = FALSE)
```

### ESC WT Combined

``` r
firre_tpm <- tpm %>%
  filter(gene_name == "Firre") %>%
  pivot_longer(3:ncol(.), names_to = "sample_id", values_to = "tpm") %>%
  left_join(rnaseq_samples, by = "sample_id") %>%
  filter(cell_type == "ESC", firre_ko == "WT")

firre_means <- firre_tpm %>%
  group_by(timepoint, firre_induced) %>%
  summarize(tpm = mean(tpm)) %>%
  mutate(timepoint = as.numeric(as.character(timepoint)))

ggplot(firre_tpm, 
       aes(x = as.numeric(as.character(timepoint)), y = tpm, color = firre_induced)) +
  geom_line(data = firre_means, aes(x = timepoint, y = tpm, color = firre_induced,
                                    group = firre_induced)) +
  geom_point() + 
  theme(legend.position = "none") +
  xlab("T (min)") +
  ggtitle("FirreOE full Firre expr.")
```

![](firre_expression_profile_files/figure-gfm/esc_wt_combined_firre_profile-1.png)<!-- -->

``` r
ggsave("figures/firre_esc_wt_profile.pdf", height = 3, width = 3, useDingbats = FALSE)
```

## Example Firre responder gene profiles

``` r
rapgef4_tpm <- tpm %>%
  filter(gene_name == "Rapgef4") %>%
  pivot_longer(3:ncol(.), names_to = "sample_id", values_to = "tpm") %>%
  left_join(rnaseq_samples, by = "sample_id") %>%
  filter(cell_type == "ESC", timecourse_length == "short")

rapgef4_means <- rapgef4_tpm %>%
  group_by(timepoint, firre_induced) %>%
  summarize(tpm = mean(tpm)) %>%
  mutate(timepoint = as.numeric(as.character(timepoint)))

## Fig. 2E
ggplot(rapgef4_tpm, 
       aes(x = as.numeric(as.character(timepoint)), y = tpm, color = firre_induced)) +
  geom_line(data = rapgef4_means, aes(x = timepoint, y = tpm, color = firre_induced,
                                      group = firre_induced)) +
  geom_point() + 
  facet_wrap(~firre_ko) +
  theme(legend.position = "none") +
  xlab("T (min)") +
  ggtitle("Rapgef4 profile short",
          subtitle = "Fig. 2E")
```

![](firre_expression_profile_files/figure-gfm/rapgef4-profile-1.png)<!-- -->

``` r
ggsave("figures/rapgef4_short.pdf", height = 1.5, width = 3, useDingbats = FALSE)
```

``` r
adgrg1_tpm <- tpm %>%
  filter(gene_name == "Adgrg1") %>%
  pivot_longer(3:ncol(.), names_to = "sample_id", values_to = "tpm") %>%
  left_join(rnaseq_samples, by = "sample_id") %>%
  filter(cell_type == "ESC", timecourse_length == "short")

adgrg1_means <- adgrg1_tpm %>%
  group_by(timepoint, firre_induced, firre_ko) %>%
  summarize(tpm = mean(tpm)) %>%
  mutate(timepoint = as.numeric(as.character(timepoint)))

## Fig. 2E
ggplot(adgrg1_tpm, 
       aes(x = as.numeric(as.character(timepoint)), y = tpm, color = firre_induced)) +
  geom_line(data = adgrg1_means, aes(x = timepoint, y = tpm, color = firre_induced,
                                     group = firre_induced)) +
  geom_point() + 
  facet_grid(~firre_ko) +
  theme(legend.position = "none") +
  xlab("T (min)") +
  ggtitle("Adgrg1 profile short",
          subtitle = "Fig. 2E")
```

![](firre_expression_profile_files/figure-gfm/adgrg1-profile-1.png)<!-- -->

``` r
ggsave("figures/adgrg1_short.pdf", height = 1.5, width = 3, useDingbats = FALSE)
```

``` r
shf_tpm <- tpm %>%
  filter(gene_name == "Shf") %>%
  pivot_longer(3:ncol(.), names_to = "sample_id", values_to = "tpm") %>%
  left_join(rnaseq_samples, by = "sample_id") %>%
  filter(cell_type == "ESC", timecourse_length == "short")

shf_means <- shf_tpm %>%
  group_by(timepoint, firre_induced, firre_ko) %>%
  summarize(tpm = mean(tpm)) %>%
  mutate(timepoint = as.numeric(as.character(timepoint)))

## Fig. 2E
ggplot(shf_tpm, 
       aes(x = as.numeric(as.character(timepoint)), y = tpm, color = firre_induced)) +
  geom_line(data = shf_means, aes(x = timepoint, y = tpm, color = firre_induced,
                                  group = firre_induced)) +
  geom_point() + 
  facet_grid(~firre_ko) +
  theme(legend.position = "none") +
  xlab("T (min)") +
  ggtitle("Shf profile short",
          subtitle = "Fig. 2E")
```

![](firre_expression_profile_files/figure-gfm/shf-profile-1.png)<!-- -->

``` r
ggsave("figures/shf_short.pdf", height = 1.5, width = 3, useDingbats = FALSE)
```

### Interpolation of Firre Rescue expression at 15 minutes

``` r
firre_0_30 <- tpm %>%
  filter(gene_name == "Firre") %>%
  pivot_longer(3:ncol(.), names_to = "sample_id", values_to = "tpm") %>%
  left_join(rnaseq_samples, by = "sample_id") %>%
  filter(cell_type == "ESC", timecourse_length == "short", timepoint %in% c(0,30, 60),
         firre_induced == "firre_induced", experiment == "ESC_KO") %>%
  group_by(gene_name, timepoint) %>%
  summarize(tpm = mean(tpm)) %>%
  mutate(timepoint = as.numeric(as.character(timepoint)))


firre_expr_at_15 <- (firre_0_30$tpm[which(firre_0_30$timepoint == 30)] - firre_0_30$tpm[which(firre_0_30$timepoint == 0)]) / 2

firre_wt_range <- tpm %>%
  filter(gene_name == "Firre") %>%
  pivot_longer(3:ncol(.), names_to = "sample_id", values_to = "tpm") %>%
  left_join(rnaseq_samples, by = "sample_id") %>%
  filter(cell_type == "ESC", timepoint %in% c(0), experiment == "ESC_WT",
         firre_induced == "control")


a <- mean(firre_wt_range$tpm)
s <-  sd(firre_wt_range$tpm)
n <- 5
error <- qnorm(0.975)*s/sqrt(n)
left <- a-error
right <- a+error
crossing_df <- tibble(gene_name = "Firre", x1 = 15, x2 = 15, y1 = 0, y2 = firre_expr_at_15)

firre_expr_at_15 / a
```

    ## [1] 2.675134

``` r
ggplot(firre_0_30 %>% filter(timepoint %in% c(0,30)), aes(x = timepoint, y = tpm)) +
  geom_hline(yintercept = a, lty = 2, color = "red") +
  geom_hline(yintercept = left, alpha = 0.3) +
  geom_hline(yintercept = right, alpha = 0.3) +
  geom_segment(mapping = aes(x = x1, xend = x2, y = y1, yend = y2), data = crossing_df, lty = 2) +
  geom_point() +
  geom_line(lty = 2, alpha = .3, aes(group = gene_name)) +
  annotate(geom = "text", x = 18, y = 70, label = paste0("y=",round(firre_expr_at_15)), size = 6) +
  ggtitle("Firre expr. in Firre RESCUE") +
  scale_x_continuous(breaks = seq(0,30,5)) +
  theme(plot.title = element_text(size = 13, face = "bold")) +
  xlab("Time (min)") +
  ylab("TPM")
```

![](firre_expression_profile_files/figure-gfm/firre-rescue-at-15-1.png)<!-- -->

``` r
ggsave("figures/firre_rescue_15min_crossing_point.pdf", height = 2.6, width = 3.5)
```