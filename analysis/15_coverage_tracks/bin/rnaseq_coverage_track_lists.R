rnaseq_dir <- "/scratch/Shares/rinn/Michael/firre_timecourse/rnaseq1/results/star_salmon/bigwig/"

esc_ko_short_all <- list(`0h R1` = list(pos = file.path(rnaseq_dir,
"ESC_KO_firre_induced_0min_dox_short_rep1.forward.bigWig"), neg =
file.path(rnaseq_dir,
"ESC_KO_firre_induced_0min_dox_short_rep1.reverse.bigWig")), `0h R2` = list(pos
= file.path(rnaseq_dir,
"ESC_KO_firre_induced_0min_dox_short_rep2.forward.bigWig"), neg =
file.path(rnaseq_dir,
"ESC_KO_firre_induced_0min_dox_short_rep2.reverse.bigWig")), `0.5h R1` =
list(pos = file.path(rnaseq_dir,
"ESC_KO_firre_induced_30min_dox_short_rep1.forward.bigWig"), neg =
file.path(rnaseq_dir,
"ESC_KO_firre_induced_30min_dox_short_rep1.reverse.bigWig")), `0.5h R2` =
list(pos = file.path(rnaseq_dir,
"ESC_KO_firre_induced_30min_dox_short_rep2.forward.bigWig"), neg =
file.path(rnaseq_dir,
"ESC_KO_firre_induced_30min_dox_short_rep2.reverse.bigWig")), `1h R1` = list(pos
= file.path(rnaseq_dir,
"ESC_KO_firre_induced_60min_dox_short_rep1.forward.bigWig"), neg =
file.path(rnaseq_dir,
"ESC_KO_firre_induced_60min_dox_short_rep1.reverse.bigWig")), `1h R2` = list(pos
= file.path(rnaseq_dir,
"ESC_KO_firre_induced_60min_dox_short_rep2.forward.bigWig"), neg =
file.path(rnaseq_dir,
"ESC_KO_firre_induced_60min_dox_short_rep2.reverse.bigWig")), `1.5h R1` =
list(pos = file.path(rnaseq_dir,
"ESC_KO_firre_induced_90min_dox_short_rep1.forward.bigWig"), neg =
file.path(rnaseq_dir,
"ESC_KO_firre_induced_90min_dox_short_rep1.reverse.bigWig")), `1.5h R2` =
list(pos = file.path(rnaseq_dir,
"ESC_KO_firre_induced_90min_dox_short_rep2.forward.bigWig"), neg =
file.path(rnaseq_dir,
"ESC_KO_firre_induced_90min_dox_short_rep2.reverse.bigWig")), `2h R1` = list(pos
= file.path(rnaseq_dir,
"ESC_KO_firre_induced_120min_dox_short_rep1.forward.bigWig"), neg =
file.path(rnaseq_dir,
"ESC_KO_firre_induced_120min_dox_short_rep1.reverse.bigWig")), `2h R2` =
list(pos = file.path(rnaseq_dir,
"ESC_KO_firre_induced_120min_dox_short_rep2.forward.bigWig"), neg =
file.path(rnaseq_dir,
"ESC_KO_firre_induced_120min_dox_short_rep2.reverse.bigWig")), `2.5h R1` =
list(pos = file.path(rnaseq_dir,
"ESC_KO_firre_induced_150min_dox_short_rep1.forward.bigWig"), neg =
file.path(rnaseq_dir,
"ESC_KO_firre_induced_150min_dox_short_rep1.reverse.bigWig")), `2.5h R2` =
list(pos = file.path(rnaseq_dir,
"ESC_KO_firre_induced_150min_dox_short_rep2.forward.bigWig"), neg =
file.path(rnaseq_dir,
"ESC_KO_firre_induced_150min_dox_short_rep2.reverse.bigWig")), `3h R1` =
list(pos = file.path(rnaseq_dir,
"ESC_KO_firre_induced_180min_dox_short_rep1.forward.bigWig"), neg =
file.path(rnaseq_dir,
"ESC_KO_firre_induced_180min_dox_short_rep1.reverse.bigWig")), `3h R2` =
list(pos = file.path(rnaseq_dir,
"ESC_KO_firre_induced_180min_dox_short_rep2.forward.bigWig"), neg =
file.path(rnaseq_dir,
"ESC_KO_firre_induced_180min_dox_short_rep2.reverse.bigWig")), `3.5h R1` =
list(pos = file.path(rnaseq_dir,
"ESC_KO_firre_induced_210min_dox_short_rep1.forward.bigWig"), neg =
file.path(rnaseq_dir,
"ESC_KO_firre_induced_210min_dox_short_rep1.reverse.bigWig")), `3.5h R2` =
list(pos = file.path(rnaseq_dir,
"ESC_KO_firre_induced_210min_dox_short_rep2.forward.bigWig"), neg =
file.path(rnaseq_dir,
"ESC_KO_firre_induced_210min_dox_short_rep2.reverse.bigWig")), `4h R1` =
list(pos = file.path(rnaseq_dir,
"ESC_KO_firre_induced_240min_dox_short_rep1.forward.bigWig"), neg =
file.path(rnaseq_dir,
"ESC_KO_firre_induced_240min_dox_short_rep1.reverse.bigWig")), `4h R2` =
list(pos = file.path(rnaseq_dir,
"ESC_KO_firre_induced_240min_dox_short_rep2.forward.bigWig"), neg =
file.path(rnaseq_dir,
"ESC_KO_firre_induced_240min_dox_short_rep2.reverse.bigWig")), `4.5h R1` =
list(pos = file.path(rnaseq_dir,
"ESC_KO_firre_induced_270min_dox_short_rep1.forward.bigWig"), neg =
file.path(rnaseq_dir,
"ESC_KO_firre_induced_270min_dox_short_rep1.reverse.bigWig")), `4.5h R2` =
list(pos = file.path(rnaseq_dir,
"ESC_KO_firre_induced_270min_dox_short_rep2.forward.bigWig"), neg =
file.path(rnaseq_dir,
"ESC_KO_firre_induced_270min_dox_short_rep2.reverse.bigWig")), `5h R1` =
list(pos = file.path(rnaseq_dir,
"ESC_KO_firre_induced_300min_dox_short_rep1.forward.bigWig"), neg =
file.path(rnaseq_dir,
"ESC_KO_firre_induced_300min_dox_short_rep1.reverse.bigWig")), `5h R2` =
list(pos = file.path(rnaseq_dir,
"ESC_KO_firre_induced_300min_dox_short_rep2.forward.bigWig"), neg =
file.path(rnaseq_dir,
"ESC_KO_firre_induced_300min_dox_short_rep2.reverse.bigWig")), `5.5h R1` =
list(pos = file.path(rnaseq_dir,
"ESC_KO_firre_induced_330min_dox_short_rep1.forward.bigWig"), neg =
file.path(rnaseq_dir,
"ESC_KO_firre_induced_330min_dox_short_rep1.reverse.bigWig")), `5.5h R2` =
list(pos = file.path(rnaseq_dir,
"ESC_KO_firre_induced_330min_dox_short_rep2.forward.bigWig"), neg =
file.path(rnaseq_dir,
"ESC_KO_firre_induced_330min_dox_short_rep2.reverse.bigWig")))

esc_ko_long_all <- list(`0h R1` = list(pos = file.path(rnaseq_dir,
"ESC_KO_firre_induced_0min_dox_long_rep1.forward.bigWig"), neg =
file.path(rnaseq_dir,
"ESC_KO_firre_induced_0min_dox_long_rep1.reverse.bigWig")), `0h R2` = list(pos =
file.path(rnaseq_dir, "ESC_KO_firre_induced_0min_dox_long_rep2.forward.bigWig"),
neg = file.path(rnaseq_dir,
"ESC_KO_firre_induced_0min_dox_long_rep2.reverse.bigWig")), `0h R3` = list(pos =
file.path(rnaseq_dir, "ESC_KO_firre_induced_0min_dox_long_rep3.forward.bigWig"),
neg = file.path(rnaseq_dir,
"ESC_KO_firre_induced_0min_dox_long_rep3.reverse.bigWig")), `12h R1` = list(pos
= file.path(rnaseq_dir,
"ESC_KO_firre_induced_720min_dox_long_rep1.forward.bigWig"), neg =
file.path(rnaseq_dir,
"ESC_KO_firre_induced_720min_dox_long_rep1.reverse.bigWig")), `12h R2` =
list(pos = file.path(rnaseq_dir,
"ESC_KO_firre_induced_720min_dox_long_rep2.forward.bigWig"), neg =
file.path(rnaseq_dir,
"ESC_KO_firre_induced_720min_dox_long_rep2.reverse.bigWig")), `12h R3` =
list(pos = file.path(rnaseq_dir,
"ESC_KO_firre_induced_720min_dox_long_rep3.forward.bigWig"), neg =
file.path(rnaseq_dir,
"ESC_KO_firre_induced_720min_dox_long_rep3.reverse.bigWig")), `24h R1` =
list(pos = file.path(rnaseq_dir,
"ESC_KO_firre_induced_1440min_dox_long_rep1.forward.bigWig"), neg =
file.path(rnaseq_dir,
"ESC_KO_firre_induced_1440min_dox_long_rep1.reverse.bigWig")), `24h R2` =
list(pos = file.path(rnaseq_dir,
"ESC_KO_firre_induced_1440min_dox_long_rep2.forward.bigWig"), neg =
file.path(rnaseq_dir,
"ESC_KO_firre_induced_1440min_dox_long_rep2.reverse.bigWig")), `24h R3` =
list(pos = file.path(rnaseq_dir,
"ESC_KO_firre_induced_1440min_dox_long_rep3.forward.bigWig"), neg =
file.path(rnaseq_dir,
"ESC_KO_firre_induced_1440min_dox_long_rep3.reverse.bigWig")), `48h R1` =
list(pos = file.path(rnaseq_dir,
"ESC_KO_firre_induced_2880min_dox_long_rep1.forward.bigWig"), neg =
file.path(rnaseq_dir,
"ESC_KO_firre_induced_2880min_dox_long_rep1.reverse.bigWig")), `48h R2` =
list(pos = file.path(rnaseq_dir,
"ESC_KO_firre_induced_2880min_dox_long_rep2.forward.bigWig"), neg =
file.path(rnaseq_dir,
"ESC_KO_firre_induced_2880min_dox_long_rep2.reverse.bigWig")), `48h R3` =
list(pos = file.path(rnaseq_dir,
"ESC_KO_firre_induced_2880min_dox_long_rep3.forward.bigWig"), neg =
file.path(rnaseq_dir,
"ESC_KO_firre_induced_2880min_dox_long_rep3.reverse.bigWig")), `96h R1` =
list(pos = file.path(rnaseq_dir,
"ESC_KO_firre_induced_5760min_dox_long_rep1.forward.bigWig"), neg =
file.path(rnaseq_dir,
"ESC_KO_firre_induced_5760min_dox_long_rep1.reverse.bigWig")), `96h R2` =
list(pos = file.path(rnaseq_dir,
"ESC_KO_firre_induced_5760min_dox_long_rep2.forward.bigWig"), neg =
file.path(rnaseq_dir,
"ESC_KO_firre_induced_5760min_dox_long_rep2.reverse.bigWig")), `96h R3` =
list(pos = file.path(rnaseq_dir,
"ESC_KO_firre_induced_5760min_dox_long_rep3.forward.bigWig"), neg =
file.path(rnaseq_dir,
"ESC_KO_firre_induced_5760min_dox_long_rep3.reverse.bigWig")))
