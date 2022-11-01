rnaseq_dir <- "/scratch/Shares/rinn/Michael/firre_timecourse/rnaseq1/results/star_salmon/bigwig/"

esc_ko_short_all_stranded <- list(`0h R1` = list(pos = file.path(rnaseq_dir,
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



esc_ko_short_all <- list(
  `0h R1` = "/scratch/Shares/rinn/Michael/firre_timecourse/rnaseq1/bw/ESC_KO_firre_induced_0min_dox_short_rep1.bigWig",
  `0h R2` = "/scratch/Shares/rinn/Michael/firre_timecourse/rnaseq1/bw/ESC_KO_firre_induced_0min_dox_short_rep2.bigWig",
  `0.5h R1` = "/scratch/Shares/rinn/Michael/firre_timecourse/rnaseq1/bw/ESC_KO_firre_induced_30min_dox_short_rep1.bigWig",
  `0.5h R2` = "/scratch/Shares/rinn/Michael/firre_timecourse/rnaseq1/bw/ESC_KO_firre_induced_30min_dox_short_rep2.bigWig",
  `1h R1` = "/scratch/Shares/rinn/Michael/firre_timecourse/rnaseq1/bw/ESC_KO_firre_induced_60min_dox_short_rep1.bigWig",
  `1h R2` = "/scratch/Shares/rinn/Michael/firre_timecourse/rnaseq1/bw/ESC_KO_firre_induced_60min_dox_short_rep2.bigWig",
  `1.5h R1` = "/scratch/Shares/rinn/Michael/firre_timecourse/rnaseq1/bw/ESC_KO_firre_induced_90min_dox_short_rep1.bigWig",
  `1.5h R2` = "/scratch/Shares/rinn/Michael/firre_timecourse/rnaseq1/bw/ESC_KO_firre_induced_90min_dox_short_rep2.bigWig",
  `2h R1` = "/scratch/Shares/rinn/Michael/firre_timecourse/rnaseq1/bw/ESC_KO_firre_induced_120min_dox_short_rep1.bigWig",
  `2h R2` = "/scratch/Shares/rinn/Michael/firre_timecourse/rnaseq1/bw/ESC_KO_firre_induced_120min_dox_short_rep2.bigWig",
  `2.5h R1` = "/scratch/Shares/rinn/Michael/firre_timecourse/rnaseq1/bw/ESC_KO_firre_induced_150min_dox_short_rep1.bigWig",
  `2.5h R2` = "/scratch/Shares/rinn/Michael/firre_timecourse/rnaseq1/bw/ESC_KO_firre_induced_150min_dox_short_rep2.bigWig",
  `3h R1` = "/scratch/Shares/rinn/Michael/firre_timecourse/rnaseq1/bw/ESC_KO_firre_induced_180min_dox_short_rep1.bigWig",
  `3h R2` = "/scratch/Shares/rinn/Michael/firre_timecourse/rnaseq1/bw/ESC_KO_firre_induced_180min_dox_short_rep2.bigWig",
  `3.5h R1` = "/scratch/Shares/rinn/Michael/firre_timecourse/rnaseq1/bw/ESC_KO_firre_induced_210min_dox_short_rep1.bigWig",
  `3.5h R2` = "/scratch/Shares/rinn/Michael/firre_timecourse/rnaseq1/bw/ESC_KO_firre_induced_210min_dox_short_rep2.bigWig",
  `4h R1` = "/scratch/Shares/rinn/Michael/firre_timecourse/rnaseq1/bw/ESC_KO_firre_induced_240min_dox_short_rep1.bigWig",
  `4h R2` = "/scratch/Shares/rinn/Michael/firre_timecourse/rnaseq1/bw/ESC_KO_firre_induced_240min_dox_short_rep2.bigWig",
  `4.5h R1` = "/scratch/Shares/rinn/Michael/firre_timecourse/rnaseq1/bw/ESC_KO_firre_induced_270min_dox_short_rep1.bigWig",
  `4.5h R2` = "/scratch/Shares/rinn/Michael/firre_timecourse/rnaseq1/bw/ESC_KO_firre_induced_270min_dox_short_rep2.bigWig",
  `5h R1` = "/scratch/Shares/rinn/Michael/firre_timecourse/rnaseq1/bw/ESC_KO_firre_induced_300min_dox_short_rep1.bigWig",
  `5h R2` = "/scratch/Shares/rinn/Michael/firre_timecourse/rnaseq1/bw/ESC_KO_firre_induced_300min_dox_short_rep2.bigWig",
  `5.5h R1` = "/scratch/Shares/rinn/Michael/firre_timecourse/rnaseq1/bw/ESC_KO_firre_induced_330min_dox_short_rep1.bigWig",
  `5.5h R2` = "/scratch/Shares/rinn/Michael/firre_timecourse/rnaseq1/bw/ESC_KO_firre_induced_330min_dox_short_rep2.bigWig"
)

proseq_dir <- "/scratch/Shares/rinn/Michael/firre_timecourse/proseq/results/mapped/rcc_bigwig"

proseq_bw_representative <- list(`0'` = list(pos = file.path(proseq_dir, "JR3116.pos.rcc.bw"),
                                             neg = file.path(proseq_dir, "JR3116.neg.rcc.bw")),
                                 `15'` = list(pos = file.path(proseq_dir, "JR3117.pos.rcc.bw"),
                                              neg = file.path(proseq_dir, "JR3117.neg.rcc.bw")),
                                 `30'` = list(pos = file.path(proseq_dir, "JR3118.pos.rcc.bw"),
                                              neg = file.path(proseq_dir, "JR3118.neg.rcc.bw")))

proseq_bw_all <- list(`0' R1` = list(pos = file.path(proseq_dir, "JR3113.pos.rcc.bw"),
                                     neg = file.path(proseq_dir, "JR3113.neg.rcc.bw")),
                      `0' R2` = list(pos = file.path(proseq_dir, "JR3116.pos.rcc.bw"),
                                     neg = file.path(proseq_dir, "JR3116.neg.rcc.bw")),
                      `0' R3` = list(pos = file.path(proseq_dir, "JR3119.pos.rcc.bw"),
                                     neg = file.path(proseq_dir, "JR3119.neg.rcc.bw")),
                      `15' R1` = list(pos = file.path(proseq_dir, "JR3114.pos.rcc.bw"),
                                      neg = file.path(proseq_dir, "JR3114.neg.rcc.bw")),
                      `15' R2` = list(pos = file.path(proseq_dir, "JR3117.pos.rcc.bw"),
                                      neg = file.path(proseq_dir, "JR3117.neg.rcc.bw")),
                      `15' R3` = list(pos = file.path(proseq_dir, "JR3120.pos.rcc.bw"),
                                      neg = file.path(proseq_dir, "JR3120.neg.rcc.bw")),
                      `30' R1` = list(pos = file.path(proseq_dir, "JR3115.pos.rcc.bw"),
                                      neg = file.path(proseq_dir, "JR3115.neg.rcc.bw")),
                      `30' R2` = list(pos = file.path(proseq_dir, "JR3118.pos.rcc.bw"),
                                      neg = file.path(proseq_dir, "JR3118.neg.rcc.bw")),
                      `30' R3` = list(pos = file.path(proseq_dir, "JR3121.pos.rcc.bw"),
                                      neg = file.path(proseq_dir, "JR3121.neg.rcc.bw")))


atacseq_bw <- list(`0'` =  "/scratch/Shares/rinn/Michael/firre_timecourse/atacseq/results/bwa/mergedLibrary/bigwig/ESC_KO_firre_induced_0_R1.mLb.clN.bigWig",
                   `30'` = "/scratch/Shares/rinn/Michael/firre_timecourse/atacseq/results/bwa/mergedLibrary/bigwig/ESC_KO_firre_induced_30_R1.mLb.clN.bigWig",
                   `60'` = "/scratch/Shares/rinn/Michael/firre_timecourse/atacseq/results/bwa/mergedLibrary/bigwig/ESC_KO_firre_induced_60_R1.mLb.clN.bigWig",
                   `90'` = "/scratch/Shares/rinn/Michael/firre_timecourse/atacseq/results/bwa/mergedLibrary/bigwig/ESC_KO_firre_induced_90_R1.mLb.clN.bigWig",
                   `120'` = "/scratch/Shares/rinn/Michael/firre_timecourse/atacseq/results/bwa/mergedLibrary/bigwig/ESC_KO_firre_induced_120_R1.mLb.clN.bigWig",
                   `150'` = "/scratch/Shares/rinn/Michael/firre_timecourse/atacseq/results/bwa/mergedLibrary/bigwig/ESC_KO_firre_induced_150_R1.mLb.clN.bigWig")

