proj <-"firre_timecourse"
user <- Sys.getenv("LOGNAME")  # identikey
sharedir <- '/scratch/Shares/rinn'
groupdir <- '/rinnlab/rinngrp/'

users <- tribble(~logname,~dirname,
                 "mism6893","Michael",
                 "laha3063","laurette/proj")

mydir <- users %>% filter(logname==user) %>% select(dirname)

