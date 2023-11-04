library(ggplot2)
library(hrbrthemes)
library(drc)
library(nls2)
library(readxl)    
library(reshape2)
library(mgcv)
library(plotly)
library(rempsyc)

set.seed(123)
rm(list = ls())

p="/media/thong/sda/RCID/CoHoai/modelling/Rshiny/scripts"

setwd(p)

dir.create("../models", recursive=T)

multiplesheets <- function(fname) {
   
  # getting info about all excel sheets
  sheets <- readxl::excel_sheets(fname)
  tibble <- lapply(sheets, function(x) readxl::read_excel(fname, sheet = x))
  data_frame <- lapply(tibble, as.data.frame)
    
  # assigning names to data frames
  names(data_frame) <- sheets
    
  # print data frame
  print(data_frame)
}

### TiO2
# specifying the path name
path <- "../data/Dilnaz Bsubt-Ecoli TiO2 OD600 96wp shaking_20230226_185953_results_statistics.xls"
xls=multiplesheets(path)

# Bacillus subtilis
dt=xls$`mean_B-subt`
dt$time=(as.numeric(rownames(dt))*15-15)/60
# normalize
dt=data.frame(apply(dt,2, function(x) x-x[1]), check.names=F)
dt=dt[dt$time%in%c(1:20),]

df=melt(dt, id.var='time')
df$variable=as.numeric(as.character(df$variable))
df$log=log2(as.numeric(df$variable))

fit=gam(value ~ te(variable, time), data = df)
saveRDS(fit, "../models/B.subt_TiO2.RDS")

# E. coli
dt=xls$`mean_E-coli`
dt$time=(as.numeric(rownames(dt))*15-15)/60
# normalize
dt=data.frame(apply(dt,2, function(x) x-x[1]), check.names=F)
dt=dt[dt$time%in%c(1:20),]

df=melt(dt, id.var='time')
df$variable=as.numeric(as.character(df$variable))
df$log=log2(as.numeric(df$variable))

fit=gam(value ~ te(variable, time), data = df)
saveRDS(fit, "../models/E.coli_TiO2.RDS")

### SNCD
# specifying the path name
path <- "../data/Meruyert Bsubt-Ecoli SNCD OD600 shaking 96wp 37d_20231009_160126_results_statistics.xls"
xls=multiplesheets(path)

# Bacillus subtilis
dt=xls$`mean_B.subt`
dt$time=(as.numeric(rownames(dt))*15-15)/60
# normalize
dt=data.frame(apply(dt,2, function(x) x-x[1]), check.names=F)
dt=dt[dt$time%in%c(1:20),]

df=melt(dt, id.var='time')
df$variable=as.numeric(as.character(df$variable))
df$log=log2(as.numeric(df$variable))

fit=gam(value ~ te(variable, time), data = df)
saveRDS(fit, "../models/B.subt_SN-CNPs.RDS")

# E. coli
dt=xls$`mean_E.coli`
dt$time=(as.numeric(rownames(dt))*15-15)/60
# normalize
dt=data.frame(apply(dt,2, function(x) x-x[1]), check.names=F)
dt=dt[dt$time%in%c(1:20),]

df=melt(dt, id.var='time')
df$variable=as.numeric(as.character(df$variable))
df$log=log2(as.numeric(df$variable))

fit=gam(value ~ te(variable, time), data = df)
saveRDS(fit, "../models/E.coli_SN-CNPs.RDS")
