library(ggplot2)
library(hrbrthemes)
library(drc)
library(nls2)
library(readxl)    
library(reshape2)
library(mgcv)
library(plotly)
library(rempsyc)

app = function(Species, NPs, data) {

data=read.table(data, header=T, sep=",")
names(data)=c("variable", "time")

n=length(unique(data[,1]))
m=length(unique(data[,2]))

# Species="B.subt"
# NPs="TiO2"

gam=readRDS(paste0("/media/thong/sda/RCID/CoHoai/modelling/Rshiny/models/", Species, "_", NPs, ".RDS"))

pred=matrix(predict(gam, data), m, n)

pred=as.data.frame(pred)
names(pred)=unique(data[,1])
pred$time=unique(data[,2])
ldata=melt(pred, id.vars="time")
ldata$log=log2(as.numeric(ldata$variable))

names(ldata)[2]="Concentration (μg/mL)"

p=ggplot(ldata, aes(x=time, y=value, fill=`Concentration (μg/mL)`, color=`Concentration (μg/mL)`)) +
geom_line() +
geom_point() +
xlab("Time (Hours)") + ylab("∆OD 600") + 
scale_colour_manual(values = rainbow(n)) +
theme_minimal() +
theme(text = element_text(size = 20))

pred=pred[,c("time", unique(data[,1]))]
names(pred)[1]="Time"
rownames(pred)=NULL
print(pred)

}