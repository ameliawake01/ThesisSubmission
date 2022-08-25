library(ggplot2)

data1 <- read.csv("/home/amelia-wake/Documents/NHMProject/PossibleDatasets/DivAbTotComm.csv")
ggplot(data=data1, aes(x=RR, y=reorder(Treatment, RR), colour = Metric, group = Metric)) +
  geom_point(size=2.5) +
  xlab("log(Response Ratio)") +
  ylab("Non-Conventional Treatment") +
  geom_errorbar(aes(xmin=Ci.lb, xmax=Ci.ub, colour = Metric), width = 0.25) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", size = 0.5) +
  theme_bw()

data2 <- read.csv("/home/amelia-wake/Documents/NHMProject/ResultsDraft/RandomEffectsDiversityResults.csv")
ggplot(data=data2, aes(x=reorder(Treatment, ResponseRatio), y=ResponseRatio, color = Group, shape = Metric, group = Group)) +
  geom_point(size=3) + 
  xlab("Treatment") +
  ylab("log(Response Ratio)") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.5) +
  theme_bw()

data3 <- read.csv("/home/amelia-wake/Documents/NHMProject/PossibleDatasets/Fixed+GroupsAbundance.csv")
ggplot(data=data3, aes(x=reorder(Soil.Taxonomy, ResponseRatio), y=ResponseRatio, color = Group, group = Group)) +
  geom_point(size=3) + 
  #geom_line(linetype = "dotted") +
  xlab("Soil Taxonomy") +
  ylab("log(Response Ratio)") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.5) +
  theme_bw()
