rm(list = ls())  # clean workspace

library(foreign)


LDage <-  read.table('LdAWeightedRegTr5.pc')
agefless <- read.table('LdAFlessAvs10ind.pc')
time <- 1:40


par(mfrow=c(1,1))
plot(time, agefless[time,5]/agefless[1,5],ylim=c(0.9,11),ylab='Employment relative to Age=1',xlab='Age',pch=21,bg="black",log='y')
points(time, LDage[time,6]/LDage[1,6],pch=21,bg="darkgray")
points(time, LDage[time,7]/LDage[1,7],pch=8,bg="black")
points(time, LDage[time,18]/LDage[1,18],pch=5,bg="darkgray")
legend(1,6,legend=c('USA   ','France   ','Ghana   ','India   '),pch=c(21,21,8,5),pt.bg=c('black','darkgray','black','darkgray'))
grid(nx = NULL, ny = NULL, col = "gray", lty = "dotted")
