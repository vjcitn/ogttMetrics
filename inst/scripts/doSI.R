
library(ogttMetrics)
data(xoverSamp)
mm = fig3tests(xoverSamp)
#ss = fig3tests(xoverSamp, type="SI")
data(SItests)
par(mfrow=c(1,2))
fig3plot(mm)
fig3plot(SItests, measTag="SI")

