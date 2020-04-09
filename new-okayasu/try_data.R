require(TDA)
require(myfs)
require(rgl)
require(ks)
require(MaskJointDensity)
require(boot)
require(doParallel)
require(foreach)
require(plotly)
require(viridis)
require(pterrace)

#原稿・スライド用のグラフを出力する
trus11.pl<-calcLandscape(torus11.aggr[["Xdiag"]][[6]])

peak11<-calc.landscape.peak(trus11.pl$Land.dim1, thresh = calcPerMeanNoDouble(torus11.aggr[["Xdiag"]][[6]]), tseq = trus11.pl$tseq, show = T)

plot(trus11.pl$tseq, trus11.pl$Land.dim1, type = "l", col=2, xlab = "(Birth + Death) / 2",ylab = "(Death - Birth) / 2", ylim=c(0, 1), cex.lab=1.6, cex.axis=1.6, lwd=2)
par(new=T)
plot(trus11.pl$tseq, trus11.pl$Land.dim2, type = "l", col=3, xlab = "(Birth + Death) / 2",ylab = "(Death - Birth) / 2", ylim=c(0, 1), cex.lab=1.6, cex.axis=1.6, lwd=2)
