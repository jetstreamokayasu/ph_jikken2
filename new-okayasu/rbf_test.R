#test of RBF
# NOT RUN {
demo(rbf_irisSnnsR)
# }
# NOT RUN {
demo(rbf_sin)
# }
# NOT RUN {
demo(rbf_sinSnnsR)
# }
# NOT RUN 
{

inputs <- as.matrix(seq(0,10,0.1))
outputs <- as.matrix(sin(inputs) + runif(inputs*0.2))
outputs <- normalizeData(outputs, "0_1")

model <- rbf(inputs, outputs, size=40, maxit=1000, 
             initFuncParams=c(0, 1, 0, 0.01, 0.01), 
             learnFuncParams=c(1e-8, 0, 1e-8, 0.1, 0.8), linOut=TRUE)

par(mfrow=c(2,1))
plotIterativeError(model)
plot(inputs, outputs)
lines(inputs, fitted(model), col="green")
# 
}

model.torus <- rbf(tseq.torus,  Land.dim1.torus7, size=40, maxit=1000, 
             initFuncParams=c(0, 1, 0, 0.01, 0.01), 
             learnFuncParams=c(1e-8, 0, 1e-8, 0.1, 0.8), linOut=TRUE)
plot(tseq.torus, Land.dim1.torus7, type = "l", col=2, xlab = "(Birth + Death) / 2",ylab = "(Death - Birth) / 2", ylim=c(0, 1.1))
lines(tseq.torus, fitted(model.torus), col="green")

