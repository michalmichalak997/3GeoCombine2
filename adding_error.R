library(dplyr)
#Input data points
input<-read.csv(file = ".txt", sep = ";", dec = ".", header = T)
input
nrow(input)

#Adding error to elevation data
error<-rnorm(nrow(input), mean=0, sd=0.05)
error
length(error)
length(input$Z)
disturbed_input<-input
disturbed_input$Z<-input$Z+error
disturbed_input$Scalar<-disturbed_input$Z
write.table(disturbed_input, file=".txt", sep=";", dec=".", row.names=F)

disturbed_input<-dplyr::select(disturbed_input, c("X", "Y", "Z"))
write.table(disturbed_input, file=".txt", sep=" ", dec=".", row.names=F, col.names = F)
