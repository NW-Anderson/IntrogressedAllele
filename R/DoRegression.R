data <- read.csv(file = 'offspring info.csv')
data <- data[!is.na(data[,1]),]

summary(m1 <- zeroinfl(Calves.sired ~ Year.of.rut | Sire.age..years., data = data))
