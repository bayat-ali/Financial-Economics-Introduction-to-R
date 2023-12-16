library("xts")
library("quantmod")

hiweb = read.csv("/Users/Shared/Files From e.localized/Desktop/Codes/R/S_HiWeb.csv")
hiweb1 = hiweb[c(3:6,8,12)]
class(hiweb$X.DTYYYYMMDD.)
x = hiweb$X.DTYYYYMMDD.
x <- data.frame(lapply(x, function(x) as.Date(as.character(x), "%Y%m%d")))
hiweb1$Date = t(x)

hiweb1$Date<-as.Date(hiweb1$Date, format="%Y-%m-%d")
hiweb1 <- xts(hiweb1[,-7], order.by=hiweb1$Date)
colnames(hiweb1) = c("Hiweb.Open","Hiweb.High","Hiweb.Low",
                     "Hiweb.Close","Hiweb.Volume","Hiweb.Adjusted")
chartSeries(hiweb1, name = "hiweb1", minor.ticks=FALSE, theme=chartTheme("white"))


foolad = read.csv("/Users/Shared/Files From e.localized/Desktop/Codes/R/foolad.csv")
foolad1 = foolad[c(3:6,8,12)]
class(foolad$X.DTYYYYMMDD.)
x = foolad$X.DTYYYYMMDD.
x <- data.frame(lapply(x, function(x) as.Date(as.character(x), "%Y%m%d")))
foolad1$Date = t(x)

foolad1$Date<-as.Date(foolad1$Date, format="%Y-%m-%d")
foolad1 <- xts(foolad1[,-7], order.by=foolad1$Date)
colnames(foolad1) = c("Hiweb.Open","Hiweb.High","Hiweb.Low",
                     "Hiweb.Close","Hiweb.Volume","Hiweb.Adjusted")
chartSeries(foolad1, name = "foolad1", minor.ticks=FALSE, theme=chartTheme("white"))

