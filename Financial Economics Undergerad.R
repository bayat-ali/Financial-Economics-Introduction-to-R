# Financial Economics
# 2021-Fall
# Introduction to R
# Ali Bayat

rm(list=ls(all=TRUE))
ls()

# Arithmetic operator and variables
x = 10
x
print(x)
show(x)


x <- 10
x += 1   # It does not work like Python
x = x + 1
x <- x + 1
x

100 + 200 

5^2

(5/3 + 1/500)^50

19 %% 3

19 %/% 3

25 %/% 4

exp(1)

options(digits = 10)
exp(1)

options(digits = 20)
exp(1)

options(digits = 6)

sin(pi)
sin(pi/6)
sin(0)

sqrt(9)

x <- 10.2546
x

floor(x)
ceiling(x)

x = 10.999
floor(x)
ceiling(x)

z <- 5*x+9
(z <- 5*x+9)

# functions
seq(from = 5, to = 50, by = 5)
seq()
seq(from = 5, to = 50)
seq(5, 50, 10)

?seq
help(seq)

seq(by = 5, to = 50, from = 5)
seq(by = -4, from = 16, to = 4)

x = 8
seq(by = -4, from = 16, to = x)


# vector
(x <- seq(5, 10, 2))
y <- seq(6, 0, 3);y
(s <- rep(5, 10))
(s <- rep(8, 3))
(n <- c(s,x))

(x <- seq(5, 20, 2))
x[5]
i <- c(5, 7, 2)
x[i]

j <- c(-3, -4, -8) #Different from Python
x[j]
x

x[1] <- 0
x

length(x)

(x <- c(5, 10, 20))
(y <- c(1, 3, 9))
x+y
x*y
x^y

c(2, 5, 6, 8, 3, 9) + c(1, 3, 5)
c(2, 5, 6, 8, 3, 9) * c(1, 3, 5)
c(2, 5, 6, 8, 3) + c(1, 3, 5)


10 + c(2, 5, 8)
10 * c(2, 5, 8)

sqrt(c(2, 8, 10))
prod(c(2, 8, 10)) #prod is abbreviation of product

1:10  # default by 1
5:15
min(100:1000)
max(100:1000)
sort(c(1200, 1960, 5, 77))

var(1:10)
(n <- c(1.5, 2.0, 2.4, 3.6))
(n.mean <- sum(n)/length(n))
n.mean - mean(n)

n.var <- sum((n-n.mean)^2)/length(n)
n.var == var(n)

n.var <- sum((n-n.mean)^2)/(length(n)-1)
n.var == var(n)

# Calculate a simple integral
dt <- 0.01
t <- seq(0, pi/6, by = dt)
f <- cos(t)
sin(pi/6)
(I <- sum(f)*dt)
I - sin(pi/6)

dt <- 0.001
t <- seq(0, pi/6, by = dt)
f <- cos(t)
sin(pi/6)
(I <- sum(f)*dt)
I - sin(pi/6)

# plots
(x <- seq(0, 150, 10))
(y <- (15 + 1/x)^(1/x))
plot(x, y)


# Working with data
n <- NA
is.na(n)
n <- c(50, NA, 200);n
(n <- c(50, NA, 200))

n
is.na(n)
any(is.na(n))
mean(n)
mean(n, na.rm = TRUE)


# Logical expressions
5 == 6
5 == 5
6 != 5
c(0, 0, 0, 1) | c(1, 0, 0, 1)
xor(c(0, 0, 0, 1) , c(1, 0, 0, 1)) 

x <- 1:20
x %% 4 == 0
(z <- x[x %% 4 == 0])
x <- c(2, NA, 3, 4)
x > 2
x[x>2]
subset(x, subset = x > 2) # does not consider NA

(x <- 10:20)
which(x %% 4 == 0)

# Approximation problem in R
2*2 == 4
sqrt(2) * sqrt(2) == 2

all.equal(sqrt(2) * sqrt(2), 2)

# sequential operator
x <- 0
x * sin(1/x) == 0
(x == 0) | (sin(1/x) == 0)
(x == 0) || (sin(1/x) == 0) # without warning


# Matrix
matrix(data = NA, nrow = 1, ncol = 1, byrow = FALSE)
(A <- matrix(data = 1:6, nrow = 2, ncol = 3, byrow = FALSE))
(B <- matrix(data = 1:6, nrow = 2, ncol = 3, byrow = TRUE))
(C <- matrix(data = 1:6, nrow = 3, ncol = 2, byrow = TRUE))

dim(A)
rbind(A, B)
cbind(A, B)
rbind(A, C)

A
A[1,2]
A[1,]
A[,2:3]

(D <- diag(c(4, 5, 6, 7, 8)))

(A <- matrix(c(3, 5, 6, 2), nrow = 2, ncol = 2))  

(B <- matrix(c(1, 2, 3, 5), nrow = 2, ncol = 2))

A + B
A * B
A %*% B

(A.inv <- solve(A))
A * A.inv
A %*% A.inv
all.equal(A %*% A.inv, matrix(c(1,0,0,1), 2, 2))

a <- solve(A, B)
B
A %*% a
all.equal(A %*% a, B)

# Fundamentals of Programming for Computing in R
# find zeros of a*x^2 + b*x + c = 0

# clear the workspace
rm(list = ls())

# input
a <- 5
b <- 8
c <- 2

# calculation
root1 <- (-b + sqrt(b^2 - 4*a*c))/(2*a)
root2 <- (-b - sqrt(b^2 - 4*a*c))/(2*a)

# output
print(c(root1,root2))


# find zeros of a*x^2 + b*x + c = 0

# clear the workspace
rm(list = ls())

# input
a <- 1
b <- 4
c <- 2

# calculation the delta
delta <- b^2 - 4*a*c

if (delta > 0){
  roots <- c((-b + sqrt(b^2 - 4*a*c))/(2*a),
             (-b - sqrt(b^2 - 4*a*c))/(2*a))
} else {
  if (delta == 0){
    roots <- -b/2*a
  } else {
    roots <- c()
  }
}

# output
show(roots)



# for loop
x_l <- seq(1, 20, 3)
sum_x <- 0
for (i in x_l){
  sum_x <- sum_x + i
  cat ("elemane konooni halghe:", i , "\n")
  cat ("        ta alan majmoo:", sum_x, "\n")
}

sum(x_l)


x_l <- 0:100
sum_x <- 0
for (i in x_l){
  sum_x <- sum_x + i
  cat ("az", i,"ta", i+1, "dar dore", i+1,"th hastim" , "\n")
}
 

# calculate n factorial

# clear the workspace
rm(list = ls())

# input
n <- 8 

# calcualtion 
n_factorial <- 1
for (i in 1:n){
  n_factorial <- n_factorial*i
}

# out
show(n_factorial)


# preallocation
n <- 10000000
x <- rep(0, n)
for (i in 1:n){
  x[i] <- i
}

# redimensioning
n <- 10000000
x <- 1
for (i in 2:n){
  x[i] <- i
}

# while loop
# calculate the first Fibonacci number greater than 150

# clear the workspace
rm(list = ls())

# Initialisie variables
F <- c(1, 1)
n <- 2

while (F[n] <= 150) {
  n <- n + 1
  F[n] <- F[n-1] + F[n-2]
  cat(n, "omin jomle fib barabar ast ba",F[n],"\n")
}
 
# output


cat("avalin adade fib > 150: F(", n, ") =", F[n], "\n")


# modat zamaen bazpardakhte vam

# clear the workspace
rm(list = ls())

# inputs
r <- 0.04 #sood
p <- 1/12 #zamane beyne har pardakht
b <- 20000000 #mablaghe vam
pardakhti <- 120000 #mablaghe pardakhti mahiane

# calculation
time <- 0
bedehi <- b
while (bedehi > 0) {
  time <- time + p
  bedehi <- bedehi*(1+r*p) - pardakhti
}

# output
cat("modat zamane pardakhte vame shoma", time, "sal\n")

# add the phyton projects

for (i in 1:3){  # Unlike Python in R, the last number itself is included.
  cat ("h", "\n")
}

for (i in 1:3){
  print ("h")
}





###################
### new session ###
###################

n <- 100
m <- 0
for (i in 1:n) {
  m <- m + i^2
}
show(m)

sum((1:n)^2)

# ifelse
ifelse(test, A, B)  

x <- c(-12, -9, 15, 14)
ifelse(x > 0, "Positive", "Negative")


# pmin & pmax
pmin(c(5, 8, 12), c(1, 29, 30), c(9, 1, 55))
pmax(c(5, 8, 12), c(1, 29, 30), c(9, 1, 55))

# program flow: What does our program do at each stage?
x <- 3
for (i in 1:3) {
  show(x)
  if(x %% 2 == 0){
    x <- -x/2
  }else {
    x <- 3*x + 1
  }
}
show(x)

# Turn Warning into an error
options(warn = 2)

# Debugging
x <- 3
for (i in 1:3) {
  show(x)
  cat("i = ", i, "\n")
  if(x %% 2 == 0){
    x <- -x/2
  }else {
    x <- 3*x + 1
  }
}
show(x)

# text
x <- "pejman"
y <- "ali"
z <- "javad"
(s <- paste(x, y, z, sep = ", "))



# format
# neshan dadane tavan haye 1 ta n asad x

# input
x <- 6
n <- 5

# neshan dadane tavan ha
cat("tavan haye", x, "\n")
cat("shomarandeh         meghdar\n\n")

meghdar <- 1
for (i in 1:n) {
  meghdar <- meghdar * x
  cat(format(i, width = 8),
      format(meghdar, width = 11),
      "\n", sep = "")
}

?format
zz <- data.frame("(row names)"= c("aaaaa", "b"), check.names = FALSE)
format(zz)
format(zz, justify = "left")

## use of nsmall
format(13.7)
format(13.7, nsmall = 3)
format(c(6.0, 13.1), digits = 2)
format(c(6.0, 13.1), digits = 2, nsmall = 1)

## use of scientific
format(2^31-1)
format(2^31-1, scientific = TRUE)
format(2^31-1, digits = 2, scientific = TRUE)

# same as above:
cat(paste(format(1: n, width = 8),
          format(x^(1:n), width = 10), "\n"), sep = "")



# input from a file;
?scan
getwd()
setwd("C:/Users/Ilia/Desktop/codes/R")

# mohasebeye miane va chrak ha
# clear the workspace
rm(list = ls())

# input
file_name = "C:/Users/Ilia/Desktop/codes/R/data1.txt"
# dade ha ya bayad ba space az ham joda shan ua new line

# read from file
data <- scan(file = file_name)

# calculations
n <- length(data)
data.sort <- sort(data)
data.1qrt <- data.sort[ceiling(n/4)]
data.med <- data.sort[ceiling(n/2)]
data.3qrt <- data.sort[ceiling(3*n/4)]

# output
cat("charake aval:", data.1qrt, "\n")
cat("miane:", data.med, "\n")
cat("charake sevom", data.3qrt, "\n")

source("C:/Users/Ilia/Desktop/codes/R/prcalmedsqrts.R")
source("C:/Users/Ilia/Desktop/codes/R/prcalmedsqrts.txt")


# Input from the keyboard;
?readline
# find the zeros of a*x^2 + b*x + c = 0

# clear the workspace
rm(list = ls())

# input
cat("peyda kardane sefr haye a*x^2 + b*x + c = 0")
a <- as.numeric(readline("a = "))#function readline data ra be soorate character mikhoone
b <- as.numeric(readline("b = "))
c <- as.numeric(readline("c = "))

# calculate delta
delta <- b^2 - 4*a*c
if (delta > 0){
  roots <- (-a + c(1, -1) * sqrt(delta))/(2*a)
}else {
  if (delta == 0){
    roots <- -b/(2*a)
  }else{
    roots <- c()
  }
}

# output
if (length(roots) == 0){
  cat("rishe haghighi nadarad")
} else if (length(roots) == 1){
  cat("rishe tekrari: ", roots, "\n")
} else {
  cat("2 rishe: ", roots[1], "va", roots[2], "\n")
}

source("C:/Users/Ilia/Desktop/codes/R/calculatetheroots.txt")
getwd()
setwd("C:/Users/Ilia/Desktop/codes/R")
getwd()
source("calculatetheroots.txt")



# output to a file;
?write


(x <- matrix(1:24, nrow = 4, ncol = 6))
write(x, file = "matrix01.txt", ncolumns = 6)
write(t(x), file = "matrix02.txt", ncolumns = 6) #t(x): transpose of x

cat(c(1:4), file = "matrix01.txt", append = TRUE) # az akhar michasboone

x <- matrix(rep(1:5, 1:5), nrow = 3, ncol = 5)
dump("x", file = "matrix03.txt")
source("matrix03.txt")
x

(a <- rep(1:5, 1:5))



# plotting
x <- seq(0, 5, by = 0.01)
y.upper <- 2*sqrt(x)
y.lower <- -2*sqrt(x)
y.max <- max(y.upper)
y.min <- min(y.lower)
plot(c(-3, 5), c(y.min, y.max), type = "n", xlab = "x", ylab = "y")
lines(x, y.upper)
lines(x, y.lower)
abline(v = -1)
points(1, 0)
text(1, 0, "focus(1, 0)")
text(1, 0, "focus(1, 0)", pos = 1)
text(1, 0, "focus(1, 0)", pos = 2)
text(1, 0, "focus(1, 0)", pos = 3)
text(1, 0, "focus(1, 0)", pos = 4)

text(-1, y.min, "x = -1")
text(-1, y.min, "x = -1", pos = 4)
title("parabola y^2 = 4*X")

par(mfrow = c(2,2)) #or mfcol
curve(x*sin(x), from = 0, to = 100, n = 1001)
curve(x*sin(x), from = 0, to = 100, n = 27)
curve(x*sin(x), from = 0, to = 1, n = 1001)
curve(x*sin(x), from = 0, to = .100, n = 1001)



###################
### new session ###
###################

# Functions;

quad <- function(a, b, c){
  #find the zeros of a*x^2 + b*x + c = 0
  if (a == 0 && b == 0 && c == 0){
    roots <- NA
  } else if (a == 0 && b == 0) {
    roots <- NULL
  } else if (a == 0) {
    roots <- -c/b
  } else {
    # calculation of delta
    delta <- b^2 - 4 * a * c
    if (delta > 0) {
      roots <- (-b + c(1, -1) * sqrt(delta))/(2*a)
    } else if (delta == 0) {
      roots <- -b/(2*a)
    } else {
      roots <- NULL
    }
  } 
  return(roots)
}
quad(1, 0, -1)
quad(1, 5, -1)
quad(0, 0, 0)
quad(1, 10, 100)
quad(0, 10, 0)


# Permutation & Combination
n_factorial <- function(n){
  #calculate n factorial
  n_fact <- prod(1:n)
  return(n_fact)
}

n_choose_r <- function(n, r){
  #calculate n choose r
  n_ch_r <- n_factorial(n)/(n_factorial(r)*n_factorial(n-r))
  return(n_ch_r)
}

n_choose_r(4, 2)
n_choose_r(6, 4)
n_choose_r(6, 44)


# trimmed mean: miangine piraste
w_mean <- function(x, k) {
  # calculate the k-th winsoried mean of the vector x
  x <- sort(x)
  n <- length(x)
  x[1:k] <- x[k+1]
  x[(n-k+1):n] <- x[k+1]
  return(mean(x))
}

x <- c(3, 45, 53, 109, 26, 61, 58, 86, 34, 77)
mean(x)
w_mean(x, 2)
x.error <- x
x.error[1] <- 1000
mean(x.error)
w_mean(x.error, 2)


# scope: motaghayer haye daroone function tasiri bar mothaghayer haye biruni nadarand
test <- function(s){
  f <- s + 1 
  return(f)
}
test(1)
f = 2
test(f)
s
f


# arguments
formals(w_mean)
formals(test)

test <- function(x = 1){
  return(x)
}
test(5)
test()

fcn <- function(argus = c("book", "car", "cat")){
  argus <- match.arg(argus)
  return(argus)
}
fcn()
fcn("car")
fcn("carcar")



# Vector based programing using function
?apply
?sapply

#estimate the density of primes
#clear the workspace
rm(list = ls())

prime <- function(n) {
  #return TRUE if n is prime
  if (n == 1) {
    is.prime <- FALSE
  } else if (n == 2) {
    is.prime <- TRUE
  } else {
    is.prime <- TRUE
    for (m in 2:(n/2)) {
      if (n %% m == 0) is.prime <- FALSE
    }
  }
  return(is.prime)
}

#input
n <- 2000

#calculate the number of primes <= m for m in 2:n
#https://drive.google.com/file/d/1UdA24jAapC9LsdFzMZW4aB0ONADd8fM-/view?usp=sharing
m.vec <- 2:n
primes <- sapply(m.vec, prime)
num.primes <- cumsum(primes)

#output
par(mfrow = c(1, 2), las = 1)   # las miad label haye mehvar haro movazi mizare ba mehvar
plot(m.vec, num.primes/m.vec, type = "l",  # num.primes/m.vec chon be donbale density hastim
     main = "prime density", xlab = "n", ylab = "")
lines(m.vec, 1/log(m.vec), col = "red")

plot(m.vec, (num.primes/m.vec)*log(m.vec), type = "l",
     main = "prime density * log(n)", xlab = "n", ylab = "")


# agar n ro ziad konim mohasebat be shedat toolani mishe ->
# behbood code ghabli:
#estimate the density of primes

#clear the workspace
rm(list = ls())

prime <- function(n) {
  #return TRUE if n is prime
  if (n == 1) {
    is.prime <- FALSE
  } else if (n == 2) {
    is.prime <- TRUE
  } else {
    is.prime <- TRUE
    m <- 2
    m.max <- sqrt(n)
    while (is.prime && m <= m.max) {
      if (n %% m == 0) is.prime <- FALSE
      m <- m + 1
    }
  }
  return(is.prime)
}

#input
n <- 200000

#calculate the number of primes <= m for m in 2:n
m.vec <- 2:n
primes <- sapply(m.vec, prime)
num.primes <- cumsum(primes)

#output
par(mfrow = c(1, 2), las = 1)
plot(m.vec, num.primes/m.vec, type = "l", 
     main = "prime density", xlab = "n", ylab = "")
lines(m.vec, 1/log(m.vec), col = "red")

plot(m.vec, num.primes/m.vec*log(m.vec), type = "l",
     main = "prime density * log(n)", xlab = "n", ylab = "")



# Recursive programing
nfact_n <- function(n) {
  #calculate n factorial
  if (n == 1){
    cat("farakhanie nfact_n(1)\n")
    return(1)
  } else {
    cat("farakhanie nfact_n(", n, ")\n", sep= "")
    return(n * nfact_n(n-1))
  }
}
nfact_n(7)



# Debugging functions
my_fun <- function(x){
  browser()
  y <- x * z
  return(y)
}

my_fun(c(1, 2, 3))



# Data structures in R
# factors
# categorical data
hair <- c("black", "gray", "gray", "brown", "blond", "brown")
is.character(hair)
is.factor(hair)
hair <- factor(hair)
is.factor(hair)
levels(hair)

hair <- factor(hair, levels = c("black", "gray", "brown", "blond", "none"))
table(hair)


# ordinal data
phys_act <- c("L", "H", "L", "L", "H", "M", "M")
phys_act <- factor(phys_act, levels = c("L", "M", "H"), ordered = TRUE)
is.ordered(phys_act)
phys_act[2] > phys_act[1]

phys_act <- factor(phys_act, levels = c("L", "M", "H"),
                   labels = c("Low", "Medium", "High" ), ordered = TRUE)
table(phys_act)


# Dataframes
?read.table
setwd("C:/Users/Ilia/Desktop/codes/R")
ufc <- read.csv("C:/Users/Ilia/Desktop/codes/R/ufc.csv")
head(ufc)
tail(ufc)
x <- ufc$height.m
x[1:6]
ufc[[5]]
ufc[["height.m"]]
head(ufc)
ufc[1:6, 5]
ufc[1:6, 4]

d_h <- ufc[4:5]
d_h[1:5,]
d_h
d.h <- ufc[c("dbh.cm", "height.m")]
d.h[1:5,]

is.data.frame(d.h)
mode(ufc)
mode(ufc[5])
mode(ufc[[5]])

ufc$volume.m3 <- pi * (ufc$dbh.cm/200)^2 * ufc$height.m # /200 baraye ine ke ye taghsim bar 2 baraye shoa va ye taghsim bar 100 baraye inke mikhaim bar asase meter mohasebe konim
# V=??*r^2*h
mean(ufc$volume.m3)

(ufc.names <- names(ufc))
names(ufc) <- c("P", "T", "S", "D", "H", "V")
names(ufc)

names(ufc) <- ufc.names
names(ufc)

head(ufc)
t.height <- subset(ufc, subset = species %in% c("DF", "GF"),
                   select = c(plot, tree, height.m))
head(t.height)
ufc

?write.table



# Lists
my.list <- list("one", TRUE, c("f", "h", "t"), c(1, 5, 9))
my.list
my.list[2]
my.list[[2]]
mode(my.list[[2]])
mode(my.list[2])
my.list[4][1]
my.list[4]
my.list[[4]][1]
my.list[[4]][2]
my.list[4][2]  # in code khoruju null mide chon 4th ozve my.list khodesh 
# list nist ke 2th misho bede, bayad daghigh begim va 2 ta [ bezarim: [[]]

mylist <-list(first = "one", second = TRUE, third = 3, tourth = c(5, 7))
names(mylist)
mylist$second

x <- list(1, c(2,5), c(3,76,7))
unlist(x)  # list migire va bordar mide

?lm
lm.xy <- lm(y ~ x, data = data.frame(x = 1:5, y = 1:5))
mode(lm.xy)
names(lm.xy)
ghahraman <- list(
  A = c(2000, 2005),
  B = c(2001, 2004, 2010),
  C = c(2002, 2003),
  D = c(2006, 2007, 2008, 2009),
  E = c())

year <- 2005
for (i in 1:length(ghahraman)) {
  if (year %in% ghahraman[[i]]) {
    winner <- names(ghahraman)[i]
  }
}
winner


# tapply
?tapply(vector, index, function)
head(ufc)
tapply(ufc$height.m, ufc$species, mean)
tapply(ufc$height.m, ufc$species, length)



###################
### new session ###
###################
# Graphical parameters: par;
# Graphical augments;

setwd("C:/Users/Ilia/Desktop/codes/R")
ufc <- read.csv("C:/Users/Ilia/Desktop/codes/R/ufc.csv")
plot(ufc$dbh.cm, ufc$height.m)
plot(ufc$dbh.cm, ufc$height.m, xlab = "Diameter (cm)", ylab = "Height(m)")

?plot
?par

opar1 <- par(las = 1, mar = c(4, 4, 3, 2)) #las = 1 yani onvane mehvare y ro ofoghi bezar
plot(ufc$dbh.cm, ufc$height.m, axes = FALSE, xlab = " ", ylab = " ", type = "n")
plot(ufc$dbh.cm, ufc$height.m, axes = FALSE, xlab = " ", ylab = " ")
plot(ufc$dbh.cm, ufc$height.m, axes = FALSE, xlab = " ", ylab = " ", type = "n")
points(ufc$dbh.cm, ufc$height.m, 
       col = ifelse(ufc$height.m > 15, "darkseagreen4", "red"),
       pch = ifelse(ufc$height.m > 15, 1, 3))

points(ufc$dbh.cm, ufc$height.m, 
       col = ifelse(ufc$height.m > 20, "darkseagreen4", "red"),
       pch = ifelse(ufc$height.m > 15, 1, 10))
axis(1)
axis(2)
opar2 <- par(las = 0)

mtext("Diameter(cm)", side = 1, line = 3)
mtext("Diameter(cm)", side = 1, line = 2)
mtext("Height(m)", side = 2, line = 2)
mtext("Height(m)", side = 2, line = 3)
box()

legend(x = 60, y = 17, c("Normal trees", "weird tree"),
       col = c("darkseagreen", "red"), 
       pch = c(1, 10),
       bty = "n")
legend(x = 60, y = 17, c("Normal trees", "weird tree"),
       col = c("darkseagreen", "red"), 
       pch = c(1, 3))

legend(x = 60, y = 17, c("Normal trees", "weird tree"),
       col = c("darkseagreen", "red"), 
       pch = c(4, 11))



# Mathematical typesetting;
curve(100 * (x^3 - x^2) + 15, from = 0, to = 1,
      xlab = expression(alpha),
      ylab = expression(100 %*% (alpha^3 - alpha^2) + 15),
      main = expression(paste("Function: (", 
                              f(alpha) == 100 %*% (alpha^3 - alpha^2) + 15)))

mtext(expression(100 %*% (alpha^3 - alpha^2) + 15), side = 2, line = 2.5)


# chejoori label mehvare y haro biarim yekam rast
curve(100 * (x^3 - x^2) + 15, from = 0, to = 1,
      xlab = expression(alpha),
      ylab = "",
      main = expression(paste("Function: (", 
                              f(alpha) == 100 %*% (alpha^3 - alpha^2) + 15)))

mtext(expression(100 %*% (alpha^3 - alpha^2) + 15), side = 2, line = 2.5)

myMu <- 0.5
mysigma <- 0.25
curve(100 * (x^3 - x^2) + 15, from = 0, to = 1,
      xlab = expression(alpha),
      ylab = expression(100 %*% (alpha^3 - alpha^2) + 15),
      main = expression(paste("Function: (", 
                              f(alpha) == 100 %*% (alpha^3 - alpha^2) + 15)))

par(usr = c(0, 1, 0, 1))
text(0.1, 0.1, bquote(sigma[alpha] == .(mysigma)), cex = 1.25)
text(0.1, .5, bquote(sigma[alpha] == .(mysigma)), cex = 1)
text(0.6, 0.6, paste("(The mean is ", myMu, ")", sep = ""), cex = 1)
text(0.5, 0.9, bquote(paste("sigma^2 = ",
                            sigma^2 == .(format(mysigma^2, 2)))))



# Permanence
setwd("C:/Users/Ilia/Desktop/codes/R")
ufc <- read.csv("C:/Users/Ilia/Desktop/codes/R/ufc.csv")
head(ufc)

pdf(file = "gragh1.pdf", width = 4, height = 3)
plot(ufc$dbh.cm, ufc$height.m, main = "UFC trees",
     xlab = "Dbh.cm", ylab = "Height (cm)")
dev.off()



# Grouped graphs;
library(lattice)
densityplot(~ dbh.cm | species, data = ufc)
bwplot(~ dbh.cm | species, data = ufc)
histogram(~ dbh.cm | species, data = ufc)
xyplot(height.m ~ dbh.cm | species, data = ufc)


# 3D plots;
setwd("C:/Users/Ilia/Desktop/codes/R")
ufc.plots <- read.csv("C:/Users/Ilia/Desktop/codes/R/ufc-plots.csv")
str(ufc.plots)
head(ufc.plots)

contourplot(vol.m3.ha ~ east * north,
            main = expression(paste("Volume (", m^3, ha^{-1}, ")", sep = "")),
            xlab = "East (m)", ylab = "North (m)",
            region = TRUE,
            aspect = "iso",
            col.regions = gray((11:1)/11),
            data = ufc.plots)


contourplot(vol.m3.ha ~ east * north,
            main = expression(paste("Volume (", m^3, ha^{-1}, ")", sep = "")),
            xlab = "East (m)", ylab = "North (m)",
            region = TRUE,
            aspect = "iso",
            col.regions = gray((11:1)/11),
            data = ufc.plots)
contourplot(vol.m3.ha ~ east * north,
            main = expression(paste("Volume (", m^3, ha^{-1}, ")", sep = "")),
            xlab = "East (m)", ylab = "North (m)",
            region = FALSE,
            aspect = "iso",
            col.regions = gray((11:1)/11),
            data = ufc.plots)
contourplot(vol.m3.ha ~ east * north,
            main = expression(paste("Volume (", m^3, ha^{-1}, ")", sep = "")),
            xlab = "East (m)", ylab = "North (m)",
            region = TRUE,
            aspect = "iso",
            col.regions = gray((11:1)/11),
            data = ufc.plots)
contourplot(vol.m3.ha ~ east * north,
            main = expression(paste("Volume (", m^3, ha^{-1}, ")", sep = "")),
            xlab = "East (m)", ylab = "North (m)",
            region = TRUE,
            aspect = "iso",
            col.regions = red((11:1)/11),
            data = ufc.plots)

wireframe(vol.m3.ha ~ east * north,
          main = expression(paste("Volume (", m^3, ha^{-1}, ")", sep = "")),
          xlab = "East (m)", ylab = "North (m)",
          data = ufc.plots)
demo(persp)
demo(image)
demo(graphics)
savehistory(file = "history5")



###################
### new session ###
###################

#data.frame
X = data.frame(x1=c(1,3,75),x2=c(3,4,2))
X$x2
fix(X)
X$x3=c(2,3.54,6)
X$x5=X$x1/X$x3
X$lx1=log(X$x1)
X$x3=NULL
Y=subset(X,select=c(x1,x2));Y
Y=subset(X,select=-c(x1,x2));Y
Y=X[,"x1"];Y
Y=X[,c("x1","lx1")];Y

B=matrix(1:15, nco=5, byrow=TRUE)
dim(B); length(B)
fix(B)
nrow(B)
ncol(B)
class(B)
colnames(B)=c("Y","G","C","T","GDP");B
colnames(B)
rownames(B)=c("Tehran","Mashhad","Kermanshah");B

B=data.frame(B) # or  B=as.data.frame(B)
class(B)
B=as.matrix(B)
class(B)
(height <- runif(20, 150, 180))   #runif: random uniform
mass <- runif(20, 50, 90)
gender <- sample(c("M", "F"), 20, rep=TRUE)
?sample
color <- c("Blue", "Green","Brown")
eyes <- sample(color, 20, rep=TRUE)
table(gender); table(eyes)  #jadvale tozie faravani
table(gender, eyes)
H <- data.frame(height, mass, gender, eyes)
H; summary(H)
head(H)
tail(H)
str(H)
H[1,]
H$height
H$gender
is.data.frame(H)
is.matrix(H)
as.matrix(H) #H as a matrix
is.matrix(H)
H1=subset(H,mass<70)

#setclass and slots as altenative for list
mo=setClass("mohsen",slots=c(x="numeric",y="numeric"))
t=mo(x=1:10,y=1:10+rnorm(10)) # or t<- new("mohsen", x=1:10,y=1:10+rnorm(10))
t
class(t)
class(mo)
t@x=8:10     
t
slot(t,"x")=1:9
t

mo1<- setClass("mohsen1",slots = c(smooth = "numeric"), contains = "mohsen")
t1 <- mo1(t, smooth = 1:10)
t1
class(t1)
mo2=setClass("mohsen2",slots=c(m="lm"),contains = "mohsen1")
fit=lm(y~x, data=data.frame(x=1:10,y=1:10+rnorm(10)))
t2=mo2(t1,m=fit)
t2
t3=list(smooth=1:10,y=1:10+rnorm(10),m=fit)
t3$smooth
t3$sm
t2@sm# error
t2@smooth
t3[3]
t2[3] # error
objects(t3)
objects(t2) # error


#################Importing and exporting Data#################
getwd()
setwd("E:/book/Master Tehran/Term 2/Financial Econometrics/data for R")
# ,-separated (Original) csv Format or Comma Separated Values
data_1st <- read.csv("E:/book/Master Tehran/Term 2/Financial Econometrics/data for R/FinData1.csv")
data_2nd <- read.table("E:/book/Master Tehran/Term 2/Financial Econometrics/data for R/FinData1.csv", header=TRUE, sep = ",")
rm(data_2nd,data_1st)
# ;-separated csv
data2_1st <- read.table("E:/book/Master Tehran/Term 2/Financial Econometrics/data for R/FinData2.csv", sep=";",header=TRUE)
data2_2nd <- read.csv2("E:/book/Master Tehran/Term 2/Financial Econometrics/data for R/FinData2.csv")
rm(data2_1st,data2_2nd)

# Text File (Tab-separated)
data4 <- read.table("E:/book/Master Tehran/Term 2/Financial Econometrics/data for R/FinData3.txt", header=TRUE, sep="\t")


getwd()
write.csv(data_1st,file="data.csv") #The result of this line will be saved in the working directory.
write.table(data_1st,file="data.txt")

A <- seq(1, 10, by=0.5)
write.table(A, "A.txt")
sink("A2.txt")
A
summary(A)
sink()
summary(A)

install.packages("xlsx",dependencies = T)
library("xlsx")
data_3nd=read.xlsx("data.xlsx",sheetName = 1)
data(Journals,package = "AER") # loading data from package

# paste: Concatenate vectors after converting to character.
paste(1:12)
paste0(1:12)# same
as.character(1:12) # same
(nth <- paste(1:12, c("st", "nd", "rd", rep("th", 9))))
(nth <- paste0(1:12, c("st", "nd", "rd", rep("th", 9))))
(nth <- paste(1:2, c("st", "nd", "rd", rep("th", 9))))
(nth <- paste0(1:12, c("st", "nd", "rd", rep("th", 9)),sep=" "))

month.abb
letters
paste(month.abb, "is the", nth, "month of the year.")
paste(month.abb, letters)
paste(month.abb, "is the", nth, "month of the year.", sep = "_*_")
paste0(nth, collapse = ", ")
paste("1st", "2nd", "3rd", collapse = ", ")
paste("1st", "2nd", "3rd", sep = ", ")



