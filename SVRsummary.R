dat<-read.csv("TcrSvCalls.csv")
G<-dat[seq(1,34,2),4:9]+dat[seq(2,34,2),4:9]
ph<-dat[seq(1,34,2),2]
loc<-dat[seq(1,34,2),3]

## sv5 distinguishing G+S refugio from all others
## sv6 perfectly distinguishes refugio stripe from green

## sv3 inv translocation, only h154
 o<-lm(as.numeric(ph[1:13]=="m") ~G[1:13,3])
 o<-lm(as.numeric(ph[1:13]=="m") ~ as.numeric(G[1:13,3]==2))
 ## having two copies of this perfectly distinguishes melanic
 
## is sv2 associated with pattern
o <-lm(as.numeric(ph[1:13]=="g") ~ G[1:13,2])
summary(o)
#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)
#(Intercept)   0.2353     0.1645   1.431    0.180
#G[1:13, 2]    0.3235     0.2097   1.543    0.151
#
#Residual standard error: 0.4795 on 11 degrees of freedom
#Multiple R-squared:  0.1779,	Adjusted R-squared:  0.1032 
#F-statistic: 2.381 on 1 and 11 DF,  p-value: 0.1511

