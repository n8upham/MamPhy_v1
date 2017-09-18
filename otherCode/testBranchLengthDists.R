# Code to test branch length distributions (sensu Venditti et al. 2010)
######


bl = rexp(100,1)

hist(bl)

dexp(bl,1,log=T)

sum(dexp(bl,1,log=T))

f = function(p,x) {-sum(dnorm(x,p[1],p[2],log=T)) }

nlm(f,p=c(1,1),bl)
