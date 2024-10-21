#ver.2019.4.07
#modified 2020.7.3
# frequency=7 
#ver.2020.09.01
# new regression model
#ver.2021.2.12
# add TLS, SLS
#ver.2023.1.11
# m1 default change
# change x-axis for Z plot
# Ver93(2023.1.30)
# sid <- c(sid,(n-20):n)
# Ver1.0
# m1, sorder default change
# Ver1.1
# smooth=T -> ver93, default -> ver92

T <- TRUE
F <- FALSE
DDAIC <- 4
Frequency <- 4
TCRATE <- 0.7^(12/Frequency)
len <- length
#dyn.load("decomp6r.dll")
len <- length
parm <- function(m,n) {
	par(mfrow=c(m,n))
}
tsplot2 <- function(x,y=NULL,...)
{
	if(is.null(y)) {
	plot.ts(x,...)
}
else {
	plot.ts(cbind(x,y),plot.type="single",...,col=1:4)
}
}
"regsiml"<-
function(data, reg = NULL, m1 = 2, period = 4, sorder = 1, log = 0,pb=2,pa=2,mtype=1,iplot=T,smooth=F)
{
    if(m1 < 3) m1 <- max(3,as.integer(length(data)/10))
    h <- 0
    pb <- pb*(sorder > 0)
    pa <- pa*(sorder > 0)    
    odata <- data
    if(log > 0) data <- log(data)
        n0 <- len(data)
    if(sorder > 0) {
        dx <- diff(data)
        
    dx <- c(0,rep(dx[1:period],pb),dx,rep(dx[(n0-period):(n0-1)],pa))
    x.tmp <- cumsum(dx)
    x.tmp <- x.tmp + (data[1]-x.tmp[period*pb+1])
    data <- x.tmp
    }
    n <- len(data)
    if(!is.null(reg)) {
        reg <- cbind(reg)
        k <- ncol(reg)
    }
    dx <- diff(data)

        
    mat <- siml.mat(n,m1,type=mtype)
    n1 <- dim(mat$K)[1]
    z.y <- mat$K %*% c(data)
#    print(dim(mat$K
#    print(length(data))
    z.y.o <- z.y
    sid <- integer(0)
    if(sorder > 0) {
        if(period > 2) {
        sid <- ceiling(2*(n1+h)/period)+(-sorder):(sorder)

if(F) {        
        if(period==4) { sid <- c(sid,n1-(( max(sorder+1-h,0) ):0))}
        if(period==12) {
            sid <- c(sid,ceiling(4*(n1+h)/period)+(-sorder):(sorder))
            sid <- c(sid,ceiling(6*(n1+h)/period)+(-sorder):(sorder))
            sid <- c(sid,ceiling(8*(n1+h)/period)+(-sorder):(sorder))
            sid <- c(sid,ceiling(10*(n1+h)/period)+(-sorder):(sorder))                                    
            sid <- c(sid,n1-(( max(sorder+1-h,0) ):0))
        }  
        if(period==7) {
            sid <- c(sid,ceiling(4*(n1+h)/period)+(-sorder):(sorder))
            sid <- c(sid,ceiling(6*(n1+h)/period)+(-sorder):(sorder))
            sid <- c(sid,n1-(( max(sorder+1-h,0) ):0))
        }  
 } else {
    if(period > 4) {
    for(i in 2:(ceiling(period/2)-1)) {
            sid <- c(sid,ceiling(2*i*(n1+h)/period)+(-sorder):(sorder))
             }
          }
                   
        }
        }
        if(!smooth) {
           if(period %% 2 == 0)  sid <- c(sid,n1-(( max(sorder+1-h,0) ):0)) #Ver.92
        } else {
           ## Ver.93 start
           if(period %% 2 == 0)  sid <- c(sid,n1-(( as.integer((n1-rev(sid)[1])/2) ):0))
            ## Ver.93 end
            }
        z.s <- z.y
        z.s[-sid] <- 0

        z.y[sid] <- 0
    }
    
    else { z.s <- rep(0,n1)
          sid <- -(1:n)
          }
        
    if(!is.null(reg)) {

        reg1 <- NULL
        if(pb > 0) {for(ii in seq(pb)) {
            reg1 <- rbind(reg1,reg[(1:period)*(ii > 0),,drop=F])
                    }
            }
        reg1 <- rbind(reg1,reg)
        if(pa > 0) {for(ii in seq(pa)) {
            reg1 <- rbind(reg1,reg[((n0-period+1):n0)*(ii > 0),,drop=F])
                    }
            }
        reg.org <- reg
        reg <- reg1
        
        z.d <- mat$K %*% reg
        
        zz <- lsfit(z.d[1:m1,],z.y[1:m1],inter=F)
      #        zz$coef <- zz$coef[-1]
        vvv2 <- mean((zz$res)^2)
 #       vvv3 <- vvv2
        res <- c(zz$res,rep(0,n1-m1))
        print(zz$coef)

        if(sorder > -1) { ## modi at 2020.12.16  
            regname <- dimnames(reg)[[2]]
            if(!is.null(regname)) {
                ireg <- substring(paste(regname,"XX",sep=""),1,2) == "ao"
                ireg <- ireg |(regname == "hol")|(regname == "sls") | (regname=="vat") |(regname=="ly")  | (regname=="vatt")
            }
#            if(sum(!ireg) > 0) {
#                z.d[(m1:n1),!ireg] <- 0
#           }
#            zddd <- c(sum(abs(z.s[sid]))*1 , sum(abs(z.s[sid]-c(z.d %*% zz$coef)[sid])))
#               if(iplot) print(zddd)
               zz.tmp <- lsfit(z.d[-sid,],z.y[-sid],inter=F)
         #      zz.tmp$coef <- zz.tmp$coef[-1]
               vvv3 <- mean((zz.tmp$res)^2)                
    print(zz.tmp$coef)
        #    if(diff(zddd) > 0 ||any(dimnames(reg)[[2]]=="ao") ||any(dimnames(reg)[[2]]=="hol") ) {
                tmp1 <- ls.print(zz,print=iplot)[[1]]
             #   print(tmp1)
                 #zz <- lsfit(z.d,z.y.o,inter=F)  
                tmp2 <- ls.print(zz.tmp,print=iplot)[[1]]
              #  print(tmp2)
         #  if(as.numeric(tmp2[3]) > as.numeric(tmp1[3]) & any(ireg)) {
           if(as.numeric(tmp2[3]) > as.numeric(tmp1[3]) | any(ireg)) { # mody at 21.2.13             
        #       zz$coef[ireg] <- zz.tmp$coef[ireg]
         #      zz$res <- zz.tmp$res
               zz <- zz.tmp
               res <- z.y
            res[-sid] <- zz$res
            # res <- zz$res
                    zz$res <- res
                    
                }
       #     }
            z.dd <- z.d
             tmtm <- substring(regname,1,3)=="reg" 
               
            if(sum(tmtm > 0)) z.dd[,tmtm] <- 0
            z.s[sid] <- z.s[sid]-c(z.dd %*% zz$coef)[sid]
       }     
        cz.d <- c(z.d %*% zz$coef)
        res <- c(zz$res[1:m1],rep(0,n1-m1))
        trade <- t( t(reg) * zz$coef)
#        para <- log(sum((res)^2)/m1)*m1+2*(k+1)
        para <- log(c(vvv2,vvv3))*n1+2*(k+1)
        
        #        res <- res + c(z.d %*% zz$coef)
        trend <- mat$inv %*% res + mat$TT %*%c(data-c(reg %*% zz$coef))
        #cumsum(c((data-c(reg %*% zz$coef))[1],mat %*% res))
#        trend <- cumsum(c(data[1],mat %*% res))

        #        browser()
#             z.s.o <- z.s                           #        print(zz$coef)

#browser()
    }
    else {
        cz.d <- rep(NA,n1)
        if(length(sid)==0) vvv3 <- mean(z.y^2)
        else vvv3 <- mean(z.y[-sid]^2)
        if(n1 > m1)  z.y[(m1+1):n1] <- 0
        trade <- cbind(rep(0,n))
#        para <- log(sum((z.y)^2)/m1)*m1+2        
        para <- log(c(sum((z.y)^2)/m1,vvv3))*n1+2        

        trend <- mat$inv %*% z.y + mat$TT%*%c(data)
        #cumsum(c(data[1],mat %*% z.y))
    }
    seasonal <- mat$inv %*% z.s
    #cumsum(c(0,c(mat %*% z.s)))
#browser()
    if(pb > 0 | pa > 0) {
        iddd4 <- pb*period + seq(n0)
        seasonal <- seasonal[iddd4]
        trend <- trend[iddd4]
        trade <- trade[iddd4,,drop=F]
        data <- data[iddd4]
        }

    z <- list(data=odata,trend=trend,seasonal=seasonal,
              noise=data-trend-apply(trade,1,sum)-seasonal,trade=trade, para=para, h=h, z=z.y.o, sid=sid, cz.d = cz.d)

	return(z)
}
"x12siml"<-
function(data, reg = NULL, trend = NULL, ilog = 0, frequency = 4, start = 
	c(1994, 1), iplot = T, sorder = NULL, mtype=1, pb=2,pa=2, smooth=F, ...)
{
	n <- length(data)
     
        if(is.null(trend)) trend <- as.integer(n/max(6,frequency))+1
        if(is.null(sorder)) { sorder <- max(1,as.integer(n/frequency/5)+1)
            if(sorder > 2) sorder <- 2+as.integer((sorder-2)/3)
            if(sorder > 5) sorder <- 6
        }
        
#ls=c(2003,1)
#tc=c(2007,4)
#ao=c(2003,2)
#rp=c(2008,3,2009,1)
#vat=c(1997,1)
    ar <- 0
    reg2 <- list(...)

    nn <- length(reg2)
    dimnames.reg <- NULL
    if(!is.null(reg)) {reg <- cbind(reg)
        dimnames.reg <- c(dimnames(reg)[[2]],rep("",dim(reg)[2]))[seq(dim(reg)[2])]
        dimnames.reg <- paste("reg",seq(dim(reg)[2]),dimnames.reg,sep="")
                       
        }
                 
	if(nn > 0) {
    if(any(names(reg2)=="osls")) {
      ey <- ceiling(n/frequency)+start[1]
      isls <- seq(nn)[names(reg2)=="osls"]
      for(i in isls) {
        sts <- reg2[[i]]
        sts <- c(rbind((sts[1]):ey,sts[2]))
        reg2[[i]] <- sts
        names(reg2)[i]<- "hol"
      }
    }
		for(i in seq(nn)) {
                    switch(names(reg2[i]),
                           "ls"={
				z <- rep(0, n)
				tt <- sum((reg2[[i]] - start) * c(frequency, 1)
				  ) + 1
				z[1:(tt - 1)] <- -1
                                z <- z+1
				reg <- cbind(reg, z)
			},
			"tc"= {
                            tcrate <- TCRATE #0.69999999999999996^(12/frequency)
				z <- rep(0, n)
				tt <- sum((reg2[[i]] - start) * c(frequency, 1)
				  ) + 1
				z[tt:n] <- tcrate^(0:(n - tt))
				reg <- cbind(reg, z)
			},
			"ao"= {
				z <- rep(0, n)
				tt <- sum((reg2[[i]] - start) * c(frequency, 1)
				  ) + 1
				z[tt] <- 1
				reg <- cbind(reg, z)
			},
			"aot"= {
				z <- rep(0, n)
				tt <- sum((reg2[[i]] - start) * c(frequency, 1)
				  ) + 1
				z[tt] <- 1
				reg <- cbind(reg, z)
			},                           
            "aao"= {
				z <- rep(0, n)
				tt <- sum((reg2[[i]] - start) * c(frequency, 1)
				  ) + 1
				z[tt] <- 1
				reg <- cbind(reg, z)
			},
			"hol"= {
				z <- rep(0, n)
                nnn2 <- length(reg2[[i]])/2
                for(j in seq(nnn2)) {
				tt <- sum((reg2[[i]][j*2-(1:0)] - start) * c(frequency, 1)
				  ) + 1
        tt <- tt[tt <= n]
				z[tt] <- 1
                    }
				reg <- cbind(reg, z)
			},
			"sls"= {
				z <- rep(0, n)
        nyy <- floor(n/frequency)+1
        yyy <- c(rep(start[1],frequency-start[2]+1),t(matrix(rep((start[1]+1):(start[1]+nyy),frequency),nrow=nyy)))[1:n]
        mmm <- c((start[2]):frequency,rep(1:frequency,nyy))[1:n]
        sty <- reg2[[i]][1]
        stm <- reg2[[i]][2]        
        z[yyy >= sty] <- -1/(frequency-1)
        z[yyy >= sty & mmm==stm] <- 1       
				reg <- cbind(reg, z)
			},
     "ly"= {
                     z <- rep(0, n)
                     if(length(reg2[[i]])==1) 
				tt <- reg2[[i]]+(0:n)*(frequency*4)
                            else
				tt <- sum((reg2[[i]] - start) * c(frequency, 1)

                         ) + 1 + (0:n)*(frequency*4)                                
				z[tt] <- 1
                                z <- z[1:n]
				reg <- cbind(reg, z)
			},

			"vat"={
				z <- rep(0, n)
				tt <- sum((reg2[[i]] - start) * c(frequency, 1)
				  ) + 1
				z[tt] <- 1
				z[tt + 1] <- -1
				reg <- cbind(reg, z)
			},
			"vatt"={
				z <- rep(0, n)
				tt <- sum((reg2[[i]] - start) * c(frequency, 1)
				  ) + 1
				z[tt] <- 1
				z[tt + 1] <- -1
				reg <- cbind(reg, z)
			},
             "rp"= {
				z <- rep(0, n)
				tt1 <- sum((reg2[[i]][1:2] - start) * c(
				  frequency, 1)) + 1
				tt2 <- sum((reg2[[i]][3:4] - start) * c(
				  frequency, 1)) + 1
				z[1:(tt1)] <- -1
				z[(tt1 + 1):(tt2 - 1)] <- (((tt1 + 1):(tt2 - 1)
                                ) - tt1)/(tt2 - tt1) - 1
                                z <- z+1
				reg <- cbind(reg, z)
			},
            "tls"= {
				z <- rep(0, n)
				tt1 <- sum((reg2[[i]][1:2] - start) * c(
				  frequency, 1)) + 1
				tt2 <- sum((reg2[[i]][3:4] - start) * c(
				  frequency, 1)) + 1
				z[(tt1):tt2] <- 1
				reg <- cbind(reg, z)
			},
			"tcrp"= {
				z <- rep(0, n)
				tt1 <- sum((reg2[[i]][1:2] - start) * c(
				  frequency, 1)) + 1
				tt2 <- sum((reg2[[i]][3:4] - start) * c(
				  frequency, 1)) + 1
				tt3 <- sum((reg2[[i]][5:6] - start) * c(
				  frequency, 1)) + 1
				tcrate <- TCRATE
#				z[1:(tt1)] <- -1
				z[(tt1 + 1):(tt2 - 1)] <- (((tt1 + 1):(tt2 - 1)
				  ) - tt1)/(tt2 - tt1)
				z[(tt2):tt3] <- tcrate^(0:(tt3 - tt2))
				z[tt3:n] <- z[tt3]
				reg <- cbind(reg, z)
			},
			"tcrp1"= {
				z <- rep(0, n)
				tt1 <- sum((reg2[[i]][1:2] - start) * c(
				  frequency, 1)) + 1
				tt2 <- sum((reg2[[i]][3:4] - start) * c(
				  frequency, 1)) + 1
				tt3 <- sum((reg2[[i]][5:6] - start) * c(
				  frequency, 1)) + 1
				tcrate <- TCRATE #0.80000000000000004	#^(12/frequency)
#				z[1:(tt1)] <- -1
				z[(tt1 + 1):(tt2 - 1)] <- (((tt1 + 1):(tt2 - 1)
				  ) - tt1)/(tt2 - tt1)
				z[tt2] <- 1
				z[(tt2 + 1):tt3] <- tcrate^(0:(tt3 - tt2 - 1))
				z[tt3:n] <- z[tt3]
				reg <- cbind(reg, z)
			},
                        {
                            nn <- nn-1
                            cat(paste("Warning:",names(reg2[i]), "is not supported\n"))
                            reg2[i] <- NULL
                        }
                        )
                  }
        }
            if(!is.null(reg)) {    dimnames(reg) <- list(NULL,
                                                       c(dimnames.reg,names(reg2)))
        nn0 <- dim(reg)[2]
                               }
#	if(!is.null(reg)) {
#		if(ncol(reg) > 6)
#			reg <- reg[, 1:6]
#	}
    
	zz <- regsiml(data, reg, m1 = trend, log = ilog, period = 
                                                             frequency, sorder = sorder,
                      mtype=mtype, pb=pb, pa=pa, iplot=iplot, smooth=smooth)
    		z.trend <- zz$trend


                                        # x11()
		if(!is.null(reg) | ar > 0)
			par(mfrow = c(3, 2))
		else par(mfrow = c(3, 2))
		if(ilog > 0)
			data <- log(data)
		z.trend <- zz$trend
		if(!is.null(reg)) {
                    zz.dumm <- seq(ncol(zz$trade))
                    zz.dumm <- zz.dumm[dimnames(zz$trade)[[2]] != "ao" &
                                       dimnames(zz$trade)[[2]] != "vat" &
                                       dimnames(zz$trade)[[2]] != "hol" &
                                       dimnames(zz$trade)[[2]] != "sls" &                                       
                                       dimnames(zz$trade)[[2]] != "ly" 
                                       ]

                    if(length(zz.dumm) >0) {                    
                        z.trend <- z.trend + c(apply(zz$trade[,zz.dumm,drop=F],1,sum)) }
                    }
    
                    zz.seasonal <- zz$seasonal
                    zz.seasonal.o <- zz.seasonal
    if(any(dimnames(reg)[[2]] == "ly"))     
                    zz.seasonal <- zz.seasonal+zz$trade[,"ly"]
    if(any(dimnames(reg)[[2]] == "hol")) {
                    ididi <- seq(nn0)[dimnames(reg)[[2]] == "hol"]
                    for(ijij in ididi) {
                    zz.seasonal <- zz.seasonal+zz$trade[,ijij]
                        }
                    }
    if(any(dimnames(reg)[[2]] == "sls")) {
                    ididii <- seq(nn0)[dimnames(reg)[[2]] == "sls"]
                    for(ijij in ididii) {
                    zz.seasonal <- zz.seasonal+zz$trade[,ijij]
                        }
                    }
    
    		z.adj <- data - zz.seasonal


	if(iplot) {    
   #                                 if(any(dimnames(reg)[[2]] == "sls")) {
   #                                 print(zz$trade[,ijij]) }
            tsplot2(ts(data, start = start, frequency = frequency), ts(
			z.trend, start = start, frequency = frequency), sub = 
                         "org&trend", main = paste("SIML",trend)) #,"h=",zz$h))
            
            tsplot2(ts(zz.seasonal, start = start, frequency = frequency), 
			sub = "seasonal",col=2)
            lines(ts(zz.seasonal.o,start=start,frequency=frequency))
                    abline(h = 0, lty = 2)

                    zz.noise <- data-(zz.seasonal+z.trend)

                    tsplot2(ts(zz.noise, start = start, frequency = frequency), sub
			 = "noise")
		abline(h = 0, lty = 2)

            if(!is.null(reg)) {
                
                if(nn0 > nn) {nn2 <- nn0-nn
                             for(jj in seq(nn2)) reg2 <- c(list(""),reg2)
                             }
            if(any(dimnames(zz$trade)[[2]] == "ao") ||
                                       any(dimnames(zz$trade)[[2]] == "vat")) {
                aoline <- apply(zz$trade[,dimnames(zz$trade)[[2]] == "ao" | dimnames(zz$trade)[[2]] == "vat",drop=F],1,sum)
                lines(ts(aoline,start=start,frequency=frequency),type="h",col=2)
                }
			tsplot2(ts(c(apply(zz$trade,1,sum)), start = start, frequency = 
				frequency), sub = "reg")
			abline(h = 0, lty = 2)
			tit <- NULL
      reg2t <- reg2          
			for(i in 1:nn0) {
        if(length(reg2t[[i]]) > 8) reg2t[[i]] <- c(reg2t[[i]][1:2],"-","-",rev(reg2t[[i]])[2:1])
				tit <- paste(tit, " ", dimnames(reg)[[2]][i], ":(", 
				  paste(reg2t[[i]], collapse = " "), ")", sep = 
				  "")
				if(i == 3)
				  tit <- paste(tit, "\n", sep = "")
			}
			title(tit, cex = 0.34999999999999998)
		}
            		if(ilog > 0) {
			data <- exp(data)
			z.adj <- exp(z.adj)
                        }
            aiccT <- zz$para[1]	
            aiccA <- zz$para[2]
	#	  if(!is.null(reg)) aicc <-aicc+2*ncol(reg)
		tsplot2(ts(data, start = start, frequency = frequency), ts(z.adj,
			start = start, frequency = frequency), sub = "org&adj", 
			main = paste("AIC=", aiccT,"(",aiccA,")"))
		if(as.integer(ar > 0) + as.integer(!is.null(reg)) == 10) {
			tsplot2(1:100, type = "n", xaxt = "n", yaxt = "n", bty
				 = "n")
			title(paste(substitute(data)), adj = 0)
		}
                if(!is.null(zz$z)) {nnz <- length(zz$z)
                                   plot(seq(nnz)/(2*nnz),zz$z,type="n",ylab="",xlab="")
                                   par(new=T)
                                   tsplot2(zz$z,main="Z",sub="freq",xaxt="n",yaxt="n",ylab="",xlab="")
                                    if(length(zz$sid) > 0) {
                                   abline(v=seq(length(zz$z))[zz$sid],col=5,lty=2)
                                         }
                                          abline(h=0,col=2,lty=1)
                                     lines(zz$z)
                                    lines(zz$cz.d,col=3)
                                   }
        } else {
            		if(ilog > 0) {
			#data <- exp(data)
			#z.adj <- exp(z.adj)
                        }
                    aiccT <- zz$para[1]	
            aiccA <- zz$para[2]

                    

        }

    
#    print(reg2)
    zz <-  list(data=data,trend=z.trend,seasonal = zz.seasonal,
                noise= zz$noise, adj=z.adj, dummy=zz$trade, aic=c(aiccT,aiccA),
                m1=trend, sorder=sorder,Z=zz$z,para=zz$para,hh=zz$h)
	return(zz)
}
"tkramp"<-
function(data, tt, reg = NULL, trend = 2, ilog = 0, frequency = 4, 
	start = c(1994, 1), iplot = T, sorder = 1,mtype=1, ...)
{
#ls=c(2003,1)
#tc=c(2007,4)
#ao=c(2003,2)
#rp=c(2008,3,2009,1)
#vat=c(1997,1)
	rpt <- tt
	n <- length(data)
	k <- rpt[3] * frequency + rpt[4] - (rpt[1] * frequency + rpt[2]) + 1
	tt <- cbind(rep(rpt[1], k), rpt[2] - 1 + seq(k))
	while(any(tt[, 2] > frequency)) {
		tt[tt[, 2] > frequency, 1] <- tt[tt[, 2] > frequency, 1] + 1
		tt[tt[, 2] > frequency, 2] <- tt[tt[, 2] > frequency, 2] - 
			frequency
	}
	z <- x12siml(data = data, reg = reg, trend = trend, ilog = 
		ilog, frequency = frequency, start = start, iplot = T, sorder
                = sorder,mtype=mtype, ...)
        hh <- z$h
	maic <- z$para[1] - DDAIC -min(as.integer(k/10),4)
	rpp <- NULL
	for(i in 1:(k - 2)) {
		for(j in (i + 2):k) {
			for(m1 in (i + 1):(j - 1)) {
				for(m2 in (m1):(min(m1 + 1, j - 1))) {
				  rp1 <- c(tt[i,  ], tt[m1,  ])
				  rp2 <- c(tt[m2,  ], tt[j,  ])
				  tmp <- x12siml(data = data, reg = reg, 
				    trend = trend, ilog = ilog,
				    frequency = frequency, start = start, iplot
				     = F, sorder = sorder,mtype=mtype, rp = rp1, rp = rp2, 
				    ...)
				  if(tmp$para[1] < maic) {
				    z <- tmp
				    maic <- tmp$para[1]
				    rpp <- c(rp1, rp2)
				    print(c(paste(maic), paste(rpp)))
				  }
				}
			}
		}
	}
	if(!is.null(rpp)) {
		tmp <- x12siml(data = data, reg = reg, trend = trend,
			ilog = ilog, frequency = frequency, start = start, 
			iplot = T, sorder = sorder,mtype=mtype, rp = rpp[1:4], rp = rpp[5:8
			], ...)
	}
	z$rpp <- rpp
	return(z)
}
"outlier"<-
function(data, tt, type = "ls", reg = NULL, trend = 2, ilog = 0, 
	frequency = 4, start = c(1994, 1), iplot = T, sorder = 1,mtype=1, ...)
{
#ls=c(2003,1)
#tc=c(2007,4)
#ao=c(2003,2)
#rp=c(2008,3,2009,1)
#vat=c(1997,1)
	n <- length(data)
	k <- tt[3] * frequency + tt[4] - (tt[1] * frequency + tt[2]) + 1
	tt <- cbind(rep(tt[1], k), tt[2] - 1 + seq(k))
	while(any(tt[, 2] > frequency)) {
		tt[tt[, 2] > frequency, 1] <- tt[tt[, 2] > frequency, 1] + 1
		tt[tt[, 2] > frequency, 2] <- tt[tt[, 2] > frequency, 2] - 
			frequency
	}
	z <- x12siml(data = data, reg = reg, trend = trend , ilog = 
		ilog, frequency = frequency, start = start, iplot = T, sorder
                = sorder,mtype=mtype, ...)
        hh <- z$h
	maic <- z$para - DDAIC-min(as.integer(k/10),4)
	rpp <- NULL
    ia <- 1
	for(i in 1:(k)) {
		switch(type,
			ao = {
				tmp <- x12siml(data = data, reg = reg, trend
				   = trend,  ilog = ilog, frequency = 
				  frequency, start = start, iplot = F, sorder
				   = sorder,mtype=mtype, ao = tt[i,  ], ...)
                ia <- 2
			}
			,
			hol = {
				tmp <- x12siml(data = data, reg = reg, trend
				   = trend,  ilog = ilog, frequency = 
				  frequency, start = start, iplot = F, sorder
				   = sorder,mtype=mtype, hol = tt[i,  ], ...)
                ia <- 2
			}
			,               
			sls = {
				tmp <- x12siml(data = data, reg = reg, trend
				   = trend,  ilog = ilog, frequency = 
				  frequency, start = start, iplot = F, sorder
				   = sorder,mtype=mtype, sls = tt[i,  ], ...)
                ia <- 2
			}
			,               
			ly = {
				tmp <- x12siml(data = data, reg = reg, trend
				   = trend,  ilog = ilog, frequency = 
				  frequency, start = start, iplot = F, sorder
				   = sorder,mtype=mtype, ly = tt[i,  ], ...)
                ia <- 2
			}
			,

                       ls = {
				tmp <- x12siml(data = data, reg = reg, trend
				   = trend,  ilog = ilog, frequency = 
				  frequency, start = start, iplot = F, sorder
				   = sorder,mtype=mtype, ls = tt[i,  ], ...)
			}
			,
			tc = {
				tmp <- x12siml(data = data, reg = reg, trend
				   = trend,  ilog = ilog, frequency = 
				  frequency, start = start, iplot = F, sorder
				   = sorder,mtype=mtype, tc = tt[i,  ], ...)
			}
			,
			rp = {
				if(i < k) {
				  for(j in (i + 1):(min(i + 10, k))) {
				    rp1 <- c(tt[i,  ], tt[j,  ])
				    tmp <- x12siml(data = data, reg = reg, 
				      trend = trend,  ilog = ilog, 
				      frequency = frequency, start = start, 
				      iplot = F, sorder = sorder,mtype=mtype, rp = rp1, ...
				      )
				    if(tmp$para[1] < maic[1]) {
				      z <- tmp
				      maic <- tmp$para
				      rpp <- c(rp1)
				      print(c(paste(maic), paste(rpp)))
				    }
				  }
				}
			}
      ,
			tls = {
				if(i < k) {
				  for(j in (i + 1):(min(i + 10, k))) {
				    rp1 <- c(tt[i,  ], tt[j,  ])
				    tmp <- x12siml(data = data, reg = reg, 
				      trend = trend,  ilog = ilog, 
				      frequency = frequency, start = start, 
				      iplot = F, sorder = sorder,mtype=mtype, tls = rp1, ...
				      )
				    if(tmp$para[1] < maic[1]) {
				      z <- tmp
				      maic <- tmp$para
				      rpp <- c(rp1)
				      print(c(paste(maic), paste(rpp)))
				    }
				  }
				}
			}
			)
		if(tmp$para[ia] < maic[ia]) {
			z <- tmp
			maic <- tmp$para
			rpp <- c(tt[i,  ])
			print(c(paste(maic), paste(rpp)))
		}
	}
	if(!is.null(rpp)) {
		switch(type,
			ao = {
				tmp <- x12siml(data = data, reg = reg, trend
				   = trend,  ilog = ilog, frequency = 
				  frequency, start = start, iplot = T, sorder
				   = sorder,mtype=mtype, ao = rpp, ...)
			}
			,
			hol = {
				tmp <- x12siml(data = data, reg = reg, trend
				   = trend,  ilog = ilog, frequency = 
				  frequency, start = start, iplot = T, sorder
				   = sorder,mtype=mtype, hol = rpp, ...)
			}
			,
      sls = {
				tmp <- x12siml(data = data, reg = reg, trend
				   = trend,  ilog = ilog, frequency = 
				  frequency, start = start, iplot = T, sorder
				   = sorder,mtype=mtype, sls = rpp, ...)
			}
			,
			ly = {
				tmp <- x12siml(data = data, reg = reg, trend
				   = trend,  ilog = ilog, frequency = 
				  frequency, start = start, iplot = T, sorder
				   = sorder,mtype=mtype, ly = rpp, ...)
			}
			,

                       ls = {
				tmp <- x12siml(data = data, reg = reg, trend
				   = trend,  ilog = ilog, frequency = 
				  frequency, start = start, iplot = T, sorder
				   = sorder,mtype=mtype, ls = rpp, ...)
			}
			,
			tc = {
				tmp <- x12siml(data = data, reg = reg, trend
				   = trend,  ilog = ilog, frequency = 
				  frequency, start = start, iplot = T, sorder
				   = sorder,mtype=mtype, tc = rpp, ...)
			}
			,
			rp = {
				tmp <- x12siml(data = data, reg = reg, trend
				   = trend,  ilog = ilog, frequency = 
				  frequency, start = start, iplot = T, sorder
				   = sorder,mtype=mtype, rp = rpp, ...)
			}
,
			tls = {
				tmp <- x12siml(data = data, reg = reg, trend
				   = trend,  ilog = ilog, frequency = 
				  frequency, start = start, iplot = T, sorder
				   = sorder,mtype=mtype, tls = rpp, ...)
			}
			)
	}
	z[[type]] <- rpp
	return(z)
}
"tcramp"<-
function(data, tt, reg = NULL, trend = 2,  ilog = 0, frequency = 4, 
	start = c(1994, 1), iplot = T, sorder = 1,mtype=1,...)
{
#ls=c(2003,1)
#tc=c(2007,4)
#ao=c(2003,2)
#rp=c(2008,3,2009,1)
#vat=c(1997,1)
	n <- length(data)
	k <- tt[3] * frequency + tt[4] - (tt[1] * frequency + tt[2]) + 1
	tt <- cbind(rep(tt[1], k), tt[2] - 1 + seq(k))
	while(any(tt[, 2] > frequency)) {
		tt[tt[, 2] > frequency, 1] <- tt[tt[, 2] > frequency, 1] + 1
		tt[tt[, 2] > frequency, 2] <- tt[tt[, 2] > frequency, 2] - 
			frequency
	}
	z <- x12siml(data = data, reg = reg, trend = trend,  ilog = 
		ilog, frequency = frequency, start = start, iplot = T, sorder
		 = sorder,mtype=mtype, ...)
    maic <- z$para[1] - DDAIC-min(as.integer(k/10),4)
    hh <- z$h
	rpp <- NULL
	for(i in 1:(k - 2)) {
		for(j in (i + 2):k) {
			for(m1 in (i + 1):(j - 1)) {
				for(m2 in (m1):(min(m1 + 1, j - 1))) {
				  rp1 <- c(tt[i,  ], tt[m1,  ])
				  rp2 <- c(tt[m2,  ], tt[j,  ])
				  rp <- c(rp1, rp2[3:4])
				  if(m1 == m2)
				    tmp <- x12siml(data = data, reg = reg, 
				      trend = trend,  ilog = ilog, 
				      frequency = frequency, start = start, 
				      iplot = F, sorder = sorder,mtype=mtype, tcrp = rp, 
				      ...)
				  else tmp <- x12siml(data = data, reg = reg, 
				      trend = trend,  ilog = ilog, 
				      frequency = frequency, start = start, 
				      iplot = F, sorder = sorder,mtype=mtype, tcrp1 = rp, 
				      ...)
				  if(tmp$para[1] < maic) {
				    z <- tmp
				    maic <- tmp$para[1]
				    rpp <- c(rp1, rp2)
				    print(c(paste(maic), paste(rpp)))
				  }
				}
			}
		}
	}
	if(!is.null(rpp)) {
		if(rpp[4] == rpp[6])
			tmp <- x12siml(data = data, reg = reg, trend = trend, 
				 ilog = ilog, frequency = frequency, 
				start = start, iplot = T, sorder = sorder,mtype=mtype, tcrp
				 = c(rpp[1:4], rpp[7:8]), ...)
		else tmp <- x12siml(data = data, reg = reg, trend = trend, ar
				 = ar, ilog = ilog, frequency = frequency, 
				start = start, iplot = T, sorder = sorder, mtype=mtype,
				tcrp1 = c(rpp[1:4], rpp[7:8]), ...)
	}
	z$rpp <- rpp
	return(z)
}
#input <- function() {
#	system("notepad tmptmp123.dat")
#	z <- scan("tmptmp123.dat")
#	system("del tmptmp123.dat")
#	return(z)
                                        #}

## n-1, y0=y0
makeMat20 <- function(n,m=NULL,h=0,type=2,trans=1,allz=1,init=F,inv=1)
# type=1 : T1m
# type=2 : T2m
# type=3 : Muller=Watson (Econometrica, 2018)
# Müller, U. K., & Watson, M. W. (2018). Long‐Run Covariability. Econometrica, 86(3), 775-804.
# type=13 : type1-type3
{
    if(is.null(m)) m <- as.integer(n/10)+3
    if(m==1) m <- 2
    J <- matrix(0,ncol=n,nrow=n)
    J1 <- J
    Jn <- J
    if(m!=n |!trans|!allz) {
    J1[,1] <- 1
    Jn[,n] <- 1
    }
    R1 <- cbind(0,diag(rep(1,n-1)))
    R1t <- t(R1)
    Rn1 <- cbind(diag(rep(1,n-1)),0)
    Rn1t <- t(Rn1)

    I1n <- array(1/n,dim=c(n,n))
    In <- diag(rep(1,n))
    Ci <- diag(rep(1,n-1))
    diag(Ci[-1,-(n-1)]) <- -1
    Cti <- t(Ci)
    C <- solve(Ci)
    Ct <- t(C)
    if(type <= 2) {
        P <- sqrt(2/(n-1+h+0.5))*
           cos((2*pi/(2*(n-1+h)+1))*
               outer(seq(n-1)-0.5,seq(n-1)-0.5)
               )
               Pr <- sqrt(2/(n-1+h+0.5))*
           sin((2*pi/(2*(n-1+h)+1))*
               outer(seq(n-1)-0.5,seq(n-1))
               )

        Pm <- P[1:m,]
        Prm <- Pr[1:m,]

    Pmt <- t(Pm)
    Prmt <- t(Prm)

        P1 <- R1t %*% C %*% Pmt %*% Pm%*% Ci %*% R1
    P2 <- Rn1t %*% Ct %*% Prmt %*% Prm %*% Cti %*% Rn1
#    T1 <- P1 + (In-P1)%*%J1
        T2 <- In
        TT1 <- NULL
        TT2 <- NULL
        if(init) {
            if(type==1) return(P1 + (In-P1)%*%J1)
            if(type==2) return(P2 + (In-P2)%*%Jn)
            }
        for(i in 1:3) {
        T1 <- P1 + (In-P1)%*%J1%*%T2            
        T2 <- P2 + (In-P2)%*%Jn%*%T1

        }
        if(trans) {
            if(allz) {
        P1 <-  P %*% Ci %*% R1
        P2 <-  Pr %*% Cti %*% Rn1            
            } else {
        P1 <-  Pm%*% Ci %*% R1
        P2 <-  Prm%*% Cti %*% Rn1
                }
        TT2 <- P2 %*%(In-Jn%*%T1)
        TT1 <- P1 %*%(In-J1%*%T2)
        
 
        }
        if(type==2) {
                     if(inv) iTmat <- Rn1t %*% Ct %*% t(Pr)
                       Tmat <- list(TT=Jn%*%T1,G=TT2,inv=iTmat)
         }     else  {
                     if(inv) iTmat <- R1t %*% C %*% t(P)
                       Tmat <- list(TT=J1%*%T2,G=TT1,inv=iTmat)
}
        }
               
         Tmat[["order"]] <- c(m,n,type)
    return(Tmat)
}
## n, y0=y1
makeMat21 <- function(n,m=NULL,h=0,type=2,trans=1,allz=1,init=F,inv=1)
# type=1 : T1m
# type=2 : T2m
# type=3 : Muller=Watson (Econometrica, 2018)
# Müller, U. K., & Watson, M. W. (2018). Long‐Run Covariability. Econometrica, 86(3), 775-804.
# type=13 : type1-type3
{
    if(is.null(m)) m <- as.integer(n/10)+3
    if(m==1) m <- 2
    J <- matrix(0,ncol=n,nrow=n)
    J1 <- J
    Jn <- J
    if(m!=n |!trans|!allz) {
    J1[,1] <- 1
    Jn[,n] <- 1
    }
    
    I1n <- array(1/n,dim=c(n,n))
    In <- diag(rep(1,n))
    Ci <- diag(rep(1,n))
    diag(Ci[-1,-n]) <- -1
    Cti <- t(Ci)
    C <- solve(Ci)
    Ct <- t(C)
    if(type <= 2) {
        P <- sqrt(2/(n+h+0.5))*
           cos((2*pi/(2*(n+h)+1))*
               outer(seq(n)-0.5,seq(n)-0.5)
               )
         Pr <- sqrt(2/(n+h+0.5))*
           sin((2*pi/(2*(n+h)+1))*
               outer(seq(n)-0.5,seq(n))
               )

        Pm <- P[1:m,]
        Prm <- Pr[1:m,]

    Pmt <- t(Pm)
    Prmt <- t(Prm)

        P1 <- C %*% Pmt %*% Pm%*% Ci
    P2 <- Ct %*% Prmt %*% Prm %*% Cti
    T2 <- P2 + (In-P2)%*%Jn
        TT1 <- NULL
        TT2 <- NULL
        if(init) {
            if(type==1) return(P1 + (In-P1)%*%J1)
            if(type==2) return(P2 + (In-P2)%*%Jn)
            }
        for(i in 1:3) {
        T1 <- P1 + (In-P1)%*%J1%*%T2
        T2 <- P2 + (In-P2)%*%Jn%*%T1
            
        }
        if(trans) {
            if(allz) {
        P1 <-  P %*% Ci
        P2 <-  Pr %*% Cti                
            } else {
        P1 <-  Pm%*% Ci
        P2 <-  Prm%*% Cti
                }
        TT1 <- P1 %*%(In-J1%*%T2)
        TT2 <- P2 %*%(In-Jn%*%T1)
        
        }
        iTmat <- NULL
        if(type==2) {
                     if(inv) iTmat <- Ct %*% t(Pr)
                     Tmat <- list(TT=Jn%*%T1,K=TT2,inv=iTmat)
        } else {
                     if(inv) iTmat <- C %*% t(P)
                       Tmat <- list(TT=J1%*%T2,K=TT1,inv=iTmat)
         }
               }
         Tmat[["order"]] <- c(m,n,type)
    return(Tmat)
}
Global.Mat.SIML <- NULL
siml.mat <- function(n,m,h=0,type=1)
    {
       if(is.null(Global.Mat.SIML)||any(Global.Mat.SIML$order!=c(m,n,type))){ 
       Pmat <- makeMat21(n,m,type=type)
#           sqrt(2/(n+h+0.5))*
#           cos(outer(seq(m)-0.5,2*pi/(2*(n+h)+1)*(seq(n)-0.5)))
       assign("Global.Mat.SIML", Pmat, envir = .GlobalEnv)
       print("make mat")
       } else { 
           Pmat <- Global.Mat.SIML
       }
        return(Pmat)
        
    }
shouhi <-
c(59107.4, 58602.2, 60976.9, 61596.2, 60350.6, 59938.8, 62395.4, 
63518, 61754.1, 61378.3, 63484.6, 64735.9, 64139.4, 61427.1, 
63619.2, 63957.6, 62419.5, 61203.2, 63555.5, 63898, 62900.2, 
61850.6, 64152.4, 64686, 64256.1, 62906.9, 64781.2, 66107.5, 
65340.3, 64751.2, 66191.8, 66695.6, 66058.9, 65559.9, 67045.4, 
67494.1, 66857.5, 65599.4, 66820.9, 68410.5, 67619.5, 66513.7, 
68306.1, 68731.2, 68164.6, 67281.5, 69101.8, 69946.6, 69537.9, 
68149.2, 68798.1, 70757.4, 70143.3, 69165.2, 69744, 70923.3, 
70479.4, 68120.3, 69292.1, 69235.8, 67912.3, 67166.9, 69183.4, 
70680, 69679.2, 68717.8, 71507.2, 71440.6, 68858, 67964, 70777.1, 
72050.6, 71257.8, 69737.9, 71322.8, 72614.4, 72442.7, 71457, 
73364.5, 74382, 74934.9, 69698.1, 71757, 73003, 72325.1, 70729.7, 
72317.2, 72895.6, 72201.8, 70307.1, 71950.2, 73066.4, 72524.9, 
71800.1, 72435.7, 74008.6, 72759.6, 71832.3, 72839.5, 74528.9
)
kakeim <-
structure(c(13416.97289, 11948.66111, 14385.75921, 13748.61782, 
12764.66629, 12629.49807, 13790.31606, 13249.80611, 13235.38625, 
13343.0109, 12815.13181, 15507.36652, 13087.14874, 11982.08271, 
14330.12401, 13462.9204, 12739.95821, 12924.48286, 13357.93597, 
13331.12116, 12620.5265, 13177.10313, 12916.38631, 15917.45961, 
13315.90963, 12645.45921, 14060.62051, 14362.15316, 13584.85767, 
12789.64726, 13723.58875, 13614.15931, 12647.35236, 13214.38231, 
12813.18009, 15491.95308, 13580.97953, 12378.48524, 14065.32461, 
13994.67978, 13172.83187, 12777.22056, 13796.18213, 13656.82341, 
12942.77283, 13846.36576, 12928.29102, 15553.62848, 13308.94365, 
12050.32813, 13937.99206, 13774.54063, 13139.43402, 12694.16979, 
12986.35515, 12890.89073, 12133.15846, 12971.966, 12789.96918, 
15406.3735, 13555.06129, 12423.47733, 14250.60223, 14173.57005, 
13211.35025, 12908.12205, 13206.09954, 13338.19822, 12710.11075, 
13267.5608, 12681.13526, 15717.95308, 13793.90157, 12393.88848, 
13957.25322, 13928.78644, 13071.65446, 12594.78074, 13211.70869, 
13028.88705, 12449.38314, 13165.37052, 13070.73433, 15580.32334, 
13276.91339, 12042.81241, 13927.79295, 13847.58725, 13105.8423, 
12765.46333, 13146.13159, 13164.48056, 12534.12601, 13277.2788, 
13215.4376, 15828.08537, 13475.234, 12019.73694, 14431.05526, 
13747.76746, 13013.57975, 12876.11048, 13390.45893, 13615.0178, 
12654.77208, 13283.13259, 13451.20846, 15191.04231, 13431.48532, 
12281.59495, 13375.88353, 13487.4195, 12937.36718, 12561.4162, 
13325.25084, 13231.03534, 12517.51458, 13415.19801, 13059.35625, 
15426.09242, 13338.94061, 12543.23909, 13918.28025, 13832.9749, 
13256.49385, 12532.21639, 13255.77738, 13486.70532, 12361.76416, 
13485.83775, 12933.16002, 15450.16491, 13669.35499, 12739.55541, 
14827.2073, 13876.13074, 13056.8627, 12591.2437, 13453.5315, 
13463.59742, 13088.70305, 13278.02731, 12971.66705, 15679.64955, 
13583.11083, 12405.92182, 15968.4445, 13446.21062, 12308.67669, 
12472.18222, 12667.34341, 12994.2442, 12508.6883, 12851.36343, 
12865.58821, 15263.68223, 13199.78727, 12127.87057, 14071.137, 
13299.75564, 12969.3052, 12275.57049, 12727.41674, 13206.5006, 
12577.23806, 12826.07399, 12380.5491, 14732.52866, 12923.01365, 
12246.05435, 13558.52944, 13557.8737, 12921.0583, 12047.14139, 
12855.0639, 12869.23623, 12268.70246, 12736.32855, 12357.90503, 
14651.94261, 10996.90083, 10370.17727, 12397.5026, 11338.15926, 
11081.19219, 11103.09278, 11802.28928, 11007.26895, 10851.08583, 
11108.80829, 11410.97308, 13571.57676, 10957.24713, 10515.18325, 
12487.4739, 11257.52856, 11135.7513, 11217.93535, 11826.68067, 
11165.08938, 10973.93118, 11439.58333, 11442.81217, 13927.74869, 
11388.24764, 10971.66842, 12617.40042, 11478.53403, 11188.08777, 
11100, 12152.47108, 11148.57745, 11042.70833, 11295.57158, 11349.79424, 
13761.95426, 11610.01043, 10626.83438, 12535.49061, 11825.86027, 
11425.59834, 11508.40336, 12219.17808, 11306.96203, 11067.92059, 
11389.75967, 11613.20755, 14030.39832, 11547.44526, 10749.21301, 
12697.3822, 11725.46973, 11392.11618, 11374.48133, 12092.99896, 
11177.8697, 10956.74562, 11252.32678, 11508.8634, 13917.62252, 
11451.41066, 10776.60695, 12633.12369, 11731.73278, 11508.31601, 
11472.36705, 11931.93717, 11380.85328, 11145.96273, 11382.47423, 
11645.66116, 13853.75901, 11525.2838, 10968.91192, 12594.65021, 
11566.42636, 11259.63489, 11102.61569, 11684.73896, 11038.07615, 
10721.27872, 10983.98398, 11384.53713, 13452.77207, 11314.37435, 
10409.56341, 12220.49689, 11294.42149, 11183.22981, 11141.51925, 
11898.32285, 11313.47962, 11021.85224, 11375.91623, 11640, 13782.70042, 
11681.67539, 10905.66038, 12890.28213, 11864.10788, 11542.61954, 
11537.01773, 12495.24815, 11982.01058, 11322.14061, 11430.96234, 
11846.79958, 13532.20697, 11809.92608, 11031.67899, 11947.36842, 
11529.90556, 11614.4958, 11874.60485, 12601.90074, 11693.12169, 
11204.42571, 11654.04837, 11679.02542, 13931.07105, 11972.60274, 
11287.21174, 12906.34755, 12001.04058, 11921.46597, 11885.8351, 
12548.4558, 11947.70544, 11294.17989, 11571.88161, 11826.78002, 
13972.39915, 11885.71429, 11180.46709, 13087.83069, 12159.4509, 
12095.89041, 12036.88093, 12341.38655, 11870.93389, 11390.67358, 
11633.81743, 12008.29016, 13952.4302, 12133.54037, 11300.20704, 
14172.34262, 11060.24096, 11395.81256, 11334.33134, 11809.19081, 
11454, 11134.65347, 11345.30938, 11692.77108, 13564.25703, 11535.14056, 
10823.4107, 12452.81124, 11498.50746, 11666.99703, 11412.93532, 
12030, 11564.12826, 11176.1194, 11539.38185, 11562.249, 13424.1206, 
11560.48387, 10852.82258, 12382.05645, 11506.53266, 11575.3012, 
11392.35412, 12138.52376, 11453.34686, 11116.81772, 11520.43868, 
11664.34263, 13468.46847, 91.6, 90.9, 109.7, 96, 95.8, 96.9, 
99.6, 97.1, 100.5, 96.5, 96.8, 104.6, 92.1, 91.9, 111, 95.8, 
96.3, 98.1, 98.9, 96.5, 102.3, 99.5, 97.1, 106.7, 94, 94.3, 112.6, 
98.5, 96.3, 100.5, 101.4, 98.5, 103.1, 99.5, 99.7, 108.5, 95.9, 
95, 114.5, 100.1, 98.1, 102.3, 101.8, 101.7, 105.7, 101.6, 102.4, 
111.2, 98.4, 97.6, 116.6, 101.7, 101.5, 104.5, 103.8, 103.1, 
106.1, 103.8, 103.2, 111.9, 98.9, 98.8, 117.1, 103.4, 103.4, 
105.5, 105.5, 104.9, 106.4, 105, 104.7, 111.5, 99.3, 100.3, 116.3, 
102.9, 102.1, 103.6, 104.7, 101.2, 105, 102.5, 99.9, 108.3, 95.4, 
93.1, 108.1, 97, 95.7, 99.1, 99.5, 97.6, 100.4, 98.3, 97.1, 105.3, 
94.7, 93.1, 109.9, 98.6, 96.5, 99.9, 100.9, 99.8, 101.7, 98.7, 
99.4, 107, 95.9, 95.1, 105.8, 96.1, 96.2, 100.6, 101.1, 100.4, 
101.8, 99.9, 99.5, 108.5, 96.4, 97.5, 110.8, 99.3, 100, 102.1, 
102.5, 101.8, 102.1, 101.5, 100.5, 108.5, 97, 96.3, 112, 100.8, 
101.4, 102.3, 103.8, 102.6, 103.7, 101.7, 101.4, 109.4, 98.9, 
97.1, 115.3, 98.8, 99.7, 101.4, 102.3, 100.5, 103.3, 101.3, 99.5, 
109.4, 98.4, 97.4, 113.2, 101.4, 100.4, 103.8, 104.2, 102.5, 
104.1, 102.9, 100.9, 109.4, 98.6, 99.8, 113.6, 102.1, 100.8, 
104.3, 104.3, 103.5, 105, 102.6, 102.4, 110.2), .Dim = c(180, 
3))
"avg12" <- function(x,...) {
    n <- length(x)
    if(all(x > 0.1E-20)) x <- log(x)
    m1 <- as.integer(n/10)+1
    m2 <- as.integer(n/60)+1
    zz <- x12siml(x,trend=m1,sorder=m2,frequency=12,...)
    avg <- zz$adj[-(1:12)]-zz$trend[-((n-11):n)]
   return(avg)
}
