library("cluster")

dropPercent <- function(x, y, percent, timesteps){
	#x is the x coordinate
	#y is the y coordinate
	#percent is what percent of data points you want to keep
	#timesteps is the number of timesteps per years-interval
	#returns in the inner x%	
	dat <- data.frame(x, y)
	dat$r <- sqrt(x^2 + y^2)
	dat$perc <- rep(0, length(dat$r))
	for(i in 1:timesteps){
		dat$perc[seq(i, length(dat$r), by = timesteps)] = rank(dat$r[seq(i, length(dat$r), by = timesteps)])/length(dat$r[seq(i, length(dat$r), by = timesteps)])
	}
	dat1 <- dat[dat$perc <= percent,]
	dat2 <- subset(dat1, select = c(x, y))
	return(dat2)
}
casesToXAndY <- function(day, cases, dayType = "m"){
	#day is a date corresponding to cases - if dayType is "m", day is 1-12, if dayType is "d", day is 1-365
	#cases is the number of cases
	#dayType is an indicator for whether day is in months or days
	if(dayType == "d"){
		cases.norm <- rep(0, length(cases))
		for(i in 1:length(cases)){
			cases.norm[i] <- cases[i]/max(cases[max(1, i-183):min(length(cases), i+183)])
		}
		exes <- cases.norm*cos((day-1)/365)
		yies <- cases.norm*sin((day-1)/365)
	}
	if(dayType == "m"){
		cases.norm <- rep(0, length(cases))
		for(i in 1:length(cases)){
			cases.norm[i] <- cases[i]/max(cases[max(1, i-6):min(length(cases), i+6)])
		}
		exes <- cases.norm*cos((day-1)/12)
		yies <- cases.norm*sin((day-1)/12)
	}
	return(list(x = exes, y = yies))
}
casesToEllipse <- function(cases, days, dayType = "m", interval, percent, years, distribution = FALSE){
	cases.norm <- casesToXAndY(cases, days, dayType)
	if(dayType == "m"){
		cases.norm <- dropPercent(cases.norm$x, cases.norm$y, 1-percent, 12)
	}
	if(dayType == "d"){
		cases.norm <- dropPercent(cases.norm$x, cases.norm$y, 1-percent, 365)
	}
	areas <- rep(0, (years - interval + 1))
	ratio <- rep(0, (years - interval + 1))
	for(i in 1:(years - interval + 1)){
		if(dayType == "m"){
			dat.x <- cases.norm$x[(1+(i-1)*12):((i+interval-1)*12)]
			dat.y <- cases.norm$y[(1+(i-1)*12):((i+interval-1)*12)]
		}
		if(dayType == "d"){
			dat.x <- cases.norm$x[(1+(i-1)*365):((i+interval-1)*365)]
			dat.y <- cases.norm$y[(1+(i-1)*365):((i+interval-1)*365)]
		}
		dat <- matrix(c(dat.x, dat.y), ncol = 2, nrow = length(dat.x), byrow = FALSE)
		exy <- predict(ellipsoidhull(as.matrix(dat)))
		me <- colMeans((exy))
		dist2center <- sqrt(rowSums((t(t(exy)-me))^2))
		width <- max(dist2center)*2
		height <- min(dist2center)*2
		areas[i] <- pi*width*height/4
		ratio[i] = height/width
	}
	if(distribution){
		areas.mean <- rep(0, (years - interval + 1))
		ratio.mean <- rep(0, (years - interval + 1))
		areas.var <- rep(0, (years - interval + 1))
		ratio.var <- rep(0, (years - interval + 1))
		
		for(i in 1:(years - interval + 1)){
			if(dayType == "m"){
				dat.x <- cases.norm$x[(1+(i-1)*12):((i+interval-1)*12)]
				dat.y <- cases.norm$y[(1+(i-1)*12):((i+interval-1)*12)]
			}
			if(dayType == "d"){
				dat.x <- cases.norm$x[(1+(i-1)*365):((i+interval-1)*365)]
				dat.y <- cases.norm$y[(1+(i-1)*365):((i+interval-1)*365)]
			}
			areas.mat <- rep(0, 100)
			ratios.mat <- rep(0, 100)
			for(j in 1:100){
				dat <- matrix(c(dat.x, dat.y), ncol = 2, nrow = length(dat.x), byrow = FALSE)
				dat <- dat[sample(1:nrow(dat), floor(nrow(dat)/2)),]
				exy <- predict(ellipsoidhull(as.matrix(dat)))
				me <- colMeans((exy))
				dist2center <- sqrt(rowSums((t(t(exy)-me))^2))
				width <- max(dist2center)*2
				height <- min(dist2center)*2
				areas.mat[j] <- pi*width*height/4
				ratios.mat[j] = height/width
			}
			areas.mean[i] = mean(areas.mat)
			areas.var[i] = var(areas.mat)
			ratio.mean[i] = mean(ratios.mat)
			ratio.var[i] = var(ratios.mat)
		}
		
		return(list(areas, ratio, areas.mean, areas.var, ratio.mean, ratio.var))
	}
	if(!distribution){
		return(list(areas, ratio))
	}
}
