# INTRODUCTION DONNEES ATTRIBUTAIRES SQL
library (RODBC)
library (plotrix)
library (reshape)		
 

ch = odbcConnect("PostgreSQL35W",uid="hug",pwd="hug")
 
sql=
"
/*
select a.gid, nom_station, azimut, angle from stations.station_meteo_skyline a
join stations.geo_station_meteofrance b on a.gid=b.gid
*/
select gid, gid nom_station, azimut, angle from stations.station_meteo_skyline
"


data=sqlQuery(ch, paste(sql, collapse=' '))
close(ch)

ch = odbcConnect("PostgreSQL35W",uid="hug",pwd="hug")
 
sql= "select gid, gid nom_station, azimut, angle from stations.station_meteo_skyline_v2"

datav2=sqlQuery(ch, paste(sql, collapse=' '))
close(ch)

mylim <- ceiling(max(data[,4])/20)*20

#data2 <- cast(data, nom_station ~ azimut, value = "angle")

for (sta in unique(data[,1])){

png(
	file=paste0("C:/Users/hugues.francois/Desktop/github/skyline/graph_",sta,".png"),
	width = 800, height=600, res=150
)
mylay<-layout(matrix(c(
1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
1,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3)
,18,2))
#layout.show(mylay)
par(mar = c(0,0,0,0))
plot.new()
mtext(unique(data[data[,1]==sta,2]), cex=1, line =-4)
par(mar=c(2,2,0,1), cex.axis=.6)
data2 <- cast(data[data[,1]==sta,2:4], nom_station ~ azimut, value = "angle")
radial.plot(
	mylim-data2[,2:length(data2)],
	labels=c("N","NE","E","SE","S","SW","W","NW"),
	rp.type="p",
	radial.lim=c(0,mylim),
	radial.labels=rev(pretty(c(0,mylim))),
	boxed.radial=F,
	grid.unit ="°",
	line.col="#648bda",
	lwd = 2,
	start=pi/2,
	clockwise = T,
	poly.col="#648bda50"
)

par(new=T)
datav22 <- cast(datav2[datav2[,1]==sta,2:4], nom_station ~ azimut, value = "angle")
radial.plot(
	mylim-datav22[,2:length(datav22)],
	labels=c(),
	rp.type="p",
	radial.lim=c(0,mylim),
	radial.labels=rev(pretty(c(0,mylim))),
	boxed.radial=F,
	grid.unit ="°",
	line.col="#2ca25f",
	lwd = 2,
	start=pi/2,
	clockwise = T,
	poly.col="#2ca25f50"
)

par(mar=c(2,4,0,1))
plot(data[data[,1]==sta,3], data[data[,1]==sta,4],
	type = "l",
	col="#648bda",
	lwd = 2,
	axes=F,
	ylim=c(0,mylim),
	xlab = NA,
	ylab = NA,
	xlim=c(0,360)
)

axis(side = 1, tck = -.02, labels = NA)
axis(side=1, line = -.8, lwd = 0, cex.axis =.7, font.lab=2)
mtext(side=1, line=1.2, cex=.6, "Azimuts (degrees)")

axis(side = 2, tck = -.02, labels = NA)
axis(side=2, line = -.6, lwd = 0, cex.axis =.7, font.lab=2)
mtext(side=2, line=1.4, cex=.6, "Skyline angle (degrees)")

#dev.off()
}

if(FALSE){
data3 <-read.csv("C:/python_script/tmp.csv",header=F,sep=";")
data3<-data3[order(data3[,1]),]
par(new=T)
plot(data3[,1], data3[,2], type = "l", col="red",
	xlim=c(0,360), ylim=c(0,mylim), xlab=NA, ylab=NA)
}
ch = odbcConnect("PostgreSQL35W",uid="hug",pwd="hug")
sql = "select azimut, angle from stations.station_meteo_skyline_v2"
data4=sqlQuery(ch, paste(sql, collapse=' '))
close(ch)
par(new=T)
plot(data4[,1], data4[,2], type = "l", col="green",
	xlim=c(0,360), ylim=c(0,mylim), xlab=NA, ylab=NA)

data3 <-read.csv("C:/Users/hugues.francois/Desktop/github/skyline/test.csv",header=F,sep=";")
par(new=T)
plot(t(data3[1,2:length(data3[1,])]), t(data3[2,2:length(data3[2,])]),
	type = "l", col="red",
	xlim=c(0,360), ylim=c(0,mylim), xlab=NA, ylab=NA)
par(new=T)
plot(t(data3[1,2:length(data3[1,])]), t(data3[3,2:length(data3[2,])]),
	type = "l", col="orange",
	xlim=c(0,360), ylim=c(0,mylim), xlab=NA, ylab=NA)

dev.off()
