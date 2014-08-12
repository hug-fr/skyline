# INTRODUCTION DONNEES ATTRIBUTAIRES SQL
library (RODBC)
library (plotrix)
library (reshape)
 

ch = odbcConnect("PostgreSQL35W",uid="hug",pwd="hug")
 
sql=
"
select nom_station, azimut, angle from stations.station_meteo_skyline a
join stations.geo_station_meteofrance b on a.gid=b.gid
"


data=sqlQuery(ch, paste(sql, collapse=' '))
close(ch)

mylim <- ceiling(max(data[,3])/20)*20

data2 <- cast(data, nom_station ~ azimut, value = "angle")

for (sta in data2[,1]){

dev.new()
mylay<-layout(matrix(c(1,2,2,2,2,2,2,2,2,2,1,3,3,3,3,3,3,3,3,3),10,2))
#layout.show(mylay)
par(mar = c(1,0,2,0))
plot.new()
mtext(sta)
par(mar=c(2,2,0,1))
radial.plot(
	mylim-data2[data2[,1]==sta,2:length(data2)],
	labels=c("N","NE","E","SE","S","SW","W","NW"),
	rp.type="p",
	#main= sta,
	radial.lim=c(0,mylim),
	line.col="#648bda",
	lwd = 2,
	start=1.56,
	clockwise = T)

par(mar=c(2,4,0,1))
plot(data[data[,1]==sta,2], data[data[,1]==sta,3],
	type = "l",
	col="#648bda",
	lwd = 2,
	xlab = "Azimuts (degrees)",
	ylab = "Skyline angle (degrees)"
)
}