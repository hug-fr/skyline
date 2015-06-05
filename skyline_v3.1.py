import numpy
import math
import gdal
import psycopg2
import csv
import time
import conn_param
from shapely import geometry
from shapely import affinity
import os

start_time = time.time()

# load complete raster
img = gdal.Open('C:\ds_test_data\\i_france_mnt_2154.tif')
band1 = img.GetRasterBand(1)

rastinit = img.GetGeoTransform()

step = int((rastinit[1]+(-rastinit[5]))/2) #for further use in line interpolation

#x,y geographic reference matrix
imgx=numpy.zeros((1,img.RasterXSize)).astype(numpy.float)
imgy=numpy.zeros((img.RasterYSize,1)).astype(numpy.float)

for i in range(0,imgx.shape[1]):
	imgx[0,i]=rastinit[0]+(i*rastinit[1])


for i in range(0,imgy.shape[0]):
	imgy[i,0]=rastinit[3]+(i*rastinit[5])

#output to csv file

csv_out = "C:\Users\hugues.francois\Desktop\github\skyline\sta_skylines.csv"
if os.path.isfile(csv_out):
	os.remove(csv_out)

csvfile = open(csv_out,"wb")
stawriter = csv.writer(csvfile)

	

#Connect to DB

myconn = psycopg2.connect("host="+conn_param.host+" dbname="+conn_param.dbname+" user="+conn_param.user+" password="+conn_param.password)

#create table in postgres
skylinetable=myconn.cursor()
query="""
drop table if exists stations.station_meteo_skyline_v2;
"""
skylinetable.execute(query)
myconn.commit()

query="""
create table stations.station_meteo_skyline_v2(
	gid varchar(50),
	azimut int4,
	angle int4,
	pt_geom geometry);
"""
skylinetable.execute(query)
myconn.commit()
print "table created"

'''
#load points and view extent
wsta=myconn.cursor()

#query to select all point from a table
query="""
with a as (
	select distinct a.gid, st_x(a.the_geom) x, st_y(a.the_geom) y, st_buffer(a.the_geom,%s) geom
	from stations.geo_station_meteofrance a, spatial.geo_departements_fra b
	where st_intersects(a.the_geom, b.the_geom)
	)

select gid, x, y, st_xmin(geom) xmin, st_ymin(geom), st_xmax(geom) xmax, st_ymax(geom) ymax
from a
where gid = '38002406'
order by gid
;
"""

query="""
#test query or to compute skyline only for a single point
with test as (
	select 'Tignes-Pierre'::varchar gid, 2160::integer alt,
	st_transform(st_geomfromtext('POINT(6.896594 45.449133)', 4326),2154) the_geom 
	--st_transform(st_geomfromtext('POINT(5.765447 45.295108)', 4326),2154) the_geom 
	),
	a as (
	select gid, alt, st_x(the_geom) x, st_y(the_geom) y, st_buffer(the_geom,%s) geom
	from test
	)

select gid, x, y, st_xmin(geom) xmin, st_ymin(geom), st_xmax(geom) xmax, st_ymax(geom) ymax, alt
from a
--where gid = '38002406'
order by gid
;
"""
'''

viewmax= 20000#3.57*math.sqrt(215)*1000
#print viewmax
#wsta.execute(query,(viewmax,))

wsta=[["Tignes",1004464.1,6490813.1],
["Autrans-Prairie",900929.9,6460461.4],
["Autrans-Retenue d'eau",901261.0,6460502.4],
["2Alpes-Coolidge",946370.1,6439200.3],
["2Alpes-Lutins",946600.5,6439619.2],
["Chamrousse-Gabourreaux",926945.6,6451060.4],
["Chamrousse-Variante",927029.1,6450992.8],
["Chamrousse-Perche",926694.0,6450415.3]]


#extract from original raster

for sta in wsta:
	final_data = []
	
	#Find row/col information et xy normalization
	xmin = rastinit[0]+((math.floor(((sta[1]-viewmax)-rastinit[0])/rastinit[1]))*rastinit[1])
	xmax = rastinit[0]+((math.floor(((sta[1]+viewmax)-rastinit[0])/rastinit[1]))*rastinit[1])
	ymin = rastinit[3]-((math.ceil((rastinit[3]-(sta[2]-viewmax))/rastinit[5]))*rastinit[5])
	ymax = rastinit[3]-((math.ceil((rastinit[3]-(sta[2]+viewmax))/rastinit[5]))*rastinit[5])
	
	stax = rastinit[0]+((math.floor((sta[1]-rastinit[0])/rastinit[1]))*rastinit[1])
	stay = rastinit[3]-(math.ceil((rastinit[3]-sta[2])/rastinit[5])*rastinit[5])
	
	if ymax >= max(imgy):
		minrow = 0
	else:
		minrow = numpy.unique(numpy.argwhere(imgy==ymax))[1]
	if ymin <= min(imgy):
		maxrow = imgy.shape[0]
	else:
		maxrow = numpy.unique(numpy.argwhere(imgy==ymin))[1]

	if xmin <= min(imgx[0,]):
		mincol=0
	else:
		mincol = numpy.unique(numpy.argwhere(imgx==xmin))[1]
	if xmax >= max(imgx[0,]):
		maxcol = imgx.shape[1]
	else:
		maxcol = numpy.unique(numpy.argwhere(imgx==xmax))[1]

	starow = maxrow-numpy.unique(numpy.argwhere(imgy==stay))[1]
	stacol = numpy.unique(numpy.argwhere(imgx==stax))[1]-mincol
	starow = starow.astype('int64')
	stacol = stacol.astype('int64')
	
	sta_xy = (stax+(rastinit[1]/2), stay+(rastinit[5]/2))
	sta_rc = (starow, stacol)
	
	#Extract array from raster
	height = band1.ReadAsArray(mincol, minrow, maxcol-mincol, maxrow-minrow)
	height = height.astype('int64')

	# get width and heigth of image
	w,h = height.shape
	
	print "raster extracted", w, h
	print height[sta_rc]
	
	#Get all intersected cells on azimuth
	
	for azimut in range (0, 360, 5):
		i = 0
		angle = numpy.zeros((1,(viewmax/step)-1)).astype(numpy.float) #initialize container for angles
		points = [] #initialize container for points
		pt_dist = []
		for dist in range (step,viewmax,step):
			ptx = sta[1]+(dist*math.sin(math.radians(azimut))) 
			pty = sta[2] + (dist * math.cos(math.radians(azimut)))
			pt = (ptx,pty)
			points.append(pt)
			pt_dist.append(dist)
						
			#get row col information
			if ptx < xmax and ptx > xmin:
				x = rastinit[0]+((math.floor((ptx-rastinit[0])/rastinit[1]))*rastinit[1])
				ptcol = numpy.unique(numpy.argwhere(imgx==x))[1]-mincol
			if pty < ymax and pty > ymin:
				y = rastinit[3]-((math.ceil((rastinit[3]-pty)/rastinit[5]))*rastinit[5])
				ptrow = numpy.unique(numpy.argwhere(imgy==y))[1]-minrow
			ptrc=(ptrow, ptcol)
			#print dist, height[ptrc]-height[sta_rc], x,y, stax, stay, ptrc, sta_rc
					
			#calculate corresponding angle to reach the height of pt
			if ptrow < w and ptcol < h:
				b = height[ptrc]-height[sta_rc] #sta[7]
				b = b.astype('float')
				#print b, b/dist, type(b), type(dist), type(b/dist)#)))*100)/100
				if b > 0:					
					angle[0,i]=math.ceil((math.degrees(math.atan(b/dist)))*100)/100
				else:
					angle[0,i]=0
			
			#print angle[0,i], max(angle[0,])
			#raw_input()
			
			i = i+1
		
		print sta[0], azimut, max(angle[0,])
		
				
		#append each azimut to final data for weather station
		data = (sta[0], azimut, max(angle[0,]), points[numpy.argwhere(angle==max(angle[0,]))[0][1]][0], points[numpy.argwhere(angle==max(angle[0,]))[0][1]][1])
		final_data.append(data)
		#print(data)

	final_data.append((final_data[0][0], 360, final_data[0][2], final_data[0][3], final_data[0][4]))

	#insert values into new table for the given weather station
	
	for values in final_data:
		
		skylinetable=myconn.cursor()
		query="""
		insert into stations.station_meteo_skyline_v2
		values(%s, %s, %s, ST_SetSRID(ST_MakePoint(%s, %s), 2154));
		"""

		skylinetable.execute(query,(values[0],values[1], values[2], values[3], values[4]))
		myconn.commit()
		
		stawriter.writerow([values[0],values[1],values[2]])		
	
	final_data = None
	print sta[0], "done"

csvfile.close()
print"done in", time.time()-start_time, "seconds"