import numpy
import math
import gdal
import psycopg2
import time
import conn_param
from shapely import geometry
from shapely import affinity

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
	pt_geom geometry,
	line_geom geometry);
"""
skylinetable.execute(query)
myconn.commit()
print "table created"

#load points and view extent
wsta=myconn.cursor()
'''
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
'''
query="""
with test as (
	select 'test'::varchar gid, st_transform(st_geomfromtext('POINT(1.8105 42.577167)', 4326),2154) the_geom
	),
	a as (
	select gid, st_x(the_geom) x, st_y(the_geom) y, st_buffer(the_geom,%s) geom
	from test
	)

select gid, x, y, st_xmin(geom) xmin, st_ymin(geom), st_xmax(geom) xmax, st_ymax(geom) ymax
from a
--where gid = '38002406'
order by gid
;
"""

viewmax= 50000#3.57*math.sqrt(215)*1000
print viewmax
wsta.execute(query,(viewmax,))

#extract from original raster

for sta in wsta:
	final_data = []
	
	#Find row/col information et xy normalization
	xmin = rastinit[0]+((math.floor((sta[3]-rastinit[0])/rastinit[1]))*rastinit[1])
	xmax = rastinit[0]+((math.floor((sta[5]-rastinit[0])/rastinit[1]))*rastinit[1])
	ymin = rastinit[3]-((math.ceil((rastinit[3]-sta[4])/rastinit[5]))*rastinit[5])
	ymax = rastinit[3]-((math.ceil((rastinit[3]-sta[6])/rastinit[5]))*rastinit[5])
	
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

	#Init ref line (north)
	
	refline=geometry.LineString([(sta[1],sta[2]), (sta[1],sta[2]+viewmax)])
	refpoint=geometry.Point(sta[1], sta[2])
	
	#Get all intersected cells on azimuth
	
	for azimut in range (0, 360, 5):
		newline = affinity.rotate(refline, -azimut, origin=(sta[1],sta [2])) #rotate line regarding azimut
		i = 0 #count for dist iterations
		angle = numpy.empty((1,viewmax/(2*step))).astype(numpy.float) #initialize container for angles
		points = [] #initialize container for points
		pt_dist = []
		
		for dist in range (0,viewmax,step*2):
			i = i+1
			pt = newline.interpolate(dist) #point on line
			
			#get row col information
			if pt.x < xmax and pt.x > xmin:
				x = rastinit[0]+((math.floor((pt.x-rastinit[0])/rastinit[1]))*rastinit[1])
				ptcol = numpy.unique(numpy.argwhere(imgx==x))[1]-mincol
			if pt.y < ymax and pt.y > ymin:
				y = rastinit[3]-((math.ceil((rastinit[3]-pt.y)/rastinit[5]))*rastinit[5])
				ptrow = maxrow-numpy.unique(numpy.argwhere(imgy==y))[1]
			ptrc=(ptrow, ptcol)				
			#calculate corresponding angle to reach the height of pt
			if ptrow < w and ptcol < h:
				b = height[ptrc]-height[sta_rc]
				a = round(refpoint.distance(pt),0)
				if b > 0:
					
					angle[0,i-1]=math.ceil((math.degrees(math.acos(a/math.sqrt(a**2+b**2))))*100)/100
					points.append(pt.wkb_hex)
					pt_dist.append(a)
				else:
					angle[0,i-1]=0
					points.append(pt.wkb_hex)
					pt_dist.append(a)
		
		print sta[0], azimut, len(points), len(angle[0,]), numpy.argwhere(angle==max(angle[0,]))[0][1], pt_dist[numpy.argwhere(angle==max(angle[0,]))[0][1]]
		
				
		#append each azimut to final data for weather station
		data = (sta[0], azimut, max(angle[0,]), points[numpy.argwhere(angle==max(angle[0,]))[0][1]], newline.wkb_hex)
		final_data.append(data)
		#print(data)

	final_data.append((final_data[0][0], 360, final_data[0][2], final_data[0][3], final_data[0][4]))

	#insert values into new table for the given weather station
	for values in final_data:
		skylinetable=myconn.cursor()
		query="""
		insert into stations.station_meteo_skyline_v2
		values(%s, %s, %s, ST_SetSRID(%s, 2154), ST_SetSRID(%s, 2154));
		"""

		skylinetable.execute(query,(values[0],values[1], values[2], values[3], values[4]))
		myconn.commit()
	
	final_data = None
	print sta[0], "done"

print"done in", time.time()-start_time, "seconds"