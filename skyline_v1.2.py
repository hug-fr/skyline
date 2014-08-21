import numpy
import math
import gdal
import psycopg2
import time
import conn_param

start_time = time.time()

# load complete raster
img = gdal.Open('C:\ds_test_data\\i_france_mnt_2154.tif')
band1 = img.GetRasterBand(1)

rastinit = img.GetGeoTransform()

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
drop table if exists stations.station_meteo_skyline;
"""
skylinetable.execute(query)
myconn.commit()

query="""
create table stations.station_meteo_skyline(
	gid varchar(50),
	azimut int4,
	angle int4);
"""
skylinetable.execute(query)
myconn.commit()
print "table created"

#load points and view extent
wsta=myconn.cursor()
'''
query="""
with a as (
	select gid, st_x(the_geom) x, st_y(the_geom) y, st_buffer(the_geom,%s) geom
	from stations.geo_station_meteofrance
	)

select gid, x, y, st_xmin(geom) xmin, st_ymin(geom), st_xmax(geom) xmax, st_ymax(geom) ymax
from a
where gid = '09139401'
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
	
	#Get all intersected cells on azimuth
	for azimut in range (0, 360, 5):
		print sta[0], azimut
		myview=0
		c= None
		#TESTS on azimut
		
######################################AZIMUT 0 / 360######################################

		if azimut == 0:
			imax = int(math.floor(viewmax/25))
			angle = []
			for i in range(1,imax):
			#incrementation depending on azimut value (take care of x,y and row, col order and sign for ptx (or pty) calcul
				ptrow = sta_rc[0]-i
				ptcol = sta_rc[1]
				ptrc=(ptrow, ptcol)
				
				#calculate corresponding angle to reach the corresponding height
				b = height[ptrc]-height[sta_rc]
				if b > 0:
					#take care of row col and order between viewpt and sta
					a = ((sta_rc[0]-ptrc[0])*25)
					angle.append(math.ceil((math.degrees(math.acos(a/math.sqrt(a**2+b**2))))*100)/100)
				else:
					angle.append(0)
			
######################################AZIMUT 180######################################

		if azimut == 180:
			imax = int(math.floor(viewmax/25))
			angle = []
			for i in range(1,imax):
			#incrementation depending on azimut value (take care of x,y and row, col order and sign for ptx (or pty) calcul
				ptrow = sta_rc[0]+i
				ptcol = sta_rc[1]
				ptrc=(ptrow, ptcol)
				
				#calculate corresponding angle to reach the corresponding height
				b = height[ptrc]-height[sta_rc]
				if b > 0:
					#take care of row col and order between viewpt and sta
					a = ((ptrc[0]-sta_rc[0])*25)
					angle.append(math.ceil((math.degrees(math.acos(a/math.sqrt(a**2+b**2))))*100)/100)
				else:
					angle.append(0)
			
######################################AZIMUT 90######################################

		if azimut == 90:
			imax = int(math.floor(viewmax/25))
			angle = []
			for i in range(1,imax):
				#incrementation depending on azimut value (take care of x,y and row, col order
				ptcol = sta_rc[1]+i
				ptrow=sta_rc[0]
				ptrc=(ptrow, ptcol)
				
				#calculate corresponding angle to reach the corresponding height
				b = height[ptrc]-height[sta_rc]
				if b > 0:
					#take care of row col and order between viewpt and sta
					a = ((ptrc[1]-sta_rc[1])*25)
					angle.append(math.ceil((math.degrees(math.acos(a/math.sqrt(a**2+b**2))))*100)/100)
				else:
					angle.append(0)

######################################AZIMUT 270######################################

		if azimut == 270:
			imax = int(math.floor(viewmax/25))
			angle = []
			for i in range(1,imax):
				#incrementation depending on azimut value (take care of x,y and row, col order
				ptcol = sta_rc[1]-i
				ptrow=sta_rc[0]
				ptrc=(ptrow, ptcol)
				
				#calculate corresponding angle to reach the corresponding height
				b = height[ptrc]-height[sta_rc]
				if b > 0:
					#take care of row col and order between viewpt and sta
					a = ((sta_rc[1]-ptrc[1])*25)
					angle.append(math.ceil((math.degrees(math.acos(a/math.sqrt(a**2+b**2))))*100)/100)
				else:
					angle.append(0)

######################################AZIMUT 45######################################

		if azimut == 45:
			#calculate corresponding angle of right-angled triangle
			alpha = math.radians(45)
			cos_alpha = math.cos(alpha)
			imax = int(math.floor(viewmax*cos_alpha/25))
			angle = []
			for i in range(1,imax):
				#incrementation depending on azimut value (take care of x,y and row, col order
				ptcol = sta_rc[1]+i
				ptrow = sta_rc[0]-i
				ptrc=(ptrow, ptcol)
				
				#calculate corresponding angle to reach the corresponding height
				b = height[ptrc]-height[sta_rc]
				if b > 0:
					#take care of row col and order between viewpt and sta
					a = ((ptrc[1]-sta_rc[1])*25)/cos_alpha
					angle.append(math.ceil((math.degrees(math.acos(a/math.sqrt(a**2+b**2))))*100)/100)
				else:
					angle.append(0)

######################################AZIMUT 135######################################

		if azimut == 135:
			#calculate corresponding angle of right-angled triangle
			alpha = math.radians(45)
			cos_alpha = math.cos(alpha)
			imax = int(math.floor(viewmax*cos_alpha/25))
			angle = []
			for i in range(1,imax):
				#incrementation depending on azimut value (take care of x,y and row, col order
				ptcol = sta_rc[1]+i
				ptrow = sta_rc[0]+i
				ptrc=(ptrow, ptcol)
				
				#calculate corresponding angle to reach the corresponding height
				b = height[ptrc]-height[sta_rc]
				if b > 0:
					#take care of row col and order between viewpt and sta
					a = ((ptrc[1]-sta_rc[1])*25)/cos_alpha
					angle.append(math.ceil((math.degrees(math.acos(a/math.sqrt(a**2+b**2))))*100)/100)
				else:
					angle.append(0)
					
######################################AZIMUT 225######################################

		if azimut == 225:
			#calculate corresponding angle of right-angled triangle
			alpha = math.radians(45)
			cos_alpha = math.cos(alpha)
			imax = int(math.floor(viewmax*cos_alpha/25))
			angle = []
			for i in range(1,imax):
				#incrementation depending on azimut value (take care of x,y and row, col order
				ptcol = sta_rc[1]-i
				ptrow = sta_rc[0]+i
				ptrc=(ptrow, ptcol)
				
				#calculate corresponding angle to reach the corresponding height
				b = height[ptrc]-height[sta_rc]
				if b > 0:
					#take care of row col and order between viewpt and sta
					a = ((ptrc[0]-sta_rc[0])*25)/cos_alpha
					angle.append(math.ceil((math.degrees(math.acos(a/math.sqrt(a**2+b**2))))*100)/100)
				else:
					angle.append(0)

######################################AZIMUT 315######################################

		if azimut == 315:
			#calculate corresponding angle of right-angled triangle
			alpha = math.radians(45)
			cos_alpha = math.cos(alpha)
			imax = int(math.floor(viewmax*cos_alpha/25))
			angle = []
			for i in range(1,imax):
				#incrementation depending on azimut value (take care of x,y and row, col order
				ptcol = sta_rc[1]-i
				ptrow = sta_rc[0]-i
				ptrc=(ptrow, ptcol)
				
				#calculate corresponding angle to reach the corresponding height
				b = height[ptrc]-height[sta_rc]
				if b > 0:
					#take care of row col and order between viewpt and sta
					a = ((sta_rc[0]-ptrc[0])*25)/cos_alpha
					angle.append(math.ceil((math.degrees(math.acos(a/math.sqrt(a**2+b**2))))*100)/100)
				else:
					angle.append(0)

######################################AZIMUTS 0 => 45######################################

		if azimut > 0 and azimut < 45:
		#calculate corresponding angle of right-angled triangle
			alpha = math.radians(azimut)
			cos_alpha = math.cos(alpha)
			sin_alpha = math.sin(alpha)
			imax = int(math.floor(viewmax*cos_alpha/25))
			angle = []
			for i in range(1,imax):
				myview = myview+25
				c = myview/cos_alpha
				#incrementation depending on azimut value (take care of x,y and row, col order and sign for ptx (or pty) calcul
				ptrow = sta_rc[0]-i
				ptx = rastinit[0]+((math.floor((sta_xy[0]+(sin_alpha*c)-rastinit[0])/rastinit[1]))*rastinit[1])
				ptcol = numpy.unique(numpy.argwhere(imgx==ptx))[1]-mincol
				ptrc=(ptrow, ptcol)
				
				#calculate corresponding angle to reach the corresponding height
				b = height[ptrc]-height[sta_rc]
				if b > 0:
					#take care of row col and order between viewpt and sta
					a = ((sta_rc[0]-ptrc[0])*25)/cos_alpha
					angle.append(math.ceil((math.degrees(math.acos(a/math.sqrt(a**2+b**2))))*100)/100)
				else:
					angle.append(0)

######################################AZIMUTS 45 => 90######################################

		if azimut > 45 and azimut < 90:
			#calculate corresponding angle of right-angled triangle
			alpha = math.radians(90-azimut)
			cos_alpha = math.cos(alpha)
			sin_alpha = math.sin(alpha)
			imax = int(math.floor(viewmax*cos_alpha/25))
			angle = []
			for i in range(1,imax):
				myview = myview+25
				c = myview/cos_alpha
				#incrementation depending on azimut value (take care of x,y and row, col order
				ptcol = sta_rc[1]+i
				pty=rastinit[3]-((math.ceil((rastinit[3]-sta_xy[1]+(sin_alpha*c))/rastinit[5]))*rastinit[5])
				ptrow=maxrow-numpy.unique(numpy.argwhere(imgy==pty))[1]
				ptrc=(ptrow, ptcol)
				
				#calculate corresponding angle to reach the corresponding height
				b = height[ptrc]-height[sta_rc]
				if b > 0:
					#take care of row col and order between viewpt and sta
					a = ((ptrc[1]-sta_rc[1])*25)/cos_alpha
					angle.append(math.ceil((math.degrees(math.acos(a/math.sqrt(a**2+b**2))))*100)/100)
				else:
					angle.append(0)

######################################AZIMUTS 90 => 135######################################

		if azimut > 90 and azimut < 135:
			#calculate corresponding angle of right-angled triangle
			alpha = math.radians(azimut-90)
			cos_alpha = math.cos(alpha)
			sin_alpha = math.sin(alpha)
			imax = int(math.floor(viewmax*cos_alpha/25))
			angle = []
			for i in range(1,imax):
				myview = myview+25
				c = myview/cos_alpha
				#incrementation depending on azimut value (take care of x,y and row, col order
				ptcol = sta_rc[1]+i
				pty=rastinit[3]-((math.ceil((rastinit[3]-sta_xy[1]-(sin_alpha*c))/rastinit[5]))*rastinit[5])
				ptrow=maxrow-numpy.unique(numpy.argwhere(imgy==pty))[1]
				ptrc=(ptrow, ptcol)
				
				#calculate corresponding angle to reach the corresponding height
				b = height[ptrc]-height[sta_rc]
				if b > 0:
					#take care of row col and order between viewpt and sta
					a = ((ptrc[1]-sta_rc[1])*25)/cos_alpha
					angle.append(math.ceil((math.degrees(math.acos(a/math.sqrt(a**2+b**2))))*100)/100)
				else:
					angle.append(0)

######################################AZIMUTS 135 => 180######################################

		if azimut > 135 and azimut < 180:
			#calculate corresponding angle of right-angled triangle
			alpha = math.radians(180-azimut)
			cos_alpha = math.cos(alpha)
			sin_alpha = math.sin(alpha)
			imax = int(math.floor(viewmax*cos_alpha/25))
			angle = []
			for i in range(1,imax):
				myview = myview+25
				c = myview/cos_alpha
				#incrementation depending on azimut value (take care of x,y and row, col order
				ptrow = sta_rc[0]+i
				ptx = rastinit[0]+((math.floor((sta_xy[0]+(sin_alpha*c)-rastinit[0])/rastinit[1]))*rastinit[1])
				ptcol = numpy.unique(numpy.argwhere(imgx==ptx))[1]-mincol
				ptrc=(ptrow, ptcol)
				
				#calculate corresponding angle to reach the corresponding height
				b = height[ptrc]-height[sta_rc]
				if b > 0:
					#take care of row col and order between viewpt and sta
					a = ((ptrc[0]-sta_rc[0])*25)/cos_alpha
					angle.append(math.ceil((math.degrees(math.acos(a/math.sqrt(a**2+b**2))))*100)/100)
				else:
					angle.append(0)
			
######################################AZIMUTS 180 => 225######################################

		if azimut > 180 and azimut < 225:
		#calculate corresponding angle of right-angled triangle
			alpha = math.radians(azimut - 180)
			cos_alpha = math.cos(alpha)
			sin_alpha = math.sin(alpha)
			imax = int(math.floor(viewmax*cos_alpha/25))
			angle = []
			for i in range(1,imax):
				myview = myview+25
				c = myview/cos_alpha
			#incrementation depending on azimut value (take care of x,y and row, col order and sign for ptx (or pty) calcul
				ptrow = sta_rc[0]+i
				ptx = rastinit[0]+((math.floor((sta_xy[0]-(sin_alpha*c)-rastinit[0])/rastinit[1]))*rastinit[1])
				ptcol = numpy.unique(numpy.argwhere(imgx==ptx))[1]-mincol
				ptrc=(ptrow, ptcol)
				
				#calculate corresponding angle to reach the corresponding height
				b = height[ptrc]-height[sta_rc]
				if b > 0:
					#take care of row col and order between viewpt and sta
					a = ((ptrc[0]-sta_rc[0])*25)/cos_alpha
					angle.append(math.ceil((math.degrees(math.acos(a/math.sqrt(a**2+b**2))))*100)/100)
				else:
					angle.append(0)
					
######################################AZIMUTS 225 => 270######################################

		if azimut > 225 and azimut < 270:
			#calculate corresponding angle of right-angled triangle
			alpha = math.radians(270 - azimut)
			cos_alpha = math.cos(alpha)
			sin_alpha = math.sin(alpha)
			imax = int(math.floor(viewmax*cos_alpha/25))
			angle = []
			for i in range(1,imax):
				myview = myview+25
				c = myview/cos_alpha
				#incrementation depending on azimut value (take care of x,y and row, col order
				ptcol = sta_rc[1]-i
				pty=rastinit[3]-((math.ceil((rastinit[3]-sta_xy[1]-(sin_alpha*c))/rastinit[5]))*rastinit[5])
				ptrow=maxrow-numpy.unique(numpy.argwhere(imgy==pty))[1]
				ptrc=(ptrow, ptcol)
				
				#calculate corresponding angle to reach the corresponding height
				b = height[ptrc]-height[sta_rc]
				if b > 0:
					#take care of row col and order between viewpt and sta
					a = ((sta_rc[1]-ptrc[1])*25)/cos_alpha
					angle.append(math.ceil((math.degrees(math.acos(a/math.sqrt(a**2+b**2))))*100)/100)
				else:
					angle.append(0)
			
######################################AZIMUTS 270 => 315######################################

		if azimut > 270 and azimut < 315:
			#calculate corresponding angle of right-angled triangle
			alpha = math.radians(azimut-270)
			cos_alpha = math.cos(alpha)
			sin_alpha = math.sin(alpha)
			imax = int(math.floor(viewmax*cos_alpha/25))
			angle = []
			for i in range(1,imax):
				myview = myview+25
				c = myview/cos_alpha
				#incrementation depending on azimut value (take care of x,y and row, col order
				ptcol = sta_rc[1]-i
				pty=rastinit[3]-((math.ceil((rastinit[3]-sta_xy[1]+(sin_alpha*c))/rastinit[5]))*rastinit[5])
				ptrow=maxrow-numpy.unique(numpy.argwhere(imgy==pty))[1]
				ptrc=(ptrow, ptcol)
				
				#calculate corresponding angle to reach the corresponding height
				b = height[ptrc]-height[sta_rc]
				if b > 0:
					#take care of row col and order between viewpt and sta
					a = ((sta_rc[1]-ptrc[1])*25)/cos_alpha
					angle.append(math.ceil((math.degrees(math.acos(a/math.sqrt(a**2+b**2))))*100)/100)
				else:
					angle.append(0)

######################################AZIMUTS 315 => 360######################################

		if azimut > 315 and azimut < 360:
		#calculate corresponding angle of right-angled triangle
			alpha = math.radians(360-azimut)
			cos_alpha = math.cos(alpha)
			sin_alpha = math.sin(alpha)
			imax = int(math.floor(viewmax*cos_alpha/25))
			angle = []
			for i in range(1,imax):
				myview = myview+25
				c = myview/cos_alpha
			#incrementation depending on azimut value (take care of x,y and row, col order and sign for ptx (or pty) calcul
				ptrow = sta_rc[0]-i
				ptx = rastinit[0]+((math.floor((sta_xy[0]-(sin_alpha*c)-rastinit[0])/rastinit[1]))*rastinit[1])
				ptcol = numpy.unique(numpy.argwhere(imgx==ptx))[1]-mincol
				ptrc=(ptrow, ptcol)
				
				#calculate corresponding angle to reach the corresponding height
				b = height[ptrc]-height[sta_rc]
				if b > 0:
					#take care of row col and order between viewpt and sta
					a = ((sta_rc[0]-ptrc[0])*25)/cos_alpha
					angle.append(math.ceil((math.degrees(math.acos(a/math.sqrt(a**2+b**2))))*100)/100)
				else:
					angle.append(0)
		
		#append each azimut to final data for weather station
		data = (sta[0], azimut, max(angle))
		final_data.append(data)
		#print(data)

	final_data.append((final_data[0][0], 360, final_data[0][2]))

	#insert values into new table for the given weather station
	for values in final_data:
		skylinetable=myconn.cursor()
		query="""
		insert into stations.station_meteo_skyline
		values(%s, %s, %s);
		"""

		skylinetable.execute(query,(values[0],values[1], values[2]))
		myconn.commit()
	
	final_data = None
	print sta[0], "done"

print"done in", time.time()-start_time, "seconds"