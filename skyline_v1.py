import numpy
import math
#import matplotlib.pyplot as plt
#import matplotlib.image
import gdal
#import osr
#import ogr
import psycopg2
#import os
#import sys

# load complete raster
img = gdal.Open('C:\ds_test_data\ign_mnt25_alpes.tif')
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
myconn = psycopg2.connect("host = psql9.grenoble.cemagref.fr dbname = xxxx user=xxxx password=xxxx")

#load points and view extent
wsta=myconn.cursor()
query="""
with a as (
	select gid, st_x(the_geom) x, st_y(the_geom) y, st_buffer(the_geom,%s) geom
	from stations.geo_station_meteofrance
	)

select gid, x, y, st_xmin(geom) xmin, st_ymin(geom), st_xmax(geom) xmax, st_ymax(geom) ymax
from a
where gid = '73024400'
order by gid
;
"""

viewmax=1000
wsta.execute(query,(viewmax,))

final_data = []

#extract from original raster

for sta in wsta:
	print sta[0]
	minrow = None
	maxrow = None
	mincol = None
	maxcol = None
	
	for i in range(0,imgx.shape[1]-1):
		if i+1 <= imgx.shape[1] and sta[3]-imgx[0,i]>0 and sta[3]-imgx[0,i+1]<=0:
			mincol = i
		if i+1 <= imgx.shape[1] and sta[5]-imgx[0,i]>0 and sta[5]-imgx[0,i+1]<=0:
			maxcol = i
	for i in range(0,imgx.shape[1]-1):
		if sta[1]-imgx[0,i]>0 and sta[1]-imgx[0,i+1]<=0:
			stacol = i-mincol
			stax = imgx[0,i]+((imgx[0,i+1]-imgx[0,i])/2)
	if minrow is None:
		minrow = 0
	if maxrow is None:
		maxrow = imgx.shape[1]
		
	for i in range(0,imgy.shape[0]-1):
		if i+1 <= imgy.shape[0] and imgy[i,0]-sta[6]>0 and imgy[i+1,0]-sta[6]<=0:
			minrow = i
		if i+1 <= imgy.shape[0] and imgy[i,0]-sta[4]>0 and imgy[i+1,0]-sta[4]<=0:
			maxrow = i
	for i in range(0,imgy.shape[0]-1):
		if imgy[i,0]-sta[2]>0 and imgy[i+1,0]-sta[2]<=0:
			starow = i-minrow
			stay = imgy[i+1,0]+((imgy[i,0]-imgy[i+1,0])/2)
	if mincol is None:
		mincol = 0
	if maxcol is None:
		maxcol = imgy.shape[0]
			
	sta_xy = (stax, stay)
	sta_rc = (starow, stacol)

	height = band1.ReadAsArray(mincol, minrow, maxcol-mincol, maxrow-minrow)

	# get width and heigth of image
	w,h = height.shape
	
	print "raster extracted"
	print w, h
	
	
	#Get all intersected cells on azimuth
	for azimut in range (0, 360, 5):
		print azimut
		myview=0
		c= None
		#TESTS on azimut
		
######################################AZIMUT 0 / 360######################################

		if azimut == 0:
			imax = viewmax/25
			viewpt_rc = []
			for i in range(1,imax):
			#incrementation depending on azimut value (take care of x,y and row, col order and sign for ptx (or pty) calcul
				ptrow = sta_rc[0]-i
				ptcol = sta_rc[1]
				ptrc=(ptrow, ptcol)
				viewpt_rc.append(ptrc)

			print "cells collected"
			
			pth = []
			for p in viewpt_rc:
				pth.append(height[p])
			maxh = max(pth)
			myindex = pth.index(maxh)
			
			#calculate distance of line from ref point to max alt point: take care of row col and order between viewpt and sta
			d = ((sta_rc[0]-viewpt_rc[myindex][0])*25)
			
######################################AZIMUT 180######################################

		if azimut == 180:
			imax = viewmax/25
			viewpt_rc = []
			for i in range(1,imax):
			#incrementation depending on azimut value (take care of x,y and row, col order and sign for ptx (or pty) calcul
				ptrow = sta_rc[0]+i
				ptcol = sta_rc[1]
				ptrc=(ptrow, ptcol)
				viewpt_rc.append(ptrc)

			print "cells collected"
			
			pth = []
			for p in viewpt_rc:
				pth.append(height[p])
			maxh = max(pth)
			myindex = pth.index(maxh)
			
			#calculate distance of line from ref point to max alt point: take care of row col and order between viewpt and sta
			d = ((viewpt_rc[myindex][0]-sta_rc[0])*25)
			
######################################AZIMUT 90######################################

		if azimut == 90:
			viewpt_rc = []
			for i in range(1,imax):
				#incrementation depending on azimut value (take care of x,y and row, col order
				ptcol = sta_rc[1]+i
				ptrow=sta_rc[0]
				ptrc=(ptrow, ptcol)
				viewpt_rc.append(ptrc)

			print "cells collected"
			
			pth = []
			for p in viewpt_rc:
				pth.append(height[p])
			maxh = max(pth)
			myindex = pth.index(maxh)
			
			#calculate distance of line from ref point to max alt point: take care of row col
			d = ((viewpt_rc[myindex][1]-sta_rc[1])*25)

######################################AZIMUT 270######################################

		if azimut == 270:
			viewpt_rc = []
			for i in range(1,imax):
				#incrementation depending on azimut value (take care of x,y and row, col order
				ptcol = sta_rc[1]-i
				ptrow=sta_rc[0]
				ptrc=(ptrow, ptcol)
				viewpt_rc.append(ptrc)

			print "cells collected"
			
			pth = []
			for p in viewpt_rc:
				pth.append(height[p])
			maxh = max(pth)
			myindex = pth.index(maxh)
			
			#calculate distance of line from ref point to max alt point: take care of row col
			d = ((sta_rc[1]-viewpt_rc[myindex][1])*25)

######################################AZIMUT 45######################################

		if azimut == 45:
			#calculate corresponding angle of right-angled triangle
			alpha = 45
			imax = int(math.floor(viewmax*math.cos(math.radians(alpha))/25))
			viewpt_rc = []
			for i in range(1,imax):
				#incrementation depending on azimut value (take care of x,y and row, col order
				ptcol = sta_rc[1]+i
				ptrow = sta_rc[0]-i
				ptrc=(ptrow, ptcol)
				viewpt_rc.append(ptrc)

			print "cells collected"
			
			pth = []
			for p in viewpt_rc:
				pth.append(height[p])
			maxh = max(pth)
			myindex = pth.index(maxh)
			
			#calculate distance of line from ref point to max alt point: take care of row col
			a = ((viewpt_rc[myindex][1]-sta_rc[1])*25)
			d = a/math.cos(math.radians(alpha))

######################################AZIMUT 135######################################

		if azimut == 135:
			#calculate corresponding angle of right-angled triangle
			alpha = 45
			imax = int(math.floor(viewmax*math.cos(math.radians(alpha))/25))
			viewpt_rc = []
			for i in range(1,imax):
				#incrementation depending on azimut value (take care of x,y and row, col order
				ptcol = sta_rc[1]+i
				ptrow = sta_rc[0]+i
				ptrc=(ptrow, ptcol)
				viewpt_rc.append(ptrc)

			print "cells collected"
			
			pth = []
			for p in viewpt_rc:
				pth.append(height[p])
			maxh = max(pth)
			myindex = pth.index(maxh)
			
			#calculate distance of line from ref point to max alt point: take care of row col
			a = ((viewpt_rc[myindex][1]-sta_rc[1])*25)
			d = a/math.cos(math.radians(alpha))

######################################AZIMUT 225######################################

		if azimut == 225:
			#calculate corresponding angle of right-angled triangle
			alpha = 45
			imax = int(math.floor(viewmax*math.cos(math.radians(alpha))/25))
			viewpt_rc = []
			for i in range(1,imax):
				#incrementation depending on azimut value (take care of x,y and row, col order
				ptcol = sta_rc[1]-i
				ptrow = sta_rc[0]+i
				ptrc=(ptrow, ptcol)
				viewpt_rc.append(ptrc)

			print "cells collected"
			
			pth = []
			for p in viewpt_rc:
				pth.append(height[p])
			maxh = max(pth)
			myindex = pth.index(maxh)
			
			#calculate distance of line from ref point to max alt point: take care of row col
			a = ((viewpt_rc[myindex][0]-sta_rc[0])*25)
			d = a/math.cos(math.radians(alpha))

######################################AZIMUT 315######################################

		if azimut == 315:
			#calculate corresponding angle of right-angled triangle
			alpha = 45
			imax = int(math.floor(viewmax*math.cos(math.radians(alpha))/25))
			viewpt_rc = []
			for i in range(1,imax):
				#incrementation depending on azimut value (take care of x,y and row, col order
				ptcol = sta_rc[1]-i
				ptrow = sta_rc[0]-i
				ptrc=(ptrow, ptcol)
				viewpt_rc.append(ptrc)

			print "cells collected"
			
			pth = []
			for p in viewpt_rc:
				pth.append(height[p])
			maxh = max(pth)
			myindex = pth.index(maxh)
			
			#calculate distance of line from ref point to max alt point: take care of row col
			a = ((sta_rc[0]-viewpt_rc[myindex][0])*25)
			d = a/math.cos(math.radians(alpha))

######################################AZIMUTS 0 => 45######################################

		if azimut > 0 and azimut < 45:
		#calculate corresponding angle of right-angled triangle
			alpha = azimut
			imax = int(math.floor(viewmax*math.cos(math.radians(alpha))/25))
			viewpt_rc = []
			for i in range(1,imax):
				myview = myview+25
				c = myview/math.cos(math.radians(alpha))
			#incrementation depending on azimut value (take care of x,y and row, col order and sign for ptx (or pty) calcul
				ptrow = sta_rc[0]-i
				ptx = sta_xy[0]+(math.sin(math.radians(alpha))*c)
				for i in range(0,imgx.shape[1]-1):
					if ptx-imgx[0,i]>0 and ptx-imgx[0,i+1]<=0:
						ptcol = i-mincol
				ptrc=(ptrow, ptcol)
				viewpt_rc.append(ptrc)

			print "cells collected"
			
			pth = []
			for p in viewpt_rc:
				pth.append(height[p])
			maxh = max(pth)
			myindex = pth.index(maxh)
			
			#calculate distance of line from ref point to max alt point: take care of row col and order between viewpt and sta
			a = ((sta_rc[0]-viewpt_rc[myindex][0])*25)
			d = a/math.cos(math.radians(alpha))
			
######################################AZIMUTS 45 => 90######################################

		if azimut > 45 and azimut < 90:
			#calculate corresponding angle of right-angled triangle
			alpha = 90-azimut
			imax = int(math.floor(viewmax*math.cos(math.radians(alpha))/25))
			viewpt_rc = []
			for i in range(1,imax):
				myview = myview+25
				c = myview/math.cos(math.radians(alpha))
				#incrementation depending on azimut value (take care of x,y and row, col order
				ptcol = sta_rc[1]+i
				pty = sta_xy[1]+(math.sin(math.radians(alpha))*c)
				for i in range(0,imgy.shape[0]-1):
					if imgy[i,0]-pty>0 and imgy[i+1,0]-pty<=0:
						ptrow = i-minrow
				ptrc=(ptrow, ptcol)
				viewpt_rc.append(ptrc)

			print "cells collected"
			
			pth = []
			for p in viewpt_rc:
				pth.append(height[p])
			maxh = max(pth)
			myindex = pth.index(maxh)
			
			#calculate distance of line from ref point to max alt point: take care of row col
			a = ((viewpt_rc[myindex][1]-sta_rc[1])*25)
			d = a/math.cos(math.radians(alpha))
			
######################################AZIMUTS 90 => 135######################################

		if azimut > 90 and azimut < 135:
			#calculate corresponding angle of right-angled triangle
			alpha = azimut-90
			imax = int(math.floor(viewmax*math.cos(math.radians(alpha))/25))
			viewpt_rc = []
			for i in range(1,imax):
				myview = myview+25
				c = myview/math.cos(math.radians(alpha))
				#incrementation depending on azimut value (take care of x,y and row, col order
				ptcol = sta_rc[1]+i
				pty = sta_xy[1]-(math.sin(math.radians(alpha))*c)
				for i in range(0,imgy.shape[0]-1):
					if imgy[i,0]-pty>0 and imgy[i+1,0]-pty<=0:
						ptrow = i-minrow
				ptrc=(ptrow, ptcol)
				viewpt_rc.append(ptrc)

			print "cells collected"
			
			pth = []
			for p in viewpt_rc:
				pth.append(height[p])
			maxh = max(pth)
			myindex = pth.index(maxh)
			
			#calculate distance of line from ref point to max alt point: take care of row col and order between viewpt and sta
			a = ((viewpt_rc[myindex][1]-sta_rc[1])*25)
			d = a/math.cos(math.radians(alpha))
			
######################################AZIMUTS 135 => 180######################################

		if azimut > 135 and azimut < 180:
			#calculate corresponding angle of right-angled triangle
			alpha = 180-azimut
			imax = int(math.floor(viewmax*math.cos(math.radians(alpha))/25))
			viewpt_rc = []
			for i in range(1,imax):
				myview = myview+25
				c = myview/math.cos(math.radians(alpha))
				#incrementation depending on azimut value (take care of x,y and row, col order
				ptrow = sta_rc[0]+i
				ptx = sta_xy[0]+(math.sin(math.radians(alpha))*c)
				for i in range(0,imgx.shape[1]-1):
					if ptx-imgx[0,i]>0 and ptx-imgx[0,i+1]<=0:
						ptcol = i-mincol
				ptrc=(ptrow, ptcol)
				viewpt_rc.append(ptrc)

			print "cells collected"
			
			pth = []
			for p in viewpt_rc:
				pth.append(height[p])
			maxh = max(pth)
			myindex = pth.index(maxh)
			
			#calculate distance of line from ref point to max alt point: take care of row col and order between viewpt and sta
			a = ((viewpt_rc[myindex][0]-sta_rc[0])*25)
			d = a/math.cos(math.radians(alpha))
			
######################################AZIMUTS 180 => 225######################################

		if azimut > 180 and azimut < 225:
		#calculate corresponding angle of right-angled triangle
			alpha = azimut - 180
			imax = int(math.floor(viewmax*math.cos(math.radians(alpha))/25))
			viewpt_rc = []
			for i in range(1,imax):
				myview = myview+25
				c = myview/math.cos(math.radians(alpha))
			#incrementation depending on azimut value (take care of x,y and row, col order and sign for ptx (or pty) calcul
				ptrow = sta_rc[0]+i
				ptx = sta_xy[0]-(math.sin(math.radians(alpha))*c)
				for i in range(0,imgx.shape[1]-1):
					if ptx-imgx[0,i]>0 and ptx-imgx[0,i+1]<=0:
						ptcol = i-mincol
				ptrc=(ptrow, ptcol)
				viewpt_rc.append(ptrc)

			print "cells collected"
			
			pth = []
			for p in viewpt_rc:
				pth.append(height[p])
			maxh = max(pth)
			myindex = pth.index(maxh)
			
			#calculate distance of line from ref point to max alt point: take care of row col and order between viewpt and sta
			a = ((viewpt_rc[myindex][0]-sta_rc[0])*25)
			d = a/math.cos(math.radians(alpha))
			

######################################AZIMUTS 225 => 270######################################

		if azimut > 225 and azimut < 270:
			#calculate corresponding angle of right-angled triangle
			alpha = 270 - azimut
			imax = int(math.floor(viewmax*math.cos(math.radians(alpha))/25))
			viewpt_rc = []
			for i in range(1,imax):
				myview = myview+25
				c = myview/math.cos(math.radians(alpha))
				#incrementation depending on azimut value (take care of x,y and row, col order
				ptx = sta_xy[0]-25
				ptcol = sta_rc[1]-i
				pty = sta_xy[1]-(math.sin(math.radians(alpha))*c)
				for i in range(0,imgy.shape[0]-1):
					if imgy[i,0]-pty>0 and imgy[i+1,0]-pty<=0:
						ptrow = i-minrow
				ptrc=(ptrow, ptcol)
				viewpt_rc.append(ptrc)

			print "cells collected"
			
			pth = []
			for p in viewpt_rc:
				pth.append(height[p])
			maxh = max(pth)
			myindex = pth.index(maxh)
			
			#calculate distance of line from ref point to max alt point: take care of row col and order between viewpt and sta
			a = ((sta_rc[1]-viewpt_rc[myindex][1])*25)
			d = a/math.cos(math.radians(alpha))
			
######################################AZIMUTS 270 => 315######################################

		if azimut > 270 and azimut < 315:
			#calculate corresponding angle of right-angled triangle
			alpha = azimut-270
			imax = int(math.floor(viewmax*math.cos(math.radians(alpha))/25))
			viewpt_rc = []
			for i in range(1,imax):
				myview = myview+25
				c = myview/math.cos(math.radians(alpha))
				#incrementation depending on azimut value (take care of x,y and row, col order
				ptx = sta_xy[0]-25
				ptcol = sta_rc[1]-i
				pty = sta_xy[1]+(math.sin(math.radians(alpha))*c)
				for i in range(0,imgy.shape[0]-1):
					if imgy[i,0]-pty>0 and imgy[i+1,0]-pty<=0:
						ptrow = i-minrow
				ptrc=(ptrow, ptcol)
				viewpt_rc.append(ptrc)

			print "cells collected"
			
			pth = []
			for p in viewpt_rc:
				pth.append(height[p])
			maxh = max(pth)
			myindex = pth.index(maxh)
			
			#calculate distance of line from ref point to max alt point: take care of row col and order between viewpt and sta
			a = ((sta_rc[1]-viewpt_rc[myindex][1])*25)
			d = a/math.cos(math.radians(alpha))

######################################AZIMUTS 315 => 360######################################

		if azimut > 315 and azimut < 360:
		#calculate corresponding angle of right-angled triangle
			alpha = 360-azimut
			imax = int(math.floor(viewmax*math.cos(math.radians(alpha))/25))
			viewpt_rc = []
			for i in range(1,imax):
				myview = myview+25
				c = myview/math.cos(math.radians(alpha))
			#incrementation depending on azimut value (take care of x,y and row, col order and sign for ptx (or pty) calcul
				ptrow = sta_rc[0]-i
				ptx = sta_xy[0]-(math.sin(math.radians(alpha))*c)
				for i in range(0,imgx.shape[1]-1):
					if ptx-imgx[0,i]>0 and ptx-imgx[0,i+1]<=0:
						ptcol = i-mincol
				ptrc=(ptrow, ptcol)
				viewpt_rc.append(ptrc)

			print "cells collected"
			
			pth = []
			for p in viewpt_rc:
				pth.append(height[p])
			maxh = max(pth)
			myindex = pth.index(maxh)
			
			#calculate distance of line from ref point to max alt point: take care of row col and order between viewpt and sta
			a = ((sta_rc[0]-viewpt_rc[myindex][0])*25)
			d = a/math.cos(math.radians(alpha))
			
		#TEST ANGLE TO SEE UPPER THAN UPPER ELEV
		for beta in range(0,90, 2):
			mybeta=beta
			h = (d*math.sin(math.radians(beta))/math.cos(math.radians(beta)))+height[sta_rc]
			if h > maxh:
				break
		data = (sta[0], azimut, mybeta)
		final_data.append(data)
		print(data)

final_data.append((final_data[0][0], 360, final_data[0][2]))

#create table in postgres

#load points and view extent
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

#insert values into new table
for values in final_data:
	skylinetable=myconn.cursor()
	query="""
	insert into stations.station_meteo_skyline
	values(%s, %s, %s);
	"""

	skylinetable.execute(query,(values[0],values[1], values[2]))
	myconn.commit()
	print "azimut", values[1], "done"

print(final_data)