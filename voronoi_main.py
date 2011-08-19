#! /usr/local/bin/python

import glob
import sys
import json
import pylab
import hashlib
from shapely.geometry import MultiPoint
from pylab import *
import voronoi_gen

PlotIt=True

#files=glob.glob('data_cities_lv3/c_*_*_*.ctr')
#files=glob.glob('data_cities_lv2/c_*_*.ctr')
files=glob.glob('data_cities_lv1/c_*.ctr')

files.sort()
PointsMap={}
for f in files:
  #print f
  
  cluster_id=f.split(".")[0]
  cluster_id=cluster_id.split("/")[1]

  #print cluster_id

  center_pnt={}
  for line in open(f):
    data=line.strip().replace("\"", "").split(" ")
    if len(data)==2:
      center_pnt[data[0]]=data[1]
  
  #print center_pnt["x"], center_pnt["y"]

  PointsMap[cluster_id]=[ center_pnt["x"], center_pnt["y"] ]

vl=voronoi_gen.GenerateVoronoi(cluster_id, PointsMap, PlotIt)      
#print vl

for vl_num, data in vl.items():

  #print data

  cluster_id=data["info"]
  geojson={}

  geojson["properties"]={}
  geojson["properties"]["_namespace"]="geo.location.foursquare"
  geojson["properties"]["_protocol"]="places"
  geojson["properties"]["_type"]="clusters"
  geojson["properties"]["_domain_id"]=cluster_id

  geojson["type"]="Feature"
  

  m = hashlib.md5()
  hash_string="geo.location.foursquare.places.clusters:"+cluster_id
  m.update(hash_string)
  geojson["id"]=m.hexdigest()

  

  
  #geojson["properties"]["geo_center"]={}
  
  #geojson["properties"]["geo_center"]["type"]="Feature"
  #geojson["properties"]["geo_center"]["geometry"]={}
  #geojson["properties"]["geo_center"]["geometry"]["type"]="Point"
  #geojson["properties"]["geo_center"]["geometry"]["coordinates"]=data["coordinate"]
  

  geojson["geometry"]={}
  geojson["geometry"]["type"]="Polygon"
  

  coords=[]
  geojson["geometry"]["coordinates"]=[list(data["obj_polygon"].exterior.coords)]

  x=[]
  y=[]
  for point in list(data["obj_polygon"].exterior.coords):
    x.append(point[0])
    y.append(point[1])
  #fill(x,y, alpha=0.6)
  #show()

  print json.dumps(geojson)#, sort_keys=True, indent=3)
  #break
  #print str_spaced.replace(" ", "")






