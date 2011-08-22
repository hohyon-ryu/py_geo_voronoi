#! /usr/local/bin/python

import voronoi
import sys
import getopt
from pylab import *
#import matplotlib
from shapely.geometry import Polygon
from shapely.ops import polygonize
import hashlib
import json

import time

import globalmaptiles

WorldRanges={}
WorldRanges["AUSTIN"]=[30.8, -98.5, 29.535, -97.031]
WorldRanges["TX"]=[36.5, -106, 25, -93]
WorldRanges["US"]=[55, -130, 23, -60]
WorldRanges["GUS"]=[60, -140, 22, -50]
WorldRanges["KR"]=[45, 120, 32, 135]
WorldRanges["W"]=[90, -180, -90, 180]

World="W"
WorldRange=WorldRanges[World]

PlotIt=False

def checkInRange(WorldRange,x,y):

  if x>=WorldRange[1] and x<=WorldRange[3]:
    if y<=WorldRange[0] and y>=WorldRange[2]:
      return True
  return False

def getExtremes(line, m_range):
  a,b,c=line

  if a==0:
    y0=(c-(a*m_range["min_x"]))/b
    y1=(c-(a*m_range["max_x"]))/b
    return [(m_range["min_x"],y0), (m_range["max_x"], y1)]

  if b==0:
    x0=(c-(b*m_range["min_y"]))/a
    x1=(c-(b*m_range["max_y"]))/a
    return [(x0, m_range["min_y"]), (x1, m_range["max_y"])]  
  
  x0=(c-(b*m_range["min_y"]))/a
  x1=(c-(b*m_range["max_y"]))/a

  y0=(c-(a*m_range["min_x"]))/b
  y1=(c-(a*m_range["max_x"]))/b

  return [(x0, m_range["min_y"]), (x1, m_range["max_y"]), (m_range["min_x"],y0), (m_range["max_x"], y1)]


def getExtreme(line, v_known, LR):
  
  global WorldRange

  global_extremes=getExtremes(line, {"min_y":WorldRange[2], "max_y":WorldRange[0], "min_x":WorldRange[1], "max_x":WorldRange[3]})

  v_unknown=None
  CurrentExtreme=None

  for extreme in global_extremes:

    if checkInRange(WorldRange,extreme[0],extreme[1])==False:
      #v_unknown=CurrentExtreme
      continue

    if v_known[0]<extreme[0] and LR==1:
      v_unknown=extreme

    if v_known[0]>extreme[0] and LR==0:
      v_unknown=extreme
    
    if CurrentExtreme:
      if CurrentExtreme[0]==extreme[0]:
        if CurrentExtreme:      

          if CurrentExtreme[1]<extreme[1] and LR==0:
            v_unknown=extreme
          else:
            v_unknown=CurrentExtreme
      
    CurrentExtreme=extreme

    #print LR, v_known, extreme, v_unknown

  return v_unknown    


def update_maxmin(m_range, x, y):
  if (m_range["max_x"]<x):
    m_range["max_x"]=x
  if (m_range["max_y"]<y):
    m_range["max_y"]=y
  if (m_range["min_x"]>x):
    m_range["min_x"]=x
  if (m_range["min_y"]>y):
    m_range["min_y"]=y
  return m_range

def linkExtremes(point1, point2, m_range):
  #print "Link Extremes"
  #print point1, point2
  global WorldRange

  if point1[1]==point2[1] or point1[0]==point2[0]:
    #print "already done case"
    return [(point1, point2)]

  if point1[1]==-point2[1] or point1[0]==-point2[0]:
    #print "worst case"
    if point1[1]==-point2[1]:
      if abs(point1[0]-m_range["min_x"])+abs(point2[0]-m_range["min_x"])<abs(point1[0]-m_range["max_x"])+abs(point2[0]-m_range["max_x"]):
        output=[(point1, (WorldRange[1], point1[1])), ((WorldRange[1], point1[1]),(WorldRange[1], point2[1])), ((WorldRange[1], point2[1]), point2)]
      else:
        output=[(point1, (WorldRange[3], point1[1])), ((WorldRange[3], point1[1]),(WorldRange[3], point2[1])), ((WorldRange[3], point2[1]), point2)]
    else:
      if abs(point1[1]-m_range["min_y"])+abs(point2[1]-m_range["min_y"])<abs(point1[1]-m_range["max_y"])+abs(point2[1]-m_range["max_y"]):
        output=[(point1, (point1[0], WorldRange[2])), ((point1[0], WorldRange[2]),(point2[0], WorldRange[2])), ((point2[0], WorldRange[2]), point2)]
      else:
        output=[(point1, (point1[0], WorldRange[0])), ((point1[0], WorldRange[0]),(point2[0], WorldRange[0])), ((point2[0], WorldRange[0]), point2)]

    return output

  #print "one corner case"
  if point1[1]==WorldRange[0] or point1[1]==WorldRange[2]:
    #print "add x:", point2[0]
    #print "add y:", point1[1]
    output=[(point1, (point2[0], point1[1])), ((point2[0], point1[1]), point2)]
  else:
    #print "add x:", point1[0]
    #print "add y:", point2[1]
    output=[(point1, (point1[0], point2[1])), ((point1[0], point2[1]), point2)]
  return output

def VoronoiLineEdges(PointsMap):
  Sitepts = []
  pts = {}

  #print CurrentDate, PointsMap[PointsMap.keys()[0]]
  for grid, stn in PointsMap.items():

    x=float(stn[0])
    y=float(stn[1])
    station=grid
    #station.extend( stn[3:])
    #print x,y,station

    pts[ (x,y) ]=station
  
  stncounts=len(pts.keys())
  #print stncounts, "points"
  
  site_points=[]
  for pt in pts.keys():
    Sitepts.append(voronoi.Site(pt[0],pt[1]))
    site_points.append( (pt[0],pt[1]) )
      

  #print "Calculating Voronoi Lattice",
  
  siteList = voronoi.SiteList(Sitepts)
  context  = voronoi.Context()
  voronoi.Edge.EDGE_NUM=0  
  voronoi.voronoi(siteList,context)

  vertices=context.vertices
  lines=context.lines
  edges=context.edges
  triangles=context.triangles
  has_edge=context.has_edge

  return vertices, lines, edges, has_edge


def VoronoiGeoJson_MultiPolygons(PointsMap, BoundingBox="W",PlotMap=False):

  vl=VoronoiPolygons(PointsMap, BoundingBox="W", PlotMap=PlotMap)      

  cluster_id="VoronoiLattice"
  geojson={}

  geojson["properties"]={}
  geojson["properties"]["_domain_id"]=cluster_id

  geojson["type"]="Feature"
  

  m = hashlib.md5()
  hash_string=cluster_id
  m.update(hash_string)
  geojson["id"]=m.hexdigest()

  geojson["geometry"]={}
  geojson["geometry"]["type"]="MultiPolygon"

  geojson["geometry"]["coordinates"]=[]

  for vl_num, data in vl.items():

    geojson["geometry"]["coordinates"].append([list(data["obj_polygon"].exterior.coords)])


  return json.dumps(geojson)


def VoronoiGeoJson_Polygons(PointsMap, BoundingBox="W",PlotMap=False):

  vl=VoronoiPolygons(PointsMap, BoundingBox="W", PlotMap=PlotMap)      

  output=""
  for vl_num, data in vl.items():

    #print data

    cluster_id=data["info"]
    geojson={}

    geojson["properties"]={}
    geojson["properties"]["_domain_id"]=cluster_id

    geojson["type"]="Feature"
    

    m = hashlib.md5()
    hash_string=cluster_id
    m.update(hash_string)
    geojson["id"]=m.hexdigest()

    geojson["geometry"]={}
    geojson["geometry"]["type"]="Polygon"

    coords=[]
    geojson["geometry"]["coordinates"]=[list(data["obj_polygon"].exterior.coords)]

    output+=json.dumps(geojson)#, sort_keys=True, indent=3)
    output+="\n"
  return output

def VoronoiPolygons(PointsMap, BoundingBox="W", PlotMap=False):
  global PlotIt
  global WorldRange
  global WorldRanges

  if type(BoundingBox)==type([]):
    if len(BoundingBox)==4:
      WorldRange=BoundingBox
    else:
      return "Error in Bounding Box"
  else:
    WorldRange=WorldRanges[BoundingBox]
    
  PlotIt=PlotMap

  currenttime=time.time()
  Sitepts = []
  pts = {}

  #print CurrentDate, PointsMap[PointsMap.keys()[0]]
  for grid, stn in PointsMap.items():

    x=float(stn[0])
    y=float(stn[1])
    station=grid
    #station.extend( stn[3:])
    #print x,y,station

    pts[ (x,y) ]=station
  
  stncounts=len(pts.keys())
  #print stncounts, "points"
  
  site_points=[]
  for pt in pts.keys():
    Sitepts.append(voronoi.Site(pt[0],pt[1]))
    site_points.append( (pt[0],pt[1]) )
      

  #print "Calculating Voronoi Lattice",
  
  siteList = voronoi.SiteList(Sitepts)
  context  = voronoi.Context()
  voronoi.Edge.EDGE_NUM=0   
  voronoi.voronoi(siteList,context)

  vertices=context.vertices
  lines=context.lines
  edges=context.edges
  
  #print edges

  #For Faster Access
  edge_dic={}
  for edge in edges:
    edge_dic[edge[0]]=edge[1:]
  
  triangles=context.triangles
  has_edge=context.has_edge

  voronoi_lattice={}

  m_range={}
  m_range["max_x"]=-9999999999
  m_range["max_y"]=-9999999999
  m_range["min_x"]=9999999999
  m_range["min_y"]=9999999999

  #Get the range!!
  for pnt in site_points:
    m_range=update_maxmin(m_range, pnt[0], pnt[1])

  #print "Getting the Polygons"

  prev_percent=0
  for station, ls in has_edge.items():


    voronoi_lattice[station]={}
    voronoi_lattice[station]["coordinate"]=site_points[station]
    voronoi_lattice[station]["info"]=pts[ site_points[station] ]

    polygon=[]
 
    prev_extreme=[]
    Verbose=True
    if Verbose: 
      current_percent=int(station/float(stncounts)*100)
      if current_percent!=prev_percent:
        #print station,"/", stncounts, current_percent, "% Done" 
        timeelapse=time.time()-currenttime
        #print station, timeelapse
        currenttime=time.time()

      prev_percent=current_percent

    #For every lines that the station owns
    for l in ls:
      e=edge_dic[l]

      v1=vertices[e[0]]
      v2=vertices[e[1]]
    
      if e[0] < 0 and checkInRange(WorldRange,v2[0],v2[1])==False: continue
      if e[1] < 0 and checkInRange(WorldRange,v1[0],v1[1])==False: continue

      if e[0] > -1 and e[1] > -1 and checkInRange(WorldRange,v1[0],v1[1])==False and checkInRange(WorldRange,v2[0],v2[1])==False:
        continue 

      if e[0] < 0 or checkInRange(WorldRange,v1[0],v1[1])==False:
        v1=getExtreme(lines[l],v2, LR=0)

        if len(prev_extreme)==0:
          prev_extreme=v1
        else:
          extreme_points=linkExtremes(prev_extreme, v1, m_range)
          for extreme_pair in extreme_points:
            polygon.append(extreme_pair)

      if e[1] < 0 or checkInRange(WorldRange,v2[0],v2[1])==False :
        v2=getExtreme(lines[l],v1, LR=1)

        if len(prev_extreme)==0:
          prev_extreme=v2
        else:
          extreme_points=linkExtremes(prev_extreme, v2, m_range)
          for extreme_pair in extreme_points:
            polygon.append(extreme_pair)

      if v1!=v2:
        polygon.append( (v1,v2) )

    if len(polygon)==0:
      sys.stderr.write ("\nThis station does not have meaningful polygon:")
      sys.stderr.write (str(pts[ site_points[station] ])+" at "+str(site_points[station])+"\n")
      for l in ls:
        e=edge_dic[l]

        v1=vertices[e[0]]
        v2=vertices[e[1]]
        
        sys.stderr.write(str(e[0])+","+str(e[1])+"\n")
        sys.stderr.write(str(v1)+","+str(v2)+"\n")

      voronoi_lattice.pop(station)
      continue
      

    #print polygon
    try:
      result = list(polygonize( polygon ))
    except ValueError:
      #print "Wrong:", polygon

      #Delete invisible short lines
      point_to_replace=[]
      for line in polygon:
        d=math.hypot( line[0][0]-line[1][0], line[0][1]-line[1][1])
        if d<1:
          polygon.remove(line)
          point_to_replace=line
      i=0
      New_Lines=[]
      for line in polygon:
        New_Lines.append(list(line))
        j=0
        for point in line:
          if point==point_to_replace[0]:
            #print line, point_to_replace
            New_Lines[i][j]=point_to_replace[1]
          j+=1
        i+=1

      polygon=tuple(New_Lines)
      
      #I could not figure out why sometimes it fails to draw a polygon. 
      try:
        result = list(polygonize( polygon ))
      except:
        voronoi_lattice.pop(station)
        continue





    finalpoly=result[0]
    #print list(finalpoly.exterior.coords)
    #voronoi_lattice[station]["polygon"]=list(finalpoly.exterior.coords)
    voronoi_lattice[station]["obj_polygon"]=finalpoly
    
    #print polygon
    


  if (PlotIt):
    for station, data in voronoi_lattice.items():

      x=[]
      y=[]
      #print data
      try:
        polygon_data=data["obj_polygon"]
        #print polygon_data
        for point in list(polygon_data.exterior.coords):
          x.append(point[0])
          y.append(point[1])
        
      except:
        print "Error", data["name"]    
      #if station == 8 :
      #plot(x,y)
      
      fill(x,y, alpha=0.6)
      
      plot(data["coordinate"][0],data["coordinate"][1])
      #text(data["coordinate"][0],data["coordinate"][1],data["info"])

    show()
    
  return voronoi_lattice



ClosestZoomLevel=7



def plot_voronoi(voronoi_lattice):
  for station, data in voronoi_lattice.items():

    x=[]
    y=[]
    polygon_data=data["obj_polygon"]
    #print polygon_data
    for point in list(polygon_data.exterior.coords):
      x.append(point[0])
      y.append(point[1])
    
  #if station == 8 :
  #plot(x,y)
    fill(x,y, alpha=0.3)

  #scatter(data["coordinate"][0], data["coordinate"][1])
#  show()


def quadGrid(grid_range):
  #print grid_range
  bix=(grid_range[0]+grid_range[2])/2
  biy=(grid_range[1]+grid_range[3])/2
  
  grid0=[grid_range[0], grid_range[1], bix, biy]
  grid1=[grid_range[0], biy, bix, grid_range[3]]
  grid2=[bix, grid_range[1], grid_range[2], biy]
  grid3=[bix, biy, grid_range[2], grid_range[3]]

  return (grid0, grid1, grid2, grid3)

def polygonize_grid(grid_range):
  x0=grid_range[1]
  y0=grid_range[0]
  x1=grid_range[3]
  y1=grid_range[2]

  polygon=( (x0, y0), (x0, y1), (x1, y1), (x1, y0), (x0, y0))
  return polygon

def get_quadkeystr(quadkey):
  for i in range(ClosestZoomLevel-len(quadkey)):
    quadkey+="X"
  return quadkey
#Zoom in!!
def GridMap(quadkey="", ploygons_to_lookat=None):

  mercator = globalmaptiles.GlobalMercator()
 
  global PlotIt

  for quadkeyadd in range(4):
    curQuadKey=quadkey+str(quadkeyadd)

    tx, ty, zl=mercator.QuadTree2TMS(curQuadKey)
    grid_latlon=mercator.TileLatLonBounds(tx, ty, zl)    
    #print grid_latlon

    polygon=polygonize_grid(grid_latlon) 

    #print polygon

    obj_polygon=Polygon(polygon)

    inter_cnt=0
    
    Grid_Weather_Stations=[]    

    ploygons_to_lookat_new=[]
    for pobj in ploygons_to_lookat:
      #print pobj, obj_polygon
      if pobj[0].intersection(obj_polygon):
        inter_cnt+=1
        ploygons_to_lookat_new.append(pobj)


    #print inter_cnt
    if inter_cnt> 2 and zl<ClosestZoomLevel:
      GridMap(curQuadKey, ploygons_to_lookat_new)

      continue

    quadkeystr=get_quadkeystr(curQuadKey)
    

    
    for pobj in  ploygons_to_lookat_new:
      if PlotIt:
        x=[]
        y=[]      
        for point in list(obj_polygon.exterior.coords):
          x.append(point[0])
          y.append(point[1])
        plot(x,y, alpha=0.3)

      print quadkeystr+"\t"+pobj[1]+"\t"+str(pobj[2][0])+"\t"+str(pobj[2][1])


def GridVoronoi(voronoi_lattice, zl, PlotMap=True):
  #print "Mapping the Grids"
  
  global PlotIt
  PlotIt=PlotMap
  global ClosestZoomLevel
  ClosestZoomLevel=zl

  if PlotMap:
    plot_voronoi(voronoi_lattice)

  ploygons_to_lookat_new=[]
  for station, data in voronoi_lattice.items():
    ploygons_to_lookat_new.append((data["obj_polygon"], data["info"], data["coordinate"] ))


  GridMap("", ploygons_to_lookat_new)

  if PlotMap:
    show()
