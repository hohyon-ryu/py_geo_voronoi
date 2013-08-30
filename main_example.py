#! /usr/local/bin/python

import sys
import voronoi_poly

#Run this using the pipe: "cat sample_city_data.csv | python main_example.py"
if __name__=="__main__":

  #Creating the PointsMap from the input data
  PointsMap={}
  for line in sys.stdin:
    data=line.strip().split(",")
    try:
      PointsMap[data[0]]=(float(data[1]),float(data[2]))
    except:
      sys.stderr.write( "(warning) Invalid Input Line: "+line)

  #1. Stations, Lines and edges
  #vl=voronoi_poly.VoronoiLineEdges(PointsMap)
  #vertices, lines, edges, station_to_edge
  #print vl

  #2. Polygons
  vl=voronoi_poly.VoronoiPolygons(PointsMap, BoundingBox="W", PlotMap=False)
  #print vl

  #3. Quadkey-based Grids on Polygons
  voronoi_poly.GridVoronoi(vl, zl=7, PlotMap=True)  

  #4. GeoJson Polygons
  #cat sample_city_data | python main_example.py > out.geojson
  #voronoi_poly.VoronoiGeoJson_Polygons(PointsMap, BoundingBox="W", PlotMap=False)

  #5. GeoJson MultiPolygons
  #voronoi_poly.VoronoiGeoJson_MultiPolygons(PointsMap, BoundingBox="W", PlotMap=True)

  
  
