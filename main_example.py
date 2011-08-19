#! /usr/local/bin/python

import sys
import voronoi_poly

#Run this using the pipe: "cat sample_city_data | python main_example.py"
if __name__=="__main__":

  #Creating the PointsMap from the input data
  PointsMap={}
  for line in sys.stdin:
    data=line.strip().split(",")
    try:
      PointsMap[data[0]]=(float(data[1]),float(data[2]))
    except:
      sys.stderr.write( "Invalid Input Line: "+line)

  #vl=voronoi_poly.VoronoiLineEdges(PointsMap)
  #vertices, lines, edges, station_to_edge
  #print vl

  #vl=voronoi_poly.VoronoiPolygons(PointsMap, BoundingBox="W", PlotMap=True)
  #print vl
  
  #cat sample_city_data | python main_example.py > out.geojson
  #vl=voronoi_poly.VoronoiGeoJson_Polygons(PointsMap, BoundingBox="W", PlotMap=False)

  vl=voronoi_poly.VoronoiGeoJson_MultiPolygons(PointsMap, BoundingBox="W", PlotMap=False)
  print vl

