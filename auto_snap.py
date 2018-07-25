# -*- coding: utf-8 -*-
'''
	auto snap les entités entre elles
'''
import geopandas as gpd
# import time
import os
# import sys
import pdb
import shapely.geometry as geom
from shapely.ops import snap, nearest_points
from fiona.crs import from_epsg

gdf = gpd.read_file('SHAPES/ZONE SRO.shp')
# gdf.crs = from_epsg(2154)
snapped_vertices = {}

def insert_vertex_at(geometry, vertex1, vertex2, inplace = False):
	point = geom.Point(vertex2)
	if isinstance(geometry, geom.Polygon):
		linear_ring = geom.LinearRing(geometry.exterior.coords)
		coords = linear_ring.coords[:]
		for i in range(0, len(coords) - 1):
			line = geom.LineString([coords[i], coords[i+1]])
			if point.distance(line) <= 1:
				if inplace:
					coords[i] = vertex1
					vertex_at = i
				else:
					coords.insert(i+1, vertex1)
					vertex_at = i + 1
				break
		# pdb.set_trace()
		return geom.Polygon(coords), vertex_at

def snap_geom(geom1, geom2, index1, index2, tolerance):
	if isinstance(geom1, geom.Polygon) and isinstance(geom2, geom.Polygon):
		vertices1 = geom1.exterior.coords[:]
		vertices2 = geom2.exterior.coords[:]
		for i in range(0, len(vertices1)):
			vertex1 = vertices1[i]
			closest_vertex = get_closest_vertex(vertex1, vertices2, tolerance)
			
			if i not in snapped_vertices[str(index1)] and closest_vertex != None:
				snapped_vertices[str(index1)].append(i)
				# if code_entite in ['NRO_08_003_044', 'NRO_08_003_157']:
				# 	print "vertex ", vertex1, " at index ", i, " with closest vertex ", closest_vertex
				vertices1[i] = closest_vertex
			elif geom.Point(vertex1).distance(geom2) <= tolerance:
				snapped_vertices[str(index1)].append(i)
				linear_ring = geom.LinearRing(geom2.exterior.coords)
				geom2, vertex_added_at = insert_vertex_at(geom2, vertex1, nearest_points(geom.Point(vertex1), linear_ring)[1].coords[0])
				gdf.loc[index2, 'geometry'] = geom2
				snapped_vertices[str(index2)].append(vertex_added_at)
		return geom.Polygon(vertices1)
		

def get_closest_vertex(vertex1, vertices2, tolerance):

	p1 = geom.Point(vertex1)
	p2 = geom.Point(vertices2[0])

	min_distance = p1.distance(p2)
	closest_vertex = vertices2[0]

	for vertex2 in vertices2[1:]:
		p2 = geom.Point(vertex2)
		distance = p1.distance(p2)

		if distance < min_distance:
			min_distance = distance
			closest_vertex = vertex2

	if min_distance <= tolerance:
		return closest_vertex
	else:
		return None


def snap_gdf(tolerance):
	for index, row in gdf.iterrows():
		snapped_vertices[str(index)] = []

	for index, row in gdf.iterrows():
		snapped_vertices[row['CODE SRO']] = []
		tmp_gdf = gdf.copy()
		tmp_gdf['distance'] = tmp_gdf.distance(row['geometry'])
		tmp_gdf = tmp_gdf[tmp_gdf['distance'] <= tolerance]
		snapped_geom = list(gdf.loc[[index]]['geometry'])[0]
		
		for index_, row_ in tmp_gdf.iterrows():
			if index != index_:
				snapped_geom = snap_geom(snapped_geom, list(gdf.loc[[index_]]['geometry'])[0], index, index_, tolerance)
			

		gdf.loc[index, 'geometry'] = snapped_geom
		gdf.to_file('SHAPES/PM_snapped.shp')

snap_gdf(100)