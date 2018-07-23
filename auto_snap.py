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
gdf.crs = from_epsg(2154)
snapped_vertices = {}

def snap_geom(geom1, geom2, tolerance, code_entite):
	if isinstance(geom1, geom.Polygon) and isinstance(geom2, geom.Polygon):
		vertices1 = geom1.exterior.coords[:]
		vertices2 = geom2.exterior.coords[:]
		for i in range(0, len(vertices1)):
			vertex1 = vertices1[i]
			closest_vertex = get_closest_vertex(vertex1, vertices2, tolerance)
			if code_entite == 'NRO_08_003_056':
				# if i == 35:
				# 	pdb.set_trace()
				print snapped_vertices[code_entite]
			if i not in snapped_vertices[code_entite]:
				if closest_vertex != None:
					snapped_vertices[code_entite].append(i)
					if code_entite == 'NRO_08_003_056':
						print "vertex ", vertex1, " at index ", i, " with closest vertex ", closest_vertex
					vertices1[i] = closest_vertex
				elif geom.Point(vertex1).distance(geom2) <= tolerance:
					snapped_vertices[code_entite].append(i)
					linear_ring = geom.LinearRing(geom2.exterior.coords)
					if code_entite == 'NRO_08_003_056':
						# if i == 35:
						# 	pdb.set_trace()
						print "vertex ", vertex1, " at index ", i, " with projected point ", nearest_points(geom.Point(vertex1), linear_ring)[1].coords[0]
					vertices1[i] = nearest_points(geom.Point(vertex1), linear_ring)[1].coords[0]
					# pdb.set_trace()
			elif code_entite == 'NRO_08_003_056':
				print "vertex ", vertex1, " at index ", i, " already snapped"
		try:
			return geom.Polygon(vertices1)
		except TypeError:
			pdb.set_trace()


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


for index, row in gdf.iterrows():
	snapped_vertices[row['CODE SRO']] = []
	tmp_gdf = gdf.copy()
	tmp_gdf['distance'] = tmp_gdf.distance(row['geometry'])
	tmp_gdf = tmp_gdf[tmp_gdf['distance'] <= 100]
	snapped_geom = row['geometry']
	
	for index_, row_ in tmp_gdf.iterrows():
		# if [row['CODE SRO'], row_['CODE SRO']] not in snap_combinations and \
		# 	[row_['CODE SRO'], row['CODE SRO']] not in snap_combinations and \
		if row['CODE SRO'] != row_['CODE SRO']:

			# snap_combinations.append([row['CODE SRO'], row_['CODE SRO']])
			print "snapping zone ", row['CODE SRO'], " with ", row_['CODE SRO']
			# snapped_geom = snap(snapped_geom, row_['geometry'], 100)
			snapped_geom = snap_geom(snapped_geom, row_['geometry'], 100, row['CODE SRO'])
		else:
			print row['CODE SRO'], " and ", row_['CODE SRO'], " already snapped"
	# closest_geom = list(tmp_gdf.sort_values('distance')['geometry'])[1]
	# print "snapping zone ", row['CODE SRO'], " with ", list(tmp_gdf.sort_values('distance')['CODE SRO'])[1]
	# I took 1 because index 0 would be the row itself
	# snapped_geom = snap(row['geometry'], closest_geom, 100)
	gdf.loc[index, 'geometry'] = snapped_geom
	# gdf.set_value(index, 'geometry', snapped_geom)



pdb.set_trace()