# -*- coding: utf-8 -*-
'''
	auto snap les entités entre elles
'''
import geopandas as gpd
import time
import os
# import sys
import pdb
import shapely.geometry as geom
from shapely.ops import nearest_points, snap
from fiona.crs import from_epsg

t0 = time.time()

gdf = gpd.read_file('SHAPES/INFRASTRUCTURES.shp')
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
				return geom.Polygon(coords), vertex_at
				# pdb.set_trace()
				# break
		# pdb.set_trace()
		try:
			return geom.Polygon(coords), None
		except UnboundLocalError:
			pdb.set_trace()

def snap_geom(gdf1, gdf2, geom1, geom2, index1, index2, tolerance):
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
				if gdf1.equals(gdf2):
					geom2, vertex_added_at = insert_vertex_at(geom2, vertex1, nearest_points(geom.Point(vertex1), linear_ring)[1].coords[0])
					gdf2.loc[index2, 'geometry'] = geom2
					snapped_vertices[str(index1)].append(vertex_added_at)
				else:
					# geom1, vertex_added_at = insert_vertex_at(geom1, vertex1, nearest_points(geom.Point(vertex1), linear_ring)[1].coords[0])
					vertices1[i] = nearest_points(geom.Point(vertex1), linear_ring)[1].coords[0]
					# gdf2.loc[index2, 'geometry'] = geom1
		return geom.Polygon(vertices1)
	

def densify_geometry(gdf2, geom1, geom2, index2, tolerance):
	if isinstance(geom1, geom.Polygon) and isinstance(geom2, geom.Polygon):
		vertices1 = geom1.exterior.coords[:]
		vertices2 = geom2.exterior.coords[:]
		for i in range(0, len(vertices1)):
			vertex1 = vertices1[i]
			closest_vertex = get_closest_vertex(vertex1, vertices2, tolerance)
			
			if i not in snapped_vertices[str(index2)] and closest_vertex != None:
					geom2, vertex_added_at = insert_vertex_at(geom2, vertex1, closest_vertex, inplace = True)
					snapped_vertices[str(index2)].append(i)
			elif geom.Point(vertex1).distance(geom2) <= tolerance:
				linear_ring = geom.LinearRing(geom2.exterior.coords)
				geom2, vertex_added_at = insert_vertex_at(geom2, vertex1, nearest_points(geom.Point(vertex1), linear_ring)[1].coords[0])
				increment_indices(vertex_added_at, index2)
				snapped_vertices[str(index2)].append(vertex_added_at)
			gdf2.loc[index2, 'geometry'] = geom2

def increment_indices(starting_index, index):
	snapped_vertices[str(index)] = [i+1 if i >= starting_index else i for i in snapped_vertices[str(index)]]

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
		tmp_gdf = gdf.copy()
		tmp_gdf['distance'] = tmp_gdf.distance(row['geometry'])
		tmp_gdf = tmp_gdf[tmp_gdf['distance'] <= tolerance]
		snapped_geom = list(gdf.loc[[index]]['geometry'])[0]
		
		for index_, row_ in tmp_gdf.iterrows():
			if index != index_:
				snapped_geom = snap_geom(gdf, gdf, snapped_geom, list(gdf.loc[[index_]]['geometry'])[0], index, index_, tolerance)
			

		gdf.loc[index, 'geometry'] = snapped_geom
		gdf.to_file('SHAPES/infra_snapped.shp')

def snap_gdf_to_gdf(gdf1, gdf2, tolerance):
	for index, row in gdf2.iterrows():
		snapped_vertices[str(index)] = []

	for index, row in gdf1.iterrows():
		gdf2['distance'] = gdf2.distance(row['geometry'])
		tmp_gdf = gdf2[gdf2['distance'] <= tolerance]
		# snapped_geom = list(gdf1.loc[[index]]['geometry'])[0]
		
		for index_, row_ in tmp_gdf.iterrows():
			print "snapping ", row['CODE SRO'], " to ", row_['CODE SRO']
			# snapped_geom = snap(list(gdf2.loc[[index_]]['geometry'])[0], list(gdf1.loc[[index]]['geometry'])[0], tolerance)
			# snapped_geom = snap_geom(gdf1, gdf2, list(gdf2.loc[[index_]]['geometry'])[0], list(gdf1.loc[[index]]['geometry'])[0], index, index_, tolerance)
			# gdf2.loc[index_, 'geometry'] = snapped_geom
			densify_geometry(gdf2, list(gdf1.loc[[index]]['geometry'])[0], list(gdf2.loc[[index_]]['geometry'])[0], index_, tolerance)
		gdf2.to_file('SHAPES/PM_snapped.shp')

def snap_geodataframe(gdf1, gdf2, tolerance):
	same_gdf = gdf1.equals(gdf2)

	for index, row in gdf1.iterrows():
		gdf2['distance'] = gdf2.distance(list(gdf1.loc[[index]]['geometry'])[0])
		gdf_to_snap = gdf2[gdf2['distance'] <= tolerance]

		for index_, row_ in gdf2.iterrows():
			if not same_gdf or index != index_:
				print "snapping ", row['CODE INFRA'], " and ", row_['CODE INFRA']
				geom1 = list(gdf2.loc[[index_]]['geometry'])[0]
				geom2 = list(gdf1.loc[[index]]['geometry'])[0]
				gdf2.loc[index_, 'geometry'] = snap(geom1, geom2, tolerance)
	gdf2.to_file('SHAPES/infra_snapped.shp')

def snap_in_centroids(gdf, tolerance):
	'''
		valid only for LineString Geometries
	'''
	for index, row in gdf.iterrows():
		print "snapping elements at head ", index
		head = geom.Point(list(gdf.loc[[index]]['geometry'])[0].coords[0])
		snap_geom_to_centroid(head, gdf, tolerance)

		print "snapping elements at tail ", index
		tail = geom.Point(list(gdf.loc[[index]]['geometry'])[0].coords[-1])
		snap_geom_to_centroid(tail, gdf, tolerance)

	gdf.to_file('SHAPES/infra_snapped.shp')

def snap_geom_to_centroid(point, gdf, tolerance):
	'''
		valid only for LineString Geometries
	'''
	gdf['distance'] = gdf.distance(point)
	indices = gdf.index[gdf['distance'] <= tolerance].tolist()

	multi_points = {}
	for index in indices:
		this_head = geom.Point(list(gdf.loc[[index]]['geometry'])[0].coords[0])
		this_tail = geom.Point(list(gdf.loc[[index]]['geometry'])[0].coords[-1])

		if point.distance(this_head) <= tolerance:
			multi_points[index] = {'point': this_head, 'place': 0}
		elif point.distance(this_tail) <= tolerance:
			multi_points[index] = {'point': this_tail, 'place': -1}
	# if len(indices) > 2:
	# 	pdb.set_trace()
	mp = geom.MultiPoint([multi_points[i]['point'] for i in multi_points.keys()])
	centroid = mp.centroid

	for i in multi_points.keys():
		geometry = list(gdf.loc[[i]]['geometry'])[0]
		coords = geometry.coords[:]
		coords[multi_points[i]['place']] = centroid.coords[0]
		gdf.loc[i, 'geometry'] = geom.LineString(coords)


# snap_gdf(100)
# gdf2 = gpd.read_file('SHAPES/INFRASTRUCTURES.shp')
# gdf1 = gpd.read_file('SHAPES/INFRASTRUCTURES.shp')
# snap_gdf_to_gdf(gdf1, gdf2, 100)
# snap_geodataframe(gdf1, gdf2, 3)

snap_in_centroids(gdf, 5)


tf = time.time()
print "execution time : ", tf-t0, " seconds" 