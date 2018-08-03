# -*- coding: utf-8 -*-
'''
	auto snap les entités entre elles
'''
import geopandas as gpd
import time
import os
import sys
import pdb
import shapely.geometry as geom
from shapely.ops import nearest_points
# from fiona.crs import from_epsg


class AutoSnapper:
	"""docstring foSnapperme"""
	def __init__(self, shapes_dir):
		os.chdir(shapes_dir)
		self.t0 = time.time()
		self.snapped_vertices = None
		self.tf = None

	def snap(self, layer, tolerance):

		self.snapped_vertices = {}
		gdf = gpd.read_file('SHAPES/' + layer + '.shp')
		# gdf.crs = from_epsg(2154)

		if isinstance(list(gdf['geometry'])[0], geom.Polygon):
			self.snap_gdf(gdf, layer, tolerance)
		elif isinstance(list(gdf['geometry'])[0], geom.LineString):
			self.snap_in_centroids(gdf, layer, tolerance)

		self.tf = time.time()
		print "execution time : ", self.tf-self.t0, " seconds" 

#------------------- SAME GDF / GEOMETRY TYPE = POLYGON ----------------------------------

	def insert_vertex_at(self, geometry, vertex1, vertex2, inplace = False):

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

	def snap_geom(self, gdf, geom1, geom2, index1, index2, tolerance):
		if isinstance(geom1, geom.Polygon) and isinstance(geom2, geom.Polygon):
			vertices1 = geom1.exterior.coords[:]
			vertices2 = geom2.exterior.coords[:]
			for i in range(0, len(vertices1)):
				vertex1 = vertices1[i]
				closest_vertex, vertex_index = self.get_closest_vertex(vertex1, vertices2, tolerance)
				
				if i not in self.snapped_vertices[str(index1)] and closest_vertex != None:
					self.snapped_vertices[str(index1)].append(i)
					# if code_entite in ['NRO_08_003_044', 'NRO_08_003_157']:
					# 	print "vertex ", vertex1, " at index ", i, " with closest vertex ", closest_vertex
					vertices1[i] = closest_vertex
				elif geom.Point(vertex1).distance(geom2) <= tolerance:
					self.snapped_vertices[str(index1)].append(i)
					linear_ring = geom.LinearRing(geom2.exterior.coords)
					
					geom2, vertex_added_at = self.insert_vertex_at(geom2, vertex1, nearest_points(geom.Point(vertex1), linear_ring)[1].coords[0])
					# if vertex_added_at > vertex_index +1:
					# 	print "different indices : get_closest_vertex index = ", vertex_index, " / insert_at index = ", vertex_added_at
					# elif vertex_added_at == vertex_index + 1:
					# 	print "same indices = ", vertex_added_at
					# else:
					# 	print "WTF case !"
					gdf.loc[index2, 'geometry'] = geom2
					self.snapped_vertices[str(index1)].append(vertex_added_at)
					
			return geom.Polygon(vertices1)	

	def get_closest_vertex(self, vertex1, vertices2, tolerance):

		p1 = geom.Point(vertex1)
		p2 = geom.Point(vertices2[0])

		min_distance = p1.distance(p2)
		closest_vertex = vertices2[0]
		vertex_index = 0

		for i in range(1, len(vertices2)):
		# for vertex2 in vertices2[1:]:
			vertex2 = vertices2[i]
			p2 = geom.Point(vertex2)
			distance = p1.distance(p2)

			if distance < min_distance:
				min_distance = distance
				closest_vertex = vertex2
				vertex_index = i

		if min_distance <= tolerance:
			return closest_vertex, vertex_index
		else:
			return None, vertex_index

	def snap_gdf(self, gdf, layer, tolerance):
		for index, row in gdf.iterrows():
			self.snapped_vertices[str(index)] = []

		for index, row in gdf.iterrows():
			tmp_gdf = gdf.copy()
			tmp_gdf['distance'] = tmp_gdf.distance(row['geometry'])
			indices = tmp_gdf.index[tmp_gdf['distance'] <= tolerance]
			snapped_geom = list(gdf.loc[[index]]['geometry'])[0]
			
			print "snapping elements at row ", index
			# for index_, row_ in tmp_gdf.iterrows():
			for index_ in indices:
				if index != index_:
					snapped_geom = self.snap_geom(gdf, snapped_geom, list(gdf.loc[[index_]]['geometry'])[0], index, index_, tolerance)
				

			gdf.loc[index, 'geometry'] = snapped_geom
			gdf.to_file('SHAPES/' + layer +'_snapped.shp')

#------------------------------------------------------------------------------------------

#------------------ SAME GDF / GEOMETRY TYPE = LINESTRING -------------------------
	def snap_in_centroids(slef, gdf, layer, tolerance):
		'''
			valid only for LineString Geometries
		'''
		for index, row in gdf.iterrows():
			print "snapping elements at head ", index
			head = geom.Point(list(gdf.loc[[index]]['geometry'])[0].coords[0])
			self.snap_geom_to_centroid(head, gdf, tolerance)

			print "snapping elements at tail ", index
			tail = geom.Point(list(gdf.loc[[index]]['geometry'])[0].coords[-1])
			self.snap_geom_to_centroid(tail, gdf, tolerance)

		gdf.to_file('SHAPES/' + layer + '_snapped.shp')

	def snap_geom_to_centroid(self, point, gdf, tolerance):
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

#---------------------------------------------------------------------------------------------


if __name__ == '__main__':
	layer = sys.argv[1]
	tolerance = float(sys.argv[2])
	process = AutoSnapper(os.getcwd())
	process.snap(layer, tolerance)

