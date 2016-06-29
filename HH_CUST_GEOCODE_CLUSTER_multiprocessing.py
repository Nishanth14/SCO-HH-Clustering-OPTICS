

import math
import json
from multiprocessing import Pool

################################################################################
# POINT
################################################################################

class Point:
    
    def __init__(self, latitude, longitude):
        
        self.latitude = latitude
        self.longitude = longitude
        self.cd = None              # core distance
        self.rd = None              # reachability distance
        self.processed = False      # has this point been processed?
        
    # --------------------------------------------------------------------------
    # calculate the distance between any two points on earth
    # --------------------------------------------------------------------------
    
    def distance(self, point):
        
        # convert coordinates to radians
        
        p1_lat, p1_lon, p2_lat, p2_lon = [math.radians(c) for c in
            self.latitude, self.longitude, point.latitude, point.longitude]
        
        numerator = math.sqrt(
            math.pow(math.cos(p2_lat) * math.sin(p2_lon - p1_lon), 2) +
            math.pow(
                math.cos(p1_lat) * math.sin(p2_lat) -
                math.sin(p1_lat) * math.cos(p2_lat) *
                math.cos(p2_lon - p1_lon), 2))

        denominator = (
            math.sin(p1_lat) * math.sin(p2_lat) +
            math.cos(p1_lat) * math.cos(p2_lat) *
            math.cos(p2_lon - p1_lon))
        
        # convert distance from radians to meters
        # note: earth's radius ~ 6372800 meters
        
        return math.atan2(numerator, denominator) * 6372800
        
    # --------------------------------------------------------------------------
    # point as GeoJSON
    # --------------------------------------------------------------------------
        
    def to_geo_json_dict(self, properties=None):
        
        return {
            'type': 'Feature',
            'geometry': {
                'type': 'Point',
                'coordinates': [
                    self.longitude,
                    self.latitude,
                ]
            },
            'properties': properties,
        }
 
    def __repr__(self):
        return '(%f, %f)' % (self.latitude, self.longitude)

################################################################################
# CLUSTER
################################################################################

class Cluster:
    
    def __init__(self, points):
        
        self.points = points
        
    # --------------------------------------------------------------------------
    # calculate the centroid for the cluster
    # --------------------------------------------------------------------------

    def centroid(self):
        
        return Point(sum([p.latitude for p in self.points])/len(self.points),
            sum([p.longitude for p in self.points])/len(self.points))
            
    # --------------------------------------------------------------------------
    # calculate the region (centroid, bounding radius) for the cluster
    # --------------------------------------------------------------------------
    
    def region(self):
        
        centroid = self.centroid()
        radius = reduce(lambda r, p: max(r, p.distance(centroid)), self.points)
        return centroid, radius
        
    # --------------------------------------------------------------------------
    # cluster as GeoJSON
    # --------------------------------------------------------------------------
        
    def to_geo_json_dict(self, user_properties=None):
        
        center, radius = self.region()
        properties = { 'radius': radius }
        if user_properties: properties.update(user_properties)
        
        return {
            'type': 'Feature',
            'geometry': {
                'type': 'Point',
                'coordinates': [
                    center.longitude,
                    center.latitude,
                ]
            },
            'properties': properties,
        }

################################################################################
# OPTICS
################################################################################

class Optics:
    
    def __init__(self, points, max_radius, min_cluster_size):
        
        self.points = points
        self.max_radius = max_radius                # maximum radius to consider
        self.min_cluster_size = min_cluster_size    # minimum points in cluster
    
    # --------------------------------------------------------------------------
    # get ready for a clustering run
    # --------------------------------------------------------------------------
    
    def _setup(self):
        
        for p in self.points:
            p.rd = None
            p.processed = False
        self.unprocessed = [p for p in self.points]
        self.ordered = []

    # --------------------------------------------------------------------------
    # distance from a point to its nth neighbor (n = min_cluser_size)
    # --------------------------------------------------------------------------
    
    def _core_distance(self, point, neighbors):

        if point.cd is not None: return point.cd
        if len(neighbors) >= self.min_cluster_size - 1:
            #print(point,neighbors,len(neighbors))
            sorted_neighbors = sorted([n.distance(point) for n in neighbors])
            point.cd = sorted_neighbors[self.min_cluster_size - 2]
            #print(point.cd)
            return point.cd
        
    # --------------------------------------------------------------------------
    # neighbors for a point within max_radius
    # --------------------------------------------------------------------------
    
    def _neighbors(self, point):
        
        return [p for p in self.points if p is not point and
            p.distance(point) <= self.max_radius]
            
    # --------------------------------------------------------------------------
    # mark a point as processed
    # --------------------------------------------------------------------------
        
    def _processed(self, point):
    
        point.processed = True
        self.unprocessed.remove(point)
        self.ordered.append(point)
    
    # --------------------------------------------------------------------------
    # update seeds if a smaller reachability distance is found
    # --------------------------------------------------------------------------

    def _update(self, neighbors, point, seeds):
        
        # for each of point's unprocessed neighbors n...

        for n in [n for n in neighbors if not n.processed]:
            
            # find new reachability distance new_rd
            # if rd is null, keep new_rd and add n to the seed list
            # otherwise if new_rd < old rd, update rd
            
            new_rd = max(point.cd, point.distance(n))
            if n.rd is None:
                n.rd = new_rd
                seeds.append(n)
            elif new_rd < n.rd:
                n.rd = new_rd
    
    # --------------------------------------------------------------------------
    # run the OPTICS algorithm
    # --------------------------------------------------------------------------

    def run(self):
        
        self._setup()
        
        # for each unprocessed point (p)...
        
        while self.unprocessed:
            point = self.unprocessed[0]
            
            # mark p as processed
            # find p's neighbors
            
            self._processed(point)
            point_neighbors = self._neighbors(point)

            # if p has a core_distance, i.e has min_cluster_size - 1 neighbors

            if self._core_distance(point, point_neighbors) is not None:
                
                # update reachability_distance for each unprocessed neighbor
                
                seeds = []
                self._update(point_neighbors, point, seeds)
                
                # as long as we have unprocessed neighbors...
                
                while(seeds):
                    
                    # find the neighbor n with smallest reachability distance
                    
                    seeds.sort(key=lambda n: n.rd)
                    n = seeds.pop(0)
                    
                    # mark n as processed
                    # find n's neighbors
                    
                    self._processed(n)
                    n_neighbors = self._neighbors(n)
                    
                    # if p has a core_distance...
                    
                    if self._core_distance(n, n_neighbors) is not None:
                        
                        # update reachability_distance for each of n's neighbors
                        
                        self._update(n_neighbors, n, seeds)
                        
        # when all points have been processed
        # return the ordered list
        #print(self.ordered)
        return self.ordered
        
    # --------------------------------------------------------------------------
    
    def cluster(self, cluster_threshold):
        
        clusters = []
        separators = []

        for i in range(len(self.ordered)):
            this_i = i
            next_i = i + 1
            this_p = self.ordered[i]
            this_rd = this_p.rd if this_p.rd else float('infinity')
            #print(this_p.rd,this_p)
            #print( cluster_threshold)
            
            # use an upper limit to separate the clusters
            
            if this_rd > cluster_threshold:
                separators.append(this_i)

        separators.append(len(self.ordered))

        for i in range(len(separators) - 1):
            start = separators[i]
            end = separators[i + 1]
            if end - start >= self.min_cluster_size:
                clusters.append(Cluster(self.ordered[start:end]))

        return clusters

# LOAD SOME POINTS

import numpy as np

import pyodbc


import pyodbc



conn1 = pyodbc.connect("DRIVER={NetezzaSQL};SERVER=npsdwh.con-way.com;PORT=5480;DATABASE=PRD_WHSE;UID=dalli.krishna;PWD=four5six;")
#conn2 = pyodbc.connect("DRIVER={NetezzaSQL};SERVER=npsdwh.con-way.com;PORT=5480;DATABASE=SANDBOX;UID=BAHUGUNAASHISH;PWD=june_2016;")
conn3 = pyodbc.connect("DRIVER={NetezzaSQL};SERVER=npsdwh.con-way.com;PORT=5480;DATABASE=PRD_WHSE;UID=YALAM.NISHANTH;PWD=nishanth519;")
curs1 = conn1.cursor()
#curs2 = conn2.cursor()
curs3 = conn3.cursor()

curs3.execute(
    "select distinct NODE_INST_ID from prd_rds..sco_trip_node_rds where LATITUDE_NBR <> 0 and LONGITUDE_NBR <> 0 and NODE_TYPE_CD <> 'SIC'  and STATUS_CD = 'COMP' and  TRIP_INST_ID in (SELECT DISTINCT TRP_INSTC_NBR   FROM PRD_WHSE..EOBR_VEHICLE_COORDINATE WHERE date(COOR_DTTM) = '2016-06-16' and  VEH_NBR in ('632 8284','632 8291','632 8297','632 8299','632 8301')) ")
rows3 = curs3.fetchall()
NODE_INST_ID = []


def cluster(NODE_INST_ID):
    curs1.execute("select A.LATITUDE_NBR as LATD,A.LONGITUDE_NBR as LNGT from prd_rds..sco_trip_node_rds A where LATITUDE_NBR <> 0 and LONGITUDE_NBR <> 0 and A.NODE_INST_ID in (?) ORDER BY  LST_UPDT_TMST DESC LIMIT 1000 ",[NODE_INST_ID])
    rows1 = curs1.fetchall()
    #print(rows1)
    print(len(rows1))
    length = len(rows1)
    no_points = int(round(length*0.03))

    lats = []
    if(no_points<=15):
       no = 2
       clust_radius = 100
    else:
        no = no_points
        clust_radius = 10
    print(no)
    for row1 in rows1:
        lats.append(Point(row1.LATD, row1.LNGT))

    points = lats
    print(clust_radius)
    optics = Optics(points, 100, no)  # 100m radius for neighbor consideration, cluster size >= 2 points
    optics.run()  # run the algorithm
    clusters = optics.cluster(clust_radius)  # 50m threshold for clustering
    myfile = open('C:\\Users\\yalam.nishanth\\Desktop\\Projects\\RPM Phase 2\\cluster_centroids.txt', 'a')
    print(clusters)
    for cluster in clusters:
        cluster_points = Cluster(cluster.points)

        centroid = cluster_points.centroid()
        print(centroid,NODE_INST_ID)

        #curs2.execute("INSERT INTO sandbox..SCO_HH_LAT_LNG_V5 (NODE_APPL_INSTC_NBR, LATD, LNGT) select ?,?,?", [NODE_INST_ID,centroid. latitude,centroid. longitude])
        '''var1 = [str(NODE_INST_ID),str(centroid. latitude),str(centroid. longitude)]
        print(var1)
        myfile.write(var1)
        print(cluster_points.region())
    myfile.close()'''
    print(len(clusters))
    #curs2.commit()

#cluster(700787595)
#958902139
if __name__ == '__main__':

  for row3 in rows3:
    NODE_INST_ID.append(row3.NODE_INST_ID)
  print(NODE_INST_ID)
  p=Pool(1)
  p.map(cluster,NODE_INST_ID)

