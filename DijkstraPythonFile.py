import xmltodict as xtd
import folium
import numpy as np
import webbrowser
import os, sys
from haversine import haversine

#Parsing raw data from the .OSM File
with open('Maps/santoMini.osm', "rb") as osm_fn:
    map_osm = xtd.parse(osm_fn)['osm']

#Parsing bounds from .OSM file
ymax = map_osm['bounds']['@maxlat']
ymin = map_osm['bounds']['@minlat']
xmax = map_osm['bounds']['@maxlon']
xmin = map_osm['bounds']['@minlon']
parsed_bounds = [xmin, xmax, ymin, ymax]

#Parsing Node
Node=map_osm['node']
Nnodes=len(Node)
Nodeid = [0]*Nnodes
xy = []
for i in range(Nnodes):
    Nodeid[i]=float(Node[i]['@id'])
    x=float(Node[i]['@lat'])
    y=float(Node[i]['@lon'])
    xy.append([x,y])
parsed_node={'id':Nodeid, 'xy':xy}

#Parsing Ways
Way=map_osm['way']
Nways=len(Way)
Wayid=[0]*Nways
nodes_in_way=[0]*Nways
tags=[0]*Nways
for i in range(Nways):
    tempWay = Way[i]
    Wayid[i] = float(tempWay['@id'])
    Nnd=len(tempWay['nd'])
    ndTemp=[0]*Nnd
    for j in range(Nnd):
        ndTemp[j]=float(tempWay['nd'][j]['@ref'])
    nodes_in_way[i] = ndTemp
    if 'tag' in tempWay.keys():
        if type(tempWay['tag']) is list:
              tags[i]=tempWay['tag']
        else:
              tags[i]=[tempWay['tag']]
    else:
        tags[i]=[]
parsed_way={'id':Wayid,'nodes':nodes_in_way, 'tags':tags}

#Parsing Relations
Relation=map_osm['relation']
Nrelation=len(Relation)
Relationid=[0]*Nrelation
for i in range(Nrelation):
    currentRelation = Relation[i]
    currentId=currentRelation['@id']
    Relationid[i]=float(currentId)
parsed_relation={'id':Relationid}

#Parsing .OSM file
parsed_osm={
    'bounds':parsed_bounds,
    'relation':parsed_relation,
    'way':parsed_way,
    'node':parsed_node,
    'attributes':map_osm.keys()
}

bounds=parsed_osm['bounds']
way=parsed_osm['way']
node=parsed_osm['node']
relation=parsed_osm['relation']

ways_num = len(way['id'])
ways_node_set=way['nodes']
node_ids = dict()
n = len(node['id'])
for i in range(n):
    node_ids[node['id'][i]] = i

road_vals = ['highway', 'motorway', 'motorway_link', 'trunk', 'trunk_link',
             'primary', 'primary_link', 'secondary', 'secondary_link',
             'tertiary', 'road', 'residential', 'living_street',
             'service', 'services', 'motorway_junction']

#Creating Connectivity between the nodes
def create_connectivity():
    connectivity_matrix = np.full((Nnodes, Nnodes), 99999)  # Establecer una distancia alta pero finita
    np.fill_diagonal(connectivity_matrix, 0)

    for currentWay in range(ways_num):
        skip = True
        for i in way['tags'][currentWay]:
            if i['@k'] in road_vals:
                skip = False
                break
        if skip:
            continue

        nodeset = ways_node_set[currentWay]
        nodes_num = len(nodeset)

        currentWayID = way['id'][currentWay]

        for firstnode_local_index in range(nodes_num - 1):
            firstnode_id = nodeset[firstnode_local_index]
            firstnode_index = node_ids.get(firstnode_id, -1)
            if firstnode_index == -1:
                continue

            secondnode_id = nodeset[firstnode_local_index + 1]
            secondnode_index = node_ids.get(secondnode_id, -1)
            if secondnode_index == -1:
                continue

            # Obtenemos las coordenadas de los nodos
            firstnode_lat, firstnode_lon = node['xy'][firstnode_index]
            secondnode_lat, secondnode_lon = node['xy'][secondnode_index]

            # Calculamos la distancia entre los nodos
            distance = haversine((firstnode_lat, firstnode_lon), (secondnode_lat, secondnode_lon))

            connectivity_matrix[firstnode_index, secondnode_index] = distance

    return connectivity_matrix

#Dijkstra Algorithm used for finding the shortest path
def dijkstra(source, connectivity_matrix, p):
    s = dict()
    s[source] = True
    p[source] = source

    v = len(connectivity_matrix)
    u = source
    d_u = float('inf')
    for i in range(v):
        if i != source and connectivity_matrix[source][i] < d_u:
            u = i
            d_u = connectivity_matrix[source][i]
    s[u] = True
    p[u] = source

    i = v-2
    while i > 0:
        u_x = source
        d_u = float('inf')

        for j in range(v):
            if s.get(j, False) == False and connectivity_matrix[source][u] != float('inf') and connectivity_matrix[u][j] != float('inf'):
                k = connectivity_matrix[source][u] + connectivity_matrix[u][j]
                connectivity_matrix[source][j] = min(connectivity_matrix[source][j], k)
                connectivity_matrix[j][source] = connectivity_matrix[source][j]

                if connectivity_matrix[source][j] == k:
                    p[j] = u
                elif connectivity_matrix[source][j] == 1:
                    p[j] = source

                if connectivity_matrix[source][j] < d_u:
                    u_x = j
                    d_u = connectivity_matrix[source][j]

        if u_x == source: break
        s[u_x] = True
        u = u_x
        i -= 1


def plot_routes(s, connectivity_matrix):
    p = dict()
    dijkstra(s, connectivity_matrix, p)

    nodes_routes_values = []
    for i in p.keys():
        if i != s and i in p:
            adder = [i, 0]
            while p[i] != i:
                adder[1] += 1
                i = p[i]
            nodes_routes_values.append(adder)

    return nodes_routes_values, p


print("Please wait while all Nodes Map is Generating...")

#Generating a map to display all the nodes
def BuildAllNodesMap():
    x1, y1 = (float(bounds[2]), float(bounds[0]))
    x2, y2 = (float(bounds[3]), float(bounds[1]))
    center = ((x1+x2)/2, (y1+y2)/2)
    map_0 = folium.Map(location = center, zoom_start = 16)

    for i in range(n):
        xy = (node['xy'][i][0], node['xy'][i][1])
        folium.CircleMarker(xy, radius=3, color="green", fill=True, fill_color="green", popup=str(i)).add_to(map_0)
    return map_0



def BuildAllNodesMap():
    x1, y1 = (float(bounds[2]), float(bounds[0]))
    x2, y2 = (float(bounds[3]), float(bounds[1]))
    center = ((x1 + x2) / 2, (y1 + y2) / 2)
    map_0 = folium.Map(location=center, zoom_start=16)

    for i in range(n):
        xy = (node['xy'][i][0], node['xy'][i][1])
        folium.CircleMarker(xy, radius=3, color="green", fill=True, fill_color="green", popup=str(i)).add_to(map_0)

    return map_0






#Generating a map to display the path between source and destination
def BuildFinalPathMap(i, p):
    if i not in p:
        print("El nodo destino no está en la ruta óptima.")
        return None

    node_cds = [(node['xy'][i][0], node['xy'][i][1])]
    while p[i] != i:
        node_cds.append((node['xy'][p[i]][0], node['xy'][p[i]][1]))
        i = p[i]

    map_0 = folium.Map(location=node_cds[-1], zoom_start=15)

    folium.CircleMarker(node_cds[-1], radius=5, color="blue", fill=True, fill_color="orange").add_to(map_0)
    folium.Marker(node_cds[0], icon=folium.Icon(color="blue", icon="circle", prefix='fa')).add_to(map_0)

    for j in range(len(node_cds) - 1):
        start = node_cds[j]
        end = node_cds[j + 1]
        folium.PolyLine(locations=[start, end], weight=5, color="blue", opacity="0.75", dash_array=10).add_to(map_0)

    return map_0

#Function to open a html file in browser
def OpenHTMLMapinBrowser(filename):
    url = "file://" + os.path.realpath(filename)
    webbrowser.open(url,new=2)

#First Map Generator to show all the Nodes
map1 = BuildAllNodesMap()
map1.save("AllNodeMap.html")
OpenHTMLMapinBrowser("AllNodeMap.html")

#Third Map Generator to show path from source to destination
while True:
    SourceNode = int(input("Ingrese el nodo de partida o 0 para salir: "))
    connectivity_matrix = create_connectivity()
    nodes_routes_values, p = plot_routes(SourceNode, connectivity_matrix)

    if SourceNode == 0:
        print("Programa finalizado")
        break


    map1 = BuildAllNodesMap()
    map1.save("AllNodeMap.html")
    OpenHTMLMapinBrowser("AllNodeMap.html")

    while True:
        DestinationNode = int(input("Ingrese el nodo de destino o -1 para seleccionar un nuevo nodo de partida o 0 para salir: "))

        if DestinationNode == -1:
            break

        if DestinationNode == 0:
            print("Programa finalizado")
            sys.exit(1)

        if DestinationNode < 0 or DestinationNode >= n:
            print("El nodo de destino seleccionado no es válido.")
            continue

        map3 = BuildFinalPathMap(DestinationNode, p)
        if map3 is not None:
            map3.save("OutputMap.html")
            OpenHTMLMapinBrowser("OutputMap.html")