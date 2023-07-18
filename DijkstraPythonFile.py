import xmltodict as xtd
import folium
import numpy as np
import webbrowser
import os, sys
from haversine import haversine

# Codigo tomado de base para el proyecto: https://github.com/shivansh1012/DijkstraOnMaps/blob/main/DijkstraPythonFile.py


# Abriendo y leyendo el archivo .OSM
with open("Maps/santoMini.osm", "rb") as osm_fn:
    map_osm = xtd.parse(osm_fn)["osm"]

# Obteniendo los limites del mapa
ymax = map_osm["bounds"]["@maxlat"]
ymin = map_osm["bounds"]["@minlat"]
xmax = map_osm["bounds"]["@maxlon"]
xmin = map_osm["bounds"]["@minlon"]
parsed_bounds = [xmin, xmax, ymin, ymax]

# Transformando los datos de los nodos en listas
Node = map_osm["node"]
Nnodes = len(Node)
Nodeid = [0] * Nnodes
xy = []
for i in range(Nnodes):
    Nodeid[i] = float(Node[i]["@id"])
    x = float(Node[i]["@lat"])
    y = float(Node[i]["@lon"])
    xy.append([x, y])
parsed_node = {"id": Nodeid, "xy": xy}

# Transformando los datos de las vias en listas
Way = map_osm["way"]
Nways = len(Way)
Wayid = [0] * Nways
nodes_in_way = [0] * Nways
tags = [0] * Nways
for i in range(Nways):
    tempWay = Way[i]
    Wayid[i] = float(tempWay["@id"])
    Nnd = len(tempWay["nd"])
    ndTemp = [0] * Nnd
    for j in range(Nnd):
        ndTemp[j] = float(tempWay["nd"][j]["@ref"])
    nodes_in_way[i] = ndTemp
    if "tag" in tempWay.keys():
        if type(tempWay["tag"]) is list:
            tags[i] = tempWay["tag"]
        else:
            tags[i] = [tempWay["tag"]]
    else:
        tags[i] = []
parsed_way = {"id": Wayid, "nodes": nodes_in_way, "tags": tags}

# Transformando los datos de las relaciones en listas
Relation = map_osm["relation"]
Nrelation = len(Relation)
Relationid = [0] * Nrelation
for i in range(Nrelation):
    currentRelation = Relation[i]
    currentId = currentRelation["@id"]
    Relationid[i] = float(currentId)
parsed_relation = {"id": Relationid}

# Creando un diccionario con los datos de los nodos, vias y relaciones
parsed_osm = {
    "bounds": parsed_bounds,
    "relation": parsed_relation,
    "way": parsed_way,
    "node": parsed_node,
    "attributes": map_osm.keys(),
}

bounds = parsed_osm["bounds"]
way = parsed_osm["way"]
node = parsed_osm["node"]
relation = parsed_osm["relation"]

ways_num = len(way["id"])
ways_node_set = way["nodes"]
node_ids = dict()
n = len(node["id"])
for i in range(n):
    node_ids[node["id"][i]] = i

road_vals = [
    "highway",
    "motorway",
    "motorway_link",
    "trunk",
    "trunk_link",
    "primary",
    "primary_link",
    "secondary",
    "secondary_link",
    "tertiary",
    "road",
    "residential",
    "living_street",
    "service",
    "services",
    "motorway_junction",
]


# Creando la conectividad entre los nodos del mapa
def create_connectivity():
    connectivity_matrix = np.full((Nnodes, Nnodes), 99999)
    np.fill_diagonal(connectivity_matrix, 0)

    for currentWay in range(ways_num):
        skip = True
        for i in way["tags"][currentWay]:
            if i["@k"] in road_vals:
                skip = False
                break
        if skip:
            continue

        nodeset = ways_node_set[currentWay]
        nodes_num = len(nodeset)


        for firstnode_local_index in range(nodes_num - 1):
            firstnode_id = nodeset[firstnode_local_index]
            firstnode_index = node_ids.get(firstnode_id, -1)
            if firstnode_index == -1:
                continue

            secondnode_id = nodeset[firstnode_local_index + 1]
            secondnode_index = node_ids.get(secondnode_id, -1)
            if secondnode_index == -1:
                continue

            # Obtener las cordenadas de los nodos
            firstnode_lat, firstnode_lon = node["xy"][firstnode_index]
            secondnode_lat, secondnode_lon = node["xy"][secondnode_index]

            # Calculamos la distancia entre los nodos gracias a la funcion haversine
            distance = haversine(
                (firstnode_lat, firstnode_lon), (secondnode_lat, secondnode_lon))

            connectivity_matrix[firstnode_index, secondnode_index] = distance

    return connectivity_matrix


# Algoritmo UCS para buscar la ruta optima
def dijkstra(source, connectivity_matrix, parents):
    visited = dict()
    visited[source] = True
    parents[source] = source

    num_nodes = len(connectivity_matrix)
    current_node = source
    min_distance = float("inf")
    for i in range(num_nodes):
        if i != source and connectivity_matrix[source][i] < min_distance:
            current_node = i
            min_distance = connectivity_matrix[source][i]
    visited[current_node] = True
    parents[current_node] = source

    remaining_nodes = num_nodes - 2
    while remaining_nodes > 0:
        next_node = source
        min_distance = float("inf")

        for j in range(num_nodes):
            if (
                visited.get(j, False) is False
                and connectivity_matrix[source][current_node] != float("inf")
                and connectivity_matrix[current_node][j] != float("inf")
            ):
                total_distance = (
                    connectivity_matrix[source][current_node]
                    + connectivity_matrix[current_node][j]
                )
                connectivity_matrix[source][j] = min(
                    connectivity_matrix[source][j], total_distance
                )
                connectivity_matrix[j][source] = connectivity_matrix[source][j]

                if connectivity_matrix[source][j] == total_distance:
                    parents[j] = current_node
                elif connectivity_matrix[source][j] == 1:
                    parents[j] = source

                if connectivity_matrix[source][j] < min_distance:
                    next_node = j
                    min_distance = connectivity_matrix[source][j]

        if next_node == source:
            break
        visited[next_node] = True
        current_node = next_node
        remaining_nodes -= 1

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

print("Espere un momento, se esta generando el mapa...")

# Generar el mapa para mostrar los nodos
def BuildAllNodesMap():
    x1, y1 = (float(bounds[2]), float(bounds[0]))
    x2, y2 = (float(bounds[3]), float(bounds[1]))
    center = ((x1 + x2) / 2, (y1 + y2) / 2)
    map_0 = folium.Map(location=center, zoom_start=16)

    for i in range(n):
        xy = (node["xy"][i][0], node["xy"][i][1])
        folium.CircleMarker(
            xy, radius=3, color="black", fill=True, fill_color="black", popup=str(i)
        ).add_to(map_0)
    return map_0


# Generar el mapa para mostrar los nodos disponibles para el destino
def BuildAllNodesMapResponse(SourceNode):
    x1, y1 = (float(bounds[2]), float(bounds[0]))
    x2, y2 = (float(bounds[3]), float(bounds[1]))
    center = ((x1 + x2) / 2, (y1 + y2) / 2)
    map_0 = folium.Map(location=center, zoom_start=16)

    for i in range(n):
        xy = (node["xy"][i][0], node["xy"][i][1])
        if i != SourceNode:
            folium.CircleMarker(xy, radius = 3, color="purple", fill=True, fill_color="purple", popup=str(i),).add_to(map_0)
        else:
            folium.CircleMarker(
                xy, radius=9, color="green", fill=True, fill_color="red", popup=str(i)
            ).add_to(map_0)

    return map_0


# Generando el mapa para mostrar la ruta optima entre el nodo origen y el nodo destino
def BuildFinalPathMap(i, p):

    node_cds = [(node["xy"][i][0], node["xy"][i][1])]
    while p[i] != i:
        node_cds.append((node["xy"][p[i]][0], node["xy"][p[i]][1]))
        i = p[i]

    map_0 = folium.Map(location=node_cds[-1], zoom_start=15)

    folium.CircleMarker(
        node_cds[-1], radius=5, color="blue", fill=True, fill_color="orange"
    ).add_to(map_0)
    folium.Marker(
        node_cds[0], icon=folium.Icon(color="blue", icon="circle", prefix="fa")
    ).add_to(map_0)

    for j in range(len(node_cds) - 1):
        start = node_cds[j]
        end = node_cds[j + 1]
        folium.PolyLine(
            locations=[start, end],
            weight=5,
            color="blue",
            opacity="0.75",
            dash_array=10,
        ).add_to(map_0)

    return map_0


# Funcion para poder abrir el mapa HTML
def OpenHTMLMapinBrowser(filename):
    url = "file://" + os.path.realpath(filename)
    webbrowser.open(url, new=2)


# Empieza la interfaz, genreando el mapa de todos los nodos
map1 = BuildAllNodesMap()
map1.save("Mapa de partida.html")
OpenHTMLMapinBrowser("Mapa de partida.html")

# Interaccion con el usuario, pidiendo datos de entrada para la ruta
while True:
    SourceNode = int(input("Ingrese el nodo de partida o 0 para salir: "))
    connectivity_matrix = create_connectivity()
    nodes_routes_values, p = plot_routes(SourceNode, connectivity_matrix)

    if SourceNode == 0:
        print("Programa finalizado")
        break

    map2 = BuildAllNodesMapResponse(SourceNode)
    map2.save("Mapa de destino.html")
    OpenHTMLMapinBrowser("Mapa de destino.html")

    while True:
        DestinationNode = int(
            input(
                "Ingrese el nodo de destino o -1 para seleccionar un nuevo nodo de partida o 0 para salir: "
            )
        )

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
            map3.save("Mapa con ruta.html")
            OpenHTMLMapinBrowser("Mapa con ruta.html")
