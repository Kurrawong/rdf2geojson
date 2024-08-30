from typing import List, Union, Optional, Tuple, Dict, AnyStr, Any

from rdflib import BNode, Graph, Literal, URIRef, DCTERMS
from rdflib.namespace import GEO, RDF, RDFS, SDO, NamespaceManager

from geojson import (
    FeatureCollection,
    Feature,
    GeometryCollection,
    MultiPolygon,
    Polygon,
    MultiLineString,
    LineString,
    MultiPoint,
    Point,
    GeoJSON,
    loads as geojson_loads,
    dumps as geojson_dumps,
)
from pyshacl import validate
from pathlib import Path
from .contrib.geomet import wkt


def get_geosparql_validator() -> Graph:
    try:
        return Graph().parse(Path(__file__).parent / "geosparql-validator.ttl")
    except FileNotFoundError:
        return Graph().parse(
            "https://raw.githubusercontent.com/opengeospatial/ogc-geosparql/master/vocabularies/validator.ttl"
        )


def make_json_key_from_iri(iri: URIRef, ns: NamespaceManager) -> Tuple[Optional[Tuple[str, str]], str]:
    """
    Returns either a tuple of ((namespace, prefix) local name)
    or (None, full_uri_str)
    :param iri:
    :return:
    """
    id_ = str(iri).split("?")[-1]
    try:
        (prefix, namespace, name) = ns.compute_qname(id_)
    except KeyError:
        # compute_qname will raise a KeyError if there is no known prefix
        (prefix, namespace, name) = None, None, id_
    if prefix is not None:
        return (namespace, prefix), name
    return None, name


def make_json_representation_of_obj(obj: Union[Literal, URIRef]) -> object:
    if isinstance(obj, URIRef):
        return str(obj)
    else:  # Literal
        return obj.toPython()


def parse_geometry(
    geom: Literal,
) -> Union[
    GeometryCollection,
    MultiPolygon,
    Polygon,
    MultiLineString,
    LineString,
    MultiPoint,
    Point,
]:
    # parse the GeoSPARQL geometry literal based on type
    # TODO: support all other GeoSPARQL 1.1 datatypes other than DGGS - GML, KML
    if geom.datatype == GEO.wktLiteral:
        return GeoJSON.to_instance(wkt.loads(geom))
    elif geom.datatype == GEO.geoJSONLiteral:
        return geojson_loads(geom)
    elif geom.datatype in [GEO.gmlLiteral, GEO.kmlLiteral]:
        raise NotImplementedError(
            "This GeoSPARQL geometry serialization format is not yet handled but "
            "eventually will be"
        )
    else:
        raise ValueError(
            "The serialization format of the supplied geometry is not one of "
            "GeoSPARQL 1.1's other than DGGS, as required"
        )

def _extract_geoms(g: Graph, pred, obj) -> List:
    geoms = []
    coords = g.value(obj, GEO.asWKT)
    if coords:
        geoms.append(parse_geometry(coords))
    coords = g.value(obj, GEO.asGeoJSON)
    if coords:
        geoms.append(parse_geometry(coords))
    else:
        # TODO handle unsupported GeosPARQL geometry serialization formats
        pass
    return geoms

def _extract_additional_property(g: Graph, pred, obj) -> Tuple[Union[str, URIRef], Any]:
    key_name = None
    value = None
    property_ids = list(g.objects(obj, SDO.propertyID))
    if len(property_ids) > 0:
        key_name = URIRef(property_ids[0])
    if key_name is None:
        text_names = list(g.objects(obj, SDO.name))
        if len(text_names) > 0:
            key_name = str(text_names[0])
    sdo_values_list = list(g.objects(obj, SDO.value))
    if len(sdo_values_list) > 0:
        value = sdo_values_list[0]
    if value is None:
        rdf_values_list = list(g.objects(obj, RDF.value))
        if len(rdf_values_list) > 0:
            value = rdf_values_list[0]
    if key_name is None:
        raise ValueError("Could not find a name or propertyID for additionalProperty")
    elif value is None:
        raise ValueError("Could not find a value for additionalProperty")
    return key_name, value

def get_converted_features(g: Graph, fc: Optional[URIRef] = None) -> List[Feature]:
    fs = []
    if fc is not None:
        feature_finder = g.objects(fc, RDFS.member)
    else:
        feature_finder = g.subjects(RDF.type, GEO.Feature)
    for f in feature_finder:
        # TODO: handle multiple Geometries per Feature
        geoms = []
        props = {}
        prop_contexts = {}
        for pred, obj in g.predicate_objects(f):
            if pred in [GEO.hasGeometry, GEO.hasDefaultGeometry]:
                geoms.extend(_extract_geoms(g, pred, obj))
            elif pred == SDO.spatial:
                # The Schema.org version of a GeoSpatial feature
                spatial_node = obj
                for p2, o2 in g.predicate_objects(spatial_node):
                    if p2 in [GEO.hasGeometry, GEO.hasDefaultGeometry]:
                        geoms.extend(_extract_geoms(g, p2, o2))

            elif pred == SDO.additionalProperty:
                # This is the Schema.org version of a Key-Value pair
                p_key, p_value =_extract_additional_property(g, pred, obj)
                if isinstance(p_key, URIRef):
                    prefix_pair, name = make_json_key_from_iri(URIRef(p_key), g.namespace_manager)
                    if prefix_pair is not None:
                        prefix_ns, prefix_name = prefix_pair
                        prop_contexts[prefix_name] = prefix_ns
                        p_key = f"{prefix_name}:{name}"
                    else:
                        p_key = str(p_key)
                props[p_key] = p_value
            elif isinstance(obj, BNode):
                # TODO: do something with Blank Nodes
                pass
            else:
                prefix_pair, name = make_json_key_from_iri(pred, g.namespace_manager)
                if prefix_pair is not None:
                    prefix_ns, prefix_name = prefix_pair
                    prop_contexts[prefix_name] = prefix_ns
                    name = f"{prefix_name}:{name}"
                props[name] = make_json_representation_of_obj(obj)
        if len(prop_contexts) > 0:
            props["@context"] = prop_contexts
        if geoms:
            fs.append(Feature(f, geometry=geoms[0], properties=props))
    return fs


def convert(g: Graph) -> GeoJSON:
    # validate the RDF data according to GeoSPARQL
    conforms, results_graph, results_text = validate(
        g,
        shacl_graph=get_geosparql_validator(),
    )
    if not conforms:
        return {}

    # TODO: consider handling multiple Feature Collections, as get_features(g, fc) allows for
    features = get_converted_features(g)
    if features:
        return FeatureCollection(features)
    else:
        return {}


def unconvert_geometry(geom: dict) -> Tuple[str, str]:
    # returns a WKT and GeoJSON representation of the unconverted geometry
    wkt_string = wkt.dumps(geom)
    geojson_string = geojson_dumps(geom)
    return wkt_string, geojson_string


def get_unconverted_features(g: Graph, gj: GeoJSON):
    type_ = gj["type"]
    if type_ == "FeatureCollection":
        features = gj["features"]
    elif type_ == "Feature":
        features = [gj]
    elif type_ == "GeometryCollection":
        features = [{"type": "Feature", "geometry": gj}]
    else:
        raise NotImplementedError(
            f"Not Implemented unconvert for GeoJSON type: {type_}"
        )
    for f in features:
        id_ = URIRef(f.id)
        unconverted_geom = unconvert_geometry(f.geometry)
        bn = BNode()
        g.add((id_, GEO.hasGeometry, bn))
        g.add((id_, RDF.type, GEO.Feature))
        g.add((bn, RDF.type, GEO.Geometry))
        g.add((bn, GEO.asWKT, Literal(unconverted_geom[0], datatype=GEO.wktLiteral)))
        g.add(
            (
                bn,
                GEO.asGeoJSON,
                Literal(unconverted_geom[1], datatype=GEO.geoJSONLiteral),
            )
        )
        properties = f.get("properties", {})
        for k, v in properties.items():
            if k == "title":
                g.add((id_, DCTERMS.title, Literal(v)))
            elif k == "label":
                g.add((id_, RDFS.label, Literal(v)))
            elif k == "description":
                g.add((id_, DCTERMS.description, Literal(v)))
            elif k == "identifier":
                g.add((id_, DCTERMS.identifier, Literal(v)))
            elif k == "sfWithin":
                g.add((id_, GEO.sfWithin, URIRef(v)))


def unconvert(gj: GeoJSON) -> Graph:
    # very basic GeoJSON to RDF conversion, for testing roundtripping
    g = Graph()
    type_ = gj["type"]
    if type_ == "FeatureCollection":
        assert "features" in gj
    elif type_ == "Feature":
        assert "geometry" in gj
    elif type_ == "GeometryCollection":
        assert "geometries" in gj
    else:
        raise NotImplementedError(
            f"Not Implemented unconvert for GeoJSON type: {type_}"
        )
    get_unconverted_features(g, gj)
    return g
