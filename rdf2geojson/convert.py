from typing import List, Union, Optional, Tuple, Dict, AnyStr, Any, Callable

from rdflib import BNode, Graph, Literal, URIRef, DCTERMS, SOSA, XSD, SKOS
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
    elif isinstance(obj, Literal):
        # Some literals cannot be represented as JSON, so we return a string
        if obj.datatype is None:
            return str(obj)
        elif obj.datatype in (XSD.date, XSD.dateTime, XSD.time):
            return str(obj)
    else:
        # The GeoJSON serializer will take care of converting this to JSON
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

def _extract_bnode(g: Graph, bn: URIRef, prop_contexts: Dict, recurse: int = 0) -> Dict:
    obs_dict = {}
    for pred, obj in g.predicate_objects(bn):
        prefix_pair, name = make_json_key_from_iri(pred, g.namespace_manager)
        if prefix_pair is not None:
            use_prefix = True
            prefix_ns, prefix_name = prefix_pair
            if prefix_name in prop_contexts:
                if prefix_ns != prop_contexts[prefix_name]:
                    # conflicting prefix with one thats already in there
                    use_prefix = False
            if use_prefix:
                prop_contexts[prefix_name] = prefix_ns
                name = f"{prefix_name}:{name}"
        if isinstance(obj, BNode):
            if recurse < 8:
                obs_dict[name] = _extract_bnode(g, obj, prop_contexts, recurse + 1)
        else:

            obs_dict[name] = make_json_representation_of_obj(obj)
    return obs_dict

def _extract_observation(g: Graph, obs: URIRef, prop_contexts: Dict) -> Dict:
    obs_dict = {"iri": str(obs)}
    members = [] # this could be an observationCollection too
    for pred, obj in g.predicate_objects(obs):
        if pred == SOSA.hasMember:
            members.append(_extract_observation(g, obj, prop_contexts))
        else:
            prefix_pair, name = make_json_key_from_iri(pred, g.namespace_manager)
            if prefix_pair is not None:
                use_prefix = True
                prefix_ns, prefix_name = prefix_pair
                if prefix_name in prop_contexts:
                    if prefix_ns != prop_contexts[prefix_name]:
                        # conflicting prefix with one thats already in there
                        use_prefix = False
                if use_prefix:
                    prop_contexts[prefix_name] = prefix_ns
                    name = f"{prefix_name}:{name}"

            if isinstance(obj, BNode):
                obs_dict[name] = _extract_bnode(g, obj, prop_contexts)
            else:
                obs_dict[name] = make_json_representation_of_obj(obj)
        if len(members) > 0:
            obs_dict["sosa:hasMember"] = members
    return obs_dict

def get_features_collections(g: Graph, iri2id: Optional[Callable[[URIRef], str]] = None) -> List[FeatureCollection]:
    feature_finder = g.subjects(RDF.type, GEO.FeatureCollection)
    fs = []
    for f in feature_finder:
        props = {}
        extras = {}
        prop_contexts = {}
        _id = None
        if iri2id is not None:
            _id = iri2id(f)
        for pred, obj in g.predicate_objects(f):
            if pred in (RDFS.label, SKOS.prefLabel) and not "title" in extras:
                extras["title"] = str(obj)
            elif pred == RDFS.member:
                # Skip the members, they are handled by get_features
                pass
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
            elif _id is None and pred == DCTERMS.identifier:
                _id = str(obj)  # the identifier becomes the ID
            else:
                prefix_pair, name = make_json_key_from_iri(pred, g.namespace_manager)
                if prefix_pair is not None:
                    prefix_ns, prefix_name = prefix_pair
                    prop_contexts[prefix_name] = prefix_ns
                    name = f"{prefix_name}:{name}"
                if isinstance(obj, BNode):
                    props[name] = _extract_bnode(g, obj, prop_contexts)
                else:
                    props[name] = make_json_representation_of_obj(obj)

        if _id is not None:
            # ID is not the same as IRI, so put iri in the properties
            props["iri"] = str(f)
        else:
            _id = str(f)

        if len(prop_contexts) > 0:
            prop_contexts["@vocab"] = "https://purl.org/geojson/vocab#"
            props["@context"] = prop_contexts

        fs.append(FeatureCollection([], id=_id, metadata=props, **extras))
    return fs

def get_converted_features(g: Graph, fc: Optional[URIRef] = None, iri2id: Optional[Callable[[URIRef], str]] = None) -> List[Feature]:
    fs = []
    if fc is not None:
        feature_finder = g.objects(fc, RDFS.member)
    else:
        feature_finder = g.subjects(RDF.type, GEO.Feature)
    for f in feature_finder:
        # TODO: handle multiple Geometries per Feature
        geoms = []
        props = {}
        extras = {}
        _id = None
        prop_contexts = {}
        associated_observations = set()
        if iri2id is not None:
            _id = iri2id(f)
        for pred, obj in g.predicate_objects(f):
            if pred in (RDFS.label, SKOS.prefLabel) and not "title" in extras:
                extras["title"] = str(obj)
            elif pred in [GEO.hasGeometry, GEO.hasDefaultGeometry]:
                geoms.extend(_extract_geoms(g, pred, obj))
            elif _id is None and pred == DCTERMS.identifier:
                _id = str(obj)  # the identifier becomes the ID
            elif pred == SDO.spatial:
                # The Schema.org version of a GeoSpatial feature
                spatial_node = obj
                spatial_geoms = []
                for p2, o2 in g.predicate_objects(spatial_node):
                    if p2 in [GEO.hasGeometry, GEO.hasDefaultGeometry]:
                        spatial_geoms.extend(_extract_geoms(g, p2, o2))
                if len(spatial_geoms) < 1:
                    # no hasGeometry in the Spatial, treat this as a Feature
                    spatial_geoms.extend(_extract_geoms(g, pred, obj))
                geoms.extend(spatial_geoms)
            elif pred == SOSA.isFeatureOfInterestOf:
                associated_observations.add(obj)
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
            else:
                prefix_pair, name = make_json_key_from_iri(pred, g.namespace_manager)
                if prefix_pair is not None:
                    prefix_ns, prefix_name = prefix_pair
                    prop_contexts[prefix_name] = prefix_ns
                    name = f"{prefix_name}:{name}"
                if isinstance(obj, BNode):
                    props[name] = _extract_bnode(g, obj, prop_contexts)
                else:
                    props[name] = make_json_representation_of_obj(obj)
        # get observations on the feature
        associated_observations = associated_observations.union(set(g.subjects(SOSA.hasFeatureOfInterest, f)))

        if len(associated_observations) > 0:
            props["sosa:isFeatureOfInterestOf"] = obs_dict_list = []
            for obs in associated_observations:
                obs_dict = _extract_observation(g, obs, prop_contexts)
                obs_dict_list.append(obs_dict)
        if _id is not None:
            # ID is not the same as IRI, so put iri in the properties
            props["iri"] = str(f)
        else:
            _id = str(f)
        if len(prop_contexts) > 0:
            prop_contexts["@vocab"] = "https://purl.org/geojson/vocab#"
            props["@context"] = prop_contexts
        if geoms:
            fs.append(Feature(_id, geometry=geoms[0], properties=props, **extras))
    return fs


def convert(g: Graph, do_validate: bool = True, iri2id: Optional[Callable[[URIRef], str]] = None) -> GeoJSON:
    if do_validate:
        # validate the RDF data according to GeoSPARQL
        conforms, results_graph, results_text = validate(
            g,
            shacl_graph=get_geosparql_validator(),
        )
        if not conforms:
            print(results_text)
            return {}
    feature_collections = get_features_collections(g, iri2id=iri2id)
    fc = None
    if len(feature_collections) > 1:
        # A GeoJSON doc can handle maximum of one Feature Collection
        fc = feature_collections[0]
    elif len(feature_collections) == 1:
        fc = feature_collections[0]

    if fc is not None:
        if "metadata" in fc and "iri" in fc["metadata"]:
            fc_iri = fc["metadata"]["iri"]
        else:
            fc_iri = fc["id"]
        features = get_converted_features(g, URIRef(fc_iri), iri2id=iri2id)
        if len(features) > 0:
            fc["features"].extend(features)
        return fc
    else:
        features = get_converted_features(g, iri2id=iri2id)
        if len(features) > 1:
            return FeatureCollection(features)
        elif len(features) == 1:
            return features[0]
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
