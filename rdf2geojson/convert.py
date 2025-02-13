from __future__ import annotations

from collections import defaultdict
from typing import List, Union, Optional, Tuple, Dict, Any, Callable

from rdflib import BNode, Graph, Literal, URIRef, DCTERMS, SOSA, XSD, SKOS, TIME
from rdflib.namespace import GEO, RDF, RDFS, SDO, Namespace, NamespaceManager

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
from .time_ont import temporal_to_string

TERN = Namespace("https://w3id.org/tern/ontologies/tern/")
PREZ = Namespace("https://prez.dev/")
SCHEMA = SDO


def get_geosparql_validator() -> Graph:
    try:
        return Graph().parse(Path(__file__).parent / "geosparql-validator.ttl")
    except FileNotFoundError:
        return Graph().parse(
            "https://raw.githubusercontent.com/opengeospatial/ogc-geosparql/master/vocabularies/validator.ttl"
        )


def make_json_key_from_iri(
    iri: URIRef, ns: NamespaceManager
) -> Tuple[Optional[Tuple[str, str]], str]:
    """
    Returns either a tuple of ((namespace, prefix), local_name)
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
        if (
            obj.datatype is None
            or obj.datatype == XSD.string
            or obj.datatype == RDF.langString
            or obj.language is not None
        ):
            return str(obj)
        elif obj.datatype in (XSD.date, XSD.dateTime, XSD.time):
            try:
                return obj.value.isoformat()
            except Exception:
                return str(obj)
        elif obj.value is not None:
            return obj.value  # This can be a number, decimal, true, false, etc
        else:
            if obj.datatype is not None:
                return {"datatype": str(obj.datatype), "value": str(obj)}
            else:
                return str(
                    obj
                )  # This will be one of our custom datatypes (eg, waMuseumID)


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
    property_ids = list(g.objects(obj, SCHEMA.propertyID))
    if len(property_ids) > 0:
        key_name = URIRef(property_ids[0])
    if key_name is None:
        text_names = list(g.objects(obj, SCHEMA.name))
        if len(text_names) > 0:
            key_name = str(text_names[0])
    sdo_values_list = list(g.objects(obj, SCHEMA.value))
    if len(sdo_values_list) > 0:
        value = sdo_values_list[0]
    if value is None:
        rdf_values_list = list(g.objects(obj, RDF.value))
        if len(rdf_values_list) > 0:
            value = rdf_values_list[0]
    if key_name is None:
        key_name = "error: Could not find a name or propertyID for additionalProperty"
    elif value is None:
        value = "error: Could not find a value for additionalProperty"
    return key_name, value


def _extract_bnode(g: Graph, bn: BNode, prop_contexts: Dict|None=None, recurse: int = 0) -> Dict:
    obs_dict = {}
    for pred, obj in g.predicate_objects(bn):
        prefix_pair, name = make_json_key_from_iri(pred, g.namespace_manager)
        if prop_contexts is not None and prefix_pair is not None:
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


def _extract_observation(g: Graph, obs: URIRef | BNode, prop_contexts: Dict) -> Dict:
    if isinstance(obs, URIRef):
        obs_dict = {"rdf:subject": str(obs)}
    else:
        obs_dict = {}
    members = []  # this could be an observationCollection too
    attribute_list = []
    attribute_list_name = "attributes"
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
                        # conflicting prefix with one that's already in there
                        use_prefix = False
                if use_prefix:
                    prop_contexts[prefix_name] = prefix_ns
                    name = f"{prefix_name}:{name}"
            if pred == TERN.hasAttribute:
                attribute_list_name = name
                attribute_list.append(_extract_attribute(g, obj, prop_contexts))
            elif isinstance(obj, BNode):
                obs_dict[name] = _extract_bnode(g, obj, prop_contexts)
            else:
                obs_dict[name] = make_json_representation_of_obj(obj)
        if len(attribute_list) > 0:
            obs_dict[attribute_list_name] = attribute_list
        if len(members) > 0:
            obs_dict["sosa:hasMember"] = members
    return obs_dict


def _extract_attribute(
    g: Graph, attr: URIRef | BNode, prop_contexts: Dict
) -> Dict | URIRef:
    pred_ob_list = list(g.predicate_objects(attr))
    if isinstance(attr, URIRef):
        if len(pred_ob_list) == 0:
            return attr
        obs_dict = {"rdf:subject": str(attr)}
    else:
        obs_dict = {}

    for pred, obj in pred_ob_list:
        prefix_pair, name = make_json_key_from_iri(pred, g.namespace_manager)
        if prefix_pair is not None:
            use_prefix = True
            prefix_ns, prefix_name = prefix_pair
            if prefix_name in prop_contexts:
                if prefix_ns != prop_contexts[prefix_name]:
                    # conflicting prefix with one that's already in there
                    use_prefix = False
            if use_prefix:
                prop_contexts[prefix_name] = prefix_ns
                name = f"{prefix_name}:{name}"
        if pred == TERN.hasValue:
            obs_dict[name] = _extract_attribute_value(g, obj, prop_contexts)
        elif isinstance(obj, BNode):
            obs_dict[name] = _extract_bnode(g, obj, prop_contexts)
        else:
            obs_dict[name] = make_json_representation_of_obj(obj)
    return obs_dict

def _extract_attribute_value(
    g: Graph, attr: URIRef | BNode, prop_contexts: Dict
) -> Dict | URIRef:
    pred_ob_list = list(g.predicate_objects(attr))
    if isinstance(attr, URIRef):
        if len(pred_ob_list) == 0:
            return attr
        obs_dict = {"rdf:subject": str(attr)}
    else:
        obs_dict = {}

    for pred, obj in pred_ob_list:
        prefix_pair, name = make_json_key_from_iri(pred, g.namespace_manager)
        if prefix_pair is not None:
            use_prefix = True
            prefix_ns, prefix_name = prefix_pair
            if prefix_name in prop_contexts:
                if prefix_ns != prop_contexts[prefix_name]:
                    # conflicting prefix with one that's already in there
                    use_prefix = False
            if use_prefix:
                prop_contexts[prefix_name] = prefix_ns
                name = f"{prefix_name}:{name}"

        if isinstance(obj, BNode):
            obs_dict[name] = _extract_bnode(g, obj, prop_contexts)
        else:
            obs_dict[name] = make_json_representation_of_obj(obj)
    return obs_dict

def _hoist_attribute(
    g: Graph, attr: URIRef | BNode
) -> Dict | URIRef:
    prez_labels = list(g.objects(attr, PREZ.label))
    use_label = None
    if len(prez_labels) > 0:
        use_label = str(prez_labels[0])
    if use_label is None:
        pref_labels = list(g.objects(attr, SKOS.prefLabel))
        if len(pref_labels) > 0:
            use_label = str(pref_labels[0])
    if use_label is None:
        dcterms_labels = list(g.objects(attr, DCTERMS.title))
        if len(dcterms_labels) > 0:
            use_label = str(dcterms_labels[0])
    if use_label is None:
        rdfs_labels = list(g.objects(attr, RDFS.label))
        if len(rdfs_labels) > 0:
            use_label = str(rdfs_labels[0])
    if use_label is None:
        # TODO: What do we actually use as the name of the attribute?
        if isinstance(attr, URIRef):
            use_label = str(attr).rsplit("/", 1)[-1]
        else:
            attrib_links = list(g.objects(attr, TERN.attribute))
            if attrib_links:
                use_label = str(attrib_links[0]).rsplit("/", 1)[-1]
            else:
                use_label = "attribute_"+str(attr).rsplit(":", 1)[-1]
    use_value = None
    prez_values = list(g.objects(attr, PREZ.value))
    if len(prez_values) > 0:
        use_value = str(prez_values[0])
    if use_value is None:
        simple_values = list(g.objects(attr, TERN.hasSimpleValue))
        if len(simple_values) > 0:
            use_value = str(simple_values[0])
    if use_value is None:
        value_links = list(g.objects(attr, TERN.hasValue))
        if len(value_links) > 0:
            tern_value = value_links[0]
            tern_value_rdf_values = list(g.objects(tern_value, RDF.value))
            if len(tern_value_rdf_values) > 0:
                use_value = str(tern_value_rdf_values[0])
            else:
                if isinstance(tern_value, URIRef):
                    use_value = str(tern_value).rsplit("/", 1)[-1]
                else:
                    use_value = "value_"+str(tern_value).rsplit(":", 1)[-1]
    return {use_label: use_value}

def _hoist_observation(
    g: Graph, observation: URIRef | BNode
) -> Dict | URIRef:
    members = list(g.objects(observation, SOSA.hasMember))
    if len(members) > 0:
        # This is an ObservationCollection
        members_results = {}
        for m in members:
            members_results.update(_hoist_observation(g, m))
        return members_results
    prez_labels = list(g.objects(observation, PREZ.label))
    use_label = None
    degraded_label = False

    if len(prez_labels) > 0:
        use_label = str(prez_labels[0])
    if use_label is None:
        pref_labels = list(g.objects(observation, SKOS.prefLabel))
        if len(pref_labels) > 0:
            use_label = str(pref_labels[0])
    if use_label is None:
        dcterms_labels = list(g.objects(observation, DCTERMS.title))
        if len(dcterms_labels) > 0:
            use_label = str(dcterms_labels[0])
    if use_label is None:
        rdfs_labels = list(g.objects(observation, RDFS.label))
        if len(rdfs_labels) > 0:
            use_label = str(rdfs_labels[0])
    if use_label is None:
        observed_properties = list(g.objects(observation, SOSA.observedProperty))
        if len(observed_properties) > 0:
            the_observed_property = observed_properties[0]
            property_labels = list(g.objects(the_observed_property, PREZ.label))
            if len(property_labels) > 0:
                use_label = str(property_labels[0])
            else:
                property_pref_labels = list(g.objects(the_observed_property, SKOS.prefLabel))
                if len(property_pref_labels) > 0:
                    use_label = str(property_pref_labels[0])
    if use_label is None:
        # TODO: What do we actually use as the name of the observation?
        degraded_label = True
        if isinstance(observation, URIRef):
            use_label = str(observation).rsplit("/", 1)[-1]
        else:
            obp_links = list(g.objects(observation, SOSA.observedPropertyr))
            if obp_links:
                use_label = str(obp_links[0]).rsplit("/", 1)[-1]
            else:
                use_label = "observation_"+str(observation).rsplit(":", 1)[-1]
    use_value = None
    prez_values = list(g.objects(observation, PREZ.value))
    if len(prez_values) > 0:
        use_value = str(prez_values[0])
    if use_value is None:
        simple_values = list(g.objects(observation, SOSA.hasSimpleResult))
        if len(simple_values) > 0:
            use_value = str(simple_values[0])
    if use_value is None:
        result_links = list(g.objects(observation, SOSA.hasResult))
        if len(result_links) > 0:
            sosa_result = result_links[0]
            sosa_result_rdfs_labels = list(g.objects(sosa_result, RDFS.label))
            if len(sosa_result_rdfs_labels) > 0:
                use_value = str(sosa_result_rdfs_labels[0])
            else:
                sosa_result_rdf_values = list(g.objects(sosa_result, RDF.value))
                if len(sosa_result_rdf_values) > 0:
                    use_value = str(sosa_result_rdf_values[0])
                else:
                    if degraded_label:
                        # label is bad, and also value is bad. Just give up on this one.
                        return {}
                    if isinstance(sosa_result, URIRef):
                        use_value = str(sosa_result).rsplit("/", 1)[-1]
                    else:
                        use_value = "result_"+str(sosa_result).rsplit(":", 1)[-1]
        else:
            if degraded_label:
                # label is bad, and also value is bad. Just give up on this one.
                return {}
            use_value = "result_"+str(observation)
    return {use_label: use_value}

def get_features_collections(
    g: Graph, iri2id: Optional[Callable[[URIRef], str]] = None
) -> List[FeatureCollection]:
    feature_finder = g.subjects(RDF.type, GEO.FeatureCollection)
    fs = []
    for f in feature_finder:
        props = {}
        extras = {}
        prop_contexts = {}
        _id = None
        if iri2id is not None:
            _id = iri2id(f)
        else:
            if "#" in str(f):
                _id = str(f).rsplit("#", 1)[-1]
            else:
                _id = str(f).rsplit("/", 1)[-1]
        for pred, obj in g.predicate_objects(f):
            if pred in (RDFS.label, SKOS.prefLabel) and "title" not in extras:
                extras["title"] = str(obj)
            elif pred == RDFS.member:
                # Skip the members, they are handled by get_features
                pass
            elif pred == SCHEMA.additionalProperty:
                # This is the Schema.org version of a Key-Value pair
                p_key, p_value = _extract_additional_property(g, pred, obj)
                if isinstance(p_key, URIRef):
                    prefix_pair, name = make_json_key_from_iri(
                        URIRef(p_key), g.namespace_manager
                    )
                    if prefix_pair is not None:
                        prefix_ns, prefix_name = prefix_pair
                        prop_contexts[prefix_name] = prefix_ns
                        p_key = f"{prefix_name}:{name}"
                    else:
                        p_key = str(p_key)
                props[p_key] = p_value
                continue
            prefix_pair, name = make_json_key_from_iri(pred, g.namespace_manager)
            if prefix_pair is not None:
                prefix_ns, prefix_name = prefix_pair
                prop_contexts[prefix_name] = prefix_ns
                name = f"{prefix_name}:{name}"
            if isinstance(obj, BNode):
                props[name] = _extract_bnode(g, obj, prop_contexts)
            else:
                props[name] = make_json_representation_of_obj(obj)

        # ID is not the same as IRI, so put iri in the properties
        props["rdf:subject"] = str(f)

        if len(prop_contexts) > 0:
            prop_contexts["@vocab"] = "https://purl.org/geojson/vocab#"
            props["@context"] = prop_contexts

        fs.append(FeatureCollection([], id=_id, metadata=props, **extras))
    return fs


def get_converted_features(
    g: Graph,
    fc: Optional[URIRef] = None,
    iri2id: Optional[Callable[[URIRef], str]] = None,
) -> List[Feature]:
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
        attribute_list = []
        _id = None
        prop_contexts = {}
        associated_observations = set()
        if iri2id is not None:
            _id = iri2id(f)
        else:
            if "#" in str(f):
                _id = str(f).rsplit("#", 1)[-1]
            else:
                _id = str(f).rsplit("/", 1)[-1]
        for pred, obj in g.predicate_objects(f):
            if pred in (RDFS.label, SKOS.prefLabel) and "title" not in extras:
                extras["title"] = str(obj)
                continue
            elif pred in [GEO.hasGeometry, GEO.hasDefaultGeometry]:
                geoms.extend(_extract_geoms(g, pred, obj))
                continue
            elif pred == SCHEMA.spatial:
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
                continue
            elif pred == SOSA.isFeatureOfInterestOf:
                associated_observations.add(obj)
                continue
            elif pred == SCHEMA.additionalProperty:
                # This is the Schema.org version of a Key-Value pair
                p_key, p_value = _extract_additional_property(g, pred, obj)
                if isinstance(p_key, URIRef):
                    prefix_pair, name = make_json_key_from_iri(
                        URIRef(p_key), g.namespace_manager
                    )
                    if prefix_pair is not None:
                        prefix_ns, prefix_name = prefix_pair
                        prop_contexts[prefix_name] = prefix_ns
                        p_key = f"{prefix_name}:{name}"
                    else:
                        p_key = str(p_key)
                props[p_key] = p_value
                continue
            elif pred == TERN.hasAttribute:
                # This is the TERN version of a Key-Value pair
                attribute_list.append(_extract_attribute(g, obj, prop_contexts))
                continue

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
        associated_observations = associated_observations.union(
            set(g.subjects(SOSA.hasFeatureOfInterest, f))
        )

        if len(associated_observations) > 0:
            props["sosa:isFeatureOfInterestOf"] = obs_dict_list = []
            for obs in associated_observations:
                obs_dict = _extract_observation(g, obs, prop_contexts)
                obs_dict_list.append(obs_dict)
        if len(attribute_list) > 0:
            props["tern:hasAttribute"] = attribute_list

        # ID is not the same as IRI, so put iri in the properties
        props["rdf:subject"] = str(f)

        if len(prop_contexts) > 0:
            prop_contexts["@vocab"] = "https://purl.org/geojson/vocab#"
            props["@context"] = prop_contexts
        if geoms:
            fs.append(Feature(_id, geometry=geoms[0], properties=props, **extras))
    return fs


def get_converted_features_for_human(
    g: Graph,
    fc: Optional[URIRef] = None,
    iri2id: Optional[Callable[[URIRef], str]] = None,
) -> List[Feature]:
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
        known_time_strings = []
        attribute_dict = defaultdict(list)
        observations_dict = defaultdict(list)
        additional_properties_dict = defaultdict(list)
        associated_observations = set()
        _id = None
        if iri2id is not None:
            _id = iri2id(f)
        else:
            if "#" in str(f):
                _id = str(f).rsplit("#", 1)[-1]
            else:
                _id = str(f).rsplit("/", 1)[-1]
        for pred, obj in g.predicate_objects(f):
            if pred in (RDFS.label, SKOS.prefLabel) and "title" not in extras:
                extras["title"] = str(obj)
                continue
            elif pred in [GEO.hasGeometry, GEO.hasDefaultGeometry]:
                geoms.extend(_extract_geoms(g, pred, obj))
                continue
            elif pred == SCHEMA.spatial:
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
                continue
            elif pred == SOSA.isFeatureOfInterestOf:
                associated_observations.add(obj)
                continue
            elif pred == SCHEMA.additionalProperty:
                # This is the Schema.org version of a Key-Value pair
                add_prop_key, add_prop_val = _extract_additional_property(g, pred, obj)
                additional_properties_dict[add_prop_key].append(str(add_prop_val))
                continue
            elif pred == TERN.hasAttribute:
                # This is the TERN version of a Key-Value pair
                _hoisted_attribute_dict = _hoist_attribute(g, obj)
                for (k, v) in _hoisted_attribute_dict.items():
                    attribute_dict[k].append(v)
                continue
            elif pred == TIME.hasTime:
                # make it a time string
                known_time_strings.append(temporal_to_string(g, obj))
                continue
            elif pred == SCHEMA.temporal:
                known_time_strings.append(temporal_to_string(g, obj))
                continue
            prefix_pair, name = make_json_key_from_iri(pred, g.namespace_manager)

            if isinstance(obj, BNode):
                props[name] = _extract_bnode(g, obj)
            elif isinstance(obj, URIRef):
                has_obj_string = None
                prez_value = list(g.objects(obj, PREZ.value))
                if len(prez_value) > 0:
                    has_obj_string = str(prez_value[0])
                if has_obj_string is None:
                    prez_label = list(g.objects(obj, PREZ.label))
                    if len(prez_label) > 0:
                        has_obj_string = str(prez_label[0])
                if has_obj_string is None:
                    skos_preflabels = list(g.objects(obj, SKOS.prefLabel))
                    if len(skos_preflabels) > 0:
                        has_obj_string = str(skos_preflabels[0])
                if has_obj_string is None:
                    dcterms_titles = list(g.objects(obj, DCTERMS.title))
                    if len(dcterms_titles) > 0:
                        has_obj_string = str(dcterms_titles[0])
                if has_obj_string is None:
                    rdfs_labels = list(g.objects(obj, RDFS.label))
                    if len(rdfs_labels) > 0:
                        has_obj_string = str(rdfs_labels[0])
                if has_obj_string is None:
                    # TODO: What do we actually use as the value of the property?
                    has_obj_string = str(obj)
                props[name] = has_obj_string
            else:
                props[name] = make_json_representation_of_obj(obj)
        if len(known_time_strings) > 1:
            props["datetime"] = "; ".join(known_time_strings)
        else:
            props["datetime"] = known_time_strings[0]
        # get observations on the feature
        associated_observations = associated_observations.union(
            set(g.subjects(SOSA.hasFeatureOfInterest, f))
        )
        if len(associated_observations) > 0:
            for obs in associated_observations:
                hoisted_dict = _hoist_observation(g, obs)
                for (k, v) in hoisted_dict.items():
                    observations_dict[k].append(v)

        for (obs_key, obs_value) in observations_dict.items():
            if obs_key not in props:
                if len(obs_value) > 1:
                    props[obs_key] = "; ".join(obs_value)
                else:
                    props[obs_key] = obs_value[0]
        for (attr_key, attr_value) in attribute_dict.items():
            if attr_key not in props:
                if len(attr_value) > 1:
                    props[attr_key] = "; ".join(attr_value)
                else:
                    props[attr_key] = attr_value[0]
        for (add_key, add_value) in additional_properties_dict.items():
            if add_key not in props:
                if len(add_value) > 1:
                    props[add_key] = "; ".join(add_value)
                else:
                    props[add_key] = add_value[0]

        if geoms:
            fs.append(Feature(_id, geometry=geoms[0], properties=props, **extras))
    return fs

def convert(
    g: Graph, do_validate: bool = True, iri2id: Optional[Callable[[URIRef], str]] = None,
    kind: str = "machine", collection_label: str|None = None,
) -> GeoJSON:
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
        if "metadata" in fc and "rdf:subject" in fc["metadata"]:
            fc_iri = fc["metadata"]["rdf:subject"]
        else:
            fc_iri = fc["id"]
        if kind == "human":
            features = get_converted_features_for_human(g, URIRef(fc_iri), iri2id=iri2id)
        else:
            features = get_converted_features(g, URIRef(fc_iri), iri2id=iri2id)
        if len(features) > 0:
            fc["features"].extend(features)
        return fc
    else:
        if kind == "human":
            features = get_converted_features_for_human(g, iri2id=iri2id)
        else:
            features = get_converted_features(g, iri2id=iri2id)
        if len(features) > 1:
            # Make a new feature collection for these Features.
            if collection_label is not None:
                return FeatureCollection(features, title=collection_label)
            else:
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
