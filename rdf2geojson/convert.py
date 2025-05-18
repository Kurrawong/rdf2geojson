from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping
from functools import lru_cache
from typing import Union, Optional, Any, Callable

from rdflib import BNode, Graph, Literal, URIRef
from rdflib.namespace import GEO, RDF, RDFS, SDO, DCTERMS, SOSA, XSD, SKOS, TIME, Namespace, NamespaceManager

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
GEO_Feature = GEO.Feature
PrezFocusNode = PREZ.FocusNode
PrezType = PREZ.type
PrezLabel = PREZ.label
RDFType = RDF.type

observation_temporal_predicates = [SCHEMA.temporal, SOSA.phenomenonTime, TERN.resultDateTime, SOSA.resultTime]

def get_geosparql_validator() -> Graph:
    try:
        return Graph().parse(Path(__file__).parent / "geosparql-validator.ttl")
    except FileNotFoundError:
        return Graph().parse(
            "https://raw.githubusercontent.com/opengeospatial/ogc-geosparql/master/vocabularies/validator.ttl"
        )


def make_json_key_from_iri(
    iri: URIRef, ns: NamespaceManager
) -> tuple[Optional[tuple[str, str]], str]:
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
    if prefix is not None and namespace is not None:
        return (namespace, prefix), name
    return None, name


def make_json_representation_of_obj(g: Graph, obj: Union[Literal, URIRef], flatten: bool = False)\
        -> Union[Mapping[Any, Any], list[Any], str, int, float]:
    if isinstance(obj, Literal):
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
            # This will be one of our custom datatypes (eg, waMuseumID)
            if obj.datatype is not None:
                if flatten:
                    dt_label = _get_annotation_label(g, obj.datatype)
                    if dt_label is None:
                        dt_label = str(obj.datatype).rsplit("/", 1)[-1]
                    return f"{str(obj)} ({dt_label})"
                else:
                    return {"datatype": str(obj.datatype), "value": str(obj)}
            else:
                return str(obj)
    return str(obj)

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


def _extract_geoms(g: Graph, pred, obj) -> list:
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


def _extract_additional_property(g: Graph, pred, obj) -> tuple[Union[str, URIRef], Any]:
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


def _extract_bnode(g: Graph, bn: BNode, prop_contexts: dict|None=None, recurse: int = 0) -> dict:
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
            obs_dict[name] = make_json_representation_of_obj(g, obj)
    return obs_dict


def _extract_observation(g: Graph, obs: URIRef | BNode, prop_contexts: dict) -> dict:
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
                obs_dict[name] = make_json_representation_of_obj(g, obj)
        if len(attribute_list) > 0:
            obs_dict[attribute_list_name] = attribute_list
        if len(members) > 0:
            obs_dict["sosa:hasMember"] = members
    return obs_dict


def _extract_attribute(
    g: Graph, attr: URIRef | BNode, prop_contexts: dict
) -> dict | URIRef:
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
            obs_dict[name] = make_json_representation_of_obj(g, obj)
    return obs_dict

def _extract_attribute_value(
    g: Graph, attr: URIRef | BNode, prop_contexts: dict
) -> dict | URIRef:
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
            obs_dict[name] = make_json_representation_of_obj(g, obj)
    return obs_dict

def _get_annotation_label(g: Graph, labelled_node: Union[URIRef,BNode]) -> Union[str, None]:
    prez_labels = list(g.objects(labelled_node, PREZ.label))
    use_label = None

    if len(prez_labels) > 0:
        use_label = str(prez_labels[0])
    if use_label is None:
        pref_labels = list(g.objects(labelled_node, SKOS.prefLabel))
        if len(pref_labels) > 0:
            use_label = str(pref_labels[0])
    if use_label is None:
        dcterms_labels = list(g.objects(labelled_node, DCTERMS.title))
        if len(dcterms_labels) > 0:
            use_label = str(dcterms_labels[0])
    if use_label is None:
        rdfs_labels = list(g.objects(labelled_node, RDFS.label))
        if len(rdfs_labels) > 0:
            use_label = str(rdfs_labels[0])
    return use_label

def _hoist_attribute(
    g: Graph, attr: URIRef | BNode
) -> dict[str, Any]:
    annotation_label = _get_annotation_label(g, attr)
    attrib_links = list(g.objects(attr, TERN.attribute))
    if len(attrib_links) > 0:
        the_attribute_link = attrib_links[0]
        attrib_labels = list(g.objects(the_attribute_link, PREZ.label))
        if len(attrib_labels) > 0:
            use_label = str(attrib_labels[0])
        else:
            attrib_pref_label = list(g.objects(the_attribute_link, SKOS.prefLabel))
            if len(attrib_pref_label) > 0:
                use_label = str(attrib_pref_label[0])
            elif annotation_label is not None:
                use_label = annotation_label
            else:
                attrib_parts = str(the_attribute_link).rsplit("/", 2)
                if len(attrib_parts) == 1:
                    use_label = attrib_parts[0]
                elif len(attrib_parts) == 2:
                    use_label = attrib_parts[-1]
                else:
                    use_label = attrib_parts[-2]+":"+attrib_parts[-1]
    else:
        use_label = annotation_label

    if use_label is None:
        # TODO: What do we actually use as the name of the attribute?
        if isinstance(attr, URIRef):
            attr_parts = str(attr).rsplit("/", 2)
            if len(attr_parts) == 1:
                use_label = attr_parts[0]
            else:
                if attr_parts[-1].isnumeric():
                    use_label = attr_parts[-2]
                else:
                    use_label = attr_parts[-1]
        else:
            if attrib_links:
                use_label = str(attrib_links[0]).rsplit("/", 1)[-1]
            else:
                use_label = "attribute_"+str(attr).rsplit(":", 1)[-1]
    use_value_node = None
    fallback_value = None
    prez_values = list(g.objects(attr, PREZ.value))
    if len(prez_values) > 0:
        use_value_node = prez_values[0]
    if use_value_node is None:
        simple_values = list(g.objects(attr, TERN.hasSimpleValue))
        if len(simple_values) > 0:
            use_value_node = simple_values[0]
    if use_value_node is None:
        value_links = list(g.objects(attr, TERN.value))
        if len(value_links) > 0:
            tern_value = value_links[0]
            tern_value_rdf_values = list(g.objects(tern_value, RDF.value))
            if len(tern_value_rdf_values) > 0:
                use_value_node = tern_value_rdf_values[0]
            else:
                if isinstance(tern_value, URIRef):
                    fallback_value = str(tern_value).rsplit("/", 1)[-1]
                else:
                    fallback_value = "value_"+str(tern_value).rsplit(":", 1)[-1]
    use_value = None
    if use_value_node is not None:
        if isinstance(use_value_node, (URIRef, BNode)):
            use_annotation_value = _get_annotation_label(g, use_value_node)
            if use_annotation_value is not None:
                use_value = use_annotation_value
            else:
                if fallback_value is None:
                    fallback_value = str(use_value_node)
        else:
            use_value = str(use_value_node)
    if use_value is None:
        if fallback_value is not None:
            use_value = fallback_value
        else:
            use_value = "error: Could not find a value for attribute"
    return {use_label: use_value}

@lru_cache(maxsize=128)
def _get_flattened_observation_collection_properties(g, observation_collection) -> dict[str, Any]:
    time_string: Optional[str] = None
    for tp in observation_temporal_predicates:
        if temporal_matches := list(g.objects(observation_collection, tp)):
            time_string = temporal_to_string(g, temporal_matches[0])
            break

    procedure_string: Optional[str] = None
    procedure_uri: Optional[URIRef|BNode] = None
    if used_procedures := list(g.objects(observation_collection, SOSA.usedProcedure)):
        for used_procedure in used_procedures:
            if procedure_method_types := list(g.objects(used_procedure, TERN.methodType)):
                procedure_uri = procedure_method_types[0]
                break
            else:
                procedure_uri = used_procedure
                break
        if procedure_uri is not None:
            if isinstance(procedure_uri, (URIRef, BNode)):
                if (check_procedure_label := _get_annotation_label(g, procedure_uri)) is not None:
                    procedure_string = check_procedure_label
                else:
                    procedure_string = str(procedure_uri)
            else:
                procedure_string = str(procedure_uri)
    flattened_attributes = {}
    if has_attributes := list(g.objects(observation_collection, TERN.hasAttribute)):
        for has_attribute in has_attributes:
            a_flat = _hoist_attribute(g, has_attribute)
            flattened_attributes[has_attribute] = a_flat

    ret = {}
    if time_string is not None:
        ret["time"] = time_string
    if procedure_uri is not None and procedure_string is not None:
        ret["procedure"] = (procedure_uri, procedure_string)
    if flattened_attributes:
        ret["attributes"] = flattened_attributes
    return ret


# Set this to True, to allow 'time' and 'procedure' flattened properties from the
# obsevation collection to be applied to the Observation, if they are different.
# If they are the same, they are still ignored.
ALLOW_SECOND_ATTRIBUTES=False
OBSERVATION_FLATTEN_EXCLUDE_TYPES = [SOSA.Sampling, TERN.Sampling]
def _hoist_and_flatten_observation(
    g: Graph, observation: URIRef | BNode,
) -> dict[str, Any]:
    if known_types := list(g.objects(observation, RDFType)):
        for kt in known_types:
            if kt in OBSERVATION_FLATTEN_EXCLUDE_TYPES:
                # Don't try to flatten samplings as if they are Samples
                return {}
    if has_children := list(g.objects(observation, SOSA.hasMember)):
        # This is an ObservationCollection, skip it, because observation_colltion properties are collected by each member
        return {"children": has_children}
    
    flattened_collection_props: dict[URIRef|BNode, dict] = {}
    if in_collections := list(g.subjects(SOSA.hasMember, observation)):
        for in_col in in_collections:
            flattened_collection_props[in_col] = _get_flattened_observation_collection_properties(g, in_col)

    degraded_label = False
    annotation_label = _get_annotation_label(g, observation)
    if observed_properties := list(g.objects(observation, SOSA.observedProperty)):
        the_observed_property = observed_properties[0]
        if property_labels := list(g.objects(the_observed_property, PREZ.label)):
            use_label = str(property_labels[0])
        else:
            if property_pref_labels := list(g.objects(the_observed_property, SKOS.prefLabel)):
                use_label = str(property_pref_labels[0])
            elif annotation_label is not None:
                use_label = annotation_label
            else:
                observed_property_parts = str(the_observed_property).rsplit("/", 2)
                if len(observed_property_parts) == 1:
                    use_label = observed_property_parts[0]
                elif len(observed_property_parts) == 2:
                    use_label = observed_property_parts[-1]
                else:
                    use_label = observed_property_parts[-2]+":"+observed_property_parts[-1]
    else:
        use_label = annotation_label

    if use_label is None:
        # TODO: What do we actually use as the name of the observation?
        degraded_label = True
        if isinstance(observation, URIRef):
            observation_parts = str(observation).rsplit("/", 2)
            if len(observation_parts) == 1:
                use_label = observation_parts[0]
            else:
                if observation_parts[-1].isnumeric():
                    use_label = observation_parts[-2]
                else:
                    use_label = observation_parts[-1]
        else:
            if observed_properties:
                use_label = str(observed_properties[0]).rsplit("/", 1)[-1]
            else:
                use_label = "observation_"+str(observation).rsplit(":", 1)[-1]
    use_value_node = None
    fallback_value = None
    if prez_values := list(g.objects(observation, PREZ.value)):
        use_value_node = prez_values[0]
    if use_value_node is None:
        if simple_values := list(g.objects(observation, SOSA.hasSimpleResult)):
            use_value_node = simple_values[0]
    if use_value_node is None:
        if result_links := list(g.objects(observation, SOSA.hasResult)):
            sosa_result = result_links[0]
            if sosa_result_rdfs_labels := list(g.objects(sosa_result, RDFS.label)):
                fallback_value = str(sosa_result_rdfs_labels[0])
            else:
                if sosa_result_rdf_values := list(g.objects(sosa_result, RDF.value)):
                    if isinstance(sosa_result_rdf_values[0], (URIRef, BNode)):
                        use_value_node = sosa_result_rdf_values[0]
                    else:
                        fallback_value = str(sosa_result_rdf_values[0])
                else:
                    if degraded_label:
                        # label is bad, and also value is bad. Just give up on this one.
                        return {}
                    if isinstance(sosa_result, URIRef):
                        fallback_value = str(sosa_result).rsplit("/", 1)[-1]
                    else:
                        fallback_value = "result_"+str(sosa_result).rsplit(":", 1)[-1]
    use_value = None
    if use_value_node is not None:
        if isinstance(use_value_node, (URIRef, BNode)):
            use_annotation_value = _get_annotation_label(g, use_value_node)
            if use_annotation_value is not None:
                use_value = use_annotation_value
            else:
                if fallback_value is None:
                    fallback_value = str(use_value_node)
        else:
            use_value = str(use_value_node)
    if use_value is None:
        if fallback_value is not None:
            use_value = fallback_value
        else:
            if degraded_label:
                # label is bad, and also value is bad. Just give up on this one.
                return {}
            use_value = "result_"+str(observation)
    this_observation_flattened_dict = {use_label: use_value}

    flattened_attributes = {}
    if has_attributes := list(g.objects(observation, TERN.hasAttribute)):
        for has_attribute in has_attributes:
            a_flat = _hoist_attribute(g, has_attribute)
            flattened_attributes[has_attribute] = a_flat
    
    time_string: Optional[str] = None
    for tp in observation_temporal_predicates:
        if temporal_matches := list(g.objects(observation, tp)):
            time_string = temporal_to_string(g, temporal_matches[0])
            break


    procedure_string: Optional[str] = None
    procedure_uri: Optional[URIRef|BNode] = None 
    if used_procedures := list(g.objects(observation, SOSA.usedProcedure)):
        for used_procedure in used_procedures:
            if procedure_method_types := list(g.objects(used_procedure, TERN.methodType)):
                procedure_uri = procedure_method_types[0]
                break
            else:
                procedure_uri = used_procedure
                break
        if procedure_uri is not None:
            if isinstance(procedure_uri, (URIRef, BNode)):
                if (check_procedure_label := _get_annotation_label(g, procedure_uri)) is not None:
                    procedure_string = check_procedure_label
                else:
                    procedure_string = str(procedure_uri)
            else:
                procedure_string = str(procedure_uri)

    if flattened_collection_props:
        for i, (coll_uri, coll_flat) in enumerate(sorted(flattened_collection_props.items())):
            if (coll_time_string := coll_flat.get("time", None)) is not None:
                if time_string is None:
                    time_string = coll_time_string
                elif ALLOW_SECOND_ATTRIBUTES and time_string == coll_time_string:
                    # same, don't duplicate it
                    pass
                elif ALLOW_SECOND_ATTRIBUTES:
                    this_observation_flattened_dict[f"{use_label} ({str(i+1)}) (datetime)"] = time_string
            if (coll_procedure_pair := coll_flat.get("procedure", None)) is not None:
                coll_proc_uri, coll_proc_string = coll_procedure_pair
                if procedure_string is None and procedure_uri is None:
                    procedure_string = coll_proc_string
                    procedure_uri = coll_proc_uri
                elif ALLOW_SECOND_ATTRIBUTES and (procedure_uri == coll_proc_uri or procedure_string == coll_proc_string):
                    # same, pass
                    pass
                elif ALLOW_SECOND_ATTRIBUTES:
                    this_observation_flattened_dict[f"{use_label} ({str(i+1)}) (procedure)"] = coll_proc_string
            if (coll_attrs_pairs := coll_flat.get("attributes", None)) is not None:
                for (coll_attr_uri, flattened_coll_attr_kv) in coll_attrs_pairs.items():
                    if coll_attr_uri not in flattened_attributes:
                        flattened_attributes[coll_attr_uri] = flattened_coll_attr_kv
                    elif ALLOW_SECOND_ATTRIBUTES:
                        for f_attr_k, v in flattened_coll_attr_kv.items():
                            this_observation_flattened_dict[f"{use_label} ({str(i+1)}) ({f_attr_k})"] = v

    if time_string is not None:
        this_observation_flattened_dict[f"{use_label} (datetime)"] = time_string
    if procedure_uri and procedure_string:
        this_observation_flattened_dict[f"{use_label} (procedure)"] = procedure_string
    if flattened_attributes:
        for attr_uri, flattened_attr_kv in flattened_attributes.items():
            for f_attr_k, v in flattened_attr_kv.items():
                this_observation_flattened_dict[f"{use_label} ({f_attr_k})"] = v
    return this_observation_flattened_dict

def get_features_collections(
    g: Graph, fc_uri: Optional[URIRef] = None,
    iri2id: Optional[Callable[[URIRef], str]] = None
) -> list[tuple[URIRef, FeatureCollection]]:
    fc_finder = g.subjects(RDFType, GEO.FeatureCollection)
    fcs = []
    for f in fc_finder:
        if fc_uri is not None:
            # Filter on a given known FeatureCollection URI
            if fc_uri != f:
                continue
        props = {}
        extras = {}
        prop_contexts = {}
        anot = None
        _id = None
        if iri2id is not None:
            _id = iri2id(f)
        else:
            if "#" in str(f):
                _id = str(f).rsplit("#", 1)[-1]
            else:
                _id = str(f).rsplit("/", 1)[-1]
        for pred, obj in g.predicate_objects(f):
            if pred == RDFType:
                if obj in [PrezFocusNode, GEO.FeatureCollection]:
                    # Don't include FocusNode or Geo.FeatureCollection in the list of RDF types, they are both implied
                    continue
            elif pred == PrezType:
                # Don't include PrezType in the list of properties, it is a hidden property
                continue
            elif pred == RDFS.member:
                # Skip the members, they are handled by get_features
                continue
            elif pred == PrezLabel:
                anot = str(obj)
                continue
            elif pred in (RDFS.label, SKOS.prefLabel) and "title" not in extras:
                extras["title"] = str(obj)
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
                props[name] = make_json_representation_of_obj(g, obj)

        # ID is not the same as IRI, so put iri in the properties
        props["rdf:subject"] = str(f)

        if len(prop_contexts) > 0:
            prop_contexts["@vocab"] = "https://purl.org/geojson/vocab#"
            props["@context"] = prop_contexts

        if "title" not in extras and anot is not None:
            extras["title"] = anot

        fcs.append((f, FeatureCollection([], id=_id, metadata=props, **extras)))
    return fcs

def get_features_collections_for_human(
    g: Graph, fc_uri: Optional[URIRef] = None,
    iri2id: Optional[Callable[[URIRef], str]] = None
) -> list[tuple[URIRef, FeatureCollection]]:
    fc_finder = g.subjects(RDFType, GEO.FeatureCollection)
    fcs = []
    for f in fc_finder:
        if fc_uri is not None:
            # Filter on a given known FeatureCollection URI
            if fc_uri != f:
                continue
        props = {}
        extras = {}
        prop_contexts = {}
        anot = None
        attribute_dict = defaultdict(list)
        props_dict_lists = defaultdict(list)
        additional_properties_dict = defaultdict(list)
        _id = None
        if iri2id is not None:
            _id = iri2id(f)
        else:
            if "#" in str(f):
                _id = str(f).rsplit("#", 1)[-1]
            else:
                _id = str(f).rsplit("/", 1)[-1]
        for pred, obj in g.predicate_objects(f):
            if pred == RDFType:
                if obj in [PrezFocusNode, GEO.FeatureCollection]:
                    # Don't include FocusNode or Geo.FeatureCollection in the list of RDF types, they are both implied
                    continue
            elif pred == PrezType:
                # Don't include PrezType in the list of properties, it is a hidden property
                continue
            elif pred == RDFS.member:
                # Skip the members, they are handled by get_features
                continue
            elif pred == PrezLabel:
                anot = str(obj)
                continue
            elif pred in (RDFS.label, SKOS.prefLabel) and "title" not in extras:
                extras["title"] = str(obj)
                continue # This is a difference between semantic and human-readable labels
                # we don't want to duplicate the "title" into the properties, in human-readable mode.
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
            prefix_pair, name = make_json_key_from_iri(pred, g.namespace_manager)
            if isinstance(obj, (URIRef, BNode)):
                bn_has_obj_string = None
                bn_prez_value = list(g.objects(obj, PREZ.value))
                if len(bn_prez_value) > 0:
                    if isinstance(bn_prez_value[0], (URIRef, BNode)):
                        bn_has_obj_string = _get_annotation_label(g, bn_prez_value[0])
                    else:
                        bn_has_obj_string = str(bn_prez_value[0])

                if bn_has_obj_string is None:
                    bn_has_obj_string = _get_annotation_label(g, obj)
                if bn_has_obj_string is None:
                    # TODO: What do we actually use as the value of the property?
                    bn_has_obj_string = str(obj)
                props_dict_lists[name].append(bn_has_obj_string)
            else:
                props_dict_lists[name].append(make_json_representation_of_obj(g, obj, flatten=True))
        for (name, values) in props_dict_lists.items():
            if len(values) > 1:
                props[name] = "; ".join(str(v) for v in values)
            else:
                props[name] = values[0]
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

        if "title" not in extras and anot is not None:
            extras["title"] = anot
        props["uri"] = str(f)
        fcs.append((f, FeatureCollection([], id=_id, metadata=props, **extras)))
    return fcs

def get_converted_features(
    g: Graph,
    fc: Optional[URIRef] = None,
    iri2id: Optional[Callable[[URIRef], str]] = None,
) -> list[Feature]:
    fs = []
    if fc is not None:
        feature_finder = g.objects(fc, RDFS.member)
    else:
        feature_finder = g.subjects(RDFType, GEO_Feature)
    for f in feature_finder:
        # TODO: handle multiple Geometries per Feature
        default_geometry = None
        centroid = None
        bounding_box = None
        geoms: list[list] = []
        props = {}
        extras = {}
        attribute_list = []
        _id = None
        prop_contexts = {}
        associated_observations = set()
        anot = None
        if iri2id is not None:
            _id = iri2id(f)
        else:
            if "#" in str(f):
                _id = str(f).rsplit("#", 1)[-1]
            else:
                _id = str(f).rsplit("/", 1)[-1]
        for pred, obj in g.predicate_objects(f):
            if pred == RDFType:
                if obj in [PrezFocusNode, GEO_Feature]:
                    # Don't include FocusNode or GeoFeature in the list of RDF types, they are both implied
                    continue
            elif pred == PrezType:
                # Don't include PrezType in the list of properties, it is a hidden property
                continue
            elif pred == PrezLabel:
                anot = str(obj)
                continue
            elif pred in (RDFS.label, SKOS.prefLabel) and "title" not in extras:
                extras["title"] = str(obj)
            elif pred == GEO.hasDefaultGeometry:
                default_geometry = _extract_geoms(g, pred, obj)
                continue
            elif pred == GEO.hasBoundingBox:
                bounding_box = _extract_geoms(g, pred, obj)
                continue
            elif pred == GEO.hasCentroid:
                centroid = _extract_geoms(g, pred, obj)
                continue
            elif pred == GEO.hasGeometry:
                geoms.append(_extract_geoms(g, pred, obj))
                continue
            elif pred == SCHEMA.spatial:
                # The Schema.org version of a GeoSpatial feature
                spatial_node = obj
                spatial_geoms: list[list] = []
                sp_default_geometry = None
                sp_bounding_box = None
                sp_centroid = None
                for p2, o2 in g.predicate_objects(spatial_node):
                    if p2 == GEO.hasDefaultGeometry:
                        sp_default_geometry = _extract_geoms(g, p2, o2)
                    elif p2 == GEO.hasBoundingBox:
                        sp_bounding_box = _extract_geoms(g, p2, o2)
                    elif p2 == GEO.hasCentroid:
                        sp_centroid = _extract_geoms(g, p2, o2)
                    elif p2 == GEO.hasGeometry:
                        spatial_geoms.append(_extract_geoms(g, p2, o2))
                if sp_default_geometry is not None and default_geometry is None:
                    default_geometry = sp_default_geometry
                if sp_bounding_box is not None and bounding_box is None:
                    bounding_box = sp_bounding_box
                if sp_centroid is not None and centroid is None:
                    centroid = sp_centroid
                if sp_default_geometry is None and sp_centroid is None and sp_bounding_box is None \
                        and len(spatial_geoms) < 1:
                    # no hasGeometry in the Spatial, treat this as a the geometry itself.
                    spatial_geoms.append(_extract_geoms(g, pred, obj))
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
                props[name] = make_json_representation_of_obj(g, obj)
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
        if "title" not in extras and anot is not None:
            extras["title"] = anot
        if len(prop_contexts) > 0:
            prop_contexts["@vocab"] = "https://purl.org/geojson/vocab#"
            props["@context"] = prop_contexts
        use_geometry = None
        if default_geometry is not None:
            use_geometry = default_geometry[0]
        elif len(geoms) > 0:
            use_geometry = geoms[0][0]
        elif bounding_box is not None:
            use_geometry = bounding_box[0]
        elif centroid is not None:
            use_geometry = centroid[0]
        if use_geometry is not None:
            fs.append(Feature(_id, geometry=use_geometry, properties=props, **extras))
    return fs


def get_converted_features_for_human(
    g: Graph,
    fc: Optional[URIRef] = None,
    iri2id: Optional[Callable[[URIRef], str]] = None,
) -> list[Feature]:
    fs = []
    if fc is not None:
        feature_finder = g.objects(fc, RDFS.member)
    else:
        feature_finder = g.subjects(RDFType, GEO_Feature)
    for f in feature_finder:
        # TODO: handle multiple Geometries per Feature
        default_geometry = None
        centroid = None
        bounding_box = None
        geoms: list[list] = []
        props = {}
        extras = {}
        known_time_strings = []
        attribute_dict = defaultdict(list)
        observations_dict = defaultdict(list)
        additional_properties_dict = defaultdict(list)
        associated_observations = set()
        props_dict_lists = defaultdict(list)
        anot = None
        _id = None
        if iri2id is not None:
            _id = iri2id(f)
        else:
            if "#" in str(f):
                _id = str(f).rsplit("#", 1)[-1]
            else:
                _id = str(f).rsplit("/", 1)[-1]
        for pred, obj in g.predicate_objects(f):
            if pred == RDFType:
                if obj in [PrezFocusNode, GEO_Feature]:
                    # Don't include FocusNode or GeoFeature in the lost of RDF types, they are both implied
                    continue
            elif pred == PrezType:
                # Don't include PrezType in the list of properties, it is a hidden property
                continue
            elif pred == PrezLabel:
                anot = str(obj)
                continue
            elif pred in (RDFS.label, SKOS.prefLabel):
                if "title" not in extras:
                    extras["title"] = str(obj)
                additional_properties_dict["label"].append(str(obj))
                continue 
            elif pred == GEO.hasDefaultGeometry:
                default_geometry = _extract_geoms(g, pred, obj)
                continue
            elif pred == GEO.hasBoundingBox:
                bounding_box = _extract_geoms(g, pred, obj)
                continue
            elif pred == GEO.hasCentroid:
                centroid = _extract_geoms(g, pred, obj)
                continue
            elif pred == GEO.hasGeometry:
                geoms.append(_extract_geoms(g, pred, obj))
                continue
            elif pred == SCHEMA.spatial:
                # The Schema.org version of a GeoSpatial feature
                spatial_node = obj
                spatial_geoms: list[list] = []
                sp_default_geometry = None
                sp_bounding_box = None
                sp_centroid = None
                for p2, o2 in g.predicate_objects(spatial_node):
                    if p2 == GEO.hasDefaultGeometry:
                        sp_default_geometry = _extract_geoms(g, p2, o2)
                    elif p2 == GEO.hasBoundingBox:
                        sp_bounding_box = _extract_geoms(g, p2, o2)
                    elif p2 == GEO.hasCentroid:
                        sp_centroid = _extract_geoms(g, p2, o2)
                    elif p2 == GEO.hasGeometry:
                        spatial_geoms.append(_extract_geoms(g, p2, o2))
                if sp_default_geometry is not None and default_geometry is None:
                    default_geometry = sp_default_geometry
                if sp_bounding_box is not None and bounding_box is None:
                    bounding_box = sp_bounding_box
                if sp_centroid is not None and centroid is None:
                    centroid = sp_centroid
                if sp_default_geometry is None and sp_centroid is None and sp_bounding_box is None \
                        and len(spatial_geoms) < 1:
                    # no hasGeometry in the Spatial, treat this as a the geometry itself.
                    spatial_geoms.append(_extract_geoms(g, pred, obj))
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

            if isinstance(obj, (URIRef, BNode)):
                bn_has_obj_string = None
                bn_prez_value = list(g.objects(obj, PREZ.value))
                if len(bn_prez_value) > 0:
                    if isinstance(bn_prez_value[0], (URIRef, BNode)):
                        bn_has_obj_string = _get_annotation_label(g, bn_prez_value[0])
                    else:
                        bn_has_obj_string = str(bn_prez_value[0])

                if bn_has_obj_string is None:
                    bn_has_obj_string = _get_annotation_label(g, obj)
                if bn_has_obj_string is None:
                    # TODO: What do we actually use as the value of the property?
                    bn_has_obj_string = str(obj)
                props_dict_lists[name].append(bn_has_obj_string)
            else:
                props_dict_lists[name].append(make_json_representation_of_obj(g, obj, flatten=True))
        for (name, values) in props_dict_lists.items():
            if len(values) > 1:
                props[name] = "; ".join(str(v) for v in values)
            else:
                props[name] = values[0]
        if "datetime" not in props and known_time_strings:
            if len(known_time_strings) > 1:
                props["datetime"] = "; ".join(known_time_strings)
            else:
                props["datetime"] = known_time_strings[0]
        # get observations on the feature
        associated_observations = associated_observations.union(
            set(g.subjects(SOSA.hasFeatureOfInterest, f))
        )
        collection_hoisted_observations: dict[URIRef|BNode, dict] = defaultdict(dict)
        hoisted_obs_have_collections: dict[URIRef|BNode, list] = defaultdict(list)
        all_hoisted_observations: dict = {}
        if associated_observations:
            for obs in associated_observations:
                if obs in all_hoisted_observations:
                    continue
                hoisted_observation_dict = _hoist_and_flatten_observation(g, obs)
                if "children" in hoisted_observation_dict:
                    # This is an observation collection.
                    all_hoisted_observations[obs] = {}
                    for obs_ch in hoisted_observation_dict["children"]:
                        if obs_ch in all_hoisted_observations:
                            hoisted_observation_dict = all_hoisted_observations[obs_ch]
                        else:
                            hoisted_observation_dict = _hoist_and_flatten_observation(g, obs_ch)
                            all_hoisted_observations[obs_ch] = hoisted_observation_dict
                        collection_hoisted_observations[obs][obs_ch] = hoisted_observation_dict
                        hoisted_obs_have_collections[obs_ch].append(obs)
                else:
                    all_hoisted_observations[obs] = hoisted_observation_dict
        # get the observations that are not part of any collections, these get added directly to the occurrence
        direct_hoisted: list = [k for k in all_hoisted_observations.keys() if k not in hoisted_obs_have_collections]
        direct_hoisted_keys: dict[str, Any] = defaultdict(list)
        _ = [direct_hoisted_keys[str_key].append(k) for k in direct_hoisted for str_key in all_hoisted_observations[k].keys()]
        for d in direct_hoisted:
            hoisted_observation_dict = all_hoisted_observations[d]
            for (k, v) in hoisted_observation_dict.items():
                observations_dict[k].append(v)
        collection_hoisted_keys: dict[str, Any] = defaultdict(list)
        _ = [collection_hoisted_keys[str_key].append(col_hoisted_child) for coll_hoisted_children in collection_hoisted_observations.values() for col_hoisted_child, coll_hoisted_kv in coll_hoisted_children.items() for str_key in coll_hoisted_kv.keys()]
        dumped_hoisted_obs = set()
        for i, coll_hoisted_k in enumerate(sorted(collection_hoisted_observations)):
            for the_hoisted_obs, the_hoisted_obs_kv in collection_hoisted_observations[coll_hoisted_k].items():
                if the_hoisted_obs in dumped_hoisted_obs:
                    # Don't add same instance of KVs
                    continue
                if len(collection_hoisted_observations) > 1:
                    # This feature is subject of multiple ObservationCollections.
                    # So we must always add grouping tags on them
                    in_grouping_tag = True
                elif not direct_hoisted_keys:
                    # There are no direct keys, so all from group-1 can be ungrouped.
                    in_grouping_tag = False
                else:
                    for (k, v) in the_hoisted_obs_kv.items():
                        if k in direct_hoisted_keys:
                            in_grouping_tag = True
                            break
                    else:
                        in_grouping_tag = False
                for (k, v) in the_hoisted_obs_kv.items():
                    if in_grouping_tag:
                        # An observation of this same flattened key already exists
                        # But we want to deliberately _not_ combine the values.
                        #if k not in direct_hoisted and
                        use_obs_dict_key = f"(collection {i+1}) {k}"
                    else:
                        use_obs_dict_key = k

                    if v not in observations_dict[use_obs_dict_key]:
                        observations_dict[use_obs_dict_key].append(v)
                dumped_hoisted_obs.add(the_hoisted_obs)

        for (obs_key, obs_values) in observations_dict.items():
            if obs_key not in props:
                if len(obs_values) > 1:
                    props[obs_key] = "; ".join(obs_values)
                else:
                    props[obs_key] = obs_values[0]
        for (attr_key, attr_value) in attribute_dict.items():
            if attr_key not in props:
                if len(attr_value) > 1:
                    props[attr_key] = "; ".join(attr_value)
                else:
                    props[attr_key] = attr_value[0]
        # additional_properties_dict is just like props_dict_lists, except its
        # added after observatons, and attributes, and only added if the key doesn't not already exist.
        for (add_key, add_value) in additional_properties_dict.items():
            if add_key not in props:
                if len(add_value) > 1:
                    props[add_key] = "; ".join(add_value)
                else:
                    props[add_key] = add_value[0]
        if "title" not in extras and anot is not None:
            extras["title"] = anot
        props["uri"] = str(f)
        use_geometry = None
        if default_geometry is not None:
            use_geometry = default_geometry[0]
        elif len(geoms) > 0:
            use_geometry = geoms[0][0]
        elif bounding_box is not None:
            use_geometry = bounding_box[0]
        elif centroid is not None:
            use_geometry = centroid[0]
        if use_geometry is not None:
            fs.append(Feature(_id, geometry=use_geometry, properties=props, **extras))

    return fs

def convert(
    g: Graph, do_validate: bool = True, iri2id: Optional[Callable[[URIRef], str]] = None,
    kind: str = "machine", fc_uri: Optional[URIRef] = None, collection_label: Optional[str] = None,
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

    if fc_uri is None and collection_label is not None:
        # When a collection_label is passed, this is a custom collection that doesn't
        # exist as a defined FeatureCollection in the graph. So don't look up FeatureCollections.
        if kind == "human":
            features = get_converted_features_for_human(g, iri2id=iri2id)
        else:
            features = get_converted_features(g, iri2id=iri2id)
        return FeatureCollection(features, title=collection_label)

    if kind == "human":
        feature_collections = get_features_collections_for_human(g, fc_uri=fc_uri, iri2id=iri2id)
    else:
        feature_collections = get_features_collections(g, fc_uri=fc_uri, iri2id=iri2id)
    fc = None
    from_fc_uri: Optional[URIRef] = None
    if len(feature_collections) > 1:
        # A GeoJSON doc can handle maximum of one Feature Collection
        from_fc_uri, fc = feature_collections[0]
    elif len(feature_collections) == 1:
        from_fc_uri, fc = feature_collections[0]
    if from_fc_uri is not None and fc is not None:
        if kind == "human":
            features = get_converted_features_for_human(g, from_fc_uri, iri2id=iri2id)
        else:
            features = get_converted_features(g, from_fc_uri, iri2id=iri2id)
        if len(features) > 0:
            fc["features"].extend(features)
        return fc
    else:
        if kind == "human":
            features = get_converted_features_for_human(g, iri2id=iri2id)
        else:
            features = get_converted_features(g, iri2id=iri2id)
        if (len(features) > 1) or (collection_label is not None):
            # Make a new feature collection for these Features.
            if collection_label is not None:
                return FeatureCollection(features, title=collection_label)
            else:
                return FeatureCollection(features)
        elif len(features) == 1:
            return features[0]
        else:
            return GeoJSON()


def unconvert_geometry(geom: dict) -> tuple[str, str]:
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
        g.add((id_, RDFType, GEO.Feature))
        g.add((bn, RDFType, GEO.Geometry))
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
