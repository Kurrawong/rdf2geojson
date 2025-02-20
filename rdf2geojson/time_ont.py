from __future__ import annotations

from rdflib.namespace import TIME, XSD, RDF
from rdflib.graph import Graph
from rdflib.term import BNode, Literal, URIRef
from typing import Union

def temporal_to_string(graph: Graph, temporal: Union[URIRef,BNode,Literal]) -> str:
    if isinstance(temporal, Literal):
        if temporal.datatype:
            if temporal.datatype == XSD.string or temporal.datatype == RDF.langString:
                return str(temporal)
            elif temporal.datatype == XSD.dateTime or temporal.datatype == XSD.dateTimeStamp:
                return str(temporal)
            elif temporal.datatype == XSD.gYear or temporal.datatype == XSD.gYearMonth:
                return str(temporal)
        else:
            return str(temporal)
    else:
        has_time = list(graph.objects(temporal, TIME.hasTime))
        if len(has_time) > 0:
            return temporal_to_string(graph, has_time[0])
        has_begins = list(graph.objects(temporal, TIME.hasBeginning))
        has_ends = list(graph.objects(temporal, TIME.hasEnd))
        if len(has_begins) > 0 and len(has_ends) < 1:
            return temporal_to_string(graph, has_begins[0])
        elif len(has_ends) > 0 and len(has_begins) < 1:
            return temporal_to_string(graph, has_ends[0])
        elif len(has_begins) > 0 and len(has_ends) > 0:
            the_begin = temporal_to_string(graph, has_begins[0])
            the_end = temporal_to_string(graph, has_ends[0])
            return f"{the_begin} - {the_end}"
    has_in_xsd_datetimestamp = list(graph.objects(temporal, TIME.inXSDDateTimeStamp))
    if len(has_in_xsd_datetimestamp) > 0:
        return temporal_to_string(graph, has_in_xsd_datetimestamp[0])
    has_in_xsd_datetime = list(graph.objects(temporal, TIME.inXSDDateTime))
    if len(has_in_xsd_datetime) > 0:
        return temporal_to_string(graph, has_in_xsd_datetime[0])
    has_in_xsd_date = list(graph.objects(temporal, TIME.inXSDDate))
    if len(has_in_xsd_date) > 0:
        return temporal_to_string(graph, has_in_xsd_date[0])
    has_in_xsd_gyearmonth = list(graph.objects(temporal, TIME.inXSDgYearMonth))
    if len(has_in_xsd_gyearmonth) > 0:
        return temporal_to_string(graph, has_in_xsd_gyearmonth[0])
    has_in_xsd_gyear = list(graph.objects(temporal, TIME.inXSDgYear))
    if len(has_in_xsd_gyear) > 0:
        return temporal_to_string(graph, has_in_xsd_gyear[0])
    return "Cannot convert time to string: " + str(temporal)


