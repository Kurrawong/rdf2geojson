# test_01.ttl from https://github.com/opengeospatial/ogc-geosparql/blob/master/examples/demo-dataset.ttl
# validate independently at https://geojson.yanzi.dev/

import pytest
from pathlib import Path
from rdf2geojson import convert
from rdflib import Graph, URIRef

from rdf2geojson.convert import unconvert

TEST_DATA_DIR = Path(__file__).parent / "data"


@pytest.mark.parametrize(
    "rdf_file,expected_valid",
    [
        (TEST_DATA_DIR / "test_01.ttl", True),
        (TEST_DATA_DIR / "test_01b_crs.ttl", True),
        (TEST_DATA_DIR / "test_01c_srid.ttl", True),
        (TEST_DATA_DIR / "test_01d_gda94.ttl", True),
        (TEST_DATA_DIR / "test_02_invalid.ttl", False),
        (TEST_DATA_DIR / "test_03_empty.ttl", False),
        (TEST_DATA_DIR / "test_04_nogeo.ttl", False),
        (TEST_DATA_DIR / "test_05_onefeature.ttl", True),
        (TEST_DATA_DIR / "test_06_sdo.ttl", True),
    ]
)
def test_convert_rdf_to_geojson(rdf_file: Path, expected_valid: bool):
    gj = convert(Graph().parse(rdf_file))

    # save outputs
    from geojson import dump
    with open(rdf_file.with_suffix(".json"), "w") as f2:
        dump(gj, f2, indent=4)

    if expected_valid:
        assert gj.is_valid
    else:
        assert (gj is None or len(gj) == 0 or
                (gj.get("type", None) == "GeoJSON" and gj.get("features", None) is None and
                 gj.get("geometry", None) is None)
                 ), "Expected invalid GeoJSON, but got valid GeoJSON"

@pytest.mark.parametrize(
    "json_file",
    [
        TEST_DATA_DIR / "test_01.json",
        TEST_DATA_DIR / "test_01b_crs.json",
        TEST_DATA_DIR / "test_01c_srid.json",
        TEST_DATA_DIR / "test_01d_gda94.json",
        TEST_DATA_DIR / "test_02_invalid.json",
        TEST_DATA_DIR / "test_03_empty.json",
        TEST_DATA_DIR / "test_04_nogeo.json",
        TEST_DATA_DIR / "test_05_onefeature.json",
    ]
)
def test_convert_geojson_to_wkt(json_file: Path):
    from geojson import load
    with open(json_file) as f:
        gj = load(f)
    if len(gj) < 1:
        pytest.skip("GeoJSON is empty")
    elif len(gj) == 1 and gj.get("type") == "GeoJSON":
         pytest.skip("GeoJSON is not a FeatureCollection or Feature, it's invalid and empty.")
    g: Graph = unconvert(gj)
    with open(json_file.with_suffix(".roundtrip.ttl"), "w") as f:
        f.write(g.serialize(format="turtle"))
        

def test_observation_collection():
        gj = convert(Graph().parse(TEST_DATA_DIR / "test_07_observation_collection.ttl"), do_validate=False, kind="human", fc_uri=URIRef("https://linked.data.gov.au/dataset/bdr/feature-collection/cfeb6552-abf6-4fc2-9fe1-d62f76c8579f"))

def test_big_observation_collection():
        rdf_file = TEST_DATA_DIR / "test_big_observation_collection.ttl"
        gj = convert(Graph().parse(rdf_file), do_validate=False, kind="human", fc_uri=URIRef("https://linked.data.gov.au/dataset/bdr/occurrence-collection/632b575b-b7eb-4804-918e-af7c65a3e4a5"))
        from geojson import dump
        with open(rdf_file.with_suffix(".json"), "w") as f2:
            dump(gj, f2, indent=4)

def test_big_observation_collection_with_attributes():
        rdf_file = TEST_DATA_DIR / "test_big_observation_collection_attributes.ttl"
        gj = convert(Graph().parse(rdf_file), do_validate=False, kind="human", fc_uri=URIRef("https://linked.data.gov.au/dataset/bdr/occurrence-collection/632b575b-b7eb-4804-918e-af7c65a3e4a5"))
        from geojson import dump
        with open(rdf_file.with_suffix(".json"), "w") as f2:
            dump(gj, f2, indent=4)

def test_big_observation_collection_with_attributes_oxigraph():
        rdf_file = TEST_DATA_DIR / "test_big_observation_collection_attributes.ttl"
        g = Graph().parse(rdf_file, format="turtle")
        from pyoxigraph import Store, RdfFormat
        store = Store()
        store.bulk_load(None, format=RdfFormat.TURTLE, path=str(rdf_file))
        gj = convert(store, do_validate=False, kind="human", fc_uri=URIRef("https://linked.data.gov.au/dataset/bdr/occurrence-collection/632b575b-b7eb-4804-918e-af7c65a3e4a5"), namespace_manager=g.namespace_manager)
        from geojson import dump
        with open(rdf_file.with_suffix(".json"), "w") as f2:
            dump(gj, f2, indent=4)

def test_to_taxon_support_human():
        rdf_file = TEST_DATA_DIR / "test_to_taxon_support.ttl"
        g = Graph().parse(rdf_file, format="turtle")
        gj = convert(g, do_validate=False, kind="human", fc_uri=URIRef("https://linked.data.gov.au/dataset/bdr/occurrence-collection/7ca6f4cb-1917-4da0-b65f-912f3d2ffbe8"))
        from geojson import dump
        with open(rdf_file.with_suffix(".human.json"), "w") as f2:
            dump(gj, f2, indent=4)

def test_to_taxon_support_machine():
        rdf_file = TEST_DATA_DIR / "test_to_taxon_support_machine.ttl"
        g = Graph().parse(rdf_file, format="turtle")
        gj = convert(g, do_validate=False, kind="machine", fc_uri=URIRef("https://linked.data.gov.au/dataset/bdr/occurrence-collection/7ca6f4cb-1917-4da0-b65f-912f3d2ffbe8"))
        from geojson import dump
        with open(rdf_file.with_suffix(".json"), "w") as f2:
            dump(gj, f2, indent=4)