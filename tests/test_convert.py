# test_01.ttl from https://github.com/opengeospatial/ogc-geosparql/blob/master/examples/demo-dataset.ttl
# validate independently at https://geojson.yanzi.dev/

import pytest
from pathlib import Path
from rdf2geojson import convert
from rdflib import Graph

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
        assert not gj

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
    g: Graph = unconvert(gj)
    with open(json_file.with_suffix(".roundtrip.ttl"), "w") as f:
        f.write(g.serialize(format="turtle"))