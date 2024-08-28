# test_01.ttl from https://github.com/opengeospatial/ogc-geosparql/blob/master/examples/demo-dataset.ttl
# validate independently at https://geojson.yanzi.dev/

import pytest
from pathlib import Path
from rdf2geojson import convert
from rdflib import Graph


TEST_DATA_DIR = Path(__file__).parent / "data"


@pytest.mark.parametrize(
    "rdf_file,expected_valid",
    [
        (TEST_DATA_DIR / "test_01.ttl", True),
        (TEST_DATA_DIR / "test_02_invalid.ttl", False),
        (TEST_DATA_DIR / "test_03_empty.ttl", False),
        (TEST_DATA_DIR / "test_04_nogeo.ttl", False),
        (TEST_DATA_DIR / "test_05_onefeature.ttl", True),
        (TEST_DATA_DIR / "test_06_sdo.ttl", True),
    ]
)
def test_convert(rdf_file, expected_valid):
    gj = convert(Graph().parse(rdf_file))

    # save outputs
    from geojson import dump
    with open(rdf_file.with_suffix(".json"), "w") as f2:
        dump(gj, f2, indent=4)

    if expected_valid:
        assert gj.is_valid
    else:
        assert not gj
