"""Unit tests for conditional-observable Kappa instrumentation (no KaSim needed)."""

import re

from kappa_smc.instrumentation import instrument_kappa_model


def write_model(tmp_path, content="%agent: A(x)\n%init: 2 A(x)\n"):
    path = tmp_path / "model.ka"
    path.write_text(content, encoding="utf-8")
    return path


def test_single_ap_creates_conditional_observable(tmp_path):
    path = write_model(tmp_path)
    model = instrument_kappa_model(str(path), [("I", "<", "2")])

    assert re.search(r"%obs:\s*'p_0'\s*\(I\s*<\s*2\)\s*\[\?\]\s*1\s*\[:\]\s*0", model)
    assert "%var: 'a_0'" not in model
    assert "$UPDATE 'a_0'" not in model


def test_multiple_aps_have_distinct_observables(tmp_path):
    path = write_model(tmp_path)
    model = instrument_kappa_model(
        str(path), [("I", "<", "2"), ("S", ">", "50"), ("R", "<", "10")]
    )

    for i, condition in enumerate(("I < 2", "S > 50", "R < 10")):
        assert model.count(f"%obs: 'p_{i}'") == 1
        assert f"({condition}) [?] 1 [:] 0" in model
    assert "'p_3'" not in model


def test_original_model_is_preserved_and_comments_removed(tmp_path):
    source = "%agent: A(x) // comment\n%init: 2 A(x)\n"
    path = write_model(tmp_path, source)
    model = instrument_kappa_model(str(path), [("I", "<", "2")])

    assert "%agent: A(x)" in model
    assert "%init: 2 A(x)" in model
    assert "// comment" not in model
    assert path.read_text(encoding="utf-8") == source


def test_empty_ap_list_does_not_add_observables(tmp_path):
    path = write_model(tmp_path)
    model = instrument_kappa_model(str(path), [])

    assert "%agent: A(x)" in model
    assert "'p_0'" not in model
