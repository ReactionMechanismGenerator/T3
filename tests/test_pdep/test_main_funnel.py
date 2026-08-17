#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_pdep test_main_funnel module

Tests for the QM-assessment funnel in ``T3.determine_species_from_pdep_network()``: the guarantee
that every P-dep network offer leaves exactly one durable record of what happened to it, including
the offers that could not be assessed at all.

That guarantee is the point of the funnel. On a real twelve-network run, seven networks -- the
majority -- left no trace, because the only record was appended deep inside the success path: a
network assessed and found not worth refining and a network never assessed at all produced the same
silence afterwards. Four of the paths that produced that silence were uncaught exceptions that ended
the whole T3 run, so the tests here are mostly about failures that used to be fatal and are now
recorded outcomes -- and about the line between the two, since a bug must never be filed away among
the networks that legitimately could not be evaluated.

Real Arkane is never invoked: ``t3.main.run_arkane_job`` is monkeypatched throughout and a synthetic
sensitivity dictionary stands in for its output, exactly as in ``test_main_wiring.py``. These tests
build a real T3 against the ``network4_2`` fixture tree, so each one costs a full T3 construction;
like ``test_main_wiring.py``, this module is slow enough to be worth excluding from a quick run.
"""

import os

import pytest

import t3.main as t3_main
from arc.common import read_yaml_file, save_yaml_file

from t3.pdep.api import load_pdep_network_assessments
from t3.pdep.cache import write_sa_cache_metadata
from t3.pdep.reason_codes import (REASON_INTERNAL_ERROR,
                                  REASON_NETWORK_DISCOVERY_FAILED,
                                  REASON_NETWORK_INPUT_WRITE_FAILED,
                                  REASON_NETWORK_PARSE_FAILED,
                                  REASON_QUALIFIED_UNCERTAIN_TS,
                                  REASON_SA_ALL_METHODS_FAILED,
                                  REASON_SA_OUTPUT_MALFORMED,
                                  REASON_SA_OUTPUT_MISSING,
                                  REASON_SA_OUTPUT_UNREADABLE,
                                  REASON_SA_STRUCTURES_MISSING,
                                  REASON_SPECIES_LABEL_MAPPING_FAILED,
                                  STATUS_INTERNAL_ERROR,
                                  STATUS_NOT_EVALUATED,
                                  STATUS_QUALIFIED,
                                  )
from t3.pdep.diagram import T3_PES_DIAGRAM_FILENAME
from t3.pdep.selector import CACHE_STATUS_GENERATED

from tests.test_pdep._wiring_helpers import (NETWORK_NAME,
                                             build_pdep_rxns_to_explore as _build_pdep_rxns_to_explore,
                                             build_sa_dict as _build_sa_dict,
                                             build_t3 as _build_t3,
                                             candidate_sa_path as _candidate_sa_path,
                                             network_path as _network_path,
                                             )


def _write_sa_output(t3, method, payload, raw=None):
    """Put a sensitivity artifact where the given ME method's Arkane job would have left one."""
    sa_path = _candidate_sa_path(t3, method=method)
    os.makedirs(os.path.dirname(sa_path), exist_ok=True)
    if raw is not None:
        with open(sa_path, 'w') as f:
            f.write(raw)
    else:
        save_yaml_file(path=sa_path, content=payload)
    return sa_path


def _arkane_writing(t3, payloads):
    """A fake ``run_arkane_job`` that reports success and leaves ``payloads[method]`` behind.

    A method absent from ``payloads`` reports failure; a method mapped to ``None`` reports SUCCESS
    and writes nothing, which is the "Arkane said it worked and produced no SA" case.
    """
    calls = list()

    def _fake_run_arkane_job(input_file, output_directory, plot, logger):
        method = os.path.basename(output_directory)
        calls.append(method)
        if method not in payloads:
            return False
        payload = payloads[method]
        if payload is not None:
            raw = payload if isinstance(payload, str) else None
            _write_sa_output(t3, method, payload if raw is None else None, raw=raw)
        return True

    return _fake_run_arkane_job, calls


def _record_path(t3):
    """Where this iteration's assessment record lives."""
    return t3.paths['PDep network assessments']


def _sole_record(t3):
    """The one assessment this pass should have produced, read back off disk rather than from memory.

    Reading the FILE is the point: an in-memory list that a crash never persisted is exactly the
    observability hole this funnel exists to close.
    """
    records = load_pdep_network_assessments(_record_path(t3))
    assert len(records) == 1, f'expected exactly one assessment, got {[r.network_id for r in records]}'
    return records[0]


class TestEveryOfferLeavesARecord(object):
    """One record per network offer, whatever became of it."""

    def test_a_network_that_qualifies_records_its_verdict_and_the_evidence_for_it(self, tmp_path, monkeypatch):
        t3 = _build_t3(tmp_path)
        pdep_rxns_to_explore = _build_pdep_rxns_to_explore(t3)
        run_arkane_job, _ = _arkane_writing(t3, {'CSE': _build_sa_dict(t3)})
        monkeypatch.setattr(t3_main, 'run_arkane_job', run_arkane_job)

        t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)

        record = _sole_record(t3)
        assert record.network_id == NETWORK_NAME
        assert record.status == STATUS_QUALIFIED
        assert record.reason_code == REASON_QUALIFIED_UNCERTAIN_TS
        assert record.iteration == t3.iteration
        assert record.final_method == 'CSE'
        assert record.cache_status == CACHE_STATUS_GENERATED
        assert record.requested_me_methods == tuple(t3.t3['sensitivity']['ME_methods'])
        # The nested selection is the evidence for the verdict. A record asserting `qualified` with
        # nothing inside it would be a conclusion that can no longer be checked.
        assert record.selection is not None and record.selection.qualified is True
        assert record.selection.network_id == record.network_id

    def test_a_network_no_method_could_analyse_is_recorded_rather_than_left_silent(self, tmp_path, monkeypatch):
        # The case that motivated the whole increment: before it, this network left the run with no
        # durable trace at all, indistinguishable from one that was assessed and found uninteresting.
        t3 = _build_t3(tmp_path)
        pdep_rxns_to_explore = _build_pdep_rxns_to_explore(t3)
        run_arkane_job, calls = _arkane_writing(t3, {})
        monkeypatch.setattr(t3_main, 'run_arkane_job', run_arkane_job)

        t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)

        assert calls == list(t3.t3['sensitivity']['ME_methods']), 'every configured method should get a turn'
        record = _sole_record(t3)
        assert record.status == STATUS_NOT_EVALUATED
        assert record.reason_code == REASON_SA_ALL_METHODS_FAILED
        assert record.selection is None, 'no selector ever ran, so there is no evidence to nest'
        assert record.sa_path is None and record.final_method is None

    def test_the_record_is_written_even_when_there_was_nothing_to_assess(self, tmp_path):
        # Absence of the file has to mean "the assessment never ran this iteration", not "it ran and
        # found no networks" -- otherwise the two are the same on disk.
        t3 = _build_t3(tmp_path)

        t3.determine_species_from_pdep_network(pdep_rxns_to_explore=[])

        assert load_pdep_network_assessments(_record_path(t3)) == []

    def test_the_lists_describe_this_iteration_only(self, tmp_path, monkeypatch):
        # Both verdict lists used to accumulate across iterations, so a network re-examined in three
        # iterations was reported three times, against a model that had changed underneath each one.
        t3 = _build_t3(tmp_path)
        pdep_rxns_to_explore = _build_pdep_rxns_to_explore(t3)
        run_arkane_job, _ = _arkane_writing(t3, {'CSE': _build_sa_dict(t3)})
        monkeypatch.setattr(t3_main, 'run_arkane_job', run_arkane_job)

        t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)
        t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)

        assert len(t3.pdep_network_assessments) == 1
        assert len(t3.pdep_network_selections) == 1
        assert len(load_pdep_network_assessments(_record_path(t3))) == 1

    def test_the_feature_being_off_clears_a_record_an_earlier_run_left_behind(self, tmp_path, monkeypatch):
        t3 = _build_t3(tmp_path)
        pdep_rxns_to_explore = _build_pdep_rxns_to_explore(t3)
        run_arkane_job, _ = _arkane_writing(t3, {'CSE': _build_sa_dict(t3)})
        monkeypatch.setattr(t3_main, 'run_arkane_job', run_arkane_job)
        t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)
        assert os.path.isfile(_record_path(t3))

        t3.t3['sensitivity']['pdep_SA_threshold'] = None
        t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)

        assert not os.path.isfile(_record_path(t3)), \
            'a record from a run where the feature was ON must not survive a run where it is OFF'


class TestFailuresThatUsedToEndTheRun(object):
    """The four crash paths, converted into recorded outcomes."""

    def test_a_network_whose_arkane_input_cannot_be_written_is_recorded(self, tmp_path, monkeypatch):
        t3 = _build_t3(tmp_path)
        pdep_rxns_to_explore = _build_pdep_rxns_to_explore(t3)

        def _refuse_to_write(**kwargs):
            raise ValueError('the network file declares no usable temperature range')

        monkeypatch.setattr(t3_main, 'write_pdep_network_file', _refuse_to_write)
        monkeypatch.setattr(t3_main, 'run_arkane_job',
                            lambda **kwargs: pytest.fail('no method should run after the write failed'))

        t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)

        record = _sole_record(t3)
        assert record.reason_code == REASON_NETWORK_INPUT_WRITE_FAILED
        assert record.status == STATUS_NOT_EVALUATED
        assert any('no usable temperature range' in warning for warning in record.warnings)

    def test_a_filesystem_failure_reading_the_network_file_is_an_outcome_not_a_crash(self, tmp_path,
                                                                                    monkeypatch):
        # The network file is this run's DATA. Whether it is missing, unreadable, or malformed is a
        # fact about the data either way, and which of the two sites that read it -- the input writer
        # or the parser -- happens to touch it first must not decide whether the campaign survives.
        t3 = _build_t3(tmp_path)
        pdep_rxns_to_explore = _build_pdep_rxns_to_explore(t3)

        def _unreadable(**kwargs):
            raise OSError('Permission denied')

        monkeypatch.setattr(t3_main, 'write_pdep_network_file', _unreadable)

        t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)

        record = _sole_record(t3)
        assert record.reason_code == REASON_NETWORK_INPUT_WRITE_FAILED
        assert any('Permission denied' in warning for warning in record.warnings)

    def test_a_programming_error_is_a_bug_recorded_under_its_own_status_and_re_raised(self, tmp_path,
                                                                                     monkeypatch):
        # The line the narrow catches draw. Filesystem trouble is a fact about the run; a TypeError
        # is a fact about the code, and filing it among the networks that could not be evaluated
        # would make a defect look like the ordinary cost of doing business.
        t3 = _build_t3(tmp_path)
        pdep_rxns_to_explore = _build_pdep_rxns_to_explore(t3)

        def _a_real_bug(**kwargs):
            raise TypeError("unsupported operand type(s) for +: 'int' and 'str'")

        monkeypatch.setattr(t3_main, 'write_pdep_network_file', _a_real_bug)

        with pytest.raises(TypeError):
            t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)

        # Recorded on the way out, so the crash is explicable afterwards -- but under its own status,
        # and the file says it is unfinished because it is.
        records = load_pdep_network_assessments(_record_path(t3), allow_incomplete=True)
        assert len(records) == 1
        assert records[0].status == STATUS_INTERNAL_ERROR
        assert records[0].reason_code == REASON_INTERNAL_ERROR
        assert records[0].network_id == NETWORK_NAME
        with pytest.raises(ValueError, match='complete=False'):
            load_pdep_network_assessments(_record_path(t3))

    def test_a_failure_to_record_the_bug_does_not_replace_the_bug(self, tmp_path, monkeypatch):
        # This handler runs precisely when something is already wrong enough that the record may be
        # unbuildable. Surfacing the complaint about writing it down, instead of the diagnosis, would
        # lose the only thing worth having.
        t3 = _build_t3(tmp_path)
        pdep_rxns_to_explore = _build_pdep_rxns_to_explore(t3)

        def _a_real_bug(**kwargs):
            raise TypeError('the original bug')

        monkeypatch.setattr(t3_main, 'write_pdep_network_file', _a_real_bug)
        real_save = t3_main.save_pdep_network_assessments

        def _fail_once_the_pass_has_started(**kwargs):
            # The write that stakes the claim on the file, before any network is tried, must still
            # succeed: if THAT cannot be written the guarantee is void and failing loudly up front is
            # right. What is under test is the write that records the crash.
            if not kwargs['assessments']:
                return real_save(**kwargs)
            raise RuntimeError('and the record failed too')

        monkeypatch.setattr(t3_main, 'save_pdep_network_assessments', _fail_once_the_pass_has_started)

        with pytest.raises(TypeError, match='the original bug'):
            t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)

    @pytest.mark.parametrize('stray_name', ['network4_backup.py',   # unparseable version
                                            'network4_99.txt',      # parseable version, not a network
                                            'network4_2.py.bak',    # the real name with a suffix
                                            ])
    def test_a_stray_file_beside_the_real_ones_is_ignored_not_fatal(self, tmp_path, monkeypatch, stray_name):
        # Each of these used to end a multi-day campaign over a file T3 never needed:
        # 'network4_backup' raised out of int(), and 'network4_99.txt' parsed as version 99 and then
        # synthesized a path to a .py file that does not exist. Neither is a network file, so neither
        # gets a say in which network is the newest.
        t3 = _build_t3(tmp_path)
        pdep_rxns_to_explore = _build_pdep_rxns_to_explore(t3)
        with open(os.path.join(t3.paths['RMG PDep'], stray_name), 'w') as f:
            f.write('# not a network file\n')
        run_arkane_job, _ = _arkane_writing(t3, {'CSE': _build_sa_dict(t3)})
        monkeypatch.setattr(t3_main, 'run_arkane_job', run_arkane_job)

        t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)

        record = _sole_record(t3)
        assert record.network_id == NETWORK_NAME, 'the real network must still have been assessed'
        assert record.status == STATUS_QUALIFIED

    def test_a_network_with_no_versioned_file_at_all_is_recorded(self, tmp_path, monkeypatch):
        # Files share the prefix and none of them names a version, so there is nothing to assess and
        # no name to assess it under -- the one case that genuinely cannot proceed.
        t3 = _build_t3(tmp_path)
        pdep_rxns_to_explore = _build_pdep_rxns_to_explore(t3)
        for name in os.listdir(t3.paths['RMG PDep']):
            if name.startswith('network4_'):
                os.rename(os.path.join(t3.paths['RMG PDep'], name),
                          os.path.join(t3.paths['RMG PDep'], 'network4_orig.py'))
        monkeypatch.setattr(t3_main, 'run_arkane_job',
                            lambda **kwargs: pytest.fail('the network was never identified'))

        t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)

        record = _sole_record(t3)
        assert record.reason_code == REASON_NETWORK_DISCOVERY_FAILED
        assert record.status == STATUS_NOT_EVALUATED
        # Identified by index alone: which revision this was about is exactly what could not be
        # determined, so a stem would be inventing the one missing fact.
        assert record.network_id == 'network4'
        assert record.network_path is None

    def test_a_network_file_that_cannot_be_parsed_is_recorded_but_still_yields_wells(self, tmp_path, monkeypatch):
        # The well analysis needs the SA output and the label map, not the parsed network. Refusing it
        # here would drop species this iteration would otherwise refine, for an unrelated reason.
        t3 = _build_t3(tmp_path)
        pdep_rxns_to_explore = _build_pdep_rxns_to_explore(t3)
        run_arkane_job, _ = _arkane_writing(t3, {'CSE': _build_sa_dict(t3)})
        monkeypatch.setattr(t3_main, 'run_arkane_job', run_arkane_job)
        wells = list()
        monkeypatch.setattr(t3_main, 'parse_pdep_network_file',
                            lambda path: (_ for _ in ()).throw(ValueError('unbalanced parentheses')))
        real_select_sensitive_wells = t3_main.select_sensitive_wells

        def _spy(**kwargs):
            result = real_select_sensitive_wells(**kwargs)
            wells.append(result)
            return result

        monkeypatch.setattr(t3_main, 'select_sensitive_wells', _spy)

        t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)

        record = _sole_record(t3)
        assert record.reason_code == REASON_NETWORK_PARSE_FAILED
        assert record.selection is None, 'the selector never ran, so there is no evidence to nest'
        assert record.network_reaction is not None, 'the labels DID map; only the network file failed'
        assert len(wells) == 1, 'the well analysis must still have run'

    def test_a_species_label_mismatch_is_recorded_and_stops_the_well_analysis(self, tmp_path, monkeypatch):
        # Deliberately stricter than the parse failure above. Without a grounded network reaction
        # string there is nothing to look up in the SA output, and falling back to an ungrounded one
        # would risk refining species for a reaction nobody asked about.
        t3 = _build_t3(tmp_path)
        pdep_rxns_to_explore = _build_pdep_rxns_to_explore(t3)
        sa_dict = _build_sa_dict(t3)
        # A structures block that is valid, and about other species entirely.
        sa_dict['structures'] = {'Ar': '1 Ar u0 p4 c0\n'}
        run_arkane_job, _ = _arkane_writing(t3, {'CSE': sa_dict})
        monkeypatch.setattr(t3_main, 'run_arkane_job', run_arkane_job)
        monkeypatch.setattr(t3_main, 'select_sensitive_wells',
                            lambda **kwargs: pytest.fail('wells must not be drawn from an ungrounded reaction'))

        species_keys = t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)

        record = _sole_record(t3)
        assert record.reason_code == REASON_SPECIES_LABEL_MAPPING_FAILED
        assert record.status == STATUS_NOT_EVALUATED
        assert record.chemkin_reaction is not None, 'the Chemkin side of the mapping was known'
        assert record.network_reaction is None, 'the network side is exactly what could not be built'
        assert species_keys == []


    @pytest.mark.skipif(hasattr(os, 'geteuid') and os.geteuid() == 0,
                        reason='root bypasses file permissions, so the unreadable file would be readable')
    def test_an_unreadable_sa_artifact_is_an_outcome_not_a_crash(self, tmp_path):
        """Test that an OSError reading the SA output is a reason code, not the end of the campaign.

        Read directly rather than through a full pass: the point is the classification, and the
        condition is one the cache check immediately before it has just proved was not true, so
        provoking it through the funnel would mean racing the two reads against each other.

        This is the sibling of the fail-closed hashing in ``t3.pdep.cache.validate_sa_cache``. The
        two run back to back on the cache-hit path, so catching the I/O error in one and letting it
        escape the other would just move the campaign-ending crash two lines later.
        """
        sa_path = str(tmp_path / 'sa_coefficients.yml')
        save_yaml_file(path=sa_path, content={'structures': {}})
        os.chmod(sa_path, 0o000)
        try:
            sa_dict, reason_code, warning = t3_main._read_pdep_sa_output(sa_path)
        finally:
            os.chmod(sa_path, 0o644)
        assert sa_dict is None
        assert reason_code == REASON_SA_OUTPUT_UNREADABLE
        assert 'could not be read' in warning


class TestTheSAArtifactIsJudgedPerMethod(object):
    """A bad artifact from one ME method must not spend the other methods' turns."""

    def test_a_method_that_produces_no_sa_does_not_deny_the_next_method_its_turn(self, tmp_path, monkeypatch):
        # The subtle one. Arkane reporting success and leaving nothing readable behind is a fact about
        # THAT method's output; until the read moved inside the loop, the first such method ended the
        # network and the later ones never ran.
        t3 = _build_t3(tmp_path)
        methods = list(t3.t3['sensitivity']['ME_methods'])
        assert len(methods) > 1, 'this test needs a fallback method to fall back to'
        pdep_rxns_to_explore = _build_pdep_rxns_to_explore(t3)
        run_arkane_job, calls = _arkane_writing(t3, {methods[0]: None, methods[1]: _build_sa_dict(t3)})
        monkeypatch.setattr(t3_main, 'run_arkane_job', run_arkane_job)

        t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)

        assert calls[:2] == methods[:2]
        record = _sole_record(t3)
        assert record.final_method == methods[1]
        assert record.status == STATUS_QUALIFIED
        # The first method's failure is still on the record; it just did not get to be the verdict.
        assert any('does not exist' in warning for warning in record.warnings)

    def test_no_method_producing_a_readable_sa_is_recorded_as_a_missing_output(self, tmp_path, monkeypatch):
        t3 = _build_t3(tmp_path)
        methods = list(t3.t3['sensitivity']['ME_methods'])
        pdep_rxns_to_explore = _build_pdep_rxns_to_explore(t3)
        # The fixture tree stages an MSC artifact, so "no method produced one" has to be arranged
        # rather than assumed: left in place, MSC would be read successfully and assessed.
        for method in methods:
            staged = _candidate_sa_path(t3, method=method)
            if os.path.isfile(staged):
                os.remove(staged)
        run_arkane_job, calls = _arkane_writing(t3, {method: None for method in methods})
        monkeypatch.setattr(t3_main, 'run_arkane_job', run_arkane_job)

        t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)

        assert calls == methods
        record = _sole_record(t3)
        # The specific diagnosis wins over the coarse one: "no method produced usable SA" is also
        # true, but it is the same failure at a coarser grain, not a second independent one.
        assert record.reason_code == REASON_SA_OUTPUT_MISSING
        assert record.secondary_reason_codes == ()

    def test_an_sa_output_that_is_not_valid_yaml_is_recorded_as_unreadable(self, tmp_path, monkeypatch):
        t3 = _build_t3(tmp_path)
        methods = list(t3.t3['sensitivity']['ME_methods'])
        pdep_rxns_to_explore = _build_pdep_rxns_to_explore(t3)
        run_arkane_job, _ = _arkane_writing(t3, {method: '{unclosed: [mapping\n' for method in methods})
        monkeypatch.setattr(t3_main, 'run_arkane_job', run_arkane_job)

        t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)

        record = _sole_record(t3)
        assert record.reason_code == REASON_SA_OUTPUT_UNREADABLE

    def test_an_sa_payload_that_is_readable_but_not_a_mapping_is_recorded_as_malformed(self, tmp_path, monkeypatch):
        # Readable, well-formed YAML that is simply not what an SA output is. This is caught before
        # the selector can ever see it, which is why it has its own code and not the selector's.
        t3 = _build_t3(tmp_path)
        methods = list(t3.t3['sensitivity']['ME_methods'])
        pdep_rxns_to_explore = _build_pdep_rxns_to_explore(t3)
        run_arkane_job, _ = _arkane_writing(t3, {method: ['not', 'a', 'mapping'] for method in methods})
        monkeypatch.setattr(t3_main, 'run_arkane_job', run_arkane_job)

        t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)

        record = _sole_record(t3)
        assert record.reason_code == REASON_SA_OUTPUT_MALFORMED
        assert record.sa_path is not None, 'the artifact WAS read; it just says nothing usable'

    def test_an_sa_output_with_no_structures_block_is_recorded(self, tmp_path, monkeypatch):
        t3 = _build_t3(tmp_path)
        methods = list(t3.t3['sensitivity']['ME_methods'])
        pdep_rxns_to_explore = _build_pdep_rxns_to_explore(t3)
        sa_dict = _build_sa_dict(t3)
        sa_dict.pop('structures')
        run_arkane_job, _ = _arkane_writing(t3, {method: sa_dict for method in methods})
        monkeypatch.setattr(t3_main, 'run_arkane_job', run_arkane_job)

        t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)

        record = _sole_record(t3)
        assert record.reason_code == REASON_SA_STRUCTURES_MISSING
        assert record.chemkin_reaction is not None


    def test_methods_that_broke_in_different_ways_all_reach_the_record(self, tmp_path, monkeypatch):
        # Which method happened to run last is not a diagnosis. Two different broken artifacts are
        # two independent defects, and a record naming only one of them undercounts the other.
        t3 = _build_t3(tmp_path)
        methods = list(t3.t3['sensitivity']['ME_methods'])
        assert len(methods) > 1, 'this test needs two methods to break differently'
        pdep_rxns_to_explore = _build_pdep_rxns_to_explore(t3)
        for method in methods:
            staged = _candidate_sa_path(t3, method=method)
            if os.path.isfile(staged):
                os.remove(staged)
        run_arkane_job, _ = _arkane_writing(t3, {methods[0]: None,                    # writes nothing
                                                 methods[1]: '{unclosed: [mapping\n'})  # writes garbage
        monkeypatch.setattr(t3_main, 'run_arkane_job', run_arkane_job)

        t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)

        record = _sole_record(t3)
        assert record.reason_code == REASON_SA_OUTPUT_MISSING
        assert record.secondary_reason_codes == (REASON_SA_OUTPUT_UNREADABLE,)

    def test_a_method_that_failed_to_run_is_named_on_the_record(self, tmp_path, monkeypatch):
        # Otherwise the codes make it look as though that method had never been configured at all.
        t3 = _build_t3(tmp_path)
        methods = list(t3.t3['sensitivity']['ME_methods'])
        pdep_rxns_to_explore = _build_pdep_rxns_to_explore(t3)
        for method in methods:
            staged = _candidate_sa_path(t3, method=method)
            if os.path.isfile(staged):
                os.remove(staged)
        run_arkane_job, _ = _arkane_writing(t3, {methods[0]: None})  # methods[1] reports failure
        monkeypatch.setattr(t3_main, 'run_arkane_job', run_arkane_job)

        t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)

        record = _sole_record(t3)
        assert record.reason_code == REASON_SA_OUTPUT_MISSING
        assert any(methods[1] in warning and 'did not complete successfully' in warning
                   for warning in record.warnings)


class TestTheRecordOnDisk(object):
    """What the artifact itself claims."""

    def test_a_finished_pass_marks_the_record_complete(self, tmp_path, monkeypatch):
        t3 = _build_t3(tmp_path)
        pdep_rxns_to_explore = _build_pdep_rxns_to_explore(t3)
        run_arkane_job, _ = _arkane_writing(t3, {'CSE': _build_sa_dict(t3)})
        monkeypatch.setattr(t3_main, 'run_arkane_job', run_arkane_job)

        t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)

        assert read_yaml_file(_record_path(t3))['complete'] is True

    def test_the_record_is_written_as_the_pass_goes_not_only_at_the_end(self, tmp_path, monkeypatch):
        # If it were written only on the way out it would be absent in exactly the case it exists to
        # explain, since the way out is what the crash prevented.
        t3 = _build_t3(tmp_path)
        pdep_rxns_to_explore = _build_pdep_rxns_to_explore(t3)

        def _a_real_bug(**kwargs):
            raise TypeError('the machine gave up')

        monkeypatch.setattr(t3_main, 'write_pdep_network_file', _a_real_bug)

        with pytest.raises(TypeError):
            t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)

        assert os.path.isfile(_record_path(t3)), 'the crash still left its own explanation behind'
        assert read_yaml_file(_record_path(t3))['complete'] is False
        assert len(load_pdep_network_assessments(_record_path(t3), allow_incomplete=True)) == 1

    def test_an_unparseable_network_with_a_valid_cache_still_reaches_the_sa(self, tmp_path, monkeypatch):
        """A network file the parser refuses must not cost a valid cached SA its well analysis.

        The cache-hit path reads the network file early, for isomer labels that feed only
        ``executed_networks``. That read must stay opportunistic. The authoritative parse happens
        much later and is deliberately non-fatal -- it records ``network_parse_failed`` and lets the
        well analysis run anyway, because that analysis needs the SA and the label map, not the
        parsed network. An early return here would silently drop the species this iteration would
        otherwise have refined, and would do it for a bookkeeping list nothing reads.

        What is asserted is narrower than that motivation: this pins that the cached SA is reached
        and read, which is the precondition the well analysis needs. Whether the analysis then
        yields anything depends on the caller's ``direction_key``, which is out of this test's scope.
        """
        t3 = _build_t3(tmp_path)
        pdep_rxns_to_explore = _build_pdep_rxns_to_explore(t3)
        sa_path = _candidate_sa_path(t3, method='CSE')
        os.makedirs(os.path.dirname(sa_path), exist_ok=True)
        save_yaml_file(sa_path, _build_sa_dict(t3))
        # Break the network file FIRST, then vouch for it: the sidecar records the hash of whatever
        # is on disk, so this is a genuinely valid cache pointing at a file the parser cannot read.
        network_path = _network_path(t3)
        with open(network_path, 'w') as f:
            f.write('network(  # unbalanced, so ast.parse refuses it\n')
        write_sa_cache_metadata(sa_path=sa_path, network_path=network_path,
                                network_id=NETWORK_NAME, method='CSE')

        def _fail_if_called(*args, **kwargs):
            pytest.fail('run_arkane_job should not be invoked when a valid cache is present.')

        monkeypatch.setattr(t3_main, 'run_arkane_job', _fail_if_called)

        t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)

        record = _sole_record(t3)
        assert record.reason_code == REASON_NETWORK_PARSE_FAILED
        # The discriminator between the two ways this reason code can be reached: these are only set
        # once the SA has actually been read, so an early return would leave them None.
        assert record.sa_path is not None, 'returned before the cached SA was ever read'
        assert record.final_method == 'CSE'

    def test_a_finished_record_from_an_earlier_run_cannot_survive_this_one(self, tmp_path, monkeypatch):
        # The window Codex found in round 58: a `complete: true` file from a previous run of the SAME
        # iteration would otherwise survive a crash that happens before the first network is
        # recorded, and answer "which networks were never evaluated?" with another run's networks.
        t3 = _build_t3(tmp_path)
        pdep_rxns_to_explore = _build_pdep_rxns_to_explore(t3)
        run_arkane_job, _ = _arkane_writing(t3, {'CSE': _build_sa_dict(t3)})
        monkeypatch.setattr(t3_main, 'run_arkane_job', run_arkane_job)
        t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)
        assert read_yaml_file(_record_path(t3))['complete'] is True

        # Injected at the cache check rather than at the writer: the first pass above left a valid
        # sidecar, so the second pass reuses it and never renders an Arkane input at all. The bug
        # this test needs is "something crashes before the first network is recorded", and the cache
        # check is the one call every method makes on both the hit and the miss path.
        monkeypatch.setattr(t3_main, 'validate_sa_cache',
                            lambda **kwargs: (_ for _ in ()).throw(TypeError('a second run, and a bug')))
        with pytest.raises(TypeError):
            t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)

        assert read_yaml_file(_record_path(t3))['complete'] is False, \
            'the stale finished record must not have survived into this pass'
        records = load_pdep_network_assessments(_record_path(t3), allow_incomplete=True)
        assert [record.reason_code for record in records] == [REASON_INTERNAL_ERROR]


class TestPESDiagram(object):
    """The normalized-E0 PES diagram is drawn for every successfully-parsed network, but a
    failure to draw it must never take the assessment down with it (see t3.main.py's
    ``_draw_pdep_pes_diagram``).
    """

    def _diagram_path(self, t3):
        return os.path.join(t3.paths['PDep SA'], NETWORK_NAME, T3_PES_DIAGRAM_FILENAME)

    def test_a_successful_assessment_draws_the_pes_diagram(self, tmp_path, monkeypatch):
        t3 = _build_t3(tmp_path)
        pdep_rxns_to_explore = _build_pdep_rxns_to_explore(t3)
        run_arkane_job, _ = _arkane_writing(t3, {'CSE': _build_sa_dict(t3)})
        monkeypatch.setattr(t3_main, 'run_arkane_job', run_arkane_job)

        t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)

        diagram_path = self._diagram_path(t3)
        assert os.path.isfile(diagram_path)
        assert os.path.getsize(diagram_path) > 0
        with open(diagram_path, 'rb') as f:
            assert f.read(8) == b'\x89PNG\r\n\x1a\n'
        # Drawing the diagram is not allowed to have changed the assessment itself.
        record = _sole_record(t3)
        assert record.status == STATUS_QUALIFIED

    def test_a_diagram_failure_is_a_logged_warning_not_a_crash(self, tmp_path, monkeypatch):
        t3 = _build_t3(tmp_path)
        pdep_rxns_to_explore = _build_pdep_rxns_to_explore(t3)
        run_arkane_job, _ = _arkane_writing(t3, {'CSE': _build_sa_dict(t3)})
        monkeypatch.setattr(t3_main, 'run_arkane_job', run_arkane_job)

        logged_warnings = []
        monkeypatch.setattr(t3.logger, 'warning', lambda message: logged_warnings.append(message))

        def _raise(*args, **kwargs):
            raise ValueError('a species with no E0, e.g.')

        monkeypatch.setattr(t3_main, 'compute_pes_diagram_data', _raise)

        t3.determine_species_from_pdep_network(pdep_rxns_to_explore=pdep_rxns_to_explore)

        assert not os.path.isfile(self._diagram_path(t3))
        assert any('PES diagram' in message and NETWORK_NAME in message for message in logged_warnings), \
            'the diagram failure should have been logged as a warning naming the network'
        # The assessment outcome itself must be unaffected by the diagram failure.
        record = _sole_record(t3)
        assert record.status == STATUS_QUALIFIED
        assert record.reason_code == REASON_QUALIFIED_UNCERTAIN_TS
