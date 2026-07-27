#!/usr/bin/env python3
# encoding: utf-8

"""
t3 tests test_pdep test_me_success

Fixtures live under ``tests/data/pdep_me/`` rather than ``tests/data/pdep_network/`` because the
latter is a ``t3.paths['PDep SA']`` target that other tests ``shutil.rmtree`` during teardown
(see e.g. ``tests/test_main.py`` lines that ``rmtree(t3.paths['PDep SA'], ...)``); a fixture
placed there would be silently deleted by an unrelated test, making this suite pass or fail
purely on test-execution ordering. ``tests/data/pdep_me/`` is not referenced by any
``t3.paths[...]`` entry or ``rmtree`` call in the test suite.
"""

import math
import os

from t3.common import TEST_DATA_BASE_PATH
from t3.pdep.me_success import MESuccessResult, check_arkane_me_success
from t3.pdep.parser import parse_pdep_network_file

PDEP_ME_DIR = os.path.join(TEST_DATA_BASE_PATH, 'pdep_me')
SUCCESS_OUTPUT = os.path.join(PDEP_ME_DIR, 'success', 'output.py')
SOFT_FAILURE_CSE_OUTPUT = os.path.join(PDEP_ME_DIR, 'soft_failure_cse', 'output.py')
HARD_FAILURE_OUTPUT = os.path.join(PDEP_ME_DIR, 'hard_failure', 'output.py')


def test_success_fixture():
    """Test that the real clean-ME-solve fixture is recognized as a success.

    Regression this guards: any change that makes ``check_arkane_me_success`` reject a
    perfectly good ME solve (e.g., an overly strict numeric check, or accidentally requiring
    positive Chebyshev coefficients) would flip ``succeeded`` to False here.
    """
    result = check_arkane_me_success(output_path=SUCCESS_OUTPUT)
    assert isinstance(result, MESuccessResult)
    assert result.succeeded is True
    assert result.reasons == tuple()
    assert len(result.reactions) == 1
    reaction = result.reactions[0]
    assert reaction.kinetics_type == 'Chebyshev'
    assert len(reaction.numeric_values) > 0
    assert all(v is not None for v in reaction.numeric_values)
    assert all(math.isfinite(v) for v in reaction.numeric_values)


def test_soft_failure_cse_fixture_rejected_with_explicit_reason():
    """Test the single most important case: a CSE-method ME solve that Arkane itself reports
    as successful (exit 0, empty stderr) but whose Chebyshev coefficients are all None must be
    rejected, and the reason must explicitly name the non-finite/None kinetics payload.

    Regression this guards: this is the entire reason this module exists. Any implementation
    that only checks exit code, stderr, or file existence/non-emptiness (all of which look fine
    for this fixture) would wrongly report ``succeeded is True`` here. A version of check 4 that
    is missing, skipped, or that only checks `bool(value)` instead of `math.isfinite` would also
    let this fixture through, since ``None`` values are falsy-but-not-caught by a naive truthy
    check in some codepaths (e.g. `if not value`) while still being non-finite here.
    """
    result = check_arkane_me_success(output_path=SOFT_FAILURE_CSE_OUTPUT)
    assert result.succeeded is False
    assert len(result.reactions) == 1
    assert result.reactions[0].kinetics_type == 'Chebyshev'
    assert any(v is None for v in result.reactions[0].numeric_values)
    joined_reasons = ' '.join(result.reasons)
    assert 'non-finite' in joined_reasons
    assert 'Chebyshev' in joined_reasons


def test_hard_failure_fixture_rejected_as_empty():
    """Test that the 0-byte output.py (Arkane pre-created it, then crashed in job.execute())
    is rejected, with a reason naming the emptiness rather than merely "does not exist".

    Regression this guards: a naive implementation that only checks ``os.path.isfile`` (as ARC's
    ``run_arkane`` does) would treat this fixture as a success, since the pre-created file does
    exist. This test would fail if the emptiness check (check 2) were removed or only logged the
    file as "not found".
    """
    assert os.path.isfile(HARD_FAILURE_OUTPUT)
    assert os.path.getsize(HARD_FAILURE_OUTPUT) == 0
    result = check_arkane_me_success(output_path=HARD_FAILURE_OUTPUT)
    assert result.succeeded is False
    assert result.reactions == tuple()
    joined_reasons = ' '.join(result.reasons)
    assert 'empty' in joined_reasons or '0 byte' in joined_reasons


def test_missing_file_rejected_without_raising():
    """Test that a nonexistent output.py path is reported as a failure, not an exception.

    Regression this guards: an implementation that lets ``open()`` raise
    ``FileNotFoundError`` straight out of ``check_arkane_me_success`` (e.g. by parsing before
    checking existence) would crash the caller instead of returning a normal result.
    """
    missing_path = os.path.join(PDEP_ME_DIR, 'does_not_exist', 'output.py')
    result = check_arkane_me_success(output_path=missing_path)
    assert result.succeeded is False
    assert result.reactions == tuple()
    assert any('does not exist' in reason for reason in result.reasons)


def test_exit_code_and_stderr_are_folded_into_reasons_on_failure():
    """Test that a caller-supplied non-zero exit code and stderr text are folded into the
    reasons alongside the payload-derived reason, on the hard-failure fixture.

    Regression this guards: if the exit_code/stderr parameters were accepted but silently
    ignored (dead parameters), this test would fail because none of the caller-supplied text
    would appear in ``reasons``.
    """
    result = check_arkane_me_success(
        output_path=HARD_FAILURE_OUTPUT,
        exit_code=1,
        stderr='rmgpy.exceptions.NetworkError: boom',
    )
    assert result.succeeded is False
    joined_reasons = ' '.join(result.reasons)
    assert 'non-zero' in joined_reasons
    assert 'boom' in joined_reasons


def test_ignorable_stderr_phrases_do_not_block_success():
    """Test that benign stderr noise (mirroring ARC's ignorable_phrases) does not, by itself,
    cause a real success fixture to be rejected.

    Regression this guards: an implementation that folds *any* non-empty stderr into a failure
    reason (without filtering benign lines) would still leave ``reasons`` non-empty here even
    though the payload is fine; since the payload check must stand on its own for ``succeeded``,
    this test also guards against a version that lets stray stderr override a clean payload.
    """
    result = check_arkane_me_success(
        output_path=SUCCESS_OUTPUT,
        exit_code=0,
        stderr='Open Babel Warning  in some_function\n',
    )
    assert result.succeeded is True
    assert result.reasons == tuple()


def test_expected_reactions_count_mismatch_is_reported():
    """Test that the optional expected_reactions count check rejects a fixture that otherwise
    looks fine when the caller expects more reactions than are present.

    Regression this guards: if the optional expected_reactions argument were unimplemented or a
    no-op, this test would fail because ``succeeded`` would remain True for a mismatched count.
    """
    result = check_arkane_me_success(output_path=SUCCESS_OUTPUT, expected_reactions=2)
    assert result.succeeded is False
    assert any('Expected 2' in reason for reason in result.reasons)


def test_expected_reactions_count_match_still_succeeds():
    """Test that the optional expected_reactions count check does not itself break a matching,
    otherwise-successful fixture.

    Regression this guards: an off-by-one or inverted comparison in the count check (e.g.
    ``>=`` instead of ``==``, or checking the wrong side) would flip this to a false failure.
    """
    result = check_arkane_me_success(output_path=SUCCESS_OUTPUT, expected_reactions=1)
    assert result.succeeded is True


def test_nonzero_exit_code_fails_even_with_a_well_formed_payload():
    """Test that a non-zero exit status fails the job even when the payload itself parses clean.

    Regression this guards: the process-level reasons are collected before the payload checks
    run, so an implementation that returns ``succeeded=True, reasons=()`` on the strength of the
    payload alone would silently discard the only record that Arkane died. Arkane can write a
    complete result for one network and then crash later in the same run, so a well-formed
    payload does not vouch for the run as a whole.
    """
    result = check_arkane_me_success(output_path=SUCCESS_OUTPUT, exit_code=1)
    assert result.succeeded is False
    assert result.reasons != tuple()
    assert any('exit' in reason.lower() for reason in result.reasons)
    # The payload is still parsed and returned, so the caller can see what was salvageable.
    assert len(result.reactions) == 1


def test_real_stderr_fails_even_with_a_well_formed_payload():
    """Test that non-ignorable stderr text fails the job even when the payload parses clean.

    Regression this guards: the same discard-on-success bug as above, reached via stderr rather
    than the exit status.
    """
    result = check_arkane_me_success(output_path=SUCCESS_OUTPUT,
                                     stderr='rmgpy.exceptions.NetworkError: boom',
                                     )
    assert result.succeeded is False
    assert any('boom' in reason for reason in result.reasons)


def test_empty_reactants_and_products_are_rejected(tmp_path):
    """Test that a pdepreaction with empty reactants and products is rejected even when its
    kinetics numbers are all finite.

    Regression this guards: ``pdepreaction(reactants=[], products=[], kinetics=Arrhenius(...))``
    has finite numeric leaves and used to be accepted, but it is not a valid k(T,P) result --
    a rate without both a source and a destination is meaningless. The failure reason must name
    the offending entry's index.
    """
    text = ("pdepreaction(\n"
            "    reactants=[], products=[],\n"
            "    kinetics=Arrhenius(A=(1,'s^-1'), n=0, Ea=(0,'kJ/mol')),\n"
            ")\n")
    output_path = str(tmp_path / 'output.py')
    with open(output_path, 'w') as f:
        f.write(text)
    result = check_arkane_me_success(output_path=output_path)
    assert result.succeeded is False
    joined_reasons = ' '.join(result.reasons)
    assert 'Reaction #0' in joined_reasons
    assert 'reactants' in joined_reasons or 'products' in joined_reasons


def test_finite_bounds_alone_do_not_satisfy_the_payload_check(tmp_path):
    """Test that finite T/P bounds can never satisfy the numeric-payload check on their own.

    Regression this guards: Tmin/Tmax/Pmin/Pmax are themselves numeric leaves, so a
    ``Chebyshev(coeffs=[], ...)`` call that also carries the (always finite) bounds used to
    yield a non-empty, all-finite ``numeric_values`` and be accepted as a solved network with
    NO rate coefficients at all. The rate payload (``coeffs``) must be checked separately from
    the bounds/metadata.
    """
    text = ("pdepreaction(\n"
            "    reactants=['A'], products=['B'],\n"
            "    kinetics=Chebyshev(coeffs=[], kunits='s^-1',\n"
            "                       Tmin=(300,'K'), Tmax=(2100,'K'),\n"
            "                       Pmin=(0.1,'bar'), Pmax=(100,'bar')),\n"
            ")\n")
    output_path = str(tmp_path / 'output.py')
    with open(output_path, 'w') as f:
        f.write(text)
    result = check_arkane_me_success(output_path=output_path)
    assert result.succeeded is False
    joined_reasons = ' '.join(result.reasons).lower()
    assert 'rate coefficient' in joined_reasons


def test_unparseable_nan_coeffs_with_finite_bounds_are_rejected(tmp_path):
    """Test that coeffs which cannot be literal-evaluated (a bare ``nan`` name) are rejected
    even when finite T/P bounds are present.

    Regression this guards: the parser omits any kinetics keyword whose value cannot be
    evaluated as a literal, so ``coeffs=[[nan]]`` used to simply vanish, leaving only the
    finite bounds behind -- and the bounds alone satisfied the old combined numeric check.
    An omitted ``coeffs`` must be distinguishable from a legitimately absent one, and must
    fail the gate.
    """
    text = ("pdepreaction(\n"
            "    reactants=['A'], products=['B'],\n"
            "    kinetics=Chebyshev(coeffs=[[nan]], kunits='s^-1',\n"
            "                       Tmin=(300,'K'), Tmax=(2100,'K'),\n"
            "                       Pmin=(0.1,'bar'), Pmax=(100,'bar')),\n"
            ")\n")
    output_path = str(tmp_path / 'output.py')
    with open(output_path, 'w') as f:
        f.write(text)
    result = check_arkane_me_success(output_path=output_path)
    assert result.succeeded is False
    assert result.reasons


def test_kinetics_with_no_numeric_parameters_is_rejected(tmp_path):
    """Test that a pdepreaction whose kinetics carries no numeric leaves at all is rejected.

    Regression this guards: the all-finite check iterates the numeric values and fails on any
    non-finite one, so an empty value list passes it VACUOUSLY -- ``all finite`` is trivially
    true of nothing. A kinetics form T3 cannot parse, or an empty ``coeffs=[]``, would then be
    read as a good solve. The guard must reject "no rate information" as firmly as "bad rate
    information".
    """
    text = ("pdepreaction(\n"
            "    reactants = ['A'],\n"
            "    products = ['B'],\n"
            "    kinetics = Chebyshev(coeffs=[], kunits='s^-1'),\n"
            ")\n")
    path = str(tmp_path / 'no_numeric_params_output.py')
    with open(path, 'w') as f:
        f.write(text)
    result = check_arkane_me_success(output_path=path)
    assert result.succeeded is False
    assert any('no numeric' in reason.lower() for reason in result.reasons)


def test_expected_reactions_derived_from_the_network_accepts_a_complete_solve():
    """Test that the count derived from the network file accepts the solve Arkane produced from it.

    ``tests/data/pdep_me/success_multi/`` holds a real Arkane MSC solve of ``network4_1.py``. The
    expected count is not restated here; it is derived from the network file itself, so this test
    fails if the derivation and the real Arkane enumeration ever drift apart.
    """
    network = parse_pdep_network_file(path=os.path.join(TEST_DATA_BASE_PATH, 'pdep_network', 'iteration_1',
                                                        'RMG', 'pdep', 'network4_1.py'))
    result = check_arkane_me_success(output_path=os.path.join(PDEP_ME_DIR, 'success_multi', 'output.py'),
                                     exit_code=0,
                                     stderr='',
                                     expected_reactions=network.expected_net_reaction_count())
    assert result.succeeded is True, result.reasons


def test_a_truncated_output_is_rejected_by_the_expected_reaction_count(tmp_path):
    """Test that a solve cut short partway through writing its net reactions is rejected.

    Every entry a truncated ``output.py`` does contain is syntactically perfect and numerically
    finite, so no per-entry check can catch this -- only the count can. This is the hole that
    ``expected_reactions`` exists to close, and it stayed open while nothing passed the argument.
    """
    network = parse_pdep_network_file(path=os.path.join(TEST_DATA_BASE_PATH, 'pdep_network', 'iteration_1',
                                                        'RMG', 'pdep', 'network4_1.py'))
    with open(os.path.join(PDEP_ME_DIR, 'success_multi', 'output.py'), 'r') as f:
        text = f.read()
    # Cut the file at the start of the last complete entry, simulating a solve that died partway.
    truncated_text = text[:text.rindex('pdepreaction(')]
    truncated_path = str(tmp_path / 'truncated_output.py')
    with open(truncated_path, 'w') as f:
        f.write(truncated_text)
    result = check_arkane_me_success(output_path=truncated_path,
                                     exit_code=0,
                                     stderr='',
                                     expected_reactions=network.expected_net_reaction_count())
    assert result.succeeded is False
    assert any('expected' in reason.lower() for reason in result.reasons), result.reasons
    # Without the count, the very same truncated file reads as a clean solve -- that is the hole.
    assert check_arkane_me_success(output_path=truncated_path, exit_code=0, stderr='').succeeded is True
