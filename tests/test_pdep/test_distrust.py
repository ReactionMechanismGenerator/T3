"""Tests for t3.pdep.distrust -- the pre-QM rank-by-distrust selection criterion (I-032)."""

from t3.pdep.parser import (PDepNetworkE0, PDepPathReaction, parse_pdep_network_e0_text,
                            parse_pdep_network_text)
from t3.pdep.pes_rounds import (CandidateSplit, QMCandidate, SKIP_NO_EVIDENCE, SKIP_OUTSIDE_WINDOW,
                                SKIP_TRUSTED_PROVENANCE, SKIP_UNMEASURABLE)
from t3.pdep.distrust import (DistrustParams, PROVENANCE_ESTIMATE, PROVENANCE_LIBRARY,
                              classify_provenance, compute_distrust, rank_candidates_by_distrust,
                              select_by_distrust)
from t3.pdep.pes_loop import _trim_candidates
from t3.schema import PESLoopConfig


# The retained r002 CHO2 round-0 surface, as objects: two [H] + O=C=O bimolecular entrance channels
# with low saddles (TS1, TS3) and one high isomerization saddle (TS2). E0 in kJ/mol, exactly as the
# reduced network file declares (E0(TS) = written value with the 191.282 kJ/mol correction removed).
_ROUND0_E0 = PDepNetworkE0(
    species={'O=[C]O(8)': -184.487, '[O]C=O(1)': -169.935, '[H](6)': 211.805, 'O=C=O(5)': -403.087},
    transition_states={'TS1': 21.347 - 191.282, 'TS2': 135.319 - 191.282, 'TS3': 30.706 - 191.282})

_FAMILY_COMMENT = ('Estimated using template [Cdd_Od;HJ] for rate rule [CO2;HJ]\n'
                   'family: R_Addition_MultipleBond')


def _rxn(label, reactants, products, ts, comment=_FAMILY_COMMENT, var=None):
    return PDepPathReaction(label=label, reactants=reactants, products=products,
                            transition_state=ts, kinetics_type='Arrhenius',
                            kinetics_comment=comment, kinetics_uncertainty_var=var)


def _cand(reaction, ts, measurable=True):
    return QMCandidate(path_reaction=reaction, ts_label=ts, family='R_Addition_MultipleBond',
                       e0_sensitivity_measurable=measurable)


def _round0_candidates():
    return (
        _cand(_rxn('reaction3', ('[H](6)', 'O=C=O(5)'), ('[O]C=O(1)',), 'TS1'), 'TS1'),
        _cand(_rxn('reaction5', ('[O]C=O(1)',), ('O=[C]O(8)',), 'TS2',
                   comment='family: intra_H_migration'), 'TS2'),
        _cand(_rxn('reaction8', ('[H](6)', 'O=C=O(5)'), ('O=[C]O(8)',), 'TS3'), 'TS3'),
    )


class TestClassifyProvenance:

    def test_family_estimate_is_distrusted(self):
        assert classify_provenance(_FAMILY_COMMENT) == PROVENANCE_ESTIMATE

    def test_reaction_library_is_trusted(self):
        comment = "Fitted to 50 data points; dA = *|/ 12.9Reaction library: 'kineticsjobs'"
        assert classify_provenance(comment) == PROVENANCE_LIBRARY

    def test_empty_comment_fails_toward_computing(self):
        # Unknown provenance is treated as distrusted, not trusted -- we do not trust by default.
        assert classify_provenance('') == PROVENANCE_ESTIMATE


class TestComputeDistrust:

    def test_barrier_is_height_over_the_lower_adjacent_configuration(self):
        lowest = min(_ROUND0_E0.transition_states.values())
        rxn = _rxn('reaction3', ('[H](6)', 'O=C=O(5)'), ('[O]C=O(1)',), 'TS1')
        score = compute_distrust(rxn, 'TS1', _ROUND0_E0, lowest, DistrustParams())
        # [H] + O=C=O asymptote = 211.805 - 403.087 = -191.282; TS1 = -169.935; barrier = 21.347.
        assert abs(score.barrier_kj - 21.347) < 1e-6
        assert abs(score.height_above_lowest_saddle_kj - 0.0) < 1e-6
        assert score.in_window is True

    def test_high_saddle_is_outside_the_window(self):
        lowest = min(_ROUND0_E0.transition_states.values())
        rxn = _rxn('reaction5', ('[O]C=O(1)',), ('O=[C]O(8)',), 'TS2',
                   comment='family: intra_H_migration')
        score = compute_distrust(rxn, 'TS2', _ROUND0_E0, lowest, DistrustParams(energy_window_kj=30.0))
        assert score.height_above_lowest_saddle_kj > 30.0
        assert score.in_window is False

    def test_missing_ts_e0_is_kept_not_dropped(self):
        rxn = _rxn('rX', ('A',), ('B',), 'TS_missing')
        score = compute_distrust(rxn, 'TS_missing', _ROUND0_E0,
                                 min(_ROUND0_E0.transition_states.values()), DistrustParams())
        assert score.barrier_kj is None
        assert score.height_above_lowest_saddle_kj is None
        assert score.in_window is True  # cannot place it -> do not rule it out

    def test_variance_raises_the_score(self):
        lowest = min(_ROUND0_E0.transition_states.values())
        params = DistrustParams()
        without = compute_distrust(_rxn('r', ('[H](6)', 'O=C=O(5)'), ('[O]C=O(1)',), 'TS1'),
                                   'TS1', _ROUND0_E0, lowest, params)
        with_var = compute_distrust(
            _rxn('r', ('[H](6)', 'O=C=O(5)'), ('[O]C=O(1)',), 'TS1', var=25.73),
            'TS1', _ROUND0_E0, lowest, params)
        assert with_var.kinetics_var == 25.73
        assert without.kinetics_var is None
        assert with_var.score > without.score

    def test_estimate_outranks_library_at_equal_geometry(self):
        lowest = min(_ROUND0_E0.transition_states.values())
        params = DistrustParams()
        estimate = compute_distrust(_rxn('r', ('[H](6)', 'O=C=O(5)'), ('[O]C=O(1)',), 'TS1'),
                                    'TS1', _ROUND0_E0, lowest, params)
        library = compute_distrust(
            _rxn('r', ('[H](6)', 'O=C=O(5)'), ('[O]C=O(1)',), 'TS1',
                 comment="Reaction library: 'kineticsjobs'"),
            'TS1', _ROUND0_E0, lowest, params)
        assert estimate.provenance == PROVENANCE_ESTIMATE
        assert library.provenance == PROVENANCE_LIBRARY
        assert estimate.score > library.score


class TestRankCandidates:

    def test_lower_barrier_ranks_first_among_equal_provenance(self):
        eligible, declined = rank_candidates_by_distrust(_round0_candidates(), _ROUND0_E0,
                                                         DistrustParams(energy_window_kj=30.0))
        order = [cand.ts_label for cand, _ in eligible]
        assert order == ['TS1', 'TS3']  # TS1 barrier 21.3 < TS3 barrier 30.7
        assert [cand.ts_label for cand, _ in declined] == ['TS2']

    def test_window_declines_the_isomerization_saddle(self):
        _, declined = rank_candidates_by_distrust(_round0_candidates(), _ROUND0_E0,
                                                  DistrustParams(energy_window_kj=30.0))
        assert [cand.ts_label for cand, _ in declined] == ['TS2']


class TestSelectByDistrust:

    def _evidence(self, ts2_coeff):
        # TS1/TS3 report a structural zero; TS2 reports whatever we pass. delta_ln_k unused by the
        # distrust selector -- present only because capture records the coefficient.
        return {'TS1': (1.27e-18, 2.9e-14), 'TS2': (ts2_coeff, abs(ts2_coeff) * 8368.0),
                'TS3': (1.27e-18, 2.9e-14)}

    def test_selects_the_entrance_channels_the_floor_cannot_reach(self):
        split = CandidateSplit(candidates=_round0_candidates())
        out = select_by_distrust(split, self._evidence(-1.43e-05), _ROUND0_E0,
                                 DistrustParams(energy_window_kj=30.0))
        assert [c.ts_label for c in out.candidates] == ['TS1', 'TS3']
        assert all(c.distrust_score is not None for c in out.candidates)
        declined = [s for s in out.skipped if s.classification == SKIP_OUTSIDE_WINDOW]
        assert [s.ts_label for s in declined] == ['TS2']

    def test_structural_zero_does_not_change_the_ranking(self):
        # The whole point: distrust ignores the sensitivity value, so flipping TS2's coefficient
        # from a real number to an exact structural zero leaves the selection identical.
        split = CandidateSplit(candidates=_round0_candidates())
        real = select_by_distrust(split, self._evidence(-1.43e-05), _ROUND0_E0, DistrustParams())
        zero = select_by_distrust(split, self._evidence(0.0), _ROUND0_E0, DistrustParams())
        assert [c.ts_label for c in real.candidates] == [c.ts_label for c in zero.candidates]
        assert [s.classification for s in real.skipped] == [s.classification for s in zero.skipped]

    def test_measured_coefficient_is_carried_for_capture(self):
        split = CandidateSplit(candidates=_round0_candidates())
        out = select_by_distrust(split, self._evidence(-1.43e-05), _ROUND0_E0, DistrustParams())
        ts1 = next(c for c in out.candidates if c.ts_label == 'TS1')
        # The structural-zero coefficient stays visible on the queued candidate (I-031), never
        # clamped, even though distrust -- not it -- justified the selection.
        assert ts1.coefficient == 1.27e-18

    def test_does_not_select_everything(self):
        # Negative control: with a high saddle on the surface, distrust declines it.
        split = CandidateSplit(candidates=_round0_candidates())
        out = select_by_distrust(split, self._evidence(-1.43e-05), _ROUND0_E0, DistrustParams())
        assert len(out.candidates) < len(split.candidates)
        assert any(s.classification == SKIP_OUTSIDE_WINDOW for s in out.skipped)

    def test_library_candidate_is_skipped_not_queued(self):
        # A trusted library value that is IN the window must still be skipped, not queued: provenance
        # gates here, it does not merely lower the rank (Copilot PR #209 finding).
        lib = _cand(_rxn('reaction3', ('[H](6)', 'O=C=O(5)'), ('[O]C=O(1)',), 'TS1',
                         comment="Reaction library: 'kineticsjobs'"), 'TS1')
        out = select_by_distrust(CandidateSplit(candidates=(lib,)), {'TS1': (1.0, 0.5)},
                                 _ROUND0_E0, DistrustParams())
        assert out.candidates == ()
        assert [s.classification for s in out.skipped] == [SKIP_TRUSTED_PROVENANCE]
        assert out.skipped[0].ts_label == 'TS1'

    def test_library_skipped_while_estimate_selected(self):
        lib = _cand(_rxn('reaction3', ('[H](6)', 'O=C=O(5)'), ('[O]C=O(1)',), 'TS1',
                         comment="Reaction library: 'kineticsjobs'"), 'TS1')
        split = CandidateSplit(candidates=(lib,) + _round0_candidates()[1:])
        out = select_by_distrust(split, self._evidence(-1.43e-05), _ROUND0_E0, DistrustParams())
        # TS1 trusted-skip, TS2 out-of-window, only TS3 (an in-window estimate) queued.
        assert [c.ts_label for c in out.candidates] == ['TS3']
        assert any(s.classification == SKIP_TRUSTED_PROVENANCE and s.ts_label == 'TS1'
                   for s in out.skipped)

    def test_candidate_without_a_finite_row_is_skipped_not_invented(self):
        split = CandidateSplit(candidates=_round0_candidates())
        evidence = {'TS1': (1.27e-18, 2.9e-14), 'TS3': (1.27e-18, 2.9e-14)}  # TS2 absent
        out = select_by_distrust(split, evidence, _ROUND0_E0, DistrustParams())
        no_row = [s for s in out.skipped if s.ts_label == 'TS2']
        assert len(no_row) == 1
        # TS2 is unimolecular on both sides, so its E0 sensitivity IS structurally measurable ->
        # a genuinely missing row is SKIP_NO_EVIDENCE, not SKIP_UNMEASURABLE.
        assert no_row[0].classification == SKIP_NO_EVIDENCE

    def test_unmeasurable_candidate_without_a_row_is_classified_unmeasurable(self):
        rxn = _rxn('reaction3', ('[H](6)', 'O=C=O(5)'), ('[O]C=O(1)',), 'TS1')
        split = CandidateSplit(candidates=(_cand(rxn, 'TS1', measurable=False),))
        out = select_by_distrust(split, {}, _ROUND0_E0, DistrustParams())
        assert len(out.candidates) == 0
        assert out.skipped[0].classification == SKIP_UNMEASURABLE


class TestTrimUnderDistrustScope:
    """`_trim_candidates` must keep the distrust order select_by_distrust imposed, not re-sort."""

    def _distrust_config(self, limit=10):
        return PESLoopConfig(pes={'network': '/abs/n.py', 'source': ['HOCHO'],
                                  'bath_gas': {'He': 1.0}},
                             qm={'scope': 'distrust', 'max_transition_states_per_round': limit})

    def test_distrust_order_is_preserved_not_resorted_by_coefficient(self):
        # Feed candidates whose distrust order (TS1 then TS3) is the REVERSE of what a
        # |coefficient| sort would produce -- so a wrongly-wired trim that re-sorts is caught.
        split = CandidateSplit(candidates=_round0_candidates())
        ranked = select_by_distrust(split, {'TS1': (1e-18, 0.0), 'TS2': (-1.43e-05, 0.12),
                                            'TS3': (5e-05, 0.42)},
                                    _ROUND0_E0, DistrustParams()).candidates
        assert [c.ts_label for c in ranked] == ['TS1', 'TS3']  # distrust order (barrier), not coeff
        trimmed = _trim_candidates(ranked, self._distrust_config(), logger=None)
        assert [c.ts_label for c in trimmed] == ['TS1', 'TS3']

    def test_limit_takes_the_most_distrusted(self):
        split = CandidateSplit(candidates=_round0_candidates())
        ranked = select_by_distrust(split, {'TS1': (1e-18, 0.0), 'TS2': (-1.43e-05, 0.12),
                                            'TS3': (1e-18, 0.0)},
                                    _ROUND0_E0, DistrustParams()).candidates
        trimmed = _trim_candidates(ranked, self._distrust_config(limit=1), logger=None)
        assert [c.ts_label for c in trimmed] == ['TS1']


class TestVarParsing:

    def test_rate_uncertainty_var_is_read(self):
        text = (
            "reaction(\n"
            "    label = 'r1', reactants = ['A'], products = ['B'], transitionState = 'TS1',\n"
            "    kinetics = Arrhenius(A=(1.0,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'),\n"
            "        uncertainty=RateUncertainty(mu=0.08, var=25.73102631834909, Tref=1000.0, N=5,\n"
            "        data_mean=0.0, correlation='x'), comment='''Estimated from node x.'''),\n)\n"
            "network(label='n', isomers=['A'], reactants=[], bathGas={'Ar': 1})\n")
        network = parse_pdep_network_text(text, network_id='n1', path='/abs/n1.py')
        assert abs(network.path_reactions[0].kinetics_uncertainty_var - 25.73102631834909) < 1e-9

    def test_absent_uncertainty_is_none(self):
        text = (
            "reaction(\n"
            "    label = 'r1', reactants = ['A'], products = ['B'], transitionState = 'TS1',\n"
            "    kinetics = Arrhenius(A=(1.0,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'),\n"
            "        comment='''family: R_Addition_MultipleBond'''),\n)\n"
            "network(label='n', isomers=['A'], reactants=[], bathGas={'Ar': 1})\n")
        network = parse_pdep_network_text(text, network_id='n1', path='/abs/n1.py')
        assert network.path_reactions[0].kinetics_uncertainty_var is None


class TestEndToEndRound0Text:
    """Parse an inline reduced-network text (the round-0 topology) and select through the real seam."""

    _NETWORK_TEXT = """
species(label='O=[C]O(8)', E0=(-184.487,'kJ/mol'))
species(label='[O]C=O(1)', E0=(-169.935,'kJ/mol'))
species(label='[H](6)', E0=(211.805,'kJ/mol'))
species(label='O=C=O(5)', E0=(-403.087,'kJ/mol'))
transitionState(label='TS1', E0=(21.347 - 191.282, 'kJ/mol'))
transitionState(label='TS2', E0=(135.319 - 191.282, 'kJ/mol'))
transitionState(label='TS3', E0=(30.706 - 191.282, 'kJ/mol'))
reaction(label='reaction3', reactants=['[H](6)', 'O=C=O(5)'], products=['[O]C=O(1)'],
         transitionState='TS1',
         kinetics=Arrhenius(A=(1.0,'m^3/(mol*s)'), n=0, Ea=(21.3,'kJ/mol'), T0=(1,'K'),
                            comment='''family: R_Addition_MultipleBond'''))
reaction(label='reaction5', reactants=['[O]C=O(1)'], products=['O=[C]O(8)'],
         transitionState='TS2',
         kinetics=Arrhenius(A=(1.0,'s^-1'), n=0, Ea=(114.0,'kJ/mol'), T0=(1,'K'),
                            comment='''family: intra_H_migration'''))
reaction(label='reaction8', reactants=['[H](6)', 'O=C=O(5)'], products=['O=[C]O(8)'],
         transitionState='TS3',
         kinetics=Arrhenius(A=(1.0,'cm^3/(mol*s)'), n=0, Ea=(30.7,'kJ/mol'), T0=(1,'K'),
                            comment='''family: R_Addition_MultipleBond'''))
network(label='PDepNetwork #1', isomers=['O=[C]O(8)', '[O]C=O(1)'], reactants=[], bathGas={'Ar': 1})
"""

    def test_parsed_network_selects_entrance_channels(self):
        network = parse_pdep_network_text(self._NETWORK_TEXT, network_id='n0', path='/abs/n0.py')
        e0 = parse_pdep_network_e0_text(self._NETWORK_TEXT, path='/abs/n0.py')
        split = CandidateSplit(candidates=tuple(
            QMCandidate(path_reaction=pr, ts_label=pr.transition_state, family=None)
            for pr in network.path_reactions))
        evidence = {'TS1': (1.27e-18, 2.9e-14), 'TS2': (-1.43e-05, 0.12), 'TS3': (1.27e-18, 2.9e-14)}
        out = select_by_distrust(split, evidence, e0, DistrustParams(energy_window_kj=30.0))
        # Both [H] + O=C=O entrance channels selected; the isomerization saddle declined.
        assert [c.ts_label for c in out.candidates] == ['TS1', 'TS3']
        assert [s.ts_label for s in out.skipped if s.classification == SKIP_OUTSIDE_WINDOW] == ['TS2']
