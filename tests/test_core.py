import numpy as np, pandas as pd, pytest
from depclaims import Atlas, build


@pytest.fixture
def cache(tmp_path):
    rng = np.random.default_rng(0)
    lines = [f'ACH-{i:06d}' for i in range(12)]
    genes = {1: 'ESS1', 2: 'ESS2', 3: 'NON1', 4: 'NON2', 5: 'SEL1', 6: 'NEVER1'}
    eff = pd.DataFrame(0.0, index=lines, columns=[f'{s} ({g})' for g, s in genes.items()])
    eff['ESS1 (1)'] = -1.0; eff['ESS2 (2)'] = -1.0
    eff['SEL1 (5)'] = [-1.2] * 3 + [0.0] * 9
    eff['NEVER1 (6)'] = rng.normal(0, 0.05, 12)
    eff.index.name = 'ModelID'
    eff.to_csv(tmp_path / 'chronos.csv'); eff.to_csv(tmp_path / 'd2.csv')
    ex = np.log2(eff * 0 + 11); ex['NEVER1 (6)'] = 0.0
    ex.reset_index().to_csv(tmp_path / 'expr.csv', index=False)
    pd.DataFrame({'ModelID': lines, 'StrippedCellLineName': [f'LINE{i}' for i in range(12)], 'CellLineName': [f'Line-{i}' for i in range(12)],
                  'OncotreeLineage': ['Breast'] * 6 + ['Lung'] * 6, 'OncotreePrimaryDisease': ['x'] * 12}).to_csv(tmp_path / 'models.csv', index=False)
    pd.DataFrame({'ENTREZ_ID': [1, 2]}).to_csv(tmp_path / 'ceg.txt', sep='\t', index=False)
    pd.DataFrame({'ENTREZ_ID': [3, 4]}).to_csv(tmp_path / 'neg.txt', sep='\t', index=False)
    return build(tmp_path / 'chronos.csv', tmp_path / 'd2.csv', tmp_path / 'expr.csv', tmp_path / 'models.csv',
                 tmp_path / 'ceg.txt', tmp_path / 'neg.txt', tmp_path / 'cache')


def test_score_dependent_line(cache):
    s = Atlas(cache).score('SEL1', 'line-0')
    assert s.depmap_model == 'ACH-000000' and s.needed_in_line and s.gene_class == 'selective'
    assert s.crispr_effect_scaled == pytest.approx(-1.2) and s.crispr_rank_in_lineage < 0.25
    assert s.expression_tpm == pytest.approx(10.0, rel=1e-4)


def test_score_unneeded_gene(cache):
    s = Atlas(cache).score(6, 'LINE7')
    assert s.gene_class == 'never' and not s.needed_in_line and s.expression_tpm == pytest.approx(0.0)


def test_unknown_inputs(cache):
    A = Atlas(cache)
    with pytest.raises(KeyError): A.score('NOPE', 'LINE1')
    with pytest.raises(KeyError): A.score('SEL1', 'NOT A LINE')


def test_lineage_profile_sorted(cache):
    p = Atlas(cache).lineage_profile('SEL1', 'Breast')
    assert list(p.cell_line[:3]) == ['LINE0', 'LINE1', 'LINE2'] and p.crispr_effect_scaled.is_monotonic_increasing
