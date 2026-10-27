from __future__ import annotations
import re


# ── Normalization helper ───────────────────────────────────────────────────────
def normalize(s: str) -> str:
    """Lowercase, strip separators/spaces for fuzzy-free key lookup."""
    return re.sub(r'[-_/\s]', '', s).lower()


# ── Canonical mapping  ────────────────────────────────────────────────────────
# Keys  : normalize(raw_name)  →  all lowercase, no separators
# Values: chosen canonical display name
#
# Guiding principles:
#   - Prefer the shortest unambiguous common abbreviation
#   - Keep phosphosite notation where it disambiguates (e.g. pRb vs Rb)
#   - Clones in parentheses (E1L3N, 28-8, SP142) are KEPT — they are
#     biologically distinct antibodies, not just naming variants
#   - R11c2 / R11c4 are uncharacterised controls in TMA-7; kept as-is
# ─────────────────────────────────────────────────────────────────────────────

CANONICAL: dict[str, str] = {

    # DAPI
    normalize('DAPI1'): 'dna(1)',
    normalize('DAPI2'): 'dna(2)',
    normalize('DAPI3'): 'dna(3)',
    normalize('DAPI4'): 'dna(4)',
    normalize('DAPI5'): 'dna(5)',
    normalize('DAPI6'): 'dna(6)',
    normalize('DAPI7'): 'dna(7)',
    normalize('DAPI8'): 'dna(8)',
    normalize('DAPI9'): 'dna(9)',
    normalize('DAPI10'): 'dna(10)',
    normalize('DAPI11'): 'dna(11)',

    # ── Nuclear / DNA-damage ─────────────────────────────────────────────────
    normalize('53BP1'):                       '53BP1',
    normalize('Histone H2AX'):               'gH2AX',      # H2AX phospho-S139
    normalize('gH2AX'):                       'gH2AX',
    normalize('hRAD51'):                      'RAD51',
    normalize('RAD51'):                       'RAD51',
    normalize('PCNA'):                        'PCNA',
    normalize('pHH3'):                        'pHH3',       # phospho-Histone H3
    normalize('Lamin-A/C'):                   'LamAC',
    normalize('Lamin-A/B/C'):                 'LamAC',
    normalize('LamAC'):                       'LamAC',
    normalize('LamB1'):                       'LamB1',      # Lamin B1 — distinct protein
    normalize('H3K4'):                        'H3K4',       # H3K4me
    normalize('H3K27'):                       'H3K27',      # H3K27me3
    normalize('MSH2'):                        'MSH2',
    normalize('Pin1'):                        'Pin1',

    # ── Cell cycle ───────────────────────────────────────────────────────────
    normalize('Antigen Ki67'):                'Ki67',
    normalize('Ki67'):                        'Ki67',
    normalize('G1/S-specific cyclin-D1'):     'Cyclin_D1',
    normalize('CCND1'):                       'Cyclin_D1',
    normalize('Rb (pS807; pS811)'):           'pRB',
    normalize('pRB'):                         'pRB',
    normalize('p21'):                         'p21',
    normalize('p53'):                         'p53',
    normalize('pRPA'):                        'pRPA',

    # ── Signalling / kinases ─────────────────────────────────────────────────
    normalize('ERK-1 (pT202; pY204)'):        'pERK',
    normalize('pERK'):                        'pERK',
    normalize('pERK-PE'):                     'pERK',       # same target, different fluorophore
    normalize('pAKT'):                        'pAKT',
    normalize('pS6RP'):                       'pS6RP',
    normalize('pMYC'):                        'pMYC',
    normalize('EGFR'):                        'EGFR',
    normalize('Epidermal growth factor receptor'): 'EGFR',
    normalize('PDGFRa'):                      'PDGFRa',
    normalize('CSF1R'):                       'CSF1R',
    normalize('RNA polymerase II CTD (pS2)'): 'RNApol2-pS2',
    normalize('BCL2'):                        'BCL2',
    normalize('CAV1'):                        'CAV1',
    normalize('BMP2'):                        'BMP2',
    normalize('BMP6'):                        'BMP6',
    normalize('ZEB1'):                        'ZEB1',

    # ── Hormone receptors / breast markers ───────────────────────────────────
    normalize('ER'):                          'ER',
    normalize('PR'):                          'PR',
    normalize('PgR'):                         'PR',
    normalize('HER2'):                        'HER2',
    normalize('AR'):                          'AR',
    normalize('Androgen Receptor'):           'AR',
    normalize('GATA3'):                       'GATA3',
    normalize('MUC1'):                        'MUC1',

    # ── Epithelial / structural ───────────────────────────────────────────────
    normalize('E-cadherin'):                  'Ecad',
    normalize('Ecad'):                        'Ecad',
    normalize('Pan-cytokeratin'):             'panCK',
    normalize('panCK'):                       'panCK',
    normalize('Vimentin'):                    'Vim',
    normalize('Vim'):                         'Vim',
    normalize('VIM'):                         'Vim',
    normalize('Aortic smooth muscle actin'):  'aSMA',
    normalize('Aortic smooth muscle actin (2)'): 'aSMA',   # duplicate channel in LSP12022
    normalize('aSMA'):                        'aSMA',
    normalize('Podoplanin'):                  'PDPN',
    normalize('PDPN'):                        'PDPN',
    normalize('Podoplanin'):                   'PDPN',
    normalize('TUBB3'):                       'TUBB3',
    normalize('Glut1'):                       'Glut1',
    normalize('CGA'):                         'CGA',        # Chromogranin A

    # ── Cytokeratins ─────────────────────────────────────────────────────────
    normalize('CK-5'):                        'CK5',
    normalize('CK5'):                         'CK5',
    normalize('CK-14'):                       'CK14',
    normalize('CK14'):                        'CK14',
    normalize('CK-17'):                       'CK17',
    normalize('CK17'):                        'CK17',
    normalize('CK-19'):                       'CK19',
    normalize('CK19'):                        'CK19',
    normalize('CK1-4'):                       'CK1-4',      # AE1/AE3 cocktail — distinct from CK5/14/17/19
    normalize('CK7'):                         'CK7',
    normalize('CK8'):                         'CK8',

    # ── Immune / T-cell ───────────────────────────────────────────────────────
    normalize('CD3'):                         'CD3',
    normalize('CD3d'):                        'CD3d',
    normalize('CD4'):                         'CD4',
    normalize('CD8'):                         'CD8',
    normalize('CD8a'):                        'CD8',        # CD8a = CD8 in most panels
    normalize('CD11b'):                       'CD11b',
    normalize('CD11c'):                       'CD11c',
    normalize('CD14'):                        'CD14',
    normalize('CD15'):                        'CD15',
    normalize('CD20'):                        'CD20',
    normalize('CD31'):                        'CD31',
    normalize('CD44'):                        'CD44',
    normalize('CD45'):                        'CD45',
    normalize('CD57'):                        'CD57',
    normalize('CD68'):                        'CD68',
    normalize('CD90'):                        'CD90',
    normalize('CD163'):                       'CD163',
    normalize('FOXP3'):                       'FOXP3',
    normalize('FoxP3'):                       'FOXP3',
    normalize('Granzyme B'):                  'GranzymeB',
    normalize('GRNZB'):                       'GranzymeB',
    normalize('LAG-3'):                       'LAG3',
    normalize('MHC class II antigen DRA'):    'HLA-DRA',
    normalize('MHC class II antigen DPB1'):   'HLA-DPB1',
    normalize('HLA-A'):                       'HLA-A',

    # ── Checkpoint / immune evasion ───────────────────────────────────────────
    normalize('PD-1'):                        'PD1',
    normalize('PD1'):                         'PD1',
    normalize('PD-L1 (E1L3N)'):              'PDL1-E1L3N',   # distinct clones!
    normalize('PD-L1 (28-8)'):               'PDL1-28-8',
    normalize('PD-L1 (SP142)'):              'PDL1-SP142',
    normalize('PDL1'):                        'PDL1',         # when clone is unknown

    # ── Apoptosis ─────────────────────────────────────────────────────────────
    normalize('CC3'):                         'CC3',          # cleaved Caspase-3
    normalize('cPARP'):                       'cPARP',        # cleaved PARP

    # ── ECM / stromal ─────────────────────────────────────────────────────────
    normalize('ColI'):                        'ColI',
    normalize('ColIV'):                       'ColIV',
    normalize('CoxIV'):                       'CoxIV',        # mitochondrial — distinct from ColIV

    # ── Misc / metabolic ──────────────────────────────────────────────────────
    normalize('R11c2'):                       'R11c2',        # uncharacterised TMA-7 control
    normalize('R11c4'):                       'R11c4',

}


# ── Public API ────────────────────────────────────────────────────────────────

def canonicalize(name: str) -> str:
    """
    Return the canonical name for a raw channel name.
    Falls back to the normalize(name) form if not in the dict,
    so unknown markers are still consistently cased/formatted.
    """
    return CANONICAL.get(normalize(name), normalize(name))


def validate_panel(channel_dict: dict, strict: bool = False) -> list[str]:
    """
    Check a {channel_name: index} dict against the canonical map.
    Returns a list of warning strings for unmapped names.
    Raises ValueError if strict=True and any are unmapped.
    """
    skip_prefixes = ('dna', 'control', 'autofluorescence')
    warnings = []
    for name in channel_dict:
        if any(normalize(name).startswith(p) for p in skip_prefixes):
            continue
        if normalize(name) not in CANONICAL:
            warnings.append(f"UNMAPPED: '{name}'  (normalized: '{normalize(name)}')")
    if strict and warnings:
        raise ValueError('\n'.join(warnings))
    return warnings


# ── Quick smoke-test ──────────────────────────────────────────────────────────
if __name__ == '__main__':
    tests = [
        ('Antigen Ki67',               'Ki67'),
        ('Ki67',                       'Ki67'),
        ('Pan-cytokeratin',            'panCK'),
        ('panCK',                      'panCK'),
        ('E-cadherin',                 'Ecad'),
        ('Ecad',                       'Ecad'),
        ('Vimentin',                   'Vim'),
        ('Vim',                        'Vim'),
        ('VIM',                        'Vim'),
        ('Aortic smooth muscle actin', 'aSMA'),
        ('aSMA',                       'aSMA'),
        ('Rb (pS807; pS811)',           'pRb'),
        ('pRB',                        'pRb'),
        ('G1/S-specific cyclin-D1',    'CyclinD1'),
        ('CCND1',                      'CyclinD1'),
        ('PD-1',                       'PD1'),
        ('PD1',                        'PD1'),
        ('PD-L1 (E1L3N)',              'PDL1-E1L3N'),
        ('PD-L1 (28-8)',               'PDL1-28-8'),
        ('FOXP3',                      'FOXP3'),
        ('FoxP3',                      'FOXP3'),
        ('Granzyme B',                 'GranzymeB'),
        ('GRNZB',                      'GranzymeB'),
        ('Epidermal growth factor receptor', 'EGFR'),
        ('ERK-1 (pT202; pY204)',       'pERK'),
        ('pERK',                       'pERK'),
        ('pERK-PE',                    'pERK'),
        ('Androgen Receptor',          'AR'),
        ('AR',                         'AR'),
        ('CD8a',                       'CD8'),
        ('CD8',                        'CD8'),
        ('CK-19',                      'CK19'),
        ('CK19',                       'CK19'),
        ('CK-5',                       'CK5'),
        ('Lamin-A/C',                  'LamAC'),
        ('Lamin-A/B/C',                'LamAC'),
        ('LamAC',                      'LamAC'),
        ('Histone H2AX',               'gH2AX'),
        ('gH2AX',                      'gH2AX'),
        ('hRAD51',                     'RAD51'),
        ('RAD51',                      'RAD51'),
        ('CoxIV',                      'CoxIV'),
        ('ColIV',                      'ColIV'),
        ('MHC class II antigen DRA',   'HLA-DRA'),
        ('Podoplanin',                 'PDPN'),
        ('PgR',                        'PR'),
    ]

    passed = failed = 0
    for raw, expected in tests:
        result = canonicalize(raw)
        status = '✓' if result == expected else '✗'
        if result != expected:
            print(f"{status} '{raw}' → '{result}'  (expected '{expected}')")
            failed += 1
        else:
            passed += 1

    print(f"\n{passed}/{passed+failed} tests passed")