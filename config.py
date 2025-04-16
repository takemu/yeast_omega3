import os

root_dir = os.path.dirname(os.path.abspath(__file__))
model_xml = os.path.join(root_dir, 'data', 'yeast-GEM-9.0.2.xml')
r_biomass = 'R_r_2111'
min_growth = 1E-6  # h^-1
MAX_FLUX = 1000

substrate_symbols = {
    'ammonium':'NH4+',
    'D-galactose':'Gal',
    'D-glucose':'Glc',
    'ethanol':'EtOH',
    'glycerol':'Glyc',
    'L-alanine':'Ala',
    'L-arginine':'Arg',
    'L-asparagine':'Asn',
    'L-aspartate':'Asp',
    'L-glutamate':'Glu',
    'L-glutamine':'Gln',
    'L-glycine':'Gly',
    'L-isoleucine':'Ile',
    'L-leucine':'Leu',
    'L-lysine':'Lys',
    'L-methionine':'Met',
    'L-phenylalanine':'Phe',
    'L-serine':'Ser',
    'L-threonine':'Thr',
    'L-tryptophan':'Trp',
    'L-tyrosine':'Tyr',
    'L-valine':'Val'
}

substrate_rxns = {
    'NH4+': 'R_r_1654',
    'Gal': 'R_r_1710',
    'Glc': 'R_r_1714',
    'EtOH': 'R_r_1761',
    'Glyc': 'R_r_1808',
    'Ala': 'R_r_1873',
    'Arg': 'R_r_1879',
    'Asn': 'R_r_1880',
    'Asp': 'R_r_1881',
    'Glu': 'R_r_1889',
    'Gln': 'R_r_1891',
    'Gly': 'R_r_1810',
    'Ile': 'R_r_1897',
    'Leu': 'R_r_1899',
    'Lys': 'R_r_1900',
    'Met': 'R_r_1902',
    'Phe': 'R_r_1903',
    'Ser': 'R_r_1906',
    'Thr': 'R_r_1911',
    'Trp': 'R_r_1912',
    'Tyr': 'R_r_1913',
    'Val': 'R_r_1914'
}

ynb_medium = [substrate_rxns[sub] for sub in ['NH4+']]
yp_medium = [substrate_rxns[sub] for sub in ['Ala', 'Arg', 'Asn', 'Asp', 'Glu', 'Gln', 'Gly', 'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val']]

