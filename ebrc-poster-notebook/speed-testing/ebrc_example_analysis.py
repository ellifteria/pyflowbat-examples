import pyflowbat as pfb
import numpy as np
import time

def main():
    my_wrkspc = pfb.pyflowbat.Workspace(full_output=False)

    my_wrkspc.calculate_beads_factors("./ebrc-poster-notebook/ebrc-poster-example-data/Beads_After.fcs", [["PE-Texas Red-A", "MEPTRs"], ["FITC-A", "MEFLs"]], 9)

    my_wrkspc.load_samples("raw", "./ebrc-poster-notebook/ebrc-poster-example-data/", lambda name: True if ('.fcs' in name and 'Beads' not in name) else False)

    my_wrkspc.apply_gate('raw', 'heks', pfb.gating.gate_heks, {'method': 'same', 'samples': ['373_C_001.fcs', '664_D_002.fcs', '373_M_001.fcs']}, gate_type=2)

    my_wrkspc.apply_gate('heks', 'singlets', pfb.gating.gate_singlets, {})

    my_wrkspc.calculate_compensation_matrix('singlets', ['Colors_DsRE2.fcs', 'Colors_mNG.fcs'], ['PE-Texas Red-A', 'FITC-A'], threshold=10**-3)

    my_wrkspc.apply_compensation_matrix("singlets", "compensated")

    my_wrkspc.apply_gate("compensated", "transfected", pfb.gating.gate_high_low, {"channel": "Alexa 750-A", "low": pfb.gating.find_percentile(my_wrkspc.sample_collections["compensated"]["Colors_Filler.fcs"], **{"channel": "Alexa 750-A", "percentile":99.5})})

    dox_conc = {
        'A':0,
        'B':10**0,
        'C':10**0.17,
        'D':10**0.33,
        'E':10**0.50,
        'F':10**0.63,
        'G':10**0.75,
        'H':10**0.88,
        'I':10**1.0,
        'J':10**1.08,
        'K':10**1.25,
        'L':10**1.33,
        'M':10**1.50,
        'N':10**1.75,
        'O':10**2.0,
        'P':10**3.0
    }

    my_wrkspc.extract_statistics('transfected', 'samples', lambda name: True if (('Beads' not in name and 'Controls' not in name and 'Colors' not in name) and '.fcs' in name) else False,
        [
            ['line', lambda file_name, fcs_data: (file_name.split("_"))[0]],
            ['sample', lambda file_name, fcs_data: (file_name.split("_"))[2]],
            ['dox_conc', lambda file_name, fcs_data: dox_conc[(file_name.split("_"))[1]]],
            ['Mean FITC-A', lambda file_name, fcs_data: np.mean(fcs_data[:, 'FITC-A'])],
            ['Mean PE-Texas Red-A', lambda file_name, fcs_data: np.mean(fcs_data[:, 'PE-Texas Red-A'])]
        ]
    )

    my_wrkspc.stats_collections["samples"]

    my_wrkspc.combine_replicates('samples', 'samples combined',
        lambda row: [row['dox_conc'], row['line']],
        [
            ['CellLine', lambda row: row['line'], False],
            ['DoxConc', lambda row: row['dox_conc'], False],
            ['Mean FITC-A', lambda row: row['Mean FITC-A'], True],
            ['Mean PE-Texas Red-A', lambda row: row['Mean PE-Texas Red-A'], True],
            ['MEFLs', lambda x: 0, True],
            ['MEPTRs', lambda x: 0, True]
        ]
    )

    my_wrkspc.apply_operation('samples combined', 'samples converted', 
        [
            ['Mean FITC-A', lambda row, inputs: max(row['Mean FITC-A'], 0)],
            ['Mean PE-Texas Red-A', lambda row, inputs: max(row['Mean PE-Texas Red-A'], 0)],
            ['MEFLs', lambda row, inputs: row['Mean FITC-A']*my_wrkspc.conversion_factors["FITC-A"]],
            ['MEFLs_stdErr', lambda row, inputs: row['Mean FITC-A']*my_wrkspc.conversion_factors["FITC-A"] * np.sqrt((row['Mean FITC-A_stdErr']/row['Mean FITC-A'])**2+(my_wrkspc.conversion_factors["FITC-A_stderr"]/my_wrkspc.conversion_factors["FITC-A"])**2)],
            ['MEPTRs', lambda row, inputs: row['Mean PE-Texas Red-A']*my_wrkspc.conversion_factors["PE-Texas Red-A"]],
            ['MEPTRs_stdErr', lambda row, inputs: row['Mean PE-Texas Red-A']*my_wrkspc.conversion_factors["PE-Texas Red-A"] * np.sqrt((row['Mean PE-Texas Red-A_stdErr']/row['Mean PE-Texas Red-A'])**2+(my_wrkspc.conversion_factors["PE-Texas Red-A_stderr"]/my_wrkspc.conversion_factors["PE-Texas Red-A"])**2)]
        ], []
    )

if __name__ == "__main__":
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))