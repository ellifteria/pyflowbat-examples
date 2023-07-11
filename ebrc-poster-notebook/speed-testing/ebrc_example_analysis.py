import pyflowbat as pfb
import numpy as np
import time

def main():
   
    import pyflowbat as pfb

    
    my_wrkspc = pfb.pyflowbat.Workspace(full_output=False)

    
    my_wrkspc.calculate_beads_factors(
        beads_file_file_path="./ebrc-poster-example-data/Beads_After.fcs",
        beads_fluorescent_channels=[
            ("PE-Texas Red-A", "MEPTRs"),
            ("FITC-A", "MEFLs")
        ],
        beads_num_pops=9
    )

    
    my_wrkspc.load_samples(
        sample_collection_name="raw",
        samples_folder_path="./ebrc-poster-example-data/",
        include=['.fcs'],
        not_include=['Beads']
    )

    
    my_wrkspc.apply_gate(
        sample_collection_name='raw',
        new_sample_collection_name='heks',
        gating_function=pfb.gating.gate_heks,
        method = 'same',
        samples = ['373_C_001.fcs', '664_D_002.fcs', '373_M_001.fcs']
    )

    
    my_wrkspc.apply_gate(
        sample_collection_name='heks',
        new_sample_collection_name='singlets',
        gating_function=pfb.gating.gate_singlets,
        gating_channel_names=["FSC-A", "FSC-H"],
        b=1*10**4
    )

    
    my_wrkspc.calculate_compensation_matrix(
        sample_collection_name="singlets",
        compensation_sample_names=['Colors_DsRE2.fcs', 'Colors_mNG.fcs'],
        compensation_channel_names=['PE-Texas Red-A', 'FITC-A'],
        threshold=10**-3
    )

    
    my_wrkspc.apply_compensation_matrix(
        sample_collection_name="singlets",
        new_sample_collection_name="compensated"
    )

    
    my_wrkspc.apply_gate(
        sample_collection_name="compensated",
        new_sample_collection_name="transfected",
        gating_function=pfb.gating.gate_high_low,
        gating_channel_name="Alexa 750-A",
        gating_channel_names=["Alexa 750-A","Alexa 750-A"],
        low=pfb.gating.find_percentile(
            workspace=my_wrkspc,
            sample_collection_name="compensated",
            sample_name="Colors_Filler.fcs",
            channel_name="Alexa 750-A",
            percentile=99
        )
    )

    
    extraction = my_wrkspc.create_statistic_extraction(
        sample_collection_name='transfected',
        statistics_collection_name='samples',
        include=['.fcs'],
        not_include=['Beads', 'Controls', 'Colors'],
        statistic_names=['line', 'dox_conc', 'sample','Mean FITC-A', 'Mean PE-Texas Red-A']
    )
    my_wrkspc.extract_statistic(
        extraction=extraction,
        statistc_name='line',
        operation=pfb.operations.split_sample_name,
        by = '_',
        index = 0
    )
    my_wrkspc.extract_statistic(
        extraction=extraction,
        statistc_name='dox_conc',
        operation=pfb.operations.split_sample_name,
        by = '_',
        index = 1
    )
    my_wrkspc.extract_statistic(
        extraction=extraction,
        statistc_name='sample',
        operation=pfb.operations.split_sample_name,
        by = '_',
        index = 2
    )
    my_wrkspc.extract_statistic(
        extraction=extraction,
        statistc_name='Mean FITC-A',
        operation=pfb.operations.channel_mean,
        channel_name = 'FITC-A' 
    )
    my_wrkspc.extract_statistic(
        extraction=extraction,
        statistc_name='Mean PE-Texas Red-A',
        operation=pfb.operations.channel_mean,
        channel_name = 'PE-Texas Red-A'
    )

    
    my_wrkspc.combine_replicates(
        statistics_collection_name='samples',
        combined_statistics_collection_name='samples combined',
        combine_by=['dox_conc', 'line'],
        combination_operations={
            'Mean FITC-A': 'mean',
            'Mean PE-Texas Red-A': 'mean'
        },
        sem_cols=[
        'Mean FITC-A',
        'Mean PE-Texas Red-A' 
        ]
    )

    
    my_wrkspc.apply_operation(
        statistics_collection_name="samples combined",
        new_statistics_collection_name="samples converted",
        statistic_name="Mean FITC-A",
        new_statistic_name="Mean FITC-A",
        operation=pfb.operations.non_negative
    )
    my_wrkspc.apply_operation(
        statistics_collection_name="samples combined",
        new_statistics_collection_name="samples converted",
        statistic_name="Mean PE-Texas Red-A",
        new_statistic_name="Mean PE-Texas Red-A",
        operation=pfb.operations.non_negative
    )
    my_wrkspc.apply_operation(
        statistics_collection_name="samples converted",
        new_statistics_collection_name="samples converted",
        statistic_name="Mean FITC-A",
        new_statistic_name="MEFLs",
        operation=pfb.operations.apply_conversion_factor,
        factor=my_wrkspc.conversion_factors["FITC-A"]
    )
    my_wrkspc.apply_operation(
        statistics_collection_name="samples converted",
        new_statistics_collection_name="samples converted",
        statistic_name=["Mean FITC-A", "Mean FITC-A_stdErr"],
        new_statistic_name="MEFLs_stdErr",
        operation=pfb.operations.compute_conversion_factor_stdErr,
        factor=my_wrkspc.conversion_factors["FITC-A"],
        factor_err=my_wrkspc.conversion_factors["FITC-A_stderr"]
    )
    my_wrkspc.apply_operation(
        statistics_collection_name="samples converted",
        new_statistics_collection_name="samples converted",
        statistic_name="Mean PE-Texas Red-A",
        new_statistic_name="MEPTRs",
        operation=pfb.operations.apply_conversion_factor,
        factor=my_wrkspc.conversion_factors["PE-Texas Red-A"]
    )
    my_wrkspc.apply_operation(
        statistics_collection_name="samples converted",
        new_statistics_collection_name="samples converted",
        statistic_name=["Mean PE-Texas Red-A", "Mean PE-Texas Red-A_stdErr"],
        new_statistic_name="MEPTRs_stdErr",
        operation=pfb.operations.compute_conversion_factor_stdErr,
        factor=my_wrkspc.conversion_factors["PE-Texas Red-A"],
        factor_err=my_wrkspc.conversion_factors["PE-Texas Red-A_stderr"]
    )



if __name__ == "__main__":
    times = []
    num_runs = 20
    for i in range(num_runs):
        start_time = time.time()
        main()
        execution_time = time.time() - start_time
        print("run: %d" %i)
        print("--- %s seconds ---" % (execution_time))
        times.append(execution_time)

    print("\nmean execution_time:")
    print("--- %s seconds ---" % (np.mean(times)))
    print("\std of execution_times:")
    print("--- %s seconds ---" % (np.std(times)))
