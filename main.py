from data_layer.communities_stats import generate_per_sample_community_statistics
from data_layer.cross_experiment_arcos_comparison import generate_cross_experiment_spatial_significance_data
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)



if __name__ == '__main__':
    experiment_types = ['mid_third_wild_type', 'mid_third_zpg_RNAi', 'cbx_inhibitor_3.125', 'cbx_inhibitor_12.5', 'cbx_inhibitor_washout']
    generate_cross_experiment_spatial_significance_data()
    generate_per_sample_community_statistics(experiment_types)


    ## RUN ALL ANALYSES FOR FIG1
    from analysis.Fig1.analysis_B import main as analysis_1A
    from analysis.Fig1.analysis_B import main as analysis_1B
    from analysis.Fig1.analysis_C import main as analysis_1C
    from analysis.Fig1.analysis_F import main as analysis_1F
    from analysis.Fig1.analysis_G import main as analysis_1G
    from analysis.Fig1.analysis_H import main as analysis_1H
    from analysis.Fig1.analysis_I import main as analysis_1I
    analysis_1A()
    analysis_1B()
    analysis_1C()
    analysis_1F()
    analysis_1G()
    analysis_1H()
    analysis_1I()

    ## RUN ALL ANALYSES FOR FIG2
    from analysis.Fig2.analysis_A import main as analysis_2A
    from analysis.Fig2.analysis_B import main as analysis_2B
    from analysis.Fig2.analysis_C_D import main as analysis_2CD
    from analysis.Fig2.analysis_E import main as analysis_2E
    from analysis.Fig2.analysis_G import main as analysis_2G
    from analysis.Fig2.analysis_H import main as analysis_2H
    analysis_2A()
    analysis_2B()
    analysis_2CD()
    analysis_2E()
    analysis_2G()
    analysis_2H()

    ## RUN ALL ANALYSES FOR FIG3
    from analysis.Fig3.analysis_B_E import main as analysis_3BE
    from analysis.Fig3.analysis_C_F import main as analysis_3CF
    from analysis.Fig3.analysis_D_G import main as analysis_3DG
    analysis_3BE()
    analysis_3CF()
    analysis_3DG()


    ## RUN ALL ANALYSES FOR SI_FIG1
    from analysis.SI1.analysis_A import main as analysis_SIA
    from analysis.SI1.analysis_B import main as analysis_SIB
    analysis_SIA()
    analysis_SIB()

    ## RUN ALL ANALYSES FOR SI_FIG2
    from analysis.SI2.analysis_A_B import main as analysis_SIA_B
    from analysis.SI2.analysis_C import main as analysis_SIC
    analysis_SIA_B()
    analysis_SIC()

    ## RUN ALL ANALYSES FOR SI_FIG3
    from analysis.SI3.analysis import main as analysis_SI3
    analysis_SI3()

    ## RUN ALL ANALYSES FOR SI_FIG4 + SI_FIG7
    from analysis.SI4_SI7.analyses import main as analysis_SI4_SI7
    analysis_SI4_SI7()

    ## RUN ALL ANALYSES FOR SI_FIG5
    from analysis.SI5.analysis import main as analysis_SI5
    analysis_SI5()

    ## RUN ALL ANALYSES FOR SI_FIG6
    from analysis.SI6.analysis import main as analysis_SI6
    analysis_SI6()




