###########
# Usage: python iAller.py
###########

import argparse
import numpy as np
from pubscripts import *
from dimreduction import pca as dimpca
from featureselection import *
from machinelearning import *





if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="it's usage tip.",
                                     description="running the pipeline")
    parser.add_argument("--config", default='config1', help="the config file")
    args = parser.parse_args()

    print('\n Step 1: Preprocessing sequence data')

    #1. reading config
    print('1. Reading config')
    print('\nYou might need to modify the path of Sequence_File in the config1 file.\n')
    # read config
    parameters = read_config.read_config(args.config)
    if parameters['Sequence_Type'] == 'Protein':
        from descproteins import *
    else:
        print("The file is not Protein!")

    # commands for encoding
    protein_cmd_coding = {
        'AAC': ['AAC.AAC(training_data, **kw)', 'AAC.AAC(testing_data, **kw)'],
        'CKSAAP': ['CKSAAP.CKSAAP(training_data, gap=%d, **kw)' % int(parameters['K_Space']), 'CKSAAP.CKSAAP(testing_data, gap=%d, **kw)' % int(parameters['K_Space'])],
        'DPC': ['DPC.DPC(training_data, **kw)', 'DPC.DPC(testing_data, **kw)'],
    }

    # Error information
    error_array = []


    #2.1. read fasta sequence and specify cmd
    print('2. Reading FASTA file of protein sequences')
    fastas = []
    cmd_coding = {}
    if parameters['Sequence_Type'] == 'Protein':
        fastas = read_fasta_sequences.read_protein_sequences(parameters['Sequence_File'])
        cmd_coding = protein_cmd_coding
    else:
        error_array.append('Sequence type can only be selected in "Protein".')

    kw = {'nclusters': 3, 'sof': 'sample', 'order': ''}
    kw['order'] = 'ACDEFGHIKLMNPQRSTVWY'

    #2.2. divide training and testing data
    training_data = []
    testing_data = []
    PSTNP_training_data = []
    PSTNP_testing_data = []
    for sequence in fastas:
        if sequence[3] == 'training':
            training_data.append(sequence)
            PSTNP_training_data.append(sequence)
            PSTNP_training_data.append([sequence[0], sequence[1], sequence[2], 'testing'])
            PSTNP_testing_data.append(sequence)
        else:
            testing_data.append(sequence)
            PSTNP_testing_data.append(sequence)


    #3. calculate descriptor 
    print('3. Calculating descriptors for training and testing data')
    training_code_dict = {}
    testing_code_dict = {}
    method_array = parameters['Method'].split(';')
    for method in method_array:
        # calculate descriptors for training data
        training_code_dict[method] = eval(cmd_coding[method][0])
        # calculate descriptors for testing data
        if len(testing_data) > 0:
            testing_code_dict[method] = eval(cmd_coding[method][1])

    training_code = np.array(training_code_dict[method_array[0]])
    testing_code = []
    if len(testing_data) > 0:
        testing_code = np.array(testing_code_dict[method_array[0]])

    for i in range(1, len(method_array)):
        if training_code_dict[method_array[i]] != 0:
            training_code = np.concatenate((training_code, np.array(training_code_dict[method_array[i]])[:, 2:]), axis=1)
            if len(testing_data) > 0:
                if testing_code_dict[method_array[i]] != 0:
                    testing_code = np.concatenate((testing_code, np.array(testing_code_dict[method_array[i]])[:, 2:]), axis=1)

    if len(testing_data) != 0 and training_code.shape[1] != testing_code.shape[1]:
        error_array.append('Descriptor(s) for testing data calculating failed.')
        testing_data = []

    training_code = training_code.tolist()
    save_file.save_file(training_code, format=parameters['Output_Format'], file='training_code.txt')
    save_file.save_file(training_code, format='tsv_1', file='training_code_1.tsv')

    if len(testing_data) > 0:
        testing_code = testing_code.tolist()
        save_file.save_file(testing_code, format=parameters['Output_Format'], file='testing_code.txt')
        save_file.save_file(testing_code, format='tsv_1', file='testing_code_1.tsv')


    #4. prepare data for feature selection and dimension reduction
    training_code_file = 'training_code.txt'
    testing_code_file = 'testing_code.txt'
    training_code_used, training_labels = read_code.read_code(training_code_file, format=parameters['Output_Format'])
    testing_code_used, testing_labels = [], []
    if len(testing_data) > 0:
        testing_code_used, testing_labels = read_code.read_code(testing_code_file, format=parameters['Output_Format'])

    #4.1. feature selection
    print('4.1. feature selection')
    training_code_selected = []
    testing_code_selected = []
    if parameters['Feature_Selection_Algorithm'] != '' and parameters['Feature_Selection_Algorithm'] in ('pearsonr'):
        training_code_used, training_labels = read_code.read_code(training_code_file,
                                                                  format=parameters['Output_Format'])
        if len(testing_data) > 0:
            testing_code_used, testing_labels = read_code.read_code(testing_code_file,
                                                                    format=parameters['Output_Format'])

        cmd = parameters['Feature_Selection_Algorithm'] + '.' + parameters[
            'Feature_Selection_Algorithm'] + '(training_code_used, training_labels)'
        result_tuple = eval(cmd)
        if isinstance(result_tuple, tuple) and len(result_tuple) > 1:
            selected_features, e = result_tuple[:2]  # only take the first two elements
        else:
            selected_features = result_tuple
            e = None

        save_file.save_FS_result(selected_features, e, parameters['Feature_Selection_Algorithm'], 'feaure_rank.txt')
        training_code_selected = select_features.select_features(training_code_used, training_labels,
                                                                 selected_features,
                                                                 parameters['Selected_Feature_Number']).tolist()
        save_file.save_file(training_code_selected, parameters['Output_Format'], 'training_code_selected.txt')
        training_code_file = 'training_code_selected.txt'

        if len(testing_data) > 0:
            testing_code_selected = select_features.select_features(testing_code_used, testing_labels,
                                                                    selected_features,
                                                                    parameters['Selected_Feature_Number']).tolist()
            save_file.save_file(testing_code_selected, parameters['Output_Format'], 'testing_code_selected.txt')
            testing_code_file = 'testing_code_selected.txt'

    #4.2. dimensionality reduction
    print('4.2. dimensionality reduction')
    if parameters['Dimension_Reduction_Algorithm'] != '' and parameters['Dimension_Reduction_Algorithm'] in ('pca'):
        training_code_used, training_labels = read_code.read_code(training_code_file,
                                                                  format=parameters['Output_Format'])
        if len(testing_data) > 0:
            testing_code_used, testing_labels = read_code.read_code(testing_code_file,
                                                                    format=parameters['Output_Format'])
        n_components = parameters['Dimension_Reduction_Number'] if parameters['Dimension_Reduction_Number'] != '' else 3
        training_code_reduced = []
        if parameters['Dimension_Reduction_Algorithm'] == 'pca':
            training_code_reduced = dimpca.pca(training_code_used, n_components=n_components)
            training_code_reduced = np.array(training_code_reduced)
            new_data = np.zeros((training_code_reduced.shape[0] , training_code_reduced.shape[1]), dtype='object')
            new_data[:, 0:] = training_code_reduced
            new_data = new_data.astype(str)
            training_code_reduced = new_data.tolist()
        save_file.save_reduction_result(training_code_reduced, file='training_code_dimension_reduction.txt')
        if n_components >= 2:
            draw_plot.plot_2d(training_code_reduced, training_labels, file='training_dimension_reduction_2d.png')
        if n_components >= 3:
            draw_plot.plot_3d(training_code_reduced, training_labels, file='training_dimension_reduction_3d.png')



    print('\n Step 2: Training and testing RF model')

    # machine learning
    ML_array = parameters['ML'].split(';')
    if parameters['ML'] != '' and parameters['Validation'] != '':
        X, y, independent = 0, 0, np.array([])                      ### X, y: training data,   independent: testing data
        X, y = read_code_ml.read_code(training_code_file, format=parameters['Output_Format'])
        classes = sorted(list(set(y)))
        if len(testing_data) > 0:
            ind_X, ind_y = read_code_ml.read_code(testing_code_file, format=parameters['Output_Format'])
            independent = np.zeros((ind_X.shape[0], ind_X.shape[1] + 1))
            independent[:, 0], independent[:, 1:] = ind_y, ind_X

        fold = parameters['Validation'] if parameters['Validation'] != '' else 5

        # RF
        n_trees = parameters['Tree_Number'] if parameters['Tree_Number'] != '' else 100

        # training model
        if len(classes) == 2:
            ML_AUCs_dict = {}
            para_info_dict = {}
            for ML in ML_array:
                para_info, cv_res, ind_res = 0, 0, 0
                if ML == 'RF':
                    para_info, cv_res, ind_res = RF.RF_Classifier(X, y, indep=independent, fold=fold, n_trees=n_trees)   ### X, y: training data,   independent: testing data
                    con = np.concatenate(tuple((f for f in cv_res)), axis=0)

                para_info_dict[ML] = para_info

                save_file.save_CV_result_binary(cv_res, '%s_CV.txt' % ML, para_info)
                ML_AUCs_dict[ML] = draw_plot.plot_roc_cv(cv_res, '%s_ROC_CV.png' % ML, label_column=0, score_column=2)
                mean_auprc = draw_plot.plot_prc_CV(cv_res, '%s_PRC_CV.png' % ML, label_column=0, score_column=2)
                cv_metrics = calculate_prediction_metrics.calculate_metrics_cv(cv_res, label_column=0, score_column=2)
                save_file.save_prediction_metrics_cv(cv_metrics, '%s_metrics_CV.txt' % ML)

                if len(testing_data) > 0:
                    save_file.save_IND_result_binary(ind_res, '%s_IND.txt' % ML, para_info)
                    ind_auc = draw_plot.plot_roc_ind(ind_res, '%s_ROC_IND.png' % ML, label_column=0, score_column=2)
                    ind_auprc = draw_plot.plot_prc_ind(ind_res, '%s_PRC_IND.png' % ML, label_column=0, score_column=2)
                    ind_metrics = calculate_prediction_metrics.calculate_metrics(ind_res[:, 0], ind_res[:, 2])
                    save_file.save_prediction_metrics_ind(ind_metrics, '%s_metrics_IND.txt' % ML)

            best_ML = ''
            best_AUC = 0
            for ml in ML_array:
                if best_AUC < ML_AUCs_dict[ml]:
                    best_AUC = ML_AUCs_dict[ml]
                    best_ML = ml
            print('\nThe model with best performance is : %s' %best_ML)
            print('The AUC value for training_data is : %s' %round(best_AUC, 2))
            print('The AUPRC value for training_data is : %s' %round(mean_auprc, 2))

            print('\nThe AUC value for testing_data is : %s' %round(ind_auc, 2))
            print('The AUPRC value for testing_data is : %s' %round(ind_auprc, 2))




