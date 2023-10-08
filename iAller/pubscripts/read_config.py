def read_config(file):
    parameters = {}
    with open(file) as f:
        records = f.readlines()
    for line in records:
        if line[0] == '#' or line.strip() == '':
            continue
        else:
            # print(line.strip())
            array = line.strip().split('=')
            parameters[array[0]] = array[1]

    for key in parameters:
        if parameters[key].isdigit():
            parameters[key] = int(parameters[key])

    if parameters['Selected_Feature_Number'] != '':
        parameters['Selected_Feature_Number'] = int(parameters['Selected_Feature_Number'])



    default = {
        'Sequence_Type': 'Protein',
        'Method':'AAC;DPC;CKSAAP',
        'Kmer_Size': 3,
        'K_space': 5,
        'Output_Format': 'svm',
        'Clustering_Algorithm': '',
        'Kmean_Cluster_Number': 2,
        'Clustering_Type': 'sample',
        'Feature_Normalization_Algorithm': '',
        'Feature_Selection_Algorithm': 'pearsonr',
        'Selected_Feature_Number': 200,
        'Dimension_Reduction_Algorithm': 'pca',
        'Dimension_Reduction_Number': 3,
        'ML': 'RF',
        'Kernel': 'rbf',
        'Cost': 1.0,
        'Gamma': '',
        'Auto_Opterimize': True,
        'Tree_Number': 100,
        'Validation': 5,
        'Ensemble': 'Yes'
    }

    for key in default:
        if key in parameters:
            if parameters[key] == '':
                parameters[key] = default[key]
    return parameters
