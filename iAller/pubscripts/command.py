protein_cmd_coding = {
    'AAC': ['AAC.AAC(training_data, **kw)', 'AAC.AAC(testing_data, **kw)'],
    'CKSAAP': ['CKSAAP.CKSAAP(training_data, gap=%d, **kw)' % int(parameters['K_Space']), 'CKSAAP.CKSAAP(testing_data, gap=%d, **kw)' % int(parameters['K_Space'])],
    'DPC': ['DPC.DPC(training_data, **kw)', 'DPC.DPC(testing_data, **kw)'],
    }