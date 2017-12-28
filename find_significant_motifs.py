def calc_model_log_likelihood(probas, labels):
    log_likelihood = 0
    Y = labels.astype(float)
    for i in range(len(Y)):
        p_true = probas[i][1]
        p_false = probas[i][0]
        y = Y[i]
        prod = ((p_true)**y) * ((p_false) ** (1-y))
        log_prod = np.log(prod)
        log_likelihood += log_prod
    return log_likelihood
    
def calc_feature_pvals(features,
                            labels,
                            test_size=0.5):
    num_motifs = features.shape[1]
    # split data into training and test sets
    training_features, test_features, training_labels, test_labels = get_split(
        features, labels, test_size = test_size)

    #  Train affinity classifier
    classifier = sklearn.linear_model.LogisticRegression(penalty='l1', n_jobs=-1)

    classifier.fit(training_features, training_labels)
    # score predictions

    probas = classifier.predict_proba(features) # [[p_false, p_true]...] 
    overall_log_likelihood = calc_model_log_likelihood(probas, labels)
    pvals = []
    for motif in features.columns.values:
        current_features = features.drop(motif, axis=1, inplace=False)
        current_training_features = training_features.drop(motif, axis=1, inplace = False)
        current_classifier = sklearn.linear_model.LogisticRegression(penalty='l1', n_jobs=-1)
        current_classifier.fit(current_training_features, training_labels)
        current_probas = current_classifier.predict_proba(current_features)
        current_log_likelihood = calc_model_log_likelihood(current_probas, labels)

        stat = -2*(current_log_likelihood - overall_log_likelihood)
        p = scipy.stats.chisqprob(stat, df=1)

        pvals.append(p)

    return pvals



if __name__ == '__main__':

num_iterations = 10

motifs = motif_score_frame.columns.values[3:]
monomer_lr_test_pvals_dict = {}
for treatment in ['veh', 'kla']:
    if treatment == 'veh':
        members = ['atf3', 'cjun', 'jund']
    else:
        members = ['atf3', 'cjun', 'fos', 'fra2', 'junb', 'jund']
    for monomer in members:
        monomer_lr_test_pvals_dict[monomer + '_' + treatment] = []
        
        for i in range(num_iterations):
            start = time.time()
            features = monomer_treatment_standardized_features_dict[monomer + '_' + treatment]
            labels = pd.Series([False if 'background' in x else True for x in features.index.values],
                               index = features.index.values)
            pvals = calc_feature_pvals(features,
                               labels,
                               test_size=0.2)
            monomer_lr_test_pvals_dict[monomer + '_' + treatment].append(pvals)
            print(monomer, treatment, time.time() - start)
    
            
