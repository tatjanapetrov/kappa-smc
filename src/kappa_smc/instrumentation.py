
def instrument_kappa_model(file_name, AP_triples):
    """
    INPUT: Kappa file and list of APs (triples) used in formula to model-check
           AP_triples = [(var, op, threshold_str), ...], e.g. ('I','<','2')

    OUTPUT: Kappa model string with extra observables p_0..p_{m-1}
            appended so they show up as the last columns in simulation results.
    """
    with open(file_name, 'r') as file:
        lines = file.read().splitlines()

    model = ""
    for l in lines:
        if '//' in l:
            model += l.split('//')[0] + "\n"  # remove line comments
        else:
            model += l + "\n"

    for i in range(len(AP_triples)):
        p = AP_triples[i]

        model += "%var: 'a_{}'".format(i) + " 1"
        model += " %obs: 'p_{}'".format(i) + " 'a_{}'".format(i)
        model += " %mod: " + p[0] + p[1] + str(p[2]) + " do $UPDATE " + "'a_{}'".format(i) + " 1; repeat " + p[0] + p[1] + str(p[2]) + "\n"

        # NOTE: we make AP false by updating a_i=0 on the "opposite side"
        # we assume integer thresholds
        if p[1] == '>':
            model += " %mod: " + p[0] + '<' + str(int(p[2]) + 1) + " do $UPDATE " + "'a_{}'".format(i) + " 0; repeat " + p[0] + '<' + str(int(p[2]) + 1) + "\n"
        elif p[1] == '<':
            model += " %mod: " + p[0] + '>' + str(int(p[2]) - 1) + " do $UPDATE " + "'a_{}'".format(i) + " 0; repeat " + p[0] + '>' + str(int(p[2]) - 1) + "\n"

    return model