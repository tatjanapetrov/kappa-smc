
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

    for i, (var, op, threshold) in enumerate(AP_triples):
        condition = f"{var} {op} {threshold}"
        model += (f"%obs: 'p_{i}' ({condition}) [?] 1 [:] 0\n")

    return model