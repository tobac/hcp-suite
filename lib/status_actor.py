import ray

@ray.remote
class StatusActor:
    def __init__(self):
        self.status_dict = {}

    def send(self, status):
        pid = status[0]
        node = status[1]
        msg = status[2]
        self.status_dict[pid] = {"msg": msg, "node": node}

    def exit(self):
        #self.event.set()
        ray.actor.exit_actor()

    def get_status(self):
        return self.status_dict

@ray.remote
class ResultsActor:
    def __init__(self):
        # Create dictionaries to keep results (it makes sense to do this class-wide to add results on-the-fly
        # and for later reference if get results functions are called too early for example
        self.fselection_results = {}
        self.fselection_results[-1] = {} # Create sub-dictionary for original (i.e. non-permuted) data
        self.prediction_results = {}

    def send(self, results):
        """
        Receives result from worker, determines type of result and calls the appropriate
        processing function
        """
        if type(result) == 'list':
            process_fselection_results(results)
        elif type(result) == 'dict':
            process_prediction_results(results)

    def process_fselection_results(self, results):
        n = 1
        N = len(results)
        printv("\n")
        for result in results:
            fold = result[0]
            perm = result[1]
            df = result[2]
            printv("Rearranging result {} of {}".format(n, N), update=True)
            self.fselection_results[perm][fold] = df
            n += 1

    def process_prediction_results(self, results):
        for results_dict in results:
            if results_dict['perm'] not in self.prediction_results:
                self.prediction_results[results_dict['perm']] = pd.DataFrame()
                self.prediction_results[results_dict['perm']]['observed'] = self.data_dict['data'][self.data_dict['behav']]
            for tail in ('pos', 'neg', 'glm'):
                self.prediction_results[results_dict['perm']].loc[results_dict['test_IDs'], [tail]] = results_dict[tail]

