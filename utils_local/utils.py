import pickle


class BatteryFig:
    def __init__(self, fig=None, path=None):
        self.fig = fig
        self.path = path

    def save_fig(self):
        pickle.dump(obj=self.fig,
                    file=open(self.path, "wb"))

    def load_fig(self):
        fig = pickle.load(
            file=open(self.path, "rb"))
        return fig

