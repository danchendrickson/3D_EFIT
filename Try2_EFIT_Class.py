import numpy as np

class EFIT_Source:
    def __init__(self):
        pass

class EFIT_model:
    def __init__(self, SimSpace, gs, source, boundaries):
        pass

class EFIT_sim:
    
    def __init__(self, model, ts, time, par = False, ):
        self.maxX = model.MaxZ
        self.maxY = model.MaxZ
        self.maxZ = model.MaxZ
        pass

