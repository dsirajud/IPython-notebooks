def trans(W):
    W = np.transpose(W)
    return W

def trans2(W):
    W.T
    return W

class test():
    def __init__(self):
        self.array = np.arange(10).reshape(2,5)*20

def transinst(inst):
    inst.array = np.transpose(inst.array)
    return inst

def transinst2(inst):
    inst.array = np.transpose(inst.array)

