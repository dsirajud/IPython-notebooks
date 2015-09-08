import sim

def testrun():
    
    w = sim.World(101)
    sim.evolve(w, 4096)
    print("done!")
