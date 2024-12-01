from simulator import *


if __name__ == '__main__':
    sample_world = nx.Graph()
    for i in range(10):
        sample_world.add_node(i, happiness=1, status=SUSCEPTIBLE, compliance=1)

    sample_world.add_weighted_edges_from([(0, 1, 0.2), (0, 3, 0.5), (2, 5, 0.2), (1, 3, 0.1)])
    nx.set_node_attributes(sample_world, {0: EXPOSED}, 'status')
    test_sim = PandemicSimulation(world=sample_world)

    for i in range(100):
        #print(test_sim.observe())
        print(test_sim.stats)
        test_sim.step()

    test_sim.update_mandates(Q=True, M=True, SD=True)
    for i in range(100):
        #print(test_sim.observe())
        print(test_sim.stats)
        test_sim.step()

    test_sim.update_mandates(Q=False, M=False, SD=False)
    for i in range(500):
        #print(test_sim.observe())
        print(test_sim.stats)
        test_sim.step()
