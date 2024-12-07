from simulator import *
from utils import *
import matplotlib.pyplot as plt


def initialize_world(N=1000, p_edge=0.1, exposed_pop=0.05):
    world = nx.Graph()
    for i in range(N):
        h = max(min(stats.norm.rvs(loc=0.5, scale=0.2), 1), 0)
        c = max(min(stats.norm.rvs(loc=0.5, scale=0.2), 1), 0)
        world.add_node(i, happiness=h, status=SUSCEPTIBLE, compliance=c)

    for node in world.nodes:
        connected = [to for (fr, to) in world.edges(node)]
        unconnected = [n for n in world.nodes() if n not in connected and n != node]
        new_edges = []
        for _ in unconnected:
            if random.random() < p_edge:
                w = max(min(stats.norm.rvs(loc=0.5, scale=0.2), 1), 0)
                new_edges.append((node, _, w))
        world.add_weighted_edges_from(new_edges)

    # Initialize subset of population as exposed
    inds = list(range(N))
    random.shuffle(inds)
    exposed = inds[:int(N*exposed_pop)]
    new_attributes = {}
    for node in exposed:
        new_attributes[node] = EXPOSED
    print(new_attributes)
    nx.set_node_attributes(world, new_attributes, 'status')
    return world


def model(o):
    return random.choice(list(range(16)))


if __name__ == '__main__':
    iters = 1000  # Number of days to simulate
    density = 0.1  # Density of connections
    init_exposed = 0.05  # Percentage of total pop that was exposed
    N = 100
    world = initialize_world(N=N, p_edge=density, exposed_pop=init_exposed)
    #nx.draw(world)
    #plt.show()
    world_sim = PandemicSimulation(world=world)
    mandates_list = list(world_sim.mandates.keys())
    print(mandates_list)
    for day in range(iters):
        mandates = world_sim.mandates
        print("Start of day", day)
        print("Mandates", mandates)
        observation = world_sim.observe()
        print("Observations")
        world_sim.print_observations()
        print("Stats")
        print(world_sim.stats)
        # Get model action
        action = model(observation)
        #print("action", action)
        # Parse action
        new_mandates = {}
        base2 = numberToBase(action, 2)
        while len(base2) < len(mandates_list):
            base2.insert(0, 0)
        print("action", base2)
        for _ in range(len(base2)):
            if base2[_]:
                new_mandates[mandates_list[_]] = True
            else:
                new_mandates[mandates_list[_]] = False
        world_sim.update_mandates(new_mandates)
        print()

    # sample_world = nx.Graph()
    # for i in range(10):
    #     sample_world.add_node(i, happiness=1, status=SUSCEPTIBLE, compliance=1)
    #
    # sample_world.add_weighted_edges_from([(0, 1, 0.2), (0, 3, 0.5), (2, 5, 0.2), (1, 3, 0.1)])
    # nx.set_node_attributes(sample_world, {0: EXPOSED}, 'status')
    # test_sim = PandemicSimulation(world=sample_world)
    #
    # for i in range(100):
    #     #print(test_sim.observe())
    #     print(test_sim.stats)
    #     test_sim.step()
    #
    # test_sim.update_mandates(Q=True, M=True, SD=True)
    # for i in range(100):
    #     #print(test_sim.observe())
    #     print(test_sim.stats)
    #     test_sim.step()
    #
    # test_sim.update_mandates(Q=False, M=False, SD=False)
    # for i in range(500):
    #     #print(test_sim.observe())
    #     print(test_sim.stats)
    #     test_sim.step()
