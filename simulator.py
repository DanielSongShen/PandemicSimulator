import networkx as nx
import random
from scipy import stats
import numpy as np


SUSCEPTIBLE = 0
EXPOSED = 1
INFECTIOUS = 2
RECOVERED = 3
DECEASED = 4

class SEIRSTransmissionModel:
    def __init__(self, trans_rate=0.1, inc_per=5, inf_per=10, fate_rate=0.0001, prot_per=14):
        self.beta = trans_rate  # Transmission rate (susceptible -> exposed)
        self.sigma = 1./inc_per  # Infectious rate (exposed -> infectious)
        self.gamma = 1./inf_per  # Recovery rate (infectious -> recovered)
        self.delta = fate_rate  # Fatality rate (infectious -> deceased)
        self.omega = 1./prot_per  # Loss of immunity rate (recovered -> susceptible)

    def step(self, status, interactions, neighbors):
        """
        :param status: Status of single node representing a person in the world
        :param interactions: List of edge weights representing interaction rate with each neighbor
        :param neighbors: List of nodes representing neighbors
        :return: Updated node state
        """
        if status == SUSCEPTIBLE:
            for i in range(len(neighbors)):
                other_status = neighbors[i]
                interaction_rate = interactions[i]
                if other_status == INFECTIOUS and random.random() < self.beta*interaction_rate:
                    return EXPOSED
        if status == EXPOSED:
            if random.random() < self.sigma:
                return INFECTIOUS
        if status == INFECTIOUS:
            if random.random() < self.gamma:
                return RECOVERED
            if random.random() < self.delta:
                return DECEASED
        if status == RECOVERED:
            if random.random() < self.omega:
                return SUSCEPTIBLE

class ErrorModel:
    def __init__(self, error_mag=0.1, error_decay=0.9, std_dev=0.1):
        self.alpha = error_mag  # Magnitude of error
        self.beta = error_decay  # Decay rate of error
        self.epsilon = stats.norm(loc=0, scale=std_dev)

    def __call__(self, *args, **kwargs):
        t = args[0]
        return self.alpha*self.epsilon.rvs()*np.exp(-self.beta*t)


class PandemicSimulation:
    def __init__(self, transmission_model=SEIRSTransmissionModel(), error_model=ErrorModel()):
        self.world = nx.Graph()
        self.t = 0
        self.stats = {"N": 0, "S": 0, "E": 0, "I": 0, "D": 0, "H": 0}
        self.mandates = {"Q": False, "M": False, "SD": False, "CT": False}
        self.transmission_model = transmission_model
        self.error_model = error_model

    def step(self):
        asdf

    def observe(self):
        result = {}
        for stat, val in self.stats.items():
            result[stat] = val + self.error_model(self.t)
        # TODO add contact tracing
        return result
