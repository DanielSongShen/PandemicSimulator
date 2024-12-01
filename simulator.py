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


class HappinessModel:
    def __init__(self, q=0.1, m=0.05, sd=0.01, reg=0.05):
        self.h_q = q  # Reduce happiness by h_q
        self.h_m = m  # Reduce happiness by h_m
        self.h_sd = sd  # Reduce happiness by h_sd
        self.h_reg = reg  # Regeneration rate of happiness

    def __call__(self, *args, **kwargs):
        # TODO test negative exponential happiness functions
        h = args[0] + self.h_reg
        mandates = args[1]
        if mandates["Q"]:
            h -= self.h_q
        if mandates["M"]:
            h -= self.h_m
        if mandates["SD"]:
            h -= self.h_sd
        return min(max(h, 0), 1)


class ComplianceModel:
    def __init__(self, influence=0.3, sensitivity=0.5):
        self.alpha = influence  # Influence of happiness on compliance
        self.k = sensitivity  # Sensitivity of compliance as a factor of happiness

    def __call__(self, *args, **kwargs):
        h = kwargs['h']
        p = kwargs['p']
        return p + self.alpha*pow(h, self.k)*(1.-p)


class PandemicSimulation:
    def __init__(self,
                 transmission_model=SEIRSTransmissionModel(),
                 error_model=ErrorModel(),
                 happiness_model=HappinessModel(),
                 compliance_model=ComplianceModel(),
                 birth_rate=0.0182,
                 death_rate=0.009841,
                 r_q=0.01,
                 r_m=0.15,
                 r_sd=0.3):
        self.world = nx.Graph()
        self.t = 0
        self.br = birth_rate
        self.dr = death_rate
        self.stats = {"N": 0, SUSCEPTIBLE: 0, EXPOSED: 0, INFECTIOUS: 0, RECOVERED: 0, DECEASED: 0, "H": 0}
        self.mandates = {"Q": False, "M": False, "SD": False, "CT": False}
        self.mandate_effects = {"Q": r_q, "M": r_m, "SD": r_sd}
        self.transmission_model = transmission_model
        self.error_model = error_model
        self.happiness_model = happiness_model
        self.compliance_model = compliance_model

    def apply_mandates(self, compliance, interactions):
        updated_interactions = interactions
        if self.mandates["Q"]:
            if random.random < compliance:
                updated_interactions = [_*self.mandate_effects["Q"] for _ in updated_interactions]
        for i in range(len(updated_interactions)):
            if self.mandates["M"]:
                if random.random < compliance:
                    updated_interactions[i] = updated_interactions[i]*self.mandate_effects["M"]
            if self.mandates["SD"]:
                if random.random < compliance:
                    updated_interactions[i] = updated_interactions[i]*self.mandate_effects["SD"]
        return updated_interactions

    def step(self):
        # TODO add birth and death rates
        for node in self.world.nodes:
            # First do behavior model
            c = self.compliance_model(c=node["compliance"], h=node["happiness"])
            interactions = self.apply_mandates(c, self.world.edges[node])
            # Transmit disease
            status = node["status"]
            self.stats[status] -= 1
            node["status"] = self.transmission_model.step(status,
                                                          interactions=interactions,
                                                          neighbors=self.world[node])
            self.stats[node["status"]] += 1
            node["happiness"] = self.happiness_model(node["happiness"], self.mandates)
        self.t += 1

    def observe(self):
        result = {}
        for stat, val in self.stats.items():
            result[stat] = val + self.error_model(self.t)
        # TODO add contact tracing
        return result
