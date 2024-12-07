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
        return status


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


class Reward:
    def __init__(self, alpha=0.5, beta=0.5, w0=1, w1=1000):
        self.alpha = alpha
        self.beta = beta
        self.w0 = w0
        self.w1 = w1

    def __call__(self, *args, **kwargs):
        s = kwargs['stats']
        h = s['H']
        n_e = s[EXPOSED]
        n_i = s[INFECTIOUS]
        n_d = s[DECEASED]
        return self.alpha*h-self.beta*(self.w0*(n_i+n_e)+self.w1*n_d)


class PandemicSimulation:
    def __init__(self,
                 world=nx.Graph(),
                 transmission_model=SEIRSTransmissionModel(),
                 error_model=ErrorModel(),
                 happiness_model=HappinessModel(),
                 compliance_model=ComplianceModel(),
                 birth_rate=0.0182,  # From United Nations
                 death_rate=0.009841,  # From United Nations
                 r_q=0.01,
                 r_m=0.15,
                 r_sd=0.3):
        self.world = world
        self.statuses = nx.get_node_attributes(self.world, "status")
        self.t = 0
        self.br = birth_rate
        self.dr = death_rate
        self.stats = {"N": len(world.nodes),
                      SUSCEPTIBLE: 0,
                      EXPOSED: 0,
                      INFECTIOUS: 0,
                      RECOVERED: 0,
                      DECEASED: 0,
                      "H": sum([world.nodes[_]["happiness"] for _ in world.nodes])
                      }
        for node in self.world.nodes:
            self.stats[self.statuses[node]] += 1
        self.mandates = {"Q": False, "M": False, "SD": False, "CT": False}
        self.mandate_effects = {"Q": r_q, "M": r_m, "SD": r_sd}
        self.transmission_model = transmission_model
        self.error_model = error_model
        self.happiness_model = happiness_model
        self.compliance_model = compliance_model

    def apply_mandates(self, compliance, interactions):
        updated_interactions = interactions
        if self.mandates["Q"]:
            if random.random() < compliance:
                updated_interactions = [_*self.mandate_effects["Q"] for _ in updated_interactions]
        for i in range(len(updated_interactions)):
            if self.mandates["M"]:
                if random.random() < compliance:
                    updated_interactions[i] = updated_interactions[i]*self.mandate_effects["M"]
            if self.mandates["SD"]:
                if random.random() < compliance:
                    updated_interactions[i] = updated_interactions[i]*self.mandate_effects["SD"]
        return updated_interactions

    def step(self):
        self.statuses = nx.get_node_attributes(self.world, "status")
        # TODO add birth and death rates
        for node in self.world.nodes:
            node_dict = self.world.nodes[node]
            # First do behavior model
            c = self.compliance_model(p=node_dict["compliance"], h=node_dict["happiness"])
            neighbors = self.world[node]
            default_interactions = [self.world.get_edge_data(node, neighbor, default=0)['weight'] for neighbor in neighbors]
            interactions = self.apply_mandates(c, default_interactions)
            # Transmit disease
            status = self.statuses[node]
            self.stats[status] -= 1
            new_status = self.transmission_model.step(status,
                                                      interactions=interactions,
                                                      neighbors=[self.statuses[_] for _ in neighbors])
            nx.set_node_attributes(self.world, {node: new_status}, 'status')
            self.stats[new_status] += 1
            self.stats["H"] -= node_dict["happiness"]
            if new_status != DECEASED:
                new_happiness = self.happiness_model(node_dict["happiness"], self.mandates)
            else:
                new_happiness = 0
            nx.set_node_attributes(self.world, {node: new_happiness}, 'happiness')
            self.stats["H"] += new_happiness

        self.t += 1
        assert self.stats["N"] == sum([val for key, val in self.stats.items() if key != "H" and key != "N"])

    def observe(self):
        result = {}
        for stat, val in self.stats.items():
            if stat != "N":
                error = self.error_model(self.t)
            else:
                error = 0
            result[stat] = val + error
        # TODO add contact tracing
        return result

    def print_observations(self):
        obs = self.observe()
        for key, val in obs.items():
            if key == 0:
                print("SUSCEPTIBLE: %.2f" % val)
            elif key == 1:
                print("EXPOSED: %.2f" % val)
            elif key == 2:
                print("INFECTIOUS: %.2f" % val)
            elif key == 3:
                print("RECOVERED: %.2f" % val)
            elif key == 4:
                print("DECEASED: %.2f" % val)
            else:
                print(f"{key}: {val: .2f}")

    def update_mandates(self, new_mandates):
        for key, val in new_mandates.items():
            self.mandates[key] = val
