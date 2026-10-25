"""Finite route computations for Repair Filters Mutational Histories by Event Count.

All normalization is conditional on arriving at a given k-hit damage class.
This is an auxiliary equilibrium-contribution law (and a leading genealogy
law only under the paper's specified branching-immigration realization).
"""
from .model import (route_partition, event_distribution, survival_factor,
                    two_hit_stats, lambda_i, route_activity)
__all__ = ['route_partition', 'event_distribution', 'survival_factor',
           'two_hit_stats', 'lambda_i', 'route_activity']
