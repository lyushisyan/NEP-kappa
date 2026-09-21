"""Machine-readable scope of the implemented temperature-renormalization route."""


def sscha_approximation(*, qha_volume=False):
    """Describe calculated quantities, not a claim of numerical convergence."""
    return {
        "label": (
            "SSCHA at QHA equilibrium volumes"
            if qha_volume else "Fixed-cell Phonopy SSCHA"
        ),
        "volume_treatment": (
            "QHA equilibrium volume; no feedback from SSCHA free energy"
            if qha_volume else "Fixed input cell at each temperature"
        ),
        "fc2_role": "Auxiliary harmonic FC2 from stochastic self-consistent fitting",
        "free_energy_hessian_calculated": False,
        "dynamic_bubble_calculated": False,
        "sscha_cell_optimization": False,
        "results_added_together": False,
        "limitations": [
            "Auxiliary frequencies are not dynamic self-energy-corrected spectral peaks.",
            "FC3/FC4 transport does not automatically add a bubble frequency shift.",
            "Completed iterations do not establish sampling or finite-size convergence.",
        ],
    }
