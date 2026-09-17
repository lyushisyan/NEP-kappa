"""Typed, read-only views over the legacy flat workflow configuration.

The views let new code depend on cohesive configuration sections while the
existing workflows continue to use their historical ``config.option`` API.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from typing import Any, Mapping, Optional, Tuple


Int3 = Tuple[int, int, int]
FloatValues = Tuple[float, ...]
StringValues = Tuple[str, ...]


def _tuple(value):
    return None if value is None else tuple(value)


def _mapping(value):
    return dict(value or {})


@dataclass(frozen=True)
class StructureConfig:
    poscar: str
    dimensionality: int
    effective_thickness: Optional[float]
    effective_area: Optional[float]
    vacuum_axis: str
    periodic_axis: str


@dataclass(frozen=True)
class WorkflowPlanConfig:
    preset: str
    steps: StringValues


@dataclass(frozen=True)
class CalculatorConfig:
    name: str
    nep_model: Optional[str]
    factory: Optional[str]
    kwargs: Mapping[str, Any]
    model_files: StringValues
    mace_model: Optional[str]
    mace_foundation: Optional[str]
    mace_device: Optional[str]
    mace_dtype: Optional[str]
    vasp_command: Optional[str]
    vasp_path: Optional[str]
    potcar_path: Optional[str]
    vasp_kwargs: Mapping[str, Any]


@dataclass(frozen=True)
class RelaxationConfig:
    enabled: bool
    workdir: str
    stages: Optional[Mapping[str, Any]]


@dataclass(frozen=True)
class ForceConstantsConfig:
    dim_fc2: Int3
    dim_fc3: Int3
    dim_fc4: Int3
    fc3_backend: str
    cutoff_fc3: float
    pair_cutoff_fc3: Optional[float]
    cutoff_fc4: float
    thirdorder_command: Optional[str]
    fourthorder_command: Optional[str]
    fc3_workdir: str
    fc4_workdir: str
    vasp_workdir: str
    use_hiphive: bool
    compact: bool
    format: str
    n_structures: int
    rattle_std: float
    cutoffs: FloatValues
    min_dist: float
    parallel: Mapping[str, Any]


@dataclass(frozen=True)
class KappaConfig:
    mesh: Int3
    temperatures: FloatValues
    method: str
    command: Optional[str]
    isotope: bool
    boundary_mfp: float
    wigner: bool
    parallel: Mapping[str, Any]


@dataclass(frozen=True)
class FourPhononConfig:
    enabled: bool
    command: str
    workdir: str
    control: Optional[str]
    fc2: Optional[str]
    fc3: Optional[str]
    fc4: Optional[str]
    harmonic_format: str
    mesh: Optional[Int3]
    temperatures: Optional[FloatValues]
    supercell: Optional[Int3]
    solver: str
    scalebroad: float
    isotopes: bool
    nonanalytic: bool
    only_harmonic: bool
    sample_3ph: int
    sample_3ph_phase_space: int
    sample_4ph: int
    sample_4ph_phase_space: int
    mpi_launcher: StringValues
    mpi_processes: int
    omp_threads: int
    omp_stacksize: str
    parallel: Mapping[str, Any]


@dataclass(frozen=True)
class PlotConfig:
    layout: str
    path: str
    path_points: Any
    path_segments: Any
    tau: str
    kappa: str
    temperature: float
    dpi: int


@dataclass(frozen=True)
class QHAConfig:
    enabled: bool
    volume_ratios: FloatValues
    dim_fc2: Optional[Int3]
    mesh: Int3
    temperatures: FloatValues
    eos: str
    pressure: float
    relax_internal: bool
    relax_fmax: float
    relax_steps: int
    displacement_distance: float
    imaginary_frequency_tolerance: float
    cutoff_frequency: float
    vasp_relax_kwargs: Mapping[str, Any]
    vasp_static_kwargs: Mapping[str, Any]


@dataclass(frozen=True)
class SCPHConfig:
    enabled: bool
    workdir: str
    temperatures: FloatValues
    initial_fc2: Optional[str]
    born: Optional[str]
    snapshots: int
    iterations: int
    transient: int
    sscha_mesh: Int3
    random_seed: int
    cutoff_frequency: float
    fc_calculator: str
    fc_calculator_options: Optional[str]
    save_datasets: bool
    run_transport: bool
    transport_fc3: Optional[str]
    transport_metadata: Optional[str]


@dataclass(frozen=True)
class QHASSCHAConfig:
    enabled: bool
    three_phonon: Optional[bool]
    four_phonon: bool


@dataclass(frozen=True)
class OutputConfig:
    result_dir: str
    progress: bool


@dataclass(frozen=True)
class WorkflowSections:
    workflow: WorkflowPlanConfig
    structure: StructureConfig
    calculator: CalculatorConfig
    relaxation: RelaxationConfig
    force_constants: ForceConstantsConfig
    kappa: KappaConfig
    fourphonon: FourPhononConfig
    plot: PlotConfig
    qha: QHAConfig
    scph: SCPHConfig
    qha_sscha: QHASSCHAConfig
    output: OutputConfig

    @classmethod
    def from_flat(cls, cfg):
        """Build typed immutable sections from a validated flat config."""
        return cls(
            workflow=WorkflowPlanConfig(
                preset=cfg.workflow_preset,
                steps=_tuple(cfg.workflow_steps),
            ),
            structure=StructureConfig(
                poscar=cfg.poscar,
                dimensionality=cfg.dimensionality,
                effective_thickness=cfg.effective_thickness,
                effective_area=cfg.effective_area,
                vacuum_axis=cfg.vacuum_axis,
                periodic_axis=cfg.periodic_axis,
            ),
            calculator=CalculatorConfig(
                name=cfg.calculator,
                nep_model=cfg.nep_model,
                factory=cfg.calculator_factory,
                kwargs=_mapping(cfg.calculator_kwargs),
                model_files=_tuple(cfg.calculator_model_files),
                mace_model=cfg.mace_model,
                mace_foundation=cfg.mace_foundation,
                mace_device=cfg.mace_device,
                mace_dtype=cfg.mace_dtype,
                vasp_command=cfg.vasp_command,
                vasp_path=cfg.vasp_path,
                potcar_path=cfg.potcar_path,
                vasp_kwargs=_mapping(cfg.vasp_kwargs),
            ),
            relaxation=RelaxationConfig(
                enabled=cfg.do_relax,
                workdir=cfg.vasp_relax_workdir,
                stages=(
                    None
                    if cfg.vasp_relax_stages is None
                    else _mapping(cfg.vasp_relax_stages)
                ),
            ),
            force_constants=ForceConstantsConfig(
                dim_fc2=_tuple(cfg.dim_fc2),
                dim_fc3=_tuple(cfg.dim_fc3),
                dim_fc4=_tuple(cfg.dim_fc4),
                fc3_backend=cfg.fc3_backend,
                cutoff_fc3=cfg.cutoff_fc3,
                pair_cutoff_fc3=cfg.pair_cutoff_fc3,
                cutoff_fc4=cfg.cutoff_fc4,
                thirdorder_command=cfg.thirdorder_command,
                fourthorder_command=cfg.fourthorder_command,
                fc3_workdir=cfg.fc3_workdir,
                fc4_workdir=cfg.fc4_workdir,
                vasp_workdir=cfg.vasp_workdir,
                use_hiphive=cfg.use_hiphive,
                compact=cfg.compact_fc,
                format=cfg.fc_format,
                n_structures=cfg.n_structures,
                rattle_std=cfg.rattle_std,
                cutoffs=_tuple(cfg.cutoffs),
                min_dist=cfg.min_dist,
                parallel=_mapping(cfg.force_parallel),
            ),
            kappa=KappaConfig(
                mesh=_tuple(cfg.mesh),
                temperatures=_tuple(cfg.temps),
                method=cfg.method,
                command=cfg.kappa_command,
                isotope=cfg.isotope,
                boundary_mfp=cfg.bfmp,
                wigner=cfg.wigner,
                parallel=_mapping(cfg.lbte_parallel),
            ),
            fourphonon=FourPhononConfig(
                enabled=cfg.fp_enabled,
                command=cfg.fp_command,
                workdir=cfg.fp_workdir,
                control=cfg.fp_control,
                fc2=cfg.fp_fc2,
                fc3=cfg.fp_fc3,
                fc4=cfg.fp_fc4,
                harmonic_format=cfg.fp_harmonic_format,
                mesh=_tuple(cfg.fp_mesh),
                temperatures=_tuple(cfg.fp_temps),
                supercell=_tuple(cfg.fp_scell),
                solver=cfg.fp_solver,
                scalebroad=cfg.fp_scalebroad,
                isotopes=cfg.fp_isotopes,
                nonanalytic=cfg.fp_nonanalytic,
                only_harmonic=cfg.fp_only_harmonic,
                sample_3ph=cfg.fp_sample_3ph,
                sample_3ph_phase_space=cfg.fp_sample_3ph_phase_space,
                sample_4ph=cfg.fp_sample_4ph,
                sample_4ph_phase_space=cfg.fp_sample_4ph_phase_space,
                mpi_launcher=_tuple(cfg.fp_mpi_launcher),
                mpi_processes=cfg.fp_mpi_processes,
                omp_threads=cfg.fp_omp_threads,
                omp_stacksize=cfg.fp_omp_stacksize,
                parallel=_mapping(cfg.fp_parallel),
            ),
            plot=PlotConfig(
                layout=cfg.plot_layout,
                path=cfg.plot_path,
                path_points=cfg.plot_path_points,
                path_segments=cfg.plot_path_segments,
                tau=cfg.plot_tau,
                kappa=cfg.plot_kappa,
                temperature=cfg.plot_temperature,
                dpi=cfg.plot_dpi,
            ),
            qha=QHAConfig(
                enabled=cfg.qha_enabled,
                volume_ratios=_tuple(cfg.qha_volume_ratios),
                dim_fc2=_tuple(cfg.qha_dim_fc2),
                mesh=_tuple(cfg.qha_mesh),
                temperatures=_tuple(cfg.qha_temps),
                eos=cfg.qha_eos,
                pressure=cfg.qha_pressure,
                relax_internal=cfg.qha_relax_internal,
                relax_fmax=cfg.qha_relax_fmax,
                relax_steps=cfg.qha_relax_steps,
                displacement_distance=cfg.qha_displacement_distance,
                imaginary_frequency_tolerance=cfg.qha_imaginary_frequency_tolerance,
                cutoff_frequency=cfg.qha_cutoff_frequency,
                vasp_relax_kwargs=_mapping(cfg.qha_vasp_relax_kwargs),
                vasp_static_kwargs=_mapping(cfg.qha_vasp_static_kwargs),
            ),
            scph=SCPHConfig(
                enabled=cfg.scph_enabled,
                workdir=cfg.scph_workdir,
                temperatures=_tuple(cfg.scph_temps),
                initial_fc2=cfg.scph_initial_fc2,
                born=cfg.scph_born,
                snapshots=cfg.scph_snapshots,
                iterations=cfg.scph_iterations,
                transient=cfg.scph_transient,
                sscha_mesh=_tuple(cfg.scph_sscha_mesh),
                random_seed=cfg.scph_random_seed,
                cutoff_frequency=cfg.scph_cutoff_frequency,
                fc_calculator=cfg.scph_fc_calculator,
                fc_calculator_options=cfg.scph_fc_calculator_options,
                save_datasets=cfg.scph_save_datasets,
                run_transport=cfg.scph_run_transport,
                transport_fc3=cfg.scph_transport_fc3,
                transport_metadata=cfg.scph_transport_metadata,
            ),
            qha_sscha=QHASSCHAConfig(
                enabled=cfg.qha_sscha_enabled,
                three_phonon=cfg.qha_sscha_three_phonon,
                four_phonon=cfg.qha_sscha_four_phonon,
            ),
            output=OutputConfig(result_dir=cfg.result_dir, progress=cfg.progress),
        )


class WorkflowConfig(argparse.Namespace):
    """Backward-compatible flat config with a typed ``sections`` view."""

    @property
    def sections(self) -> WorkflowSections:
        return WorkflowSections.from_flat(self)
