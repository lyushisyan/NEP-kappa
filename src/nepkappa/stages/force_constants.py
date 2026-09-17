"""Independent FC2, FC3, and FC4 numerical workflow stages."""

from __future__ import annotations

import shlex
import shutil

import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read, write
from hiphive import (
    ClusterSpace,
    ForceConstantPotential,
    StructureContainer,
    enforce_rotational_sum_rules,
)
from hiphive.structure_generation import generate_mc_rattled_structures
from hiphive.utilities import prepare_structures
from phonopy import Phonopy
from trainstation import Optimizer

from nepkappa.runtime import (
    ase_to_phonopy,
    phonopy_to_ase,
    progress_iter,
    run_activity_task,
)
from nepkappa.stages.base import HostStage


class HiPhiveStage(HostStage):
    """Fit and export FC2/FC3 with HiPhive."""

    name = "hiphive"

    def __init__(self, host, *, include_fc3=True):
        super().__init__(host)
        self.include_fc3 = include_fc3

    def execute(self):
        workflow = self.host
        print("\n[Step 2 - HiPhive] Generating Training Data & Fitting")
        cfg = workflow.cfg
        np.random.seed(42)
        workflow._save_phono3py_metadata()

        nx, ny, nz = cfg.dim_fc3 if self.include_fc3 else cfg.dim_fc2
        atoms_ideal = workflow.prim.repeat((nx, ny, nz))
        print(
            f"  - Generating {cfg.n_structures} rattled structures "
            f"(std={cfg.rattle_std}, min_dist={cfg.min_dist})"
        )
        structures = generate_mc_rattled_structures(
            atoms_ideal,
            cfg.n_structures,
            cfg.rattle_std,
            cfg.min_dist,
        )

        print(
            f"  - Computing forces with {workflow._calculator_name().upper()} "
            f"for {len(structures)} structures..."
        )
        for index, atoms in enumerate(
            progress_iter(
                structures,
                enabled=workflow.show_progress,
                total=len(structures),
                desc=workflow._force_desc("HiPhive"),
                unit="structure",
            ),
            1,
        ):
            forces = workflow._calculate_forces_cached(atoms, "hiphive", index)
            atoms.calc = SinglePointCalculator(atoms, forces=forces)
            atoms.get_forces()
            if index % 50 == 0:
                print(f"    Processed {index}/{len(structures)}")

        cutoffs = cfg.cutoffs if self.include_fc3 else cfg.cutoffs[:1]
        print(f"  - Fitting Force Constants (Cutoffs: {cutoffs})")
        cluster_space = ClusterSpace(workflow.prim, cutoffs)
        container = StructureContainer(cluster_space)
        for structure in prepare_structures(structures, atoms_ideal):
            container.add_structure(structure)

        optimizer = Optimizer(container.get_fit_data())
        optimizer.train()
        print(f"    RMSE: {optimizer.rmse_train:.6f}")

        parameters = enforce_rotational_sum_rules(
            cluster_space,
            optimizer.parameters,
            ["Huang", "Born-Huang"],
        )
        potential = ForceConstantPotential(cluster_space, parameters)
        potential.write(str(workflow.hiphive_model_path))
        print(f"  - HiPhive model saved to {workflow.hiphive_model_path}")

        export_label = "FC2 and FC3" if self.include_fc3 else "FC2"
        print(f"  - Exporting {export_label} from HiPhive model")
        phonopy_fc2 = Phonopy(
            ase_to_phonopy(workflow.prim),
            supercell_matrix=np.diag(cfg.dim_fc2),
        )
        fc2_supercell = phonopy_to_ase(phonopy_fc2.supercell)
        fc2_full = potential.get_force_constants(fc2_supercell).get_fc_array(2)
        workflow._write_fc2_exports(fc2_full, fc2_full=fc2_full)

        if self.include_fc3:
            phonopy_fc3 = Phonopy(
                ase_to_phonopy(workflow.prim),
                supercell_matrix=np.diag(cfg.dim_fc3),
            )
            fc3_supercell = phonopy_to_ase(phonopy_fc3.supercell)
            fc3_full = potential.get_force_constants(fc3_supercell).get_fc_array(3)
            workflow._write_fc3_exports(
                fc3_full,
                fc3_full=fc3_full,
                supercell_atoms=fc3_supercell,
                primitive_atoms=workflow.prim,
            )


class FiniteDisplacementStage(HostStage):
    """Generate FC2 and optional FC3 with phono3py finite differences."""

    name = "finite-displacement"

    def __init__(self, host, *, include_fc3=True):
        super().__init__(host)
        self.include_fc3 = include_fc3

    def execute(self):
        workflow = self.host
        print("\n[Step 2 - FiniteDisp] Phono3py Finite Displacement Method")
        cfg = workflow.cfg
        ph3 = workflow._make_phono3py()

        ph3.generate_fc2_displacements()
        ph3.save(str(workflow.disp_path))
        fc2_supercells = ph3.phonon_supercells_with_displacements
        print(f"  - Generated {len(fc2_supercells)} FC2 supercells")
        print("  - Computing forces for FC2...")
        forces_fc2 = []
        for index, supercell in enumerate(
            progress_iter(
                fc2_supercells,
                enabled=workflow.show_progress,
                total=len(fc2_supercells),
                desc=workflow._force_desc("FC2"),
                unit="structure",
            ),
            1,
        ):
            forces_fc2.append(
                workflow._calculate_forces_cached(
                    phonopy_to_ase(supercell), "fc2", index
                )
            )
        ph3.phonon_forces = np.array(forces_fc2)

        print("  - Producing force constants with phono3py finite differences")
        print("  - Producing FC2...")
        compact = getattr(cfg, "compact_fc", True)
        print(f"  - Force-constant layout: {'compact' if compact else 'full'}")
        run_activity_task(
            "Producing FC2",
            lambda: ph3.produce_fc2(is_compact_fc=compact),
            enabled=workflow.show_progress,
        )
        print("  - Symmetrizing FC2...")
        ph3.symmetrize_fc2()
        workflow._write_fc2_exports(
            ph3.fc2,
            p2s_map=ph3.phonon_primitive.p2s_map if compact else None,
            phonopy_primitive=ph3.phonon_primitive,
        )

        if not self.include_fc3:
            return None

        pair_cutoff = getattr(cfg, "pair_cutoff_fc3", None)
        if pair_cutoff is None:
            ph3.generate_displacements()
        else:
            ph3.generate_displacements(cutoff_pair_distance=pair_cutoff)
        ph3.save(str(workflow.disp_path))
        fc3_supercells = ph3.supercells_with_displacements
        print(f"  - Generated {len(fc3_supercells)} FC3 supercells")
        print("  - Computing forces for FC3...")
        forces_fc3 = []
        for index, supercell in enumerate(
            progress_iter(
                fc3_supercells,
                enabled=workflow.show_progress,
                total=len(fc3_supercells),
                desc=workflow._force_desc("FC3"),
                unit="structure",
            ),
            1,
        ):
            forces_fc3.append(
                workflow._calculate_forces_cached(
                    phonopy_to_ase(supercell), "fc3", index
                )
            )
        ph3.forces = np.array(forces_fc3)

        print("  - Producing FC3...")
        run_activity_task(
            "Producing FC3",
            lambda: ph3.produce_fc3(is_compact_fc=compact),
            enabled=workflow.show_progress,
        )
        print("  - Symmetrizing FC3...")
        ph3.symmetrize_fc3()
        export_kwargs = {}
        if workflow._fc_format() in {"shengbte", "both"}:
            export_kwargs = {
                "supercell_atoms": phonopy_to_ase(ph3.supercell),
                "primitive_atoms": phonopy_to_ase(ph3.primitive),
                "phono3py_primitive": ph3.primitive,
            }
        workflow._write_fc3_exports(
            ph3.fc3,
            fc3_nonzero_indices=ph3.fc3_nonzero_indices if compact else None,
            p2s_map=ph3.primitive.p2s_map if compact else None,
            **export_kwargs,
        )
        return None


class ThirdOrderStage(HostStage):
    """Generate FC3 through thirdorder sow, force evaluation, and reap."""

    name = "thirdorder-fc3"

    def execute(self):
        workflow = self.host
        print("\n[Step 2 - FC3] ShengBTE thirdorder finite displacements")
        workflow.fc3_root.mkdir(parents=True, exist_ok=True)
        write(
            str(workflow.fc3_root / "POSCAR"),
            workflow.prim,
            format="vasp",
            direct=True,
            vasp5=True,
        )
        print(f"  - thirdorder working directory: {workflow.fc3_root}")
        print(f"  - FC3 supercell dimension: {workflow.cfg.dim_fc3}")
        print(
            "  - FC3 cutoff: "
            f"{workflow._format_fourthorder_value(workflow.cfg.cutoff_fc3)}"
        )

        sow_command = workflow._thirdorder_args("sow")
        workflow._clear_displacement_patterns(workflow.fc3_root, "3RD.POSCAR.*")
        print(f"  - Running command: {shlex.join(sow_command)}")
        return_code = workflow._run_command(sow_command, cwd=workflow.fc3_root)
        if return_code != 0:
            raise RuntimeError(
                f"thirdorder sow failed with return code {return_code}"
            )

        patterns = workflow._thirdorder_pattern_paths()
        if not patterns:
            raise RuntimeError(
                "thirdorder sow did not generate 3RD.POSCAR.* in "
                f"{workflow.fc3_root}"
            )
        print(f"  - Generated {len(patterns)} FC3 displaced structures")
        paths = _calculate_order_forces(
            workflow,
            patterns,
            backend="thirdorder",
            label="FC3-thirdorder",
        )

        reap_command = workflow._thirdorder_args("reap")
        reap_input = "\n".join(str(path) for path in paths) + "\n"
        print(f"  - Running command: {shlex.join(reap_command)}")
        return_code = workflow._run_command_with_input(
            reap_command, reap_input, cwd=workflow.fc3_root
        )
        if return_code != 0:
            raise RuntimeError(
                f"thirdorder reap failed with return code {return_code}"
            )

        generated = workflow.fc3_root / "FORCE_CONSTANTS_3RD"
        if not generated.is_file() or generated.stat().st_size == 0:
            raise RuntimeError(
                f"thirdorder reap finished but {generated} was not found."
            )
        workflow._export_thirdorder_fc3(generated)


class FourthOrderStage(HostStage):
    """Generate FC4 through Fourthorder sow, force evaluation, and reap."""

    name = "fourthorder-fc4"

    def execute(self):
        workflow = self.host
        print("\n[Step 2 - FC4] FourPhonon Fourthorder finite displacements")
        workflow.fc4_root.mkdir(parents=True, exist_ok=True)
        write(
            str(workflow.fc4_root / "POSCAR"),
            workflow.prim,
            format="vasp",
            direct=True,
            vasp5=True,
        )
        print(f"  - Fourthorder working directory: {workflow.fc4_root}")
        print(f"  - FC4 supercell dimension: {workflow.cfg.dim_fc4}")
        print(
            "  - FC4 cutoff: "
            f"{workflow._format_fourthorder_value(workflow.cfg.cutoff_fc4)}"
        )

        sow_command = workflow._fourthorder_args("sow")
        workflow._clear_displacement_patterns(workflow.fc4_root, "4TH.POSCAR.*")
        print(f"  - Running command: {shlex.join(sow_command)}")
        return_code = workflow._run_command(sow_command, cwd=workflow.fc4_root)
        if return_code != 0:
            raise RuntimeError(
                f"Fourthorder sow failed with return code {return_code}"
            )

        patterns = workflow._fc4_pattern_paths()
        if not patterns:
            raise RuntimeError(
                "Fourthorder sow did not generate 4TH.POSCAR.* in "
                f"{workflow.fc4_root}"
            )
        print(f"  - Generated {len(patterns)} FC4 displaced structures")
        paths = _calculate_order_forces(
            workflow,
            patterns,
            backend="fourthorder",
            label="FC4",
        )

        reap_command = workflow._fourthorder_args("reap")
        reap_input = "\n".join(str(path) for path in paths) + "\n"
        print(f"  - Running command: {shlex.join(reap_command)}")
        return_code = workflow._run_command_with_input(
            reap_command, reap_input, cwd=workflow.fc4_root
        )
        if return_code != 0:
            raise RuntimeError(
                f"Fourthorder reap failed with return code {return_code}"
            )

        generated = workflow.fc4_root / "FORCE_CONSTANTS_4TH"
        if not generated.is_file():
            raise RuntimeError(
                f"Fourthorder reap finished but {generated} was not found."
            )
        if generated.resolve() != workflow.fc4_path.resolve():
            shutil.copyfile(generated, workflow.fc4_path)
        print(f"  - Generated: {workflow.fc4_path}")


def _calculate_order_forces(workflow, patterns, *, backend, label):
    """Evaluate one third-/fourth-order displacement set with cache reuse."""
    root = workflow.fc3_root if backend == "thirdorder" else workflow.fc4_root
    calculator = workflow._calculator_name()
    paths = []
    for index, pattern in enumerate(
        progress_iter(
            patterns,
            enabled=workflow.show_progress,
            total=len(patterns),
            desc=workflow._force_desc(label),
            unit="structure",
        ),
        1,
    ):
        job_dir = root / f"job-{index:05d}"
        vasprun_path = job_dir / "vasprun.xml"
        digest, payload = workflow._job_input_fingerprint(pattern, backend)
        if (
            vasprun_path.is_file()
            and vasprun_path.stat().st_size > 0
            and workflow._job_matches_inputs(job_dir, digest)
        ):
            print(f"    Skipping completed {label} #{index}: {vasprun_path}")
        else:
            atoms = read(str(pattern), format="vasp")
            if backend == "thirdorder" and calculator == "vasp":
                vasprun_path = workflow._run_vasp_thirdorder_job(
                    atoms, job_dir, index
                )
            elif backend == "fourthorder" and calculator == "vasp":
                vasprun_path = workflow._run_vasp_fc4_job(atoms, job_dir, index)
            elif backend == "thirdorder" and calculator == "nep":
                vasprun_path = workflow._run_nep_thirdorder_job(
                    atoms, job_dir, index
                )
            elif backend == "fourthorder" and calculator == "nep":
                vasprun_path = workflow._run_nep_fc4_job(atoms, job_dir, index)
            else:
                vasprun_path = workflow._run_external_order_job(atoms, job_dir)
            workflow._record_job_inputs(job_dir, digest, payload)
        paths.append(vasprun_path.resolve())
    return paths
