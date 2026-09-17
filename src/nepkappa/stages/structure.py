"""Input-structure loading and optional relaxation stage."""

from __future__ import annotations

from ase.io import read, write

from nepkappa.stages.base import HostStage


class StructureRelaxationStage(HostStage):
    """Load the primitive structure and apply the configured relaxation backend."""

    name = "relax-structure"

    def __init__(self, host, *, ase_relax=None):
        super().__init__(host)
        self.ase_relax = ase_relax

    def validate(self):
        if not getattr(self.host.cfg, "poscar", None):
            raise ValueError("structure.poscar is required for relaxation")

    def execute(self):
        workflow = self.host
        print("\n[Step 1] Relax Structure")
        print(f"  - Reading from {workflow.cfg.poscar}")
        primitive = read(workflow.cfg.poscar)

        if not workflow.cfg.do_relax:
            write(
                str(workflow.relaxed_poscar_path),
                primitive,
                format="vasp",
                direct=True,
                vasp5=True,
            )
            print("  - Relaxation disabled; copied input structure for downstream steps")
            print(f"  - Structure saved to {workflow.relaxed_poscar_path}")
            workflow.prim = primitive
            return primitive

        calculator = workflow._calculator_name()
        if calculator == "vasp":
            primitive = workflow._run_vasp_relaxation(primitive)
        elif calculator == "nep":
            print(f"  - Relaxing structure using NEP: {workflow.cfg.nep_model}")
            primitive.calc = workflow._make_nep_calculator()
            primitive = self._relax_ase(primitive)
        else:
            print(f"  - Relaxing structure using {calculator.upper()} ASE calculator")
            primitive.calc = workflow._make_external_backend().calculator()
            primitive = self._relax_ase(primitive)
        workflow.prim = primitive
        return primitive

    def _relax_ase(self, atoms):
        relax = self.ase_relax
        if relax is None:
            from calorine.tools import relax_structure

            relax = relax_structure
        relax(atoms, fmax=1.0e-3)
        relax(atoms, fmax=1.0e-5)
        write(str(self.host.relaxed_poscar_path), atoms)
        print(f"  - Relaxed structure saved to {self.host.relaxed_poscar_path}")
        return atoms
