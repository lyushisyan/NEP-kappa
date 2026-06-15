from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from nepkappa.config import format_config, parse_workflow_args


def test_bulk_yaml_uses_finite_displacement():
    args = parse_workflow_args("examples-input/input_1-bulk-nep-rta.yaml")

    assert args.poscar == "examples-input/POSCAR_bulk"
    assert args.nep_model == "potentials/Si_Bulk_Fan.txt"
    assert args.calculator == "nep"
    assert args.do_relax is True
    assert args.dim == [3, 3, 3]
    assert args.mesh == [21, 21, 21]
    assert args.fc_calculator == "symfc"
    assert args.use_hiphive is False


def test_film_yaml_uses_hiphive():
    args = parse_workflow_args("examples-input/input_7-film-nep-hiphive-rta.yaml")

    assert args.poscar == "examples-input/POSCAR_film"
    assert args.nep_model == "potentials/Si_NWs_XuKe.txt"
    assert args.calculator == "nep"
    assert args.do_relax is False
    assert args.use_hiphive is True
    assert args.n_structures == 500
    assert args.cutoffs == [5.0, 4.0]


def test_vasp_yaml_uses_vasp_calculator():
    args = parse_workflow_args("examples-input/input_5-bulk-vasp-rta.yaml")

    assert args.calculator == "vasp"
    assert args.nep_model is None
    assert args.do_relax is True
    assert args.vasp_command == (
        "env -u DISPLAY -u XAUTHORITY mpirun -np 24 "
        "/mnt/DATA/home/sxliu98/vasp/vasp.6.4.3/bin/vasp_std"
    )
    assert args.vasp_path == "/mnt/DATA/home/sxliu98/vasp/vasp.6.4.3/bin/vasp_std"
    assert args.potcar_path == "/mnt/DATA/home/sxliu98/vasp/potpaw_PBE.64"
    assert args.vasp_workdir == "vasp-runs"
    assert args.vasp_relax_workdir == "vasp-relax"
    assert args.vasp_kwargs["encut"] == 650
    assert args.vasp_kwargs["kspacing"] == 0.2
    assert args.vasp_kwargs["kgamma"] is True
    assert args.vasp_kwargs["kpar"] == 2
    assert args.vasp_kwargs["ncore"] == 3
    assert args.vasp_kwargs["lreal"] is False
    assert args.vasp_relax_stages["coarse"]["ediffg"] == -0.05
    assert args.vasp_relax_stages["fine"]["ediffg"] == -0.01


def test_finite_displacement_display_hides_hiphive_defaults():
    args = parse_workflow_args("examples-input/input_1-bulk-nep-rta.yaml")
    summary = format_config(args)

    assert "use_hiphive" in summary
    assert "n_structures" not in summary
    assert "rattle_std" not in summary
    assert "cutoffs" not in summary
    assert "min_dist" not in summary
