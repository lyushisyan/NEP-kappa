from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from nepkappa.config import format_compare_config, format_config, parse_compare_args, parse_workflow_args


def test_bulk_yaml_uses_finite_displacement():
    args = parse_workflow_args("examples/1-bulk-nep-rta.yaml")

    assert args.poscar == "examples/POSCAR_bulk"
    assert args.nep_model == "potentials/Si_Bulk_Fan.txt"
    assert args.calculator == "nep"
    assert args.dimensionality == 3
    assert args.do_relax is True
    assert args.dim is None
    assert args.dim_fc2 == [3, 3, 3]
    assert args.dim_fc3 == [3, 3, 3]
    assert args.mesh == [21, 21, 21]
    assert args.use_hiphive is False
    assert args.isotope is False
    assert args.bfmp == 1.0e6
    assert args.plot_layout == "both"
    assert args.plot_path == "seekpath"
    assert args.plot_tau == "total"
    assert args.plot_kappa == "all"
    assert args.plot_temperature == 300


def test_legacy_dim_sets_fc2_and_fc3(tmp_path):
    config = tmp_path / "legacy-dim.yaml"
    config.write_text(
        """
structure:
  poscar: examples/POSCAR_bulk
calculator:
  name: nep
  nep_model: potentials/Si_Bulk_Fan.txt
force-constant:
  dim: [2, 3, 4]
""",
        encoding="utf-8",
    )

    args = parse_workflow_args(config)

    assert args.dim == [2, 3, 4]
    assert args.dim_fc2 == [2, 3, 4]
    assert args.dim_fc3 == [2, 3, 4]


def test_film_yaml_uses_hiphive():
    args = parse_workflow_args("examples/8-film-nep-hiphive-rta.yaml")

    assert args.poscar == "examples/POSCAR_film"
    assert args.dimensionality == 2
    assert args.effective_thickness == 10.0
    assert args.vacuum_axis == "z"
    assert args.nep_model == "potentials/Si_NWs_XuKe.txt"
    assert args.calculator == "nep"
    assert args.do_relax is False
    assert args.use_hiphive is True
    assert args.n_structures == 500
    assert args.cutoffs == [5.0, 4.0]


def test_nanowire_geometry_requires_effective_area(tmp_path):
    config = tmp_path / "wire.yaml"
    config.write_text(
        """
structure:
  poscar: examples/POSCAR_bulk
  dimensionality: 1
calculator:
  name: nep
  nep_model: potentials/Si_Bulk_Fan.txt
""",
        encoding="utf-8",
    )

    try:
        parse_workflow_args(config)
    except SystemExit as exc:
        assert exc.code == 2
    else:
        raise AssertionError("Expected parser to reject missing effective_area")


def test_vasp_yaml_uses_vasp_calculator():
    args = parse_workflow_args("examples/5-bulk-vasp-rta.yaml")

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
    args = parse_workflow_args("examples/1-bulk-nep-rta.yaml")
    summary = format_config(args)

    assert "use_hiphive" in summary
    assert "n_structures" not in summary
    assert "rattle_std" not in summary
    assert "cutoffs" not in summary
    assert "min_dist" not in summary


def test_compare_yaml_parses_dft_and_nep_dirs():
    args = parse_compare_args("examples/compare.yaml")

    assert args.dft_dir == "results/dft"
    assert args.nep_dir == "results/nep"
    assert args.compare_dir == "comparison"
    assert args.reference_label == "DFT"
    assert args.candidate_label == "NEP"
    assert args.plot_layout == "both"
    assert args.plot_path == "seekpath"
    assert args.plot_tau == "total"
    assert args.plot_kappa == "all"

    summary = format_compare_config(args)
    assert "dft_dir" in summary
    assert "nep_dir" in summary
