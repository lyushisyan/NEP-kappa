from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from nepkappa.cli import main


def test_info_command(capsys):
    assert main(["info", "examples-input/input_7-film-nep-hiphive-rta.yaml"]) == 0

    output = capsys.readouterr().out
    assert "potentials/Si_NWs_XuKe.txt" in output
    assert "n_structures" in output
