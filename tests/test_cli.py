from pathlib import Path
from io import StringIO
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from nepkappa.cli import Tee, build_parser, main


def test_force_constant_commands_parse():
    parser = build_parser()

    assert parser.parse_args(["fc2", "input.yaml"]).command == "fc2"
    assert parser.parse_args(["fc2fc3", "input.yaml"]).command == "fc2fc3"
    assert parser.parse_args(["fc", "input.yaml"]).command == "fc"


def test_info_command(capsys):
    assert main(["info", "examples/8-film-nep-hiphive-rta.yaml"]) == 0

    output = capsys.readouterr().out
    assert "potentials/Si_NWs_XuKe.txt" in output
    assert "n_structures" in output


def test_tee_compacts_carriage_return_updates_in_log():
    terminal = StringIO()
    log = StringIO()
    tee = Tee(terminal, log)

    tee.write("start\n")
    tee.write("\rProducing FC3: 00:00 elapsed")
    tee.write("\rProducing FC3: 00:01 elapsed")
    tee.write("\rProducing FC3: 00:02 elapsed\n")
    tee.write("done\n")

    assert terminal.getvalue() == (
        "start\n"
        "\rProducing FC3: 00:00 elapsed"
        "\rProducing FC3: 00:01 elapsed"
        "\rProducing FC3: 00:02 elapsed\n"
        "done\n"
    )
    assert log.getvalue() == "start\nProducing FC3: 00:02 elapsed\ndone\n"
