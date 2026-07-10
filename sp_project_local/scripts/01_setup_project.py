#!/usr/bin/env python3
from pathlib import Path

from sp_utils import PROJECT_HOME, PROJECT_SCRATCH, PROJECT_WORK, RAW_ROOT, ensure_project_dirs, write_markdown


def main() -> None:
    ensure_project_dirs()
    env_path = PROJECT_HOME / "config" / "paths.env"
    env_path.write_text(
        "\n".join(
            [
                f"SP_PROJECT_HOME={PROJECT_HOME}",
                f"SP_PROJECT_WORK={PROJECT_WORK}",
                f"SP_PROJECT_SCRATCH={PROJECT_SCRATCH}",
                f"SP_RAW={RAW_ROOT}",
                "SP_LSF_PROJECT=acc_DiseaseGeneCell",
                "SP_CONDA_ENV=omicverse",
                "SP_DEFAULT_RESOLUTION=square_008um",
                "SP_PILOT_SAMPLE=SC000895-R4",
                "",
            ]
        )
    )
    checks = []
    for path in [PROJECT_HOME, PROJECT_WORK, PROJECT_SCRATCH]:
        test = path / ".write_test"
        test.write_text("ok\n")
        test.unlink()
        checks.append(f"- writable: `{path}`")
    checks.append(f"- raw root exists: `{RAW_ROOT}` = {RAW_ROOT.exists()}")
    write_markdown(PROJECT_HOME / "results" / "setup_report.md", ["# Setup Report", "", *checks])
    print(f"Project initialized at {PROJECT_HOME}")
    print(f"Config written to {env_path}")


if __name__ == "__main__":
    main()
